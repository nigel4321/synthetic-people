"""Carriers spill-to-disk for the streaming-cohort safe-yield heap.

At n=1M, even packed carriers (post-PR-#71, ~8 bytes per row)
accumulate to ~60 KB per site averaging over the SFS. The safe-yield
heap in :func:`coalescent._stream_cohort_pass2` holds sites pending
position-sort; under realistic overlay-position distributions (real
dbSNP / ClinVar coordinates span the whole chromosome) the heap depth
approaches ``O(N_sites)``, giving ~60 GB parent peak RSS at WGS-scale.

This module spills the heavy ``carriers`` payload to a per-chromosome
sidecar file at heap-push time. Heap entries store only a small
``(offset, length)`` reference (~16 bytes plus dict overhead). On
heap-pop the carriers are read back and re-attached to the site dict
before it's yielded. Per-site bytes-in-RAM falls from ~60 KB to
~200 bytes — parent peak RSS becomes ``O(heap_depth × ~200 B)``
regardless of n.

Disk cost: roughly doubles the per-chrom Arrow scratch budget. The
:func:`cli._preflight_arrow_disk_check` accounts for it.

Output remains byte-identical to the no-spill path — the cross-mode
parity tests (``CohortModeArrowParityTest``,
``CohortModeArrowStreamingParityTest``) lock that in.

See ``PERFORMANCE_BUDGETS.md`` § "Known scaling ceiling" for the
empirical 2026-05-12 n=1M OOM that motivated this fix (Fix B.1 in
the carriers-scaling investigation).
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np


# Each packed carrier row is two int32 scalars (haplotype index +
# allele index) = 8 bytes. Hardcoded because changing the packed
# shape would break Arrow-file byte-parity across releases (the
# shape was established in PR #71 and is part of the on-disk
# contract).
_BYTES_PER_CARRIER_ROW = 8


class CarriersSidecar:
    """Per-chromosome append-only scratch file for spilled carriers.

    Lifecycle: one sidecar per chromosome. The file is opened on
    construction, written-to during the streaming-pass-2 walk, and
    unlinked unconditionally on ``close()`` — even if the streaming
    pass raised mid-chrom. The sidecar is transient scratch; nothing
    downstream depends on it once the chrom's Arrow file is fully
    written.

    Two file descriptors are held: one append-only write fd (Python
    file object) and one read-only fd (``os.open``). Reads use
    ``os.pread`` so they don't move the write fd's position. The
    write fd is unbuffered (``buffering=0``) so reads of just-written
    data via ``pread`` see them immediately — Python's default 8 KB
    write buffer would otherwise hide recent writes from the kernel.

    Single-writer, single-reader by construction (the streaming pass
    is parent-only; workers consume the Arrow file via mmap after
    the parent is done). No locking needed.

    Use as a context manager so cleanup runs on the exception path:

    ::

        with CarriersSidecar(scratch_dir / "carriers.chr22.spill") as sc:
            offset, n_bytes = sc.write(carriers_array)
            carriers_back = sc.read(offset, n_bytes)
            # ... use ...
        # File unlinked on exit.
    """

    def __init__(self, path: Path):
        self._path = Path(path)
        self._path.parent.mkdir(parents=True, exist_ok=True)
        self._write_fh = None
        self._read_fd = -1
        self._offset = 0
        self._closed = False
        # Exception-safe construction: if the second fd-open fails
        # after the first succeeds, the partially-acquired write
        # handle would leak AND the file would linger on disk.
        # Catch the failure, undo prior state, then re-raise.
        try:
            # Unbuffered writes (``buffering=0``) so reads via
            # ``os.pread`` see them immediately — Python's default
            # 8 KB write buffer would otherwise hide recent writes.
            self._write_fh = open(self._path, "wb", buffering=0)
            # Separate read fd: ``os.pread`` is positional and does
            # not touch the read fd's offset, so writes (which
            # advance the write fd's offset) and reads don't
            # interfere.
            self._read_fd = os.open(str(self._path), os.O_RDONLY)
        except OSError:
            if self._write_fh is not None:
                try:
                    self._write_fh.close()
                except OSError:
                    pass
                self._write_fh = None
            try:
                self._path.unlink(missing_ok=True)
            except OSError:
                pass
            raise

    def write(self, carriers: np.ndarray) -> tuple[int, int]:
        """Append packed carriers; return ``(offset, n_bytes)``.

        ``carriers`` must be a 2-D ``np.int32`` array of shape
        ``(n_rows, 2)`` — the post-PR-#71 packed shape. Empty
        carriers (shape ``(0, 2)``) write 0 bytes and return
        ``(current_offset, 0)`` so the round-trip is symmetric.

        Raises ``TypeError`` for wrong dtype, ``ValueError`` for
        wrong shape, and ``OSError`` on a short write. These are
        real exceptions (not ``assert``) so they fire under
        ``python -O`` too — the sidecar's byte layout depends on
        the input shape and a silent-corruption mode is not OK.
        """
        if carriers.dtype != np.int32:
            raise TypeError(
                f"carriers must be np.int32; got {carriers.dtype}",
            )
        if carriers.ndim != 2 or carriers.shape[1] != 2:
            raise ValueError(
                f"carriers must be shape (N, 2); got {carriers.shape}",
            )
        if carriers.shape[0] == 0:
            return (self._offset, 0)
        # Ensure C-contiguous so ``tobytes`` is a single block.
        if not carriers.flags["C_CONTIGUOUS"]:
            carriers = np.ascontiguousarray(carriers)
        data = carriers.tobytes()
        n_expected = len(data)
        # Unbuffered file.write(bytes) on a raw binary file returns
        # the bytes actually written, which may be less than
        # requested on EINTR / disk-full / signal interruption.
        # Loop until everything's flushed, or raise on a true
        # zero-byte write (which would loop forever).
        view = memoryview(data)
        n_written = 0
        while n_written < n_expected:
            chunk = self._write_fh.write(view[n_written:])
            if chunk is None:
                # ``buffering=0`` should never return None per
                # PEP 3116, but defensive.
                raise OSError("sidecar write returned None")
            if chunk == 0:
                raise OSError(
                    f"sidecar write made no progress: "
                    f"wrote {n_written}/{n_expected} bytes "
                    f"before stalling",
                )
            n_written += chunk
        offset = self._offset
        self._offset += n_expected
        return (offset, n_expected)

    def read(self, offset: int, n_bytes: int) -> np.ndarray:
        """Read packed carriers at ``offset``; return a fresh
        ``(n_rows, 2)`` int32 array.

        Returns an empty ``(0, 2)`` array when ``n_bytes == 0``
        (the no-carriers case, where the original site had no
        non-zero haplotypes).

        Raises ``ValueError`` if ``n_bytes`` isn't a multiple of
        the per-row size (8 = 2 × int32). Raises ``OSError`` on
        a short read (would indicate file corruption or a
        programmer bug in offset/length tracking).
        """
        if n_bytes == 0:
            return np.zeros((0, 2), dtype=np.int32)
        if n_bytes % _BYTES_PER_CARRIER_ROW != 0:
            raise ValueError(
                f"sidecar read length {n_bytes} is not a multiple "
                f"of the per-row size ({_BYTES_PER_CARRIER_ROW} "
                f"bytes) — would corrupt the (N, 2) reshape",
            )
        # ``os.pread`` on Linux can return fewer bytes than
        # requested on EINTR / short read; loop until satisfied
        # or detect a real EOF (zero progress) as corruption.
        chunks: list = []
        n_remaining = n_bytes
        cur_offset = offset
        while n_remaining > 0:
            chunk = os.pread(self._read_fd, n_remaining, cur_offset)
            if not chunk:
                raise OSError(
                    f"sidecar short read at offset={offset}: "
                    f"requested {n_bytes} bytes, got "
                    f"{n_bytes - n_remaining}",
                )
            chunks.append(chunk)
            n_remaining -= len(chunk)
            cur_offset += len(chunk)
        data = chunks[0] if len(chunks) == 1 else b"".join(chunks)
        n_rows = n_bytes // _BYTES_PER_CARRIER_ROW
        # ``frombuffer`` returns a read-only view of the underlying
        # bytes; ``.copy()`` so the caller can mutate (e.g. workers
        # zero out portions during fan-out).
        return np.frombuffer(
            data, dtype=np.int32,
        ).reshape(n_rows, 2).copy()

    def close(self) -> None:
        """Close both fds and unlink the sidecar file.

        Idempotent — safe to call from a ``finally`` block even
        after an earlier ``close()`` or after a partial-
        construction failure (in which case one or both fds may
        not exist). Swallows ``OSError`` from partial cleanup
        paths; the goal is to leave the filesystem clean, not to
        fail loudly on a sidecar we were about to discard anyway.
        """
        if self._closed:
            return
        self._closed = True
        if self._write_fh is not None:
            try:
                self._write_fh.close()
            except OSError:
                pass
            self._write_fh = None
        if self._read_fd != -1:
            try:
                os.close(self._read_fd)
            except OSError:
                pass
            self._read_fd = -1
        try:
            self._path.unlink(missing_ok=True)
        except OSError:
            pass

    @property
    def path(self) -> Path:
        return self._path

    @property
    def n_bytes_written(self) -> int:
        return self._offset

    def __enter__(self) -> "CarriersSidecar":
        return self

    def __exit__(self, *exc) -> bool:
        self.close()
        return False
