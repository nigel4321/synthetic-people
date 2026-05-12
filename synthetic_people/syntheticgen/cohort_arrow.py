"""Arrow IPC streaming intermediate for the cohort phase (Phase 5d.1).

Replaces the in-memory sites-list hand-off between parent and BCF
workers with a memory-mapped Arrow IPC file. The hand-off used to fail
at n=3000+ through CPython refcount-COW divergence (see
``PERFORMANCE_PLAN.md`` §5e); this module is the load-bearing change
that unblocks the path to n=1M.

Flow at a glance:

* Parent walks the tree once, applies overlays, and feeds site dicts
  to :func:`write_arrow_file`. Each batch (default
  ``batch_size=256``) is densified into a ``(batch_size, 2 *
  n_samples)`` int8 haplotype matrix and flushed via
  ``pyarrow.ipc.new_file().write_batch(...)``. Parent peak RSS during
  write is bounded by the empirical fit
  ``~9.5 bytes/element × batch_size × n_haplotypes`` (Spike 2b).
* Workers ``pa.memory_map`` the file and iterate via
  :func:`read_arrow_slice`, slicing ``[2*sample_lo : 2*sample_hi]``
  columns from each batch's genotype matrix as a zero-copy numpy view.
  No refcount-COW divergence because numpy holds one PyObject wrapper
  for the whole mmap'd array, not one per element (validated end-to-
  end by Spike 2 on 2026-05-09).

Schema:

* Per-row: ``pos`` int64 / ``id`` string / ``ref`` string / ``alts``
  list<string> / ``acs`` list<int32> / ``afs`` list<float32> /
  ``genotypes`` fixed_size_list<int8, 2*n_samples> / nullable
  INFO-overlay fields (``clnsig``, ``clndn``, ``cosmic_id``,
  ``cosmic_gene``, ``svtype``, ``svlen``, ``end``, ``cipos_lo``,
  ``cipos_hi``).
* Schema metadata (constant per file): ``chrom``, ``n_samples``,
  ``n_haplotypes``, ``format_version``.

Conditional pyarrow dep. Top-level functions raise
``ImportError`` with an install hint if pyarrow is missing — the
``cohort-mode sites_list`` (Phase B) path remains available without
pyarrow for n ≤ 30k.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Iterator

import numpy as np


DEFAULT_BATCH_SIZE = 256


META_KEY_CHROM = b"synthetic_people:chrom"
META_KEY_N_SAMPLES = b"synthetic_people:n_samples"
META_KEY_N_HAPLOTYPES = b"synthetic_people:n_haplotypes"
META_KEY_FORMAT_VERSION = b"synthetic_people:format_version"
FORMAT_VERSION = "1"


def _require_pyarrow():
    try:
        import pyarrow  # noqa: F401
        import pyarrow.ipc  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "Phase 5d.1 cohort-arrow path requires pyarrow. Install with "
            "`pip install pyarrow`, or pass `--cohort-mode sites_list` "
            "to stay on the Phase B path (supported up to n ~= 30000)."
        ) from exc


def cohort_schema(n_samples: int, chrom: str):
    """Build the Arrow IPC schema for one chromosome's cohort intermediate.

    Schema metadata captures ``chrom`` / ``n_samples`` /
    ``n_haplotypes`` / ``format_version`` so workers don't need
    those passed alongside the file path.
    """
    _require_pyarrow()
    import pyarrow as pa

    n_haplotypes = 2 * n_samples
    fields = [
        pa.field("pos", pa.int64(), nullable=False),
        pa.field("id", pa.string(), nullable=False),
        pa.field("ref", pa.string(), nullable=False),
        pa.field("alts", pa.list_(pa.string()), nullable=False),
        pa.field("acs", pa.list_(pa.int32()), nullable=False),
        pa.field("afs", pa.list_(pa.float32()), nullable=False),
        pa.field(
            "genotypes",
            pa.list_(pa.int8(), n_haplotypes),
            nullable=False,
        ),
        pa.field("clnsig", pa.string(), nullable=True),
        pa.field("clndn", pa.string(), nullable=True),
        pa.field("cosmic_id", pa.string(), nullable=True),
        pa.field("cosmic_gene", pa.string(), nullable=True),
        pa.field("svtype", pa.string(), nullable=True),
        pa.field("svlen", pa.int64(), nullable=True),
        pa.field("end", pa.int64(), nullable=True),
        pa.field("cipos_lo", pa.int32(), nullable=True),
        pa.field("cipos_hi", pa.int32(), nullable=True),
    ]
    metadata = {
        META_KEY_CHROM: chrom.encode(),
        META_KEY_N_SAMPLES: str(n_samples).encode(),
        META_KEY_N_HAPLOTYPES: str(n_haplotypes).encode(),
        META_KEY_FORMAT_VERSION: FORMAT_VERSION.encode(),
    }
    return pa.schema(fields, metadata=metadata)


def _densify_carriers_to_matrix(
    sites: list,
    n_haplotypes: int,
) -> np.ndarray:
    """Build a ``(len(sites), n_haplotypes)`` int8 haplotype matrix
    from a list of site dicts. Accepts either sparse ``carriers`` or
    dense ``gts`` shape — matches the dual support already in
    ``CohortBcfWriter.write_site``.
    """
    n = len(sites)
    matrix = np.zeros((n, n_haplotypes), dtype=np.int8)
    for i, s in enumerate(sites):
        carriers = s.get("carriers")
        if carriers is not None:
            for hap_idx, allele_idx in carriers:
                matrix[i, hap_idx] = allele_idx
            continue
        gts = s.get("gts")
        if gts is not None:
            from .cohort_sites import carriers_from_dense_gts
            for hap_idx, allele_idx in carriers_from_dense_gts(gts):
                matrix[i, hap_idx] = allele_idx
        # else: all-zero row (valid — site recorded without genotypes)
    return matrix


def _build_batch(sites: list, n_samples: int, schema):
    """Construct one ``pyarrow.RecordBatch`` from a list of site dicts.

    ``schema`` must be the writer's declared schema (from
    :func:`cohort_schema`). The batch is built with this schema
    explicitly — without it, ``RecordBatch.from_pydict`` infers a
    schema from the data that doesn't match the writer's (different
    nullability, missing schema metadata, possibly different integer
    widths), and the IPC writer rejects the batch with
    ``ArrowInvalid: Tried to write record batch with different
    schema``.
    """
    import pyarrow as pa

    n_haplotypes = 2 * n_samples
    gt_matrix = _densify_carriers_to_matrix(sites, n_haplotypes)
    gt_flat = gt_matrix.reshape(-1)

    gt_values = pa.array(gt_flat, type=pa.int8())
    gt_arr = pa.FixedSizeListArray.from_arrays(gt_values, n_haplotypes)

    columns = {
        "pos": pa.array([s["pos"] for s in sites], type=pa.int64()),
        "id": pa.array([s.get("id") or "." for s in sites], type=pa.string()),
        "ref": pa.array([s["ref"] for s in sites], type=pa.string()),
        "alts": pa.array(
            [list(s["alts"]) for s in sites], type=pa.list_(pa.string())
        ),
        "acs": pa.array(
            [list(s.get("acs") or [0] * len(s["alts"])) for s in sites],
            type=pa.list_(pa.int32()),
        ),
        "afs": pa.array(
            [list(s.get("afs") or [0.0] * len(s["alts"])) for s in sites],
            type=pa.list_(pa.float32()),
        ),
        "genotypes": gt_arr,
        "clnsig": pa.array(
            [s.get("clnsig") for s in sites], type=pa.string()
        ),
        "clndn": pa.array(
            [s.get("clndn") for s in sites], type=pa.string()
        ),
        "cosmic_id": pa.array(
            [s.get("cosmic_id") for s in sites], type=pa.string()
        ),
        "cosmic_gene": pa.array(
            [s.get("cosmic_gene") for s in sites], type=pa.string()
        ),
        "svtype": pa.array(
            [s.get("svtype") for s in sites], type=pa.string()
        ),
        "svlen": pa.array(
            [s.get("svlen") for s in sites], type=pa.int64()
        ),
        "end": pa.array(
            [s.get("end") for s in sites], type=pa.int64()
        ),
        "cipos_lo": pa.array(
            [s["cipos"][0] if s.get("cipos") else None for s in sites],
            type=pa.int32(),
        ),
        "cipos_hi": pa.array(
            [s["cipos"][1] if s.get("cipos") else None for s in sites],
            type=pa.int32(),
        ),
    }
    return pa.RecordBatch.from_pydict(columns, schema=schema)


def stream_sites_to_arrow_batches(
    sites_iter: Iterable[dict],
    n_samples: int,
    batch_size: int = DEFAULT_BATCH_SIZE,
    schema=None,
) -> Iterator:
    """Yield one ``pyarrow.RecordBatch`` per ``batch_size`` sites.

    Sparse-to-dense haplotype densification happens here, per-batch.
    Memory bound is governed by the per-batch matrix:
    ``~9.5 bytes/element × batch_size × n_haplotypes`` (Spike 2b
    empirical fit). At ``batch_size=256``, ``n_samples=1_000_000`` the
    predicted parent peak is ~5 GB.

    ``schema`` should be the writer's declared schema (from
    :func:`cohort_schema`). When omitted, a schema is built locally
    with an empty ``chrom`` — fine for callers that don't write the
    batches to an IPC file. When passed to a downstream
    ``ipc.new_file(...).write_batch(...)`` writer the same schema
    must be used both there and here, otherwise the writer rejects
    the batch with ``ArrowInvalid: Tried to write record batch with
    different schema``.
    """
    _require_pyarrow()
    if batch_size <= 0:
        raise ValueError(f"batch_size must be positive, got {batch_size}")

    if schema is None:
        schema = cohort_schema(n_samples, chrom="")
    buffer: list = []
    for site in sites_iter:
        buffer.append(site)
        if len(buffer) >= batch_size:
            yield _build_batch(buffer, n_samples, schema)
            buffer = []
    if buffer:
        yield _build_batch(buffer, n_samples, schema)


def write_arrow_file(
    arrow_path: Path,
    chrom: str,
    n_samples: int,
    sites_iter: Iterable[dict],
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> int:
    """Stream ``sites_iter`` to an Arrow IPC file at ``arrow_path``.

    Returns the number of sites written. Existing files at
    ``arrow_path`` are overwritten (Arrow IPC writers create fresh
    files; there's no append mode for the file format).
    """
    _require_pyarrow()
    import pyarrow.ipc as ipc

    schema = cohort_schema(n_samples, chrom)
    n_written = 0
    arrow_path = Path(arrow_path)
    arrow_path.parent.mkdir(parents=True, exist_ok=True)

    with ipc.new_file(str(arrow_path), schema) as writer:
        for batch in stream_sites_to_arrow_batches(
            sites_iter, n_samples, batch_size=batch_size, schema=schema,
        ):
            writer.write_batch(batch)
            n_written += batch.num_rows
    return n_written


def read_arrow_metadata(arrow_path: Path) -> dict:
    """Read a cohort Arrow file's schema metadata + batch shape.

    Returns a dict with: ``chrom``, ``n_samples``, ``n_haplotypes``,
    ``format_version``, ``num_record_batches``, ``num_rows``.
    """
    _require_pyarrow()
    import pyarrow as pa
    import pyarrow.ipc as ipc

    with pa.memory_map(str(arrow_path), "r") as mm:
        reader = ipc.open_file(mm)
        meta = reader.schema.metadata or {}
        num_batches = reader.num_record_batches
        num_rows = sum(
            reader.get_batch(i).num_rows for i in range(num_batches)
        )
        return {
            "chrom": (meta.get(META_KEY_CHROM) or b"").decode(),
            "n_samples": int(
                (meta.get(META_KEY_N_SAMPLES) or b"0").decode()
            ),
            "n_haplotypes": int(
                (meta.get(META_KEY_N_HAPLOTYPES) or b"0").decode()
            ),
            "format_version": (
                meta.get(META_KEY_FORMAT_VERSION) or b""
            ).decode(),
            "num_record_batches": num_batches,
            "num_rows": num_rows,
        }


def read_arrow_slice(
    arrow_path: Path,
    sample_lo: int,
    sample_hi: int,
) -> Iterator[dict]:
    """Iterate site dicts for sample slice ``[sample_lo, sample_hi)``.

    Yields dicts shaped exactly like the in-memory site dict (so they
    plug into ``CohortBcfWriter.write_site`` without changes), but
    populated with ``gts: ["a|b", ...]`` of length
    ``sample_hi - sample_lo`` instead of sparse ``carriers``. The
    BCF writer's existing legacy ``gts`` path consumes them directly.

    Genotype reads are zero-copy numpy views into the mmap'd file —
    that's the design point, not an optimisation. The per-site
    Python loop is the realistic worker shape Spike 2 validated.
    """
    _require_pyarrow()
    import pyarrow as pa
    import pyarrow.ipc as ipc

    if sample_hi <= sample_lo:
        return

    with pa.memory_map(str(arrow_path), "r") as mm:
        reader = ipc.open_file(mm)
        meta = reader.schema.metadata or {}
        chrom = (meta.get(META_KEY_CHROM) or b"").decode()
        n_samples = int((meta.get(META_KEY_N_SAMPLES) or b"0").decode())
        n_haplotypes = 2 * n_samples

        if sample_hi > n_samples:
            raise ValueError(
                f"sample_hi {sample_hi} exceeds n_samples {n_samples} "
                f"in {arrow_path}"
            )
        if sample_lo < 0:
            raise ValueError(f"sample_lo must be >= 0, got {sample_lo}")

        hap_lo = 2 * sample_lo
        hap_hi = 2 * sample_hi
        slice_n = sample_hi - sample_lo

        for batch_idx in range(reader.num_record_batches):
            batch = reader.get_batch(batch_idx)
            n_rows = batch.num_rows

            gt_flat = batch.column("genotypes").values.to_numpy(
                zero_copy_only=True
            )
            gt_matrix = gt_flat.reshape(n_rows, n_haplotypes)

            pos_arr = batch.column("pos").to_numpy(zero_copy_only=True)
            ids = batch.column("id").to_pylist()
            refs = batch.column("ref").to_pylist()
            alts = batch.column("alts").to_pylist()
            acs = batch.column("acs").to_pylist()
            afs = batch.column("afs").to_pylist()
            clnsig = batch.column("clnsig").to_pylist()
            clndn = batch.column("clndn").to_pylist()
            cosmic_id = batch.column("cosmic_id").to_pylist()
            cosmic_gene = batch.column("cosmic_gene").to_pylist()
            svtype = batch.column("svtype").to_pylist()
            svlen = batch.column("svlen").to_pylist()
            end = batch.column("end").to_pylist()
            cipos_lo = batch.column("cipos_lo").to_pylist()
            cipos_hi = batch.column("cipos_hi").to_pylist()

            for i in range(n_rows):
                hap_slice = gt_matrix[i, hap_lo:hap_hi]
                gts = [
                    f"{hap_slice[2 * j]}|{hap_slice[2 * j + 1]}"
                    for j in range(slice_n)
                ]
                site: dict[str, Any] = {
                    "chrom": chrom,
                    "pos": int(pos_arr[i]),
                    "id": ids[i],
                    "ref": refs[i],
                    "alts": list(alts[i]),
                    "acs": list(acs[i]),
                    "afs": list(afs[i]),
                    "n_haplotypes": n_haplotypes,
                    "gts": gts,
                }
                if clnsig[i] is not None:
                    site["clnsig"] = clnsig[i]
                if clndn[i] is not None:
                    site["clndn"] = clndn[i]
                if cosmic_id[i] is not None:
                    site["cosmic_id"] = cosmic_id[i]
                if cosmic_gene[i] is not None:
                    site["cosmic_gene"] = cosmic_gene[i]
                if svtype[i] is not None:
                    site["svtype"] = svtype[i]
                if svlen[i] is not None:
                    site["svlen"] = svlen[i]
                if end[i] is not None:
                    site["end"] = end[i]
                if cipos_lo[i] is not None and cipos_hi[i] is not None:
                    site["cipos"] = (cipos_lo[i], cipos_hi[i])
                yield site


def read_arrow_carriers(arrow_path: Path, pos: int) -> np.ndarray:
    """Diagnostic helper: re-derive packed sparse carriers for one position.

    Returns ``np.ndarray`` of shape ``(n_carriers, 2)``, dtype
    ``np.int32`` — the same packed shape every other producer in
    the codebase emits. See ``cohort_sites.py`` module docstring.

    Replaces the post-hoc ``pickle.dump(sites_list)`` introspection
    that was possible under the in-memory sites-list path (called out
    in the §5d persistent-regressions tracker as item #2 — the on-disk
    Arrow file is the new source-of-truth).

    Scans batches in order; returns the first match (positions are
    expected to be unique within a chromosome's cohort file). Returns
    an empty ``(0, 2)`` array if the position is not found.
    """
    _require_pyarrow()
    import pyarrow as pa
    import pyarrow.ipc as ipc

    with pa.memory_map(str(arrow_path), "r") as mm:
        reader = ipc.open_file(mm)
        meta = reader.schema.metadata or {}
        n_haplotypes = int(
            (meta.get(META_KEY_N_HAPLOTYPES) or b"0").decode()
        )

        for batch_idx in range(reader.num_record_batches):
            batch = reader.get_batch(batch_idx)
            pos_arr = batch.column("pos").to_numpy(zero_copy_only=True)
            matches = np.where(pos_arr == pos)[0]
            if len(matches) == 0:
                continue
            row_idx = int(matches[0])
            gt_flat = batch.column("genotypes").values.to_numpy(
                zero_copy_only=True
            )
            gt_matrix = gt_flat.reshape(batch.num_rows, n_haplotypes)
            row = gt_matrix[row_idx]
            nonzero = np.flatnonzero(row)
            if not nonzero.size:
                return np.zeros((0, 2), dtype=np.int32)
            return np.column_stack(
                (nonzero, row[nonzero])
            ).astype(np.int32, copy=False)
    return np.zeros((0, 2), dtype=np.int32)
