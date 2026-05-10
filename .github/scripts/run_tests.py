"""Unittest driver that emits a GitHub Actions step-summary.

Runs ``unittest.TestLoader().discover(...)`` against a tests
directory, prints the usual verbose output to stderr so CI logs
remain useful, and writes a markdown summary to
``$GITHUB_STEP_SUMMARY`` when the workflow exports that variable.
Locally (no ``$GITHUB_STEP_SUMMARY`` set) the same summary is
printed to stderr after the run so a developer running the
script directly sees the same at-a-glance status.

Used by ``.github/workflows/tests.yml`` to make the workflow run
page show a per-project pass / fail / skip table with collapsible
details for failures, errors, and skips. Replaces a plain
``python -m unittest discover -s tests -v`` invocation which
landed only in the raw log.

Usage::

    python .github/scripts/run_tests.py <project-dir> [test-spec]

``<project-dir>`` is the directory containing the project's
``tests/`` subdirectory. ``[test-spec]`` is an optional dotted-name
restriction (e.g. ``tests.test_config.PydanticModelTest``); when
omitted the full test discovery runs.

Exit code: 0 if every test passed (skips are fine), 1 otherwise —
matches ``python -m unittest`` semantics.
"""

from __future__ import annotations

import os
import sys
import time
import unittest
from pathlib import Path


def _md_escape(s: str) -> str:
    """Escape characters that confuse the markdown table renderer."""
    return s.replace("|", "\\|").replace("\n", " ")


def _details_block(summary: str, body: str) -> list[str]:
    """Render a ``<details><summary>``...``</summary>...</details>`` block.

    GitHub's markdown renderer recognises raw HTML inside markdown
    files; the workflow summary page does too. Collapsible sections
    keep the summary scannable when there are dozens of skips.
    """
    return [
        "<details>",
        f"<summary>{summary}</summary>",
        "",
        body,
        "</details>",
        "",
    ]


def render_summary(
    result: unittest.TestResult,
    duration_s: float,
    project_label: str,
) -> str:
    """Build the markdown summary for one unittest run."""
    n_total = result.testsRun
    n_fail = len(result.failures)
    n_err = len(result.errors)
    n_skip = len(result.skipped)
    n_pass = n_total - n_fail - n_err - n_skip
    ok_icon = "✅" if (n_fail == 0 and n_err == 0) else "❌"

    lines: list[str] = [
        f"## {ok_icon} Test summary — `{project_label}`",
        "",
        "| Status | Count |",
        "|---|---:|",
        f"| Passed | {n_pass} |",
        f"| Failed | {n_fail} |",
        f"| Errored | {n_err} |",
        f"| Skipped | {n_skip} |",
        f"| **Total** | **{n_total}** |",
        "",
        f"Duration: {duration_s:.2f} s",
        "",
    ]

    if result.failures:
        lines.append("### Failures")
        lines.append("")
        for test, tb in result.failures:
            lines.extend(_details_block(
                f"<code>{_md_escape(str(test))}</code>",
                "```\n" + tb.rstrip() + "\n```",
            ))

    if result.errors:
        lines.append("### Errors")
        lines.append("")
        for test, tb in result.errors:
            lines.extend(_details_block(
                f"<code>{_md_escape(str(test))}</code>",
                "```\n" + tb.rstrip() + "\n```",
            ))

    if result.skipped:
        body = "\n".join(
            f"- `{_md_escape(str(test))}` — {_md_escape(reason)}"
            for test, reason in result.skipped
        )
        lines.extend(_details_block(
            f"Skipped tests ({len(result.skipped)})",
            body,
        ))

    return "\n".join(lines)


def _build_suite(project_dir: Path, test_spec: str | None):
    """Load the test suite. With a spec, load that dotted name; without,
    discover everything under ``project_dir/tests``.

    ``top_level_dir`` is left as the default (= ``start_dir``) so the
    discovery works whether or not ``tests/`` has an ``__init__.py``
    — synthetic_people's does, nextflow_pipeline's doesn't.
    The project's import root is already on ``sys.path`` by this
    point, so test-side ``from <package>.x import y`` still resolves.
    """
    loader = unittest.TestLoader()
    if test_spec:
        return loader.loadTestsFromName(test_spec)
    tests_dir = project_dir / "tests"
    if not tests_dir.is_dir():
        sys.exit(f"no tests directory at {tests_dir}")
    return loader.discover(start_dir=str(tests_dir))


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        sys.exit(
            "usage: run_tests.py <project-dir> [test-spec]\n"
            "  project-dir: directory containing tests/ "
            "(e.g. synthetic_people)\n"
            "  test-spec:   optional unittest dotted name "
            "(e.g. tests.test_config.PydanticModelTest)"
        )

    project_dir = Path(argv[1]).resolve()
    if not project_dir.is_dir():
        sys.exit(f"project dir does not exist: {project_dir}")
    test_spec = argv[2] if len(argv) >= 3 else None

    # Make the project importable; tests historically rely on
    # ``sys.path`` containing the project root.
    sys.path.insert(0, str(project_dir))

    suite = _build_suite(project_dir, test_spec)
    runner = unittest.TextTestRunner(verbosity=2, stream=sys.stderr)
    t0 = time.monotonic()
    result = runner.run(suite)
    duration_s = time.monotonic() - t0

    summary = render_summary(result, duration_s, project_dir.name)
    step_summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if step_summary:
        with open(step_summary, "a", encoding="utf-8") as fh:
            fh.write(summary)
            fh.write("\n")
    else:
        sys.stderr.write("\n")
        sys.stderr.write(summary)
        sys.stderr.write("\n")

    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
