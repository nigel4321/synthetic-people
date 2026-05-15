"""Tests for the YAML config-file layer.

Covers:

- The Pydantic models accept / reject inputs with clear error
  messages.
- ``discover_config_file`` finds the well-known filename in cwd only.
- ``load_and_validate_config`` repackages Pydantic errors into a
  single SystemExit with every problem listed at once.
- The CLI > config > defaults precedence in
  ``merge_config_into_args``.
- ``format_effective_values`` correctly tags sources and skips
  default-equal values.
- The committed JSON Schema file matches what
  ``generate_json_schema()`` produces (sync test).

Skipped cleanly if pydantic / PyYAML aren't installed.
"""

from __future__ import annotations

import io
import sys
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

try:
    import pydantic  # noqa: F401
    import yaml  # noqa: F401
    HAS_DEPS = True
except ImportError:
    HAS_DEPS = False


REPO_SCHEMA_PATH = (
    Path(__file__).resolve().parent.parent
    / "generate_people_config.schema.json"
)


@unittest.skipUnless(HAS_DEPS, "pydantic + PyYAML not installed")
class PydanticModelTest(unittest.TestCase):
    """Type / bounds / enum / cross-field validators on the models."""

    def setUp(self):
        from syntheticgen.config import _models
        self.Config = _models()

    def test_minimal_config_is_valid(self):
        c = self.Config(schema_version=1)
        # All defaults apply.
        self.assertEqual(c.cohort.n, 10)
        self.assertEqual(c.cohort.build, "GRCh38")
        self.assertEqual(c.performance.cohort_mode, "auto")
        self.assertEqual(c.performance.cohort_arrow_batch_size, 256)

    def test_schema_version_required(self):
        with self.assertRaises(pydantic.ValidationError) as ctx:
            self.Config()
        self.assertIn("schema_version", str(ctx.exception))

    def test_schema_version_must_match_supported(self):
        with self.assertRaises(pydantic.ValidationError) as ctx:
            self.Config(schema_version=999)
        msg = str(ctx.exception)
        self.assertIn("schema_version 999", msg)
        self.assertIn("expected 1", msg)

    def test_unknown_field_rejected(self):
        with self.assertRaises(pydantic.ValidationError) as ctx:
            self.Config(schema_version=1, cohort={"n": 10, "bogus": 5})
        self.assertIn("bogus", str(ctx.exception))

    def test_n_must_be_positive(self):
        with self.assertRaises(pydantic.ValidationError) as ctx:
            self.Config(schema_version=1, cohort={"n": 0})
        self.assertIn("greater than or equal to 1", str(ctx.exception))

    def test_cohort_mode_enum_validated(self):
        with self.assertRaises(pydantic.ValidationError) as ctx:
            self.Config(
                schema_version=1,
                performance={"cohort_mode": "arow"},  # typo
            )
        self.assertIn("cohort_mode", str(ctx.exception))

    def test_cohort_mode_accepts_every_cli_choice(self):
        # Regression: the pydantic Literal must stay in sync with
        # the cli's ``--cohort-mode`` choices. The original Literal
        # was missing ``arrow-streaming`` (added to the cli in
        # PR #52 but not the schema), which made every config file
        # carrying that value crash with a validation error before
        # the run could start.
        from syntheticgen.cli import _parser
        parser = _parser(Path(__file__).resolve().parent.parent)
        cli_choices = None
        for action in parser._actions:
            if action.dest == "cohort_mode":
                cli_choices = tuple(action.choices)
                break
        self.assertIsNotNone(
            cli_choices, "cli has no --cohort-mode action",
        )
        # Every cli choice must validate cleanly through pydantic.
        for choice in cli_choices:
            c = self.Config(
                schema_version=1,
                performance={"cohort_mode": choice},
            )
            self.assertEqual(c.performance.cohort_mode, choice)

    def test_inject_density_bounded_zero_to_one(self):
        with self.assertRaises(pydantic.ValidationError):
            self.Config(
                schema_version=1,
                overlays={"clinvar": {"inject_density": 1.5}},
            )

    def test_sv_length_max_must_be_geq_min(self):
        with self.assertRaises(pydantic.ValidationError) as ctx:
            self.Config(
                schema_version=1,
                structural_variants={"length_min": 100, "length_max": 50},
            )
        self.assertIn("length_max", str(ctx.exception))

    def test_admixture_fractions_must_sum_to_one_when_enabled(self):
        with self.assertRaises(pydantic.ValidationError) as ctx:
            self.Config(
                schema_version=1,
                admixture={
                    "enabled": True,
                    "eur_frac": 0.5,
                    "sas_frac": 0.3,
                    "afr_frac": 0.3,  # sum = 1.1
                },
            )
        self.assertIn("sum to 1.0", str(ctx.exception))

    def test_admixture_fractions_unchecked_when_disabled(self):
        # If admixture is off the fractions are ignored — useful so a
        # user with admixture: enabled=false can still ship a config
        # with leftover ancestry values from a previous run.
        c = self.Config(
            schema_version=1,
            admixture={
                "enabled": False,
                "eur_frac": 0.5, "sas_frac": 0.3, "afr_frac": 0.3,
            },
        )
        self.assertFalse(c.admixture.enabled)


@unittest.skipUnless(HAS_DEPS, "pydantic + PyYAML not installed")
class DiscoverConfigFileTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_returns_none_when_absent(self):
        from syntheticgen.config import discover_config_file
        self.assertIsNone(discover_config_file(self.dir))

    def test_returns_path_when_present(self):
        from syntheticgen.config import (
            DEFAULT_CONFIG_FILENAME,
            discover_config_file,
        )
        target = self.dir / DEFAULT_CONFIG_FILENAME
        target.write_text("schema_version: 1\n")
        result = discover_config_file(self.dir)
        self.assertEqual(result, target)

    def test_does_not_walk_up_to_parent(self):
        from syntheticgen.config import (
            DEFAULT_CONFIG_FILENAME,
            discover_config_file,
        )
        # Place the config in the parent directory; discovery from a
        # child directory must NOT find it. Predictability for users
        # running multiple jobs from different cwds.
        (self.dir / DEFAULT_CONFIG_FILENAME).write_text(
            "schema_version: 1\n",
        )
        child = self.dir / "subdir"
        child.mkdir()
        self.assertIsNone(discover_config_file(child))


@unittest.skipUnless(HAS_DEPS, "pydantic + PyYAML not installed")
class LoadAndValidateConfigTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.path = self.dir / "test_config.yaml"

    def tearDown(self):
        self.tmp.cleanup()

    def test_valid_yaml_loads(self):
        from syntheticgen.config import load_and_validate_config
        self.path.write_text(
            "schema_version: 1\n"
            "cohort:\n  n: 100\n  seed: 7\n"
        )
        c = load_and_validate_config(self.path)
        self.assertEqual(c.cohort.n, 100)
        self.assertEqual(c.cohort.seed, 7)

    def test_invalid_yaml_syntax_raises_systemexit_with_path(self):
        from syntheticgen.config import load_and_validate_config
        self.path.write_text("schema_version: 1\ncohort:\n  n: [unterminated")
        with self.assertRaises(SystemExit) as ctx:
            load_and_validate_config(self.path)
        msg = str(ctx.exception)
        self.assertIn(str(self.path), msg)
        self.assertIn("not valid YAML", msg)

    def test_non_mapping_top_level_raises(self):
        from syntheticgen.config import load_and_validate_config
        self.path.write_text("- just a list")
        with self.assertRaises(SystemExit) as ctx:
            load_and_validate_config(self.path)
        self.assertIn("must be a YAML mapping", str(ctx.exception))

    def test_validation_errors_listed_together(self):
        from syntheticgen.config import load_and_validate_config
        # Two distinct violations: bad n and bad cohort_mode.
        self.path.write_text(
            "schema_version: 1\n"
            "cohort:\n  n: -10\n"
            "performance:\n  cohort_mode: arow\n"
        )
        with self.assertRaises(SystemExit) as ctx:
            load_and_validate_config(self.path)
        msg = str(ctx.exception)
        self.assertIn("Config validation failed", msg)
        self.assertIn("cohort.n", msg)
        self.assertIn("performance.cohort_mode", msg)


@unittest.skipUnless(HAS_DEPS, "pydantic + PyYAML not installed")
class MergePrecedenceTest(unittest.TestCase):
    """The contract: CLI > config > defaults. Verified end-to-end via
    the real parser, since the merge logic depends on argparse dest
    names and the SUPPRESS-shadow-parse trick."""

    def setUp(self):
        from syntheticgen.cli import _parser
        self.parser = _parser(
            Path(__file__).resolve().parent.parent,
        )
        from syntheticgen.config import (
            load_and_validate_config,
            merge_config_into_args,
            parse_explicit_cli_args,
        )
        self.load = load_and_validate_config
        self.merge = merge_config_into_args
        self.explicit = parse_explicit_cli_args
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.path = self.dir / "c.yaml"
        self.path.write_text(
            "schema_version: 1\n"
            "cohort:\n  n: 3000\n  seed: 42\n"
            "performance:\n  cohort_mode: arrow\n  workers: 8\n"
        )

    def tearDown(self):
        self.tmp.cleanup()

    def _run(self, argv):
        args = self.parser.parse_args(argv)
        config = self.load(self.path)
        explicit = self.explicit(self.parser, argv)
        # Pass parser=parser so the merge applies argparse type
        # coercion to config-loaded values — without it, Path-typed
        # flags would stay as str (the bug this test class guards).
        merged = self.merge(args, config, explicit, parser=self.parser)
        return merged, explicit

    def test_config_value_used_when_no_cli(self):
        merged, _ = self._run([])
        self.assertEqual(merged.n, 3000)
        self.assertEqual(merged.seed, 42)
        self.assertEqual(merged.cohort_mode, "arrow")
        self.assertEqual(merged.workers, 8)

    def test_cli_overrides_config(self):
        merged, explicit = self._run(["--workers", "4", "--n", "100"])
        self.assertEqual(merged.workers, 4)
        self.assertEqual(merged.n, 100)
        # config-only values still applied
        self.assertEqual(merged.cohort_mode, "arrow")
        self.assertEqual(merged.seed, 42)
        # explicit_cli reports the dests, not the flags
        self.assertIn("workers", explicit)
        self.assertIn("n", explicit)
        self.assertNotIn("cohort_mode", explicit)

    def test_default_when_neither_cli_nor_config(self):
        # No --build in argv, no `cohort.build` in config -> argparse
        # default applies.
        merged, _ = self._run([])
        self.assertEqual(merged.build, "GRCh38")

    def test_male_fraction_default_when_neither_cli_nor_config(self):
        # M13.1: --male-fraction defaults to 0.5 when absent from
        # both the CLI and the config. Pins the cli > config >
        # defaults precedence for this field.
        merged, _ = self._run([])
        self.assertEqual(merged.male_fraction, 0.5)

    def test_male_fraction_cli_flag_sets_args(self):
        # The CLI flag should land in args.male_fraction as a float.
        merged, explicit = self._run(["--male-fraction", "0.2"])
        self.assertAlmostEqual(merged.male_fraction, 0.2)
        self.assertIn("male_fraction", explicit)

    def test_male_fraction_cli_overrides_config(self):
        # Write a config that sets male_fraction, then pass a
        # different value on the CLI — CLI wins.
        self.path.write_text(
            "schema_version: 1\n"
            "cohort:\n  n: 3000\n  seed: 42\n  male_fraction: 0.8\n"
        )
        merged, _ = self._run(["--male-fraction", "0.2"])
        self.assertAlmostEqual(merged.male_fraction, 0.2)

    def test_male_fraction_config_used_when_no_cli(self):
        self.path.write_text(
            "schema_version: 1\n"
            "cohort:\n  n: 3000\n  seed: 42\n  male_fraction: 0.8\n"
        )
        merged, _ = self._run([])
        self.assertAlmostEqual(merged.male_fraction, 0.8)

    def test_male_fraction_cli_rejects_out_of_range(self):
        # PR #98 review (Copilot): Pydantic enforces ge=0.0/le=1.0
        # on the YAML field, but the CLI's bare ``type=float`` was
        # accepting nonsense values that flow into _draw_sexes and
        # produce degenerate all-male / all-female / NaN cohorts.
        # The custom ``type=_male_fraction_arg`` callable mirrors
        # the Pydantic constraint.
        for bad in ("-0.1", "1.5", "2", "-1", "nan", "inf", "-inf"):
            with self.subTest(value=bad):
                with self.assertRaises(SystemExit):
                    self._run(["--male-fraction", bad])

    def test_male_fraction_cli_accepts_boundary_values(self):
        # Endpoints 0.0 and 1.0 are valid (ge / le, not gt / lt).
        for ok in ("0.0", "0", "1.0", "1", "0.5"):
            with self.subTest(value=ok):
                merged, _ = self._run(["--male-fraction", ok])
                self.assertAlmostEqual(merged.male_fraction, float(ok))

    def test_no_config_path_leaves_args_unchanged(self):
        # merge with config=None must be a no-op.
        from syntheticgen.config import merge_config_into_args
        args = self.parser.parse_args(["--n", "55"])
        merged = merge_config_into_args(args, None, {"n"})
        self.assertEqual(merged.n, 55)
        self.assertEqual(merged.workers, 0)  # argparse default

    def test_path_typed_flags_coerced_from_config(self):
        # Regression: previously --cache-dir / --output-dir /
        # --dbsnp-vcf / --cosmic-vcf / --profile-memory all have
        # argparse ``type=Path``, but the config-merge path stored
        # their values as ``str`` (YAML's natural type for a string).
        # Downstream code that called ``cache_dir.mkdir(...)`` then
        # crashed with ``AttributeError: 'str' object has no
        # attribute 'mkdir'``. Lock in the parser-aware coercion
        # so all five flags become Path objects regardless of source.
        self.path.write_text(
            "schema_version: 1\n"
            "output:\n"
            "  dir: /tmp/synthpeople_test_out\n"
            "  cache_dir: /tmp/synthpeople_test_cache\n"
            "overlays:\n"
            "  rsid:\n"
            "    vcf: /tmp/dbsnp_override.vcf.gz\n"
            "  cosmic:\n"
            "    vcf: /tmp/cosmic_override.vcf.gz\n"
            "performance:\n"
            "  profile_memory: /tmp/memprof.tsv\n"
        )
        merged, _ = self._run([])
        self.assertIsInstance(merged.output_dir, Path)
        self.assertIsInstance(merged.cache_dir, Path)
        self.assertIsInstance(merged.dbsnp_vcf, Path)
        self.assertIsInstance(merged.cosmic_vcf, Path)
        self.assertIsInstance(merged.profile_memory, Path)
        # Spot-check one value to confirm the string content survived.
        self.assertEqual(
            str(merged.cache_dir), "/tmp/synthpeople_test_cache",
        )

    def test_null_optional_path_stays_none(self):
        # Optional path flags (--dbsnp-vcf, --cosmic-vcf,
        # --profile-memory) all default to None. If the config
        # explicitly sets one to ``null``, the merge must NOT call
        # ``Path(None)`` (which raises). None should pass through.
        self.path.write_text(
            "schema_version: 1\n"
            "overlays:\n"
            "  rsid:\n"
            "    vcf: null\n"
            "performance:\n"
            "  profile_memory: null\n"
        )
        merged, _ = self._run([])
        self.assertIsNone(merged.dbsnp_vcf)
        self.assertIsNone(merged.profile_memory)

    def test_cli_path_still_wins_over_config_path(self):
        # CLI > config precedence must still hold for path-typed
        # flags after the coercion change. argparse already parsed
        # the CLI value as Path; coercion must not override it.
        self.path.write_text(
            "schema_version: 1\n"
            "output:\n"
            "  cache_dir: /tmp/from_config\n"
        )
        merged, explicit = self._run(["--cache-dir", "/tmp/from_cli"])
        self.assertEqual(str(merged.cache_dir), "/tmp/from_cli")
        self.assertIsInstance(merged.cache_dir, Path)
        self.assertIn("cache_dir", explicit)

    def test_realistic_mixed_cli_and_config_across_sections(self):
        # End-to-end mixed-sources scenario: a "stable baseline"
        # config covering several sections (cohort, simulation,
        # overlays, performance, output, sequencing_errors), with
        # a handful of CLI overrides for a one-off run. Verifies:
        #
        # 1. Every config-only field reaches the merged Namespace
        #    with the expected value AND the correct Python type
        #    (Path for path-typed flags, int/float/bool/str
        #    otherwise).
        # 2. CLI flags override their config counterparts cleanly,
        #    regardless of which section the field lives in.
        # 3. Fields touched by neither CLI nor config fall through
        #    to argparse defaults.
        # 4. ``explicit_cli`` reports the dests the user actually
        #    typed and nothing else.
        self.path.write_text(
            "schema_version: 1\n"
            "cohort:\n"
            "  n: 200\n"
            "  seed: 7\n"
            "  build: GRCh37\n"
            "  chromosomes: '21,22'\n"
            "  chr_length_mb: 5.0\n"
            "simulation:\n"
            "  demo_model: OutOfAfrica_3G09\n"
            "  population: YRI\n"
            "overlays:\n"
            "  clinvar:\n"
            "    inject_density: 0.02\n"
            "  rsid:\n"
            "    density: 0.3\n"
            "    vcf: /tmp/dbsnp.vcf.gz\n"
            "sequencing_errors:\n"
            "  gt_flip_rate: 0.002\n"
            "  dropout_rate: 0.001\n"
            "performance:\n"
            "  workers: 4\n"
            "  cohort_mode: arrow\n"
            "  fanout_batch_size: 8\n"
            "output:\n"
            "  dir: /tmp/baseline_out\n"
            "  mode: cohort\n"
        )
        merged, explicit = self._run([
            # Overrides from three different sections — cohort,
            # performance, output — plus one path-typed flag to
            # exercise CLI > config precedence on Paths.
            "--n", "500",
            "--workers", "16",
            "--output-dir", "/tmp/one_off_out",
        ])

        # CLI wins for the three overridden dests.
        self.assertEqual(merged.n, 500)
        self.assertEqual(merged.workers, 16)
        self.assertEqual(str(merged.output_dir), "/tmp/one_off_out")
        self.assertIsInstance(merged.output_dir, Path)

        # Config-only values land with the right values + types
        # across several sections.
        # cohort
        self.assertEqual(merged.seed, 7)
        self.assertEqual(merged.build, "GRCh37")
        self.assertEqual(merged.chromosomes, "21,22")
        self.assertEqual(merged.chr_length_mb, 5.0)
        # simulation
        self.assertEqual(merged.demo_model, "OutOfAfrica_3G09")
        self.assertEqual(merged.population, "YRI")
        # overlays
        self.assertEqual(merged.clinvar_inject_density, 0.02)
        self.assertEqual(merged.rsid_density, 0.3)
        self.assertIsInstance(merged.dbsnp_vcf, Path)
        self.assertEqual(str(merged.dbsnp_vcf), "/tmp/dbsnp.vcf.gz")
        # sequencing_errors
        self.assertEqual(merged.error_rate, 0.002)
        self.assertEqual(merged.dropout_rate, 0.001)
        # performance (workers overridden above — these are the
        # config-only siblings)
        self.assertEqual(merged.cohort_mode, "arrow")
        self.assertEqual(merged.fanout_batch_size, 8)
        # output
        self.assertEqual(merged.mode, "cohort")

        # Untouched fields fall through to argparse defaults.
        # ``somatic`` and ``svs_per_person`` are both absent from
        # this test's config and not passed on the CLI, so they
        # should land at the argparse defaults (False and 3
        # respectively, matching the pydantic model defaults too).
        self.assertFalse(merged.somatic)
        self.assertEqual(merged.svs_per_person, 3)

        # explicit_cli reports exactly the dests the user typed.
        self.assertEqual(
            explicit, {"n", "workers", "output_dir"},
        )

    def test_int_typed_flags_still_int_after_coercion(self):
        # Defensive: the coercion is generic over any callable
        # ``action.type``. Confirm an int-typed flag stays int
        # (pydantic delivers int already; coerce is the identity
        # int() call). Catches a regression where we accidentally
        # narrowed coercion to only Path.
        merged, _ = self._run([])
        self.assertIsInstance(merged.n, int)
        self.assertEqual(merged.n, 3000)
        self.assertIsInstance(merged.workers, int)

    def test_parse_explicit_cli_args_does_not_mutate_parser(self):
        # The shadow-parse must restore parser defaults so the parser
        # remains usable afterwards. Round-trip check.
        before = {a.dest: a.default for a in self.parser._actions}
        before_set_defaults = dict(self.parser._defaults)
        _ = self.explicit(self.parser, ["--n", "10"])
        after = {a.dest: a.default for a in self.parser._actions}
        after_set_defaults = dict(self.parser._defaults)
        self.assertEqual(before, after)
        self.assertEqual(before_set_defaults, after_set_defaults)
        # And the real parse still produces the same defaults.
        args = self.parser.parse_args([])
        self.assertEqual(args.workers, 0)
        self.assertEqual(args.cohort_mode, "auto")


@unittest.skipUnless(HAS_DEPS, "pydantic + PyYAML not installed")
class FormatEffectiveValuesTest(unittest.TestCase):
    def setUp(self):
        from syntheticgen.cli import _parser
        from syntheticgen.config import (
            load_and_validate_config,
            merge_config_into_args,
            parse_explicit_cli_args,
            format_effective_values,
        )
        self.parser = _parser(
            Path(__file__).resolve().parent.parent,
        )
        self.load = load_and_validate_config
        self.merge = merge_config_into_args
        self.explicit = parse_explicit_cli_args
        self.format = format_effective_values
        self.parser_defaults = {
            a.dest: a.default for a in self.parser._actions
            if a.dest not in ("help", "config", "no_config")
        }
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "c.yaml"
        self.path.write_text(
            "schema_version: 1\n"
            "cohort:\n  n: 3000\n"
            "performance:\n  cohort_mode: arrow\n  workers: 8\n"
        )

    def tearDown(self):
        self.tmp.cleanup()

    def test_tags_config_only_value(self):
        argv = []
        args = self.parser.parse_args(argv)
        cfg = self.load(self.path)
        explicit = self.explicit(self.parser, argv)
        merged = self.merge(args, cfg, explicit)
        lines = self.format(merged, self.parser_defaults, cfg, explicit)
        text = "\n".join(lines)
        self.assertIn("[config]", text)
        self.assertIn("n", text)
        self.assertIn("cohort_mode", text)

    def test_tags_cli_override(self):
        argv = ["--workers", "4"]
        args = self.parser.parse_args(argv)
        cfg = self.load(self.path)
        explicit = self.explicit(self.parser, argv)
        merged = self.merge(args, cfg, explicit)
        lines = self.format(merged, self.parser_defaults, cfg, explicit)
        text = "\n".join(lines)
        self.assertTrue(any(
            "workers" in line and "overrides config value 8" in line
            for line in lines
        ), text)

    def test_default_equal_values_suppressed(self):
        # A key whose value matches the default AND wasn't on the CLI
        # AND isn't in the config must NOT appear.
        argv = []
        args = self.parser.parse_args(argv)
        cfg = self.load(self.path)
        explicit = self.explicit(self.parser, argv)
        merged = self.merge(args, cfg, explicit)
        lines = self.format(merged, self.parser_defaults, cfg, explicit)
        text = "\n".join(lines)
        # ``build`` is not in the config and not on CLI; default GRCh38
        # should be silent.
        self.assertNotIn("build", text)


class SchemaSyncTest(unittest.TestCase):
    """The committed JSON Schema file at the repo root must match
    what generate_json_schema() produces from the current Pydantic
    models. If a developer changes a model but forgets to refresh the
    schema file, this test fails the PR.

    To refresh::

        .venv/bin/python -m syntheticgen.config \
            synthetic_people/generate_people_config.schema.json
    """

    @unittest.skipUnless(HAS_DEPS, "pydantic + PyYAML not installed")
    def test_committed_schema_matches_models(self):
        from syntheticgen.config import (
            generate_json_schema,
            serialize_schema,
        )
        if not REPO_SCHEMA_PATH.is_file():
            self.fail(
                f"Schema file missing at {REPO_SCHEMA_PATH}. Run: "
                "python -m syntheticgen.config "
                "synthetic_people/generate_people_config.schema.json"
            )
        on_disk = REPO_SCHEMA_PATH.read_text()
        expected = serialize_schema(generate_json_schema())
        if on_disk != expected:
            self.fail(
                f"{REPO_SCHEMA_PATH} is out of sync with the Pydantic "
                f"models in syntheticgen/config.py. Regenerate with: "
                f"python -m syntheticgen.config "
                f"synthetic_people/generate_people_config.schema.json"
            )


@unittest.skipUnless(HAS_DEPS, "pydantic + PyYAML not installed")
class RenderDefaultConfigYamlTest(unittest.TestCase):
    """The starter YAML emitted by ``--print-config`` must:

    - Be valid YAML the loader accepts.
    - Load as a no-op (every field equal to its built-in default,
      so a fresh user can edit incrementally).
    - Carry every top-level section the models define, so future
      additions don't silently fail to appear in the starter.
    - Carry a leading ``# description`` comment per field so the
      file documents itself.
    - Be deterministic so the output can be diffed across builds.
    """

    def setUp(self):
        from syntheticgen.config import render_default_config_yaml
        self.rendered = render_default_config_yaml()

    def test_round_trips_through_loader(self):
        # Write to a temp file and feed through the real loader the
        # cli would use. Round-trips imply every value rendered is
        # acceptable to the pydantic validators (including bounded
        # ranges, enum values, and cross-field rules).
        from syntheticgen.config import load_and_validate_config
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False,
        ) as fh:
            fh.write(self.rendered)
            path = Path(fh.name)
        try:
            cfg = load_and_validate_config(path)
        finally:
            path.unlink(missing_ok=True)
        self.assertEqual(cfg.schema_version, 1)

    def test_loaded_config_equals_minimal_config(self):
        # The starter is a no-op: every key matches the value the
        # loader would assign from defaults given ``schema_version:
        # 1`` alone. So both objects' dict-dumps must match.
        from syntheticgen.config import _models, load_and_validate_config
        Config = _models()
        minimal = Config(schema_version=1)
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False,
        ) as fh:
            fh.write(self.rendered)
            path = Path(fh.name)
        try:
            from_yaml = load_and_validate_config(path)
        finally:
            path.unlink(missing_ok=True)
        self.assertEqual(
            minimal.model_dump(), from_yaml.model_dump(),
            "Rendered starter must produce a no-op config",
        )

    def test_every_top_level_section_present(self):
        # Sentinel against silent omissions: every nested section
        # in the model must appear at the start of a line in the
        # rendered output.
        for section in (
            "cohort:", "simulation:", "overlays:",
            "structural_variants:", "sequencing_errors:",
            "performance:", "output:", "admixture:",
            "legacy_background:",
        ):
            self.assertIn(
                f"\n{section}\n", self.rendered,
                f"missing top-level section: {section}",
            )

    def test_overlays_subsections_present(self):
        # The nested overlays.{clinvar,rsid,cosmic} sections render
        # under ``overlays:`` with two-space indent.
        for sub in ("  clinvar:", "  rsid:", "  cosmic:"):
            self.assertIn(
                f"\n{sub}\n", self.rendered,
                f"missing overlays subsection: {sub}",
            )

    def test_carries_schema_version_line(self):
        self.assertIn("schema_version: 1\n", self.rendered)

    def test_carries_yaml_language_server_pragma(self):
        # The IDE-integration comment must be near the top so VS
        # Code / IntelliJ pick up the schema for autocomplete.
        self.assertIn(
            "# yaml-language-server: $schema=", self.rendered,
        )

    def test_every_leaf_field_has_a_description_comment(self):
        # For every ``key: value`` line that isn't a section header
        # (``key:`` alone) and isn't ``schema_version`` (which has
        # its own multi-line preamble), the *immediately preceding*
        # non-blank line must be a ``# ...`` comment. This is the
        # property that makes the file self-documenting and is
        # easy to regress when new fields are added.
        import re
        lines = self.rendered.split("\n")
        leaf_pat = re.compile(r"^(\s*)([A-Za-z_]\w*):\s+\S")
        for idx, line in enumerate(lines):
            m = leaf_pat.match(line)
            if not m:
                continue
            key = m.group(2)
            if key == "schema_version":
                continue
            # Walk backwards over blank lines to the previous line.
            j = idx - 1
            while j >= 0 and lines[j].strip() == "":
                j -= 1
            self.assertTrue(
                j >= 0 and lines[j].lstrip().startswith("#"),
                f"leaf field {key!r} on line {idx + 1} has no "
                f"preceding # description comment "
                f"(previous non-blank line: {lines[j]!r})",
            )

    def test_deterministic_across_calls(self):
        # Two consecutive renderings must be byte-identical so the
        # output is diff-stable for users who commit it.
        from syntheticgen.config import render_default_config_yaml
        self.assertEqual(self.rendered, render_default_config_yaml())

    def test_includes_carrier_field_descriptions(self):
        # Spot-check a handful of fields by their description text
        # so we catch the case where the renderer accidentally
        # stops emitting the ``description=`` from pydantic.
        for snippet in (
            "Cohort size (number of person VCFs).",
            "Reference build assembly.",
            "Cohort intermediate between simulation and BCF write.",
            "Output directory for per-person VCFs and cohort BCFs.",
            "Run M6 EUR + SAS + AFR -> UK admixture and emit ancestry BEDs.",
        ):
            self.assertIn(snippet, self.rendered)


@unittest.skipUnless(HAS_DEPS, "pydantic + PyYAML not installed")
class PrintConfigCliTest(unittest.TestCase):
    """The ``--print-config`` flag short-circuits cli.main and writes
    the starter YAML to stdout. This test invokes the real main
    function with stdout captured so the integration is exercised
    end-to-end."""

    def test_print_config_writes_to_stdout_and_exits_zero(self):
        from syntheticgen.cli import main
        from syntheticgen.config import render_default_config_yaml
        buf = io.StringIO()
        real_stdout = sys.stdout
        sys.stdout = buf
        try:
            rc = main(["--print-config"])
        finally:
            sys.stdout = real_stdout
        self.assertEqual(rc, 0)
        self.assertEqual(buf.getvalue(), render_default_config_yaml())

    def test_print_config_output_is_valid_loader_input(self):
        # End-to-end: capture the cli's stdout, feed it back into
        # the loader, expect a no-op config (every field == default).
        from syntheticgen.cli import main
        from syntheticgen.config import (
            _models, load_and_validate_config,
        )
        buf = io.StringIO()
        real_stdout = sys.stdout
        sys.stdout = buf
        try:
            main(["--print-config"])
        finally:
            sys.stdout = real_stdout
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False,
        ) as fh:
            fh.write(buf.getvalue())
            path = Path(fh.name)
        try:
            cfg = load_and_validate_config(path)
        finally:
            path.unlink(missing_ok=True)
        Config = _models()
        self.assertEqual(
            cfg.model_dump(), Config(schema_version=1).model_dump(),
        )


if __name__ == "__main__":
    unittest.main()
