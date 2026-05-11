"""Config-file support for generate_people.

Optional alongside the CLI: the user can drop a YAML file named
``generate_people_config.yaml`` in the working directory and the
tool picks it up automatically. CLI flags still win over config
values; config values still win over built-in defaults. Discovery
is cwd-only by design — predictability for scientists running
multiple jobs from different directories beats convenience.

This module owns:

- The Pydantic models that define the config schema (every field has
  a default matching the existing argparse defaults, so a config
  file containing only ``schema_version: 1`` is valid and changes
  nothing).
- The discovery + load helpers ``discover_config_file`` and
  ``load_and_validate_config``.
- The merge helper ``merge_config_into_args`` that resolves
  ``cli > config > defaults`` precedence, returning the dest-name
  argparse Namespace the rest of cli.py expects.
- The "effective values" formatter that explains, line by line,
  which non-default value came from which source — so a scientist
  surprised by the run's behaviour can read the cause from stderr
  rather than guess at it.
- A JSON Schema export for IDE integration (VS Code / IntelliJ YAML
  language server). The committed schema file is regenerated from
  the Pydantic models; a sync test in ``tests/test_config.py``
  fails CI if the committed file drifts.

Conditional deps: ``pydantic>=2.0`` + ``PyYAML>=6.0``. Both are now
in ``requirements.txt`` — config is a core UX feature, not a
perf-scale opt-in like ``pyarrow``.
"""

from __future__ import annotations

import argparse
import json
import sys
from argparse import Namespace
from pathlib import Path
from typing import Any, Literal, Optional


CURRENT_SCHEMA_VERSION = 1
DEFAULT_CONFIG_FILENAME = "generate_people_config.yaml"
SCHEMA_FILENAME = "generate_people_config.schema.json"


def _require_deps():
    try:
        import pydantic  # noqa: F401
        import yaml  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "Config-file support requires pydantic + PyYAML. Reinstall "
            "with `pip install -r synthetic_people/requirements.txt`."
        ) from exc


# ---------------------------------------------------------------------------
# Pydantic models — one per config-file section
# ---------------------------------------------------------------------------
#
# Every field has a default matching the existing argparse defaults so an
# under-specified config is valid (and a no-op). The ``argparse_dest`` in
# each field's ``json_schema_extra`` is the dest-name the merge helper
# uses to map config keys back onto argparse args; it's also visible to
# downstream tooling that introspects the schema.


def _models():
    """Lazy model construction — defers the pydantic import so this
    module is importable without pydantic for syntax-level tools."""
    _require_deps()
    from pydantic import BaseModel, ConfigDict, Field, model_validator

    class _Strict(BaseModel):
        model_config = ConfigDict(
            extra="forbid",
            validate_default=True,
            populate_by_name=True,
        )

    class CohortConfig(_Strict):
        n: int = Field(
            default=10, ge=1,
            description="Cohort size (number of person VCFs).",
            json_schema_extra={"argparse_dest": "n"},
        )
        build: Literal["GRCh37", "GRCh38"] = Field(
            default="GRCh38",
            description="Reference build assembly.",
            json_schema_extra={"argparse_dest": "build"},
        )
        seed: Optional[int] = Field(
            default=None,
            description="Master RNG seed; omit for fresh randomness each run.",
            json_schema_extra={"argparse_dest": "seed"},
        )
        chromosomes: str = Field(
            default="22",
            description="Chromosomes spec: list / range / mix (e.g. '22', '19-22,X').",
            json_schema_extra={"argparse_dest": "chromosomes"},
        )
        chr_length_mb: float = Field(
            default=5.0, ge=0.0,
            description="Simulated prefix per chromosome in Mb; 0 = full length.",
            json_schema_extra={"argparse_dest": "chr_length_mb"},
        )

    class SimulationConfig(_Strict):
        demo_model: str = Field(
            default="OutOfAfrica_3G09",
            description="stdpopsim model id, or 'none' for constant-size Ne.",
            json_schema_extra={"argparse_dest": "demo_model"},
        )
        population: str = Field(
            default="CEU",
            description="Sampling population for the demographic model.",
            json_schema_extra={"argparse_dest": "population"},
        )
        rec_rate: float = Field(
            default=1e-8, ge=0.0,
            description="Recombination rate, consulted only when demo_model='none'.",
            json_schema_extra={"argparse_dest": "rec_rate"},
        )
        mu: float = Field(
            default=1.29e-8, ge=0.0,
            description="Mutation rate, consulted only when demo_model='none'.",
            json_schema_extra={"argparse_dest": "mu"},
        )

    class ClinVarConfig(_Strict):
        inject_density: float = Field(
            default=0.01, ge=0.0, le=1.0,
            description="Fraction of cohort sites overwritten with ClinVar records.",
            json_schema_extra={"argparse_dest": "clinvar_inject_density"},
        )
        sig_filter: str = Field(
            default="Pathogenic,Likely_pathogenic,Pathogenic/Likely_pathogenic",
            description="Comma-separated CLNSIG values to keep from the ClinVar VCF.",
            json_schema_extra={"argparse_dest": "clinvar_sig"},
        )

    class RsidConfig(_Strict):
        density: float = Field(
            default=0.20, ge=0.0, le=1.0,
            description="Fraction of cohort sites rewritten to a known dbSNP variant.",
            json_schema_extra={"argparse_dest": "rsid_density"},
        )
        vcf: Optional[str] = Field(
            default=None,
            description="Override the rsID source VCF; default uses the cached ClinVar VCF's INFO/RS.",
            json_schema_extra={"argparse_dest": "dbsnp_vcf"},
        )

    class CosmicConfig(_Strict):
        enabled: bool = Field(
            default=False,
            description="Enable COSMIC overlay (--somatic on CLI).",
            json_schema_extra={"argparse_dest": "somatic"},
        )
        vcf: Optional[str] = Field(
            default=None,
            description="Path to COSMIC-format VCF (registration required).",
            json_schema_extra={"argparse_dest": "cosmic_vcf"},
        )
        inject_density: float = Field(
            default=0.005, ge=0.0, le=1.0,
            description="Fraction of cohort sites overwritten with COSMIC records when enabled.",
            json_schema_extra={"argparse_dest": "cosmic_inject_density"},
        )

    class OverlaysConfig(_Strict):
        clinvar: ClinVarConfig = ClinVarConfig()
        rsid: RsidConfig = RsidConfig()
        cosmic: CosmicConfig = CosmicConfig()

    class StructuralVariantsConfig(_Strict):
        per_person: int = Field(
            default=3, ge=0,
            description="Structural variants (DEL/DUP/INV) per person.",
            json_schema_extra={"argparse_dest": "svs_per_person"},
        )
        length_min: int = Field(
            default=50, ge=1,
            description="Minimum SV length in bp (log-uniform draw).",
            json_schema_extra={"argparse_dest": "sv_length_min"},
        )
        length_max: int = Field(
            default=10000, ge=1,
            description="Maximum SV length in bp.",
            json_schema_extra={"argparse_dest": "sv_length_max"},
        )

        @model_validator(mode="after")
        def _length_bounds(self):
            if self.length_max < self.length_min:
                raise ValueError(
                    f"structural_variants.length_max ({self.length_max}) "
                    f"must be >= length_min ({self.length_min})"
                )
            return self

    class SequencingErrorsConfig(_Strict):
        gt_flip_rate: float = Field(
            default=0.001, ge=0.0, le=1.0,
            description="Per-call probability of a GT flip (lightweight noise model).",
            json_schema_extra={"argparse_dest": "error_rate"},
        )
        dropout_rate: float = Field(
            default=0.0005, ge=0.0, le=1.0,
            description="Per-call probability of a coverage dropout.",
            json_schema_extra={"argparse_dest": "dropout_rate"},
        )

    class PerformanceConfig(_Strict):
        workers: int = Field(
            default=0, ge=0,
            description="Worker processes; 0=auto (cpu_count), 1=serial.",
            json_schema_extra={"argparse_dest": "workers"},
        )
        cohort_mode: Literal["auto", "sites_list", "arrow"] = Field(
            default="auto",
            description="Cohort intermediate between simulation and BCF write.",
            json_schema_extra={"argparse_dest": "cohort_mode"},
        )
        cohort_arrow_batch_size: int = Field(
            default=256, ge=1,
            description="Sites per Arrow record batch when cohort_mode='arrow'.",
            json_schema_extra={"argparse_dest": "cohort_arrow_batch_size"},
        )
        fanout_batch_size: int = Field(
            default=4, ge=1,
            description="Persons grouped per bcftools-query invocation during fan-out.",
            json_schema_extra={"argparse_dest": "fanout_batch_size"},
        )
        chr_chunk_mb: float = Field(
            default=0.0, ge=0.0,
            description="Per-chromosome msprime chunk size in Mb; 0=auto-pick.",
            json_schema_extra={"argparse_dest": "chr_chunk_mb"},
        )
        no_resume: bool = Field(
            default=False,
            description="If true, ignore any existing cohort.meta.json + cohort BCFs.",
            json_schema_extra={"argparse_dest": "no_resume"},
        )
        profile_memory: Optional[str] = Field(
            default=None,
            description="Path to write a memory-profile TSV (~1s cadence + phase marks).",
            json_schema_extra={"argparse_dest": "profile_memory"},
        )

    class OutputConfig(_Strict):
        dir: str = Field(
            default="./out",
            description="Output directory for per-person VCFs and cohort BCFs.",
            json_schema_extra={"argparse_dest": "output_dir"},
        )
        cache_dir: str = Field(
            default="./cache",
            description="Cache directory for ClinVar / dbSNP / COSMIC downloads.",
            json_schema_extra={"argparse_dest": "cache_dir"},
        )
        mode: Literal["per-person", "cohort", "both"] = Field(
            default="per-person",
            description="What the run writes: per-person VCFs, cohort BCFs only, or both.",
            json_schema_extra={"argparse_dest": "mode"},
        )

    class AdmixtureConfig(_Strict):
        enabled: bool = Field(
            default=False,
            description="Run M6 EUR + SAS + AFR -> UK admixture and emit ancestry BEDs.",
            json_schema_extra={"argparse_dest": "admixture"},
        )
        eur_frac: float = Field(
            default=0.60, ge=0.0, le=1.0,
            description="European ancestry proportion.",
            json_schema_extra={"argparse_dest": "eur_frac"},
        )
        sas_frac: float = Field(
            default=0.25, ge=0.0, le=1.0,
            description="South Asian ancestry proportion.",
            json_schema_extra={"argparse_dest": "sas_frac"},
        )
        afr_frac: float = Field(
            default=0.15, ge=0.0, le=1.0,
            description="African ancestry proportion (eur + sas + afr must sum to 1.0).",
            json_schema_extra={"argparse_dest": "afr_frac"},
        )

        @model_validator(mode="after")
        def _fractions_sum(self):
            if self.enabled:
                total = self.eur_frac + self.sas_frac + self.afr_frac
                if abs(total - 1.0) > 1e-6:
                    raise ValueError(
                        f"admixture.eur_frac + sas_frac + afr_frac must sum to "
                        f"1.0; got {total} (eur={self.eur_frac}, "
                        f"sas={self.sas_frac}, afr={self.afr_frac})"
                    )
            return self

    class LegacyBackgroundConfig(_Strict):
        enabled: bool = Field(
            default=False,
            description="Use M4 1000G-pool + power-law SFS sampler (--legacy-background).",
            json_schema_extra={"argparse_dest": "legacy_background"},
        )
        background_glob: Optional[str] = Field(
            default=None,
            description="Source glob(s) for common variants (legacy path only).",
            json_schema_extra={"argparse_dest": "background_glob"},
        )
        n_background: int = Field(
            default=500, ge=1,
            description="Shared background site count (legacy path only).",
            json_schema_extra={"argparse_dest": "n_background"},
        )
        af_min: float = Field(
            default=0.05, ge=0.0, le=1.0,
            description="Minimum AF when loading the legacy pool.",
            json_schema_extra={"argparse_dest": "af_min"},
        )
        sfs_alpha: float = Field(
            default=2.0, gt=0.0,
            description="Power-law exponent for the legacy SFS.",
            json_schema_extra={"argparse_dest": "sfs_alpha"},
        )

    class Config(_Strict):
        """Top-level config-file schema."""
        schema_version: int = Field(
            ...,
            description=(
                "Config-file schema version. Required so the loader "
                "can reject incompatibly-old configs with a clear "
                "message rather than silently mis-interpreting fields."
            ),
        )
        cohort: CohortConfig = CohortConfig()
        simulation: SimulationConfig = SimulationConfig()
        overlays: OverlaysConfig = OverlaysConfig()
        structural_variants: StructuralVariantsConfig = StructuralVariantsConfig()
        sequencing_errors: SequencingErrorsConfig = SequencingErrorsConfig()
        performance: PerformanceConfig = PerformanceConfig()
        output: OutputConfig = OutputConfig()
        admixture: AdmixtureConfig = AdmixtureConfig()
        legacy_background: LegacyBackgroundConfig = LegacyBackgroundConfig()

        @model_validator(mode="after")
        def _schema_version_supported(self):
            if self.schema_version != CURRENT_SCHEMA_VERSION:
                raise ValueError(
                    f"schema_version {self.schema_version} is not supported "
                    f"by this build (expected {CURRENT_SCHEMA_VERSION}). "
                    f"Update the config to schema_version: "
                    f"{CURRENT_SCHEMA_VERSION} and consult the changelog "
                    f"for any field renames."
                )
            return self

    return Config


def discover_config_file(cwd: Path) -> Optional[Path]:
    """Return ``cwd / generate_people_config.yaml`` if it exists, else
    ``None``. Discovery is cwd-only by design — see module docstring."""
    candidate = Path(cwd) / DEFAULT_CONFIG_FILENAME
    return candidate if candidate.is_file() else None


def load_and_validate_config(path: Path):
    """Open a YAML config, run Pydantic validation, return the parsed
    Config object. Validation errors are repackaged into a single
    ``SystemExit`` with every problem listed at once (citing the field
    path) — scientists shouldn't have to fix one error, re-run, see
    the next, fix it, re-run, etc.
    """
    _require_deps()
    import yaml
    from pydantic import ValidationError

    Config = _models()
    path = Path(path)
    try:
        with open(path, "r", encoding="utf-8") as fh:
            data = yaml.safe_load(fh)
    except yaml.YAMLError as exc:
        raise SystemExit(
            f"Config file at {path} is not valid YAML: {exc}"
        )
    if data is None:
        # Empty file. Treat as "config present but empty" — the
        # schema_version requirement still bites, which is the right
        # message to surface.
        data = {}
    if not isinstance(data, dict):
        raise SystemExit(
            f"Config file at {path} must be a YAML mapping at top level; "
            f"got {type(data).__name__}."
        )
    try:
        return Config.model_validate(data)
    except ValidationError as exc:
        lines = [f"Config validation failed in {path}:"]
        for err in exc.errors():
            loc = ".".join(str(p) for p in err["loc"])
            lines.append(f"  - {loc}: {err['msg']}")
        raise SystemExit("\n".join(lines))


def _flatten_config_to_argparse_dests(config) -> dict:
    """Walk the config model and produce a flat
    ``{argparse_dest: value}`` mapping. Only fields the user
    explicitly set survive — Pydantic-fill defaults are dropped.

    Subtle: ``model_fields_set`` only reports keys at one level; we
    recurse into nested BaseModels (the overlay / admixture / etc.
    sections) so a config that sets ``overlays.clinvar.inject_density``
    surfaces just that field as explicit, not its sibling defaults.
    """
    out: dict = {}
    if config is None:
        return out

    from pydantic import BaseModel
    for section_name in config.model_fields_set:
        if section_name == "schema_version":
            continue
        section = getattr(config, section_name)
        if isinstance(section, BaseModel):
            _flatten_section(section, out)
        else:
            # Top-level scalar (none today, but future-proof).
            _add_with_dest(config, section_name, section, out)
    return out


def _flatten_section(section, out: dict) -> None:
    from pydantic import BaseModel
    for field_name in section.model_fields_set:
        value = getattr(section, field_name)
        if isinstance(value, BaseModel):
            _flatten_section(value, out)
        else:
            _add_with_dest(section, field_name, value, out)


def _add_with_dest(model, field_name: str, value, out: dict) -> None:
    field = type(model).model_fields[field_name]
    extra = field.json_schema_extra or {}
    dest = extra.get("argparse_dest") if isinstance(extra, dict) else None
    if dest is None:
        return
    out[dest] = value


def parse_explicit_cli_args(parser, argv) -> set:
    """Return the set of argparse ``dest`` names the user explicitly
    passed in ``argv``. The shadow parse uses
    ``default=argparse.SUPPRESS`` on every action AND clears
    ``parser._defaults`` (the dict ``parser.set_defaults(...)``
    populates outside of any action) so the resulting Namespace
    contains only user-typed flags. The original parser is not
    mutated permanently — both maps are saved and restored.
    """
    saved_action_defaults = {a.dest: a.default for a in parser._actions}
    saved_parser_defaults = dict(parser._defaults)
    try:
        for action in parser._actions:
            action.default = argparse.SUPPRESS
        parser._defaults.clear()
        explicit_ns = parser.parse_args(argv)
    finally:
        for action in parser._actions:
            action.default = saved_action_defaults.get(
                action.dest, action.default,
            )
        parser._defaults.clear()
        parser._defaults.update(saved_parser_defaults)
    return set(vars(explicit_ns).keys())


def merge_config_into_args(
    args: Namespace,
    config,
    explicit_cli: set,
) -> Namespace:
    """Apply ``cli > config > defaults`` precedence.

    ``args`` is the regular argparse Namespace (defaults filled in).
    ``config`` is the loaded Pydantic Config or ``None``.
    ``explicit_cli`` is the set of dest names the user actually typed.

    Returns a fresh Namespace with the merged values. The original
    ``args`` is left untouched so callers can compare before/after if
    they want.
    """
    merged = Namespace(**vars(args))
    if config is None:
        return merged

    flat = _flatten_config_to_argparse_dests(config)
    for dest, value in flat.items():
        if dest in explicit_cli:
            continue  # CLI wins
        if not hasattr(merged, dest):
            continue  # config has a field with no CLI counterpart; skip
        setattr(merged, dest, value)
    return merged


def format_effective_values(
    merged: Namespace,
    parser_defaults: dict,
    config,
    explicit_cli: set,
) -> list:
    """Produce per-line stderr output documenting where each non-
    default value came from.

    Output looks like::

        cohort.n            = 3000           [config]
        performance.workers = 8              [cli, overrides config value 4]

    Only keys whose resolved value differs from the argparse default
    are emitted — the goal is to surface the run's *effective*
    behaviour, not list everything verbosely.
    """
    if config is None:
        config_flat: dict = {}
    else:
        config_flat = _flatten_config_to_argparse_dests(config)

    lines = []
    for dest, default in parser_defaults.items():
        value = getattr(merged, dest, default)
        if value == default and dest not in explicit_cli and dest not in config_flat:
            continue

        if dest in explicit_cli:
            if dest in config_flat and config_flat[dest] != value:
                source = f"cli, overrides config value {config_flat[dest]!r}"
            else:
                source = "cli"
        elif dest in config_flat:
            source = "config"
        else:
            # Resolved value differs from default but neither CLI nor config
            # claimed it — should not happen in practice, but be defensive.
            source = "default"

        lines.append(f"  {dest:<28} = {value!r:<20} [{source}]")
    return lines


# ---------------------------------------------------------------------------
# JSON Schema export (for VS Code / IntelliJ YAML language server)
# ---------------------------------------------------------------------------


def render_default_config_yaml() -> str:
    """Render a starter ``generate_people_config.yaml`` with every
    field set to its built-in default and a leading ``# description``
    comment for each key.

    The output is a fully valid config — saving it as
    ``generate_people_config.yaml`` and running with no flags
    behaves identically to running with no config at all. Intended
    use is ``generate_people --print-config > generate_people_config.yaml``
    as a starting point a new user can then edit down.

    Determinism: field iteration order is fixed by the pydantic
    model definitions, and ``yaml.safe_dump`` is used per-leaf so
    booleans, ``None`` (rendered ``null``), and strings escape the
    way the loader expects. Re-emission is byte-identical for a
    given build of the models.
    """
    _require_deps()
    Config = _models()

    lines: list[str] = [
        "# generate_people_config.yaml",
        "#",
        "# Starter config emitted by `generate_people --print-config`.",
        "# Every field is set to its built-in default — running with",
        "# this file and no CLI flags behaves identically to running",
        "# with no config at all. Edit the values you want to change;",
        "# CLI flags still override config values, config values",
        "# still override built-in defaults.",
        "#",
        "# Schema reference (for IDE auto-complete / validation):",
        "# https://github.com/nigel4321/synthetic-people/blob/main/"
        "synthetic_people/generate_people_config.schema.json",
        "# Full documentation: TUTORIAL.md §10.",
        "",
        "# yaml-language-server: $schema=./"
        f"{SCHEMA_FILENAME}",
        "",
        "# Schema version. Required. The loader rejects unknown",
        "# schema versions with a clear message so an incompatible",
        f"# field rename can never silently mis-interpret an old config.",
        f"schema_version: {CURRENT_SCHEMA_VERSION}",
        "",
    ]
    _render_model_fields(
        Config, lines, indent=0, skip=("schema_version",),
    )
    # Strip any trailing blank lines, then ensure exactly one
    # terminal newline.
    while lines and lines[-1] == "":
        lines.pop()
    return "\n".join(lines) + "\n"


def _render_model_fields(
    model_class,
    lines: list,
    indent: int,
    skip: tuple = (),
) -> None:
    """Walk a pydantic model's fields and append YAML lines.

    Sub-models (nested ``BaseModel`` subclasses) become indented
    sections. Leaf fields become ``key: value`` lines preceded by
    a ``# description`` comment line drawn from the field's
    pydantic ``description=``.
    """
    import yaml
    from pydantic import BaseModel

    pad = "  " * indent
    fields = list(model_class.model_fields.items())
    for idx, (name, field) in enumerate(fields):
        if name in skip:
            continue
        annotation = field.annotation
        description = field.description or ""

        if isinstance(annotation, type) and issubclass(
            annotation, BaseModel
        ):
            # Nested section: descend with one extra indent level.
            # A blank line between top-level sections makes the
            # output readable when redirected to a file.
            if indent == 0 and idx > 0:
                lines.append("")
            if description:
                lines.append(f"{pad}# {description}")
            lines.append(f"{pad}{name}:")
            _render_model_fields(annotation, lines, indent + 1)
            continue

        if description:
            lines.append(f"{pad}# {description}")
        default_value = field.default
        # safe_dump produces compact one-line values (``42``,
        # ``true``, ``null``, ``'1-22'``) but appends a ``...`` YAML
        # end-of-document marker on scalar values. Take just the
        # first line so the marker is dropped.
        rendered = yaml.safe_dump(
            default_value, default_flow_style=True,
        ).split("\n", 1)[0]
        lines.append(f"{pad}{name}: {rendered}")


def generate_json_schema() -> dict:
    """Return the JSON Schema for the config file, suitable for
    publishing as ``generate_people_config.schema.json``."""
    Config = _models()
    schema = Config.model_json_schema()
    schema["$schema"] = "https://json-schema.org/draft/2020-12/schema"
    schema["title"] = "generate_people_config"
    schema["description"] = (
        "Schema for synthetic_people's generate_people_config.yaml. "
        "Documented at: https://github.com/nigel4321/synthetic-people/"
        "blob/main/synthetic_people/TUTORIAL.md"
    )
    return schema


def serialize_schema(schema: dict) -> str:
    """Deterministic serialisation of the JSON Schema (sorted keys,
    2-space indent, trailing newline) so the committed file is
    diff-stable across regenerations."""
    return json.dumps(schema, indent=2, sort_keys=True) + "\n"


def regenerate_schema_file(schema_path: Path) -> None:
    """Write a freshly-generated schema to disk. Used by the
    developer-facing regenerate script and the sync test."""
    schema_path = Path(schema_path)
    schema_path.write_text(serialize_schema(generate_json_schema()))


if __name__ == "__main__":
    # Allow `python -m syntheticgen.config <out_path>` to dump the
    # JSON Schema to disk. Convenient when a developer changes a
    # model and needs to refresh the committed file.
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(SCHEMA_FILENAME)
    regenerate_schema_file(out)
    print(f"Wrote {out}")
