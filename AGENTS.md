# genmod: agent guide

genmod is a click CLI that annotates and ranks variants in VCF files: genetic inheritance models
per family, region/frequency annotation, rank scores from .ini models, compound scoring.
Python >= 3.11 (CI runs 3.13 only), managed with uv. Entry point `genmod.commands.base:cli`.

## Commands and how they chain

Typical pipeline: `annotate` → `models` → `score` → `compound` → `sort`.

| Command | Module | Adds / does |
|---|---|---|
| `annotate` | commands/annotate_variant.py | `Annotation` (genes from bundled genmod/annotations/ensembl_genes_{37,38}.txt.gz, region start padded 4 kb), tabix scores/frequencies (1000G, ExAC, CADD, SPIDEX) via pysam |
| `models` | commands/annotate_models.py | `GeneticModels`, `Compounds`, `ModelScore` per family (PED file via ped_parser). Genes from INFO `Annotation` (`-k`), or `--vep` (CSQ) / `--snpeff` (ANN) |
| `score` | commands/score_variants.py | `RankScore`, `RankScoreNormalized`, `RankScoreMinMax`, `RankResult` from a rank model .ini (`-c`). Single process |
| `compound` | commands/score_compounds.py | Rescores compound pairs; needs `RankScoreMinMax` from `score` and `-o` |
| `sort` | commands/genmod_sort.py | Sorts by rank score or position (shells out to unix `sort -V`) |
| `filter` | commands/filter_variants.py | Threshold on one INFO float (default `1000GAF` < 0.05) |

`commands/analyze.py` and `commands/summarize_variants.py` are dead (not registered).

## Core mechanics

- **VCF handling** (genmod/vcf_tools): `HeaderParser` parses `##` lines (`vep_columns` from CSQ
  `Format:`, `ann_columns` from the SnpEff ANN description). Variants are plain dicts:
  `get_variant_dict` + `get_info_dict` (raw strings; flags map to `[]`).
  `parse_annotations(raw, columns)` turns CSQ/ANN into a list of dicts.
- **models / compound are multiprocess**: `utils/get_batches.py` (main process) groups consecutive
  variants whose genes overlap into a batch → `JoinableQueue` → workers (`VariantAnnotator`,
  `CompoundScorer`) → results queue → `VariantPrinter` process. One `None` per worker stops it.
  Workers only read `variant["annotation"]` (a set of genes); CSQ/ANN is resolved in `get_batches`
  via `ConsequenceSource` (`VEP`/`SNPEFF` in utils/get_features.py) and `INTERESTING_SO_TERMS`.
- **Genetic models** (annotate_models/genetic_models.py, models/*.py): AD, AD_dn, AR_hom, AR_hom_dn,
  AR_comp, AR_comp_dn, XR, XR_dn, XD, XD_dn. Only X gets X-linked models; Y/MT are treated as
  autosomal. `--strict`, `--phased` (PS tags), `-r` reduced-penetrance genes (IDs must match the
  annotation's gene IDs) change the rules.
- **Rank models** (.ini via configobj; bundled in genmod/configs/): `[Version]`, `[Categories]`
  (`category_aggregation` min/max/sum), then one section per plugin with sub-sections as rules.
  `config_parser.py` builds a `Plugin` (genmod/plugin/plugin.py, vendored from extract_vcf) and a
  `ScoreFunction` per plugin.
  - Numeric rules are intervals `[lower, upper)`; string rules need `string` + `priority`
    (priority picks the winner, score scores it); `[[not_reported]]` sets the missing score.
  - String rules match whole tokens after splitting on `separators` (`exact_match = false` = legacy
    substring matching).
  - `info_key = CSQ` or `ANN` + `csq_key = <column>`: genmod itself splits transcripts on `,` and
    columns on `|`; `separators` only split inside the column (use `'&'`, never `,`/`|`).

## Development

```bash
uv sync
uv run pytest                          # run from repo root: fixture paths are relative
uv run ruff check .                    # lint (CI: lint.yml)
uv run ruff check --select I .         # import order, NOT covered by `ruff check .` (CI: check_format.yml)
uv run ruff format --check .           # CI runs `ruff format --diff .`
uv run ruff check --select I --fix . && uv run ruff format .   # fix both
uv run genmod --help
```

- CI on every PR: tests in Docker (`make test`, `make test_dist`), lint, import order + format, and a
  changelog check: **every PR must edit CHANGELOG.md** (or carry the `Skip-Changelog` label).
- Version lives only in `pyproject.toml` (+ `uv lock`); `__version__` is read from package metadata.
  Release steps: RELEASE.md. Stale, ignore: `.travis.yml`, `setup.cfg`, `MANIFEST.in`, `CITATION.cff`.
- No type checker is configured; new code gets type hints anyway (builtin generics, `Optional`).

## Conventions

- CHANGELOG.md: Keep a Changelog, entries under `## [unreleased]` → Added/Changed/Fixed, one line
  each, ending with the **PR** link `([#N](https://github.com/Clinical-Genomics/genmod/pull/N))`.
  The issue number goes in the PR description/commit, not the changelog. Flag scoring behaviour
  changes in bold.
- Branches: `issue-<N>-<slug>`, `fix-<slug>`, `docs-<slug>`. Commits: short, imperative, capitalised.
- Logging: `logger = logging.getLogger(__name__)`, lazy `%s` args, not f-strings/.format.
- Tests: `tests/<area>/test_*.py`; new tests use `# GIVEN / # WHEN / # THEN` comments and type hints.
  CLI tests use `CliRunner().invoke(cmd, args)` and parse `result.stdout_bytes`; helpers in
  `tests/functionality/test_utils.py` (imported as `from test_utils import ...`, which works only
  because that folder has no `__init__.py`; don't add one). Fixtures: `tests/fixtures/` (SnpEff/VEP
  twins in `tests/fixtures/snpeff/`, ClinVar CSQ in `score_variants/`).

## Gotchas

- `models -p 1` prints from the `VariantPrinter` child process, so `CliRunner` captures nothing:
  don't pass `-p 1` in CLI tests.
- Compound pair order in `Compounds=` comes from a set: compare as sets in tests.
- `compound` writes via `outfile.name`: it needs `-o`. Its batching must use the same gene source
  (`Annotation` / `--vep` / `--snpeff`) as `models` did.
- `str.lstrip("chr")` / `.strip("RankScore")` are used as prefix strips in places; they strip
  characters. Don't copy that pattern.
- A plugin value of `0`/`0.0` is scored as not reported (`if not value`), and scores are truncated
  with `int()`.
