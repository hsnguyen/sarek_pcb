# Copilot Instructions for nf-core/sarek

## Project and architecture

`nf-core/sarek` is a Nextflow DSL2 pipeline for germline and somatic variant
analysis of whole-genome, exome, and targeted sequencing data.

- `main.nf` is the entry point: it resolves genome attributes, runs pipeline
  initialization and completion utilities, and invokes `NFCORE_SAREK`.
- `NFCORE_SAREK` prepares references, intervals, optional annotation caches,
  and SnpSift databases before passing channels to `workflows/sarek/main.nf`.
- `workflows/sarek/main.nf` dispatches input by type, performs preprocessing,
  separates normal, tumor-only, and matched tumor/normal paths, invokes the
  variant-calling dispatcher subworkflows, then post-processes, annotates, and
  aggregates reports and software versions for MultiQC.
- Atomic processes live in `modules/`; imported nf-core components are in
  `modules/nf-core/`, while pipeline-owned components are in `modules/local/`.
  Compose modules in `subworkflows/local/`; do not put process logic directly
  in the top-level workflows unless it is pipeline orchestration.
- `nextflow.config` owns parameter defaults, profiles, plugins, and includes
  `conf/base.config` plus the tool-specific files under `conf/modules/`.
  Put process arguments, prefixes, `publishDir` rules, and resource overrides
  in the relevant module config. Add user-facing parameter defaults to
  `nextflow.config` and update `nextflow_schema.json` with
  `nf-core pipelines schema build`.
- `conf/test.config` provides the minimal CI dataset. Tests use `tests/` and
  shared scenario/snapshot helpers in `tests/lib/UTILS.groovy`; imported
  nf-core module tests are deliberately ignored by `nf-test.config`.

## Environment, testing, and linting

Run nf-core commands from the project environment:

```bash
conda activate nf-core
```

Use Docker, Singularity, or Conda as the container profile; profiles are
comma-separated. The smallest useful pipeline smoke test is:

```bash
nextflow run . -profile test,docker --outdir results
nextflow run . -profile debug,test,docker --outdir results
```

Run nf-test tests and snapshots as follows:

```bash
# Full suite
nf-test test --profile debug,test,docker --verbose

# One test file
nf-test test tests/variant_calling_haplotypecaller.nf.test --profile debug,test,docker

# Fast plumbing check
nf-test test tests/default.nf.test --profile debug,test,docker -stub

# Update an intentional snapshot change
nf-test test tests/<name>.nf.test --profile debug,test,docker --update-snapshot
```

Run formatting/hooks and pipeline lint before a PR:

```bash
pre-commit run --all-files
nf-core pipelines lint
```

CI shards changed nf-test cases across Docker, Conda, and Singularity, and
runs against Nextflow `25.10.2` plus `latest-everything`. Use the `debug`
profile to surface process-selector warnings locally.

## Nextflow conventions

- Read `docs/DEVELOPER_GUIDELINES.md` before changing pipeline code. Use
  4-space indentation and Harshil alignment for `include`, `take:`, and
  `emit:` blocks.
- Any code touched in a PR must follow strict syntax: give closures explicit
  parameters, use underscore-prefixed names for deliberately dropped tuple
  values, and prefer explicit types where applicable.
- Keep sample metadata in the leading `meta` map of tuples. Add fields with
  `meta + [key: value]` and remove temporary fields with
  `meta - meta.subMap(...)` before emitting. Key fields include `patient`,
  `sample`, `status`, `lane`, `id`, `data_type`, `num_intervals`, and
  `variantcaller`.
- Name channels `ch_output_from_<process>` for initial outputs and
  `ch_<previousprocess>_for_<nextprocess>` for intermediate or terminal
  channels.
- Prefer dataflow controls (`filter`, `branch`, `mix`) to `if` statements and
  never add `ext.when`. Existing `ext.when` configuration is legacy and should
  be replaced with channel-driven flow when touched.
- A keyed `join` must specify `failOnDuplicate: true` and
  `failOnMismatch: true` unless unmatched items are intentionally handled with
  `remainder: true`. Use `combine` only for cartesian products. When a group
  size is known, use `groupKey(...)/groupTuple()` or `groupTuple(size: ...)`
  to avoid blocking.
- Let modules that publish `versions` or `multiqc` through topic channels do
  so; collect them at the top-level workflow with `channel.topic(...)` rather
  than adding duplicate explicit `.mix()` wiring.
- Use nf-core resource labels from `conf/base.config` (`process_single`,
  `process_low`, `process_medium`, `process_high`, etc.) instead of inventing
  per-process defaults. Module-specific behavior belongs in
  `conf/modules/<tool>.config`.

## Change coupling

- A new tool needs its module/subworkflow wiring, `conf/modules/<tool>.config`,
  an nf-test scenario, and any relevant MultiQC configuration. Update
  `README.md`, `docs/usage.md`, `docs/output.md`, `CITATIONS.md`,
  `CHANGELOG.md`, and the SVG/PNG metro map when its functionality or outputs
  require it.
- A new variant caller must be registered in the applicable germline, somatic,
  or tumor-only dispatcher; emit `meta.variantcaller`; and be classified in
  `subworkflows/local/post_variantcalling/main.nf` so it is not silently
  excluded from normalization, filtering, or consensus handling.
- Output changes require `docs/output.md` and `CHANGELOG.md`. Changelog
  entries use PR numbers and remain in ascending PR-number order within a
  section.

## Git workflow

Develop from `origin/dev` on `fix/issue-XXXX` or `feat/issue-XXXX` branches,
and target PRs to `dev`. Keep branches local and do not push or amend commits
unless explicitly requested.
