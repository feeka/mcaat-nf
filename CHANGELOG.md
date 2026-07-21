# RNABioInfo/mcaat-nf: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of RNABioInfo/mcaat-nf, created with the [nf-core](https://nf-co.re/)
template (nf-core/tools 4.0.2). This pipeline is built **with** the nf-core template
but is **not** an nf-core pipeline.

### `Added`

- Read QC subworkflow: FastQC (raw), fastp adapter/quality trimming in strict paired
  mode, optional Bowtie2 host depletion, optional BBDuk phiX scrub, optional
  deterministic `seqtk sample` subsampling, FastQC (trimmed) and `seqkit stats`.
  The phiX174 reference (NC_001422.1, public domain) is shipped as
  `assets/phix174.fasta.gz`, so `--remove_phix` works out of the box.
- **R1/R2 pair-sync gate**: samples whose mates disagree on record count are excluded
  before dispatch with a named warning. MCAAT treats two input files as a strictly
  synchronised paired-end library and the vendored megahit aborts on an odd combined
  read count with a completely empty stderr, so this gate converts the single most
  likely real-world failure from undebuggable into explicit.
- `MCAAT_RUN`: CRISPR array detection in un-assembled reads. Hardened wrapper with
  `< /dev/null` on every invocation, bash-derived `--threads`/`--ram`, a dedicated
  `--output-folder`, `--input-files` placed last, a log-content exit-code remap,
  post-condition assertions on `parameters.json` and `graph.sdbg_info`, zero-result
  normalisation, a verified megahit-intermediate prune, and toggle-invariant graph
  retention.
- `MCAAT_PARSEARRAYS`: parses `CRISPR_Arrays_*.txt` into a tidy TSV, a spacer FASTA
  and a per-sample `*.provenance.tsv`, and reads `parameters.json` back to assert that
  every flag the pipeline set was actually honoured (MCAAT's argument parser silently
  ignores unknown flags).
- Optional, **experimental**, off-by-default `--run_phage_curation` stage
  (`MCAAT_PHAGECURATE`), including a FASTA header rewrite so records are addressable
  in a cohort.
- Graph re-entry: an optional `sdbg` samplesheet column that skips array detection and
  feeds a previously published SDBG straight to the optional phage stage.
- Cohort layer: `cohort_arrays.tsv`, `cohort_summary.tsv`, `cohort_spacer_sharing.tsv`,
  `cohort_spacer_redundancy.tsv`, `cohort_repeat_families.tsv` and a de-duplicated
  cohort spacer FASTA. `COHORT_SUMMARISE` consumes the per-sample `*.arrays.tsv` and
  `*.provenance.tsv` tables directly.
- MultiQC report with zero-config custom content: a cohort summary table, per-sample
  array/spacer counts, spacer- and repeat-length distributions, and the cohort sharing
  heatmap.
- Pipeline-owned MCAAT container built without `-march=native` and without an
  `ENTRYPOINT`, with `libbz2` and `procps` present, and pinned to a release tag.
- Test profiles `test`, `test_phage`, `test_graph` and `test_full`; one pipeline-level
  nf-test suite (`tests/pipeline/test_default.nf.test`) — per-module suites for the four
  local modules are **not yet written**; CI on Nextflow 25.10.4 exercising **both** the
  docker and singularity profiles.

### `Changed`

- **The pipeline is named `mcaat-nf`**, everywhere and without exception. The manifest is
  `RNABioInfo/mcaat-nf`, the home page is `https://github.com/RNABioInfo/mcaat-nf`, the
  MultiQC section ids are `mcaat_*`, the environment-variable prefix is `MCAAT_`, and the
  test-data path segment is `mcaat-nf`. The DSL2 workflow itself remains `MCAAT` in
  `workflows/mcaat.nf`, and the `bin/mcaat_*.py` script names are unchanged. Every
  reference to the earlier working name of this repository has been removed.
- **The pipeline is licensed GPL-3.0**, matching MCAAT. The nf-core template portions
  it reuses remain MIT-licensed by their authors. See the licence section of
  `README.md`.
- `params.pipelines_testdata_base_path` now defaults to
  `https://raw.githubusercontent.com/RNABioInfo/mcaat-nf-testdata/main/`.

### `Removed`

### `Fixed`

- `COHORT_SUMMARISE` was called with three mismatched inputs — a pre-concatenated
  cohort table and every sample's raw `parameters.json`. The last of
  those collides on staging as soon as a cohort has two samples, because every sample
  emits a file with that exact name. It now takes the per-sample `*.arrays.tsv` and
  `*.provenance.tsv` collections, and `MCAAT_PARSEARRAYS.out.provenance` is finally
  consumed.
- `params.cohort_min_spacer_redundancy` was read by `COHORT_SUMMARISE` but declared in
  neither `nextflow.config` nor `nextflow_schema.json`, so it silently resolved to
  `null`. It is now declared in both, defaulting to `2`.
- Publish-path collision in `<outdir>/reports`: `CSVTK_CONCAT` and `COHORT_SUMMARISE`
  both wrote a file called `cohort_arrays.tsv`. They now publish under distinct names.
- `conf/modules.config` carried `withName:` selectors for processes that do not exist
  (`MCAAT_DETECT`, `MCAAT_CAS`, `MCAAT_PARSE`, `MCAAT_AGGREGATE`). Removed; every
  remaining selector matches a real process.
- The default-profile nf-test asserted spacer FASTA headers against a pattern anchored
  with `$` immediately after the identifier, which cannot match — `bin/mcaat_parse_arrays.py`
  writes trailing ` repeat_len=.. spacer_len=..` attributes on every record. The
  assertion was wrong, not the parser; the assertion was fixed.
- `conf/test_full.config` claimed `assets/phix174.fasta.gz` was not shipped. It is
  shipped, and the comment and the `remove_phix` setting have been revisited
  accordingly.

### `Dependencies`

| Dependency | Version   |
| ---------- | --------- |
| `Nextflow` | `25.10.4` |
| `MCAAT`    | `1.0.0`   |
| `bbmap`    | `39.18`   |
| `bowtie2`  | `2.5.4`   |
| `csvtk`    | `0.37.0`  |
| `fastp`    | `1.1.0`   |
| `fastqc`   | `0.12.1`  |
| `multiqc`  | `1.35`    |
| `python`   | `3.12`    |
| `seqkit`   | `2.13.0`  |
| `seqtk`    | `1.4`     |

Every row above is the version pinned in the corresponding module's `environment.yml`; the
container images are the matching Wave builds. Two of those images are multi-tool: the `bbmap`
image also carries `pigz`, and the `bowtie2` image also carries `pigz`, `htslib` and `samtools`.
Those co-packaged versions are fixed by the image digest rather than declared in the recipe, and
the modules do not report them, so they are not listed here. If you need them for a methods
section, read them out of the image itself (`pigz --version`, `samtools --version`).

### `Deprecated`

Not applicable — initial release.

### `Known issues`

- **The test-data repository does not exist yet.** Every test profile (`test`,
  `test_phage`, `test_graph`, `test_full`) resolves its samplesheet from
  `params.pipelines_testdata_base_path`, which defaults to
  `https://raw.githubusercontent.com/RNABioInfo/mcaat-nf-testdata/main/`. That
  repository **must be created and populated** — with `mcaat-nf/samplesheet/*.csv`,
  `mcaat-nf/reads/*.fastq.gz` and `mcaat-nf/graph/` — before `-profile test` can run at
  all. Until then the test profiles fail at input resolution. This is stated here
  rather than papered over; see `docs/usage.md` for the required layout. Pointing
  `--pipelines_testdata_base_path` at a local directory is a valid interim workaround.
- **Module-level test coverage is missing.** The only hand-written suite in this repository
  is `tests/pipeline/test_default.nf.test`; the five other `*.nf.test` files are vendored
  nf-core subworkflow tests, and `nf-test.config`'s `ignore` list excludes them. None of
  `MCAAT_RUN`, `MCAAT_PARSEARRAYS`, `MCAAT_PHAGECURATE` or `COHORT_SUMMARISE` has a
  `tests/` directory, so the workarounds those modules carry — the exit-code remap, the
  `< /dev/null` stdin guard, the megahit intermediate prune and the zero-result
  normalisation — are exercised only indirectly, through the one end-to-end run. Writing
  per-module suites for those four is the highest-value test work outstanding.
- The `conda` profile cannot run the array or phage stages: MCAAT has no Bioconda
  recipe, so those processes declare no `conda` directive.
  The `conda` profile is consequently excluded from the CI matrix.
- `linux/amd64` only. MCAAT's vendored megahit rank/select relies on
  `__builtin_popcount` lowering and spoa pulls `<immintrin.h>`; there is no arm64 build.
- `--ram` was functionally inert in MCAAT 1.0.0 (a gigabyte figure was passed into
  megahit's byte-valued `--host_mem`). This is fixed in the pinned 1.0.1 container. The
  Nextflow `memory` directive should still be sized empirically; see `docs/usage.md`.
- `--mcaat_keep_graph` is the one documented cache-busting toggle. Changing it
  invalidates every cached `MCAAT_RUN` task. Set it once per campaign.
- The phage / protospacer stage is experimental (MCAAT v2.0.0 track), off by default,
  and configured to be ignored after two failed attempts.
