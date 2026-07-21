# RNABioInfo/mcaat-nf: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of RNABioInfo/mcaat-nf, created with the [nf-core](https://nf-co.re/) template
(nf-core/tools 4.0.2). This pipeline is built with the nf-core template and is not an nf-core
pipeline.

### `Added`

- Read QC subworkflow: FastQC (raw), fastp adapter and quality trimming in strict paired mode,
  optional Bowtie2 host depletion, optional BBDuk phiX scrub, optional deterministic
  `seqtk sample` subsampling, FastQC (trimmed) and `seqkit stats`. The phiX174 reference
  (NC_001422.1, public domain) is shipped as `assets/phix174.fasta.gz`, so `--remove_phix`
  requires no additional download.
- R1/R2 pair-sync gate. Samples whose mates disagree on record count are excluded before dispatch
  and reported in a named warning. MCAAT treats two input files as a strictly synchronised
  paired-end library; the vendored megahit aborts on an odd combined read count and writes nothing
  to stderr.
- `MCAAT_RUN`: CRISPR array detection in un-assembled reads. The wrapper applies:
  - `< /dev/null` on every invocation;
  - bash-derived `--threads`;
  - a dedicated `--output-folder`;
  - `--input-files` placed last on the command line;
  - a log-content exit-code remap;
  - post-condition assertions on `parameters.json` and `graph.sdbg_info`;
  - zero-result normalisation;
  - a verified megahit-intermediate prune;
  - graph retention controlled by `--mcaat_keep_graph` alone, independent of the optional-stage
    toggles.
- `MCAAT_PARSEARRAYS`: parses `CRISPR_Arrays_*.txt` into a tidy TSV, a spacer FASTA and a
  per-sample `*.provenance.tsv`. It reads `parameters.json` back and asserts that every flag the
  pipeline set was honoured. MCAAT's argument parser ignores unknown flags without reporting them.
- Optional, experimental, off-by-default `--run_phage_curation` stage (`MCAAT_PHAGECURATE`),
  including a FASTA header rewrite that makes records addressable across a cohort.
- Graph re-entry: an optional `sdbg` samplesheet column that skips array detection and feeds a
  previously published SDBG to the optional phage stage.
- Cohort layer: `cohort_arrays.tsv`, `cohort_summary.tsv`, `cohort_spacer_sharing.tsv`,
  `cohort_spacer_redundancy.tsv`, `cohort_repeat_families.tsv` and a de-duplicated cohort spacer
  FASTA. `COHORT_SUMMARISE` consumes the per-sample `*.arrays.tsv` and `*.provenance.tsv` tables
  directly.
- MultiQC report with zero-config custom content: cohort summary table, per-sample array and spacer
  counts, spacer- and repeat-length distributions, and the cohort sharing heatmap.
- Pipeline-owned MCAAT container, built without `-march=native` and without an `ENTRYPOINT`, with
  `libbz2` and `procps` present, pinned to a release tag.
- Test profiles `test`, `test_phage`, `test_graph` and `test_full`. Hand-written nf-test suites:
  `tests/pipeline/test_default.nf.test`, `tests/pipeline/test_phage.nf.test` and
  `tests/pipeline/test_graph_reentry.nf.test`. CI runs on Nextflow 25.10.4 and exercises both the
  docker and singularity profiles.

### `Changed`

- The pipeline is named `mcaat-nf` throughout. The manifest is `RNABioInfo/mcaat-nf`, the home page
  is `https://github.com/RNABioInfo/mcaat-nf`, the MultiQC section ids are `mcaat_*`, the
  environment-variable prefix is `MCAAT_`, and the test-data path segment is `mcaat-nf`. The DSL2
  workflow remains `MCAAT` in `workflows/mcaat.nf` and the `bin/mcaat_*.py` script names are
  unchanged. References to the earlier working name of this repository have been removed.
- The pipeline is licensed GPL-3.0, matching MCAAT. The nf-core template portions it reuses remain
  MIT-licensed by their authors. See the licence section of `README.md`.
- `params.pipelines_testdata_base_path` now defaults to
  `https://raw.githubusercontent.com/RNABioInfo/mcaat-nf-testdata/main/`.

### `Removed`

### `Fixed`

- `COHORT_SUMMARISE` was called with three mismatched inputs, including a pre-concatenated cohort
  table and every sample's raw `parameters.json`. Every sample emits a file with that exact name,
  so the inputs collide on staging as soon as a cohort has two samples. The process now takes the
  per-sample `*.arrays.tsv` and `*.provenance.tsv` collections, and
  `MCAAT_PARSEARRAYS.out.provenance` is consumed.
- `params.cohort_min_spacer_redundancy` was read by `COHORT_SUMMARISE` but declared in neither
  `nextflow.config` nor `nextflow_schema.json`, and resolved to `null` without a warning. It is now
  declared in both, with a default of `2`.
- Publish-path collision in `<outdir>/reports`: `CSVTK_CONCAT` and `COHORT_SUMMARISE` both wrote a
  file called `cohort_arrays.tsv`. They now publish under distinct names.
- `conf/modules.config` carried `withName:` selectors for processes that do not exist
  (`MCAAT_DETECT`, `MCAAT_PARSE`, `MCAAT_AGGREGATE`). They were removed; every remaining selector
  matches a real process.
- The default-profile nf-test asserted spacer FASTA headers against a pattern anchored with `$`
  immediately after the identifier, which cannot match: `bin/mcaat_parse_arrays.py` writes trailing
  ` repeat_len=.. spacer_len=..` attributes on every record. The assertion was corrected; the
  parser was not changed.
- `conf/test_full.config` stated that `assets/phix174.fasta.gz` was not shipped. The file is
  shipped; the comment and the `remove_phix` setting were corrected.

### `Dependencies`

| Dependency  | Version   |
| ----------- | --------- |
| `Nextflow`  | `25.10.4` |
| `MCAAT`     | `1.0.1`   |
| `bbmap`     | `39.18`   |
| `bowtie2`   | `2.5.4`   |
| `coreutils` | `9.5`     |
| `csvtk`     | `0.37.0`  |
| `fastp`     | `1.1.0`   |
| `fastqc`    | `0.12.1`  |
| `multiqc`   | `1.35`    |
| `python`    | `3.12`    |
| `seqkit`    | `2.13.0`  |
| `seqtk`     | `1.4`     |

The MCAAT version is the tag of the pinned image in `params.mcaat_container`
(`docker.io/feeka94/mcaat:1.0.1`); MCAAT has no `--version` flag, and the reported version comes
from the `MCAAT_VERSION` environment variable baked into that image. The `coreutils` row is the
version `CAT_FASTQ` reports as `cat`. Every row other than `MCAAT` and `python` is the version
pinned in the corresponding nf-core module's `environment.yml`, and the container images are the
matching Wave builds. `python` is pinned in the `container` directives of the two local
Python modules, `modules/local/mcaat/parsearrays/main.nf` and `modules/local/cohort_summarise/main.nf`.

Some images are multi-tool. The `bbmap` image also carries `pigz` 2.8. The `bowtie2` image also
carries `pigz` 2.8, `htslib` 1.21 and `samtools` 1.21. The `cat/fastq` image carries `coreutils` 9.5,
`grep` 3.11, `gzip` 1.13, `lbzip2` 2.5, `sed` 4.8 and `tar` 1.34. The modules do not report these
co-packaged versions. To read them from the image directly, run `pigz --version`,
`samtools --version` or `gzip --version`.

### `Deprecated`

Not applicable to the initial release.

### `Known issues`

- The test-data repository does not exist yet. Every test profile (`test`, `test_phage`,
  `test_graph`, `test_full`) resolves its samplesheet from `params.pipelines_testdata_base_path`,
  which defaults to `https://raw.githubusercontent.com/RNABioInfo/mcaat-nf-testdata/main/`. That
  repository must be created and populated with `mcaat-nf/samplesheet/*.csv`,
  `mcaat-nf/reads/*.fastq.gz` and `mcaat-nf/graph/` before `-profile test` can run. Until then the
  test profiles fail at input resolution. See `docs/usage.md` for the required layout. Pointing
  `--pipelines_testdata_base_path` at a local directory is a valid interim workaround.
- Module-level test coverage is missing. The hand-written suites in this repository are the three
  files under `tests/pipeline/`. The remaining `*.nf.test` files are vendored nf-core module and
  subworkflow tests, excluded by the `ignore` list in `nf-test.config`. None of `MCAAT_RUN`,
  `MCAAT_PARSEARRAYS`, `MCAAT_PHAGECURATE` or `COHORT_SUMMARISE` has a `tests/` directory, so the
  workarounds those modules carry — the exit-code remap, the `< /dev/null` stdin guard, the megahit
  intermediate prune and the zero-result normalisation — are exercised only through the end-to-end
  runs.
- `linux/amd64` only. MCAAT's vendored megahit rank/select relies on `__builtin_popcount` lowering
  and spoa pulls `<immintrin.h>`; there is no arm64 build. The `arm` profile runs the images under
  emulation.
- `--ram` had no effect in MCAAT 1.0.0: a gigabyte figure was passed into megahit's byte-valued
  `--host_mem`. This is fixed in the pinned 1.0.1 container. Size the Nextflow `memory` directive
  empirically; see `docs/usage.md`.
- Changing `--mcaat_keep_graph` invalidates every cached `MCAAT_RUN` task. Set it once per campaign.
- The phage / protospacer stage is experimental (MCAAT v2.0.0 track), off by default, and
  configured to retry twice and then be ignored.
