# RNABioInfo/mcaat-nf: Usage

## Table of contents

- [Introduction](#introduction)
- [Samplesheet input](#samplesheet-input)
- [Running the pipeline](#running-the-pipeline)
- [Parameters](#parameters)
  - [Input/output](#inputoutput)
  - [Read QC](#read-qc)
  - [MCAAT array detection](#mcaat-array-detection)
  - [Stage activation](#stage-activation)
  - [Phage / protospacer curation](#phage--protospacer-curation)
  - [Reporting](#reporting)
  - [Fixed MCAAT settings](#fixed-mcaat-settings)
- [Graph re-entry](#graph-re-entry)
- [Test profiles and test data](#test-profiles-and-test-data)
- [Resources and tuning for large metagenomes](#resources-and-tuning-for-large-metagenomes)
- [Resume behaviour and caching](#resume-behaviour-and-caching)
- [Exit codes and troubleshooting](#exit-codes-and-troubleshooting)
- [Core Nextflow arguments](#core-nextflow-arguments)
- [Custom configuration](#custom-configuration)

## Introduction

This pipeline runs [MCAAT](https://github.com/RNABioInfo/mcaat) over quality-controlled metagenomic
reads to find CRISPR repeat–spacer arrays without assembly, then parses and aggregates the results
across a cohort.

Requirements:

| Component        | Requirement                                                             |
| ---------------- | ----------------------------------------------------------------------- |
| Nextflow         | `>=25.10.4`                                                             |
| Container engine | Docker, Singularity / Apptainer, Podman or Charliecloud                 |

MCAAT has no Bioconda package and its processes declare no Conda environment. The `conda` profile
inherited from the nf-core template runs the QC and reporting layers only; it cannot run the
array-detection or phage stages.

## Samplesheet input

The samplesheet is a comma-separated file with a header row, passed with `--input`.

```bash
--input '[full path to samplesheet file]'
```

### Columns

| Column    | Required | Description                                                                                                                                             |
| --------- | -------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | yes      | Sample identifier. No spaces. Rows sharing a value are concatenated before analysis. Used verbatim in every output path and filename.                    |
| `run`     | no       | Sequencing run or lane identifier. No spaces. Keeps multiple rows of the same sample distinct; dropped once the rows have been merged.                   |
| `fastq_1` | yes\*    | Path to the R1 FASTQ (or the single-end FASTQ). Extension must be `.fastq.gz`, `.fq.gz`, `.fastq` or `.fq`.                                              |
| `fastq_2` | no       | Path to the R2 FASTQ. Same extensions. An empty value marks the sample as single-end. Cannot be supplied without `fastq_1`.                              |
| `sdbg`    | no       | Path to a pre-built MCAAT graph directory, for re-analysis without rebuilding it. Mutually exclusive with `fastq_1`. See [Graph re-entry](#graph-re-entry). |

\* Each row must supply either `fastq_1` or `sdbg`.

The pair `sample` + `run` must be unique across the samplesheet.

### Multiple runs of the same sample

Give the rows the same `sample` value and different `run` values. R1 files are concatenated with
each other and R2 files with each other, in samplesheet order, before any QC step:

```csv
sample,run,fastq_1,fastq_2,sdbg
CONTROL_REP1,L001,AEG588A1_S1_L001_R1_001.fastq.gz,AEG588A1_S1_L001_R2_001.fastq.gz,
CONTROL_REP1,L002,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,
CONTROL_REP1,L003,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,
```

All rows for one sample must agree on single- vs paired-endedness. A mixed set stops the run at
launch with an error naming the sample.

### Full example

```csv
sample,run,fastq_1,fastq_2,sdbg
gut_A,L001,/data/gut_A_L001_R1.fastq.gz,/data/gut_A_L001_R2.fastq.gz,
gut_A,L002,/data/gut_A_L002_R1.fastq.gz,/data/gut_A_L002_R2.fastq.gz,
soil_B,,/data/soil_B.fastq.gz,,
rerun_C,,,,/archive/mcaat_run1/arrays/rerun_C/graph
```

### Paired-end input handling

MCAAT accepts one or two input files and treats two files as a strictly synchronised paired-end
library: record _n_ of file 1 is taken as the mate of record _n_ of file 2. If the two files hold
different numbers of records, the graph builder aborts. Its stderr is redirected to `/dev/null`
inside MCAAT, so `.command.err` is empty for that failure.

The pipeline therefore:

- runs `fastp` in strict paired mode with `--detect_adapter_for_pe`, and never routes failed or
  orphan reads onward (`--fastp_save_trimmed_fail` writes them to a separate published file only);
- never produces merged reads;
- subsamples both mates with the same `seqtk` seed;
- asserts, from the `num_seqs` column of `seqkit stats`, that R1 and R2 record counts are identical
  before MCAAT is dispatched. Samples failing this check are excluded from the run with a
  `MCAAT-NF: EXCLUDING sample` message, and the rest of the cohort continues.

`--skip_pairsync_check` disables the assertion. `--skip_qc` bypasses the whole QC subworkflow,
including the assertion.

## Running the pipeline

```bash
nextflow run RNABioInfo/mcaat-nf \
   --input ./samplesheet.csv \
   --outdir ./results \
   -profile docker
```

This writes results to `./results`, plus `work/` (intermediates) and `.nextflow*` (logs and cache)
in the launch directory.

To supply parameters from a file:

```bash
nextflow run RNABioInfo/mcaat-nf -profile docker -params-file params.yaml
```

```yaml title="params.yaml"
input: "./samplesheet.csv"
outdir: "./results"
subsample_reads: 20000000
remove_phix: true
```

> [!WARNING]
> Pipeline parameters set in a `-c` config file are not validated and interact with profile
> precedence. Pass parameters on the CLI or through `-params-file`.

### Updating and version pinning

```bash
nextflow pull RNABioInfo/mcaat-nf                    # update to the latest version
nextflow run RNABioInfo/mcaat-nf -r 1.0.0 ...        # pin an exact release
```

Record the `-r` value used for a run. The pipeline writes the full resolved parameter set to
`results/pipeline_info/params_<suffix>.json` on every run.

## Parameters

Parameters are validated against `nextflow_schema.json` at launch. Run
`nextflow run RNABioInfo/mcaat-nf --help` for the short parameter list, `--help_full` for all
parameters including hidden ones.

### Input/output

| Parameter       | Type   | Default | Description                                |
| --------------- | ------ | ------- | ------------------------------------------ |
| `input`         | path   | —       | **Required.** Samplesheet CSV.             |
| `outdir`        | path   | —       | **Required.** Results directory.           |
| `email`         | string | `null`  | Address to send the completion summary to. |
| `multiqc_title` | string | `null`  | Title used for the MultiQC report.         |

### Read QC

All QC steps are on by default and are switched off with `skip_*` flags.

| Parameter                 | Type              | Default                                 | Description                                                                                                                                                            |
| ------------------------- | ----------------- | --------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `skip_qc`                 | boolean           | `false`                                 | Bypass the whole QC subworkflow and feed samplesheet reads straight to MCAAT. The R1/R2 pair-sync gate is part of this subworkflow and does not run when it is set.     |
| `skip_fastqc`             | boolean           | `false`                                 | Skip both FastQC instances (raw and trimmed).                                                                                                                          |
| `skip_fastp`              | boolean           | `false`                                 | Skip adapter and quality trimming.                                                                                                                                     |
| `fastp_save_trimmed_fail` | boolean           | `false`                                 | Write reads that failed the fastp filters to a separate published file. These reads are not fed to MCAAT.                                                               |
| `save_trimmed_reads`      | boolean           | `false`                                 | Publish the trimmed FASTQ files.                                                                                                                                       |
| `host_fasta`              | path              | `null`                                  | Host genome FASTA (local path or `https://` URL). Setting it enables host depletion with Bowtie2. The index is built once per run and reused for every sample.          |
| `host_bowtie2_index`      | path              | `null`                                  | Pre-built Bowtie2 index directory. Takes precedence over `host_fasta` and skips the index build. Supplying both emits a warning at launch.                              |
| `save_hostremoved_reads`  | boolean           | `false`                                 | Publish host-depleted FASTQ files.                                                                                                                                     |
| `remove_phix`             | boolean           | `false`                                 | Enable BBDuk k-mer scrubbing (`k=31 hdist=1`) against `phix_reference`.                                                                                                |
| `phix_reference`          | path              | `${projectDir}/assets/phix174.fasta.gz` | Contaminant FASTA for BBDuk. The phiX174 genome (NCBI NC_001422.1, 5386 bp) is shipped in `assets/`, so `--remove_phix` needs no download. BBDuk builds an in-memory k-mer table over the whole reference; a mammalian genome does not fit. |
| `subsample_reads`         | integer or number | `null`                                  | Read budget per mate, passed to `seqtk sample`. An integer takes that many reads; a float `0 < x < 1` takes that fraction. `null` disables subsampling.                 |
| `subsample_seed`          | integer           | `100`                                   | `seqtk sample` seed (`-s`). Applied identically to both mates, which preserves pairing through subsampling. Minimum 1.                                                  |
| `skip_pairsync_check`     | boolean           | `false`                                 | Disable the R1/R2 record-count assertion. Without it, a de-synchronised pair reaches MCAAT and the task fails with exit 1 and an empty `.command.err`.                  |

### MCAAT array detection

| Parameter                      | Type    | Default                          | Description                                                                                                                                                              |
| ------------------------------ | ------- | -------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `mcaat_min_spacer_len`         | integer | `23`                             | Minimum spacer length (`--min-spacer-len`), applied as a post-processing filter. Minimum 1.                                                                              |
| `mcaat_max_spacer_len`         | integer | `50`                             | Maximum spacer length (`--max-spacer-len`). Minimum 1.                                                                                                                   |
| `mcaat_min_repeat_len`         | integer | `23`                             | Minimum repeat length (`--min-repeat-len`). Schema minimum 23: the de Bruijn graph is built with a hard-coded _k_ of 23, so shorter repeats are not representable.        |
| `mcaat_max_repeat_len`         | integer | `55`                             | Maximum repeat length (`--max-repeat-len`). Minimum 23.                                                                                                                  |
| `mcaat_cycle_min_length`       | integer | `27`                             | Minimum multicycle length in graph nodes (`--cycle-min-length`). Minimum 1.                                                                                              |
| `mcaat_cycle_max_length`       | integer | `77`                             | Maximum multicycle length in graph nodes (`--cycle-max-length`). Minimum 1.                                                                                              |
| `mcaat_threshold_multiplicity` | integer | `20`                             | Multiplicity a candidate cycle must reach to be reported (`--threshold-multiplicity`). Lower values recover rarer array-bearing organisms and increase runtime. Minimum 1. |
| `mcaat_low_abundance`          | boolean | `true`                           | Low-abundance mode (`--low-abundance`). Rendered as the lowercase string `true` or `false`; MCAAT accepts only `1`, `true` or `yes` and reads anything else as false.     |
| `mcaat_keep_graph`             | boolean | `true`                           | Keep the graph in the task work directory after the run. Required by the phage stage and by `save_mcaat_graph`. Changing it invalidates cached `MCAAT_RUN` tasks — see [Resume behaviour](#resume-behaviour-and-caching). |
| `prune_graph`                  | boolean | `true`                           | Delete graph-builder intermediates that MCAAT never reads back (`outfile_prefix.bin`, `outfile_prefix.lib_info`, `data.lib`, `graph.mercy_cand.*`). Set `false` for a directory other MEGAHIT-derived tools can consume; it is roughly 3–4× larger. |
| `save_mcaat_graph`             | boolean | `false`                          | Publish the (pruned) graph under `arrays/<sample>/graph/`. Requires `mcaat_keep_graph`.                                                                                   |
| `graph_publish_mode`           | string  | `symlink`                        | `publishDir` mode used for the graph directory only. One of `symlink`, `copy`, `link`, `rellink`.                                                                        |
| `save_mcaat_cycles`            | boolean | `false`                          | Publish `cycles.txt`. One line per cycle, 23 bp per node; the file can be enormous, which is why publishing is off by default. Its line order comes from an unordered-map iteration and is not byte-reproducible between runs. |
| `mcaat_scratch`                | string  | `null`                           | Value for the Nextflow `scratch` directive on all MCAAT processes, e.g. `'$TMPDIR'` or an absolute path. Relocates every MCAAT intermediate to that storage.              |
| `mcaat_args`                   | string  | `null`                           | Raw arguments appended to the `mcaat` command line before `--input-files`. MCAAT ignores unrecognised arguments without a warning, so a typo here has no visible effect; the `parameters.json` readback below is the detector. Setting it emits a warning at launch. |
| `mcaat_container`              | string  | `docker.io/feeka94/mcaat:1.0.1`  | Container image providing the `mcaat` binary. Built from this repository's `Dockerfile.local` with `-march=x86-64-v2`, `spoa_optimize_for_native=OFF`, no `ENTRYPOINT`, and `ENV MCAAT_VERSION` set. |

#### The `parameters.json` readback

MCAAT's argument parser has no unknown-argument branch: a misspelled flag is discarded without a
warning and the run proceeds with defaults. The array-parsing step reads MCAAT's own
`parameters.json` back and asserts that the resolved thread count, multiplicity threshold, cycle
bounds, spacer/repeat bounds and autoclean setting match the values the pipeline set. A mismatch
fails the task and prints expected versus observed. The readback can only check flags the pipeline
itself sets.

The same step writes the values it read into `<sample>.provenance.tsv`, which the cohort layer
consumes so that the resolved MCAAT settings become columns of `cohort_summary.tsv`.

### Stage activation

| Parameter            | Type    | Default | Description                                                                             |
| -------------------- | ------- | ------- | ----------------------------------------------------------------------------------------- |
| `run_phage_curation` | boolean | `false` | Enable the protospacer / phage curation stage. Experimental. Requires `mcaat_keep_graph`. |

### Phage / protospacer curation

The stage belongs to MCAAT's unreleased v2.0.0 track and is off by default.

| Parameter        | Type    | Default | Description                                    |
| ---------------- | ------- | ------- | ---------------------------------------------- |
| `save_phage_log` | boolean | `false` | Publish the gzipped `phage-curate` stdout log. |

Every other setting in this MCAAT subcommand is hard-coded (search depth 50, overridden per cycle;
contig length bounds 500–2500 bp) or accepted and never read (`--beam-width`).

The underlying search is a recursive depth-first traversal with no candidate cap and no wall-clock
budget, followed by two O(n²) `std::search` subpath filters, so wall time is the binding constraint
rather than memory. `MCAAT_PHAGECURATE` retries twice and is then ignored, so one sample cannot fail
the cohort.

The reconstructed paths are graph walks around a spacer match. The pipeline provides no independent
evidence that a walk corresponds to a contiguous element.

### Reporting

| Parameter                      | Type    | Default | Description                                                                                             |
| ------------------------------ | ------- | ------- | --------------------------------------------------------------------------------------------------------- |
| `skip_cohort_aggregation`      | boolean | `false` | Skip the cross-sample layer (`CSVTK_CONCAT`, `SEQKIT_RMDUP`, `COHORT_SUMMARISE`).                        |
| `skip_multiqc`                 | boolean | `false` | Skip the MultiQC report.                                                                                |
| `cohort_min_spacer_sharing`    | integer | `2`     | Minimum number of shared spacers for a sample pair to appear in the sharing matrix. Minimum 0; `0` emits every pair, giving an N² table. |
| `cohort_min_spacer_redundancy` | integer | `2`     | Minimum number of samples a spacer sequence must occur in to appear in `cohort_spacer_redundancy.tsv`. Minimum 0; `1` emits every spacer in the study. |
| `multiqc_config`               | path    | `null`  | Extra MultiQC config, merged with the pipeline's own.                                                   |
| `multiqc_logo`                 | path    | `null`  | Custom logo for the report.                                                                             |
| `multiqc_methods_description`  | path    | `null`  | Custom methods-description YAML.                                                                        |

Both cohort thresholds are rendered onto the `COHORT_SUMMARISE` command line as
`--min-spacer-sharing` / `--min-spacer-redundancy`. `ext.args` is appended after them, so a config
can override either value.

### Fixed MCAAT settings

| Setting                                    | Behaviour                                                                                                                                                                                              |
| ------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| k-mer size                                 | Hard-coded to 23 in the graph build. `mcaat_min_repeat_len` cannot go below 23 for the same reason.                                                                                                     |
| minimum k-mer frequency                    | Hard-coded to 1, so no error-k-mer filtering takes place. Resource use tracks total sequenced bases.                                                                                                    |
| `--ram`                                    | Derived inside the task at runtime as 90 % of the cgroup memory limit, with a 2 GB floor, falling back to `/proc/meminfo` when no cgroup limit is set. It is not interpolated from `task.memory`, so retuning the `memory` directive does not change the task hash. The clamp keeps MCAAT's `ram <= 1.0` and `ram > total system RAM` validation branches unreachable; both route into an interactive stdin prompt. |
| `--threads`                                | Derived from `task.cpus` and clamped to `nproc`. MCAAT rejects a thread count above `std::thread::hardware_concurrency()` and routes that failure into the interactive prompt.                          |
| `--autoclean`                              | Always `false`. The default is `true` and recursively deletes `<out>/graph`. Cleanup is done by the pipeline in shell, so the resolved command text does not change when a stage toggle is flipped.     |
| `--output-folder`                          | Always `mcaat`, a subdirectory inside the task directory. The MCAAT default is a localtime-derived `mcaat_run_<timestamp>` name.                                                                        |
| `--settings`                               | Not used. A settings file's `input-files` key concatenates with the command-line one instead of replacing it, which can turn a single-end run into a mis-paired paired-end run.                          |
| `--benchmark`                              | Accepted and stored, but its consumer is never called.                                                                                                                                                 |
| phage `--beam-width`, depth, length bounds | Accepted but never read, or hard-coded.                                                                                                                                                                |
| stdin                                      | All MCAAT invocations read from `/dev/null`. On a settings-validation failure MCAAT prompts `Remove it? (y/n):` and reads an uninitialised char; one branch calls `fs::remove_all()` on the output folder. With a tty or a pipe on stdin the job blocks holding a scheduler slot. |

Launch-time validation rejects `mcaat_min_spacer_len > mcaat_max_spacer_len`,
`mcaat_min_repeat_len > mcaat_max_repeat_len`, `mcaat_cycle_min_length > mcaat_cycle_max_length`,
`run_phage_curation` without `mcaat_keep_graph`, and `save_mcaat_graph` without `mcaat_keep_graph`.

## Graph re-entry

Building the de Bruijn graph is the most expensive step. An existing graph from a previous published
run can be re-analysed with the phage stage without rebuilding:

```csv
sample,run,fastq_1,fastq_2,sdbg
rerun_C,,,,/archive/mcaat_run1/arrays/rerun_C/graph
```

```bash
nextflow run RNABioInfo/mcaat-nf \
   --input samplesheet_graph.csv \
   --outdir results_rerun \
   --run_phage_curation \
   -profile singularity
```

Rules:

- The row must supply `sdbg` and must not supply `fastq_1`.
- The directory must contain `graph.sdbg_info` and at least one `graph.sdbg.<N>` shard. MCAAT loads
  the graph without validating the path, so the pipeline checks it and fails with a named message
  when the file is missing or empty.
- The stage also needs that sample's arrays. They are located as the sibling of the graph directory:
  `<sdbg>/../CRISPR_Arrays_*.txt`. That is the layout a previous run publishes under
  `arrays/<sample>/`, so pointing `sdbg` at `…/arrays/<sample>/graph` resolves them. When no array
  files are found there, the pipeline stops at launch and names the sample.
- Graph re-entry applies to `--run_phage_curation`. It does not re-run array detection, and a
  re-entry row contributes no rows to the cohort tables: it carries no `parameters.json` and no
  parsed array table of its own.

To produce re-usable graphs, run with `--save_mcaat_graph --graph_publish_mode copy`. The default
graph publish mode `symlink` points back into `work/` and does not survive `nextflow clean`.

## Test profiles and test data

> [!IMPORTANT]
> The test-data repository does not exist yet, so `-profile test` cannot run.

The four test profiles are `test`, `test_phage`, `test_graph` and `test_full`. Each builds its
`--input` from `params.pipelines_testdata_base_path`, whose default is:

```
https://raw.githubusercontent.com/RNABioInfo/mcaat-nf-testdata/main/
```

Until that repository is created and populated, launching any test profile fails while resolving the
samplesheet URL, before MCAAT is reached.

Profile settings:

| Profile      | Samplesheet                          | Notable parameters                                                                       |
| ------------ | ------------------------------------ | ---------------------------------------------------------------------------------------- |
| `test`       | `mcaat-nf/samplesheet/samplesheet.csv`       | `subsample_reads 2000`, `subsample_seed 100`; MCAAT defaults; resource limits 4 cpus / 15 GB / 1 h. |
| `test_phage` | `mcaat-nf/samplesheet/samplesheet.csv`       | Same input and MCAAT settings as `test`, plus `run_phage_curation true`, `mcaat_keep_graph true`, `save_phage_log true`; resource limits 4 cpus / 15 GB / 1 h. |
| `test_graph` | `mcaat-nf/samplesheet/samplesheet_graph.csv` | `run_phage_curation true`, `mcaat_keep_graph true`, `save_phage_log true`, `subsample_reads 2000`, `subsample_seed 100`; resource limits 4 cpus / 15 GB / 1 h. Samplesheet rows carry `sdbg` and no `fastq_1`. |
| `test_full`  | `mcaat-nf/samplesheet/samplesheet_full.csv`  | `subsample_reads null`, `remove_phix true`, `mcaat_keep_graph true`, `prune_graph true`, `save_mcaat_graph false`, `save_mcaat_cycles false`, `run_phage_curation false`, `multiqc_title 'mcaat-nf full-size test'`. Not run in PR CI; dispatched manually to obtain resource measurements. |

The layout the profiles expect, relative to the base path:

```
mcaat-nf/
├── samplesheet/
│   ├── samplesheet.csv          # used by -profile test and -profile test_phage
│   ├── samplesheet_graph.csv    # used by -profile test_graph (sdbg column populated)
│   └── samplesheet_full.csv     # used by -profile test_full
├── reads/
│   └── *.fastq.gz               # the FASTQ files the samplesheets point at
└── graph/
    └── <sample>/                # a pre-built SDBG plus its CRISPR_Arrays_*.txt sibling,
                                 # for the graph re-entry profile
```

Constraints on the `-profile test` reads: a few thousand read pairs, small enough for the GitHub
Actions budget, containing at least one CRISPR array so the run exercises the parser and the cohort
layer as well as the zero-result path. The spacer-sharing matrix and the MultiQC heatmap require at
least two samples; with a single sample they are skipped.

The profiles can be run against local files:

```bash
nextflow run RNABioInfo/mcaat-nf -profile test,docker \
   --pipelines_testdata_base_path /path/to/local/testdata/ \
   --outdir results_test
```

The trailing slash is required: the profiles concatenate the base path with
`mcaat-nf/samplesheet/…` directly.

## Resources and tuning for large metagenomes

### Read budget

MCAAT builds its graph with a minimum k-mer frequency of 1, so every sequencing error contributes
edges. Graph size, and therefore RAM, disk and runtime, scales with total sequenced bases rather
than with community complexity. Deeper sequencing of the same community costs proportionally more.

`--subsample_reads` controls that budget:

```bash
--subsample_reads 20000000     # 20 M reads per mate
--subsample_seed  100          # same seed for both mates
--subsample_reads 0.25         # 25 % of reads
```

Halving depth roughly halves the graph.

### Default resource labels

| Process             | Label                | cpus | memory (attempt 1) | time (attempt 1) |
| ------------------- | -------------------- | ---- | ------------------ | ---------------- |
| `MCAAT_RUN`         | `process_mcaat_run`  | 8    | 72 GB              | 16 h             |
| `MCAAT_PHAGECURATE` | `process_phage`      | 4    | 48 GB              | 24 h             |

`memory` and `time` are multiplied by `task.attempt` on each retry; `cpus` is fixed and does not
change on retry. The cycle finder allocates one visited bitset of `(nodes+63)/64` words per thread,
so memory grows linearly with core count and a higher core count on an out-of-memory retry increases
peak memory. Fixed `cpus` also pins the number of graph shard files, since the count of
`graph.sdbg.<N>` files equals `--num_cpu_threads`, which keeps `-resume` hashes and output
declarations stable across profiles.

Three phases are single-threaded regardless of `cpus`: `RecursiveReduction` over all tips, per-bucket
bitset zeroing, and the SPOA post-processing stage. Measured throughput plateaus at around 24 cores.

### Calibrating memory

The `memory` directive is the point at which the scheduler or cgroup kills the job. Peak resident
memory is approximately the larger of the two phases, which do not overlap:

```
build   ~  total_bases / 4  +  8 * n_reads          bytes
detect  ~  2 * n_edges  +  cpus * n_edges / 8       bytes
```

Procedure:

1. Run ten representative samples with `-with-report -with-trace`.
2. For each, run `grep 'Memory estimate' work/<hash>/.command.out` to read MCAAT's own visited-bitset
   figure.
3. Read `peak_rss` from `results/pipeline_info/execution_trace_*.txt`.
4. Set `memory` to the observed maximum plus about 30 %.

Blank trace columns indicate a container image without `ps` (the `procps` package). The
pipeline-owned image includes it; older MCAAT images do not.

```groovy title="custom_resources.config"
process {
    withLabel: process_mcaat_run {
        memory = { 128.GB * task.attempt }
        time   = { 24.h   * task.attempt }
    }
}
```

```bash
nextflow run RNABioInfo/mcaat-nf -profile singularity -c custom_resources.config ...
```

### Disk

`MCAAT_RUN` requests `50 GB + 4 × (input GB)`, scaled by `task.attempt`. The multiplier covers the
2-bit re-encoding of every input read written by the graph builder (`outfile_prefix.bin`, about
0.29 bytes per base, transient), the graph shards (about 2 bytes per edge), `cycles.txt`, and the
retained pruned graph.

Retained-graph sizing is about 2 GB of pruned graph per gigabase of input. A 100-sample cohort at
10 Gbp each occupies around 2 TB in `work/`. Setting `--mcaat_keep_graph false` removes that residue,
makes the phage stage unavailable, and invalidates cached `MCAAT_RUN` tasks — see
[Resume behaviour](#resume-behaviour-and-caching).

### Node-local scratch

Set `--mcaat_scratch '$TMPDIR'` to use node-local storage. MCAAT writes its `--output-folder`
relative to the task working directory, so relocating the task directory relocates every
intermediate. It is off by default because not every site provisions node-local storage.

### CPU architecture

The default image `docker.io/feeka94/mcaat:1.0.1` is compiled with `-march=x86-64-v2` and runs on
heterogeneous clusters. The MCAAT authors' `docker.io/feeka94/mcaat:1.0.0` image is compiled with
`-march=native` and can abort with `Illegal instruction` (exit 132) on a CPU older than its build
host. Exit 132 is not retried, since a retry onto the same node class fails identically.

The 1.0.0 image also declares `ENTRYPOINT ["mcaat"]`. Under Docker and Podman the MCAAT processes
pass `--entrypoint ""`, which keeps that image working; Singularity and Apptainer ignore
`ENTRYPOINT`.

MCAAT is linux/amd64 only: the vendored megahit rank/select code relies on `__builtin_popcount`
lowering and spoa pulls `<immintrin.h>`. The `arm` profile runs the images under emulation.

## Resume behaviour and caching

`-resume` works as usual. The specifics below apply to the MCAAT stages.

**Stage toggles do not change the array-detection cache.** The `MCAAT_RUN` command text contains no
reference to `run_phage_curation`, `save_mcaat_graph` or `save_mcaat_cycles`; those steer publishing
only. The following sequence costs one graph build:

```bash
# month 1
nextflow run RNABioInfo/mcaat-nf --input s.csv --outdir results -profile singularity

# month 6 — resumes every MCAAT_RUN from cache, runs only the phage chain
nextflow run RNABioInfo/mcaat-nf --input s.csv --outdir results -profile singularity \
    --run_phage_curation -resume
```

**`--mcaat_keep_graph` changes the cache.** Graph retention is performed in shell inside the task, so
the toggle changes the resolved command text and invalidates every cached `MCAAT_RUN`. Set it once
per campaign.

**Changing `cpus` changes the cache** and changes the number of graph shard files. Retuning
`process_mcaat_run` cpus partway through a cohort forces a rebuild.

**Other inputs to the task hash:** any parameter in the `mcaat_options` group, `mcaat_args`, the input
files themselves, and the container image tag or digest. The `memory` directive is not part of the
hash, because `--ram` is derived at runtime rather than interpolated from `task.memory`.

**`work/` and published graphs.** With the default `--graph_publish_mode symlink`,
`arrays/<sample>/graph/` is a symlink into `work/`, and `nextflow clean` leaves dangling links. Use
`--graph_publish_mode copy` for a published graph that outlives the work directory.

## Exit codes and troubleshooting

The MCAAT processes replace the nf-core default retry window (130–145 plus 104) with a narrower one:
`104`, `137`, `140`, `143` and `175`–`177` are retried, `maxRetries = 2`. MCAAT's release `main()`
aborts with SIGABRT (134) on every CLI and validation error, and a `CycleFinder` `std::bad_alloc`
aborts with the same status, so the module scripts disambiguate the two by grepping MCAAT's own log
and remap them before falling through to the raw status.

| Code                             | Meaning                                                | Retried | Action                                                                                                                    |
| -------------------------------- | ------------------------------------------------------ | ------- | ----------------------------------------------------------------------------------------------------------------------- |
| `64`                             | Deterministic MCAAT argument or validation error, MCAAT exited 0 without writing `parameters.json`, or the phage stage found no `graph.sdbg_info` | no | Read `mcaat.log`, published as `arrays/<sample>/<sample>.mcaat.log`. The message is there, not in `.command.err`.         |
| `65`                             | Graph build finished but produced no `graph.sdbg_info` | no      | Check FASTQ integrity and R1/R2 record parity. The graph builder's stderr is redirected to `/dev/null`, so `.command.err` is empty. |
| `137`                            | Allocation failure (`std::bad_alloc`) or OOM kill      | yes     | Memory doubles on retry. If it still fails, subsample the reads.                                                          |
| `132`                            | Illegal instruction                                    | no      | The binary was built for a newer CPU than the node. Use the default image, or pin the queue.                              |
| `134`                            | Abort matching neither log pattern                     | no      | Open an issue with `mcaat.log` attached.                                                                                  |
| `1`                              | Graph-builder fatal error                              | no      | `.command.err` is empty; inspect the inputs.                                                                              |
| `104`, `140`, `143`, `175`–`177` | Scheduler and walltime kills                           | yes     | Time doubles on retry.                                                                                                    |

`MCAAT_PHAGECURATE` retries twice and its third failure is ignored rather than failing the run.

**A sample is missing from the results.** Look for a `MCAAT-NF: EXCLUDING sample` line in the
Nextflow log: the pair-sync gate dropped it because R1 and R2 had different record counts.

**A sample reports zero arrays.** Zero arrays is a reported result. The sample gets a header-only
`CRISPR_Arrays_1.txt`, a zero row in `reports/cohort_summary.tsv`, and a `0` in the MultiQC General
Statistics table. MCAAT expresses the zero case three ways — no `CRISPR_Arrays_*.txt` when no spacers
were extracted, none when `cycles.txt` is unreadable, and a header-only file when all candidates were
filtered out — and the module normalises the first two into the third.

**`-profile test` fails immediately.** The test-data repository does not exist yet; see
[Test profiles and test data](#test-profiles-and-test-data).

**Locating the error message.** For the array-detection step it is in `mcaat.log` and
`.command.out`. MCAAT redirects the graph builder's stderr to `/dev/null`, so `.command.err` is empty
for graph-build failures.

## Core Nextflow arguments

> [!NOTE]
> These are Nextflow's own options and take a single hyphen.

### `-profile`

Container profiles: `docker`, `singularity`, `podman`, `apptainer`, `charliecloud`, `shifter`,
`wave`. The `arm` profile adds `--platform=linux/amd64` under Docker for emulated runs. The `conda`
and `mamba` profiles cannot run the MCAAT stages.

Test profiles: `test`, `test_phage`, `test_graph`, `test_full`. See
[Test profiles and test data](#test-profiles-and-test-data).

Profiles are applied left to right, so the container profile goes last: `-profile test,docker`.

### `-resume`

Restart from cached results. See [Resume behaviour](#resume-behaviour-and-caching). A specific
earlier run can be resumed by name or session ID; `nextflow log` lists them.

### `-c`

Supply extra configuration, not parameters.

## Custom configuration

### Per-process resources

Resource requests are set per label. To change one process:

```groovy
process {
    withName: 'MCAAT_RUN' {
        memory = 200.GB
        time   = 48.h
    }
}
```

Use this mechanism for steps failing at exit 137 or 140. See the nf-core guide on
[tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources).

### Running in the background

```bash
nextflow run RNABioInfo/mcaat-nf -profile singularity --input s.csv --outdir results -bg
```

Alternatives: run under `screen` or `tmux`, or submit the Nextflow process itself as a job.

### Nextflow memory

Cap Nextflow's own JVM:

```bash
export NXF_OPTS='-Xms1g -Xmx4g'
```
