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
  - [Phage / protospacer curation (experimental)](#phage--protospacer-curation-experimental)
  - [Reporting](#reporting)
  - [Parameters that deliberately do not exist](#parameters-that-deliberately-do-not-exist)
- [Graph re-entry](#graph-re-entry)
- [Test profiles and test data](#test-profiles-and-test-data)
- [Resources and tuning for large metagenomes](#resources-and-tuning-for-large-metagenomes)
- [Resume behaviour and caching](#resume-behaviour-and-caching)
- [Exit codes and troubleshooting](#exit-codes-and-troubleshooting)
- [Core Nextflow arguments](#core-nextflow-arguments)
- [Custom configuration](#custom-configuration)

## Introduction

This pipeline runs [MCAAT](https://github.com/RNABioInfo/mcaat) over quality-controlled
metagenomic reads to find CRISPR repeat–spacer arrays without assembling anything, then parses and
aggregates the results across a cohort.

Everything below assumes you have Nextflow `>=25.10.4` and one container engine (Docker, Singularity
/ Apptainer, Podman or Charliecloud). A `conda` profile is inherited from the nf-core template, but
**it cannot run the detection step**: MCAAT has no Bioconda package, so its processes declare no
Conda environment. Use a container engine.

## Samplesheet input

You need a comma-separated file with a header row. Pass it with `--input`.

```bash
--input '[full path to samplesheet file]'
```

### Columns

| Column    | Required | Description                                                                                                                                          |
| --------- | -------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | yes      | Sample identifier. No spaces. Rows sharing a value are concatenated before analysis. This string is used verbatim in every output path and filename.  |
| `run`     | no       | Sequencing run or lane identifier. No spaces. Only used to keep multiple rows of the same sample distinct; it is dropped once the rows have been merged. |
| `fastq_1` | yes\*    | Full path to the R1 FASTQ (or the single-end FASTQ). Must end `.fastq.gz`, `.fq.gz`, `.fastq` or `.fq`.                                               |
| `fastq_2` | no       | Full path to the R2 FASTQ. Leaving it empty marks the sample as single-end.                                                                           |
| `sdbg`    | no       | Full path to a pre-built MCAAT graph directory, for re-analysis without rebuilding it. Mutually exclusive with `fastq_1`. See [Graph re-entry](#graph-re-entry). |

\* `fastq_1` is required unless the row supplies `sdbg` instead.

### Multiple runs of the same sample

Give the rows the same `sample` value and different `run` values. R1 files are concatenated with
each other and R2 files with each other, in samplesheet order, before any QC step:

```csv
sample,run,fastq_1,fastq_2,sdbg
CONTROL_REP1,L001,AEG588A1_S1_L001_R1_001.fastq.gz,AEG588A1_S1_L001_R2_001.fastq.gz,
CONTROL_REP1,L002,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,
CONTROL_REP1,L003,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,
```

All rows for one sample must agree on single- vs paired-endedness. Mixing them is an error and the
pipeline stops at launch, naming the offending sample.

### Full example

```csv
sample,run,fastq_1,fastq_2,sdbg
gut_A,L001,/data/gut_A_L001_R1.fastq.gz,/data/gut_A_L001_R2.fastq.gz,
gut_A,L002,/data/gut_A_L002_R1.fastq.gz,/data/gut_A_L002_R2.fastq.gz,
soil_B,,/data/soil_B.fastq.gz,,
rerun_C,,,,/archive/mcaat_run1/arrays/rerun_C/graph
```

### Why paired input is handled so strictly

MCAAT accepts exactly one or two input files, and it treats two files as a **strictly synchronised**
paired-end library — record _n_ of file 1 is assumed to be the mate of record _n_ of file 2. If the
two files hold different numbers of records, the underlying graph builder aborts with an empty error
log, which is close to undebuggable.

The pipeline therefore:

- runs `fastp` in strict paired mode and never routes failed/orphan reads onward
  (`--fastp_save_trimmed_fail` writes them to a separate published file only);
- never produces merged reads;
- subsamples both mates with the same `seqtk` seed;
- and asserts, from `seqkit stats`, that R1 and R2 record counts are identical **before** MCAAT is
  dispatched. Samples that fail this check are excluded from the run with an explicit message and
  the rest of the cohort continues.

You can disable that assertion with `--skip_pairsync_check`.

## Running the pipeline

```bash
nextflow run RNABioInfo/mcaat-nf \
   --input ./samplesheet.csv \
   --outdir ./results \
   -profile docker
```

This writes the results to `./results`, plus `work/` (intermediates) and `.nextflow*` (logs and
cache) in the launch directory.

To supply parameters from a file instead of the command line:

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
> Do not put pipeline parameters in a `-c` config file. Use the CLI or `-params-file`. Parameters
> set in `-c` files are not validated and interact badly with profile precedence.

### Updating and reproducibility

```bash
nextflow pull RNABioInfo/mcaat-nf                    # update to the latest version
nextflow run RNABioInfo/mcaat-nf -r 1.0.0 ...        # pin an exact release
```

Always record the `-r` you used. The pipeline writes the full resolved parameter set to
`results/pipeline_info/params_<suffix>.json` on every run.

## Parameters

Every parameter below is validated against `nextflow_schema.json` at launch. Run
`nextflow run RNABioInfo/mcaat-nf --help` for the short form, `--help_full` for everything.

Every parameter listed here reaches a tool. If you find one that does not, that is a bug — please
report it.

### Input/output

| Parameter       | Type   | Default | Description                                |
| --------------- | ------ | ------- | ------------------------------------------ |
| `input`         | path   | —       | **Required.** Samplesheet CSV.             |
| `outdir`        | path   | —       | **Required.** Results directory.           |
| `email`         | string | `null`  | Address to send the completion summary to. |
| `multiqc_title` | string | `null`  | Title used for the MultiQC report.         |

### Read QC

All QC steps are on by default, so they are switched off with `skip_*` flags.

| Parameter                 | Type              | Default                                 | Description                                                                                                                                                                     |
| ------------------------- | ----------------- | --------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `skip_qc`                 | boolean           | `false`                                 | Bypass the whole QC subworkflow and feed samplesheet reads straight to MCAAT. **Unsafe for paired input** — the pair-sync gate lives inside the QC subworkflow and will not run. |
| `skip_fastqc`             | boolean           | `false`                                 | Skip both FastQC instances.                                                                                                                                                     |
| `skip_fastp`              | boolean           | `false`                                 | Skip adapter/quality trimming.                                                                                                                                                  |
| `fastp_save_trimmed_fail` | boolean           | `false`                                 | Write reads that failed fastp filters to a separate published file. They are never fed to MCAAT.                                                                                |
| `save_trimmed_reads`      | boolean           | `false`                                 | Publish the trimmed FASTQ files.                                                                                                                                                |
| `host_fasta`              | path              | `null`                                  | Host genome FASTA (local path or `https://` URL). Setting it enables host depletion with Bowtie2.                                                                               |
| `host_bowtie2_index`      | path              | `null`                                  | Pre-built Bowtie2 index directory. Takes precedence over `host_fasta` and skips the index build.                                                                                |
| `save_hostremoved_reads`  | boolean           | `false`                                 | Publish host-depleted FASTQ files.                                                                                                                                              |
| `remove_phix`             | boolean           | `false`                                 | Enable BBDuk k-mer scrubbing against `phix_reference`.                                                                                                                          |
| `phix_reference`          | path              | `${projectDir}/assets/phix174.fasta.gz` | Contaminant FASTA for BBDuk. The phiX174 genome **is shipped** in `assets/`, so `--remove_phix` works with no extra download. Must stay small — BBDuk builds an in-memory k-mer table, so **never** point this at a mammalian genome. |
| `subsample_reads`         | integer or number | `null`                                  | Read budget per mate. An integer takes that many reads; a float `0 < x < 1` takes that fraction. `null` disables subsampling.                                                    |
| `subsample_seed`          | integer           | `100`                                   | `seqtk sample` seed. Applied identically to both mates, which is what preserves pairing.                                                                                        |
| `skip_pairsync_check`     | boolean           | `false`                                 | Disable the R1/R2 record-count assertion. Do not set this.                                                                                                                      |

### MCAAT array detection

| Parameter                      | Type    | Default   | Description                                                                                                                                                                                                                         |
| ------------------------------ | ------- | --------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `mcaat_min_spacer_len`         | integer | `23`      | Minimum spacer length retained during post-processing.                                                                                                                                                                              |
| `mcaat_max_spacer_len`         | integer | `50`      | Maximum spacer length.                                                                                                                                                                                                              |
| `mcaat_min_repeat_len`         | integer | `23`      | Minimum repeat length. **Cannot be lower than 23**: the de Bruijn graph is built with a hard-coded _k_ of 23, so shorter repeats are not representable.                                                                              |
| `mcaat_max_repeat_len`         | integer | `55`      | Maximum repeat length.                                                                                                                                                                                                              |
| `mcaat_cycle_min_length`       | integer | `27`      | Minimum multicycle length searched in the graph.                                                                                                                                                                                    |
| `mcaat_cycle_max_length`       | integer | `77`      | Maximum multicycle length searched.                                                                                                                                                                                                 |
| `mcaat_threshold_multiplicity` | integer | `20`      | Minimum edge multiplicity for a node to seed a cycle search. Lower it to catch rarer array-bearing organisms, at a steep runtime cost.                                                                                               |
| `mcaat_low_abundance`          | boolean | `true`    | Low-abundance mode.                                                                                                                                                                                                                 |
| `mcaat_keep_graph`             | boolean | `true`    | Keep the graph in the task work directory after the run. Required by the phage stage and by `save_mcaat_graph`. **This is the one documented cache-busting switch** — see [Resume behaviour](#resume-behaviour-and-caching).          |
| `prune_graph`                  | boolean | `true`    | Delete graph-builder intermediates that MCAAT never reads back (`outfile_prefix.bin`, `outfile_prefix.lib_info`, `data.lib`, `graph.mercy_cand.*`). Set `false` only if you need a directory that other MEGAHIT-derived tools can consume; it is roughly 3–4× larger. |
| `save_mcaat_graph`             | boolean | `false`   | Publish the (pruned) graph under `arrays/<sample>/graph/`. Requires `mcaat_keep_graph`.                                                                                                                                              |
| `graph_publish_mode`           | string  | `symlink` | `publishDir` mode used for the graph only, so publishing does not copy tens of gigabytes. One of `symlink`, `copy`, `link`, `rellink`.                                                                                               |
| `save_mcaat_cycles`            | boolean | `false`   | Publish `cycles.txt`. It can be enormous and its line order is not deterministic.                                                                                                                                                   |
| `mcaat_scratch`                | string  | `null`    | If set (e.g. `'$TMPDIR'` or an absolute path), used as the Nextflow `scratch` directive for all MCAAT processes, relocating every intermediate to node-local storage.                                                                |
| `mcaat_args`                   | string  | `null`    | Raw arguments appended to the `mcaat` command line. **Hazard:** MCAAT silently ignores unrecognised arguments, so a typo here has no visible effect. The only thing that catches it is the `parameters.json` readback described below. |

#### The `parameters.json` readback

MCAAT's argument parser has no unknown-argument branch: a misspelled flag is discarded without a
warning and the run proceeds with defaults. To make that detectable, the array-parsing step reads
MCAAT's own `parameters.json` back and asserts that the resolved thread count, multiplicity
threshold, cycle bounds, spacer/repeat bounds and autoclean setting match what the pipeline meant to
set. A mismatch fails the task and prints expected vs. observed. If you use `mcaat_args`, this is
your safety net; check the task log rather than assuming the flag took effect.

The same step writes the values it read into `<sample>.provenance.tsv`, which the cohort layer
consumes so that the resolved MCAAT settings end up as columns of `cohort_summary.tsv`.

### Stage activation

| Parameter            | Type    | Default | Description                                                                        |
| -------------------- | ------- | ------- | ---------------------------------------------------------------------------------- |
| `run_phage_curation` | boolean | `false` | Enable protospacer / phage curation. **Experimental.** Requires `mcaat_keep_graph`. |

This is the only stage toggle.

### Phage / protospacer curation (experimental)

| Parameter        | Type    | Default | Description                             |
| ---------------- | ------- | ------- | --------------------------------------- |
| `save_phage_log` | boolean | `false` | Publish the gzipped `phage-curate` log. |

There is nothing else to tune: every other knob in this MCAAT subcommand is hard-coded (search
depth, and the 500–2500 bp contig length bounds) or accepted-but-unread. The stage is compute-bound
and can be pathologically slow on a dense graph, so it is configured to retry twice and then be
ignored, rather than fail the whole cohort.

Treat its output as a lead rather than a result: the paths it reconstructs are graph walks around a
spacer match, with no independent evidence that the walk corresponds to a real contiguous element.

### Reporting

| Parameter                      | Type    | Default | Description                                                                         |
| ------------------------------ | ------- | ------- | ----------------------------------------------------------------------------------- |
| `skip_cohort_aggregation`      | boolean | `false` | Skip the cohort concatenation, deduplication and summary steps.                     |
| `skip_multiqc`                 | boolean | `false` | Skip the MultiQC report.                                                            |
| `cohort_min_spacer_sharing`    | integer | `2`     | Minimum number of shared spacers for a sample pair to appear in the sharing matrix. |
| `cohort_min_spacer_redundancy` | integer | `2`     | Minimum number of samples a spacer sequence must occur in to appear in `cohort_spacer_redundancy.tsv`. |
| `multiqc_config`               | path    | `null`  | Extra MultiQC config, merged with the pipeline's own.                               |
| `multiqc_logo`                 | path    | `null`  | Custom logo for the report.                                                         |
| `multiqc_methods_description`  | path    | `null`  | Custom methods-description YAML.                                                    |

Both cohort thresholds are rendered into the `COHORT_SUMMARISE` command line as
`--min-spacer-sharing` / `--min-spacer-redundancy`. `ext.args` is appended after them, so a config
can still override either value without the module knowing about it.

### Settings that are not configurable

| Setting                                    | Why                                                                                                                                                                                                                                                        |
| ------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| k-mer size                                 | Hard-coded to 23 in the graph build. This is also why `mcaat_min_repeat_len` cannot go below 23.                                                                                                                                                            |
| minimum k-mer frequency                    | Hard-coded to 1, i.e. no error-k-mer filtering. This is why resource use tracks total sequenced bases.                                                                                                                                                      |
| `--ram`                                    | The flag exists but has no effect (a gigabyte figure is passed into a byte-valued field). The pipeline derives a value internally purely to keep MCAAT's own validation from tripping; `memory` is an OOM ceiling, not a cap MCAAT respects.                 |
| `--threads`                                | Derived from `task.cpus` and clamped to `nproc`.                                                                                                                                                                                                           |
| `--autoclean`                              | Always `false`. Cleanup is done by the pipeline in shell, so that the command text does not change when you flip a stage toggle.                                                                                                                            |
| `--output-folder`                          | Always `mcaat`, inside the task directory.                                                                                                                                                                                                                 |
| `--settings`                               | A settings file's `input-files` key concatenates with the command-line one instead of replacing it, which can silently turn a single-end run into a mis-paired paired-end run.                                                                              |
| `--benchmark`                              | Accepted and stored, but its consumer is never called.                                                                                                                                                                                                     |
| phage `--beam-width`, depth, length bounds | Accepted but never read, or hard-coded.                                                                                                                                                                                                                    |

## Graph re-entry

Building the de Bruijn graph is by far the most expensive step. If you already have a graph from a
previous published run, you can re-analyse it with the phage stage without rebuilding:

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

- The row must supply `sdbg` and must **not** supply `fastq_1`.
- The directory must contain `graph.sdbg_info` and at least one `graph.sdbg.<N>` shard. MCAAT
  loads the graph without validating the path, so the pipeline checks it for you and fails with a
  clear message if the file is missing or empty.
- The stage also needs that sample's arrays. They are located by convention, as the sibling of the
  graph directory: `<sdbg>/../CRISPR_Arrays_*.txt`. That is exactly the layout a previous run
  published under `arrays/<sample>/`, so pointing `sdbg` at `…/arrays/<sample>/graph` works
  unchanged. If no array files are found there, the pipeline stops at launch and names the sample.
- Graph re-entry is only meaningful with `--run_phage_curation`; it does not re-run array
  detection, and a re-entry row therefore contributes no rows to the cohort tables (it carries no
  `parameters.json` and no parsed array table of its own).

To produce re-usable graphs in the first place, run with `--save_mcaat_graph --graph_publish_mode copy`.
The default publish mode for the graph is `symlink`, which points back into `work/` and therefore
does **not** survive `nextflow clean`.

## Test profiles and test data

> [!IMPORTANT]
> **The test-data repository does not exist yet, so `-profile test` cannot run today.** This is
> stated plainly rather than papered over.

The four test profiles are `test`, `test_phage`, `test_graph` and `test_full`. Every one of them
builds its `--input` from `params.pipelines_testdata_base_path`, whose default is:

```
https://raw.githubusercontent.com/RNABioInfo/mcaat-nf-testdata/main/
```

**That repository has to be created and populated before any of the test profiles will run.**
Until it exists, launching `-profile test` fails while resolving the samplesheet URL — the failure
is immediate and has nothing to do with MCAAT.

The layout the profiles expect, relative to that base path:

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

The reads for `-profile test` should be small enough for CI: a few thousand read pairs spiked with
at least one real CRISPR array, so that the run exercises the parser and the cohort layer rather
than only the zero-result path. At least two samples are needed for the spacer-sharing matrix and
the MultiQC heatmap to be produced at all — with a single sample they are deliberately skipped.

Until the repository exists, you can still exercise the profiles against local files:

```bash
nextflow run RNABioInfo/mcaat-nf -profile test,docker \
   --pipelines_testdata_base_path /path/to/local/testdata/ \
   --outdir results_test
```

The trailing slash is required: the profiles concatenate the base path with
`mcaat-nf/samplesheet/…` directly.

## Resources and tuning for large metagenomes

### The one lever that actually works

MCAAT builds its graph with no minimum k-mer frequency filter, so every sequencing error contributes
edges. Graph size — and therefore RAM, disk and runtime — scales with **total sequenced bases**, not
with community complexity. Deeper sequencing of the same community costs proportionally more.

`--subsample_reads` is consequently the pipeline's principal cost control:

```bash
--subsample_reads 20000000     # 20 M reads per mate
--subsample_reads 0.25         # 25 % of reads
```

Halving depth roughly halves the graph. Do this before reaching for larger nodes.

### Default resource labels

| Process             | cpus | memory (attempt 1) | time (attempt 1) |
| ------------------- | ---- | ------------------ | ---------------- |
| `MCAAT_RUN`         | 8    | 72 GB              | 16 h             |
| `MCAAT_PHAGECURATE` | 4    | 48 GB              | 24 h             |

`memory` and `time` double on each retry; **`cpus` never changes**. That is deliberate: the cycle
finder allocates one visited-bitset per thread, sized from the number of graph edges, so adding
cores adds memory linearly. Escalating cores after an out-of-memory kill makes the next attempt
worse, not better. Fixing `cpus` also pins the number of graph shard files, which keeps `-resume`
hashes stable across profiles.

### Calibrating memory

Because `--ram` has no effect inside MCAAT, the `memory` directive is only the point at which the
scheduler or cgroup kills the job. It has to be measured rather than derived. Peak resident memory
is approximately the larger of the two phases (they do not overlap):

```
build   ~  total_bases / 4  +  8 * n_reads          bytes
detect  ~  2 * n_edges  +  cpus * n_edges / 8       bytes
```

Practical procedure:

1. Run ten representative samples with `-with-report -with-trace`.
2. For each, `grep 'Memory estimate' work/<hash>/.command.out` — MCAAT prints its own visited-bitset
   figure.
3. Read `peak_rss` from `results/pipeline_info/execution_trace_*.txt`.
4. Set `memory` to the observed maximum plus about 30 %.

If the trace columns are blank, your container image lacks `ps` (the `procps` package). The
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

`MCAAT_RUN` requests `50 GB + 4 × (input GB)`, scaled by attempt. The multiplier covers the 2-bit
re-encoding of every input read that the graph builder writes (transient, but large), the graph
shards themselves, `cycles.txt`, and the retained pruned graph.

Rough retained-graph sizing: about **2 GB of pruned graph per gigabase of input**. For a 100-sample
cohort at 10 Gbp each that is around 2 TB sitting in `work/`. If that is your binding constraint,
set `--mcaat_keep_graph false` for the campaign — but read the caching note below first, and be
aware it makes the phage stage unavailable.

### Node-local scratch

If your cluster provides fast node-local storage, set `--mcaat_scratch '$TMPDIR'`. MCAAT writes
everything relative to its working directory, so relocating the task directory relocates all the
intermediates. It is off by default because not every site provisions such storage.

### Heterogeneous clusters

Use the pipeline's own container image (`quay.io/rnabioinfo/mcaat:1.0.0-portable`). The upstream
`docker.io/feeka94/mcaat:1.0.0` image is compiled with `-march=native`, so it can die with
`Illegal instruction` (exit 132) on nodes older than the machine that built it. Exit 132 is
deliberately never retried, because retrying onto the same node class fails identically.

## Resume behaviour and caching

`-resume` works normally, with a few specifics worth knowing.

**Stage toggles do not bust the array-detection cache.** The `MCAAT_RUN` command text deliberately
contains no reference to `run_phage_curation`, `save_mcaat_graph` or `save_mcaat_cycles`. Those only
steer publishing. So this sequence costs one graph build, not two:

```bash
# month 1
nextflow run RNABioInfo/mcaat-nf --input s.csv --outdir results -profile singularity

# month 6 — resumes every MCAAT_RUN from cache, pays only for the phage chain
nextflow run RNABioInfo/mcaat-nf --input s.csv --outdir results -profile singularity \
    --run_phage_curation -resume
```

**`--mcaat_keep_graph` does bust it.** It is the single toggle that changes the resolved command
text, because graph retention is done in shell inside the task. Decide it once per campaign and do
not change it mid-study. It is separated from the stage toggle precisely so that everything else
stays cache-invariant.

**Changing `cpus` busts it too**, and additionally changes the number of graph shard files. Do not
retune `process_mcaat_run` cpus partway through a cohort.

**Other things that change the hash:** any MCAAT parameter in the `mcaat_options` group, `mcaat_args`,
the input files themselves, and the container image tag or digest.

**`work/` is not disposable while you rely on published graphs.** With the default
`--graph_publish_mode symlink`, `arrays/<sample>/graph/` is a symlink into `work/`. `nextflow clean`
will leave you with dangling links. Use `--graph_publish_mode copy` if the published graph needs to
outlive the work directory.

## Exit codes and troubleshooting

The pipeline remaps some of MCAAT's exit statuses so that deterministic failures are not retried
three times before being reported.

| Code                             | Meaning                                                | Retried? | What to do                                                                                                              |
| -------------------------------- | ------------------------------------------------------ | -------- | ------------------------------------------------------------------------------------------------------------------------- |
| `64`                             | Deterministic MCAAT argument or validation error       | no       | Read `mcaat.log` (published as `arrays/<sample>/<sample>.mcaat.log`) — the real message is there, not in `.command.err`. |
| `65`                             | Graph build finished but produced no `graph.sdbg_info` | no       | Check FASTQ integrity and R1/R2 record parity. The underlying builder's stderr is suppressed, so `.command.err` will be empty by design. |
| `137`                            | Allocation failure or OOM kill                         | yes      | Memory doubles on retry. If it still fails, subsample the reads.                                                        |
| `132`                            | Illegal instruction                                    | no       | The binary was built for a newer CPU than the node. Use the portable image, or pin the queue.                           |
| `134`                            | Abort that matched no known pattern                    | no       | Unknown deterministic failure — open an issue with `mcaat.log` attached.                                                |
| `1`                              | Graph-builder fatal error                              | no       | Deterministic. `.command.err` will be empty; inspect the inputs.                                                        |
| `104`, `140`, `143`, `175`–`177` | Scheduler / walltime kills                             | yes      | Time doubles on retry.                                                                                                  |

**"My sample disappeared from the results."** Look for a `MCAAT-NF: EXCLUDING sample` line in the
Nextflow log: the pair-sync gate dropped it because R1 and R2 had different record counts.

**"A sample reports zero arrays."** That is a legitimate result and is reported explicitly — the
sample gets a header-only `CRISPR_Arrays_1.txt`, a zero row in `reports/cohort_summary.tsv`, and a
`0` in the MultiQC General Statistics table. It is never silently dropped.

**"`-profile test` fails immediately."** The test-data repository does not exist yet; see
[Test profiles and test data](#test-profiles-and-test-data).

**"Where is the real error message?"** For the array-detection step, in `mcaat.log` / `.command.out`.
MCAAT redirects the graph builder's stderr to `/dev/null`, so an empty `.command.err` is expected
even for genuine failures.

## Core Nextflow arguments

> [!NOTE]
> These are Nextflow's own options and use a _single_ hyphen.

### `-profile`

Use one of the container profiles. The template's `conda` profile exists but cannot run the
detection step, because MCAAT has no Bioconda package.

- `docker`
- `singularity`
- `podman`
- `apptainer`
- `charliecloud`
- `wave`

Test profiles: `test`, `test_phage`, `test_graph`, `test_full`. None of them can run until the
test-data repository exists — see [above](#test-profiles-and-test-data).

Profiles are applied left to right, so put the container profile last:
`-profile test,docker`.

### `-resume`

Restart from cached results. See [Resume behaviour](#resume-behaviour-and-caching). You can also
resume a specific earlier run by name or session ID (`nextflow log` lists them).

### `-c`

Supply extra _configuration_ (not parameters).

## Custom configuration

### Per-process resources

Requests are per-process, using labels. To change one process rather than a whole label:

```groovy
process {
    withName: 'MCAAT_RUN' {
        memory = 200.GB
        time   = 48.h
    }
}
```

If a step is failing at exit 137 or 140, that is the mechanism to reach for. See the nf-core guide
on [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources).

### Running in the background

```bash
nextflow run RNABioInfo/mcaat-nf -profile singularity --input s.csv --outdir results -bg
```

Or run it under `screen`/`tmux`, or submit the Nextflow process itself as a job.

### Nextflow memory

Nextflow's own JVM can be capped:

```bash
export NXF_OPTS='-Xms1g -Xmx4g'
```
