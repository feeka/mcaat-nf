# RNABioInfo/mcaat-nf: Output

## Introduction

This document describes every file the pipeline publishes into `--outdir`. Nothing is listed here
that the pipeline does not actually produce.

The cross-sample products in [`reports/`](#reports) are the primary result: they are what this
pipeline produces that running MCAAT sample-by-sample does not. Per-sample outputs live under
`arrays/` and (optionally) `phage/`, and the raw MCAAT reports are always published alongside the
parsed tables so you are never dependent on the parser having got it right.

## Directory layout

```
<outdir>/
├── qc/
│   ├── fastqc/raw/<sample>/
│   ├── fastqc/trimmed/<sample>/
│   ├── fastp/<sample>/
│   ├── bowtie2/<sample>/
│   ├── bbduk/<sample>/
│   └── seqkit/<sample>/
├── arrays/<sample>/
├── phage/<sample>/               (only with --run_phage_curation)
├── reports/
├── multiqc/
└── pipeline_info/
```

## Reports

<details markdown="1">
<summary>Output files</summary>

- `reports/`
  - `cohort_arrays.tsv` — every spacer from every sample, one row each. Written by the cohort
    summariser, which validates and sanitises each row as it streams it.
  - `cohort_arrays_concat.tsv` — the same rows, produced independently by `csvtk concat` straight
    from the per-sample tables. Kept as an unprocessed cross-check; it carries no extra columns.
  - `cohort_summary.tsv` — one row per sample.
  - `cohort_spacer_sharing.tsv` — one row per sample pair that shares spacers.
  - `cohort_spacer_redundancy.tsv` — one row per spacer sequence seen in several samples.
  - `cohort_repeat_families.tsv` — one row per distinct consensus repeat.
  - `cohort_spacers.unique.fasta` — the pooled spacer set, deduplicated by sequence.

</details>

### `cohort_arrays.tsv`

Every per-sample `<sample>.arrays.tsv` concatenated under a single header, with the column set
asserted row by row. Column meanings are identical to the per-sample file, documented
[below](#samplearraystsv).

`cohort_arrays_concat.tsv` holds the same content from a plain `csvtk concat` of the same inputs.
The two files exist under distinct names on purpose: they used to collide on the publish path, and
overwriting one with the other silently would have made it impossible to tell which tool wrote the
file you were reading. If they ever differ, the summariser's version is the one that has been
validated.

### `cohort_summary.tsv`

One row per sample. Samples that legitimately found nothing appear here with zeros; they are never
omitted. The sample roster is built from the staged filenames as well as the data rows, so a
zero-array sample that contributed no rows at all still gets a row here.

The final three columns are provenance echoed from MCAAT's own `parameters.json`, by way of the
per-sample `<sample>.provenance.tsv`.

| Column                         | Description                                                                      |
| ------------------------------ | -------------------------------------------------------------------------------- |
| `sample`                       | Sample identifier.                                                               |
| `n_arrays`                     | Number of distinct arrays.                                                       |
| `n_spacers`                    | Total spacer occurrences.                                                        |
| `n_spacers_unique`             | Distinct spacer sequences within the sample.                                     |
| `n_consensus_repeats`          | Distinct consensus repeats.                                                      |
| `median_repeat_len`            | Median consensus repeat length.                                                  |
| `median_spacer_len`            | Median spacer length.                                                            |
| `min_spacer_len`               | Shortest spacer.                                                                 |
| `max_spacer_len`               | Longest spacer.                                                                  |
| `n_spacers_shared`             | Unique spacers also seen in at least one other sample.                           |
| `pct_spacers_shared`           | `100 × n_spacers_shared / n_spacers_unique`, or `0.0` when there are no spacers.  |
| `n_repeats_shared`             | Consensus repeats also seen in at least one other sample.                        |
| `mcaat_threads`                | Thread count MCAAT actually resolved (from `parameters.json`).                   |
| `mcaat_ram_gb`                 | RAM figure MCAAT recorded. Provenance only — it does not bound anything.         |
| `mcaat_threshold_multiplicity` | Multiplicity threshold MCAAT actually resolved.                                  |

```tsv
sample	n_arrays	n_spacers	n_spacers_unique	n_consensus_repeats	median_repeat_len	median_spacer_len	min_spacer_len	max_spacer_len	n_spacers_shared	pct_spacers_shared	n_repeats_shared	mcaat_threads	mcaat_ram_gb	mcaat_threshold_multiplicity
gut_A	12	97	95	11	31.0	34.0	26	48	7	7.37	2	8	64	20
soil_B	4	22	22	4	29.0	32.0	28	41	7	31.82	2	8	64	20
empty_C	0	0	0	0	0.0	0.0	0	0	0	0.0	0	8	64	20
```

### `cohort_spacer_sharing.tsv`

Long format, one row per unordered sample pair with `sample_a < sample_b`. Only pairs sharing at
least `--cohort_min_spacer_sharing` spacers are written.

| Column             | Description                                    |
| ------------------ | ---------------------------------------------- |
| `sample_a`         | First sample of the pair.                      |
| `sample_b`         | Second sample of the pair.                     |
| `n_shared_spacers` | Exact-sequence spacer matches.                 |
| `jaccard`          | Jaccard index over the two unique spacer sets. |
| `n_shared_repeats` | Shared consensus repeats.                      |

```tsv
sample_a	sample_b	n_shared_spacers	jaccard	n_shared_repeats
gut_A	soil_B	7	0.0637	2
```

Matches are exact sequence identity, in the orientation MCAAT emitted. There is no reverse-complement
folding and no clustering, so this is a conservative lower bound on sharing.

### `cohort_spacer_redundancy.tsv`

The same sharing information organised by spacer rather than by sample pair: one row per spacer
sequence that occurs in at least `--cohort_min_spacer_redundancy` samples. This is where to look for
individual spacers that recur across the cohort — usually the interesting ones, because a spacer
shared between communities implies a shared mobile-element pressure.

| Column                | Description                                                                         |
| --------------------- | ----------------------------------------------------------------------------------- |
| `spacer`              | The spacer sequence.                                                                |
| `spacer_len`          | Its length.                                                                         |
| `n_samples`           | Number of samples carrying it.                                                      |
| `samples`             | Semicolon-separated sample identifiers.                                             |
| `n_occurrences`       | Total occurrences across the cohort (a spacer can occur more than once per sample). |
| `n_consensus_repeats` | Number of distinct consensus repeats it was found next to.                          |
| `consensus_repeats`   | Semicolon-separated list of those repeats.                                          |

```tsv
spacer	spacer_len	n_samples	samples	n_occurrences	n_consensus_repeats	consensus_repeats
AGGCTTAAGCCTGAGGCTAAGCTAGTTCAGGCTA	34	2	gut_A;soil_B	3	1	GTTTCAATCCACGCACCCGTGGGGGTGATCC
TTGACCTGAAGCTTAGCCTGAAGCTTGACCTGAA	34	2	gut_A;soil_B	2	2	GTTTCAATCCACGCACCCGTGGGGGTGATCC;CGGTTTATCCCCGCTGGCGCGGGGAACAC
```

### `cohort_repeat_families.tsv`

One row per distinct consensus repeat across the cohort. This is the table to start from when asking
which CRISPR systems recur across your samples.

| Column              | Description                                   |
| ------------------- | --------------------------------------------- |
| `consensus_repeat`  | The repeat sequence.                          |
| `repeat_len`        | Its length.                                   |
| `n_samples`         | Number of samples carrying it.                |
| `samples`           | Semicolon-separated sample identifiers.       |
| `n_arrays_total`    | Arrays with this consensus across the cohort. |
| `n_spacers_total`   | Spacer occurrences belonging to those arrays. |
| `n_spacers_unique`  | Distinct spacer sequences within the family.  |
| `median_spacer_len` | Median spacer length within the family.       |

```tsv
consensus_repeat	repeat_len	n_samples	samples	n_arrays_total	n_spacers_total	n_spacers_unique	median_spacer_len
GTTTCAATCCACGCACCCGTGGGGGTGATCC	31	2	gut_A;soil_B	3	26	24	34.0
CGGTTTATCCCCGCTGGCGCGGGGAACAC	29	1	gut_A	1	2	2	32.0
```

### `cohort_spacers.unique.fasta`

The pooled per-sample spacer FASTAs, deduplicated on sequence with `seqkit rmdup -s`. Each retained
record keeps the header of its first occurrence, so one representative sample and array remain
visible; the full sample list for a spacer is in `cohort_spacer_redundancy.tsv`. Record order
follows the concatenation order of the per-sample files, which is not guaranteed stable between
runs — do not checksum this file.

This is the file to hand to BLAST, `bowtie2` or a viral-genome search if you want to identify
protospacer matches against known mobile elements.

## Arrays (per sample)

<details markdown="1">
<summary>Output files</summary>

- `arrays/<sample>/`
  - `CRISPR_Arrays_*.txt` — MCAAT's own human-readable report, unmodified.
  - `<sample>.arrays.tsv` — the parsed, tidy version.
  - `<sample>.spacers.fasta` — every spacer as FASTA.
  - `parameters.json` — MCAAT's record of the parameters it actually resolved.
  - `<sample>.mcaat.log` — MCAAT's stdout, including the progress banners.
  - `graph/` — the succinct de Bruijn graph (only with `--save_mcaat_graph`).
  - `cycles.txt` — raw multicycles (only with `--save_mcaat_cycles`).

</details>

### `CRISPR_Arrays_*.txt`

MCAAT's native output, published verbatim. Large runs are split across numbered files; note that
`# Arrays` and `# Spacers` in the header are **run-wide totals repeated in every split file**, so do
not sum them.

```text
# MCAAT — CRISPR Array Output
# Generated : 2026-07-21 11:04:33
# Arrays    : 12
# Spacers   : 97

>Array_7  spacers=3
GTTTCAATCCACGCACCCGTGGGGGTGATCC
        -------------------------------	AGGCTTAAGCCTGAGGCTAAGCTAGTTCAGGCTA
        -------------------------------	TTGACCTGAAGCTTAGCCTGAAGCTTGACCTGAA
        -------------------------------	GGAATTCCAATCGATCGAATACCCACTTGACCTG

>Array_8  spacers=2
CGGTTTATCCCCGCTGGCGCGGGGAACAC
        CGGTTTATCCCCGCTGGCGCGGGGAAC	TTGCCAGGTCAATCGCGTTAACCGTTAGCCAT
        CGGTTTATCCCCGCTGGCGCGGGGAAC	ACCTGGTTACGCCAATTGGCACCGTTAACGGA
```

The line after `>Array_N` is the consensus repeat. Each following line is eight spaces, then the
repeat variant observed for that spacer, then a **single tab**, then the spacer. A run of dashes
means "identical to the consensus"; a literal sequence means the observed repeat differed, and it
may be shorter than the consensus.

If a sample finds nothing, this file still exists, containing only the header with `# Arrays : 0`.
That is a real result, not an error.

### `<sample>.arrays.tsv`

One row per spacer. This is the file to load into R or pandas.

| Column                        | Type   | Description                                                                           |
| ----------------------------- | ------ | ------------------------------------------------------------------------------------- |
| `sample`                      | string | Sample identifier.                                                                    |
| `array_id`                    | int    | Array number as MCAAT assigned it, globally within the run.                           |
| `source_file`                 | string | Basename of the `CRISPR_Arrays_*.txt` the row came from.                              |
| `consensus_repeat`            | string | Consensus repeat of the array.                                                        |
| `consensus_repeat_len`        | int    | Its length.                                                                           |
| `n_spacers`                   | int    | Number of spacers in the array, as declared on the `>Array_N` line.                   |
| `spacer_index`                | int    | 1-based position of this spacer within the array, in file order.                      |
| `repeat_variant`              | string | The repeat instance adjacent to this spacer. Dash runs are expanded to the consensus. |
| `repeat_variant_is_consensus` | bool   | `true` when the source column was a dash run.                                         |
| `spacer`                      | string | Spacer sequence.                                                                      |
| `spacer_len`                  | int    | Its length.                                                                           |

```tsv
sample	array_id	source_file	consensus_repeat	consensus_repeat_len	n_spacers	spacer_index	repeat_variant	repeat_variant_is_consensus	spacer	spacer_len
gut_A	7	CRISPR_Arrays_1.txt	GTTTCAATCCACGCACCCGTGGGGGTGATCC	31	3	1	GTTTCAATCCACGCACCCGTGGGGGTGATCC	true	AGGCTTAAGCCTGAGGCTAAGCTAGTTCAGGCTA	34
gut_A	7	CRISPR_Arrays_1.txt	GTTTCAATCCACGCACCCGTGGGGGTGATCC	31	3	2	GTTTCAATCCACGCACCCGTGGGGGTGATCC	true	TTGACCTGAAGCTTAGCCTGAAGCTTGACCTGAA	34
gut_A	7	CRISPR_Arrays_1.txt	GTTTCAATCCACGCACCCGTGGGGGTGATCC	31	3	3	GTTTCAATCCACGCACCCGTGGGGGTGATCC	true	GGAATTCCAATCGATCGAATACCCACTTGACCTG	34
gut_A	8	CRISPR_Arrays_1.txt	CGGTTTATCCCCGCTGGCGCGGGGAACAC	29	2	1	CGGTTTATCCCCGCTGGCGCGGGGAAC	false	TTGCCAGGTCAATCGCGTTAACCGTTAGCCAT	32
gut_A	8	CRISPR_Arrays_1.txt	CGGTTTATCCCCGCTGGCGCGGGGAACAC	29	2	2	CGGTTTATCCCCGCTGGCGCGGGGAAC	false	ACCTGGTTACGCCAATTGGCACCGTTAACGGA	32
```

Notes on interpretation:

- `array_id` is not a coordinate. Arrays are graph structures reconstructed from reads; there is no
  contig, no genome position and no strand.
- The repeat variant is constant within an array block — MCAAT records one repeat variant per array,
  not one per spacer. `repeat_variant_is_consensus` is therefore also block-constant.
- A non-consensus variant is frequently **shorter** than the consensus. The consensus is assembled
  from a tail and a front segment; the observed variant is the untailed repeat.
- A sample with zero arrays produces this file with a header row and no data rows.

The same step also writes `<sample>.provenance.tsv` — one row, `sample` plus the eleven MCAAT
settings read back out of `parameters.json`. It is an internal handoff to the cohort layer rather
than a published deliverable: its content reaches you as the `mcaat_*` columns of
`cohort_summary.tsv`, and the underlying `parameters.json` is published in full next to it.

### `<sample>.spacers.fasta`

Every spacer, unwrapped, one sequence per line. The identifier encodes sample, array and position —
which is what makes the pooled cohort FASTA addressable — and the description carries the length of
**this spacer's repeat variant, not the consensus**, together with the spacer length. The consensus
figure lives in `consensus_repeat_len` in `<sample>.arrays.tsv`; for an array whose observed variant
is a truncation of the consensus the two differ, and `repeat_len=` reports the shorter one.

```text
>gut_A|array7|sp1 repeat_len=31 spacer_len=34
AGGCTTAAGCCTGAGGCTAAGCTAGTTCAGGCTA
>gut_A|array7|sp2 repeat_len=31 spacer_len=34
TTGACCTGAAGCTTAGCCTGAAGCTTGACCTGAA
>gut_A|array7|sp3 repeat_len=31 spacer_len=34
GGAATTCCAATCGATCGAATACCCACTTGACCTG
>gut_A|array8|sp1 repeat_len=27 spacer_len=32
TTGCCAGGTCAATCGCGTTAACCGTTAGCCAT
>gut_A|array8|sp2 repeat_len=27 spacer_len=32
ACCTGGTTACGCCAATTGGCACCGTTAACGGA
```

Array 7 above is the case where variant and consensus coincide, so `repeat_len=31` equals
`consensus_repeat_len`. Array 8 is the case where they do not: its `consensus_repeat_len` is 29 but
its `repeat_variant` is the 27 bp `CGGTTTATCCCCGCTGGCGCGGGGAAC`, and 27 is what the header carries.
Joining `repeat_len=` against `consensus_repeat_len` will therefore mis-join; join on the
identifier, or on `repeat_variant`.

The header is `>{sample}|array{N}|sp{i}` followed by **a space and two `key=value` attributes**.
Anything matching against these headers has to allow for that trailing description — a pattern
anchored immediately after `sp{i}` matches no record at all.

Note that the identifier stops at the first space, so `gut_A|array7|sp1` is what most tools will
carry through. These are the sequences to BLAST or map against viral databases if you want
protospacer matches against known phages, and they are the input to the cohort sharing analysis.

### `parameters.json`

MCAAT's own record of the settings it resolved. Published as provenance. The pipeline reads it back
and fails the sample if the values disagree with what it intended to set, which is the only way to
detect a flag that MCAAT accepted and then ignored.

Two fields in it are misleading and are not used by anything: `output_file` names a
`CRISPR_Arrays.txt` that the release build never writes, and `cycles_folder` names a directory that
is created but never written to.

### `<sample>.mcaat.log`

MCAAT's stdout for the run: the `[1/4] … [4/4]` progress banners, the memory estimate it printed for
its own visited bitsets, and any error message. **Read this file, not `.command.err`.** MCAAT
suppresses the graph builder's stderr entirely, so a failed run typically has an empty `.command.err`
and the actual explanation here.

### `graph/`

Published only with `--save_mcaat_graph`, and by default as a **symlink** into the work directory
(`--graph_publish_mode`). Contents:

- `graph.sdbg_info` — graph metadata.
- `graph.sdbg.<N>` — the graph shards. The number of shards equals the thread count used, so this
  file set changes if you change `cpus`.

With `--prune_graph true` (default) the graph-builder intermediates that MCAAT never reads back have
been removed; the directory is still a valid value for the samplesheet's `sdbg` column
([graph re-entry](usage.md#graph-re-entry)), but it is no longer a complete MEGAHIT working
directory. Set `--prune_graph false` if you need that.

### `cycles.txt`

Published only with `--save_mcaat_cycles`. One line per multicycle found in the graph. It can be very
large, and its line order derives from hash-map iteration, so two runs on identical input can produce
byte-different files. Do not checksum it.

## Phage / protospacer curation (experimental)

<details markdown="1">
<summary>Output files</summary>

- `phage/<sample>/`
  - `<sample>.phage_contigs.fasta` — reconstructed contiguous paths carrying protospacer matches.
  - `phage_curate.log.gz` — the run log (only with `--save_phage_log`).

</details>

Produced only with `--run_phage_curation`, which is off by default and experimental.

Every record is between 501 and 2500 bp by construction; those bounds are hard-coded in MCAAT.

MCAAT emits these records using the reconstructed spacer sequence itself as the FASTA identifier,
which is neither unique nor informative in a cohort, so the pipeline rewrites the headers:

```text
>gut_A|phage1 spacer=AGGCTTAAGCCTGAGGCTAAGCTAGTTCAGGCTA len=1183
TTGACAGCCTGAAGGTCTTAACGCCATTGCGGATCAGCTTTACCGGCAAAGTCGATCAGCT
...
```

`spacer=` is the array spacer whose protospacer this path was reconstructed around, so it joins back
to `<sample>.spacers.fasta` and `cohort_arrays.tsv` on the sequence. `len=` is the contig length.

The stage is compute-bound and unbounded in the pathological case. It retries twice and is then
ignored, so a single difficult sample cannot fail the cohort — if a sample has no
`phage_contigs.fasta`, check the Nextflow log for an ignored task. These contigs are graph walks
around a spacer match and carry no independent support; treat them as leads.

## Read QC

<details markdown="1">
<summary>Output files</summary>

- `qc/fastqc/raw/<sample>/` and `qc/fastqc/trimmed/<sample>/` — `*.html`, `*.zip`.
- `qc/fastp/<sample>/` — `*.json`, `*.html`, `*.log`, and `*.fastp.fastq.gz` with
  `--save_trimmed_reads`, plus failed reads with `--fastp_save_trimmed_fail`.
- `qc/bowtie2/<sample>/` — alignment `*.log`, and host-depleted `*.fastq.gz` with
  `--save_hostremoved_reads`.
- `qc/bbduk/<sample>/` — `*.log` from phiX removal.
- `qc/seqkit/<sample>/` — `*.tsv` read statistics computed on the reads that are handed to MCAAT.

</details>

The `seqkit` table is not only for the report: it is the source of the R1/R2 record-count assertion
that keeps mis-paired libraries away from MCAAT. If a sample was excluded, its `seqkit` TSV is where
you can see the two counts that disagreed.

The "trimmed" FastQC instance runs on the reads _after_ every enabled QC step including subsampling,
so it describes exactly what MCAAT received.

## MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html` — the report.
  - `multiqc_data/` — parsed data behind every plot, as flat files.
  - `multiqc_plots/` — exported static plot images.

</details>

Beyond the standard QC sections, the report contains exactly five MCAAT contributions. They are
emitted as `*_mqc.tsv` / `*_mqc.yaml` files carrying their own configuration, and their order is
fixed by `assets/multiqc_config.yml`:

- **General Statistics** columns `CRISPR arrays` and `Spacers` per sample (id `mcaat_genstats`,
  contributed per sample rather than as a section of its own). A sample that found nothing shows
  `0`, not a blank.
- **MCAAT cohort summary** (`mcaat_cohort_summary`) — the `cohort_summary.tsv` table.
- **MCAAT arrays and spacers** (`mcaat_arrays_bargraph`) — bar graph of arrays and unique spacers
  per sample.
- **MCAAT spacer length distribution** (`mcaat_spacerlen`) — line graph per sample.
- **MCAAT repeat length distribution** (`mcaat_repeatlen`) — line graph per sample, counting each
  array's consensus repeat once.
- **MCAAT spacer sharing** (`mcaat_sharing`) — sample × sample heatmap of exact spacer matches.

The sharing heatmap is skipped, with a note on stderr, when the cohort has fewer than two samples:
MultiQC renders a degenerate single-cell plot otherwise.

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report_<suffix>.html`, `execution_timeline_<suffix>.html`,
    `execution_trace_<suffix>.txt`, `pipeline_dag_<suffix>.html` — Nextflow's own run reports.
  - `mcaat_software_versions.yml` — versions of every tool used.
  - `params_<suffix>.json` — the fully resolved parameter set for the run.

</details>

`execution_trace_<suffix>.txt` carries the `peak_rss` column used by the memory-calibration
procedure in [`usage.md`](usage.md#calibrating-memory). If that column is empty, the container image
in use lacks `ps`.
