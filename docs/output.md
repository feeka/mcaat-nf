# RNABioInfo/mcaat-nf: Output

## Introduction

This document describes every file the pipeline publishes into `--outdir`. Per-sample outputs
live under `arrays/` and, when phage curation is enabled, `phage/`. Cross-sample outputs live
under [`reports/`](#reports). MCAAT's own reports are published alongside the parsed tables.

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

| Path                                   | Format | Contents                                                       |
| -------------------------------------- | ------ | -------------------------------------------------------------- |
| `reports/cohort_arrays.tsv`            | TSV    | Every spacer from every sample, one row each.                  |
| `reports/cohort_arrays_concat.tsv`     | TSV    | The same rows, produced by `csvtk concat` from the per-sample tables. |
| `reports/cohort_summary.tsv`           | TSV    | One row per sample.                                            |
| `reports/cohort_spacer_sharing.tsv`    | TSV    | One row per sample pair that shares spacers.                   |
| `reports/cohort_spacer_redundancy.tsv` | TSV    | One row per spacer sequence seen in several samples.           |
| `reports/cohort_repeat_families.tsv`   | TSV    | One row per distinct consensus repeat.                         |
| `reports/cohort_spacers.unique.fasta`  | FASTA  | The pooled spacer set, deduplicated by sequence.               |

All tables are UTF-8, LF-terminated, with a single header row. Every table is written even when
the cohort contributed no arrays; in that case only the header row is present.

### `cohort_arrays.tsv`

Every per-sample `<sample>.arrays.tsv` concatenated under a single header. The column set is
asserted row by row: a file missing any of the eleven columns is skipped with a warning on
stderr. Column meanings are identical to the per-sample file, documented
[below](#samplearraystsv).

`cohort_arrays_concat.tsv` holds the same content from a plain `csvtk concat` of the same
inputs, without the row-by-row column assertion and without extra columns. The two files carry
distinct names so that neither overwrites the other on the shared publish path. If they differ,
`cohort_arrays.tsv` is the validated version.

### `cohort_summary.tsv`

One row per sample, sorted by `sample`. Samples that found nothing appear with zeros. The sample
roster is built from three sources: the staged per-sample filenames, the `sample` column of the
per-sample `<sample>.provenance.tsv` tables, and the data rows themselves. A sample that
contributed no data rows still gets a row.

The final three columns are provenance echoed from MCAAT's `parameters.json` by way of
`<sample>.provenance.tsv`. They are `NA` when `parameters.json` was missing or unparseable.

| Column                         | Type   | Description                                                                     |
| ------------------------------ | ------ | ------------------------------------------------------------------------------- |
| `sample`                       | string | Sample identifier.                                                              |
| `n_arrays`                     | int    | Number of distinct arrays.                                                      |
| `n_spacers`                    | int    | Total spacer occurrences.                                                       |
| `n_spacers_unique`             | int    | Distinct spacer sequences within the sample.                                    |
| `n_consensus_repeats`          | int    | Distinct consensus repeats.                                                     |
| `median_repeat_len`            | float  | Median consensus repeat length, one decimal place.                              |
| `median_spacer_len`            | float  | Median spacer length, one decimal place.                                        |
| `min_spacer_len`               | int    | Shortest spacer. `0` when the sample has no spacers.                            |
| `max_spacer_len`               | int    | Longest spacer. `0` when the sample has no spacers.                             |
| `n_spacers_shared`             | int    | Unique spacers also seen in at least one other sample.                          |
| `pct_spacers_shared`           | float  | `100 × n_spacers_shared / n_spacers_unique`, two decimal places; `0.00` when there are no spacers. |
| `n_repeats_shared`             | int    | Consensus repeats also seen in at least one other sample.                       |
| `mcaat_threads`                | string | Thread count MCAAT resolved, from `parameters.json`.                            |
| `mcaat_ram_gb`                 | string | The memory ceiling MCAAT used, its own default of 95 % of system RAM.           |
| `mcaat_threshold_multiplicity` | string | Multiplicity threshold MCAAT resolved.                                          |

```tsv
sample	n_arrays	n_spacers	n_spacers_unique	n_consensus_repeats	median_repeat_len	median_spacer_len	min_spacer_len	max_spacer_len	n_spacers_shared	pct_spacers_shared	n_repeats_shared	mcaat_threads	mcaat_ram_gb	mcaat_threshold_multiplicity
gut_A	12	97	95	11	31.0	34.0	26	48	7	7.37	2	8	64	20
soil_B	4	22	22	4	29.0	32.0	28	41	7	31.82	2	8	64	20
empty_C	0	0	0	0	0.0	0.0	0	0	0	0.00	0	8	64	20
```

### `cohort_spacer_sharing.tsv`

Long format, one row per unordered sample pair with `sample_a < sample_b`. Only pairs sharing at
least `--cohort_min_spacer_sharing` spacers (default `2`) are written. Rows are sorted by
`sample_a`, then `sample_b`.

Matches are exact sequence identity in the orientation MCAAT emitted. No reverse-complement
folding and no clustering are applied, so the counts are a lower bound on sharing.

| Column             | Type   | Description                                                                     |
| ------------------ | ------ | ------------------------------------------------------------------------------- |
| `sample_a`         | string | First sample of the pair.                                                       |
| `sample_b`         | string | Second sample of the pair.                                                      |
| `n_shared_spacers` | int    | Exact-sequence spacer matches.                                                  |
| `jaccard`          | float  | Jaccard index over the two unique spacer sets, six decimal places.              |
| `n_shared_repeats` | int    | Shared consensus repeats.                                                       |

```tsv
sample_a	sample_b	n_shared_spacers	jaccard	n_shared_repeats
gut_A	soil_B	7	0.063700	2
```

### `cohort_spacer_redundancy.tsv`

One row per spacer sequence that occurs in at least `--cohort_min_spacer_redundancy` samples
(default `2`). Rows are sorted by descending sample count, then descending occurrence count,
then spacer sequence.

| Column                | Type   | Description                                                                        |
| --------------------- | ------ | ----------------------------------------------------------------------------------- |
| `spacer`              | string | The spacer sequence.                                                                |
| `spacer_len`          | int    | Its length.                                                                         |
| `n_samples`           | int    | Number of samples carrying it.                                                      |
| `samples`             | string | Semicolon-separated sample identifiers, sorted.                                     |
| `n_occurrences`       | int    | Total occurrences across the cohort. A spacer can occur more than once per sample.  |
| `n_consensus_repeats` | int    | Number of distinct consensus repeats it was found next to.                          |
| `consensus_repeats`   | string | Semicolon-separated list of those repeats, sorted.                                  |

```tsv
spacer	spacer_len	n_samples	samples	n_occurrences	n_consensus_repeats	consensus_repeats
AGGCTTAAGCCTGAGGCTAAGCTAGTTCAGGCTA	34	2	gut_A;soil_B	3	1	GTTTCAATCCACGCACCCGTGGGGGTGATCC
TTGACCTGAAGCTTAGCCTGAAGCTTGACCTGAA	34	2	gut_A;soil_B	2	2	GTTTCAATCCACGCACCCGTGGGGGTGATCC;CGGTTTATCCCCGCTGGCGCGGGGAACAC
```

### `cohort_repeat_families.tsv`

One row per distinct consensus repeat across the cohort. Rows are sorted by descending sample
count, then descending spacer count, then repeat sequence.

| Column              | Type   | Description                                   |
| ------------------- | ------ | --------------------------------------------- |
| `consensus_repeat`  | string | The repeat sequence.                          |
| `repeat_len`        | int    | Its length.                                   |
| `n_samples`         | int    | Number of samples carrying it.                |
| `samples`           | string | Semicolon-separated sample identifiers, sorted. |
| `n_arrays_total`    | int    | Arrays with this consensus across the cohort. |
| `n_spacers_total`   | int    | Spacer occurrences belonging to those arrays. |
| `n_spacers_unique`  | int    | Distinct spacer sequences within the family.  |
| `median_spacer_len` | float  | Median spacer length within the family, one decimal place. |

```tsv
consensus_repeat	repeat_len	n_samples	samples	n_arrays_total	n_spacers_total	n_spacers_unique	median_spacer_len
GTTTCAATCCACGCACCCGTGGGGGTGATCC	31	2	gut_A;soil_B	3	26	24	34.0
CGGTTTATCCCCGCTGGCGCGGGGAACAC	29	1	gut_A	1	2	2	32.0
```

### `cohort_spacers.unique.fasta`

The pooled per-sample spacer FASTAs, deduplicated on sequence with `seqkit rmdup -s`. Each
retained record keeps the header of its first occurrence, so one representative sample and array
remain visible. The full sample list for a spacer is in `cohort_spacer_redundancy.tsv`. Header
format is the per-sample format documented under
[`<sample>.spacers.fasta`](#samplespacersfasta).

Record order follows the concatenation order of the per-sample files, which is not stable between
runs. Checksums of this file are not reproducible.

```text
>gut_A|array7|sp1 repeat_len=31 spacer_len=34
AGGCTTAAGCCTGAGGCTAAGCTAGTTCAGGCTA
>gut_A|array8|sp1 repeat_len=27 spacer_len=32
TTGCCAGGTCAATCGCGTTAACCGTTAGCCAT
```

The records are the input to the cohort sharing analysis and can be used as query sequences for
BLAST, `bowtie2` or a viral-genome search to find protospacer matches against known mobile
elements.

## Arrays

| Path                                     | Format | Contents                                                        |
| ---------------------------------------- | ------ | ---------------------------------------------------------------- |
| `arrays/<sample>/CRISPR_Arrays_*.txt`    | text   | MCAAT's human-readable report, unmodified.                       |
| `arrays/<sample>/<sample>.arrays.tsv`    | TSV    | One row per spacer.                                              |
| `arrays/<sample>/<sample>.spacers.fasta` | FASTA  | Every spacer, one record each.                                   |
| `arrays/<sample>/parameters.json`        | JSON   | MCAAT's record of the parameters it resolved.                    |
| `arrays/<sample>/<sample>.mcaat.log`     | text   | MCAAT's stdout and stderr for the run.                           |
| `arrays/<sample>/graph/`                 | dir    | The succinct de Bruijn graph. Only with `--save_mcaat_graph`.    |
| `arrays/<sample>/cycles.txt`             | text   | Raw multicycles. Only with `--save_mcaat_cycles`.                |

`<sample>.provenance.tsv` is written by the same step and is not published. It holds one header
row and one data row: `sample` plus eleven `mcaat_*` columns echoed from `parameters.json` —
`mcaat_threads`, `mcaat_ram_gb`, `mcaat_autoclean`, `mcaat_threshold_multiplicity`,
`mcaat_cycle_min_length`, `mcaat_cycle_max_length`, `mcaat_low_abundance`,
`mcaat_spacer_min_length`, `mcaat_spacer_max_length`, `mcaat_repeat_min_length` and
`mcaat_repeat_max_length`. Each is `NA` when the key is absent. Three of them appear in
`cohort_summary.tsv`; `parameters.json` is published in full.

### `CRISPR_Arrays_*.txt`

MCAAT's native output, published verbatim. Large runs are split across numbered files. `# Arrays`
and `# Spacers` in the header are run-wide totals repeated in every split file; they are not
per-file counts and must not be summed.

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

Grammar:

- The preamble is five lines: the title line (containing a literal UTF-8 em dash), `# Generated`,
  `# Arrays`, `# Spacers`, and a blank line.
- `>Array_N  spacers=<count>` opens a block, with two spaces before `spacers=`.
- The next line, at column 0, is the consensus repeat.
- Each following line is eight spaces, then the repeat variant observed for that spacer, then a
  single tab, then the spacer.
- A run of dashes in the variant column means the repeat is identical to the consensus. A literal
  sequence means the observed repeat differed; it may be shorter than the consensus.
- A blank line closes the block.

When a sample finds nothing, this file exists containing only the header, with `# Arrays : 0`.
MCAAT expresses a zero result in three ways, all exiting 0: no `CRISPR_Arrays_*.txt` because no
spacers were extracted, no `CRISPR_Arrays_*.txt` because `cycles.txt` was unreadable, or a
header-only `CRISPR_Arrays_1.txt` because all candidates were filtered out. The pipeline
synthesises a header-only `CRISPR_Arrays_1.txt` in the first two cases, so all three produce the
third shape.

### `<sample>.arrays.tsv`

One row per spacer.

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

Properties of the table:

- `array_id` is not a coordinate. Arrays are graph structures reconstructed from reads; there is
  no contig, no genome position and no strand.
- The repeat variant is block-constant: MCAAT records one repeat variant per array, not one per
  spacer, so `repeat_variant_is_consensus` is block-constant too. A block that mixes dashed and
  explicit variants produces a warning on stderr and the rows are still emitted.
- A non-consensus variant is frequently shorter than the consensus. The consensus is assembled
  from a tail and a front segment; the observed variant is the untailed repeat.
- Input files are read in numeric order of the `_<N>` filename suffix.
- Rows whose array block declares a `spacers=N` that disagrees with the number of parsed spacer
  lines are emitted, with a warning on stderr.
- A sample with zero arrays produces this file with a header row and no data rows.
- Tabs, carriage returns and newlines inside a field value are replaced with spaces.

### `<sample>.spacers.fasta`

Every spacer, unwrapped, one sequence per line. The header is:

```text
>{sample}|array{N}|sp{i} repeat_len={bp} spacer_len={bp}
```

The identifier is `{sample}|array{N}|sp{i}`, followed by a space and two `key=value` description
attributes. `repeat_len=` is the length of this spacer's repeat variant, not of the consensus.
The consensus length is `consensus_repeat_len` in `<sample>.arrays.tsv`. For an array whose
observed variant is a truncation of the consensus, the two differ.

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

In the example, array 7's variant equals its consensus, so `repeat_len=31` equals
`consensus_repeat_len`. Array 8 has `consensus_repeat_len` 29 and a 27 bp `repeat_variant`
(`CGGTTTATCCCCGCTGGCGCGGGGAAC`), and the header carries 27. Joining `repeat_len=` against
`consensus_repeat_len` therefore mis-joins. Join on the identifier or on `repeat_variant`.

A pattern matching against these headers must allow for the trailing description; a pattern
anchored at end-of-line immediately after `sp{i}` matches no record. The identifier ends at the
first space, so most tools carry through `gut_A|array7|sp1`.

The identifier is stable across reruns for fixed input: `array_id` comes from MCAAT's global
numbering and `spacer_index` is the 1-based position inside the block.

### `parameters.json`

MCAAT's record of the settings it resolved, published as provenance. The pipeline reads it back
and fails the sample when the recorded values disagree with the values the pipeline set. MCAAT's
argument parse loop has no unknown-argument branch, so a flag it does not understand is ignored
and the run proceeds with defaults; this readback is what detects that. Nine keys are compared:
`autoclean`, `cycle_finder.threshold_multiplicity`, `cycle_finder.cycle_min_length`,
`cycle_finder.cycle_max_length`, `cycle_finder.low_abundance`, `dna_sequence.spacer_min_length`,
`dna_sequence.spacer_max_length`, `dna_sequence.repeat_min_length` and
`dna_sequence.repeat_max_length`. The `threads` value is echoed into `<sample>.provenance.tsv`
but is not compared.

Two fields in the file are not used by anything: `output_file` names a `CRISPR_Arrays.txt` that
the release build never writes, and `cycles_folder` names a directory that is created but never
written to.

### `<sample>.mcaat.log`

MCAAT's stdout and stderr for the run: the `[1/4] … [4/4]` progress banners, the memory estimate
MCAAT printed for its own visited bitsets, and any error message. MCAAT redirects the graph
builder's stderr to `/dev/null`, so a failed run has an empty `.command.err` and the explanation
is in this file.

### `graph/`

Published only with `--save_mcaat_graph`. The publish mode is `--graph_publish_mode`, default
`symlink`, so by default the directory is a symlink into the work directory.

| File               | Contents                                                                          |
| ------------------ | --------------------------------------------------------------------------------- |
| `graph.sdbg_info`  | Graph metadata, naming the shards.                                                 |
| `graph.sdbg.<N>`   | Graph shards. The number of shards equals the thread count used, so the file set changes with `cpus`. |

With `--prune_graph true` (default), the graph-builder intermediates that MCAAT never reads back
are removed: `outfile_prefix.bin`, `outfile_prefix.lib_info`, `data.lib` and
`graph.mercy_cand.*`. The directory remains a valid value for the samplesheet's `sdbg` column
([graph re-entry](usage.md#graph-re-entry)) but is no longer a complete MEGAHIT working
directory. Set `--prune_graph false` to retain the intermediates.

### `cycles.txt`

Published only with `--save_mcaat_cycles`. One line per multicycle found in the graph. The file
can be very large. Line order derives from hash-map iteration, so two runs on identical input can
produce byte-different files. Checksums of this file are not reproducible.

## Phage curation

| Path                                            | Format | Contents                                                    |
| ----------------------------------------------- | ------ | ------------------------------------------------------------ |
| `phage/<sample>/<sample>.phage_contigs.fasta`   | FASTA  | Reconstructed contiguous paths carrying protospacer matches. |
| `phage/<sample>/phage_curate.log.gz`            | gzip   | The `mcaat phage-curate` run log. Only with `--save_phage_log`. |

Produced only with `--run_phage_curation` (default `false`). The stage wraps
`mcaat phage-curate`, a subcommand on the unreleased MCAAT v2.0.0 track. It walks the SDBG
outward from each spacer, so it requires `--mcaat_keep_graph true`; a missing or empty
`graph.sdbg_info` fails the task with exit code 64.

Every record is between 501 and 2500 bp. Both bounds are hard-coded in MCAAT.

MCAAT writes each record as `>phage_excerpt_<spacer sequence>`, which is neither unique nor
sample-identified. The pipeline rewrites the headers to:

```text
>gut_A|phage1 spacer=AGGCTTAAGCCTGAGGCTAAGCTAGTTCAGGCTA len=1183
TTGACAGCCTGAAGGTCTTAACGCCATTGCGGATCAGCTTTACCGGCAAAGTCGATCAGCT
...
```

| Attribute | Description                                                                                       |
| --------- | ------------------------------------------------------------------------------------------------- |
| `spacer=` | The bare nucleotide sequence of the array spacer this path was reconstructed around. Joins to `<sample>.spacers.fasta` and to the `spacer` column of `cohort_arrays.tsv`. |
| `len=`    | Contig length in base pairs.                                                                       |

A run that mapped zero arrays exits 0 and writes no FASTA, so `<sample>.phage_contigs.fasta` is
absent for that sample. The stage retries twice and is then ignored, so one sample cannot fail
the cohort; an ignored task appears in the Nextflow log. The contigs are graph walks around a
spacer match and carry no independent support.

## Read QC

| Path                                          | Format    | Contents                                                                   |
| --------------------------------------------- | --------- | --------------------------------------------------------------------------- |
| `qc/fastqc/raw/<sample>/*.html`, `*.zip`      | HTML, zip | FastQC on the reads as they arrive from the samplesheet.                   |
| `qc/fastqc/trimmed/<sample>/*.html`, `*.zip`  | HTML, zip | FastQC on the reads handed to MCAAT, after every enabled QC step including subsampling. |
| `qc/fastp/<sample>/*.json`, `*.html`, `*.log` | JSON, HTML, text | fastp trimming report.                                              |
| `qc/fastp/<sample>/*.fastp.fastq.gz`          | FASTQ.gz  | Trimmed reads. Only with `--save_trimmed_reads`.                           |
| `qc/fastp/<sample>/*.fail.fastq.gz`           | FASTQ.gz  | Reads that failed fastp filters. Only with `--fastp_save_trimmed_fail`.    |
| `qc/bowtie2/<sample>/*.log`                   | text      | Host-depletion alignment log.                                              |
| `qc/bowtie2/<sample>/*.fastq.gz`              | FASTQ.gz  | Host-depleted reads. Only with `--save_hostremoved_reads`.                 |
| `qc/bbduk/<sample>/*.log`                     | text      | BBDuk log from phiX removal.                                               |
| `qc/seqkit/<sample>/*.tsv`                    | TSV       | `seqkit stats --all` on the reads handed to MCAAT.                         |

The `seqkit` table also supplies the R1/R2 record-count assertion that keeps mis-paired libraries
from reaching MCAAT. Its `num_seqs` column is the value compared. When a sample was excluded, the
two disagreeing counts are in that file.

## MultiQC

| Path                            | Format | Contents                                    |
| ------------------------------- | ------ | -------------------------------------------- |
| `multiqc/multiqc_report.html`   | HTML   | The report.                                  |
| `multiqc/multiqc_data/`         | dir    | Parsed data behind every plot, as flat files. |
| `multiqc/multiqc_plots/`        | dir    | Exported static plot images.                 |

Alongside the standard QC sections, the report carries six MCAAT contributions. They are emitted
as `*_mqc.tsv` / `*_mqc.yaml` files carrying their own configuration, and the order of the five
sections is fixed by `custom_content.order` in `assets/multiqc_config.yml`.

| Id                     | Section name                      | Type               | Contents                                                                          |
| ---------------------- | --------------------------------- | ------------------ | ---------------------------------------------------------------------------------- |
| `mcaat_genstats`       | MCAAT CRISPR arrays               | General Statistics | Per-sample columns `CRISPR arrays` and `Spacers`. A sample that found nothing shows `0`, not a blank. Contributed per sample rather than as a section of its own, and configured by `table_columns_visible` / `table_columns_placement`. The file also carries `Unique spacers`, `Median spacer bp` and `Median repeat bp`, hidden by default. |
| `mcaat_cohort_summary` | MCAAT cohort summary              | table              | The `cohort_summary.tsv` table.                                                    |
| `mcaat_arrays_bargraph`| MCAAT arrays and spacers          | bargraph           | Arrays and unique spacers per sample.                                              |
| `mcaat_spacerlen`      | MCAAT spacer length distribution  | linegraph          | Spacer length distribution per sample.                                             |
| `mcaat_repeatlen`      | MCAAT repeat length distribution  | linegraph          | Repeat length distribution per sample, counting each array's consensus repeat once. |
| `mcaat_sharing`        | MCAAT spacer sharing              | heatmap            | Sample × sample matrix of exact spacer matches. The diagonal is set to zero.       |

`mcaat_cohort_summary` is written for every run, including an empty cohort. The length
distributions are skipped when no cohort arrays table was supplied. The sharing heatmap is
skipped, with a note on stderr, when the cohort has fewer than two samples; MultiQC renders a
degenerate single-cell plot otherwise.

## Pipeline information

| Path                                                | Format | Contents                                        |
| --------------------------------------------------- | ------ | ------------------------------------------------ |
| `pipeline_info/execution_report_<suffix>.html`      | HTML   | Nextflow execution report.                       |
| `pipeline_info/execution_timeline_<suffix>.html`    | HTML   | Nextflow timeline.                               |
| `pipeline_info/execution_trace_<suffix>.txt`        | text   | Nextflow trace, one row per task.                |
| `pipeline_info/pipeline_dag_<suffix>.html`          | HTML   | Nextflow DAG.                                    |
| `pipeline_info/mcaat_software_versions.yml`         | YAML   | Versions of every tool used.                     |
| `pipeline_info/params_<suffix>.json`                | JSON   | The fully resolved parameter set for the run.    |

`<suffix>` is `params.trace_report_suffix`, which defaults to the run start time formatted as
`yyyy-MM-dd_HH-mm-ss`.

`execution_trace_<suffix>.txt` carries the `peak_rss` column used by the memory-calibration
procedure in [`usage.md`](usage.md#calibrating-memory). The column is empty when the container
image in use lacks `ps`.
