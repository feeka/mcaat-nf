# RNABioInfo/mcaat-nf

[![GitHub Actions CI Status](https://github.com/RNABioInfo/mcaat-nf/actions/workflows/nf-test.yml/badge.svg)](https://github.com/RNABioInfo/mcaat-nf/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/RNABioInfo/mcaat-nf/actions/workflows/linting.yml/badge.svg)](https://github.com/RNABioInfo/mcaat-nf/actions/workflows/linting.yml)
[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

> **This is not an nf-core pipeline.** It is built _with_ the nf-core template and follows
> nf-core conventions, but it is developed and maintained outside the nf-core organisation.
> Please refer to it as `RNABioInfo/mcaat-nf`, never as `nf-core/mcaat-nf`.

## Introduction

**RNABioInfo/mcaat-nf** is a bioinformatics pipeline for discovering **CRISPR arrays directly
in raw, unassembled metagenomic reads**, and for aggregating those discoveries across a cohort of
samples.

Detection is performed by [MCAAT](https://github.com/RNABioInfo/mcaat), which builds a succinct de
Bruijn graph over the reads and identifies _multicycles_ — the structural signature that a
repeat–spacer array leaves in a de Bruijn graph. No assembly, no binning, and no reference
database are involved, so arrays are recovered from low-abundance community members that would
never appear in an assembly.

The pipeline wraps that detection step in the things a real study needs around it: read quality
control, per-sample parsing into machine-readable tables, cross-sample spacer and repeat
aggregation, and a single MultiQC report.

### Scope

- **Protospacer / phage curation is experimental and off by default.** It is enabled only with
  `--run_phage_curation`, it wraps an unreleased MCAAT subcommand, and it is configured to be
  retried twice and then ignored rather than fail a cohort. Treat its output as a lead, not a
  result.
- **No spacer–protospacer database search.** The pipeline hands you a clean, addressable spacer
  FASTA; matching it against viral or plasmid databases is left to you.

## Pipeline summary

1. **Merge sequencing runs** ([`cat`](https://www.gnu.org/software/coreutils/)) — rows in the
   samplesheet that share a `sample` value are concatenated per mate before anything else runs.
2. **Raw read QC** ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)).
3. **Adapter and quality trimming** ([`fastp`](https://github.com/OpenGene/fastp)) — strictly
   paired, merged output disabled.
4. **Host read removal**, optional ([`Bowtie2`](https://bowtie-bio.sourceforge.net/bowtie2/)) —
   enabled by supplying `--host_fasta` or `--host_bowtie2_index`.
5. **phiX / contaminant removal**, optional ([`BBDuk`](https://jgi.doe.gov/data-and-tools/bbtools/)) —
   the phiX174 reference is shipped as `assets/phix174.fasta.gz`.
6. **Read subsampling**, optional ([`seqtk sample`](https://github.com/lh3/seqtk)) — the pipeline's
   principal lever on MCAAT's runtime and memory envelope.
7. **Post-QC read statistics and pair-synchrony gate**
   ([`SeqKit`](https://bioinf.shenwei.me/seqkit/)) — samples whose R1 and R2 record counts disagree
   are excluded with a named message _before_ MCAAT is dispatched.
8. **Trimmed read QC** ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)).
9. **CRISPR array detection** ([`MCAAT`](https://github.com/RNABioInfo/mcaat)) — succinct de Bruijn
   graph construction followed by multicycle detection and repeat/spacer resolution.
10. **Array parsing** — MCAAT's human-readable report is converted to a tidy TSV plus a spacer
    FASTA, and the run's `parameters.json` is read back and asserted against the parameters the
    pipeline intended to set.
11. **Protospacer / phage curation**, optional and **experimental** (`--run_phage_curation`) —
    reconstruction of contiguous protospacer-bearing paths from the retained graph.
12. **Cohort aggregation** ([`csvtk`](https://bioinf.shenwei.me/csvtk/),
    [`SeqKit`](https://bioinf.shenwei.me/seqkit/) and a bundled Python summariser) — cohort array
    table, unique spacer set, per-sample summary, pairwise spacer-sharing matrix, spacer
    redundancy table and repeat families.
13. **Report** ([`MultiQC`](http://multiqc.info/)).

### The one experimental stage

Step **11** wraps a MCAAT subcommand that is part of MCAAT's unreleased v2.0.0 track. It is
**disabled by default** and must be switched on explicitly with `--run_phage_curation`. It
consumes the succinct de Bruijn graph, so it additionally requires `--mcaat_keep_graph true`
(the default).

Its search is an uncapped recursive traversal with no wall-clock budget, so on a dense graph it
can be pathologically slow. The stage is therefore configured to retry twice and then be ignored:
one difficult sample cannot fail a cohort. If a sample has no `phage_contigs.fasta`, look for an
ignored task in the Nextflow log.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core-style pipelines, please refer to
> [this page](https://nf-co.re/docs/usage/installation) on how to set up Nextflow.

> [!IMPORTANT]
> `-profile test` is **not yet runnable.** Every test profile resolves its samplesheet from
> `params.pipelines_testdata_base_path`, which now points at
> `https://raw.githubusercontent.com/RNABioInfo/mcaat-nf-testdata/main/`. **That repository does
> not exist yet and has to be created and populated before any `test*` profile can run.**
> See [`docs/usage.md`](docs/usage.md#test-profiles-and-test-data) for the exact layout it needs.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,run,fastq_1,fastq_2,sdbg
gut_A,L001,/data/gut_A_L001_R1.fastq.gz,/data/gut_A_L001_R2.fastq.gz,
gut_A,L002,/data/gut_A_L002_R1.fastq.gz,/data/gut_A_L002_R2.fastq.gz,
soil_B,,/data/soil_B.fastq.gz,,
```

Each row is one FASTQ pair (or one single-end FASTQ). Rows that share a `sample` value are merged.
Omit `fastq_2` for single-end data. The `sdbg` column is only used for graph re-entry — see
[`docs/usage.md`](docs/usage.md#graph-re-entry).

Now, you can run the pipeline using:

```bash
nextflow run RNABioInfo/mcaat-nf \
   -profile <docker/singularity/podman/apptainer/institute> \
   --input samplesheet.csv \
   --outdir results
```

A complete run with QC tuning and the experimental stage enabled:

```bash
nextflow run RNABioInfo/mcaat-nf \
   -profile singularity \
   --input samplesheet.csv \
   --outdir results \
   --host_fasta /refs/GRCh38.fa.gz \
   --remove_phix \
   --subsample_reads 20000000 \
   --run_phage_curation \
   -resume
```

> [!WARNING]
> Provide pipeline parameters via the CLI or a Nextflow `-params-file`. Custom config files
> including those provided by the `-c` Nextflow option can be used to provide any configuration
> _**except for parameters**_; see
> [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](docs/usage.md).

## Pipeline output

The pipeline writes per-sample arrays, spacers and logs under `arrays/`, all cross-sample products
under `reports/`, and a single MultiQC report under `multiqc/`. For details of every published
file, see the [output documentation](docs/output.md).

## Credits

RNABioInfo/mcaat-nf was written by Fikrat Talibli (RNA Bioinformatics, University of Stuttgart),
who is also the author of MCAAT.

## Contributions and Support

Bugs and feature requests belong on the
[issue tracker](https://github.com/RNABioInfo/mcaat-nf/issues). There is no `CONTRIBUTING.md`
yet; until there is, please open an issue before starting work on a pull request, and run
`nf-test test` and the pre-commit hooks before submitting it.

Bugs in the **detection tool itself** (rather than in the workflow around it) belong on the
[MCAAT issue tracker](https://github.com/RNABioInfo/mcaat/issues).

## Citations

If you use RNABioInfo/mcaat-nf for your analysis, please cite the MCAAT paper:

> **Metagenomic CRISPR array analysis directly from unassembled reads.**
> _MicroLife_, 2025. doi: [10.1093/femsml/uqaf016](https://doi.org/10.1093/femsml/uqaf016)

An extensive list of references for the tools used by the pipeline can be found in the
[`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the
[nf-core](https://nf-co.re) community, reused here under the MIT licence.

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm,
> Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x)

## Licence

The licensing situation has three distinct parts, and it is worth stating each of them plainly:

- **This pipeline** — everything in this repository — is released under the
  **GNU General Public License v3.0**. The full text is in [`LICENSE`](LICENSE).
- **MCAAT** is licensed under **GPL-3.0** and is invoked as an external tool from a container
  image. No MCAAT source code is vendored in this repository. Because the pipeline is a GPL-3.0
  work that orchestrates a GPL-3.0 tool, the two are compatible with no further conditions.
- **The nf-core template portions** — the repository skeleton, the `utils_nfcore_*`
  subworkflows, the nf-core modules under `modules/nf-core/` and the schema/lint infrastructure —
  remain **MIT-licensed by their authors** (the nf-core community). MIT is compatible with
  GPL-3.0, so redistributing them as part of this GPL-3.0 work is permitted; their MIT copyright
  notices are retained where they appear.

If you redistribute this pipeline, redistribute it under GPL-3.0 and keep all three notices
intact.
