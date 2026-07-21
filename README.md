# RNABioInfo/mcaat-nf

[![GitHub Actions CI Status](https://github.com/RNABioInfo/mcaat-nf/actions/workflows/nf-test.yml/badge.svg)](https://github.com/RNABioInfo/mcaat-nf/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/RNABioInfo/mcaat-nf/actions/workflows/linting.yml/badge.svg)](https://github.com/RNABioInfo/mcaat-nf/actions/workflows/linting.yml)
[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

A Nextflow pipeline that runs [MCAAT](https://github.com/RNABioInfo/mcaat) over a set of samples.

MCAAT finds CRISPR arrays in raw, un-assembled metagenomic reads. It builds a succinct de Bruijn
graph and detects multicycles - the structural signature of CRISPR repeat-spacer arrays - without
any prior assembly step.

The pipeline adds read quality control, parses MCAAT's output into TSV and FASTA files, aggregates
spacers and repeats across samples, and produces a MultiQC report.

## Requirements

| Requirement      | Value                                                                        |
| ---------------- | ---------------------------------------------------------------------------- |
| Nextflow         | `>=25.10.4`                                                                   |
| Container engine | Docker, Singularity, Podman, Apptainer, Shifter or Charliecloud                |
| MCAAT image      | `docker.io/feeka94/mcaat:1.0.1` (`--mcaat_container`)                          |
| Architecture     | The MCAAT image is `linux/amd64` only; on ARM hosts use `-profile arm`         |

## Pipeline summary

1. **Merge sequencing runs** ([`cat`](https://www.gnu.org/software/coreutils/)) — samplesheet rows
   that share a `sample` value are concatenated per mate before any other step runs.
2. **Raw read QC** ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)).
3. **Adapter and quality trimming** ([`fastp`](https://github.com/OpenGene/fastp)) — strictly
   paired, merged output disabled.
4. **Host read removal**, optional ([`Bowtie2`](https://bowtie-bio.sourceforge.net/bowtie2/)) —
   enabled by `--host_fasta` or `--host_bowtie2_index`.
5. **phiX and contaminant removal**, optional
   ([`BBDuk`](https://jgi.doe.gov/data-and-tools/bbtools/)) — the phiX174 reference ships as
   `assets/phix174.fasta.gz`.
6. **Read subsampling**, optional ([`seqtk sample`](https://github.com/lh3/seqtk)) — controls
   MCAAT's runtime and memory requirements.
7. **Post-QC read statistics and pair-synchrony gate**
   ([`SeqKit`](https://bioinf.shenwei.me/seqkit/)) — samples whose R1 and R2 record counts disagree
   are excluded with a named message before MCAAT is dispatched.
8. **Trimmed read QC** ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)).
9. **CRISPR array detection** ([`MCAAT`](https://github.com/RNABioInfo/mcaat)) — succinct de Bruijn
   graph construction, multicycle detection, repeat and spacer resolution.
10. **Array parsing** — MCAAT's human-readable report is converted to a TSV plus a spacer FASTA, and
    the run's `parameters.json` is read back and asserted against the parameters the pipeline set.
11. **Protospacer and phage curation**, optional and experimental (`--run_phage_curation`) —
    reconstruction of contiguous protospacer-bearing paths from the retained graph.
12. **Cohort aggregation** ([`csvtk`](https://bioinf.shenwei.me/csvtk/),
    [`SeqKit`](https://bioinf.shenwei.me/seqkit/) and a bundled Python summariser) — cohort array
    table, unique spacer set, per-sample summary, pairwise spacer-sharing matrix, spacer redundancy
    table and repeat families.
13. **Report** ([`MultiQC`](http://multiqc.info/)).

### Phage curation stage

Step 11 wraps `mcaat phage-curate`, a subcommand on MCAAT's unreleased v2.0.0 track.

| Property             | Value                                                                                  |
| -------------------- | -------------------------------------------------------------------------------------- |
| Default state        | Disabled; enable with `--run_phage_curation`                                             |
| Additional required  | `--mcaat_keep_graph true` (the default). The stage takes the succinct graph as input.      |
| Search               | Recursive graph traversal with no depth cap and no wall-clock limit; runtime scales with graph density |
| Error strategy       | Retry twice, then ignore. An ignored task does not fail the run.                          |
| Output               | `phage/<sample>/<sample>.phage_contigs.fasta`, optional                                   |

When a sample has no `phage_contigs.fasta`, check the Nextflow log for an ignored task for that
sample.

## Usage

> [!NOTE]
> For Nextflow installation, see [this page](https://nf-co.re/docs/usage/installation).

> [!IMPORTANT]
> `-profile test` does not work yet. It downloads its test data from
> `https://raw.githubusercontent.com/RNABioInfo/mcaat-nf-testdata/main/`, which does not exist.
> The required files are listed in [`docs/usage.md`](docs/usage.md#test-profiles-and-test-data).
> Running on your own data does not need this.

### Samplesheet

`samplesheet.csv`:

```csv
sample,run,fastq_1,fastq_2,sdbg
gut_A,L001,/data/gut_A_L001_R1.fastq.gz,/data/gut_A_L001_R2.fastq.gz,
gut_A,L002,/data/gut_A_L002_R1.fastq.gz,/data/gut_A_L002_R2.fastq.gz,
soil_B,,/data/soil_B.fastq.gz,,
```

| Column    | Required | Description                                                                    |
| --------- | -------- | ------------------------------------------------------------------------------ |
| `sample`  | Yes      | Sample identifier. Rows sharing a value are merged per mate.                     |
| `run`     | No       | Run or lane identifier, used to distinguish rows of the same sample.             |
| `fastq_1` | Yes      | Path to the R1 FASTQ, or to the single-end FASTQ.                                |
| `fastq_2` | No       | Path to the R2 FASTQ. Omit for single-end data.                                  |
| `sdbg`    | No       | Path to an existing succinct de Bruijn graph directory, for graph re-entry.      |

Graph re-entry is described in [`docs/usage.md`](docs/usage.md#graph-re-entry).

### Running the pipeline

```bash
nextflow run RNABioInfo/mcaat-nf \
   -profile <docker/singularity/podman/apptainer/institute> \
   --input samplesheet.csv \
   --outdir results
```

A run with QC tuning and the phage curation stage enabled:

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
> Provide pipeline parameters via the CLI or a Nextflow `-params-file`. Custom config files,
> including those provided by the `-c` Nextflow option, can be used to provide any configuration
> _**except for parameters**_; see
> [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

Full parameter reference and tuning guidance: [usage documentation](docs/usage.md).

## Pipeline output

| Path                    | Contents                                                                       |
| ----------------------- | ------------------------------------------------------------------------------ |
| `qc/`                   | FastQC, fastp, Bowtie2, BBDuk and SeqKit outputs, one subdirectory per sample    |
| `arrays/<sample>/`      | MCAAT report, parsed array TSV, spacer FASTA, `parameters.json`, MCAAT log       |
| `phage/<sample>/`       | `<sample>.phage_contigs.fasta`; `phage_curate.log.gz` with `--save_phage_log`. Written only with `--run_phage_curation`. |
| `reports/`              | Cohort array table, per-sample summary, spacer sharing, spacer redundancy, repeat families, unique spacer FASTA |
| `multiqc/`              | MultiQC report                                                                   |
| `pipeline_info/`        | Execution timeline, report, trace and DAG                                        |

Every published file, its format and its columns are documented in the
[output documentation](docs/output.md).

## Credits

RNABioInfo/mcaat-nf was written by Fikrat Talibli (RNA Bioinformatics, University of Stuttgart), who
is also the author of MCAAT.

## Contributions and support

Report bugs and feature requests on the
[issue tracker](https://github.com/RNABioInfo/mcaat-nf/issues). There is no `CONTRIBUTING.md` yet;
until there is, open an issue before starting work on a pull request, and run `nf-test test` and the
pre-commit hooks before submitting it.

Bugs in the detection tool itself, rather than in the workflow around it, belong on the
[MCAAT issue tracker](https://github.com/RNABioInfo/mcaat/issues).

## Citations

If you use RNABioInfo/mcaat-nf for your analysis, cite the MCAAT paper:

> Talibli F, Voß B. **Metagenomic CRISPR Array Analysis Tool: a novel graph-based approach to finding CRISPR arrays in metagenomic datasets.** _microLife_, 2025. doi: [10.1093/femsml/uqaf016](https://doi.org/10.1093/femsml/uqaf016)

References for the tools used by the pipeline are listed in
[`CITATIONS.md`](CITATIONS.md).

This pipeline uses code and infrastructure developed and maintained by the
[nf-core](https://nf-co.re) community, reused here under the MIT licence.

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm,
> Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x)

## Licence

| Component                | Licence  | Terms                                                                           |
| ------------------------ | -------- | ------------------------------------------------------------------------------- |
| This pipeline            | GPL-3.0  | Everything in this repository. Full text in [`LICENSE`](LICENSE).                 |
| MCAAT                    | GPL-3.0  | Invoked as an external tool from a container image; no MCAAT source is vendored here. A GPL-3.0 work orchestrating a GPL-3.0 tool is compatible with no further conditions. |
| nf-core template portions | MIT     | The repository skeleton, the `utils_nfcore_*` subworkflows, the modules under `modules/nf-core/` and the schema and lint infrastructure remain MIT-licensed by the nf-core community. MIT is compatible with GPL-3.0, so redistributing them as part of this GPL-3.0 work is permitted; their MIT copyright notices are retained where they appear. |

Redistribute this pipeline under GPL-3.0 and keep all three notices intact.
