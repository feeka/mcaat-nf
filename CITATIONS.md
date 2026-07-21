# RNABioInfo/mcaat-nf: Citations

## [MCAAT](https://github.com/RNABioInfo/mcaat)

> Talibli F, Voß B. MCAAT: Metagenomic CRISPR Array Analysis Tool. _microLife_. 2025. doi: [10.1093/femsml/uqaf016](https://doi.org/10.1093/femsml/uqaf016).

MCAAT is the tool this pipeline orchestrates. It detects CRISPR arrays directly in
un-assembled metagenomic reads by building a succinct de Bruijn graph and finding
multicycles in it. **If you use this pipeline, cite MCAAT.**

MCAAT vendors and links the following components, which should be cited alongside it
when the array-detection stage is used:

- **MEGAHIT** (succinct de Bruijn graph construction)

  > Li D, Liu CM, Luo R, Sadakane K, Lam TW. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. _Bioinformatics_. 2015 May 15;31(10):1674-6. doi: [10.1093/bioinformatics/btv033](https://doi.org/10.1093/bioinformatics/btv033). PubMed PMID: [25609793](https://pubmed.ncbi.nlm.nih.gov/25609793/).

- **SPOA** (partial order alignment, used for repeat/spacer consensus in post-processing)

  > Vaser R, Sović I, Nagarajan N, Šikić M. Fast and accurate de novo genome assembly from long uncorrected reads. _Genome Research_. 2017 May;27(5):737-746. doi: [10.1101/gr.214270.116](https://doi.org/10.1101/gr.214270.116). PubMed PMID: [28100585](https://pubmed.ncbi.nlm.nih.gov/28100585/).

## Pipeline framework

- [Nextflow](https://www.nextflow.io/)

  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. _Nature Biotechnology_. 2017 Apr 11;35(4):316-319. doi: [10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820). PubMed PMID: [28398311](https://pubmed.ncbi.nlm.nih.gov/28398311/).

- [nf-core](https://nf-co.re/)

  > Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. _Nature Biotechnology_. 2020 Feb 13;38(3):276-278. doi: [10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x). PubMed PMID: [32055031](https://pubmed.ncbi.nlm.nih.gov/32055031/).

  `mcaat-nf` is built **with** the nf-core template and tooling but is **not** an
  nf-core pipeline. Please do not refer to it as `nf-core/mcaat-nf`.

## Pipeline tools

- [BBMap / BBDuk](https://sourceforge.net/projects/bbmap/) — optional phiX and contaminant k-mer scrub.

  > Bushnell B. BBMap: A Fast, Accurate, Splice-Aware Aligner. _Lawrence Berkeley National Laboratory_. LBNL Report #: LBNL-7065E. 2014.

- [Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/) — optional host read depletion.

  > Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. _Nature Methods_. 2012 Mar 4;9(4):357-9. doi: [10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923). PubMed PMID: [22388286](https://pubmed.ncbi.nlm.nih.gov/22388286/).

- [csvtk](https://bioinf.shenwei.me/csvtk/) — cohort table concatenation.

  > Shen W, Sipos B, Zhao L. SeqKit2: A Swiss army knife for sequence and alignment processing. _iMeta_. 2024. doi: [10.1002/imt2.191](https://doi.org/10.1002/imt2.191).

- [fastp](https://github.com/OpenGene/fastp) — adapter and quality trimming.

  > Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. _Bioinformatics_. 2018 Sep 1;34(17):i884-i890. doi: [10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560). PubMed PMID: [30423086](https://pubmed.ncbi.nlm.nih.gov/30423086/).

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) — raw and trimmed read QC.

  > Andrews S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

- [MultiQC](https://multiqc.info/) — aggregated QC and CRISPR summary report.

  > Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. _Bioinformatics_. 2016 Oct 1;32(19):3047-8. doi: [10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354). PubMed PMID: [27312411](https://pubmed.ncbi.nlm.nih.gov/27312411/).

- [SeqKit](https://bioinf.shenwei.me/seqkit/) — read statistics (`seqkit stats`, used by the R1/R2 pair-sync gate) and cohort spacer de-duplication (`seqkit rmdup`).

  > Shen W, Le S, Li Y, Hu F. SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. _PLOS ONE_. 2016;11(10):e0163962. doi: [10.1371/journal.pone.0163962](https://doi.org/10.1371/journal.pone.0163962). PubMed PMID: [27706213](https://pubmed.ncbi.nlm.nih.gov/27706213/).

- [seqtk](https://github.com/lh3/seqtk) — deterministic paired subsampling (`--subsample_reads`).

- [Python](https://www.python.org/) — the three helper scripts this pipeline ships in
  `bin/`. All three are standard-library only; no third-party Python package is
  installed or required.

  > Van Rossum G, Drake FL. Python 3 Reference Manual. Scotts Valley, CA: CreateSpace; 2009.

  | Script                        | Role                                                                                                                     |
  | ----------------------------- | ------------------------------------------------------------------------------------------------------------------------ |
  | `bin/mcaat_parse_arrays.py`   | Parses MCAAT's `CRISPR_Arrays_*.txt` into `<sample>.arrays.tsv`, `<sample>.spacers.fasta` and `<sample>.provenance.tsv`, and performs the `parameters.json` readback assertion. |
  | `bin/mcaat_aggregate.py`      | Builds the cohort tables: cohort arrays, per-sample summary, pairwise spacer sharing, spacer redundancy and repeat families. |
  | `bin/mcaat_multiqc_sections.py` | Renders those cohort tables as MultiQC custom content (`mcaat_*_mqc.tsv` / `.yaml`).                                     |

## Reference data distributed with this pipeline

- `assets/phix174.fasta.gz` is the Enterobacteria phage phiX174 reference genome, NCBI accession [NC_001422.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1) (public domain). It is the default `--phix_reference` for the optional BBDuk scrub.

That is the only reference dataset shipped here.

## Software packaging and containerisation tools

- [Anaconda](https://anaconda.com)

  > Anaconda Software Distribution. _Computer software_. Vers. 2-2.4.0. Anaconda, 2016.

- [Bioconda](https://bioconda.github.io/)

  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. _Nature Methods_. 2018 Jul;15(7):475-476. doi: [10.1038/s41592-018-0046-7](https://doi.org/10.1038/s41592-018-0046-7). PubMed PMID: [29967506](https://pubmed.ncbi.nlm.nih.gov/29967506/).

- [BioContainers](https://biocontainers.pro/)

  > da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Vera Alvarez R, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. _Bioinformatics_. 2017 Aug 15;33(16):2580-2582. doi: [10.1093/bioinformatics/btx192](https://doi.org/10.1093/bioinformatics/btx192). PubMed PMID: [28379341](https://pubmed.ncbi.nlm.nih.gov/28379341/).

- [Docker](https://www.docker.com/)

  > Merkel D. Docker: lightweight Linux containers for consistent development and deployment. _Linux Journal_. 2014 Mar;2014(239):2.

- [Singularity / Apptainer](https://apptainer.org/)

  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. _PLOS ONE_. 2017;12(5):e0177459. doi: [10.1371/journal.pone.0177459](https://doi.org/10.1371/journal.pone.0177459). PubMed PMID: [28494014](https://pubmed.ncbi.nlm.nih.gov/28494014/).
