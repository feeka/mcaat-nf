# RNABioInfo/mcaat-nf: Documentation

The RNABioInfo/mcaat-nf documentation is split into the following pages:

- [Usage](usage.md)
  - An overview of how the pipeline works, how to run it, and a description of every parameter,
    including the samplesheet format, the state of the test data, resource tuning for large
    metagenomes, and `-resume` behaviour.
- [Output](output.md)
  - An overview of every file the pipeline publishes, with worked examples of the array table,
    the spacer FASTA and the cohort tables.

The rationale for the way the pipeline is shaped — the module boundaries and the reason each
workaround exists — lives next to the code it justifies, in the comment blocks of
`modules/local/**/main.nf`.

There is no `CONTRIBUTING.md` yet. Until there is, open an issue on the
[repository tracker](https://github.com/RNABioInfo/mcaat-nf/issues) before starting work on a
pull request, and run `nf-test test` plus the pre-commit hooks before submitting it.

## A note on nf-core

This pipeline is built **with** the [nf-core](https://nf-co.re) template and follows nf-core
conventions — module structure, `nextflow_schema.json`, nf-test, linting, the standard container
profiles — but it is **not an nf-core pipeline**. It is developed and maintained in the
`RNABioInfo` organisation, outside the nf-core organisation and outside its review process.

Concretely:

- Refer to it as `RNABioInfo/mcaat-nf`, never as `nf-core/mcaat-nf`.
- It does not use the nf-core logo or branding.
- Issues and support requests belong on
  [this repository's tracker](https://github.com/RNABioInfo/mcaat-nf/issues), not on nf-core
  Slack or the nf-core issue trackers.
- `.nf-core.yml` sets `is_nfcore: false`, which is why some nf-core lint checks are configured
  differently from a pipeline inside the organisation.

## A note on scope

The protospacer / phage curation stage is present, but it is experimental and off by default
(`--run_phage_curation`).

## Licensing

This pipeline is licensed **GPL-3.0** (see [`LICENSE`](../LICENSE)). MCAAT is also GPL-3.0 and is
invoked as an external tool from a container image; no MCAAT source code is vendored here. The
nf-core template portions this repository reuses remain MIT-licensed by their authors, and their
notices are retained. We gratefully acknowledge the nf-core community for that template and
tooling.

The only reference file shipped in this repository is `assets/phix174.fasta.gz` (NC_001422.1,
public domain), the default `--phix_reference` for the optional BBDuk scrub.
