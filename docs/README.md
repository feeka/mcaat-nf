# RNABioInfo/mcaat-nf: Documentation

## Pages

| Page                     | Contents                                                                                                                                                                             |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| [Usage](usage.md)        | How the pipeline runs, the samplesheet format, every parameter, the state of the test data, resource settings for large metagenomes, and `-resume` behaviour.                          |
| [Output](output.md)      | Every file the pipeline publishes, with examples of the array table, the spacer FASTA and the cohort tables.                                                                           |

## Design notes in the source

Module boundaries and the reason each workaround exists are documented in the comment blocks of
`modules/local/**/main.nf`, next to the code they describe.

## Contributing

This repository has no `CONTRIBUTING.md`. Before starting work on a pull request:

1. Open an issue on the [repository tracker](https://github.com/RNABioInfo/mcaat-nf/issues).
2. Run `nf-test test`.
3. Run the pre-commit hooks.

## Relationship to nf-core

This pipeline is built with the [nf-core](https://nf-co.re) template and follows nf-core
conventions: module structure, `nextflow_schema.json`, nf-test, linting and the standard container
profiles. It is not an nf-core pipeline. It is developed and maintained in the `RNABioInfo`
organisation, outside the nf-core organisation and outside its review process.

- Refer to it as `RNABioInfo/mcaat-nf`, not as `nf-core/mcaat-nf`.
- It does not use the nf-core logo or branding.
- Send issues and support requests to
  [this repository's tracker](https://github.com/RNABioInfo/mcaat-nf/issues), not to nf-core Slack
  or the nf-core issue trackers.
- `.nf-core.yml` sets `is_nfcore: false`. Some nf-core lint checks are therefore configured
  differently from a pipeline inside the organisation.

## Licensing

| Component                            | Licence     | Distribution                                                          |
| ------------------------------------ | ----------- | --------------------------------------------------------------------- |
| This pipeline                        | GPL-3.0     | See [`LICENSE`](../LICENSE).                                          |
| MCAAT                                | GPL-3.0     | Invoked as an external tool from a container image; no MCAAT source code is vendored in this repository. |
| Reused nf-core template portions     | MIT         | Remain MIT-licensed by their authors; their notices are retained.     |

The nf-core community produced the template and tooling this repository reuses.

## Reference data in this repository

| File                        | Accession                                                                   | Licence       | Use                                              |
| --------------------------- | --------------------------------------------------------------------------- | ------------- | ------------------------------------------------ |
| `assets/phix174.fasta.gz`   | [NC_001422.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1)             | Public domain | Default `--phix_reference` for the BBDuk scrub.  |

This is the only reference file shipped in this repository.
