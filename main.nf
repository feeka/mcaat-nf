#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RNABioInfo/mcaat-nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/RNABioInfo/mcaat-nf

    CRISPR array and protospacer discovery directly in unassembled metagenomic
    reads, powered by MCAAT.

    This pipeline is built WITH the nf-core template but is NOT an nf-core
    pipeline. Please do not refer to it as `nf-core/mcaat-nf`.
----------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MCAAT                   } from './workflows/mcaat'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_mcaat_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_mcaat_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline
//
workflow RNABIOINFO_MCAAT {

    take:
    ch_samplesheet    // channel: [ val(meta), [ reads ] ]                    -- rows that carry FASTQ
    ch_graph_reentry  // channel: [ val(meta), path(graph), [ arrays ] ]      -- `sdbg` rows

    main:

    //
    // WORKFLOW: Run pipeline
    //
    MCAAT(
        ch_samplesheet,
        ch_graph_reentry
    )

    emit:
    multiqc_report = MCAAT.out.multiqc_report // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    // Validates the parameters and the samplesheet, prints the summary, and
    // splits the samplesheet into read-bearing rows and `sdbg` graph re-entry
    // rows. All launch-time guards (mcaat_keep_graph implied by the optional
    // phage-curation stage and by --save_mcaat_graph, min-repeat-len >= k, the
    // min/max ordering checks) fire here, BEFORE any process is dispatched.
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    RNABIOINFO_MCAAT(
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.graph_reentry
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        RNABIOINFO_MCAAT.out.multiqc_report
    )
}
