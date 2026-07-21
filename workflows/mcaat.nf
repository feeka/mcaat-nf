/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MCAAT -- CRISPR array and protospacer discovery in unassembled metagenomic
             reads.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    read QC -> MCAAT array detection -> [phage curation]
            -> cohort aggregation -> MultiQC

    The stage in brackets is OPTIONAL and DISABLED BY DEFAULT. It is gated by a
    plain `if (params.run_phage_curation)` in the workflow body -- no dynamic
    dispatch, no lazy includes -- so when it is off, not one task, container pull
    or channel operator belonging to it is ever created. Its one prerequisite (a
    retained SDBG) is rejected at LAUNCH TIME by validateInputParameters() in
    subworkflows/local/utils_nfcore_mcaat_pipeline, before a single read is
    staged.

    This pipeline wraps MCAAT's array-detection and phage-curation subcommands.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { READ_QC                } from '../subworkflows/local/read_qc'

include { MCAAT_RUN              } from '../modules/local/mcaat/run/main'
include { MCAAT_PARSEARRAYS      } from '../modules/local/mcaat/parsearrays/main'
include { MCAAT_PHAGECURATE      } from '../modules/local/mcaat/phagecurate/main'
include { COHORT_SUMMARISE       } from '../modules/local/cohort_summarise/main'

include { CAT_FASTQ              } from '../modules/nf-core/cat/fastq/main'
include { CSVTK_CONCAT           } from '../modules/nf-core/csvtk/concat/main'
include { SEQKIT_RMDUP           } from '../modules/nf-core/seqkit/rmdup/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_mcaat_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MCAAT {
    take:
    ch_samplesheet   // channel: [ val(meta), [ [fastq_1, fastq_2], ... ] ]  grouped by sample
    ch_graph_reentry // channel: [ val(meta), path(graph), [ CRISPR_Arrays_*.txt ] ]

    main:
    def ch_multiqc_files = channel.empty()
    def ch_versions = channel.empty()

    //
    // MODULE: merge multi-run / multi-lane samples
    //
    // CAT_FASTQ runs BEFORE MCAAT, never after: MCAAT accepts exactly one or two
    // input files and unconditionally treats two as a synchronised paired-end
    // library (src/core/sdbg_build.cpp:66-77). There is nowhere for a third file
    // to go, so lanes must be concatenated first.
    //
    def ch_fastq = ch_samplesheet.branch { _meta, fastqs ->
        single: fastqs.size() == 1
        multiple: true
    }

    CAT_FASTQ(ch_fastq.multiple.map { meta, fastqs -> [meta, fastqs.flatten()] })

    def ch_reads_raw = ch_fastq.single
        .map { meta, fastqs -> [meta, fastqs.flatten()] }
        .mix(CAT_FASTQ.out.reads)

    //
    // SUBWORKFLOW: read QC
    //
    // Bypassed entirely by --skip_qc, which also removes the R1/R2 pair-sync
    // gate. That is documented as unsafe for two-file input and warned about at
    // launch time.
    //
    def ch_reads_for_mcaat = ch_reads_raw
    if (!params.skip_qc) {
        READ_QC(ch_reads_raw)
        ch_reads_for_mcaat = READ_QC.out.reads
        ch_multiqc_files = ch_multiqc_files.mix(READ_QC.out.multiqc_files)
    }

    //
    // MODULE: MCAAT array detection
    //
    // The script text of MCAAT_RUN is deliberately INVARIANT to the optional
    // stage toggle: --autoclean is hard-wired false, cleanup is done in bash,
    // and graph retention is controlled solely by params.mcaat_keep_graph. A
    // user who runs arrays-only today and re-runs with --run_phage_curation six
    // months from now resumes MCAAT_RUN from cache and pays only for the
    // protospacer chain.
    //
    MCAAT_RUN(ch_reads_for_mcaat)

    //
    // Unify freshly built graphs with samplesheet graph re-entry rows.
    //
    // A re-entry row skips MCAAT_RUN completely and brings its own previously
    // published CRISPR_Arrays_*.txt (resolved and validated at launch time).
    // Such rows feed the OPTIONAL stage only: they carry no parameters.json, so
    // they cannot take part in the readback assertion or the cohort tables.
    //
    def ch_graph_all = MCAAT_RUN.out.graph
        .mix(ch_graph_reentry.map { meta, graph, _arrays -> [meta, graph] })

    def ch_arrays_all = MCAAT_RUN.out.arrays
        .mix(ch_graph_reentry.map { meta, _graph, arrays -> [meta, arrays] })

    //
    // MODULE: parse the MCAAT array report
    //
    // `expected_params` is the intent side of the parameters.json readback.
    // MCAAT's argument parser has NO unknown-argument branch
    // (src/main.cpp:611-612), so a flag that was accepted and dropped is
    // otherwise completely undetectable. `threads` is deliberately absent here:
    // the module fills it in from its own task.cpus, which is the only place
    // that value is actually known after resourceLimits have been applied.
    //
    def ch_expected_params = channel.value([
        autoclean: false,
        threshold_multiplicity: params.mcaat_threshold_multiplicity,
        cycle_min_length: params.mcaat_cycle_min_length,
        cycle_max_length: params.mcaat_cycle_max_length,
        low_abundance: params.mcaat_low_abundance,
        spacer_min_length: params.mcaat_min_spacer_len,
        spacer_max_length: params.mcaat_max_spacer_len,
        repeat_min_length: params.mcaat_min_repeat_len,
        repeat_max_length: params.mcaat_max_repeat_len,
    ])

    MCAAT_PARSEARRAYS(
        MCAAT_RUN.out.arrays.join(MCAAT_RUN.out.params_json, by: 0),
        ch_expected_params
    )

    ch_multiqc_files = ch_multiqc_files.mix(MCAAT_PARSEARRAYS.out.mqc)

    //
    // MODULE: optional protospacer / phage curation (EXPERIMENTAL, off by default)
    //
    // Kept as a bare module rather than a subworkflow because it is a single
    // process with no fan-out. Its errorStrategy degrades to 'ignore' after two
    // attempts: DepthLimitedSearch is an uncapped recursive DFS followed by two
    // O(n^2) subpath filters with no candidate cap and no wall-clock budget, and
    // one pathological sample must never fail a cohort run that sits behind
    // .collect().
    //
    if (params.run_phage_curation) {
        MCAAT_PHAGECURATE(ch_arrays_all.join(ch_graph_all, by: 0))
    }

    //
    // Cohort aggregation
    //
    if (!params.skip_cohort_aggregation) {
        //
        // MODULE: concatenate every per-sample array table into one cohort table.
        //
        // This is a convenience deliverable for downstream users, NOT an input to
        // COHORT_SUMMARISE -- see the comment on that call below. It is published
        // under its own ext.prefix (conf/modules.config) so that it cannot collide
        // in <outdir>/reports with the cohort_arrays.tsv that COHORT_SUMMARISE
        // writes from the same source tables.
        //
        CSVTK_CONCAT(
            MCAAT_PARSEARRAYS.out.tsv
                .map { _meta, tsv -> tsv }
                .collect()
                .map { tsvs -> [[id: 'cohort'], tsvs] },
            'tsv',
            'tsv'
        )

        //
        // MODULE: pooled, sequence-deduplicated cohort spacer set.
        //
        // -s (dedup by sequence, not by name) is set via ext.args. The sharing
        // matrix is NOT derived from this -- it is computed in COHORT_SUMMARISE
        // from the per-sample array tables, so the analysis does not depend on an
        // output format we do not control.
        //
        SEQKIT_RMDUP(
            MCAAT_PARSEARRAYS.out.fasta
                .map { _meta, fasta -> fasta }
                .collectFile(name: 'cohort_spacers.fasta')
                .map { fasta -> [[id: 'cohort'], fasta] }
        )

        //
        // MODULE: cohort summary, spacer sharing matrix, repeat families.
        //
        // COHORT_SUMMARISE takes exactly TWO collected file lists, in this order:
        //
        //   1. every per-sample *.arrays.tsv      (MCAAT_PARSEARRAYS.out.tsv)
        //   2. every per-sample *.provenance.tsv  (MCAAT_PARSEARRAYS.out.provenance)
        //
        // Both are the PER-SAMPLE tables, not the concatenated cohort table: the
        // module rebuilds its sample roster from the staged FILENAMES so that a
        // sample in which MCAAT legitimately found zero arrays still gets an
        // explicit zero row. A pre-concatenated table has no filename per sample
        // and such a sample would silently vanish from cohort_summary.tsv.
        //
        // MCAAT's raw parameters.json must NEVER be passed here. Every sample
        // emits a file with that exact name, so a cohort of two or more samples
        // aborts at staging with "input file name collision"; *.provenance.tsv is
        // the per-sample, sample-identified replacement and is consumed here and
        // nowhere else.
        //
        COHORT_SUMMARISE(
            MCAAT_PARSEARRAYS.out.tsv.map { _meta, tsv -> tsv }.collect(),
            MCAAT_PARSEARRAYS.out.provenance.map { _meta, tsv -> tsv }.collect()
        )

        ch_multiqc_files = ch_multiqc_files.mix(COHORT_SUMMARISE.out.mqc.flatten())
    }

    //
    // Collate software versions.
    //
    // Version reporting is topic-channel based. `MODULE.out.versions` does not
    // exist on current nf-core/modules -- never mix it.
    //
    def ch_versions_topic = channel.topic('versions').distinct()

    def ch_versions_branched = ch_versions_topic.branch { entry ->
        versions_file: entry instanceof Path
        versions_tuple: true
    }

    // eval-based topic entries arrive as (process, tool, version) triples.
    def ch_versions_from_tuples = ch_versions_branched.versions_tuple.map { process, tool, version ->
        def process_name = process.toString().tokenize(':').last()
        "${process_name}:\n  ${tool}: ${version.toString().trim()}\n"
    }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(ch_versions_branched.versions_file))
        .mix(ch_versions_from_tuples)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'mcaat_software_versions.yml',
            sort: true,
            newLine: true,
        )

    //
    // MODULE: MultiQC
    //
    def ch_multiqc_report = channel.empty()
    if (!params.skip_multiqc) {
        // The pipeline's OWN config is always FIRST in the config slot. It carries
        // the report_section_order and the `custom_content: order:` list that place
        // the MCAAT sections after the QC sections; losing it silently reorders the
        // report.
        //
        // A user-supplied --multiqc_config is APPENDED after ours rather than
        // replacing it (see the MULTIQC call below for how the two are folded into
        // the module's list-capable config slot). MultiQC applies both and the
        // later --config wins key-by-key, which is exactly the "Merged with the
        // pipeline's own config" behaviour nextflow_schema.json documents.
        // Replacing ours outright would make that description a lie.
        def multiqc_config = file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)

        def extra_multiqc_config = params.multiqc_config
            ? file(params.multiqc_config, checkIfExists: true)
            : []

        def multiqc_logo = params.multiqc_logo ? file(params.multiqc_logo, checkIfExists: true) : []

        def summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        def ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))

        def methods_description_file = params.multiqc_methods_description
            ? file(params.multiqc_methods_description, checkIfExists: true)
            : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
        def ch_methods_description = channel.value(methodsDescriptionText(methods_description_file))

        ch_multiqc_files = ch_multiqc_files
            .mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
            .mix(ch_collated_versions)
            .mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

        // Current nf-core/multiqc takes a SINGLE tuple, not five positional
        // channels. Any five-argument example found in a released pipeline is
        // pinned to an older SHA.
        //
        // The tuple has SIX elements and the order is load-bearing
        // (modules/nf-core/multiqc/main.nf:11):
        //
        //   tuple val(meta), path(multiqc_files, stageAs: "?/*"),
        //         path(multiqc_config, stageAs: "?/*"),
        //         path(multiqc_logo), path(replace_names), path(sample_names)
        //
        // There is NO separate `extra_multiqc_config` slot in this module
        // version. Instead the config slot is itself list-capable -- the module
        // renders `--config <a> --config <b>` when it receives a List -- which is
        // how a user config is merged with ours rather than replacing it.
        //
        // Order inside that list is the whole point: the PIPELINE's own config
        // must come first and the user's second, because MultiQC lets the later
        // --config win key-by-key. Ours carries the report_section_order and the
        // `custom_content: order:` list; the user's is an override on top of it.
        // That is exactly the "Merged with the pipeline's own config" behaviour
        // nextflow_schema.json documents.
        //
        // Emitting the WRONG number of elements does not fail loudly -- Nextflow
        // destructures by position and logs only a cardinality warning. A seventh
        // element would shift the logo into `replace_names`, and MultiQC would
        // then try to parse a PNG as a two-column TSV and fail the task. Re-read
        // modules/nf-core/multiqc/main.nf after `nf-core modules update multiqc`
        // and re-check this arity before changing anything here.
        def multiqc_config_slot = extra_multiqc_config
            ? [multiqc_config, extra_multiqc_config]
            : multiqc_config

        MULTIQC(
            ch_multiqc_files.collect().map { files ->
                [
                    [id: 'mcaat'],
                    files,
                    multiqc_config_slot,
                    multiqc_logo,
                    [],
                    [],
                ]
            }
        )

        ch_multiqc_report = MULTIQC.out.report.toList()
    }

    emit:
    multiqc_report = ch_multiqc_report // channel: [ path(multiqc_report.html) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
