//
// Subworkflow with functionality specific to the RNABioInfo/mcaat-nf pipeline
//
// This file owns:
//   * pipeline initialisation (version print, schema validation, config checks)
//   * ALL launch-time fail-fast guards (spec section 7.2) -- these run BEFORE a
//     single process is dispatched, so a misconfigured optional stage costs
//     seconds, not three hours of SDBG construction
//   * samplesheet parsing, single/paired-end derivation and lane grouping
//   * graph re-entry row resolution
//   * completion email / summary
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN   } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { samplesheetToList       } from 'plugin/nf-schema'
include { completionEmail         } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary       } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification          } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE   } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters.
    //
    // Everything that can be known without touching a read file is checked HERE.
    // The optional MCAAT stage is gated on an asset (a retained SDBG) whose
    // absence would otherwise only surface after the graph build -- i.e. hours
    // in, at full cost.
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //
    def ch_rows = channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map { meta, fastq_1, fastq_2, sdbg ->
            // nf-schema emits [] (falsy) for absent optional file/directory columns.
            tuple(meta, fastq_1 ?: [], fastq_2 ?: [], sdbg ?: [])
        }
        .branch { _meta, _fastq_1, _fastq_2, sdbg ->
            graph_reentry: sdbg as boolean
            reads: true
        }

    //
    // Graph re-entry rows (spec section 6.7).
    //
    // A row carrying `sdbg` skips MCAAT_RUN entirely. It is only coherent as a
    // RE-ANALYSIS of a previously published run, so the CRISPR_Arrays_*.txt that
    // the optional stage consumes are located by convention as siblings of the
    // graph/ directory, exactly as `<outdir>/arrays/<sample>/` lays them out.
    //
    def ch_graph_reentry = ch_rows.graph_reentry.map { meta, fastq_1, _fastq_2, sdbg ->
        //
        // Graph re-entry feeds the OPTIONAL phage-curation stage and nothing else.
        //
        // MCAAT_PARSEARRAYS is deliberately fed from MCAAT_RUN.out alone, because
        // a re-entry row carries no parameters.json and therefore cannot take part
        // in the readback assertion that is the pipeline's only detector for a
        // silently-dropped MCAAT flag. Consequently, with the optional stage off
        // there is no consumer for these rows at all: the graph is validated, the
        // sibling arrays are located, and then the sample is dropped on the floor.
        // The run would exit 0 with no <outdir>/arrays/<sample>/ and no row in
        // cohort_summary.tsv. Silent per-sample loss under a success exit code is
        // the worst failure shape available in a cohort study, so refuse it here.
        //
        if (!params.run_phage_curation) {
            error(
                "Samplesheet row '${meta.id}' supplies an 'sdbg' graph re-entry path, but\n" +
                "  --run_phage_curation is not enabled.\n" +
                "\n" +
                "  Graph re-entry only feeds the optional phage-curation stage: it skips MCAAT_RUN\n" +
                "  and carries no parameters.json, so it cannot enter array parsing or the cohort\n" +
                "  tables. With that stage off this row is a no-op and the sample would silently\n" +
                "  vanish from the results.\n" +
                "\n" +
                "  Either enable the optional stage:\n" +
                "      --run_phage_curation\n" +
                "  or remove the 'sdbg' column for this row and supply 'fastq_1'/'fastq_2' instead."
            )
        }

        if (fastq_1) {
            error(
                "Samplesheet row '${meta.id}' supplies BOTH an 'sdbg' directory and 'fastq_1'.\n" +
                "  Graph re-entry skips MCAAT_RUN entirely, so reads for that row would be silently ignored.\n" +
                "  Remove one of the two columns for this row."
            )
        }

        def graph_dir = file(sdbg, checkIfExists: true)

        if (!file("${graph_dir}/graph.sdbg_info").exists()) {
            error(
                "Samplesheet row '${meta.id}': '${graph_dir}' is not an MCAAT SDBG directory.\n" +
                "  Expected to find '${graph_dir}/graph.sdbg_info'.\n" +
                "  'mcaat phage-curate' does not validate the --graph path: SDBG::LoadFromFile is called\n" +
                "  with no metadata check (src/main.cpp:309-316), so a wrong path yields undefined\n" +
                "  behaviour rather than a clean error. This guard exists so that never happens."
            )
        }

        //
        // Locate the sibling CRISPR_Arrays_*.txt.
        //
        // The published layout is <outdir>/arrays/<sample>/{graph,CRISPR_Arrays_*.txt}:
        // conf/modules.config flattens BOTH the graph directory and the array
        // tables out of the task-relative `mcaat/` prefix with the same
        // `saveAs: { fn -> fn.tokenize('/')[-1] }`, so the arrays sit in the
        // graph's PARENT.
        //
        // The GRANDparent is searched as a fallback for output trees published
        // before the graph block in conf/modules.config grew its `saveAs`: those
        // have <sample>/mcaat/graph (parent = <sample>/mcaat, which holds ONLY the
        // graph) while the array tables were already flattened to <sample>/. The
        // parent alone also covers an sdbg pointed straight at a raw MCAAT task
        // work directory, where graph/ and CRISPR_Arrays_*.txt share one mcaat/
        // directory. Both candidates lie inside the SAME sample's tree, so
        // widening the search cannot pull in another sample's arrays. Nearest
        // candidate first; the first non-empty one wins.
        //
        def candidate_dirs = [graph_dir.getParent(), graph_dir.getParent()?.getParent()]
            .findAll { dir -> dir }
            .unique()

        // `file()` on a glob returns a List when it matches 0 or 2+ paths and a
        // bare Path when it matches exactly one, so every result is normalised
        // before it is tested for emptiness.
        def arrays_dir = candidate_dirs.find { dir ->
            def found = file("${dir}/CRISPR_Arrays_*.txt")
            (found instanceof List ? found : [found]) as boolean
        }

        def matched = arrays_dir ? file("${arrays_dir}/CRISPR_Arrays_*.txt") : []
        def arrays = matched instanceof List ? matched : [matched]

        if (!arrays) {
            error(
                "Samplesheet row '${meta.id}': no CRISPR_Arrays_*.txt found next to the graph directory.\n" +
                "  Looked in: ${candidate_dirs.join(', ')}\n" +
                "  Graph re-entry is a re-analysis path over a prior PUBLISHED run. A run launched with\n" +
                "  --save_mcaat_graph publishes the graph as <outdir>/arrays/<sample>/graph, alongside\n" +
                "  that sample's CRISPR_Arrays_*.txt -- point 'sdbg' at that directory."
            )
        }

        tuple([id: meta.id, single_end: true, from_graph: true], graph_dir, arrays)
    }

    //
    // Read rows: derive single/paired-endedness, strip `run`, merge lanes by sample.
    //
    def ch_samplesheet = ch_rows.reads
        .map { meta, fastq_1, fastq_2, _sdbg ->
            def sample_meta = [id: meta.id, single_end: !fastq_2]
            def fastqs = fastq_2 ? [fastq_1, fastq_2] : [fastq_1]
            // Group key is the sample id; `run` is deliberately dropped here and
            // never reaches MCAAT_RUN (spec section 4.4).
            tuple(meta.id, sample_meta, fastqs)
        }
        .groupTuple()
        .map { samplesheet -> validateInputSamplesheet(samplesheet) }

    emit:
    samplesheet   = ch_samplesheet
    graph_reentry = ch_graph_reentry
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    def summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal()
            )
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error("Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Check and validate pipeline parameters.
//
// Every check here is a LAUNCH-TIME check: it must be answerable from params
// alone (plus cheap filesystem stats) and it must fail before any task is
// submitted. Anything that can only be known from the data belongs downstream.
//
def validateInputParameters() {

    //
    // Graph retention. `mcaat_keep_graph` is the single documented cache-busting
    // toggle; the optional phage-curation stage consumes the SDBG that MCAAT_RUN
    // leaves behind.
    //
    if (params.run_phage_curation && !params.mcaat_keep_graph) {
        error("--run_phage_curation requires --mcaat_keep_graph true (the phage stage consumes the SDBG).")
    }

    if (params.save_mcaat_graph && !params.mcaat_keep_graph) {
        error("--save_mcaat_graph requires --mcaat_keep_graph true.")
    }

    //
    // k is hard-coded to 23 in the graph build (src/core/sdbg_build.cpp:224), so a
    // repeat shorter than k cannot be represented, let alone found.
    //
    if (params.mcaat_min_repeat_len < 23) {
        error("--mcaat_min_repeat_len cannot be below 23: k is hard-coded to 23 in the graph build.")
    }

    //
    // Interval sanity. MCAAT accepts inverted intervals without complaint and
    // then quietly returns nothing.
    //
    if (params.mcaat_min_spacer_len > params.mcaat_max_spacer_len) {
        error("--mcaat_min_spacer_len (${params.mcaat_min_spacer_len}) is greater than --mcaat_max_spacer_len (${params.mcaat_max_spacer_len}).")
    }

    if (params.mcaat_min_repeat_len > params.mcaat_max_repeat_len) {
        error("--mcaat_min_repeat_len (${params.mcaat_min_repeat_len}) is greater than --mcaat_max_repeat_len (${params.mcaat_max_repeat_len}).")
    }

    if (params.mcaat_cycle_min_length > params.mcaat_cycle_max_length) {
        error("--mcaat_cycle_min_length (${params.mcaat_cycle_min_length}) is greater than --mcaat_cycle_max_length (${params.mcaat_cycle_max_length}).")
    }

    //
    // Host depletion assets.
    //
    if (params.host_bowtie2_index && params.host_fasta) {
        log.warn("Both --host_bowtie2_index and --host_fasta were supplied; the pre-built index takes precedence and BOWTIE2_BUILD will be skipped.")
    }

    //
    // --skip_qc removes the pair-sync gate, which is the pipeline's only defence
    // against the single most likely real-world MCAAT failure.
    //
    if (params.skip_qc) {
        log.warn(
            "--skip_qc bypasses the entire read-QC subworkflow, INCLUDING the R1/R2 record-count gate.\n" +
            "  MCAAT treats two input files as a strictly synchronised paired-end library and megahit aborts\n" +
            "  on an odd combined read count with its stderr already redirected to /dev/null -- i.e. exit 1\n" +
            "  and a completely empty .command.err. Only use --skip_qc with single-end input or with reads\n" +
            "  you have already verified to be strictly paired."
        )
    }

    if (params.skip_pairsync_check && !params.skip_qc) {
        log.warn("--skip_pairsync_check disables the R1/R2 record-count assertion. See the --skip_qc warning above for what this costs you.")
    }

    //
    // Escape hatch. MCAAT's parse loop has NO unknown-argument branch
    // (src/main.cpp:611-612), so a typo here is accepted and dropped in silence.
    //
    if (params.mcaat_args) {
        log.warn(
            "--mcaat_args is a raw escape hatch: '${params.mcaat_args}'.\n" +
            "  MCAAT silently ignores unrecognised arguments, so a typo will NOT fail the run.\n" +
            "  The parameters.json readback in MCAAT_PARSEARRAYS is the only detector -- read it."
        )
    }
}

//
// Validate channels from the input samplesheet.
//
// Input is the grouped tuple [ id, [meta, ...], [fastqs, ...] ].
// All rows for one sample must agree on single/paired-endedness: MCAAT accepts
// exactly one or two files and unconditionally treats two as a synchronised
// paired-end library, so there is nowhere for a mixed set to go.
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    def endedness_ok = metas.collect { meta -> meta.single_end }.unique().size() == 1
    if (!endedness_ok) {
        error(
            "Please check input samplesheet -> Multiple runs of a sample must be of the same datatype, " +
            "i.e. all single-end or all paired-end: ${metas[0].id}"
        )
    }

    return [metas[0], fastqs]
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def citation_text = [
        "Tools used in the workflow included:",
        "MCAAT (Talibli et al.),",
        "FastQC (Andrews 2010),",
        params.skip_fastp ? "" : "fastp (Chen et al. 2018),",
        params.host_fasta || params.host_bowtie2_index ? "Bowtie2 (Langmead and Salzberg 2012)," : "",
        params.remove_phix ? "BBTools BBDuk," : "",
        params.subsample_reads ? "seqtk," : "",
        "SeqKit (Shen et al. 2016),",
        "csvtk,",
        "MultiQC (Ewels et al. 2016)",
        "."
    ].findAll { entry -> entry }.join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
        "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</li>",
        params.skip_fastp ? "" : "<li>Chen S, Zhou Y, Chen Y, Gu J. (2018) fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17):i884-i890. doi: 10.1093/bioinformatics/bty560</li>",
        params.host_fasta || params.host_bowtie2_index ? "<li>Langmead B, Salzberg SL. (2012) Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4):357-359. doi: 10.1038/nmeth.1923</li>" : "",
        "<li>Shen W, Le S, Li Y, Hu F. (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLoS ONE 11(10):e0163962. doi: 10.1371/journal.pone.0163962</li>",
        "<li>Ewels P, Magnusson M, Lundin S, Kaller M. (2016) MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19):3047-3048. doi: 10.1093/bioinformatics/btw354</li>"
    ].findAll { entry -> entry }.join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert to a named map so it can be used with the familiar NXF ${workflow} syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    def doi_ref = ""
    def manifest_doi = meta.manifest_map.doi
    if (manifest_doi) {
        def doi_list = manifest_doi instanceof List ? manifest_doi : manifest_doi.toString().tokenize(",")
        doi_ref = doi_list.collect { doi -> "https://doi.org/${doi.toString().trim().replace('https://doi.org/', '')}" }.join(', ')
    }
    meta["doi_text"] = doi_ref ? "(doi: <a href=\'${doi_ref}\'>${doi_ref}</a>)" : ""
    meta["nodoi_text"] = doi_ref ? "" : "<li>If available, make sure to update the text to reflect the actual DOI of the pipeline.</li>"

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()

    def methods_text = mqc_methods_yaml.text

    def engine = new groovy.text.GStringTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
