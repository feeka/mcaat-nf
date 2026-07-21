// COHORT_SUMMARISE
//
// Cohort layer. Merges every per-sample table emitted by MCAAT_PARSEARRAYS into the
// cross-sample deliverables, then renders them as MultiQC custom content.
//
// This is the step that produces science MCAAT cannot produce alone: spacer sharing
// between samples, repeat families spanning a cohort, and a redundancy table naming the
// exact spacer sequences that recur across the study.
//
// The process name and the module path are load-bearing and must stay in lock-step:
// `workflows/mcaat.nf` includes this file as `../modules/local/cohort_summarise/main`,
// and `conf/modules.config` keys the publishDir block off the RESOLVED process name
// `COHORT_SUMMARISE`. An `include { X as COHORT_SUMMARISE }` alias would NOT be enough --
// `withName:` matches the declared name, not the alias -- and the reports would silently
// never reach `<outdir>/reports/`.
//
// INPUT CONTRACT (exactly two collected file lists, in this order):
//
//   1. per-sample `<sample>.arrays.tsv`      from MCAAT_PARSEARRAYS.out.tsv.collect()
//   2. per-sample `<sample>.provenance.tsv`  from MCAAT_PARSEARRAYS.out.provenance.collect()
//
// BOTH slots take the PER-SAMPLE files, never a pre-concatenated cohort table: the
// per-sample filenames are the only surviving evidence of a sample that legitimately
// found zero arrays and therefore contributes no rows at all. The filenames themselves
// are written by `bin/mcaat_parse_arrays.py` via `--out-tsv`/`--out-provenance`, which
// MCAAT_PARSEARRAYS renders as `${prefix}.arrays.tsv` and `${prefix}.provenance.tsv`;
// the roster glob below matches exactly those two suffixes.
//
// NEVER pass MCAAT's raw `parameters.json` files here. Every sample emits a file with
// that exact name, so a cohort of more than one sample aborts at staging with
// "input file name collision"; and even for a single sample the file carries no sample
// identity. `*.provenance.tsv` exists precisely to solve both problems.
//
// Slot order is nonetheless not load-bearing for correctness: `bin/mcaat_aggregate.py`
// classifies everything it is given by FILENAME SUFFIX, so a mis-ordered call site
// produces correct tables rather than silently-`NA` provenance columns. The declared
// order is still the contract.
//
// Notes that matter:
//
//   * Every sample that reached this stage appears in `cohort_summary.tsv`, including
//     samples where MCAAT legitimately found zero arrays. A zero-array metagenome is a
//     common, valid result and must be an explicit zero row, never a missing row. The
//     sample roster is therefore established from three independent sources -- the
//     staged filenames (`--samples`), the provenance tables, and the data rows -- so a
//     zero-array sample survives even though it contributes no data rows anywhere.
//   * `cohort_spacers.unique.fasta` is declared `optional: true` because the shipped
//     `conf/modules.config` SETS `ext.args = '--no-unique-fasta'` on this process, so on
//     a default run this file is not written at all. That is deliberate and is half of a
//     two-part contract: SEQKIT_RMDUP publishes a file of the SAME NAME into the SAME
//     `<outdir>/reports` directory, and two processes copying different content to one
//     publish path is a race whose winner depends on task completion order. SEQKIT_RMDUP
//     owns the name because that is what docs/output.md describes. Suppressing this copy
//     is safe: the output is optional, so its absence is not a "Missing output file"
//     failure, and nothing in `workflows/mcaat.nf` consumes `out.unique_spacers`.
//     If you ever re-enable it, RENAME one of the two first.

process COHORT_SUMMARISE {
    tag "cohort"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'quay.io/biocontainers/python:3.12'}"

    input:
    path arrays_tsv, stageAs: 'arrays/*'
    path provenance_tsv, stageAs: 'provenance/*'

    output:
    path "cohort_arrays.tsv", emit: arrays
    path "cohort_summary.tsv", emit: summary
    path "cohort_spacer_sharing.tsv", emit: sharing
    path "cohort_spacer_redundancy.tsv", emit: redundancy
    path "cohort_repeat_families.tsv", emit: families
    path "cohort_spacers.unique.fasta", emit: unique_spacers, optional: true
    path "*_mqc.{tsv,yaml}", emit: mqc
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    // The two cohort thresholds are pipeline params, so they are rendered here rather
    // than left to conf/modules.config. `ext.args` is appended AFTER them: argparse
    // takes the last occurrence of an option, so an operator can still override either
    // value from a config without this module having to know about it.
    // `!= null` rather than the elvis operator: 0 is a legitimate (if unusual) threshold
    // and must not silently fall back to the default. Both params ARE declared in
    // nextflow.config and nextflow_schema.json with these same defaults -- the fallback
    // here is a belt-and-braces guard for a stripped-down `-params-file`, not the
    // primary source of the value.
    def min_sharing = params.cohort_min_spacer_sharing != null ? params.cohort_min_spacer_sharing : 2
    def min_redundancy = params.cohort_min_spacer_redundancy != null ? params.cohort_min_spacer_redundancy : 2
    """
    export LC_ALL=C.UTF-8
    export PYTHONIOENCODING=utf-8

    # Nextflow stages nothing for an empty input list, so these directories may not
    # exist. Both scripts treat an absent directory as empty, but creating them keeps
    # the command line uniform and the failure modes obvious.
    mkdir -p arrays provenance

    # Guard against the one wiring mistake this module cannot absorb: MCAAT's raw
    # parameters.json instead of the per-sample *.provenance.tsv. With more than one
    # sample Nextflow never even reaches this script (every file is called
    # parameters.json, so staging aborts with an input file name collision); with
    # exactly one sample it would run to completion and write NA into every provenance
    # column. Fail loudly instead, and say what to change.
    for stray in arrays/*.json provenance/*.json; do
        [ -e "\$stray" ] || continue
        echo "COHORT_SUMMARISE: '\$stray' was staged as a cohort input." >&2
        echo "COHORT_SUMMARISE: this module consumes exactly two collected channels," >&2
        echo "COHORT_SUMMARISE: MCAAT_PARSEARRAYS.out.tsv and" >&2
        echo "COHORT_SUMMARISE: MCAAT_PARSEARRAYS.out.provenance -- never MCAAT's raw" >&2
        echo "COHORT_SUMMARISE: parameters.json, which carries no sample identity and" >&2
        echo "COHORT_SUMMARISE: collides on staging as soon as there are two samples." >&2
        exit 64
    done

    # Sample roster from the STAGED FILENAMES. This is deliberately independent of the
    # file contents: a sample where MCAAT legitimately found zero arrays contributes a
    # header-only arrays table, and its identity would otherwise be lost before the
    # summary is written. The two globs are the exact filenames MCAAT_PARSEARRAYS asks
    # bin/mcaat_parse_arrays.py to write -- '<prefix>.arrays.tsv' and
    # '<prefix>.provenance.tsv' -- so stripping those two suffixes recovers meta.id.
    ROSTER=""
    for f in arrays/*.arrays.tsv provenance/*.provenance.tsv; do
        [ -e "\$f" ] || continue
        b=\$(basename "\$f")
        b=\${b%.arrays.tsv}
        b=\${b%.provenance.tsv}
        ROSTER="\${ROSTER}\${b},"
    done
    SAMPLES_ARG=""
    if [ -n "\$ROSTER" ]; then
        SAMPLES_ARG="--samples \${ROSTER%,}"
    fi

    if command -v mcaat_aggregate.py > /dev/null 2>&1; then
        AGGREGATE="mcaat_aggregate.py"
        SECTIONS="mcaat_multiqc_sections.py"
    else
        AGGREGATE="python3 ${projectDir}/bin/mcaat_aggregate.py"
        SECTIONS="python3 ${projectDir}/bin/mcaat_multiqc_sections.py"
    fi

    \$AGGREGATE \\
        --arrays-dir arrays \\
        --provenance-dir provenance \\
        \$SAMPLES_ARG \\
        --outdir . \\
        --min-spacer-sharing ${min_sharing} \\
        --min-spacer-redundancy ${min_redundancy} \\
        ${args}

    \$SECTIONS \\
        --summary cohort_summary.tsv \\
        --cohort-tsv cohort_arrays.tsv \\
        --sharing cohort_spacer_sharing.tsv \\
        --outdir . \\
        ${args2}
    """

    stub:
    """
    printf 'sample\\tarray_id\\tsource_file\\tconsensus_repeat\\tconsensus_repeat_len\\tn_spacers\\tspacer_index\\trepeat_variant\\trepeat_variant_is_consensus\\tspacer\\tspacer_len\\n' > cohort_arrays.tsv
    printf 'sample\\tn_arrays\\tn_spacers\\tn_spacers_unique\\tn_consensus_repeats\\tmedian_repeat_len\\tmedian_spacer_len\\tmin_spacer_len\\tmax_spacer_len\\tn_spacers_shared\\tpct_spacers_shared\\tn_repeats_shared\\tmcaat_threads\\tmcaat_ram_gb\\tmcaat_threshold_multiplicity\\n' > cohort_summary.tsv
    printf 'sample_a\\tsample_b\\tn_shared_spacers\\tjaccard\\tn_shared_repeats\\n' > cohort_spacer_sharing.tsv
    printf 'spacer\\tspacer_len\\tn_samples\\tsamples\\tn_occurrences\\tn_consensus_repeats\\tconsensus_repeats\\n' > cohort_spacer_redundancy.tsv
    printf 'consensus_repeat\\trepeat_len\\tn_samples\\tsamples\\tn_arrays_total\\tn_spacers_total\\tn_spacers_unique\\tmedian_spacer_len\\n' > cohort_repeat_families.tsv
    : > cohort_spacers.unique.fasta
    {
        printf "# id: 'mcaat_arrays_bargraph'\\n"
        printf "# plot_type: 'bargraph'\\n"
        printf 'Sample\\tCRISPR arrays\\tUnique spacers\\n'
    } > mcaat_arrays_bargraph_mqc.tsv
    {
        printf 'id: "mcaat_cohort_summary"\\n'
        printf 'plot_type: "table"\\n'
        printf 'data: {}\\n'
    } > mcaat_cohort_summary_mqc.yaml
    """
}
