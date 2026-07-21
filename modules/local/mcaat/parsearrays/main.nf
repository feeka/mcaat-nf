// MCAAT_PARSEARRAYS
//
// Turns MCAAT's human-readable `CRISPR_Arrays_*.txt` into machine-readable per-sample
// tables, and asserts that the parameters MCAAT actually recorded match the ones this
// pipeline intended.
//
// All format knowledge lives in `bin/mcaat_parse_arrays.py`. Nothing in this module --
// and nothing in any workflow or config -- may hard-code a column name or a grammar
// assumption about the MCAAT array output.
//
// The process name and the module path are load-bearing and must stay in lock-step:
// `workflows/mcaat.nf` includes this file as `../modules/local/mcaat/parsearrays/main`,
// and `conf/modules.config` keys `ext.prefix` and every publishDir block off the
// RESOLVED process name `MCAAT_PARSEARRAYS`. Renaming the process without renaming the
// selector silently drops the publishDir blocks and nothing lands in
// `<outdir>/arrays/<sample>/`; an `include { X as MCAAT_PARSEARRAYS }` alias does NOT
// fix that, because `withName:` matches the declared name, not the alias.
//
// Two behaviours are load-bearing and must not be "simplified" away:
//
//   1. The `arrays` input is staged into a DEDICATED `arrays_in/` directory. MCAAT's own
//      downstream stages glob any filename *containing* `CRISPR_Arrays_`, so array files
//      from two samples sharing one directory would be silently merged. Keeping the
//      convention here means array files are never renamed with a sample prefix.
//
//   2. The `expected_params` readback is the ONLY detector for a flag MCAAT accepted and
//      silently dropped: MCAAT's argument parse loop has no unknown-argument branch, so a
//      typo in `--mcaat_args` proceeds with defaults and no warning. A mismatch fails the
//      task on purpose.

process MCAAT_PARSEARRAYS {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'quay.io/biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(arrays, stageAs: 'arrays_in/*'), path(params_json)
    val expected_params

    output:
    tuple val(meta), path("*.arrays.tsv"),     emit: tsv
    tuple val(meta), path("*.spacers.fasta"),  emit: fasta
    tuple val(meta), path("*.provenance.tsv"), emit: provenance
    path "*_mcaat_genstats_mqc.tsv",           emit: mqc
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def expected_json = groovy.json.JsonOutput.toJson(expected_params ?: [:])
    def params_arg = params_json ? "--parameters ${params_json}" : ''
    """
    export LC_ALL=C.UTF-8
    export PYTHONIOENCODING=utf-8

    # The array files carry a literal UTF-8 em-dash in their preamble; a C locale would
    # make Python's default stdio encoding ASCII and mangle every warning we print.

    cat > expected_params.json <<'EXPECTED_JSON'
${expected_json}
EXPECTED_JSON

    # Nextflow puts \$projectDir/bin on PATH and bind-mounts it into the container. If the
    # checkout lost the executable bit (a Windows clone, or a zip download) fall back to
    # invoking the interpreter directly against the same bind-mounted path.
    if command -v mcaat_parse_arrays.py > /dev/null 2>&1; then
        PARSE="mcaat_parse_arrays.py"
    else
        PARSE="python3 ${projectDir}/bin/mcaat_parse_arrays.py"
    fi

    \$PARSE \\
        --sample '${prefix}' \\
        --arrays arrays_in \\
        ${params_arg} \\
        --expected expected_params.json \\
        --out-tsv '${prefix}.arrays.tsv' \\
        --out-fasta '${prefix}.spacers.fasta' \\
        --out-mqc '${prefix}_mcaat_genstats_mqc.tsv' \\
        --out-provenance '${prefix}.provenance.tsv' \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tarray_id\\tsource_file\\tconsensus_repeat\\tconsensus_repeat_len\\tn_spacers\\tspacer_index\\trepeat_variant\\trepeat_variant_is_consensus\\tspacer\\tspacer_len\\n' \\
        > '${prefix}.arrays.tsv'
    : > '${prefix}.spacers.fasta'
    printf 'sample\\tmcaat_threads\\tmcaat_ram_gb\\tmcaat_autoclean\\tmcaat_threshold_multiplicity\\tmcaat_cycle_min_length\\tmcaat_cycle_max_length\\tmcaat_low_abundance\\tmcaat_spacer_min_length\\tmcaat_spacer_max_length\\tmcaat_repeat_min_length\\tmcaat_repeat_max_length\\n%s\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\n' '${prefix}' \\
        > '${prefix}.provenance.tsv'
    {
        printf "# id: 'mcaat_genstats'\\n"
        printf "# plot_type: 'generalstats'\\n"
        printf 'Sample\\tCRISPR_arrays\\tCRISPR_spacers\\tCRISPR_spacers_unique\\tCRISPR_median_spacer_len\\tCRISPR_median_repeat_len\\n'
        printf '%s\\t0\\t0\\t0\\t0.0\\t0.0\\n' '${prefix}'
    } > '${prefix}_mcaat_genstats_mqc.tsv'
    """
}
