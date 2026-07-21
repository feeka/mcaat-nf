// ============================================================================
//  MCAAT_PHAGECURATE -- optional protospacer / phage contig curation
// ============================================================================
//  Wraps `mcaat phage-curate`. EXPERIMENTAL: the subcommand is on the
//  unreleased MCAAT v2.0.0 track and is disabled by default in this pipeline
//  (`params.run_phage_curation = false`).
//
//  Walks the SDBG outward from each spacer to reconstruct the contiguous
//  mobile-element context the spacer was acquired from.
// ============================================================================

process MCAAT_PHAGECURATE {
    tag "${meta.id}"
    label 'process_phage'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://' + params.mcaat_container
        : params.mcaat_container}"

    // ------------------------------------------------------------------
    // DO NOT "CLEAN THIS UP" -- see the long comment in
    // modules/local/mcaat/run/main.nf. The upstream published image sets
    // ENTRYPOINT ["mcaat"], which turns Nextflow's `/bin/bash -ue .command.sh`
    // into arguments to mcaat; they are silently discarded and the process
    // aborts with 134 without running the script at all. Fatal under
    // docker/podman, entirely invisible under singularity.
    // ------------------------------------------------------------------
    containerOptions { workflow.containerEngine in ['docker', 'podman'] ? '--entrypoint ""' : null }

    input:
    tuple val(meta), path(arrays, stageAs: 'arrays_in/*'), path(graph)

    output:
    // Glob rather than "${meta.id}.phage_contigs.fasta" so that a task.ext.prefix
    // override in conf/modules.config cannot desynchronise the declaration from
    // the filename the script actually writes.
    tuple val(meta), path("*.phage_contigs.fasta"), emit: fasta, optional: true
    tuple val(meta), path("phage_curate.log.gz"), emit: log
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # phage-curate exposes no --threads flag either, and in --graph mode
    # omp_set_num_threads() is never called at all, so OMP_NUM_THREADS is the
    # only thing standing between this task and every core on the node.
    export OMP_NUM_THREADS=${task.cpus}
    export OMP_PROC_BIND=false
    export LC_ALL=C.UTF-8

    # Written first so an early exit fails on its own code rather than on a
    # missing output. MCAAT has no --version flag, so the version is baked
    # into the image as ENV MCAAT_VERSION.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mcaat: \$(echo "\${MCAAT_VERSION:-1.0.0}")
    END_VERSIONS

    # Exactly this sample's arrays into a dedicated directory: array discovery
    # globs any filename merely CONTAINING 'CRISPR_Arrays_', non-recursively.
    mkdir -p arrays phage
    cp arrays_in/*.txt arrays/

    # phage-curate calls SDBG::LoadFromFile with no .sdbg_info existence check,
    # so a stale or wrong path gives undefined behaviour on a default-
    # constructed metadata object instead of a clean error.
    if [ ! -s ${graph}/graph.sdbg_info ]; then
        echo "MCAAT_PHAGECURATE: ${graph}/graph.sdbg_info is missing or empty." >&2
        echo "MCAAT_PHAGECURATE: the phage stage consumes the SDBG built by MCAAT_RUN;" >&2
        echo "MCAAT_PHAGECURATE: --run_phage_curation requires --mcaat_keep_graph true." >&2
        exit 64
    fi

    # ------------------------------------------------------------------
    # ContiguousPaths.fasta is opened with std::ios::app and is NEVER
    # truncated (src/protospacer_detection/phage_curator.cpp:419-420). A retry
    # into a reused directory therefore silently DOUBLES every record rather
    # than overwriting. Belt and braces, even though `phage/` is fresh here.
    # ------------------------------------------------------------------
    rm -f phage/ContiguousPaths.fasta

    # --output is always explicit: the subcommand defaults its output directory
    # to the --arrays directory, which would write ContiguousPaths.fasta into
    # the same folder the array stage publishes from.
    # stdin from /dev/null for the same reason as everywhere else.
    set +e
    mcaat phage-curate \\
        --arrays arrays \\
        --graph ${graph} \\
        --output phage \\
        < /dev/null > phage_curate.log 2>&1
    STATUS=\$?
    set -e

    # ------------------------------------------------------------------
    # The std::bad_alloc probe MUST run on the UNCOMPRESSED log, BEFORE the
    # gzip below. DO NOT rewrite it as `zcat phage_curate.log.gz | grep -q`:
    # process.shell sets `-o pipefail`, and `grep -q` exits 0 at the first
    # match, which kills zcat with SIGPIPE (141). pipefail then promotes 141
    # to the pipeline's status, so the test would read FALSE precisely when
    # the pattern IS present, and the task would fall through to the raw
    # status -- which conf/modules.config does not retry. The failure is
    # SIZE-DEPENDENT: a tiny test log fits the 64 KB pipe buffer and zcat
    # exits 0, so it passes in -profile test_phage and only bites on real data.
    # ------------------------------------------------------------------
    BAD_ALLOC=no
    if grep -qa 'std::bad_alloc' phage_curate.log; then
        BAD_ALLOC=yes
    fi

    gzip -f phage_curate.log

    if [ "\$STATUS" -ne 0 ]; then
        if [ "\$BAD_ALLOC" = "yes" ]; then
            echo "MCAAT_PHAGECURATE: phage-curate failed with std::bad_alloc -- remapping to 137." >&2
            exit 137
        fi
        exit "\$STATUS"
    fi

    # ------------------------------------------------------------------
    # Header rewrite.
    #
    # MCAAT writes each record as
    #     >phage_excerpt_<reconstructed spacer nucleotide sequence>
    # (src/protospacer_detection/phage_curator.cpp:570) -- a FIXED
    # 'phage_excerpt_' prefix followed by the spacer sequence itself, with no
    # array id, no coordinates and no uniqueness guarantee, which makes records
    # unaddressable in a cohort. Rewrite to
    #     >{sample}|phage{n} spacer={seq} len={bp}
    # keeping the original spacer as a searchable attribute. Every record is
    # 501-2500 bp by construction (both bounds are hard-coded upstream).
    #
    # THE `sub()` BELOW IS LOAD-BEARING, DO NOT DROP IT. `spacer=` is the only
    # linkage between a curated contig and the spacer it was reconstructed
    # from, so its value must be the bare nucleotide sequence and nothing else.
    # Leave the 'phage_excerpt_' prefix in and every join against
    # <sample>.spacers.fasta or the `spacer` column of cohort_arrays.tsv
    # matches zero records -- silently, because a 14-character prefix still
    # looks like a plausible identifier. `sub()` is a no-op if a future MCAAT
    # drops the prefix, so this stays correct either way.
    #
    # A run that mapped zero arrays exits 0 and writes nothing at all, which is
    # a legitimate negative -- hence `fasta` is declared optional.
    # ------------------------------------------------------------------
    if [ -s phage/ContiguousPaths.fasta ]; then
        awk -v s='${prefix}' '
            /^>/ { n++; sp = substr(\$0, 2); sub(/^phage_excerpt_/, "", sp); next }
            { print ">" s "|phage" n " spacer=" sp " len=" length(\$0); print \$0 }
        ' phage/ContiguousPaths.fasta > ${prefix}.phage_contigs.fasta
    else
        echo "MCAAT_PHAGECURATE: phage-curate produced no contiguous paths for ${prefix}." >&2
    fi
    """

    stub:
    // `prefix` must be re-declared here: `script:` and `stub:` are separate
    // closures, so the `def prefix` above is NOT in scope and referencing it
    // would fail every stub run -- including the offline stub run that
    // `nf-core pipelines download` performs in download_pipeline.yml.
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '>${prefix}|phage1 spacer=ACGTACGTACGTACGTACGTACGTACGTACGT len=8\\nACGTACGT\\n' \\
        > ${prefix}.phage_contigs.fasta
    echo "stub" | gzip -c > phage_curate.log.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mcaat: \$(echo "\${MCAAT_VERSION:-1.0.0}")
    END_VERSIONS
    """
}
