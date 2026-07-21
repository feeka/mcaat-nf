// ============================================================================
//  MCAAT_RUN -- CRISPR array detection directly in unassembled reads
// ============================================================================
//  Builds a succinct de Bruijn graph (vendored megahit SDBG) from one or two
//  FASTQ files and detects CRISPR arrays as multicycles in that graph.
//
//  This module is deliberately defensive. Every guard below exists because of
//  a specific, verified behaviour of the MCAAT v1.0.0 binary. Read the comment
//  before removing anything -- none of it is boilerplate.
// ============================================================================

process MCAAT_RUN {
    tag "${meta.id}"
    label 'process_mcaat_run'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://' + params.mcaat_container
        : params.mcaat_container}"

    // ------------------------------------------------------------------
    // DO NOT "CLEAN THIS UP".
    //
    // The upstream published image (docker.io/feeka94/mcaat:1.0.0) declares
    // `ENTRYPOINT ["mcaat"]`. Nextflow's Docker executor runs
    //     docker run ... <image> /bin/bash -ue .command.sh
    // and Docker PREPENDS the entrypoint, so the container actually executes
    //     mcaat /bin/bash -ue .command.sh
    // MCAAT's argument loop has no unknown-argument branch, so all three
    // tokens are silently discarded; with neither --input-files nor --graph
    // set it throws "No input provided" out of an unguarded main(), which
    // becomes std::terminate -> SIGABRT -> exit status 134. .command.sh is
    // never executed and no outputs are produced, so this presents as a tool
    // failure rather than a container failure.
    //
    // Critically, Singularity/Apptainer's `exec` ignores ENTRYPOINT, so this
    // bug is FATAL under -profile docker and COMPLETELY INVISIBLE under
    // -profile singularity. CI must exercise both.
    //
    // The pipeline-owned image sets no ENTRYPOINT, so this is a no-op there.
    // It is kept so that a user who overrides `container` with the upstream
    // image still gets a working run.
    // ------------------------------------------------------------------
    containerOptions { workflow.containerEngine in ['docker', 'podman'] ? '--entrypoint ""' : null }

    input:
    tuple val(meta), path(reads, stageAs: 'reads/*')

    output:
    tuple val(meta), path("mcaat/CRISPR_Arrays_*.txt"), emit: arrays
    tuple val(meta), path("mcaat/parameters.json"), emit: params_json
    tuple val(meta), path("mcaat.log"), emit: log
    tuple val(meta), path("mcaat/graph"), emit: graph, optional: true
    tuple val(meta), path("mcaat/cycles.txt"), emit: cycles, optional: true
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: (params.mcaat_args ?: '')
    // MCAAT accepts only the literal lowercase '1' / 'true' / 'yes' and
    // silently means FALSE for anything else -- including Groovy's 'True'.
    // Never interpolate a Boolean's toString() here.
    def low_abundance = (params.mcaat_low_abundance == null ? true : params.mcaat_low_abundance) ? 'true' : 'false'
    // Graph retention is controlled by ONE dedicated param, independent of the
    // optional-stage toggles, so that flipping --run_phage_curation later does
    // not change the resolved script text and bust every cached task.
    def keep_graph = (params.mcaat_keep_graph == null ? true : params.mcaat_keep_graph) ? 'yes' : 'no'
    def prune_graph = (params.prune_graph == null ? true : params.prune_graph) ? 'yes' : 'no'
    def thr_mult = params.mcaat_threshold_multiplicity ?: 20
    def cycle_min = params.mcaat_cycle_min_length ?: 27
    def cycle_max = params.mcaat_cycle_max_length ?: 77
    def spacer_min = params.mcaat_min_spacer_len ?: 23
    def spacer_max = params.mcaat_max_spacer_len ?: 50
    def repeat_min = params.mcaat_min_repeat_len ?: 23
    def repeat_max = params.mcaat_max_repeat_len ?: 55
    """
    export OMP_NUM_THREADS=${task.cpus}
    export OMP_PROC_BIND=false
    export LC_ALL=C.UTF-8

    # ------------------------------------------------------------------
    # versions.yml is written FIRST, before any step that can exit early,
    # so a task that fails a post-condition still fails on its own exit
    # code rather than on a confusing "Missing output file(s)".
    #
    # MCAAT has NO --version flag (verified: the version constant lives in
    # the never-called print_usage() at src/main.cpp:371-373 and is not
    # wired into the argument loop). The version is therefore read from
    # MCAAT_VERSION, which the pipeline-owned image bakes in with ENV.
    # Once `mcaat --version` exists upstream, replace the echo with a real
    # invocation.
    # ------------------------------------------------------------------
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mcaat: \$(echo "\${MCAAT_VERSION:-1.0.0}")
    END_VERSIONS

    # ------------------------------------------------------------------
    # Resource derivation is done HERE, in bash, and deliberately NOT
    # interpolated from task.memory.
    #
    # WHY NOT task.memory: interpolating it would put the memory directive
    # into the resolved script text, and hence into the task hash. Retuning
    # the process_mcaat_run label -- or any retry that scales memory -- would
    # then invalidate every cached MCAAT_RUN and force a full SDBG rebuild of
    # the entire cohort. Reading the limit at RUNTIME gives the same number
    # without touching the hash. DO NOT reintroduce a task.memory
    # interpolation here.
    #
    # WHY THE CGROUP AND NOT /proc/meminfo: --ram used to be inert (a gigabyte
    # figure was passed to megahit's byte-valued --host_mem), so the value we
    # sent did not matter and host memory was a safe thing to derive from.
    # As of the pinned 1.0.1 container that bug is FIXED and megahit now
    # honours the number. /proc/meminfo still reports HOST memory inside a
    # container -- cgroups are not virtualised there -- so deriving from it
    # would now hand megahit a budget many times the container's actual limit
    # and get the task OOM-killed. Read the cgroup limit instead, falling back
    # to host memory only when unconstrained.
    #
    # The clamp also keeps MCAAT's two validation traps unreachable: `ram <=
    # 1.0` and `ram > get_total_system_ram()` (src/main.cpp:513-520) both route
    # into the interactive stdin prompt or an uncaught throw. A cgroup limit is
    # always <= host total, so the upper trap cannot fire.
    # ------------------------------------------------------------------
    LIMIT_B=""
    if [ -r /sys/fs/cgroup/memory.max ]; then
        V=\$(cat /sys/fs/cgroup/memory.max 2>/dev/null || echo max)
        if [ "\$V" != "max" ]; then LIMIT_B="\$V"; fi
    fi
    if [ -z "\$LIMIT_B" ] && [ -r /sys/fs/cgroup/memory/memory.limit_in_bytes ]; then
        V=\$(cat /sys/fs/cgroup/memory/memory.limit_in_bytes 2>/dev/null || echo 0)
        # cgroup v1 reports a near-INT64_MAX sentinel when unlimited.
        if [ "\$V" -gt 0 ] 2>/dev/null && [ "\$V" -lt 4611686018427387904 ] 2>/dev/null; then
            LIMIT_B="\$V"
        fi
    fi
    if [ -z "\$LIMIT_B" ]; then
        LIMIT_B=\$(awk '/^MemTotal:/ {print \$2 * 1024}' /proc/meminfo 2>/dev/null || echo 0)
    fi
    RAM_G=\$(awk -v b="\$LIMIT_B" 'BEGIN { v = int(b / 1073741824 * 0.9); if (v < 2) v = 2; print v }')
    echo "MCAAT_RUN: memory limit \${LIMIT_B} B -> --ram \${RAM_G}G" >&2

    # MCAAT REJECTS --threads greater than std::thread::hardware_concurrency()
    # and routes that failure into the interactive prompt. nproc honours the
    # CPU affinity mask, so clamping here makes that branch unreachable.
    NPROC=\$(nproc 2>/dev/null || echo ${task.cpus})
    THREADS=${task.cpus}
    if [ "\$THREADS" -gt "\$NPROC" ]; then
        THREADS="\$NPROC"
    fi
    if [ "\$THREADS" -lt 1 ]; then
        THREADS=1
    fi

    # ------------------------------------------------------------------
    # The invocation itself.
    #
    #  * --output-folder is ALWAYS explicit. The default is a localtime-derived
    #    'mcaat_run_<timestamp>' relative to the CWD -- non-deterministic and
    #    timezone-dependent, which makes output globs and -resume unusable.
    #    It is also a dedicated SUBDIRECTORY, never the task work dir root, so
    #    that the stray fs::remove_all() below cannot reach staged inputs.
    #
    #  * < /dev/null is LOAD-BEARING. On any settings-validation failure the
    #    release main() prints "Remove it? (y/n): " and does `cin >> answer`
    #    on an UNINITIALISED char (src/main.cpp:919-930); one branch calls
    #    fs::remove_all() on the output folder. With a tty or a pipe on stdin
    #    the job blocks forever holding a scheduler slot; with /dev/null the
    #    extraction fails immediately at EOF.
    #
    #  * --autoclean is ALWAYS false. The default is true and recursively
    #    deletes <out>/graph, which the optional phage-curation stage needs.
    #    Cleanup is done in bash below instead, so the resolved script text
    #    stays invariant to the optional-stage toggles.
    #
    #  * --input-files is LAST. It greedily consumes every following token
    #    that does not begin with '-' (src/main.cpp:460-466), so anything
    #    placed after it would be swallowed into the file list and trip the
    #    "one or two input files" throw -- which ABORTS, it does not exit 1.
    # ------------------------------------------------------------------
    set +e
    mcaat \\
        --output-folder mcaat \\
        --threads "\$THREADS" \\
        --ram "\${RAM_G}G" \\
        --autoclean false \\
        --threshold-multiplicity ${thr_mult} \\
        --cycle-min-length ${cycle_min} \\
        --cycle-max-length ${cycle_max} \\
        --low-abundance ${low_abundance} \\
        --min-spacer-len ${spacer_min} \\
        --max-spacer-len ${spacer_max} \\
        --min-repeat-len ${repeat_min} \\
        --max-repeat-len ${repeat_max} \\
        ${args} \\
        --input-files ${reads} \\
        < /dev/null > mcaat.log 2>&1
    STATUS=\$?
    set -e

    cat mcaat.log

    # ------------------------------------------------------------------
    # Exit-code remapping.
    #
    # The release main() calls parse_arguments() WITHOUT a try/catch, so every
    # CLI and validation error escapes as std::terminate -> SIGABRT -> shell
    # status 134 -- which sits inside the nf-core default retry window and
    # would burn 4x walltime on a deterministic typo. A CycleFinder
    # std::bad_alloc ALSO aborts with 134 and IS worth retrying with more
    # memory, so a plain status filter must get one of the two wrong. We
    # disambiguate on the terminate handler's own message and fall through to
    # the raw status when neither pattern matches.
    #
    # NOTE these greps read a PLAIN FILE, never a `zcat | grep` pipeline:
    # process.shell sets `-o pipefail`, and `grep -q` exits at the first match,
    # killing the decompressor with SIGPIPE (141) -- which pipefail would then
    # promote to the pipeline's status, making the test read FALSE exactly when
    # the pattern IS present.
    #
    #    64 = deterministic MCAAT error    (never retried)
    #    65 = graph build produced nothing (never retried)
    #   137 = allocation failure           (retried, memory doubles)
    # ------------------------------------------------------------------
    if [ "\$STATUS" -ne 0 ]; then
        if grep -q 'std::bad_alloc' mcaat.log; then
            echo "MCAAT_RUN: MCAAT failed with std::bad_alloc -- remapping to 137 for retry." >&2
            exit 137
        fi
        if grep -qE 'std::runtime_error|std::invalid_argument|std::out_of_range|Unknown argument|^Error:' mcaat.log; then
            echo "MCAAT_RUN: MCAAT failed with a deterministic argument/validation error." >&2
            echo "MCAAT_RUN: read mcaat.log (published) or .command.out, NOT .command.err --" >&2
            echo "MCAAT_RUN: megahit stderr is dup2'd to /dev/null during the graph build." >&2
            exit 64
        fi
        exit "\$STATUS"
    fi

    # --- post-conditions ----------------------------------------------------
    if [ ! -f mcaat/parameters.json ]; then
        echo "MCAAT_RUN: MCAAT exited 0 but wrote no parameters.json." >&2
        exit 64
    fi
    if [ ! -s mcaat/graph/graph.sdbg_info ]; then
        echo "MCAAT_RUN: SDBG build produced no graph.sdbg_info." >&2
        echo "MCAAT_RUN: megahit stderr is redirected to /dev/null inside MCAAT, so" >&2
        echo "MCAAT_RUN: .command.err is EXPECTED to be empty. Check input FASTQ" >&2
        echo "MCAAT_RUN: integrity and R1/R2 record-count parity -- two input files are" >&2
        echo "MCAAT_RUN: unconditionally treated as a synchronised paired-end library" >&2
        echo "MCAAT_RUN: and megahit xfatal()s on an odd combined read count." >&2
        exit 65
    fi

    # ------------------------------------------------------------------
    # Zero-result normalisation.
    #
    # Finding no CRISPR arrays is a common and entirely legitimate metagenomic
    # result, but MCAAT expresses it in THREE different shapes, all exiting 0:
    #   (a) no spacers extracted        -> no CRISPR_Arrays_*.txt at all
    #   (b) cycles.txt unreadable       -> no CRISPR_Arrays_*.txt at all
    #   (c) all candidates filtered out -> a header-only CRISPR_Arrays_1.txt
    # We collapse (a) and (b) into (c) so the `arrays` output can be
    # non-optional, the downstream parser has exactly one code path, and a
    # true negative appears as an explicit zero row in the cohort table
    # instead of silently vanishing from the run.
    #
    # The em-dash is emitted as its raw UTF-8 bytes so the synthesised file is
    # byte-identical in shape to MCAAT's own output.
    # ------------------------------------------------------------------
    if ! ls mcaat/CRISPR_Arrays_*.txt >/dev/null 2>&1; then
        printf '# MCAAT \\xe2\\x80\\x94 CRISPR Array Output\\n# Generated : NA\\n# Arrays    : 0\\n# Spacers   : 0\\n\\n' \\
            > mcaat/CRISPR_Arrays_1.txt
        echo "MCAAT_RUN: MCAAT found no arrays; synthesised a header-only CRISPR_Arrays_1.txt." >&2
    fi

    # ------------------------------------------------------------------
    # Prune megahit intermediates.
    #
    # Verified against libs/megahit/src/sdbg/sdbg_raw_content.cpp:18-55 --
    # LoadSdbgRawContent opens ONLY <prefix>.sdbg_info and the <prefix>.sdbg.<N>
    # shards named inside that metadata. Everything removed below is read by no
    # MCAAT code path. outfile_prefix.bin in particular is a full 2-bit
    # re-encoding of every input read and is the one component that scales with
    # RAW INPUT rather than with graph complexity.
    #
    # The number of graph.sdbg.<N> shards equals --num_cpu_threads, which is
    # why the graph is declared as a DIRECTORY above and never as a shard glob.
    # ------------------------------------------------------------------
    if [ "${prune_graph}" = "yes" ]; then
        rm -f mcaat/graph/outfile_prefix.bin \\
              mcaat/graph/outfile_prefix.lib_info \\
              mcaat/graph/data.lib \\
              mcaat/graph/graph.mercy_cand.*
    fi

    # MCAAT creates <out>/cycles/ (src/main.cpp:695) and never writes into it.
    rmdir mcaat/cycles 2>/dev/null || true

    # --- graph retention: the ONE documented cache-busting toggle -----------
    if [ "${keep_graph}" != "yes" ]; then
        rm -rf mcaat/graph
    fi
    """

    stub:
    """
    # A stub: block is mandatory in every process. download_pipeline.yml runs
    # `nf-core pipelines download` with Singularity, attempts a STUB RUN first,
    # and fails the build if the Singularity image count increases during
    # execution ("The pipeline has no support for offline runs!").
    mkdir -p mcaat/graph
    printf '# MCAAT \\xe2\\x80\\x94 CRISPR Array Output\\n# Generated : NA\\n# Arrays    : 0\\n# Spacers   : 0\\n\\n' > mcaat/CRISPR_Arrays_1.txt
    echo '{}' > mcaat/parameters.json
    touch mcaat/graph/graph.sdbg_info mcaat/graph/graph.sdbg.0 mcaat/cycles.txt mcaat.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mcaat: \$(echo "\${MCAAT_VERSION:-1.0.0}")
    END_VERSIONS
    """
}
