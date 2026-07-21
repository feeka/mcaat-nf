//
// READ_QC -- read quality control ahead of MCAAT array detection.
//
// Chain (every step individually skippable, all ON by default):
//
//   FASTQC_RAW -> FASTP -> [host depletion] -> [phiX scrub] -> [subsample]
//              -> FASTQC_TRIMMED + SEQKIT_STATS -> PAIR-SYNC GATE -> MCAAT
//
// The pair-sync gate at the end is the highest-value defensive construct in the
// whole pipeline. MCAAT hands two input files to megahit as a strictly
// SYNCHRONISED paired-end library (src/core/sdbg_build.cpp:66-77). megahit calls
// xfatal() -> exit(1) on an odd combined read count, and MCAAT has already
// dup2'd megahit's stderr to /dev/null for the entire build phase
// (src/core/sdbg_build.cpp:8-19). So the single most likely real-world failure
// -- desynchronised mates -- presents to the user as "exit status 1, empty
// .command.err", after the full cost of a graph build. We refuse to dispatch
// such a sample at all.
//

include { FASTQC as FASTQC_RAW     } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED } from '../../../modules/nf-core/fastqc/main'
include { FASTP                    } from '../../../modules/nf-core/fastp/main'
include { BOWTIE2_BUILD            } from '../../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN            } from '../../../modules/nf-core/bowtie2/align/main'
include { BBMAP_BBDUK              } from '../../../modules/nf-core/bbmap/bbduk/main'
include { SEQTK_SAMPLE             } from '../../../modules/nf-core/seqtk/sample/main'
include { SEQKIT_STATS             } from '../../../modules/nf-core/seqkit/stats/main'

workflow READ_QC {
    take:
    ch_reads // channel: [ val(meta), [ reads ] ]

    main:
    def ch_multiqc_files = channel.empty()

    //
    // MODULE: FastQC on the raw reads
    //
    if (!params.skip_fastqc) {
        FASTQC_RAW(ch_reads)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.map { _meta, zip -> zip })
    }

    //
    // MODULE: fastp -- adapter and quality trimming
    //
    // save_merged is HARD-WIRED false. A merged FASTQ would be a third file with
    // nowhere to go: MCAAT accepts exactly one or two inputs and infers
    // paired-endedness from the token count alone. Likewise reads that fail the
    // fastp filters are never routed onward, only optionally published.
    //
    def ch_after_fastp = ch_reads
    if (!params.skip_fastp) {
        FASTP(
            ch_reads.map { meta, reads -> [meta, reads, []] },
            false,
            params.fastp_save_trimmed_fail,
            false
        )
        ch_after_fastp = FASTP.out.reads
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map { _meta, json -> json })
    }

    //
    // MODULES: optional host depletion with Bowtie2
    //
    // A pre-built index takes precedence over the FASTA and skips BOWTIE2_BUILD.
    // save_unaligned = true, sort_bam = false: we consume .fastq, never the BAM.
    //
    def ch_after_host = ch_after_fastp
    if (params.host_fasta || params.host_bowtie2_index) {
        def ch_host_fasta = params.host_fasta
            ? channel.value([[id: 'host'], file(params.host_fasta, checkIfExists: true)])
            : channel.value([[id: 'host'], []])

        def ch_host_index = params.host_bowtie2_index
            ? channel.value([[id: 'host'], file(params.host_bowtie2_index, checkIfExists: true)])
            : BOWTIE2_BUILD(ch_host_fasta).index.first()

        BOWTIE2_ALIGN(
            ch_after_fastp,
            ch_host_index,
            ch_host_fasta,
            true,
            false
        )
        ch_after_host = BOWTIE2_ALIGN.out.fastq
        ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log.map { _meta, log_file -> log_file })
    }

    //
    // MODULE: optional phiX / small-contaminant scrub with BBDuk
    //
    // The reference must stay SMALL. BBDuk builds an in-memory k-mer table of the
    // reference; a mammalian genome here does not fit in the task's memory.
    //
    def ch_after_phix = ch_after_host
    if (params.remove_phix) {
        def ch_phix = channel.value(file(params.phix_reference, checkIfExists: true))
        BBMAP_BBDUK(ch_after_host, ch_phix)
        ch_after_phix = BBMAP_BBDUK.out.reads
        ch_multiqc_files = ch_multiqc_files.mix(BBMAP_BBDUK.out.log.map { _meta, log_file -> log_file })
    }

    //
    // MODULE: optional subsampling with seqtk
    //
    // This is the pipeline's ONLY lever on MCAAT's resource envelope. MCAAT
    // builds its SDBG with megahit's -m 1, so every error k-mer is retained and
    // graph size -- hence RAM, disk and the per-thread visited bitsets -- tracks
    // TOTAL SEQUENCED BASES rather than community complexity.
    //
    // The seed is applied identically to both mates (conf/modules.config sets
    // ext.args = "-s${params.subsample_seed}"), which is what preserves pairing.
    //
    def ch_reads_qc = ch_after_phix
    if (params.subsample_reads) {
        SEQTK_SAMPLE(ch_after_phix.map { meta, reads -> [meta, reads, params.subsample_reads] })
        ch_reads_qc = SEQTK_SAMPLE.out.reads
    }

    //
    // MODULE: FastQC on the processed reads
    //
    if (!params.skip_fastqc) {
        FASTQC_TRIMMED(ch_reads_qc)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.map { _meta, zip -> zip })
    }

    //
    // MODULE: SeqKit stats -- also the input to the pair-sync gate
    //
    SEQKIT_STATS(ch_reads_qc)
    ch_multiqc_files = ch_multiqc_files.mix(SEQKIT_STATS.out.stats.map { _meta, stats -> stats })

    //
    // PAIR-SYNC GATE
    //
    // SEQKIT_STATS emits one TSV row per input file. A paired sample must yield
    // exactly two rows with identical num_seqs; a single-end sample exactly one.
    // Anything else is excluded from dispatch with a named, actionable message
    // rather than being allowed to fail opaquely inside megahit.
    //
    def ch_gate = SEQKIT_STATS.out.stats
        .map { meta, tsv ->
            def rows = tsv.splitCsv(sep: '\t', header: true)
            def counts = rows.collect { row -> (row.num_seqs.toString().replaceAll(',', '')) as long }
            [meta, counts]
        }
        .branch { meta, counts ->
            pass: meta.single_end ? counts.size() == 1 : (counts.size() == 2 && counts[0] == counts[1])
            fail: true
        }

    ch_gate.fail.view { meta, counts ->
        "MCAAT-NF: EXCLUDING sample '${meta.id}' -- R1/R2 record counts differ (${counts}). " +
        "MCAAT treats two files as a strictly synchronised paired-end library and megahit " +
        "aborts on an odd total read count with an empty stderr. Re-run QC in strict paired mode."
    }

    def ch_reads_for_mcaat = params.skip_pairsync_check
        ? ch_reads_qc
        : ch_reads_qc.join(ch_gate.pass, by: 0).map { meta, reads, _counts -> [meta, reads] }

    emit:
    reads         = ch_reads_for_mcaat // channel: [ val(meta), [ reads ] ]
    multiqc_files = ch_multiqc_files   // channel: path
}
