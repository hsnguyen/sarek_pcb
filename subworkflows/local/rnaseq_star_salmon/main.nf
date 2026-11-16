include { FASTQC } from '../../../modules/nf-core/fastqc'
include { STAR_ALIGN } from '../../../modules/nf-core/star/align'
include { STAR_GENOMEGENERATE } from '../../../modules/nf-core/star/genomegenerate'
include { SALMON_INDEX } from '../../../modules/nf-core/salmon/index'
include { SALMON_QUANT } from '../../../modules/nf-core/salmon/quant'
include { RNASEQ_COUNTMATRIX } from '../../../modules/local/rnaseq/countmatrix'

workflow RNASEQ_STAR_SALMON {
    take:
    samples
    fasta
    gtf
    transcript_fasta
    star_index
    salmon_index

    main:
    ch_fasta = channel.fromPath(fasta, checkIfExists: true)
    ch_gtf = channel.fromPath(gtf, checkIfExists: true)
    ch_transcript_fasta = channel.fromPath(transcript_fasta, checkIfExists: true)

    if (star_index) {
        ch_star_index = channel.fromPath(star_index, checkIfExists: true).map { index -> [[:], index] }
    }
    else {
        STAR_GENOMEGENERATE(
            ch_fasta.map { reference -> [[:], reference] },
            ch_gtf.map { annotation -> [[:], annotation] },
        )
        ch_star_index = STAR_GENOMEGENERATE.out.index
    }

    if (salmon_index) {
        ch_salmon_index = channel.fromPath(salmon_index, checkIfExists: true).map { index -> [[:], index] }
    }
    else {
        SALMON_INDEX(
            ch_transcript_fasta.combine(ch_fasta).map { transcripts, reference -> [[:], transcripts, reference] },
        )
        ch_salmon_index = SALMON_INDEX.out.index
    }

    FASTQC(samples)

    STAR_ALIGN(
        samples,
        ch_star_index,
        ch_gtf.map { annotation -> [[:], annotation] },
        false,
    )

    SALMON_QUANT(
        samples,
        ch_salmon_index
            .combine(ch_gtf)
            .combine(ch_transcript_fasta)
            .map { index_meta, index, annotation, transcripts -> [index_meta, index, annotation, transcripts] }
            .first(),
    )

    ch_count_matrix_input = STAR_ALIGN.out.read_per_gene_tab
        .collect()
        .map { entries ->
            [
                [id: 'rnaseq'],
                entries.collect { entry -> entry[0] },
                entries.collect { entry -> entry[1] },
            ]
        }

    RNASEQ_COUNTMATRIX(ch_count_matrix_input)

    emit:
    fastqc = FASTQC.out.zip
    star_log = STAR_ALIGN.out.log_final
    star_bam = STAR_ALIGN.out.bam
    salmon = SALMON_QUANT.out.results
    samplesheet = RNASEQ_COUNTMATRIX.out.samplesheet
    counts = RNASEQ_COUNTMATRIX.out.counts
    features = RNASEQ_COUNTMATRIX.out.features
}
