include { MULTIQC } from '../modules/nf-core/multiqc'
include { RNASEQ_STAR_SALMON } from '../subworkflows/local/rnaseq_star_salmon'
include { RNASEQ_DESEQ2_GPROFILER2 } from '../subworkflows/local/rnaseq_deseq2_gprofiler2'

workflow RNASEQ {
    take:
    samples

    main:
    if (!params.rnaseq_fasta || !params.rnaseq_gtf || !params.rnaseq_transcript_fasta) {
        error('RNA-seq mode requires --rnaseq_fasta, --rnaseq_gtf, and --rnaseq_transcript_fasta.')
    }
    if (!params.rnaseq_contrasts) {
        error('RNA-seq mode requires --rnaseq_contrasts.')
    }
    RNASEQ_STAR_SALMON(
        samples,
        params.rnaseq_fasta,
        params.rnaseq_gtf,
        params.rnaseq_transcript_fasta,
        params.rnaseq_star_index,
        params.rnaseq_salmon_index,
    )

    RNASEQ_DESEQ2_GPROFILER2(
        RNASEQ_STAR_SALMON.out.samplesheet,
        RNASEQ_STAR_SALMON.out.counts,
        RNASEQ_STAR_SALMON.out.features,
        params.rnaseq_contrasts,
    )

    ch_multiqc_files = RNASEQ_STAR_SALMON.out.fastqc
        .mix(RNASEQ_STAR_SALMON.out.star_log)
        .collect()
        .map { files -> [[id: 'rnaseq'], files, [], [], [], []] }

    MULTIQC(ch_multiqc_files)

    emit:
    multiqc_publish = MULTIQC.out.report
    multiqc_report = MULTIQC.out.report
}
