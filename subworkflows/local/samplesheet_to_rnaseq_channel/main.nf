workflow SAMPLESHEET_TO_RNASEQ_CHANNEL {
    take:
    ch_from_samplesheet

    main:
    ch_rnaseq_samples = ch_from_samplesheet
        .map { meta, fastq_1, fastq_2, _spring_1, _spring_2, _table, _cram, _crai, _bam, _bai, _contamination, _vcf, _variantcaller ->
            if (!fastq_1) {
                error("RNA sample '${meta.sample}' requires fastq_1.")
            }
            if (fastq_2 && !meta.lane) {
                error("Paired-end RNA sample '${meta.sample}' requires a lane value.")
            }
            [meta.sample, meta + [id: meta.sample.toString(), single_end: !fastq_2], fastq_1, fastq_2]
        }
        .groupTuple()
        .map { _sample, metas, fastq_1, fastq_2 ->
            def meta = metas.first() - metas.first().subMap('lane')
            def paired_end = !meta.single_end
            if (metas.any { item -> item.single_end != !fastq_2.first() }) {
                error("RNA sample '${meta.id}' mixes single-end and paired-end libraries.")
            }
            [meta, paired_end
                ? fastq_1.indices.collectMany { index -> [fastq_1[index], fastq_2[index]] }
                : fastq_1]
        }
        .ifEmpty {
            error('No RNA samples found in the input samplesheet.')
        }

    emit:
    ch_rnaseq_samples
}
