include { PURECN_COVERAGE } from '../modules/nf-core/purecn/coverage'
include { PURECN_INTERVALFILE } from '../modules/nf-core/purecn/intervalfile'
include { PURECN_RUN } from '../modules/nf-core/purecn/run'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM } from '../modules/nf-core/samtools/convert'
include { SAREK as SAREK_NO_MUTECT2 } from './sarek'
include { SAREK as SAREK_MUTECT2 } from './sarek'

workflow PURITY {
    take:
    input_sample
    aligner
    skip_tools
    step
    ascat_alleles
    ascat_loci
    ascat_loci_gc
    ascat_loci_rt
    bbsplit_index
    bcftools_annotations
    bcftools_annotations_tbi
    bcftools_columns
    bcftools_header_lines
    cf_chrom_len
    chr_files
    cnvkit_reference
    dbsnp
    dbsnp_tbi
    dbsnp_vqsr
    dict
    fasta
    fasta_fai
    germline_resource
    germline_resource_tbi
    index_alignment
    intervals_and_num_intervals
    intervals_bed_combined
    intervals_bed_combined_for_variant_calling
    intervals_bed_gz_tbi_and_num_intervals
    intervals_bed_gz_tbi_combined
    intervals_for_preprocessing
    known_indels_vqsr
    known_sites_indels
    known_sites_indels_tbi
    known_sites_snps
    known_sites_snps_tbi
    known_snps_vqsr
    mappability
    msisensor2_models
    msisensorpro_scan
    ngscheckmate_bed
    mutect2_pon
    mutect2_pon_tbi
    sentieon_dnascope_model
    varlociraptor_chunk_size
    varlociraptor_events_germline
    varlociraptor_events_somatic
    varlociraptor_events_tumor_only
    varlociraptor_fdr
    varlociraptor_scenario_germline
    varlociraptor_scenario_somatic
    varlociraptor_scenario_tumor_only
    snpeff_cache
    snpeff_db
    vep_cache
    vep_cache_version
    vep_extra_files
    vep_fasta
    vep_genome
    vep_species
    snpsift_db

    main:
    if (!params.purecn_normaldb) {
        error('Purity mode requires --purecn_normaldb.')
    }

    ch_tumour_samples_for_purity = input_sample
        .filter { meta, _input, _input_index -> meta.status == 1 }
        .ifEmpty {
            error('Purity mode requires at least one tumour sample.')
        }

    ch_tumour_with_vcf = ch_tumour_samples_for_purity
        .filter { meta, _input, _input_index -> meta.vcf }

    ch_tumour_without_vcf = ch_tumour_samples_for_purity
        .filter { meta, _input, _input_index -> !meta.vcf }

    ch_vcf_from_samplesheet = ch_tumour_with_vcf
        .map { meta, _input, _input_index -> [meta.id, meta.vcf] }

    SAREK_NO_MUTECT2(
        ch_tumour_with_vcf.map { meta, input, input_index -> [meta - meta.subMap('vcf'), input, input_index] },
        aligner,
        skip_tools,
        step,
        'no_tools',
        ascat_alleles,
        ascat_loci,
        ascat_loci_gc,
        ascat_loci_rt,
        bbsplit_index,
        bcftools_annotations,
        bcftools_annotations_tbi,
        bcftools_columns,
        bcftools_header_lines,
        cf_chrom_len,
        chr_files,
        cnvkit_reference,
        dbsnp,
        dbsnp_tbi,
        dbsnp_vqsr,
        dict,
        fasta,
        fasta_fai,
        germline_resource,
        germline_resource_tbi,
        index_alignment,
        intervals_and_num_intervals,
        intervals_bed_combined,
        intervals_bed_combined_for_variant_calling,
        intervals_bed_gz_tbi_and_num_intervals,
        intervals_bed_gz_tbi_combined,
        intervals_for_preprocessing,
        known_indels_vqsr,
        known_sites_indels,
        known_sites_indels_tbi,
        known_sites_snps,
        known_sites_snps_tbi,
        known_snps_vqsr,
        mappability,
        msisensor2_models,
        msisensorpro_scan,
        ngscheckmate_bed,
        mutect2_pon,
        mutect2_pon_tbi,
        sentieon_dnascope_model,
        varlociraptor_chunk_size,
        varlociraptor_events_germline,
        varlociraptor_events_somatic,
        varlociraptor_events_tumor_only,
        varlociraptor_fdr,
        varlociraptor_scenario_germline,
        varlociraptor_scenario_somatic,
        varlociraptor_scenario_tumor_only,
        snpeff_cache,
        snpeff_db,
        vep_cache,
        vep_cache_version,
        vep_extra_files,
        vep_fasta,
        vep_genome,
        vep_species,
        snpsift_db,
    )

    SAREK_MUTECT2(
        ch_tumour_without_vcf,
        aligner,
        skip_tools,
        step,
        'mutect2',
        ascat_alleles,
        ascat_loci,
        ascat_loci_gc,
        ascat_loci_rt,
        bbsplit_index,
        bcftools_annotations,
        bcftools_annotations_tbi,
        bcftools_columns,
        bcftools_header_lines,
        cf_chrom_len,
        chr_files,
        cnvkit_reference,
        dbsnp,
        dbsnp_tbi,
        dbsnp_vqsr,
        dict,
        fasta,
        fasta_fai,
        germline_resource,
        germline_resource_tbi,
        index_alignment,
        intervals_and_num_intervals,
        intervals_bed_combined,
        intervals_bed_combined_for_variant_calling,
        intervals_bed_gz_tbi_and_num_intervals,
        intervals_bed_gz_tbi_combined,
        intervals_for_preprocessing,
        known_indels_vqsr,
        known_sites_indels,
        known_sites_indels_tbi,
        known_sites_snps,
        known_sites_snps_tbi,
        known_snps_vqsr,
        mappability,
        msisensor2_models,
        msisensorpro_scan,
        ngscheckmate_bed,
        mutect2_pon,
        mutect2_pon_tbi,
        sentieon_dnascope_model,
        varlociraptor_chunk_size,
        varlociraptor_events_germline,
        varlociraptor_events_somatic,
        varlociraptor_events_tumor_only,
        varlociraptor_fdr,
        varlociraptor_scenario_germline,
        varlociraptor_scenario_somatic,
        varlociraptor_scenario_tumor_only,
        snpeff_cache,
        snpeff_db,
        vep_cache,
        vep_cache_version,
        vep_extra_files,
        vep_fasta,
        vep_genome,
        vep_species,
        snpsift_db,
    )

    ch_sarek_for_purecn = SAREK_NO_MUTECT2.out.cram_variant_calling
        .mix(SAREK_MUTECT2.out.cram_variant_calling)
        .branch { _meta, alignment, _index ->
        bam: alignment.toString().endsWith('.bam')
        cram: alignment.toString().endsWith('.cram')
    }

    CRAM_TO_BAM(
        ch_sarek_for_purecn.cram,
        fasta.combine(fasta_fai).map { meta, fasta_file, _meta_fai, fai -> [meta, fasta_file, fai] }.collect(),
    )

    ch_cram_to_bam_for_purecn = CRAM_TO_BAM.out.bam
        .join(CRAM_TO_BAM.out.bai, failOnDuplicate: true, failOnMismatch: true)
        .mix(ch_sarek_for_purecn.bam)

    PURECN_INTERVALFILE(
        intervals_bed_combined.map { bed -> [[id: 'purecn_intervals'], bed] },
        fasta,
        params.purecn_genome ?: 'hg38',
    )

    PURECN_COVERAGE(
        ch_cram_to_bam_for_purecn,
        PURECN_INTERVALFILE.out.txt.collect { _meta, intervals -> intervals },
    )

    ch_vcf_for_purecn = ch_vcf_from_samplesheet
        .concat(SAREK_MUTECT2.out.mutect2_vcf.map { meta, vcf -> [meta.id, vcf] })

    ch_purecn_run_input = PURECN_COVERAGE.out.loess_txt
        .map { meta, coverage -> [meta.id, meta, coverage] }
        .join(ch_vcf_for_purecn, failOnDuplicate: true, failOnMismatch: true)
        .combine(PURECN_INTERVALFILE.out.txt.collect { _meta, intervals -> intervals })
        .map { _id, meta, coverage, vcf, intervals -> [meta, intervals, coverage, vcf] }

    PURECN_RUN(
        ch_purecn_run_input,
        channel.fromPath(params.purecn_normaldb, checkIfExists: true).collect(),
        params.purecn_mapping_bias
            ? channel.fromPath(params.purecn_mapping_bias, checkIfExists: true).collect()
            : channel.value([]),
        params.purecn_genome ?: 'hg38',
    )

    emit:
    multiqc_publish = channel.empty()
    multiqc_report = channel.empty()
    purity = PURECN_RUN.out.csv
}
