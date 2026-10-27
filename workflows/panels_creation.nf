include { BAM_CREATE_SOM_PON_GATK } from '../subworkflows/nf-core/bam_create_som_pon_gatk'
include { SAREK } from './sarek'

workflow PANELS_CREATION {
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
    ch_normal_samples_for_mutect2_pon = input_sample
        .filter { meta, _input, _input_index ->
            meta.status == 0 && meta.pon_db?.tokenize(';')?.contains('mutect2')
        }
        .ifEmpty {
            error("No normal samples request a Mutect2 panel of normals. Set 'pon_db' to 'mutect2' for at least one normal sample.")
        }

    SAREK(
        ch_normal_samples_for_mutect2_pon,
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

    ch_mutect2_pon_input = SAREK.out.cram_variant_calling.map { meta, cram, crai ->
        [
            meta + [pon_type: 'mutect2', pon_workspace: "${params.outdir}/pon/mutect2/mutect2"],
            cram,
            crai,
        ]
    }

    BAM_CREATE_SOM_PON_GATK(
        ch_mutect2_pon_input,
        fasta,
        fasta_fai.map { meta, fai, _gzi -> [meta, fai, []] },
        dict,
        intervals_bed_combined_for_variant_calling,
        intervals_and_num_intervals,
    )

    emit:
    multiqc_publish = SAREK.out.multiqc_publish
    multiqc_report = SAREK.out.multiqc_report
    mutect2_pon = BAM_CREATE_SOM_PON_GATK.out.pon_vcf
    mutect2_pon_tbi = BAM_CREATE_SOM_PON_GATK.out.pon_index
}
