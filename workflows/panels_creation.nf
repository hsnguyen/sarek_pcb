include { BAM_CREATE_SOM_PON_GATK } from '../subworkflows/nf-core/bam_create_som_pon_gatk'
include { CNVKIT_BATCH } from '../modules/nf-core/cnvkit/batch'
include { PURECN_COVERAGE } from '../modules/nf-core/purecn/coverage'
include { PURECN_INTERVALFILE } from '../modules/nf-core/purecn/intervalfile'
include { PURECN_NORMALDB } from '../modules/nf-core/purecn/normaldb'
//include { SAMTOOLS_VIEW } from '../modules/nf-core/samtools/view'
// CRAM_TO_BAM conversion
include { SAMTOOLS_CONVERT as CRAM_TO_BAM } from '../modules/nf-core/samtools/convert'
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
    ch_normal_samples_for_panels = input_sample.filter { sample -> sample[0].status == 0 }
    //call sarek but just for reads pre-processing and alignment here
    SAREK(
        ch_normal_samples_for_panels,
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
    ch_panel_alignments_cram = SAREK.out.cram_variant_calling
    //1. Build PoN for GATK mutect2
    BAM_CREATE_SOM_PON_GATK(
        ch_panel_alignments_cram,
        fasta,
        fasta_fai.map { meta, fai -> [meta, fai, []] },
        dict,
        'pcb_mutect2',
        intervals_bed_combined_for_variant_calling,
        intervals_and_num_intervals,
    )

    //2. Build PoN for cnvkit: need BAM input!!!
    // Below is piece of code copied from sarek.nf since don't want to change it:
    // Differentiate between bam and cram files
    cram_variant_calling_status_tmp = ch_panel_alignments_cram.branch { _meta, file, _index ->
        bam: file.toString().endsWith('.bam')
        cram: file.toString().endsWith('.cram')
    }

    // convert cram files
    CRAM_TO_BAM(cram_variant_calling_status_tmp.cram, fasta.combine(fasta_fai).map { meta, fasta_, _meta_fai, fai -> [meta, fasta_, fai] }.collect())

    // gather all bam files
    ch_panel_alignments_bam = CRAM_TO_BAM.out.bam
        .join(CRAM_TO_BAM.out.bai, by: [0])
        .mix(cram_variant_calling_status_tmp.bam)
        .map { meta, bam, bai ->
            [meta + [data_type: 'bam'], bam, bai]
    }

    CNVKIT_BATCH(
        ch_panel_alignments_bam.map { _meta, bam, _bai -> [[id: 'cnvkit'], bam] }.groupTuple().map { meta, bams -> [meta, [], [], bams, []] },
        fasta.map { meta, fasta_file -> [meta, fasta_file, []] },
        intervals_bed_combined.map { bed -> [[id: 'cnvkit'], bed] },
        [[:], []],
        true,
      )

    //3. Build PoN for purecn
    PURECN_INTERVALFILE(
        intervals_bed_combined.map { bed -> [[id: 'purecn_intervals'], bed] },
        fasta,
        params.purecn_genome ?: 'hg38',
    )
    
    PURECN_COVERAGE(
        //ch_panel_alignments_bam.map { _meta, bam, bai -> [[id: _meta.id], bam, bai] },
        ch_panel_alignments_bam,
        PURECN_INTERVALFILE.out.txt.collect{ it[1] },
    )

    ch_purecn_normaldb = PURECN_COVERAGE.out.txt.collect{ it[1] }.map { coverage_files -> [[id: 'purecn_normaldb'], coverage_files] }
        .combine(BAM_CREATE_SOM_PON_GATK.out.genomicsdb.collect { _meta, gdb -> [gdb] }.toList())
        .map { meta, cov, gdb -> [meta, cov, gdb, []]}
        //.combine(BAM_CREATE_SOM_PON_GATK.out.pon_vcf.collect {_meta, vcf -> [vcf] }.toList())
        //.combine(BAM_CREATE_SOM_PON_GATK.out.pon_index.collect {_meta, tbi -> [tbi] }.toList())
        //.map { meta, cov, vcf, tbi -> [meta, cov, vcf, tbi] }
    PURECN_NORMALDB(
        ch_purecn_normaldb,
        params.purecn_genome ?: 'hg38',
        params.purecn_assay ?: 'pcb_wes',
    )

    emit:
    multiqc_publish = SAREK.out.multiqc_publish
    multiqc_report = SAREK.out.multiqc_report
//    cnvkit_pon = CNVKIT_BATCH.out.cnn
//    mutect2_pon = BAM_CREATE_SOM_PON_GATK.out.pon_vcf
//    mutect2_pon_tbi = BAM_CREATE_SOM_PON_GATK.out.pon_index
//   purecn_normaldb = PURECN_NORMALDB.out.rds
}
