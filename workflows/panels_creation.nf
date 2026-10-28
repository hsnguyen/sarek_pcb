include { BAM_CREATE_SOM_PON_GATK } from '../subworkflows/nf-core/bam_create_som_pon_gatk'
include { CNVKIT_BATCH } from '../modules/nf-core/cnvkit/batch'
include { PURECN_COVERAGE } from '../modules/nf-core/purecn/coverage'
include { PURECN_INTERVALFILE } from '../modules/nf-core/purecn/intervalfile'
include { PURECN_NORMALDB } from '../modules/nf-core/purecn/normaldb'
include { SAMTOOLS_VIEW } from '../modules/nf-core/samtools/view'
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
    ch_normal_samples_for_panels = input_sample
        .filter { meta, _input, _input_index ->
            meta.status == 0 && meta.pon_db
        }
        .ifEmpty {
            error("No normal samples request a panel of normals. Set 'pon_db' for at least one normal sample.")
        }

    ch_panel_requests_by_id = ch_normal_samples_for_panels
        .map { meta, _input, _input_index -> [meta.id, meta] }

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

    ch_panel_alignments = SAREK.out.cram_variant_calling
        .map { meta, alignment, index -> [meta.id, meta, alignment, index] }
        .join(ch_panel_requests_by_id, failOnDuplicate: true, failOnMismatch: true)
        .map { _id, meta, alignment, index, request_meta ->
            [meta + request_meta.subMap('pon_db'), alignment, index]
        }

    ch_panel_alignments_by_format = ch_panel_alignments.branch { _meta, alignment, _index ->
        bam: alignment.toString().endsWith('.bam')
        cram: alignment.toString().endsWith('.cram')
    }

    SAMTOOLS_VIEW(
        ch_panel_alignments_by_format.cram,
        fasta.combine(fasta_fai).map { meta, fasta_file, _meta_fai, fai -> [meta, fasta_file, fai] }.collect(),
        [[:], []],
        [[:], []],
        'bai',
    )

    ch_panel_bam = ch_panel_alignments_by_format.bam.mix(SAMTOOLS_VIEW.out.bam.join(SAMTOOLS_VIEW.out.bai, failOnDuplicate: true, failOnMismatch: true))

    ch_mutect2_pon_input = ch_panel_alignments
        .filter { meta, _cram, _crai -> hasPonType(meta, 'mutect2') }
        .map { meta, cram, crai ->
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

    CNVKIT_BATCH(
        ch_panel_bam.filter { meta, _bam, _bai -> hasPonType(meta, 'cnvkit') }.map { _meta, bam, _bai -> [[id: 'cnvkit'], bam] }.groupTuple().map { meta, bams -> [meta, [], [], bams, []] },
        fasta.map { meta, fasta_file -> [meta, fasta_file, []] },
        intervals_bed_combined.map { bed -> [[id: 'cnvkit'], bed] },
        [[:], []],
        true,
    )

    PURECN_INTERVALFILE(
        intervals_bed_combined.map { bed -> [[id: 'purecn'], bed] },
        fasta,
        params.purecn_genome ?: params.ascat_genome,
    )

    PURECN_COVERAGE(
        ch_panel_bam.filter { meta, _bam, _bai -> hasPonType(meta, 'purecn') },
        PURECN_INTERVALFILE.out.txt.collect(),
    )

    PURECN_NORMALDB(
        PURECN_COVERAGE.out.txt.collect().map { coverage_files -> [[id: 'purecn'], coverage_files, [], []] },
        params.purecn_genome ?: params.ascat_genome,
        params.purecn_assay,
    )

    emit:
    multiqc_publish = SAREK.out.multiqc_publish
    multiqc_report = SAREK.out.multiqc_report
    cnvkit_pon = CNVKIT_BATCH.out.cnn
    mutect2_pon = BAM_CREATE_SOM_PON_GATK.out.pon_vcf
    mutect2_pon_tbi = BAM_CREATE_SOM_PON_GATK.out.pon_index
    purecn_normaldb = PURECN_NORMALDB.out.rds
}

def hasPonType(meta, ponType) {
    meta.pon_db
        ?.tokenize(';')
        ?.any { requestedPonType ->
            requestedPonType.equalsIgnoreCase(ponType)
        }
}
