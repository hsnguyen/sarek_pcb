//
// Run GATK mutect2, genomicsdbimport and createsomaticpanelofnormals
//

include { GATK4_CREATESOMATICPANELOFNORMALS } from '../../../modules/nf-core/gatk4/createsomaticpanelofnormals'
include { GATK4_GENOMICSDBIMPORT            } from '../../../modules/nf-core/gatk4/genomicsdbimport'
include { GATK4_MERGEMUTECTSTATS            } from '../../../modules/nf-core/gatk4/mergemutectstats'
include { GATK4_MERGEVCFS                   } from '../../../modules/nf-core/gatk4/mergevcfs'
include { GATK4_MUTECT2                     } from '../../../modules/nf-core/gatk4/mutect2'

workflow BAM_CREATE_SOM_PON_GATK {
    take:
    ch_input // channel: [ val(meta), path(input), path(input_index) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_fai // channel: [ val(meta), path(fai), path(gzi) ]
    ch_dict // channel: [ val(meta), path(dict) ]
    ch_intervals_gendb // channel: [ path(interval_file) ]
    ch_intervals_num // channel: [ path(intervals), val(num_intervals) ] or [ [], 0 ] if no intervals

    main:
    // Combine input and intervals for scatter and gather strategy
    // Move num_intervals to meta map for GATK4_MUTECT2
    ch_input_intervals = ch_input
        .combine(ch_intervals_num)
        .map { meta, input, input_index, intervals_, num_intervals ->
            [meta + [num_intervals: num_intervals], input, input_index, intervals_]
        }

    // Perform variant calling for each sample using mutect2 module in panel of normals mode
    GATK4_MUTECT2(
        ch_input_intervals,
        ch_fasta,
        ch_fai,
        ch_dict,
        [],
        [],
        [],
        [],
    )

    // Branch outputs based on whether intervals were used
    vcf_branch = GATK4_MUTECT2.out.vcf.branch { meta, _vcf ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    tbi_branch = GATK4_MUTECT2.out.tbi.branch { meta, _tbi ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    stats_branch = GATK4_MUTECT2.out.stats.branch { meta, _stats ->
        intervals: meta.num_intervals > 1
        no_intervals: meta.num_intervals <= 1
    }

    // Only when using intervals: group outputs by sample and merge
    GATK4_MERGEMUTECTSTATS(stats_branch.intervals.map { meta, stats -> [groupKey(meta, meta.num_intervals), stats] }.groupTuple())
    GATK4_MERGEVCFS(
        vcf_branch.intervals.map { meta, vcf -> [groupKey(meta, meta.num_intervals), vcf] }.groupTuple(),
        ch_dict,
    )

    // Mix merged and non-interval outputs, remove num_intervals from meta
    ch_vcf = vcf_branch.no_intervals
        .mix(GATK4_MERGEVCFS.out.vcf)
        .map { meta, vcf -> [meta - meta.subMap('num_intervals'), vcf] }

    ch_tbi = tbi_branch.no_intervals
        .mix(GATK4_MERGEVCFS.out.tbi)
        .map { meta, tbi -> [meta - meta.subMap('num_intervals'), tbi] }

    ch_stats = stats_branch.no_intervals
        .mix(GATK4_MERGEMUTECTSTATS.out.stats)
        .map { meta, stats -> [meta - meta.subMap('num_intervals'), stats] }

    ch_gendb_input = ch_vcf
        .join(ch_tbi, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf, tbi -> [meta.pon_db, vcf, tbi] }
        .groupTuple()
        .map { pon_db, vcfs, tbis ->
            java.nio.file.Path pon_db_path = file(pon_db).toAbsolutePath()
            [
                [
                    id: pon_db_path.fileName.toString(),
                    pon_db: pon_db_path.toString(),
                    pon_db_exists: pon_db_path.exists(),
                    pon_db_parent: pon_db_path.parent.toString(),
                ],
                vcfs,
                tbis,
            ]
        }
        .combine(ch_intervals_gendb)
        .combine(ch_dict.map { _meta, dict -> [dict] })
        .map { meta, vcf, tbi, interval, dict ->
            [
                meta,
                vcf,
                tbi,
                interval,
                meta.pon_db_exists ? file(meta.pon_db, checkIfExists: true) : [],
                dict,
            ]
        }

    // Create a new workspace or append to the existing workspace for each pon_db.
    GATK4_GENOMICSDBIMPORT(
        ch_gendb_input,
        false,
        ch_gendb_input.map { meta, _vcf, _tbi, _interval, _workspace, _dict -> meta.pon_db_exists },
        false,
    )

    // Create PON from the genomicsdb workspace using createsomaticpanelofnormals
    genomicsdb = GATK4_GENOMICSDBIMPORT.out.genomicsdb.mix(GATK4_GENOMICSDBIMPORT.out.updatedb)
    GATK4_CREATESOMATICPANELOFNORMALS(genomicsdb, ch_fasta, ch_fai.map { meta, fai, _gzi -> [meta, fai] }, ch_dict)

    versions = channel.empty()
        .mix(
            GATK4_MUTECT2.out.versions,
            GATK4_MERGEMUTECTSTATS.out.versions,
            GATK4_GENOMICSDBIMPORT.out.versions,
        )

    emit:
    genomicsdb // channel: [ val(meta), path(genomicsdb) ]
    mutect2_index = ch_tbi // channel: [ val(meta), path(tbi) ]
    mutect2_stats = ch_stats // channel: [ val(meta), path(stats) ]
    mutect2_vcf   = ch_vcf // channel: [ val(meta), path(vcf) ]
    pon_index     = GATK4_CREATESOMATICPANELOFNORMALS.out.tbi // channel: [ val(meta), path(tbi) ]
    pon_vcf       = GATK4_CREATESOMATICPANELOFNORMALS.out.vcf // channel: [ val(meta), path(vcf) ]
    versions      // channel: [ versions.yml ]
}
