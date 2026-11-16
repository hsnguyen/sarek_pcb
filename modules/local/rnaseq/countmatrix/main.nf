process RNASEQ_COUNTMATRIX {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), val(sample_meta), path(count_files)

    output:
    tuple val(meta), path('rnaseq_samplesheet.tsv'), emit: samplesheet
    tuple val(meta), path('rnaseq_counts.tsv'), emit: counts
    tuple val(meta), path('rnaseq_features.tsv'), emit: features

    script:
    def samples = sample_meta.collect { item ->
        "${item.id}\t${item.patient}\t${item.status}\t${item.assay}\t${item.strandness}\t${item.seq_platform}"
    }.join('\n')
    def count_columns = sample_meta.collect { item ->
        item.strandness == 'forward' ? 3 : item.strandness == 'reverse' ? 4 : 2
    }.join('\n')
    def count_file_paths = count_files.collect { item -> item.toString() }.join('\n')
    def sample_ids = sample_meta.collect { item -> item.id }.join('\t')
    """
    cat <<-END_SAMPLESHEET > rnaseq_samplesheet.tsv
    sample	patient	status	assay	strandness	seq_platform
    ${samples}
    END_SAMPLESHEET

    cat <<-END_COLUMNS > columns.tsv
    ${count_columns}
    END_COLUMNS
    cat <<-END_COUNT_FILES > count_files.tsv
    ${count_file_paths}
    END_COUNT_FILES

    printf 'gene_id\\t%s\\n' '${sample_ids}' > rnaseq_counts.tsv
    awk 'NR > 4 { print \$1 }' ${count_files[0]} > gene_ids.tsv
    index=0
    while IFS= read -r column && IFS= read -r count_file <&3; do
        awk -v column="\$column" 'NR > 4 { print \$column }' "\$count_file" > "counts_\${index}.tsv"
        index=\$((index + 1))
    done < columns.tsv 3< count_files.tsv
    paste gene_ids.tsv counts_*.tsv >> rnaseq_counts.tsv
    awk 'BEGIN { OFS="\\t"; print "gene_id", "gene_name" } NR > 4 { print \$1, \$1 }' ${count_files[0]} > rnaseq_features.tsv
    """

    stub:
    """
    printf 'sample\\tpatient\\tstatus\\tassay\\tstrandness\\tseq_platform\\n' > rnaseq_samplesheet.tsv
    printf 'gene_id\\n' > rnaseq_counts.tsv
    printf 'gene_id\\tgene_name\\n' > rnaseq_features.tsv
    """
}
