include { DESEQ2_DIFFERENTIAL } from '../../../modules/nf-core/deseq2/differential'
include { GPROFILER2_GOST } from '../../../modules/nf-core/gprofiler2/gost'
include { SHINYNGS_APP } from '../../../modules/nf-core/shinyngs/app'
include { SHINYNGS_STATICEXPLORATORY } from '../../../modules/nf-core/shinyngs/staticexploratory'
include { SHINYNGS_STATICDIFFERENTIAL } from '../../../modules/nf-core/shinyngs/staticdifferential'
include { QUARTO_NOTEBOOK } from '../../../modules/nf-core/quarto/notebook'

workflow RNASEQ_DESEQ2_GPROFILER2 {
    take:
    samplesheet
    counts
    features
    contrasts

    main:
    ch_experiment = samplesheet
        .combine(counts)
        .map { meta, sample_file, _counts_meta, count_file -> [meta, sample_file, count_file] }

    ch_contrasts = channel
        .fromPath(contrasts, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def id = "${row.variable}-${row.reference}-${row.target}"
            [[id: id], row.variable, row.reference, row.target, row.formula ?: '', row.comparison ?: '']
        }

    DESEQ2_DIFFERENTIAL(
        ch_contrasts,
        ch_experiment,
        channel.value([[:], []]),
        channel.value([[:], []]),
    )

    GPROFILER2_GOST(
        DESEQ2_DIFFERENTIAL.out.results,
        channel.value([[:], []]),
        channel.value([[:], []]),
    )

    ch_exploratory_input = ch_experiment
        .combine(features)
        .map { meta, sample_file, count_file, _feature_meta, feature_file ->
            [meta, sample_file, feature_file, [count_file], 'status']
        }

    SHINYNGS_STATICEXPLORATORY(ch_exploratory_input)

    ch_differential_input = DESEQ2_DIFFERENTIAL.out.results
        .combine(ch_experiment.combine(features))
        .map { contrast_meta, result, _experiment_meta, sample_file, count_file, _feature_meta, feature_file ->
            [[id: contrast_meta.id], result, sample_file, feature_file, count_file]
        }

    SHINYNGS_STATICDIFFERENTIAL(
        ch_differential_input.map { meta, result, _sample_file, _feature_file, _count_file -> [meta, result] },
        ch_differential_input.map { _meta, _result, _sample_file, _feature_file, _count_file -> [[:], [], [], []] },
    )

    ch_shiny_app_input = DESEQ2_DIFFERENTIAL.out.results
        .collect()
        .combine(ch_experiment.combine(features))
        .map { results, _experiment_meta, sample_file, count_file, _feature_meta, feature_file ->
            [[id: 'rnaseq'], sample_file, feature_file, [count_file], results.collect { result -> result[1] }]
        }

    SHINYNGS_APP(
        ch_shiny_app_input.map { meta, sample_file, feature_file, assay_files, _results -> [meta, sample_file, feature_file, assay_files] },
        channel.fromPath(contrasts, checkIfExists: true).map { contrast_file -> [[id: 'rnaseq'], contrast_file, []] },
        1,
        channel.value([]),
        ch_shiny_app_input.map { _meta, _sample_file, _feature_file, _assay_files, results -> results },
    )

    QUARTO_NOTEBOOK(
        ch_experiment.map { meta, _sample_file, _count_file -> [meta, file("${projectDir}/assets/rnaseq_report.qmd")] },
        [gprofiler2_organism: 'hsapiens'],
        DESEQ2_DIFFERENTIAL.out.results
            .mix(GPROFILER2_GOST.out.all_enrich)
            .collect(),
        channel.value([]),
    )

    emit:
    deseq2 = DESEQ2_DIFFERENTIAL.out.results
    gprofiler2 = GPROFILER2_GOST.out.all_enrich
    shiny = SHINYNGS_APP.out.app
    exploratory_plots = SHINYNGS_STATICEXPLORATORY.out.boxplots_png
        .mix(SHINYNGS_STATICEXPLORATORY.out.densities_png)
        .mix(SHINYNGS_STATICEXPLORATORY.out.pca2d_png)
        .mix(SHINYNGS_STATICEXPLORATORY.out.pca3d_png)
        .mix(SHINYNGS_STATICEXPLORATORY.out.mad_png)
        .mix(SHINYNGS_STATICEXPLORATORY.out.dendro)
    differential_plots = SHINYNGS_STATICDIFFERENTIAL.out.volcanos_png
    report = QUARTO_NOTEBOOK.out.html
}
