process CELLCONSENSUS_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jeffquinnmsk/cellconsensus:latest' :
        'docker.io/jeffquinnmsk/cellconsensus:latest' }"

    input:
    tuple val(meta), path(adata)

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_annotate"
    def args = task.ext.args ?: ''
    def cluster_key_arg = params.cluster_key ? "--cluster-key ${params.cluster_key}" : ''
    def cancer_arg = params.include_cancer ? "--include-cancer" : ''
    def cancer_types_arg = params.cancer_types ? "--cancer-types ${params.cancer_types}" : ''
    def save_model_arg = params.save_model ? "--save-model ${prefix}/model.pkl" : ''
    """
    mkdir -p ${prefix}

    cellconsensus-annotate \\
        --output ${prefix}/predictions.csv \\
        --output-h5ad ${prefix}/annotated.h5ad \\
        --clustering ${params.clustering} \\
        --levels ${params.levels} \\
        --output-format ${params.output_format} \\
        --n-neighbors ${params.n_neighbors} \\
        --n-neighbors-lvl2 ${params.n_neighbors_lvl2} \\
        --n-neighbors-lvl3 ${params.n_neighbors_lvl3} \\
        --n-smooth ${params.n_smooth} \\
        --ref-top-k ${params.ref_top_k} \\
        --graph-level ${params.graph_level} \\
        ${cluster_key_arg} \\
        ${cancer_arg} \\
        ${cancer_types_arg} \\
        ${save_model_arg} \\
        ${args} \\
        ${adata}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellconsensus: \$(python -c "import cellconsensus; print(cellconsensus.__version__)")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_annotate"
    """
    mkdir -p ${prefix}
    touch ${prefix}/predictions.csv
    touch ${prefix}/annotated.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellconsensus: 1.2.0
    END_VERSIONS
    """
}
