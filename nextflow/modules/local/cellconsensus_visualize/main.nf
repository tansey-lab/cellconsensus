process CELLCONSENSUS_VISUALIZE {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jeffquinnmsk/cellconsensus:latest' :
        'docker.io/jeffquinnmsk/cellconsensus:latest' }"

    input:
    tuple val(meta), path(annotate_dir)

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_visualize"
    def args = task.ext.args ?: ''
    def spatial_key_arg = params.visualize_spatial_key ? "--spatial-key ${params.visualize_spatial_key}" : ''
    def point_size_arg = params.visualize_point_size ? "--point-size ${params.visualize_point_size}" : ''
    """
    mkdir -p ${prefix}

    cellconsensus-visualize \\
        --output ${prefix}/plots.pdf \\
        ${spatial_key_arg} \\
        ${point_size_arg} \\
        ${args} \\
        ${annotate_dir}/annotated.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellconsensus: \$(python -c "import cellconsensus; print(cellconsensus.__version__)")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_visualize"
    """
    mkdir -p ${prefix}
    touch ${prefix}/plots.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellconsensus: 1.2.0
    END_VERSIONS
    """
}
