process MTBSEQ_TBFULL {
    tag "${params.project}"
    label 'process_medium'

    conda "bioconda::mtbseq=1.0.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mtbseq:1.0.4--hdfd78af_2':
        'quay.io/biocontainers/mtbseq:1.0.4--hdfd78af_2' }"

    input:
    path(*)
    env(USER)
    tuple path(ref_resistance_list),path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)

    output:
    path("Called")
    path("Position_Tables")
    path("Classification")
    path("Statistics")
    path("Called/*_${params.library_name}*gatk_position_variants*.tab"), emit: position_variants
    path("Position_Tables/*_${params.library_name}*.gatk_position_table.tab"), emit: position_tables

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    MTBseq --step TBfull \
       --thread ${task.cpus} \
       --project ${params.project} \
       --minbqual ${params.minbqual} \
       --mincovf ${params.mincovf} \
       --mincovr ${params.mincovr} \
       --minphred ${params.minphred} \
       --minfreq ${params.minfreq} \
       --resilist ${ref_resistance_list} \
       --unambig ${params.unambig} \
       --window ${params.window} \
       --distance ${params.distance} \
       --intregions ${ref_interesting_regions} \
       --categories ${ref_gene_categories} \
       --basecalib ${ref_base_quality_recalibration} \
       $args \
        1>>.command.out \
        2>>.command.err \
        || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mtbseq: \$(echo \$(MTBseq --version 2>&1) | sed 's/^.*MTBseq //' ))
    END_VERSIONS
    """
}
