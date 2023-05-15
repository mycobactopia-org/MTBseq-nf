process MTBSEQ_TBBWA {
    tag "$meta.id"
    label 'process_medium'
    conda "bioconda::mtbseq=1.0.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mtbseq:1.0.4--hdfd78af_2':
        'quay.io/biocontainers/mtbseq:1.0.4--hdfd78af_2' }"

    input:
    tuple val(meta), path("${meta.id}_${params.library_name}_R?.fastq.gz")
    env(USER)
    tuple path(ref_resistance_list), path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)

    output:
    path("Bam/${meta.id}_${params.library_name}*.{bam,bai,bamlog}")
    tuple val(meta), path("Bam/${genomeFileName}_${params.library_name}*.{bam,bai}"), emit: bam_tuple
    path("Bam/${meta.id}_${params.library_name}*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir Bam

    MTBseq --step TBbwa \
            --threads ${task.cpus} \
            --project ${params.project} \
            --resilist ${ref_resistance_list} \
            --intregions ${ref_interesting_regions} \
            --categories ${ref_gene_categories} \
            --basecalib ${ref_base_quality_recalibration} \
        1>>.command.out \
        2>>.command.err \
        || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mtbseq: \$(echo \$(MTBseq --version 2>&1) | sed 's/^.*MTBseq //' ))
    END_VERSIONS
    """
    """
}
