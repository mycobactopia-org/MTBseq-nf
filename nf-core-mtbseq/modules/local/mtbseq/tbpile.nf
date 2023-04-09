process MTBSEQ_TBPILE {
    tag "$meta.id - $params.project"
    label 'process_medium'
    conda "bioconda::mtbseq=1.0.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mtbseq:1.0.4--hdfd78af_2':
        'quay.io/biocontainers/mtbseq:1.0.4--hdfd78af_2' }"

    input:
    tuple val(meta), path("GATK_BAM/*")
    env(USER)
    tuple path(ref_resistance_list), path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)

    output:
    path("Mpileup/${genomeFileName}_${params.library_name}*.gatk.{mpileup,mpileuplog}")
    tuple val(meta), path("Mpileup/${genomeFileName}_${params.library_name}*.gatk.mpileup"), emit: mpileup
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir Mpileup
    MTBseq --step TBpile \\
    --threads ${task.cpus} \\
    --project ${params.project} \\
    --resilist ${ref_resistance_list} \\
    --intregions ${ref_interesting_regions} \\
    --categories ${ref_gene_categories} \\
    --basecalib ${ref_base_quality_recalibration} \\
    $args \\
    1>>.command.out \\
    2>>.command.err \\
    || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mtbseq: \$(echo \$(MTBseq --version 2>&1) | sed 's/^.*MTBseq //' ))
    END_VERSIONS
    """
}
