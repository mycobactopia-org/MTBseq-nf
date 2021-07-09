nextflow.enable.dsl = 2

params.resultsDir = "${params.outdir}/tbbwa"
params.saveMode = 'copy'
params.shouldPublish = true

process TBBWA {
    tag "${genomeFileName}"
    publishDir params.resultsDir, mode: params.saveMode, enabled: params.shouldPublish

    input:
    tuple val(genomeFileName), path("${genomeFileName}_${params.library_name}_R?.fastq.gz")
    path(gatk_jar)
    env USER

    output:
    path("${genomeFileName}/Bam/${genomeFileName}_${params.library_name}*.{bam,bai,bamlog}")
    tuple val(genomeFileName), path("${genomeFileName}/Bam/${genomeFileName}_${params.library_name}*.bam"), emit: bam
    val(genomeFileName), emit:genomes_names

    script:

    """

    gatk-register ${gatk_jar}

    mkdir ${genomeFileName}
    MTBseq --step TBbwa --threads ${task.cpus}
    mv  Bam ./${genomeFileName}/
    """

    stub:

    """
    mkdir ${genomeFileName}
    mkdir ${genomeFileName}/Bam
    mkdir ${genomeFileName}/Bam/${genomeFileName}
    touch ${genomeFileName}/Bam/${genomeFileName}_${params.library_name}.bam
    touch ${genomeFileName}/Bam/${genomeFileName}_${params.library_name}.bai
    touch ${genomeFileName}/Bam/${genomeFileName}_${params.library_name}.bamlog
    echo "MTBseq --step TBbwa --threads ${task.cpus}"
    """

}

