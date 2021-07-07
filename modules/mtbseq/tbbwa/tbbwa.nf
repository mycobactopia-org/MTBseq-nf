nextflow.enable.dsl = 2

params.resultsDir = "${params.outdir}/tbbwa"
params.saveMode = 'copy'
params.shouldPublish = true

process TBBWA {
    tag "${genomeFileName}"
    publishDir params.resultsDir, mode: params.saveMode, enabled: params.shouldPublish
    container 'quay.io/biocontainers/mtbseq:1.0.3--pl526_1'
    cpus 8
    memory "15 GB"

    input:
    tuple val(genomeFileName), path("${genomeFileName}_${params.library_name}_R?.fastq.gz")
    path(gatk_jar)
    env USER

    output:
    path("${genomeFileName}/Bam/${genomeFileName}_${params.library_name}*.{bam,bai,bamlog}")

    script:

    """

    gatk-register ${gatk_jar}

    mkdir ${genomeFileName}
    MTBseq --step TBbwa --thread ${task.cpus}
    mv  Bam ./${genomeFileName}/
    """

    stub:

    """
    mkdir ${genomeFileName}
    mkdir ${genomeFileName}/Bam
    mkdir ${genomeFileName}/Bam/${genomeFileName}
    touch ${genomeFileName}/Bam/${genomeFileName}/${genomeFileName}_${params.library_name}.bam
    touch ${genomeFileName}/Bam/${genomeFileName}/${genomeFileName}_${params.library_name}.bai
    touch ${genomeFileName}/Bam/${genomeFileName}/${genomeFileName}_${params.library_name}.bamlog
    echo "MTBseq --step TBbwa --thread ${task.cpus}"
    """

}

