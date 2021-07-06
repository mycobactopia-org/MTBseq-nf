nextflow.enable.dsl = 2

params.resultsDir = "${params.outdir}/tb_bwa"
params.saveMode = 'copy'
params.shouldPublish = true

process TB_BWA {
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
    touch ${genomeFileName}/bam
    echo "MTBseq --step TBbwa --thread ${task.cpus}"
    """

}

