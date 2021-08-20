nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/tbbwa"
params.save_mode = 'copy'
params.should_publish = true

process TBBWA {
    tag "${genomeFileName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

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


    MTBseq --step TBbwa \
    --threads ${task.cpus} \
    1>>.command.out \
    2>>.command.err \
    || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq


    mkdir ${genomeFileName}
    mv  Bam ./${genomeFileName}/
    """

    stub:

    """
    echo "MTBseq --step TBbwa --threads ${task.cpus}"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    touch ${task.process}_${genomeFileName}_out.log
    touch ${task.process}_${genomeFileName}_err.log

    mkdir ${genomeFileName}
    mkdir ${genomeFileName}/Bam
    mkdir ${genomeFileName}/Bam/${genomeFileName}
    touch ${genomeFileName}/Bam/${genomeFileName}_${params.library_name}.bam
    touch ${genomeFileName}/Bam/${genomeFileName}_${params.library_name}.bai
    touch ${genomeFileName}/Bam/${genomeFileName}_${params.library_name}.bamlog

    """

}

