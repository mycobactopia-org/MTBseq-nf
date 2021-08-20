nextflow.enable.dsl = 2

process TBREFINE {
    tag "${genomeFileName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(genomeFileName), path("Bam/${genomeFileName}_${params.library_name}*.bam")
    path(gatk_jar)
    env USER

    output:
    path("${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}*.gatk.{bam,bai,bamlog,grp,intervals}")
    tuple val(genomeFileName), path("${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}*gatk.bam"), emit: gatk_bam

    script:

    """

    gatk-register ${gatk_jar}

    mkdir ${genomeFileName}
    MTBseq --step TBrefine --threads ${task.cpus} 2>${task.process}_${genomeFileName}_err.log 1>${task.process}_${genomeFileName}_out.log
    mv  GATK_Bam ./${genomeFileName}/
    """

    stub:

    """
    echo "MTBseq --step TBrefine --threads ${task.cpus}"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir ${genomeFileName}
    mkdir ${genomeFileName}/GATK_Bam
    touch ${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}.gatk.bam
    touch ${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}.gatk.bai
    touch ${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}.gatk.bamlog
    touch ${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}.gatk.grp
    touch ${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}.gatk.intervals
    """
}
