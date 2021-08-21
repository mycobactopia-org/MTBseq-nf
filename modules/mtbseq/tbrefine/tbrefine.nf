nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/tbrefine"
params.save_mode = 'copy'
params.should_publish = true

process TBREFINE {
    tag "${genomeFileName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(genomeFileName), path("Bam/${genomeFileName}_${params.library_name}*.bam")
    path(gatk_jar)
    env(USER)

    output:
    tuple val(genomeFileName), path("${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}*gatk.{bam,bai,bamlog,grp,intervals}"), emit: gatk_bam

    script:

    """

    gatk-register ${gatk_jar}

    mkdir GATK_Bam

    MTBseq --step TBrefine \
    --threads ${task.cpus} \
    1>>.command.out \
    2>>.command.err \
    || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

    mkdir ${genomeFileName}
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
