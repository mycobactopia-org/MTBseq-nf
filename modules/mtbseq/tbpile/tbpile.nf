nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/tbpile"
params.save_mode = 'copy'
params.should_publish = true

process TBPILE {
    tag "${genomeFileName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish
    stageInMode 'copy'

    input:
    tuple val(genomeFileName), path("GATK_Bam/${genomeFileName}_${params.library_name}_*gatk.bam")
    path(gatk_jar)
    env(USER)

    output:
    path("${genomeFileName}/Mpileup/${genomeFileName}_${params.library_name}*.gatk.{mpileup,mpileuplog}")
    tuple val(genomeFileName), path("${genomeFileName}/Mpileup/${genomeFileName}_${params.library_name}*.gatk.mpileup"), emit: mpileup

    script:

    """

    gatk-register ${gatk_jar}

    mkdir Mpileup

    MTBseq --step TBpile \
    --threads ${task.cpus} \
    1>>.command.out \
    2>>.command.err \
    || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

    mkdir ${genomeFileName}
    mv  Mpileup ./${genomeFileName}/
    """

    stub:

    """
    echo "MTBseq --step TBpile --threads ${task.cpus}"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir ${genomeFileName}
    mkdir ${genomeFileName}/Mpileup
    touch ${genomeFileName}/Mpileup/${genomeFileName}_${params.library_name}.gatk.mpileup
    touch ${genomeFileName}/Mpileup/${genomeFileName}_${params.library_name}.gatk.mpileuplog

    """

}
