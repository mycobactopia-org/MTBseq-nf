nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/tbstats"
params.save_mode = 'copy'
params.should_publish = true

process TBSTATS {
    tag "${genomeFileName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(genomeFileName), path("Bam/*")
    path(gatk_jar)
    env(USER)

    output:
    path("${genomeFileName}/Statistics/Mapping_and_Variant_Statistics.tab")

    script:

    """
    gatk-register ${gatk_jar}

    mkdir Statistics

    MTBseq --step TBstats \
    --threads ${task.cpus} \
    1>>.command.out \
    2>>.command.err \
    || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq


    mkdir ${genomeFileName}
    mv  Statistics ./${genomeFileName}/
    """

    stub:

    """
    echo "MTBseq --step TBstats"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir ${genomeFileName}
    mkdir ${genomeFileName}/Statistics
    touch ${genomeFileName}/Statistics/Mapping_and_Variant_Statistics.tab

    """

}
