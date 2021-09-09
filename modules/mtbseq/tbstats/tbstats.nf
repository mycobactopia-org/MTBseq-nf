nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/tbstats"
params.save_mode = 'copy'
params.should_publish = true

process TBSTATS {
    tag "${params.project_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    path("Bam/*")
    path("Position_Tables/*")
    path(gatk_jar)
    env(USER)

    output:
    path("Statistics/Mapping_and_Variant_Statistics.tab")

    script:

    """
    gatk-register ${gatk_jar}

    mkdir Statistics

    MTBseq --step TBstats \
    --threads ${task.cpus} \
    1>>.command.out \
    2>>.command.err \
    || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

    """

    stub:

    """
    echo "MTBseq --step TBstats --threads ${task.cpus}"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir Statistics
    touch Statistics/Mapping_and_Variant_Statistics.tab
    """

}
