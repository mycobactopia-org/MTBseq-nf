nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/tbvariants"
params.save_mode = 'copy'
params.should_publish = true
params.mincovf = 4
params.mincovr = 4
params.minphred = 4
params.minfreq = 75

process TBVARIANTS {
    tag "${genomeFileName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(genomeFileName), path("Position_Tables/${genomeFileName}_${params.library_name}*.gatk_position_table.tab")
    path(gatk_jar)
    env(USER)

    output:
    path("${genomeFileName}/Called/${genomeFileName}_${params.library_name}*gatk_position_{uncovered,variants}*.tab")
    path("${genomeFileName}/Called/${genomeFileName}_${params.library_name}*gatk_position_variants*.tab"), emit: tbjoin_input

    script:

    """

    gatk-register ${gatk_jar}

    MTBseq --step TBvariants \
    --threads ${task.cpus} \
    --mincovf ${params.mincovf} \
    --mincovr ${params.mincovr} \
    --minphred ${params.minphred} \
    --minfreq ${params.minfreq} \
    1>>.command.out \
    2>>.command.err \
    || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

    mkdir ${genomeFileName}
    mv  Called ./${genomeFileName}
    """

    stub:

    """
    echo "MTBseq --step TBvariants --mincovf ${params.mincovf} \
        --mincovr ${params.mincovr} \
        --minphred ${params.minphred} \
        --minfreq ${params.minfreq}"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir ${genomeFileName}
    mkdir ${genomeFileName}/Called
    touch ${genomeFileName}/Called/${genomeFileName}_${params.library_name}.gatk_position_uncovered_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_outmode000.tab
    touch ${genomeFileName}/Called/${genomeFileName}_${params.library_name}.gatk_position_variants_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_outmode000.tab

    """

}
