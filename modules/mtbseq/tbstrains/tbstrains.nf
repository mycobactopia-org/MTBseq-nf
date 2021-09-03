nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/tbstrains"
params.save_mode = 'copy'
params.should_publish = true
params.mincovf = 4
params.mincovr = 4
params.minphred = 4
params.minfreq = 75


process TBSTRAINS {
    tag "${params.project_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    path("Position_Tables/*")
    path(gatk_jar)
    env(USER)

    output:
    path("Classification/Strain_Classification.tab")

    script:

    """

    gatk-register ${gatk_jar}

    mkdir Classification

    MTBseq --step TBstrains \
    --threads ${task.cpus} \
    --mincovf ${params.mincovf} \
    --mincovr ${params.mincovr} \
    --minphred ${params.minphred} \
    --minfreq ${params.minfreq} \
    1>>.command.out \
    2>>.command.err \
    || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq


    """

    stub:

    """
    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    echo "MTBseq --step TBstrains --mincovf ${params.mincovf} \
        --mincovr ${params.mincovr} \
        --minphred ${params.minphred} \
        --minfreq ${params.minfreq}"

    mkdir ${genomeFileName}
    mkdir ${genomeFileName}/Classification
    touch ${genomeFileName}/Classification/Strain_Classification.tab

    """

}
