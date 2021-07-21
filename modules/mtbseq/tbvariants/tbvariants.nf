nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar


params.resultsDir = "${params.outdir}/tbvariants"
params.saveMode = 'copy'
params.shouldPublish = true
params.mincovf = 4
params.mincovr = 4
params.minphred = 4
params.minfreq = 75

process TBVARIANTS {
    tag "${genomeFileName}"
    publishDir params.resultsDir, mode: params.saveMode, enabled: params.shouldPublish

    input:
    tuple val(genomeFileName), path("Position_Tables/${genomeFileName}_${params.library_name}*.gatk_position_table.tab")
    path(gatk_jar)
    env USER

    output:
    path("${genomeFileName}/Called/${genomeFileName}_${params.library_name}*gatk_position_{uncovered,variants}*.tab")
    path("${genomeFileName}/Called/${genomeFileName}_${params.library_name}*gatk_position_variants*.tab"), emit: tbjoin_input

    script:

    """
    echo "MTBseq --step TBvariants --mincovf ${params.mincovf} \
        --mincovr ${params.mincovr} \
        --minphred ${params.minphred} \
        --minfreq ${params.minfreq}"

    gatk-register ${gatk_jar}

    mkdir ${genomeFileName}
    MTBseq --step TBvariants \
        --mincovf ${params.mincovf} \
        --mincovr ${params.mincovr} \
        --minphred ${params.minphred} \
        --minfreq ${params.minfreq}
    mv  Called ./${genomeFileName}
    """

    stub:

    """
    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir ${genomeFileName}
    mkdir ${genomeFileName}/Called
    touch ${genomeFileName}/Called/${genomeFileName}_${params.library_name}.gatk_position_uncovered_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_outmode000.tab
    touch ${genomeFileName}/Called/${genomeFileName}_${params.library_name}.gatk_position_variants_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_outmode000.tab

    """

}
