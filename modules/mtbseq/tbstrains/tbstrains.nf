nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar


params.results_dir = "${params.outdir}/tbstrains"
params.save_mode = 'copy'
params.should_publish = true
params.mincovf = 4
params.mincovr = 4
params.minphred = 4
params.minfreq = 75


// TODO: Add the tbjoin workflow
process TBSTRAINS {
    tag "${genomeFileName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(genomeFileName), path("Position_Tables/${genomeFileName}_${params.library_name}*.gatk_position_table.tab")
    path(gatk_jar)
    env USER

    output:
    path("${genomeFileName}/Classification/Strain_Classification.tab")

    script:

    """

    gatk-register ${gatk_jar}

    mkdir ${genomeFileName}
    MTBseq --step TBstrains --mincovf ${params.mincovf} \
        --mincovr ${params.mincovr} \
        --minphred ${params.minphred} \
        --minfreq ${params.minfreq} \
        2>${task.process}_${genomeFileName}_err.log 1>${task.process}_${genomeFileName}_out.log
    mv  Classification ./${genomeFileName}/
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
