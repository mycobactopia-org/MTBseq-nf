nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar


params.resultsDir = "${params.outdir}/tblist"
params.saveMode = 'copy'
params.shouldPublish = true
params.minbqual = 13

// TODO: Add the tbjoin workflow
process TBLIST {
    tag "${genomeFileName}"
    publishDir params.resultsDir, mode: params.saveMode, enabled: params.shouldPublish

    input:
    tuple val(genomeFileName), path("Mpileup/${genomeFileName}_${params.library_name}*.gatk.mpileup")
    path(gatk_jar)
    env USER

    output:
    path("${genomeFileName}/Position_Tables/${genomeFileName}_${params.library_name}*.gatk_position_table.tab"), emit: tbjoin_input
    tuple val(genomeFileName), path("${genomeFileName}/Position_Tables/${genomeFileName}_${params.library_name}*.gatk_position_table.tab"), emit: position_table

    script:

    """

    gatk-register ${gatk_jar}

    mkdir ${genomeFileName}
    MTBseq --step TBlist --threads ${task.cpus} --minbqual ${params.minbqual} 2>${task.process}_${genomeFileName}_err.log 1>${task.process}_${genomeFileName}_out.log
    mv  Position_Tables ./${genomeFileName}/
    """

    stub:

    """
    echo "MTBseq --step TBlist --threads ${task.cpus} --minbqual ${params.minbqual}"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir ${genomeFileName}
    mkdir ${genomeFileName}/Position_Tables
    touch ${genomeFileName}/Position_Tables/${genomeFileName}_${params.library_name}.gatk_position_table.tab

    """

}
