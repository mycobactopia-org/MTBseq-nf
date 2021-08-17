nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar


params.resultsDir = "${params.outdir}/tbstats"
params.saveMode = 'copy'
params.shouldPublish = true

// TODO: Add the tbjoin workflow
process TBSTATS {
    tag "${genomeFileName}"
    publishDir params.resultsDir, mode: params.saveMode, enabled: params.shouldPublish

    input:
    tuple val(genomeFileName), path("Bam/${genomeFileName}_${params.library_name}*.bam"), path("Position_Tables/${genomeFileName}_${params.library_name}*.gatk_position_table.tab")
    path(gatk_jar)
    env USER

    output:
    path("${genomeFileName}/Statistics/Mapping_and_Variant_Statistics.tab")

    script:

    """

    gatk-register ${gatk_jar}

    mkdir ${genomeFileName}
    MTBseq --step TBstats 2>${task.process}_${genomeFileName}_err.log 1>${task.process}_${genomeFileName}_out.log
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
