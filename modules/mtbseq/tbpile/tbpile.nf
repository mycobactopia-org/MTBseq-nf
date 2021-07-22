nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar


params.resultsDir = "${params.outdir}/tbpile"
params.saveMode = 'copy'
params.shouldPublish = true

// TODO: Add the tbjoin workflow
process TBPILE {
    tag "${genomeFileName}"
    publishDir params.resultsDir, mode: params.saveMode, enabled: params.shouldPublish

    input:
    tuple val(genomeFileName), path("GATK_Bam/${genomeFileName}_${params.library_name}_*gatk.bam")
    path(gatk_jar)
    env USER

    output:
    path("${genomeFileName}/Mpileup/${genomeFileName}_${params.library_name}*.gatk.{mpileup,mpileuplog}")
    tuple val(genomeFileName), path("${genomeFileName}/Mpileup/${genomeFileName}_${params.library_name}*.gatk.mpileup"), emit: mpileup
    script:

    """

    gatk-register ${gatk_jar}

    mkdir ${genomeFileName}
    MTBseq --step TBpile --threads ${task.cpus}
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
