nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar


params.resultsDir = "${params.outdir}/mtbseq"
params.saveMode = 'copy'
params.shouldPublish = true

// TODO: Add the tbjoin workflow
process MTBSEQ {
    tag "${genomeFileName}"
    publishDir params.resultsDir, mode: params.saveMode, enabled: params.shouldPublish
    container 'quay.io/biocontainers/mtbseq:1.0.3--pl526_1'
    cpus 8
    memory "15 GB"

    input:
    tuple val(genomeFileName), path("${genomeFileName}_somelib_R?.fastq.gz")
    path(gatk_jar)
    env USER

    output:
    path("${genomeFileName}")

    script:

    """
    #ignore the silly exit 1 status even after successful execution
    set +e

    gatk-register ${gatk_jar}

    mkdir ${genomeFileName}
   
    MTBseq --step TBfull --thread ${task.cpus}
    
    mv  Amend ./${genomeFileName}/
    mv  Bam ./${genomeFileName}/
    mv  Called ./${genomeFileName}/
    mv  Classification ./${genomeFileName}/
    mv  GATK_Bam ./${genomeFileName}/
    mv  Groups ./${genomeFileName}/
    mv  Joint ./${genomeFileName}/
    mv  Mpileup ./${genomeFileName}/
    mv  Position_Tables ./${genomeFileName}/
    mv  Statistics ./${genomeFileName}/
    """


}

