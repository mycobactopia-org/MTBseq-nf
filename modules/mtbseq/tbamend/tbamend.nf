nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar


params.resultsDir = "${params.outdir}/tbamend"
params.saveMode = 'copy'
params.shouldPublish = true

process TBAMEND {
    publishDir params.resultsDir, mode: params.saveMode, enabled: params.shouldPublish
    container 'quay.io/biocontainers/mtbseq:1.0.3--pl526_1'
    cpus 8
    memory "15 GB"

    input:
    tuple path(samples), val(project_name), path(Joint/[PROJECT]_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples.tab)
    path(gatk_jar)

    output:
    path("Amend")

    script:

    """

    gatk-register ${gatk_jar}

    MTBseq --step TBamend --thread ${task.cpus}

    """
    stub:
    """
    mkdir Amend
    touch Amend/${project_name}
    """
}
