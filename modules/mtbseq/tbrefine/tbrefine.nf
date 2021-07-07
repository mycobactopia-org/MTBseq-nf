nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar


params.resultsDir = "${params.outdir}/tbrefine"
params.saveMode = 'copy'
params.shouldPublish = true

// TODO: Add the tbjoin workflow
process TBREFINE {
    tag "${genomeFileName}"
    publishDir params.resultsDir, mode: params.saveMode, enabled: params.shouldPublish

    input:
    tuple val(genomeFileName), path("${genomeFileName}_${params.library_name}_*.bam")
    path(gatk_jar)
    env USER

    output:
    path("${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}_*.gatk.{bam,bai,bamlog,grp,intervals}")

    script:

    """

    gatk-register ${gatk_jar}

    mkdir ${genomeFileName}
    MTBseq --step TBrefine --thread ${task.cpus}
    mv  GATK_bam ./${genomeFileName}/
    """

    stub:

    """
    mkdir ${genomeFileName}
    mkdir ${genomeFileName}/GATK_bam
    mkdir ${genomeFileName}/GATK_bam/${genomeFileName}
    touch ${genomeFileName}/GATK_bam/${genomeFileName}/${genomeFileName}_${params.library_name}.gatk.bam
    touch ${genomeFileName}/GATK_bam/${genomeFileName}/${genomeFileName}_${params.library_name}.gatk.bai
    touch ${genomeFileName}/GATK_bam/${genomeFileName}/${genomeFileName}_${params.library_name}.gatk.bamlog
    touch ${genomeFileName}/GATK_bam/${genomeFileName}/${genomeFileName}_${params.library_name}.gatk.grp
    touch ${genomeFileName}/GATK_bam/${genomeFileName}/${genomeFileName}_${params.library_name}.gatk.intervals
    echo "MTBseq --step TBbwa --thread ${task.cpus}"
    """
}
