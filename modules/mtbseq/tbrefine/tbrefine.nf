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
    tuple val(genomeFileName), path("Bam/${genomeFileName}_${params.library_name}*.bam")
    path(gatk_jar)
    env USER

    output:
    path("${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}*.gatk.{bam,bai,bamlog,grp,intervals}")
    tuple val(genomeFileName), path("${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}*gatk.bam"), emit: gatk_bam

    script:

    """

    gatk-register ${gatk_jar}

    mkdir ${genomeFileName}
    MTBseq --step TBrefine --threads ${task.cpus}
    mv  GATK_Bam ./${genomeFileName}/
    """

    stub:

    """
    echo "MTBseq --step TBrefine --threads ${task.cpus}"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir ${genomeFileName}
    mkdir ${genomeFileName}/GATK_Bam
    touch ${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}.gatk.bam
    touch ${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}.gatk.bai
    touch ${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}.gatk.bamlog
    touch ${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}.gatk.grp
    touch ${genomeFileName}/GATK_Bam/${genomeFileName}_${params.library_name}.gatk.intervals
    """
}
