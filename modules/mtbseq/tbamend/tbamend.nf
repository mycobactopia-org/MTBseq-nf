nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar


params.resultsDir = "${params.outdir}/tbamend"
params.saveMode = 'copy'
params.shouldPublish = true
params.project_name = "mtbseq"

process TBAMEND {
    tag "${project_name}"
    publishDir params.resultsDir, mode: params.saveMode, enabled: params.shouldPublish

    input:
    tuple path(samples_file), path("${params.mtbseq_project_name}_*_samples.tab")
    path(gatk_jar)
    env USER

    output:
    path("Amend/*")
    path  ("Amend/*_joint_*_samples_amended_*_phylo_*.tab"), emit: samples_amended

    script:

    """

    gatk-register ${gatk_jar}

    mkdir Amend
    MTBseq --step TBamend --samples ${samples_file} --project ${params.project_name}

    """
    stub:
    """
    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir Amend
    touch Amend/${project_name}_joint_cf4_cr4_fr75_ph4_samples5_amended.tab
    touch Amend/${project_name}_joint_cf4_cr4_fr75_ph4_samples5_amended_u95_phylo.tab
    touch Amend/${project_name}_joint_cf4_cr4_fr75_ph4_samples5_amended_u95_phylo.fasta
    touch Amend/${project_name}_joint_cf4_cr4_fr75_ph4_samples5_amended_u95_phylo.plainIDs.fasta
    touch Amend/${project_name}_joint_cf4_cr4_fr75_ph4_samples5_amended_u95_phylo_w12.tab
    touch Amend/${project_name}_joint_cf4_cr4_fr75_ph4_samples5_amended_u95_phylo_w12.fasta
    touch Amend/${project_name}_joint_cf4_cr4_fr75_ph4_samples5_amended_u95_phylo_w12.plainIDs.fasta
    touch Amend/${project_name}_joint_cf4_cr4_fr75_ph4_samples5_amended_u95_phylo_w12_removed.tab

    echo "MTBseq --step TBamend --samples ${samples_file} --project ${params.project_name}"
    """
}
