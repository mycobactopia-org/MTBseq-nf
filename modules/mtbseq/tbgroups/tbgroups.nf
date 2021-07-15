nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar


params.resultsDir = "${params.outdir}/tbgroups"
params.saveMode = 'copy'
params.shouldPublish = true

// TODO: Add the tbjoin workflow
process TBGROUPS {
    tag "${project_name}"
    publishDir params.resultsDir, mode: params.saveMode, enabled: params.shouldPublish

    input:
    tuple val(project_name), path("Amend/${params.mtbseq_project_name}*_joint_*_samples_amended_*_phylo_*.tab")
    path(gatk_jar)
    env USER

    output:

    path("Groups/${params.mtbseq_project_name}_joint_*_samples_amended_*_phylo_*.{matrix,groups}")

    script:
    """
    gatk-register ${gatk_jar}
    mkdir Groups
    MTBseq --step TBgroups --project ${project_name}
    """

    stub:
    """
    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir Groups
    touch Groups/${params.mtbseq_project_name}_joint_cf4_cr4_fr75_ph4_samples5_amended_u95_phylo_w12.matrix
    touch Groups/${params.mtbseq_project_name}_joint_cf4_cr4_fr75_ph4_samples35_amended_u95_phylo_w12_d12.groups
    echo "MTBseq --step TBgroups  --project ${project_name}"
    """
}
