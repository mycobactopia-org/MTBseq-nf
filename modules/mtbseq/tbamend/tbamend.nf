nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar


params.resultsDir = "${params.outdir}/tbamend"
params.saveMode = 'copy'
params.shouldPublish = true

process TBAMEND {
    tag "${project_name}"
    publishDir params.resultsDir, mode: params.saveMode, enabled: params.shouldPublish

    input:
    tuple path(samples_file), val(project_name), path("${params.mtbseq_project_name}_*_samples.tab")
    path(gatk_jar)
    env USER

    output:
    path("Amend/*")
    tuple val(project_name), path  ("Amend/*_joint_*_samples_amended_*_phylo_*.tab"), emit: samples_amended

    script:

    """

    gatk-register ${gatk_jar}

    mkdir Amend
    MTBseq --step TBamend --samples ${samples_file} --project ${project_name}

    """
    stub:
    """
    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir Amend
    touch Amend/${params.mtbseq_project_name}_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended.tab
    touch Amend/${params.mtbseq_project_name}_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended_[unambig]_phylo.tab
    touch Amend/${params.mtbseq_project_name}_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended_[unambig]_phylo.fasta
    touch Amend/${params.mtbseq_project_name}_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended_[unambig]_phylo_plainIDs.fasta
    touch Amend/${params.mtbseq_project_name}_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended_[unambig]_phylo_[window].tab
    touch Amend/${params.mtbseq_project_name}_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended_[unambig]_phylo_[window].fasta
    touch Amend/${params.mtbseq_project_name}_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended_[unambig]_phylo_[window]_plainIDs.fasta
    touch Amend/${params.mtbseq_project_name}_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended_[unambig]_phylo_[window]_removed.tab

    echo "MTBseq --step TBamend --samples ${samples_file} --project ${project_name}"
    """
}
