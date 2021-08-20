nextflow.enable.dsl = 2
// NOTE: To properly setup the gatk inside the docker image
// - Download the gatk-3.8.0 tar file from here https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
// - tar -xvf GATK_TAR_FILE
// - gatk-register gatk_folder/gatk_jar


params.results_dir = "${params.outdir}/tbjoin"
params.save_mode = 'copy'
params.should_publish = true
params.project_name = "mtbseq"
params.mincovf = 4
params.mincovr = 4
params.minphred = 4
params.minfreq = 75


// TODO: Add the tbjoin workflow
process TBJOIN {
    tag "${params.project_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    path(samples_file)
    path("Called/*")
    path("Position_Tables/*")
    path(gatk_jar)
    env USER

    output:
    path ("Joint/${params.project_name}_joint*samples.{tab,log}")
    tuple path(samples_file), path("Joint/${params.project_name}_joint*samples.tab"), emit: joint_samples

    script:
    """
    gatk-register ${gatk_jar}

    mkdir Joint
    MTBseq --step TBjoin --samples ${samples} \
        --project ${params.mtbseq_project_name} \
        --mincovf ${params.mincovf} \
        --mincovr ${params.mincovr} \
        --minphred ${params.minphred} \
        --minfreq ${params.minfreq} \
        2>${task.process}_${params.project_name}_err.log 1>${task.process}_${mtbseq_project_name}_out.log
    """
    stub:

    """
    echo "MTBseq --step TBjoin --samples ${samples_file} \
        --project ${params.project_name} \
        --mincovf ${params.mincovf} \
        --mincovr ${params.mincovr} \
        --minphred ${params.minphred} \
        --minfreq ${params.minfreq}"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir Joint
    touch Joint/${params.mtbseq_project_name}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5.tab
    touch Joint/${params.mtbseq_project_name}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5.log

    """

}
