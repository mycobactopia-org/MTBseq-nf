nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/tbgroups"
params.save_mode = 'copy'
params.should_publish = true
params.project_name = "mtbseq"

process TBGROUPS {
    tag "${params.project_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    path("Amend/*")
    path(gatk_jar)
    env(USER)

    output:
    path("Groups/${params.project_name}_joint_*_samples_amended_*_phylo_*.{matrix,groups}")

    script:
    """
    gatk-register ${gatk_jar}

    mkdir Groups

    MTBseq --step TBgroups \
    --threads ${task.cpus} \
    --project ${params.project_name} \
    1>>.command.out \
    2>>.command.err \
    || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

    """

    stub:
    """
    echo "MTBseq --step TBgroups  --project ${params.project_name}"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir Groups
    touch Groups/${params.project_name}_joint_cf4_cr4_fr75_ph4_samples5_amended_u95_phylo_w12.matrix
    touch Groups/${params.project_name}_joint_cf4_cr4_fr75_ph4_samples35_amended_u95_phylo_w12_d12.groups

    """
}
