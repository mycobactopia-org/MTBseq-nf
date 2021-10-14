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
    path(samplesheet_tsv)
    path(gatk_jar)
    tuple path("${params.mtb_ref_name}.*"), path(resilist), path(intregions), path(categories), path(basecalib)
    env(USER)

    output:
    path("Groups/*")

    script:
    """
    gatk-register ${gatk_jar}

    # setting up the references as requested by MTBseq manual
    mv ${params.mtb_ref_name}.* /MTBseq_source/var/ref/.

    mkdir Groups

    MTBseq --step TBgroups \
    --threads ${task.cpus} \
    --samples ${samplesheet_tsv} \
    --project ${params.project_name} \
    --ref ${params.mtb_ref_name} \
    --resilist ${resilist} \
    --intregions ${intregions} \
    --categories ${categories} \
    1>>.command.out \
    2>>.command.err \
    || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

    """

    stub:
    """
    echo "MTBseq --step TBgroups \
    --threads ${task.cpus} \
    --samples ${samplesheet_tsv} \
    --project ${params.project_name}"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir Groups
    touch Groups/${params.project_name}_joint_cf4_cr4_fr75_ph4_samples5_amended_u95_phylo_w12.matrix
    touch Groups/${params.project_name}_joint_cf4_cr4_fr75_ph4_samples35_amended_u95_phylo_w12_d12.groups

    """
}
