nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/tbjoin"
params.save_mode = 'copy'
params.should_publish = true
params.project_name = "mtbseq"
params.mincovf = 4
params.mincovr = 4
params.minphred = 4
params.minfreq = 75


process TBJOIN {
    tag "${params.project_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    path("Called/*")
    path("Position_Tables/*")
    path(samplesheet_tsv)
    path(gatk_jar)
    tuple path("${params.mtb_ref_name}.*"), path(resilist), path(intregions), path(categories), path(basecalib)
    env(USER)

    output:
    path("Joint/${params.project_name}_joint*samples*.{tab,log}")
    path("Joint/${params.project_name}_joint*samples*.tab"), emit: joint_samples

    script:
    """
    gatk-register ${gatk_jar}

    # setting up the references as requested by MTBseq manual
    mv ${params.mtb_ref_name}.* /MTBseq_source/var/ref/.

    mkdir Joint

    MTBseq --step TBjoin \
    --threads ${task.cpus} \
    --samples ${samplesheet_tsv} \
    --project ${params.project_name} \
    --mincovf ${params.mincovf} \
    --mincovr ${params.mincovr} \
    --minphred ${params.minphred} \
    --minfreq ${params.minfreq} \
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
    echo "MTBseq --step TBjoin \
    --threads ${task.cpus} \
    --samples ${samplesheet_tsv} \
    --project ${params.project_name} \
    --mincovf ${params.mincovf} \
    --mincovr ${params.mincovr} \
    --minphred ${params.minphred} \
    --minfreq ${params.minfreq}"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir Joint
    touch Joint/${params.project_name}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5.tab
    touch Joint/${params.project_name}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5.log

    """

}
