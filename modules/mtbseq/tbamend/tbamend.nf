nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/tbamend"
params.save_mode = 'copy'
params.should_publish = true
params.project_name = "mtbseq"
params.mincovf = 4
params.mincovr = 4
params.minphred = 4
params.minfreq = 75
params.unambig = 95
params.window = 12
params.distance = 12

process TBAMEND {
    tag "${project_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple path(samples_file), path("Joint/${params.project_name}_*_samples.tab")
    path(gatk_jar)
    env(USER)

    output:
    path("Amend/*")
    path  ("Amend/*_joint_*_samples_amended_*_phylo_*.tab"), emit: samples_amended

    script:

    """
    gatk-register ${gatk_jar}

    mkdir Amend

    MTBseq --step TBamend \
    --threads ${task.cpus} \
    --samples ${samples_file} \
    --project ${params.project_name} \
    --mincovf ${params.mincovf} \
    --mincovr ${params.mincovr} \
    --minphred ${params.minphred} \
    --minfreq ${params.minfreq} \
    --unambig ${params.unambig} \
    --window ${params.window} \
    --distance ${params.distance} \
    1>>.command.out \
    2>>.command.err \
    || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq


    """
    stub:
    """
    echo "MTBseq --step TBamend --samples ${samples_file} \
    --project ${params.project_name} \
    --mincovf ${params.mincovf} \
    --mincovr ${params.mincovr} \
    --minphred ${params.minphred} \
    --minfreq ${params.minfreq} \
    --unambig ${params.unambig} \
    --window ${params.window} \
    --distance ${params.distance} "

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir Amend
    touch Amend/${project_name}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5_amended.tab
    touch Amend/${project_name}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5_amended_u95_phylo.tab
    touch Amend/${project_name}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5_amended_u95_phylo.fasta
    touch Amend/${project_name}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5_amended_u95_phylo.plainIDs.fasta
    touch Amend/${project_name}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5_amended_u95_phylo_w12.tab
    touch Amend/${project_name}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5_amended_u95_phylo_w12.fasta
    touch Amend/${project_name}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5_amended_u95_phylo_w12.plainIDs.fasta
    touch Amend/${project_name}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5_amended_u95_phylo_w12_removed.tab

    """

}
