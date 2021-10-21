process TBAMEND {
    tag "${params.project}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    path("Joint/*")
    path(samplesheet_tsv)
    path(gatk_jar)
    tuple path("${ref_reference_genome_name}.*"), path(ref_resistance_list), path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)
    env(USER)

    output:
    path("Amend/*"), emit: samples_amended

    script:

    """
    gatk-register ${gatk_jar}

    mkdir Amend

    MTBseq --step TBamend \
        --threads ${task.cpus} \
        --samples ${samplesheet_tsv} \
        --project ${params.project_name} \
        --mincovf ${params.mincovf} \
        --mincovr ${params.mincovr} \
        --minphred ${params.minphred} \
        --minfreq ${params.minfreq} \
        --unambig ${params.unambig} \
        --window ${params.window} \
        --distance ${params.distance} \
        --ref ${ref_reference_genome_name} \
        --resilist ${ref_resistance_list} \
        --intregions ${ref_interesting_regions} \
        --categories ${ref_gene_categories} \
    1>>.command.out \
    2>>.command.err \
    || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq


    """
    stub:
    """
    echo " MTBseq --step TBamend \
    --threads ${task.cpus} \
    --samples ${samplesheet_tsv} \
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
