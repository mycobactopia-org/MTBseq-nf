process TBVARIANTS {
    tag "${genomeFileName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(genomeFileName), path("Position_Tables/*")
    path(gatk_jar)
    env(USER)

    output:
    path("Called/${genomeFileName}_${params.library_name}*gatk_position_{uncovered,variants}*.tab")
    path("Called/${genomeFileName}_${params.library_name}*gatk_position_variants*.tab"), emit: tbjoin_input

    script:

    """

    gatk-register ${gatk_jar}

    # setting up the references as requested by MTBseq manual
    mv ${ref_reference_genome_name}.* /MTBseq_source/var/ref/.

    mkdir Called

    MTBseq --step TBvariants \
    --threads ${task.cpus} \
    --mincovf ${params.mincovf} \
    --mincovr ${params.mincovr} \
    --minphred ${params.minphred} \
    --minfreq ${params.minfreq} \
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
    echo "MTBseq --step TBvariants \
    --threads ${task.cpus} \
    --mincovf ${params.mincovf} \
    --mincovr ${params.mincovr} \
    --minphred ${params.minphred} \
    --minfreq ${params.minfreq}"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir Called
    touch Called/${genomeFileName}_${params.library_name}.gatk_position_uncovered_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_outmode000.tab
    touch Called/${genomeFileName}_${params.library_name}.gatk_position_variants_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_outmode000.tab

    """

}
