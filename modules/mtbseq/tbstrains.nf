process TBSTRAINS {
    tag "${params.project}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    path("Position_Tables/*")
    path(gatk_jar)
    env(USER)

    output:
    path("Classification/Strain_Classification.tab")

    script:

    """

    gatk-register ${gatk_jar}


    mkdir Classification

    MTBseq --step TBstrains \
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
    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    echo "MTBseq --step TBstrains \
        --threads ${task.cpus} \
        --mincovf ${params.mincovf} \
        --mincovr ${params.mincovr} \
        --minphred ${params.minphred} \
        --minfreq ${params.minfreq}"

    mkdir Classification
    touch Classification/Strain_Classification.tab

    """

}
