process TBVARIANTS {
    tag "${genomeFileName} - ${params.project}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(genomeFileName), path("Position_Tables/*")
    path(gatk_jar)
    env(USER)
    tuple path(ref_resistance_list), path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)

    output:
    path("Called/${genomeFileName}_${params.library_name}*gatk_position_{uncovered,variants}*.tab")
    path("Called/${genomeFileName}_${params.library_name}*gatk_position_variants*.tab"), emit: tbjoin_input

    script:

    """

    gatk-register ${gatk_jar}


    mkdir Called

    ${params.mtbseq_path} --step TBvariants \
        --threads ${task.cpus} \
        --project ${params.project} \
        --mincovf ${params.mincovf} \
        --mincovr ${params.mincovr} \
        --minphred ${params.minphred} \
        --minfreq ${params.minfreq} \
        --resilist ${ref_resistance_list} \
        --intregions ${ref_interesting_regions} \
        --categories ${ref_gene_categories} \
        --basecalib ${ref_base_quality_recalibration} \
    1>>.command.out \
    2>>.command.err \
    || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

    check_logs.sh
    """

    stub:

    """
    echo "${params.mtbseq_path} --step TBvariants \
        --threads ${task.cpus} \
        --project ${params.project} \
        --mincovf ${params.mincovf} \
        --mincovr ${params.mincovr} \
        --minphred ${params.minphred} \
        --minfreq ${params.minfreq} \
        --resilist ${ref_resistance_list} \
        --intregions ${ref_interesting_regions} \
        --categories ${ref_gene_categories} \
        --basecalib ${ref_base_quality_recalibration}"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir Called
    touch Called/${genomeFileName}_${params.library_name}.gatk_position_uncovered_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_outmode000.tab
    touch Called/${genomeFileName}_${params.library_name}.gatk_position_variants_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_outmode000.tab

    """

}
