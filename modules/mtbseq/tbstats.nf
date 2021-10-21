process TBSTATS {
    tag "${params.project}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    path("Bam/*")
    path("Position_Tables/*")
    path(gatk_jar)
    tuple path("${ref_reference_genome_name}.*"), path(ref_resistance_list), path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)
    env(USER)

    output:
    path("Statistics/Mapping_and_Variant_Statistics.tab")

    script:

    """
    gatk-register ${gatk_jar}


    mkdir Statistics

    MTBseq --step TBstats \
        --threads ${task.cpus} \
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
    echo "MTBseq --step TBstats --threads ${task.cpus}"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir Statistics
    touch Statistics/Mapping_and_Variant_Statistics.tab
    """

}
