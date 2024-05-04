process TBSTATS {
    tag "cohort"
    label 'process_max'
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("Bam/*")
        path("Position_Tables/*")
        env(USER)
        tuple path(ref_resistance_list), path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)

    output:
        path("Statistics/Mapping_and_Variant_Statistics.tab")

    script:

        """
        mkdir Statistics

        ${params.mtbseq_path} --step TBstats \\
            --threads ${task.cpus} \\
            --project ${params.project} \\
            --resilist ${ref_resistance_list} \\
            --intregions ${ref_interesting_regions} \\
            --categories ${ref_gene_categories} \\
            --basecalib ${ref_base_quality_recalibration} \\
        1>>.command.out \\
        2>>.command.err \\
        || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

        """

    stub:

        """
        sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

        echo "${params.mtbseq_path} --step TBstats \
            --threads ${task.cpus} \
            --project ${params.project} \
            --resilist ${ref_resistance_list} \
            --intregions ${ref_interesting_regions} \
            --categories ${ref_gene_categories} \
            --basecalib ${ref_base_quality_recalibration}"

        mkdir Statistics
        touch Statistics/Mapping_and_Variant_Statistics.tab
        """

}
