process TBSTATS {
    tag "cohort"
    label 'process_single_high_memory'
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    conda "bioconda::mtbseq=1.1.0"
    container "${'quay.io/biocontainers/mtbseq:1.1.0--hdfd78af_0'}"
    input:
        path("Bam/*")
        path("Position_Tables/*")
        env(USER)
        tuple path(ref_resistance_list), path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)

    output:
        path("Statistics/Mapping_and_Variant_Statistics.tab"), emit: statistics

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
