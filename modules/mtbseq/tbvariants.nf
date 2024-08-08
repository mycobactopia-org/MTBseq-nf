process TBVARIANTS {
    tag "${meta.id}"
    label 'process_single'
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    conda "bioconda::mtbseq=1.1.0"
    container "${'bquay.io/biocontainers/mtbseq:1.1.0--hdfd78af_0'}"
    input:
        tuple val(meta), path("Position_Tables/*")
        env(USER)
        tuple path(ref_resistance_list), path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)

    output:
        path("Called/${meta.id}_${params.library_name}*gatk_position_{uncovered,variants}*.tab")
        path("Called/${meta.id}_${params.library_name}*gatk_position_variants*.tab"), emit: tbjoin_input

    script:
        def args = task.ext.args ?: "--mincovf ${params.mincovf} --mincovr ${params.mincovr} --minphred ${params.minphred} --minfreq ${params.minfreq}"

        """
        mkdir Called

        ${params.mtbseq_path} --step TBvariants \\
            --threads ${task.cpus} \\
            --project ${params.project} \\
            --resilist ${ref_resistance_list} \\
            --intregions ${ref_interesting_regions} \\
            --categories ${ref_gene_categories} \\
            --basecalib ${ref_base_quality_recalibration} \\
            ${args} \\
        1>>.command.out \\
        2>>.command.err \\
        || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

        """

    stub:

        """
        sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

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


        mkdir Called
        touch Called/${meta.id}_${params.library_name}.gatk_position_uncovered_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_outmode000.tab
        touch Called/${meta.id}_${params.library_name}.gatk_position_variants_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_outmode000.tab

        """

}
