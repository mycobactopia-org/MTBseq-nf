process TBJOIN {
    tag "${params.project}"
    label 'process_high'
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("Called/*")
        path("Position_Tables/*")
        path(samplesheet_tsv)
        env(USER)
        tuple path(ref_resistance_list), path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)

    output:
        path("Joint/${params.project}_joint*samples*.{tab,log}")
        path("Joint/${params.project}_joint*samples*.tab"), emit: joint_samples

    script:
        def args = task.ext.args ?: " --project ${params.project} --mincovf ${params.mincovf} --mincovr ${params.mincovr} --minphred ${params.minphred} --minfreq ${params.minfreq} --distance ${params.distance}"
        """
        mkdir Joint

        ${params.mtbseq_path} --step TBjoin \\
            --threads ${task.cpus} \\
            --samples ${samplesheet_tsv} \\
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

        echo "${params.mtbseq_path} --step TBjoin \
            --threads ${task.cpus} \
            --samples ${samplesheet_tsv} \
            --project ${params.project} \
            --mincovf ${params.mincovf} \
            --mincovr ${params.mincovr} \
            --minphred ${params.minphred} \
            --minfreq ${params.minfreq} \
            --resilist ${ref_resistance_list} \
            --intregions ${ref_interesting_regions} \
            --categories ${ref_gene_categories} \
            --basecalib ${ref_base_quality_recalibration}"


        mkdir Joint
        touch Joint/${params.project}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5.tab
        touch Joint/${params.project}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5.log

        """

}
