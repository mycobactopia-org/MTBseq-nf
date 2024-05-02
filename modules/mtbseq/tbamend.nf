process TBAMEND {
    tag "cohort"
    label 'process_single'
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("Joint/*")
        path(samplesheet_tsv)
        env(USER)
        tuple path(ref_resistance_list), path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)

    output:
        path("Amend/*"), emit: samples_amended

    script:

    def args = task.ext.args ?: "--project ${params.project} --mincovf ${params.mincovf} --mincovr ${params.mincovr} --minphred ${params.minphred} --minfreq ${params.minfreq} --unambig ${params.unambig} --window ${params.window} --distance ${params.distance}"

        """
        mkdir Amend

        ${params.mtbseq_path} --step TBamend \\
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

        echo " ${params.mtbseq_path} --step TBamend \
            --threads ${task.cpus} \
            --samples ${samplesheet_tsv} \
            --project ${params.project} \
            --mincovf ${params.mincovf} \
            --mincovr ${params.mincovr} \
            --minphred ${params.minphred} \
            --minfreq ${params.minfreq} \
            --unambig ${params.unambig} \
            --window ${params.window} \
            --distance ${params.distance} \
            --resilist ${ref_resistance_list} \
            --intregions ${ref_interesting_regions} \
            --categories ${ref_gene_categories} \
            --basecalib ${ref_base_quality_recalibration}  "


        mkdir Amend
        touch Amend/${params.project}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5_amended.tab
        touch Amend/${params.project}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5_amended_u95_phylo.tab
        touch Amend/${params.project}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5_amended_u95_phylo.fasta
        touch Amend/${params.project}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5_amended_u95_phylo.plainIDs.fasta
        touch Amend/${params.project}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5_amended_u95_phylo_w12.tab
        touch Amend/${params.project}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5_amended_u95_phylo_w12.fasta
        touch Amend/${params.project}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5_amended_u95_phylo_w12.plainIDs.fasta
        touch Amend/${params.project}_joint_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_samples5_amended_u95_phylo_w12_removed.tab

        """

}
