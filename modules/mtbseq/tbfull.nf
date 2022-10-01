process TBFULL {
    tag "${params.project}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("*")
        path(gatk_jar)
        env(USER)
        tuple path(ref_resistance_list), path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)

    output:
        path("Called")
        path("Position_Tables")
        path("Classification")
        path("Statistics")
        path("Called/*_${params.library_name}*gatk_position_variants*.tab"), emit: position_variants
        path("Position_Tables/*_${params.library_name}*.gatk_position_table.tab"), emit: position_tables

    script:

        """

        ${ params.load_gatk38_jar ? "gatk-register ${gatk_jar}" : ""}


        ${params.mtbseq_path} --step TBfull \
            --thread ${task.cpus} \
            --project ${params.project} \
            --minbqual ${params.minbqual} \
            --mincovf ${params.mincovf} \
            --mincovr ${params.mincovr} \
            --minphred ${params.minphred} \
            --minfreq ${params.minfreq} \
            --resilist ${ref_resistance_list} \
            --unambig ${params.unambig} \
            --window ${params.window} \
            --distance ${params.distance} \
            --intregions ${ref_interesting_regions} \
            --categories ${ref_gene_categories} \
            --basecalib ${ref_base_quality_recalibration} \
        1>>.command.out \
        2>>.command.err \
        || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

        """

    stub:
        """
        sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

        echo " ${params.mtbseq_path} --step TBfull \
            --thread ${task.cpus} \
            --project ${params.project} \
            --minbqual ${params.minbqual} \
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
            --basecalib ${ref_base_quality_recalibration} "


        mkdir GATK_Bam
        touch GATK_Bam/${genomeFileName}_${params.library_name}.gatk.bam
        touch GATK_Bam/${genomeFileName}_${params.library_name}.gatk.bai
        touch GATK_Bam/${genomeFileName}_${params.library_name}.gatk.bamlog
        touch GATK_Bam/${genomeFileName}_${params.library_name}.gatk.grp
        touch GATK_Bam/${genomeFileName}_${params.library_name}.gatk.intervals
        mkdir Bam
        mkdir Bam/${genomeFileName}
        touch Bam/${genomeFileName}_${params.library_name}.bam
        touch Bam/${genomeFileName}_${params.library_name}.bai
        touch Bam/${genomeFileName}_${params.library_name}.bamlog
        mkdir Called
        touch Called/${genomeFileName}_${params.library_name}.gatk_position_uncovered_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_outmode000.tab
        touch Called/${genomeFileName}_${params.library_name}.gatk_position_variants_cf${params.mincovf}_cr${params.mincovr}_fr${params.minfreq}_ph${params.minphred}_outmode000.tab
        mkdir Mpileup
        touch Mpileup/${genomeFileName}_${params.library_name}.gatk.mpileup
        touch Mpileup/${genomeFileName}_${params.library_name}.gatk.mpileuplog
        mkdir Classification
        touch Classification/Strain_Classification.tab
        mkdir Position_Tables
        touch Position_Tables/${genomeFileName}_${params.library_name}.gatk_position_table.tab
        mkdir Statistics
        touch Statistics/Mapping_and_Variant_Statistics.tab

        """

}
