process TBREFINE {
    tag "${meta.id}"
    label 'process_medium'
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(meta), path("Bam/")
        env(USER)
        tuple path(ref_resistance_list), path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)

    output:
        tuple val(meta), path("GATK_Bam/${meta.id}_${params.library_name}*gatk.{bam,bai,bamlog,grp,intervals}"), emit: gatk_bam

    script:

        """
        mkdir GATK_Bam

        ${params.mtbseq_path} --step TBrefine \\
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

        echo " ${params.mtbseq_path} --step TBrefine \
            --threads ${task.cpus} \
            --project ${params.project} \
            --resilist ${ref_resistance_list} \
            --intregions ${ref_interesting_regions} \
            --categories ${ref_gene_categories} \
            --basecalib ${ref_base_quality_recalibration}"


        mkdir GATK_Bam
        touch GATK_Bam/${meta.id}_${params.library_name}.gatk.bam
        touch GATK_Bam/${meta.id}_${params.library_name}.gatk.bai
        touch GATK_Bam/${meta.id}_${params.library_name}.gatk.bamlog
        touch GATK_Bam/${meta.id}_${params.library_name}.gatk.grp
        touch GATK_Bam/${meta.id}_${params.library_name}.gatk.intervals
        """
}
