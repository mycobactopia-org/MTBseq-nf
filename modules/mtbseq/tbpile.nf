process TBPILE {
    tag "${meta.id} - ${params.project}"
    label 'process_medium'
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish
    stageInMode 'copy'

    input:
        tuple val(meta.id), path("GATK_Bam/*")
        env(USER)
        tuple path(ref_resistance_list), path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)

    output:
        path("Mpileup/${meta.id}_${params.library_name}*.gatk.{mpileup,mpileuplog}")
        tuple val(meta), path("Mpileup/${meta.id}_${params.library_name}*.gatk.mpileup"), emit: mpileup

    script:

        """
        mkdir Mpileup

        ${params.mtbseq_path} --step TBpile \\
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
        echo "${params.mtbseq_path} --step TBpile \
            --threads ${task.cpus} \
            --project ${params.project} \
            --resilist ${ref_resistance_list} \
            --intregions ${ref_interesting_regions} \
            --categories ${ref_gene_categories} \
            --basecalib ${ref_base_quality_recalibration}"

        sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

        mkdir Mpileup
        touch Mpileup/${meta.id}_${params.library_name}.gatk.mpileup
        touch Mpileup/${meta.id}_${params.library_name}.gatk.mpileuplog

        """

}
