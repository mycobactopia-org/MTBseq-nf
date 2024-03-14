process TBBWA {
    tag "${meta.id} - ${params.project}"
    label 'process_medium'
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(meta), path("${meta.id}_${params.library_name}_R?.fastq.gz")
        env(USER)
        tuple path(ref_resistance_list), path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)

    output:
        path("Bam/${meta.id}_${params.library_name}*.{bam,bai,bamlog}")
        tuple val(meta), path("Bam/${meta.id}_${params.library_name}*.{bam,bai}"), emit: bam_tuple
        path("Bam/${meta.id}_${params.library_name}*.bam"), emit: bam

    script:

        """
        mkdir Bam

        ${params.mtbseq_path} --step TBbwa \\
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
        echo " ${params.mtbseq_path} --step TBbwa \
            --threads ${task.cpus} \
            --project ${params.project} \
            --resilist ${ref_resistance_list} \
            --intregions ${ref_interesting_regions} \
            --categories ${ref_gene_categories} \
            --basecalib ${ref_base_quality_recalibration} "

        sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

        touch ${task.process}_${meta.id}_out.log
        touch ${task.process}_${meta.id}_err.log

        mkdir Bam
        touch Bam/${meta.id}_${params.library_name}.bam
        touch Bam/${meta.id}_${params.library_name}.bai
        touch Bam/${meta.id}_${params.library_name}.bamlog

        """

}

