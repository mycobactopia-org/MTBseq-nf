process TBBWA {
    tag "${genomeFileName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(genomeFileName), path("${genomeFileName}_${params.library_name}_R?.fastq.gz")
    path(gatk_jar)
    tuple path("${ref_reference_genome_name}.*"), path(ref_resistance_list), path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)
    env(USER)

    output:
    path("Bam/${genomeFileName}_${params.library_name}*.{bam,bai,bamlog}")
    tuple val(genomeFileName), path("Bam/${genomeFileName}_${params.library_name}*.bam"), emit: bam_tuple
    path("Bam/${genomeFileName}_${params.library_name}*.bam"), emit: bam

    script:

    """

    gatk-register ${gatk_jar}

    mkdir Bam

    # setting up the references as requested by MTBseq manual
    mv ${ref_reference_genome_name}.* /MTBseq_source/var/ref/.

    MTBseq --step TBbwa \
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
    echo "MTBseq --step TBbwa \
    --threads ${task.cpus}"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    touch ${task.process}_${genomeFileName}_out.log
    touch ${task.process}_${genomeFileName}_err.log

    mkdir Bam
    touch Bam/${genomeFileName}_${params.library_name}.bam
    touch Bam/${genomeFileName}_${params.library_name}.bai
    touch Bam/${genomeFileName}_${params.library_name}.bamlog

    """

}

