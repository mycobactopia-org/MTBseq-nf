process TBLIST {
    tag "${genomeFileName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(genomeFileName), path("Mpileup/${genomeFileName}_${params.library_name}*.gatk.mpileup")
    path(gatk_jar)
    tuple path("${ref_reference_genome_name}.*"), path(ref_resistance_list), path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)
    env(USER)

    output:
    path("Position_Tables/${genomeFileName}_${params.library_name}*.gatk_position_table.tab"), emit: tbjoin_input
    tuple val(genomeFileName), path("Position_Tables/${genomeFileName}_${params.library_name}*.gatk_position_table.tab"), emit: position_table_tuple
    path("Position_Tables/${genomeFileName}_${params.library_name}*.gatk_position_table.tab"), emit: position_table

    script:

    """

    gatk-register ${gatk_jar}


    mkdir Position_Tables

    MTBseq --step TBlist \
    --threads ${task.cpus} \
    --minbqual ${params.minbqual} \
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
    echo "MTBseq --step TBlist --threads ${task.cpus} --minbqual ${params.minbqual}"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    touch ${task.process}_${genomeFileName}_out.log
    touch ${task.process}_${genomeFileName}_err.log

    mkdir Position_Tables
    touch Position_Tables/${genomeFileName}_${params.library_name}.gatk_position_table.tab

    """

}
