process TBLIST {
    tag "${genomeFileName} - ${params.project}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(genomeFileName), path("Mpileup/${genomeFileName}_${params.library_name}*.gatk.mpileup")
    path(gatk_jar)
    env(USER)
    tuple path(ref_resistance_list), path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)

    output:
    path("Position_Tables/${genomeFileName}_${params.library_name}*.gatk_position_table.tab"), emit: tbjoin_input
    tuple val(genomeFileName), path("Position_Tables/${genomeFileName}_${params.library_name}*.gatk_position_table.tab"), emit: position_table_tuple
    path("Position_Tables/${genomeFileName}_${params.library_name}*.gatk_position_table.tab"), emit: position_table

    script:

    """

    gatk-register ${gatk_jar}


    mkdir Position_Tables

    ${params.mtbseq_path} --step TBlist \
        --threads ${task.cpus} \
        --project ${params.project} \
        --minbqual ${params.minbqual} \
        --resilist ${ref_resistance_list} \
        --intregions ${ref_interesting_regions} \
        --categories ${ref_gene_categories} \
        --basecalib ${ref_base_quality_recalibration} \
    1>>.command.out \
    2>>.command.err \
    || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq



    """

    stub:

    """
    echo "${params.mtbseq_path} --step TBlist \
        --threads ${task.cpus} \
        --project ${params.project} \
        --minbqual ${params.minbqual} \
        --resilist ${ref_resistance_list} \
        --intregions ${ref_interesting_regions} \
        --categories ${ref_gene_categories} \
        --basecalib ${ref_base_quality_recalibration}"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    touch ${task.process}_${genomeFileName}_out.log
    touch ${task.process}_${genomeFileName}_err.log

    mkdir Position_Tables
    touch Position_Tables/${genomeFileName}_${params.library_name}.gatk_position_table.tab

    """

}
