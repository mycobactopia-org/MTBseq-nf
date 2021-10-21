process TBSTRAINS {
    tag "${params.project}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    path("Position_Tables/*")
    path(gatk_jar)
    env(USER)
    tuple path(ref_resistance_list), path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)

    output:
    path("Classification/Strain_Classification.tab")

    script:

    """

    gatk-register ${gatk_jar}


    mkdir Classification

    ${params.mtbseq_path} --step TBstrains \
        --threads ${task.cpus} \
        --project ${params.project} \
        --mincovf ${params.mincovf} \
        --mincovr ${params.mincovr} \
        --minphred ${params.minphred} \
        --minfreq ${params.minfreq} \
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
    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    echo "${params.mtbseq_path} --step TBstrains \
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

    mkdir Classification
    touch Classification/Strain_Classification.tab

    """

}
