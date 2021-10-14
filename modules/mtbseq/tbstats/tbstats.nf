nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/tbstats"
params.save_mode = 'copy'
params.should_publish = true

process TBSTATS {
    tag "${params.project_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    path("Bam/*")
    path("Position_Tables/*")
    path(gatk_jar)
    tuple path("${params.mtb_ref_name}.*"), path(resilist), path(intregions), path(categories), path(basecalib)
    env(USER)

    output:
    path("Statistics/Mapping_and_Variant_Statistics.tab")

    script:

    """
    gatk-register ${gatk_jar}

    # setting up the references as requested by MTBseq manual
    mv ${params.mtb_ref_name}.* /MTBseq_source/var/ref/.

    mkdir Statistics

    MTBseq --step TBstats \
    --threads ${task.cpus} \
    --ref ${params.mtb_ref_name} \
    --resilist ${resilist} \
    --intregions ${intregions} \
    --categories ${categories} \
    1>>.command.out \
    2>>.command.err \
    || true               # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

    """

    stub:

    """
    echo "MTBseq --step TBstats --threads ${task.cpus}"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir Statistics
    touch Statistics/Mapping_and_Variant_Statistics.tab
    """

}
