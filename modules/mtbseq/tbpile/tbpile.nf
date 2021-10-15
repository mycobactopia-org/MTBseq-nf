nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/tbpile"
params.save_mode = 'copy'
params.should_publish = true

process TBPILE {
    tag "${genomeFileName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish
    stageInMode 'copy'

    input:
    tuple val(genomeFileName), path("GATK_Bam/*")
    path(gatk_jar)
    tuple path("${ref_reference_genome_name}.*"), path(ref_resistance_list), path(ref_interesting_regions), path(ref_gene_categories), path(ref_base_quality_recalibration)
    env(USER)

    output:
    path("Mpileup/${genomeFileName}_${params.library_name}*.gatk.{mpileup,mpileuplog}")
    tuple val(genomeFileName), path("Mpileup/${genomeFileName}_${params.library_name}*.gatk.mpileup"), emit: mpileup

    script:

    """

    gatk-register ${gatk_jar}

    # setting up the references as requested by MTBseq manual
    mv ${ref_reference_genome_name}.* /MTBseq_source/var/ref/.

    mkdir Mpileup

    MTBseq --step TBpile \
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
    echo "MTBseq --step TBpile --threads ${task.cpus}"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir Mpileup
    touch Mpileup/${genomeFileName}_${params.library_name}.gatk.mpileup
    touch Mpileup/${genomeFileName}_${params.library_name}.gatk.mpileuplog

    """

}
