nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/tbrefine"
params.save_mode = 'copy'
params.should_publish = true

process TBREFINE {
    tag "${genomeFileName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(genomeFileName), path("Bam/")
    path(gatk_jar)
    tuple path("${params.mtb_ref_name}.*"), path(resilist), path(intregions), path(categories), path(basecalib)
    env(USER)

    output:
    tuple val(genomeFileName), path("GATK_Bam/${genomeFileName}_${params.library_name}*gatk.{bam,bai,bamlog,grp,intervals}"), emit: gatk_bam

    script:

    """

    gatk-register ${gatk_jar}

    # setting up the references as requested by MTBseq manual
    mv ${params.mtb_ref_name}.* /MTBseq_source/var/ref/.

    mkdir GATK_Bam

    MTBseq --step TBrefine \
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
    echo "MTBseq --step TBrefine --threads ${task.cpus}"

    sleep \$[ ( \$RANDOM % 10 )  + 1 ]s

    mkdir GATK_Bam
    touch GATK_Bam/${genomeFileName}_${params.library_name}.gatk.bam
    touch GATK_Bam/${genomeFileName}_${params.library_name}.gatk.bai
    touch GATK_Bam/${genomeFileName}_${params.library_name}.gatk.bamlog
    touch GATK_Bam/${genomeFileName}_${params.library_name}.gatk.grp
    touch GATK_Bam/${genomeFileName}_${params.library_name}.gatk.intervals
    """
}
