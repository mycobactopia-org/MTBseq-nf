nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/tbbwa"
params.save_mode = 'copy'
params.should_publish = true

process TBBWA {
    tag "${genomeFileName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(genomeFileName), path("${genomeFileName}_${params.library_name}_R?.fastq.gz")
    path(gatk_jar)
    tuple path(${params.mtb_ref_name}.*), path(resilist), path(intregions), path(categories), path(basecalib)
    env(USER)


    output:
    path("Bam/${genomeFileName}_${params.library_name}*.{bam,bai,bamlog}")
    tuple val(genomeFileName), path("Bam/${genomeFileName}_${params.library_name}*.bam"), emit: bam_tuple
    path("Bam/${genomeFileName}_${params.library_name}*.bam"), emit: bam

    script:

    """

    gatk-register ${gatk_jar}

    mkdir Bam

    MTBseq --step TBbwa \
    --threads ${task.cpus} \
    --ref mtb_ref ${params.mtb_ref_name} \
    --resilist ${resilist} \
    --intregions ${intregions} \
    --categories ${categories} \
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

