nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/tblist"
params.save_mode = 'copy'
params.should_publish = true
params.minbqual = 13

process TBLIST {
    tag "${genomeFileName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(genomeFileName), path("Mpileup/${genomeFileName}_${params.library_name}*.gatk.mpileup")
    path(gatk_jar)
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

    mkdir ${genomeFileName}
    mkdir ${genomeFileName}/Position_Tables
    touch ${genomeFileName}/Position_Tables/${genomeFileName}_${params.library_name}.gatk_position_table.tab

    """

}
