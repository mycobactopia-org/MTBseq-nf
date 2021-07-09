nextflow.enable.dsl = 2

workflow PER_SAMPLE_ANALYSIS {
    reads_ch = Channel.fromFilePairs(params.reads)

    TBBWA(reads_ch, params.gatk38_jar, params.user)
    TBREFINE(TBBWA.out.bam, params.gatk38_jar, params.user)
    TBPILE(TBREFINE.out.gatk_bam, params.gatk38_jar, params.user)
    TBLIST(TBPILE.out.mpileup, params.gatk38_jar, params.user)
    TBVARIANTS(TBLIST.out.position_table, params.gatk38_jar, params.user)
    TBSTATS(TBBWA.out.bam.join(TBLIST.out.position_table), params.gatk38_jar, params.user)
    TBSTRAINS(TBLIST.out.position_table, params.gatk38_jar, params.user)

}

