nextflow.enable.dsl = 2

include { MTBSEQ_NF } from "./workflows/mtbseq.nf"



workflow {

//--------------------------------------------
// Source the input reads
//--------------------------------------------


    //NOTE: Expected structure of input CSV samplesheet
    // sampleName,read1,read2
    // ERX1933967,/full/path/to/ERX1933967_R1.fastq.gz,/full/path/to/ERX1933967_R2.fastq.gz

    reads_ch = Channel.fromPath(params.input_samplesheet)
            .splitCsv(header: false, skip: 1)
            .map { row -> {
                        sampleName        = row[0]
                        read1             = row[1]
                        read2             = row[2]
                    }

                return tuple(sampleName, tuple(file(read1), file(read2)))
            }


    MTBSEQ_NF(reads_ch)

}
