nextflow.enable.dsl = 2

include { MTBSEQ_NF } from "./workflows/mtbseq.nf"

//--------------------------------------------
// Source the input reads
//--------------------------------------------


if( params.run_type == "folder" ) {

    reads_ch = Channel.fromFilePairs(params.input_folder)

} else {
    // Default to reading from a samplesheet

    //NOTE: Expected structure of input CSV samplesheet
    // sampleName    read1    read2
    // ERX1933967    R1,      R2

    reads_ch = Channel.fromPath(params.input_samplesheet)
            .splitCsv(header: false, skip: 1)
            .map { row -> {
                        read1             = row[1]
                        read2             = row[2]
                    }

                return tuple(tuple(file(read1), file(read2)))
            }
        }
}



workflow {

    MTBSEQ_NF(reads_ch)

}
