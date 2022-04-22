#Checking MTBseq
mtbseq_logfile="MTBseq*.log"

if ls $mtbseq_logfile 1> /dev/null 2>&1; then
    echo "$mtbseq_logfile exists."
    #Count Errors in MTBseq*log
    mtbseq_error_quantity=`cat $mtbseq_logfile | grep -c "ERROR"`

    if [ "$(($error_quantity))" > 0 ]
    then
        echo "[LOG] Found Error in $mtbseq_logfile"
        grep "ERROR" $mtbseq_logfile
    fi
fi

# Checking GATK_Bam
gatk_bam_logfile="GATK_Bam/*gatk.bamlog"
if ls $gatk_bam_logfile 1> /dev/null 2>&1; then
    echo "$gatk_bam_logfile exists."
    #Count Errors in GATK_Bam
    gatk_bam_error_quantity=`cat $gatk_bam_logfile | grep -c "ERROR"`

    if [ "$(($gatk_bam_error_quantity))" > 0 ]
    then
        echo "[LOG] Found Error in $gatk_bam_logfile"
        grep "ERROR" $gatk_bam_logfile
    fi
fi

#Checking BWA_Bam
bwa_bam_logfile="Bam/*bamlog"
if ls $bwa_bam_logfile 1> /dev/null 2>&1; then
    echo "$bwa_bam_logfile exists."
    #Count Errors in Bam/*bamlog
    bwa_bam_error_quantity=`cat  $bwa_bam_logfile | grep -c --ignore-case "error"`

    if [ "$(($error_quantity))" > 0 ]
    then
        echo "[LOG] Found Error in $bwa_bam_logfile"
        grep --ignore-case "error" $bwa_bam_logfile
    fi
fi
