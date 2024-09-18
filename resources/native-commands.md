
# Table of Contents

1.  [Native commands used by mtbseq based on the log of tbfull](#org63d26ae)
    1.  [bwa mem  - TBbwa](#org2fd46ba)
    2.  [samtools view - TBbwa](#org781fe6c)
    3.  [samtools sort - TBbwa](#org8ae18b1)
    4.  [samtools index - TBbwa](#org04f2c55)
    5.  [samtools rmdup - TBbwa](#orge57e785)
    6.  [samtools index - TBbwa](#org94f77e7)
    7.  [bwa mapping? - TBbwa](#org214b6a8)
        1.  [GATK RealignerTargetCreator](#orgd4bd11c)
        2.  [GATK IndelRealigner](#org5851d40)
        3.  [GATK BaseRecalibrator](#orge7101fe)
        4.  [GATK PrintReads](#orgaea9fea)
    8.  [mpileup creation - TBpile](#org36f0bc5)
        1.  [samtools mpileup](#org969b80a)
    9.  [Creating position lists - TBlist](#orgad724d5)
    10. [Variant calling - TBVariants](#orgf67d99f)
    11. [Refining variants - TBrefine](#orgb3ba568)
        1.  [gatk RealignerTargetCreator](#org67e9383)
        2.  [gatk IndelRealigner](#orgb083c6a)
        3.  [gatk BaseRecalibrator](#org0f401c9)
        4.  [gatk PrintReads](#orgcf088a9)
    12. [Statistics for mapping - TBstats](#orgfb28061)
        1.  [samtools](#org607b5a0)
    13. [Calling Strains - TBstrains](#org3abd486)



<a id="org63d26ae"></a>

# Native commands used by mtbseq based on the log of tbfull

This are the commands found on the mtbseq log related


<a id="org2fd46ba"></a>

## bwa mem  - TBbwa

bwa for mapping

    mem -t 2 -R t 2 -R '@RG\tID:5765_somelib\tSM:5765\tPL:Illumina\tLB:somelib' /root/miniconda3/envs/mtbseq/share/mtbseq-1.0.3-1/var/ref/M._tuberculosis_H37Rv_2015-11-13.fasta  /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/5765_somelib_R1.fastq.gz /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/5765_somelib_R2.fastq.gz > /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/Bam/5765_somelib.sam 2>> /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/Bam/5765_somelib.bamlog


<a id="org781fe6c"></a>

## samtools view - TBbwa

samtools for converting from sam to bam

    samtools view -@ 2 -b -T /root/miniconda3/envs/mtbseq/share/mtbseq-1.0.3-1/var/ref/M._tuberculosis_H37Rv_2015-11-13.fasta -o /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/Bam/5765_somelib.bam /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/Bam/5765_somelib.sam 2>> /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/Bam/5765_somelib.bamlog


<a id="org8ae18b1"></a>

## samtools sort - TBbwa

samtools for sortimg the bam file

    samtools sort -@ 2 -T /tmp/5765_somelib.sorted -o /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/Bam/5765_somelib.sorted.bam /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/Bam/5765_somelib.bam 2>> /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/Bam/5765_somelib.bamlog


<a id="org04f2c55"></a>

## samtools index - TBbwa

samtools for indexing the sorted bam

    samtools index -b /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/Bam/5765_somelib.sorted.bam 2>> /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/Bam/5765_somelib.bamlog


<a id="orge57e785"></a>

## samtools rmdup - TBbwa

samtools for removing putative pcr duplicates

    samtools rmdup /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/Bam/5765_somelib.sorted.bam /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/Bam/5765_somelib.nodup.bam 2>> /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/Bam/5765_somelib.bamlog


<a id="org94f77e7"></a>

## samtools index - TBbwa

recreating the index after removing the pcr duplicates

    samtools index -b /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/Bam/5765_somelib.nodup.bam 2>> /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/Bam/5765_somelib.bamlog


<a id="org214b6a8"></a>

## bwa mapping? - TBbwa

refining mappings:


<a id="orgd4bd11c"></a>

### GATK RealignerTargetCreator

    gatk --analysis_type RealignerTargetCreator --reference_sequence /root/miniconda3/envs/mtbseq/share/mtbseq-1.0.3-1/var/ref/M._tuberculosis_H37Rv_2015-11-13.fasta --input_file /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/Bam/5765_somelib.bam --downsample_to_coverage 10000 --num_threads 2 --out /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/GATK_Bam/5765_somelib.gatk.intervals 2>> /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/GATK_Bam/5765_somelib.gatk.bamlog


<a id="org5851d40"></a>

### GATK IndelRealigner

    gatk --analysis_type IndelRealigner --reference_sequence /root/miniconda3/envs/mtbseq/share/mtbseq-1.0.3-1/var/ref/M._tuberculosis_H37Rv_2015-11-13.fasta --input_file /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/Bam/5765_somelib.bam --defaultBaseQualities 12 --targetIntervals /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/GATK_Bam/5765_somelib.gatk.intervals --noOriginalAlignmentTags --out /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/GATK_Bam/5765_somelib.realigned.bam 2>> /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/GATK_Bam/5765_somelib.gatk.bamlog


<a id="orge7101fe"></a>

### GATK BaseRecalibrator

    gatk --analysis_type BaseRecalibrator --reference_sequence /root/miniconda3/envs/mtbseq/share/mtbseq-1.0.3-1/var/ref/M._tuberculosis_H37Rv_2015-11-13.fasta --input_file /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/GATK_Bam/5765_somelib.realigned.bam --knownSites /root/miniconda3/envs/mtbseq/share/mtbseq-1.0.3-1/var/res/MTB_Base_Calibration_List.vcf --maximum_cycle_value 600 --num_cpu_threads_per_data_thread 2 --out /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/GATK_Bam/5765_somelib.gatk.grp 2>> /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/GATK_Bam/5765_somelib.gatk.bamlog


<a id="orgaea9fea"></a>

### GATK PrintReads

    gatk -T --analysis_type PrintReads --reference_sequence /root/miniconda3/envs/mtbseq/share/mtbseq-1.0.3-1/var/ref/M._tuberculosis_H37Rv_2015-11-13.fasta --input_file /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/GATK_Bam/5765_somelib.realigned.bam --BQSR /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/GATK_Bam/5765_somelib.gatk.grp --num_cpu_threads_per_data_thread 2 --out /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/GATK_Bam/5765_somelib.gatk.bam  2>> /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/GATK_Bam/5765_somelib.gatk.bamlog


<a id="org36f0bc5"></a>

## mpileup creation - TBpile


<a id="org969b80a"></a>

### samtools mpileup

    samtools mpileup -B -A -x -f /root/miniconda3/envs/mtbseq/share/mtbseq-1.0.3-1/var/ref/M._tuberculosis_H37Rv_2015-11-13.fasta /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/GATK_Bam/5765_somelib.gatk.bam > /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/Mpileup/5765_somelib.gatk.mpileup 2>> /var/lib/docker/volumes/minikube/_data/lib/docker/volumes/nftowervolume/_data/data/mtbseqtest_conda/Mpileup/5765_somelib.gatk.mpileuplog


<a id="orgad724d5"></a>

## Creating position lists - TBlist


<a id="orgf67d99f"></a>

## Variant calling - TBVariants


<a id="orgb3ba568"></a>

## Refining variants - TBrefine


<a id="org67e9383"></a>

### gatk RealignerTargetCreator

    gatk --analysis_type RealignerTargetCreator --reference_sequence /usr/local/share/mtbseq-1.0.3-1/var/ref/M._tuberculosis_H37Rv_2015-11-13.fasta --input_file /home/minikubeuser/mtbseq-parallel-test/mtbseq-nf/work/d5/9dd48017a20a050eb8884b09488603/Bam/ERR841438_1_somelib.bam --downsample_to_coverage 10000 --num_threads 2 --out /home/minikubeuser/mtbseq-parallel-test/mtbseq-nf/work/d5/9dd48017a20a050eb8884b09488603/GATK_Bam/ERR841438_1_somelib.gatk.intervals 2>> /home/minikubeuser/mtbseq-parallel-test/mtbseq-nf/work/d5/9dd48017a20a050eb8884b09488603/GATK_Bam/ERR841438_1_somelib.gatk.bamlog


<a id="orgb083c6a"></a>

### gatk IndelRealigner

    
    gatk --analysis_type IndelRealigner --reference_sequence /usr/local/share/mtbseq-1.0.3-1/var/ref/M._tuberculosis_H37Rv_2015-11-13.fasta --input_file /home/minikubeuser/mtbseq-parallel-test/mtbseq-nf/work/d5/9dd48017a20a050eb8884b09488603/Bam/ERR841438_1_somelib.bam --defaultBaseQualities 12 --targetIntervals /home/minikubeuser/mtbseq-parallel-test/mtbseq-nf/work/d5/9dd48017a20a050eb8884b09488603/GATK_Bam/ERR841438_1_somelib.gatk.intervals --noOriginalAlignmentTags --out /home/minikubeuser/mtbseq-parallel-test/mtbseq-nf/work/d5/9dd48017a20a050eb8884b09488603/GATK_Bam/ERR841438_1_somelib.realigned.bam 2>> /home/minikubeuser/mtbseq-parallel-test/mtbseq-nf/work/d5/9dd48017a20a050eb8884b09488603/GATK_Bam/ERR841438_1_somelib.gatk.bamlog


<a id="org0f401c9"></a>

### gatk BaseRecalibrator

    gatk --analysis_type BaseRecalibrator --reference_sequence /usr/local/share/mtbseq-1.0.3-1/var/ref/M._tuberculosis_H37Rv_2015-11-13.fasta --input_file /home/minikubeuser/mtbseq-parallel-test/mtbseq-nf/work/d5/9dd48017a20a050eb8884b09488603/GATK_Bam/ERR841438_1_somelib.realigned.bam --knownSites /usr/local/share/mtbseq-1.0.3-1/var/res/MTB_Base_Calibration_List.vcf --maximum_cycle_value 600 --num_cpu_threads_per_data_thread 2 --out /home/minikubeuser/mtbseq-parallel-test/mtbseq-nf/work/d5/9dd48017a20a050eb8884b09488603/GATK_Bam/ERR841438_1_somelib.gatk.grp 2>> /home/minikubeuser/mtbseq-parallel-test/mtbseq-nf/work/d5/9dd48017a20a050eb8884b09488603/GATK_Bam/ERR841438_1_somelib.gatk.bamlog


<a id="orgcf088a9"></a>

### gatk PrintReads

    gatk -T --analysis_type PrintReads --reference_sequence /usr/local/share/mtbseq-1.0.3-1/var/ref/M._tuberculosis_H37Rv_2015-11-13.fasta --input_file /home/minikubeuser/mtbseq-parallel-test/mtbseq-nf/work/d5/9dd48017a20a050eb8884b09488603/GATK_Bam/ERR841438_1_somelib.realigned.bam --BQSR /home/minikubeuser/mtbseq-parallel-test/mtbseq-nf/work/d5/9dd48017a20a050eb8884b09488603/GATK_Bam/ERR841438_1_somelib.gatk.grp --num_cpu_threads_per_data_thread 2 --out /home/minikubeuser/mtbseq-parallel-test/mtbseq-nf/work/d5/9dd48017a20a050eb8884b09488603/GATK_Bam/ERR841438_1_somelib.gatk.bam  2>> /home/minikubeuser/mtbseq-parallel-test/mtbseq-nf/work/d5/9dd48017a20a050eb8884b09488603/GATK_Bam/ERR841438_1_somelib.gatk.bamlog


<a id="orgfb28061"></a>

## Statistics for mapping - TBstats


<a id="org607b5a0"></a>

### samtools


<a id="org3abd486"></a>

## Calling Strains - TBstrains

