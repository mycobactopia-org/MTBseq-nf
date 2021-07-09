All credits to this manual goes to the developers of MTBseq, please check the original [MANUAL](https://github.com/ngs-fzb/MTBseq_source/blob/master/MANUAL.md)

## [OPTIONS & VALUES]

Using OPTIONS and their respective VALUES, a user can modify the execution of MTBseq by directly calling
functional modules and supplying parameters used in the analysis.


 [OPTIONS]
     [VALUES]  
 
### Steps
 **--step**  The OPTION **--step**> is essential and requires a user specified VALUE. MTBseq has a modular architecture with functional modules encapsulating distinct steps in the analysis pipeline. Using the **--step** OPTION, you can specify whether the full pipeline is executed, call a specific pipeline step only or as starting point for the pipeline. During execution, MTBseq will only create output files not already present in the working environment. If you choose a specific pipeline step, make sure that the output files of previous steps needed as input are present. In the following, all possible VALUES for this OPTION are described:
 #### TBfull
Description: 
      For executing the whole pipeline. If you choose this VALUE, make sure that you specify all analysis parameters for which you do not want to use the preset defaults using the appropriate OPTIONS.    
 
####   TBbwa
Description:
For mapping next generation sequencing read files to a reference genome, using BWA mem. This is the first step within the pipeline. Depending on the library type (single-end or paired-end), this module will align single end or paired end FASTQ files to a reference genome. Please ensure that FASTQ files used as input follow the required naming scheme described above.

By default, MTBseq uses the Mycobacterium tuberculosis H37Rv genome (NC_000962.3) as a reference and maps reads with the BWA mem program with default settings. After mapping, files (.sam) are converted into binary mapping files (.bam). The mapping is sorted, indexed and putative PCR duplicates are removed, using the program SAMTOOLS. 

 Wherever possible, multi-threading is activated by the user provided VALUES for the **--threads** OPTION.
 
                 Input:
                 [SampleID]_[LibID]_[*]_[Direction].fastq.gz
                 Output:
                 Bam/[SampleID]_[LibID]_[*].bam
                 Bam/[SampleID]_[LibID]_[*].bai
                 Bam/[SampleID]_[LibID]_[*].bamlog
#### TBrefine
Description:
for realignment around insertions and deletions and base call recalibration, using the program GATK. The GATK program is invoked with the parameters:
                 `--downsample_to_coverage 10000 --defaultBaseQualities 12 --maximum_cycle_value 600 --noOriginalAlignmentTags`
                 
For base call recalibration, MTBseq employs a set of known resistance associated variants if _Mycobacterium tuberculosis_ H37Rv was used as reference genome. The calibration list is stored in the directory "var/res/MTB_Base_Calibration_List.vcf" of the package. For other reference genomes, the file needs to be specified with the **--basecalib** OPTION or this step will be skipped.

                 Input:
                 Bam/[SampleID]_[LibID]_[*].bam
                 Output:
                 GATK_Bam/[SampleID]_[LibID]_[*].gatk.bam
                 GATK_Bam/[SampleID]_[LibID]_[*].gatk.bai
                 GATK_Bam/[SampleID]_[LibID]_[*].gatk.bamlog
                 GATK_Bam/[SampleID]_[LibID]_[*].gatk.grp
                 GATK_Bam/[SampleID]_[LibID]_[*].gatk.intervals
#### TBpile
Description:
For creating pileup files (.mpileup) from refined mapping files (.gatk.bam), using the program SAMTOOLS. The SAMTOOLS program is executed with the parameters - B - A

                 Input:
                 GATK_Bam/[SampleID]_[LibID]_[*].gatk.bam
                 Output:
                 Mpileup/[SampleID]_[LibID]_[*].gatk.mpileup
                 Mpileup/[SampleID]_[LibID]_[*].gatk.mpileuplog
                 
#### TBlist
Description:
For creating position lists from pileup files (.gatk.mpileup). Position list files capture the essential data from a mapping in a table format with 21 columns, containing the nucleotide counts in mapped reads for each position of the reference genome. This step will take into account the **--threads** and the **--minbqual**> OPTIONS. The columns in the position list file are listed in table 1.

                 Input:
                 Mpileup/[SampleID]_[LibID]_[*].gatk.mpileup
                 Output:
                 Position_Tables/[SampleID]_[LibID]_[*].gatk_position_table.tab

#### TBvariants
Description: 
For variant detection from position lists. Variant calling can be run in different detection modes using the OPTIONS **--all_vars**, **--snp_vars**, and **--lowfreq_vars** (with the respective changes to the detection algorithm explained in detail in the description for these OPTIONS). In the default mode, a variant is called if it is supported by a minimum number of 4 reads in both forward and reverse direction, by a minimum number of 4 reads indicating the allele with a phred score of at least 20, and if the allele is indicated by a minimum frequency of 75% in mapped reads. These thresholds can all be set using the respective modifiers ( **--mincovf**, **--mincovr**, **--minphred20**, **--minfreq**). In addition, positions fulfilling these thresholds will be counted as 'unambigous' in the calculation of quality values for the dataset. The parameters used by the pipeline will be recorded in the file names by dedicated fields, and the detection mode will be recorded in the output files as a binary string (e.g. "outmode100" indicates **--all_vars** but not **--snp_vars** and **--lowfreq_vars** active).

        Input:
        Position_Tables/[SampleID]_[LibID]_[*].gatk_position_table.tab
        Output:
        Called/[SampleID]_[LibID]_[*].gatk_position_uncovered_[mincovf]_[mincovr]_[minfreq]_[minphred20]_[all_vars][snp_vars][lowfreq_vars].tab
        Called/[SampleID]_[LibID]_[*].gatk_position_variants_[mincovf]_[mincovr]_[minfreq]_[minphred20]_[all_vars][snp_vars][lowfreq_vars].tab

#### TBstats
Description:
For calculating an overview of mapping quality and detected variants for a dataset, using the SAMTOOLS flagstat program. This step creates or updates a tabular delimited file "Mapping_and_Variant_Statistics.tab". This file stores all sample statistics for the analyzed datasets present in the working environment. This module employs the same thresholds to discern unambiguously covered positions as the TBvariants module (set by the **--mincovf**, **--mincovr**, **--minphred20**, **--minfreq**).The columns of the output are shown in table 2.*

        Input:**
        Bam/[SampleID]_[LibID]_[*].bam
        Position_Tables/[SampleID]_[LibID]_[*].gatk_position_table.tab
        Output:
        Statistics/Mapping_and_Variant_Statistics.tab
        
#### TBstrains
Description:
For lineage classification based on a set of phylogenetic SNPs (Homolka et al., 2012; Coll et al., 2014; Merker et al., 2015). This module creates a tabular delimited file within the "Classification" directory. Within this file, the majority lineage is reported for each dataset. This entry also gives an indication of the data quality for the positions used to infer the phylogenetic classification. The column indicating the quality will contain the label _good_ if all phylogenetic positions contained in the set are covered at least 10fold and show a frequency of at least 75%, _bad_ if any phylogenetic position does not meet these, and _ugly_ if any phylogenetic position does not have a clear base call.

        Input:
        Position_Tables/[SampleID]_[LibID]_[*].gatk_position_table.tab
        Output:
        Classification/Strain_Classification.tab
#### TBjoin
Description:
For creating a comparative SNP analysis of a set of samples specified by the user with the **--samples** OPTION. The comparative analysis can be run in different detection modes using the OPTIONS **--snp_vars**, and **--lowfreq_vars** (with the respective changes to the detection algorithm explained in detail in the description for these OPTIONS). For the joint analysis, first a scaffold of all variant positions is built from the individual variant files. Second, for all positions with a variant detected for any of the input samples, the allele information is recalculated from the original mappings to produce a comprehensive table. Output files will be created starting with the project name specified with the **--project** OPTION, otherwise the default file name will start with “NONE”.

The first line within the tabular delimited file joint variant table is a sample header line. The second line starts with position specific columns and continuing with sample specific entries shown in table 3.
         
        Input:
        samples.txt
        --project
        Called/[SampleID]_[LibID]_[*].gatk_position_variants_[mincovf]_[mincovr]_[minfreq]_[minphred20]_[all_vars][snp_vars][lowfreq_vars].tab
        Position_Tables/[SampleID]_[LibID]_[*].gatk_position_table.tab
        Output
        Joint/[PROJECT]_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples.tab
        Joint/[PROJECT]_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples.log
#### TBamend

Description:
For post-processing of joint variant tables. If you call this step directly, you need to specify the correct **--project** OPTION. This step will produce a comprehensive variant table including calculated summary values for each position. In addition, the set of positions will be processed in consecutive filtering steps. In the first step, positions will be excluded if less than a minimum percentage of samples have unambiguous base calls at this position, with this threshold set by the **--unambig** OPTION.
In addition, all samples need to have either a SNP or wild type base at the position, and positions within repetitive regions of the reference genome or within a resistance associated genes are excluded.
This filtering step results in output files carrying the "\_amended_[unambig]\_phylo" ending; a full table (ending in .tab), a FASTA file containing the aligned alleles of all samples for the given position (.fasta), and a corresponding FASTA file with the headers consisting solely of the respective sample ID (\_plainIDs.fasta).
The second filtering step removes positions that are located within a maximum distance to each other in the same sample, with this threshold set by the **--window** OPTION. 
Output files from this filtering step are generated with the same naming scheme as for the first step and carry the selected window threshold as additional field. In addition, positions not passing the window criteria are reported in the output file carrying the tag "removed".

        Input
        samples.txt
        --project           
        Joint/[PROJECT]_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples.tab
        Output
        Amend/[PROJECT]_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended.tab
        Amend/[PROJECT]_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended_[unambig]_phylo.tab
        Amend/[PROJECT]_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended_[unambig]_phylo.fasta
        Amend/[PROJECT]_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended_[unambig]_phylo_plainIDs.fasta
        Amend/[PROJECT]_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended_[unambig]_phylo_[window].tab
        Amend/[PROJECT]_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended_[unambig]_phylo_[window].fasta
        Amend/[PROJECT]_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended_[unambig]_phylo_[window]_plainIDs.fasta
        Amend/[PROJECT]_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended_[unambig]_phylo_[window]_removed.tab

#### TBgroups
Description: for inferring likely related isolates based on pairwise distance of distinct SNP positions. In an agglomerative process, samples are grouped together if they are within the threshold set by the **--distance** OPTION. The output files consist of a text file listing the detected groups and ungrouped isolates, and the calculated pairwise distance matrix.

        Input:
        Amend/[PROJECT]_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended_[unambig]_phylo_[window].tab
        Output:
        Groups/[PROJECT]_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended_[unambig]_phylo_[window].matrix
        Groups/[PROJECT]_joint_[mincovf]_[mincovr]_[minfreq]_[minphred20]_samples_amended_[unambig]_phylo_[window]_[distance].groups

#### Aditional configuration

 **--continue** If a module was chosen with the **--step** OPTION, the **--continue** OPTION ensures that the pipeline will continue the analysis with downstream modules. This is automatically set if the **--step** OPTION is set to the VALUE **TBfull**.

 **--samples** This OPTION requires a user supplied file specifying a set of datasets (e.g. samples.txt) as VALUE. The file must be a two-column, tab-separated file. Column 1 has to be your **[SampleID]**. Column 2 has to be your **[LibID]**. **TBjoin** requires this OPTION to be set.

 **--project** This OPTION allows the user to set a project name for the steps **TBjoin**, **TBamend** and **TBgroups**. If you do not support a project name, **[NONE]** is used as a default value.

 **--ref** This OPTION sets the reference genome for the read mapping. By default, the genome of <i>Mycobacterium tuberculosis</i> H37Rv (NC_000962.3) is set as reference. User supplied FASTA files for other reference genomes should be placed in the directory /MTBseq_source/var/ref/, and the respective name given without .fasta extension. Please be aware that for other reference genomes, you need to provide the respective annotation files as well or annotations will be skipped.

 **--resilist** This OPTION sets a list of known variant positions associated to drug resistance for resistance prediction. Give the full path to the file. The required structure of the file can be seen here: /MTBseq_source/var/res/MTB_Resistance_Mediating.txt
 
 **--intregions** This OPTION sets a list of interesting regions to be used for annotation of detected variants. Give the full path to the file. The required structure of the file can be seen here: /MTBseq_source/var/res/MTB_Extended_Resistance_Mediating.txt

 **--categories** This OPTION specifies a gene categories file to annotate essential and non-essential genes as well as repetitive regions. SNPs in repetitive regions will be excluded for phylogenetic analysis. Give the full path to the file. The required structure of the file can be seen here: /MTBseq_source/var/cat/MTB_Gene_Categories.txt

 **--basecalib** This OPTION specifies a file for base quality recalibration. The list must be in VCF format and should contain known SNPs. Give the full path to the file. The required structure of the file can be seen here: /MTBseq_source/var/res/MTB_Base_Calibration_List.vcf

 **--all_vars** This OPTION is used in the modules **TBvariants**, **TBstats**, **TBjoin**, and **TBstrains**. By default, the OPTION is not active. Setting this OPTION will skip all filtering steps and report the calculated information for all positions in the input file.

 **--snp_vars** This OPTION is used in **TBvariants**, **TBstats**, **TBjoin**, and **TBstrains**. By default, the OPTION is not active. Setting this OPTION will add an additional filter that excludes all variants except SNPs.

 **--lowfreq_vars** This OPTION is used in **TBvariants**, **TBstats**, **TBjoin**, and **TBstrains**. By default, the OPTION is not active. Setting this OPTION has major implications on how the mapping data for each position is processed. By default, the majority allele is called and taken for further calculations. If the **--lowfreq_vars** OPTION is set, MTBseq will consider the majority allele distinct from wild type, if such an allele is present. This means that only in this detection mode, MTBseq will report variants present only in subpopulations, i.e. low frequency mutations. Of course, OPTIONS **--mincovf**, **--mincovr**, **--minphred20**, and **--minfreq** need to be set accordingly. Please be aware that output generated in this detection mode should not be used for phylogenetic analysis.

 **--minbqual** This OPTION is used in **TBlist**. By default, the OPTION is set to 13. The OPTION sets a threshold for the sequence data quality to be used for the mpileup creation.

 **--mincovf** This OPTION is used in **TBvariants**, **TBjoin**, **TBamend**, and **TBstrains**. By default, the OPTION is set to 4. The OPTION sets a minimum forward read coverage threshold. Alleles must have a forward coverage of this VALUE or higher to be considered.

 **--mincovr** This OPTION is used in **TBvariants**, **TBjoin**, **TBamend**, and **TBstrains**. By default, the OPTION is set to 4. The OPTION sets a minimum reverse read coverage threshold. Alleles must have a reverse coverage of this VALUE or higher to be considered.


 **--minphred** This OPTION is used in **TBvariants**, **TBjoin**, **TBamend**, and **TBstrains**. By default, the OPTION is set to 4. The OPTION sets a minimum number of reads indicating an allele with a phred score of at least 20.

 **--minfreq** This OPTION is used in **TBvariants**, **TBjoin**, **TBamend**, and **TBstrains**. By default, the OPTION is set to 75. The OPTION sets a minimum frequency for an allele.

 **--unambig** This OPTION is used in **TBamend**. By default, the OPTION is set to 95. The option sets a minimum percentage of samples with unambiguous information for position.

 **--window** This OPTION is used in **TBamend**. By default, the OPTION is set to 12. The OPTION sets a window size in which the algorithm scans for the occurrence of multiple variants within the same sample. If more than one variant occurs within this window in the same sample, the positions will be excluded.

 **--distance** This OPTION is used in **TBgroups**. By default, the OPTION is set to 12. The OPTION sets a SNP distance that is used to classify samples into groups of samples, using agglomerative clustering. If SNP distances between samples are less or equal this VALUE, they are grouped together.

 **--quiet** This OPTION turns off the display logging function and will report the logging only in a file, called "MTBseq_[DATE]_[USER].log".

 **--threads** This OPTION is used in **TBbwa**, **TBmerge**, **TBrefine**, **TBpile** and **TBlist**. By default, the OPTION is set to 1. The OPTION sets the maximum number of CPUs to use within the pipeline. You can use more than one core in order to execute the pipeline faster. 8 is the current maximum.

 **--help** This OPTION will show you all available OPTIONs and corresponding VALUEs used by MTBseq.

 **--version** This OPTION will show you the current version of MTBseq.
 **--check** This OPTION will check the dependencies of MTBseq.
