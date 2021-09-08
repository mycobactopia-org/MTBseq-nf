# Mtbseq-nf

**NOTE**: This still a work in progress, the code is functional however the usage doc needs some love!

We would like to Thank the developers of MTBseq for putting in the intial effort!

Project with aim to create a nextflow wrapper for mtbseq workflow.

Heavily based on [ngs-fzb/MTBseq_source](https://github.com/ngs-fzb/MTBseq_source)

Contributions are warmly accepted!


# Improvements over the original code

- Reliance of bioconda for installing packages
- Reproducibility via containers 
- Deployability on a range of infratructure (cloud/on-prem clusters/local machine)


# License


The insipiration for this project itself gentb-snakemake doesn't have a license as of [db4d8fc637711f31544047ed43d9098923efd8f6](https://github.com/farhat-lab/gentb-snakemake/tree/db4d8fc637711f31544047ed43d9098923efd8f6)

However, the gentb-snakemake project builds upon other tools which have a variety of licenses https://github.com/farhat-lab/gentb-snakemake/search?q=LICENSE.

The components related to gentb-nf project itself (the Nextflow code) are licensed under the liberal MPL-2.0 license.

