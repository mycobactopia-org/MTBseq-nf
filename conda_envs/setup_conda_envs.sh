#!/usr/bin/env bash

set -xue

# NOTE: Please replace `conda` with `mamba` if it is installed for faster installs.

# NOTE: The conda environments are expected by the `conda_local` profile to be created within `conda_envs` directory

conda env create -p mtbseq-nf-env --file mtbseq-nf-env.yml
