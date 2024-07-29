#!/usr/bin/env bash

source $(conda info --base)/etc/profile.d/conda.sh
conda activate mchip-c

./Helpers/MChIPC_interaction_calling.R
