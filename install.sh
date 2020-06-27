#!/usr/bin/env bash

conda env create -n mtoolbox -f envs/mtoolbox.yaml

# create alias for env activation and running
echo 'alias mtoolbox-activate="conda activate mtoolbox; export PATH='`pwd`':'`pwd`'/scripts:$(conda run -n mtoolbox echo $PATH)"' >> ~/.bash_profile
