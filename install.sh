#!/usr/bin/env bash

# conda env create -n mtoolbox -f envs/mtoolbox.yaml

# create alias for env activation and running
echo 'alias mtoolbox-activate="export PATH='`pwd`':'`pwd`'/scripts:$(conda run -n mtoolbox echo $PATH); conda activate mtoolbox"' >> ~/.bash_profile
