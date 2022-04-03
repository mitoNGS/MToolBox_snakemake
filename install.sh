#!/usr/bin/env bash

mamba env create -n mtoolbox -f envs/mtoolbox.yaml

# create alias for env activation and running
unalias mtoolbox-activate 2> /dev/null
echo 'alias mtoolbox-activate="export PATH='`pwd`':'`pwd`'/scripts:$PATH && conda activate mtoolbox"' >> ~/.bash_profile

