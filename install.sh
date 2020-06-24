#!/usr/bin/env bash

# conda env create -n mtoolbox -f envs/mtoolbox.yaml

conda activate mtoolbox

# install required packages not available through conda
pip install mtoolnote

conda deactivate

# create alias for env activation and running
echo 'alias mtoolbox_activate="conda activate mtoolbox; export PATH='`pwd`':$PATH; export PATH='`pwd`'/scripts:$PATH"' >> ~/.bash_profile
