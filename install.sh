#!/usr/bin/env bash

conda env create -n mtoolbox -f envs/mtoolbox.yaml

# install bamUtils
conda activate mtoolbox
git clone https://github.com/statgen/bamUtil.git
cd bamUtil
make cloneLib
make
make install INSTALLDIR=$(dirname $(which python))
conda deactivate

# create alias for env activation and running
echo 'alias mtoolbox_activate="conda activate mtoolbox; export PATH='`pwd`':$PATH; export PATH='`pwd`'/scripts:$PATH"' >> ~/.bash_profile
