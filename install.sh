#!/usr/bin/env bash

#conda env create -n mtoolbox -f envs/mtoolbox.yaml

# create alias for env activation and running
#echo 'alias mtoolbox-activate="export PATH='`pwd`':'`pwd`'/scripts:$(conda run -n mtoolbox python -c "import os; print(os.environ.get(\"PATH\"))"); conda activate mtoolbox"' >> ~/.bash_profile
unalias mtoolbox-activate 2> /dev/null
echo 'alias mtoolbox-activate="export PATH='`pwd`':'`pwd`'/scripts:$PATH && conda activate mtoolbox"' >> ~/.bash_profile

