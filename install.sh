#!/bin/bash

# directory of script
nctc="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# will return error if env called nctc-tools already installed
conda env create -f ${nctc}/env/nctc_tools.yaml

# add script location to path
echo 'export PATH="''$PATH'':'$nctc'"' >> ~/.bashrc
source ~/.bashrc

# make script executable
chmod +x $nctc/nctc_tools.py
