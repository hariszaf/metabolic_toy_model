#!/bin/bash

# Update
sudo apt-get update

# Install Miniconda silently
MINICONDA_DIR="/opt/miniconda"
if [ ! -d "$MINICONDA_DIR" ]; then
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
    sudo bash /tmp/miniconda.sh -b -p "$MINICONDA_DIR"
    rm /tmp/miniconda.sh
fi

# Add conda to PATH
echo 'export PATH="$MINICONDA_DIR/bin:$PATH"' >> ~/.bashrc
export PATH="$MINICONDA_DIR/bin:$PATH"

# Initialize conda for bash
eval "$(conda shell.bash hook)"
conda init bash

# Accept conda terms of service
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

conda config --remove-key channels
conda config --add channels conda-forge
conda config --set channel_priority strict



# Get microbetag
git clone https://github.com/msysbio/microbetag.git
cd microbetag
sudo bash setup_environment.sh
pip install -e .
cd ..
