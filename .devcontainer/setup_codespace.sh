#!/bin/bash

# Update
# -----
sudo apt-get update


# Install Miniconda
# -----
MINICONDA_DIR="/opt/miniconda"
CONDA_BIN="/opt/miniconda/bin/conda"

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
echo "Ensuring Anaconda TOS acceptance..."
$CONDA_BIN tos accept --all


# Get microbetag
# -----
# git clone https://github.com/hariszaf/microbetag.git
git clone -b codespace https://github.com/hariszaf/microbetag.git --single-branch
cd microbetag
sudo bash setup_environment.sh
pip install -e .
cd ..
