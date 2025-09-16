#!/bin/bash

# Update
# -----
sudo apt-get update


# Install Miniconda
# -----
MINICONDA_DIR="/opt/miniconda"
CONDA_BIN="$MINICONDA_DIR/bin/conda"

if [ ! -d "$MINICONDA_DIR" ]; then
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
    sudo bash /tmp/miniconda.sh -b -p "$MINICONDA_DIR"
    rm /tmp/miniconda.sh
fi

# Add conda to PATH
echo 'export PATH="$MINICONDA_DIR/bin:$PATH"' >> ~/.bashrc
export PATH="$MINICONDA_DIR/bin:$PATH"

# Initialize conda for bash (without sudo)
eval "$($CONDA_BIN shell.bash hook)"

# Force base to use Python 3.10
$CONDA_BIN install -n base python=3.10 -y

# Accept conda terms of service
echo "Ensuring Anaconda TOS acceptance..."
$CONDA_BIN config --system --set always_yes true
$CONDA_BIN config --system --set auto_update_conda false


# Get microbetag
# -----
git clone -b codespace https://github.com/hariszaf/microbetag.git --single-branch
cd microbetag
sudo bash setup_environment.sh
pip install -e .
cd ..
