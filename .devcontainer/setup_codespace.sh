#!/bin/bash

# Update
# -----
sudo apt-get update


# Install Miniconda
# -----
MINICONDA_DIR="$HOME/miniconda"
CONDA_BIN="$MINICONDA_DIR/bin/conda"

if [ ! -d "$MINICONDA_DIR" ]; then
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p "$MINICONDA_DIR"  # no sudo anymore
    rm /tmp/miniconda.sh
fi

# Add conda to PATH
echo 'export PATH="$MINICONDA_DIR/bin:$PATH"' >> ~/.bashrc
export PATH="$MINICONDA_DIR/bin:$PATH"

# Initialize conda for bash (without sudo)
eval "$($CONDA_BIN shell.bash hook)"


# Accept conda TOS
echo "Ensuring Anaconda TOS acceptance..."
$CONDA_BIN tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
$CONDA_BIN tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Force base to use Python 3.10
echo "Having Python 3.10 in base environment..."
$CONDA_BIN install -n base python=3.10 -y

# Accept conda terms of service
$CONDA_BIN config --system --set always_yes true
$CONDA_BIN config --system --set auto_update_conda false


# Get microbetag
# -----
git clone -b codespace https://github.com/hariszaf/microbetag.git --single-branch
cd microbetag
sudo bash setup_environment.sh
pip install -e .
cd ..
