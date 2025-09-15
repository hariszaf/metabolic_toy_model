#!/bin/bash

echo "hello friend"

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

conda activate base