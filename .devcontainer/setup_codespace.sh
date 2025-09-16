#!/bin/bash

# Update
# -----
sudo apt-get update


# Install Miniconda
# -----
MINICONDA_DIR="$HOME/miniconda"
CONDA_BIN="$MINICONDA_DIR/bin/conda"

# --- Install Miniconda if missing ---
if [ ! -d "$MINICONDA_DIR" ]; then
    echo "Installing Miniconda..."
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p "$MINICONDA_DIR"
    rm /tmp/miniconda.sh
fi

# --- Make conda available in the current script ---
export PATH="$MINICONDA_DIR/bin:$PATH"
eval "$($CONDA_BIN shell.bash hook)"

# --- Accept TOS for main and R channels ---
$CONDA_BIN tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
$CONDA_BIN tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# --- Force Python 3.10 in base environment ---
$CONDA_BIN install -n base python=3.10 -y

# --- Make conda available in all new terminals ---
BASHRC_LINE="export PATH=\"$MINICONDA_DIR/bin:\$PATH\""
CONDA_SH_LINE="source $MINICONDA_DIR/etc/profile.d/conda.sh"
BASE_ACTIVATE_LINE="conda activate base"

grep -qxF "$BASHRC_LINE" ~/.bashrc || echo "$BASHRC_LINE" >> ~/.bashrc
grep -qxF "$CONDA_SH_LINE" ~/.bashrc || echo "$CONDA_SH_LINE" >> ~/.bashrc
grep -qxF "$BASE_ACTIVATE_LINE" ~/.bashrc || echo "$BASE_ACTIVATE_LINE" >> ~/.bashrc

# --- Set user-level defaults ---
$CONDA_BIN config --set always_yes true
$CONDA_BIN config --set auto_update_conda false

echo -e "-- Miniconda setup on the Codespace has been completed. \U0001F389 \n\n"


# Get microbetag
# -----
git clone -b codespace https://github.com/hariszaf/microbetag.git --single-branch
cd microbetag
bash setup_environment.sh --phenotrex # no sudo
# pip install -e .
cd ..
