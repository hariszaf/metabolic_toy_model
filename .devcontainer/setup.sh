#!/bin/bash
set -e

sudo apt-get update
sudo apt-get install -y libsuitesparse-dev

python -m pip install --upgrade pip
python -m pip install -r requirements.txt

# Required for carveme
wget http://github.com/bbuchfink/diamond/releases/download/v2.1.25/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
sudo mv diamond /usr/local/bin/