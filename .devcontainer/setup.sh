#!/bin/bash
set -e

sudo apt-get update
sudo apt-get install -y libsuitesparse-dev

python -m pip install --upgrade pip
python -m pip install -r requirements.txt