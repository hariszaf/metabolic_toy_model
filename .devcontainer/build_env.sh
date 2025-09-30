

sudo apt-get update
# Requirement for dingo
sudo apt-get install -y libsuitesparse-dev

# In codesapces, the current working directory will be the /workspaces/metabolic_toy_model
pip install -r requirements.txt

# Install PolyRound 
CWD=$(pwd)
cd .. 
git clone https://gitlab.com/csb.ethz/PolyRound.git
cd PolyRound 
pip install -e .
cd CWD=$(pwd)

