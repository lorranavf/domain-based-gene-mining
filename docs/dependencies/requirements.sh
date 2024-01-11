#!/bin/sh

# from repo apt (ubuntu store)
apt install mafft -y
apt install iqtree -y
apt install hmmer -y
apt install ncbi-blast+ -y
apt install emboss -y
apt install docker.io -y

# from repo pip (python store)
pip install pandas
pip install biopython
pip install pybiolib
pip install cialign
pip install gdown

# from other sources
mkdir programs
cd programs

curl -O -L https://github.com/UCLOrengoGroup/cath-tools/releases/download/v0.16.10/cath-resolve-hits.ubuntu-20.04
mv cath-resolve-hits.ubuntu-20.04 /usr/bin/cath-resolve-hits
chmod 755 /usr/bin/cath-resolve-hits

gdown --id 1gbTfAi_i6HnBGkLtvbo6_nocKUltOeo9 --output signalp6.tar.gz
tar -xzvf signalp6.tar.gz
pip install signalp6_fast/signalp-6-package/
SIGNALP_DIR=$(python3 -c "import signalp; import os; print(os.path.dirname(signalp.__file__))" )
cp -r signalp6_fast/signalp-6-package/models/* $SIGNALP_DIR/model_weights/

gdown --id 1msW1RpyZHW9pUOANH4phI6tYenh4gQV5 --output deeploc2.tar.gz
tar -xzvf deeploc2.tar.gz
pip install deeploc2_package/

cd ..
rm -r programs

# another settings
# pip install torch torchvision torchaudio

# nano /workspace/.pyenv_mirror/user/current/lib/python3.11/site-packages/lightning_fabric/utilities/cloud_io.py
# editar linha 52
# return torch.load(f, map_location=torch.device('cpu'))  # type: ignore[arg-type]

# docker 
# https://docs.docker.com/engine/install/ubuntu/#install-using-the-convenience-script
# curl -fsSL https://get.docker.com -o get-docker.sh
# sh get-docker.sh