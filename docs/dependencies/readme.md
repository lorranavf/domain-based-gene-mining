# GeneMining Environment Setup

Before proceeding with the GeneMining setup, ensure that Docker and Mamba are installed on your system.

You can find installation instructions at:
- [Docker](https://docs.docker.com/get-docker/)
- [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)

Then, create the mamba environment:

```bash
# create mamba env
mamba create -n genemining python=3.9
mamba activate genemining
```
Now, proceed with the installation of the required packages and tools for GeneMining:

First, download the following programs:
- [Deeploc2](https://services.healthtech.dtu.dk/services/DeepLoc-2.0/)
- [SignalP 6.0](https://services.healthtech.dtu.dk/services/SignalP-6.0/)

```bash
# install python packages
mamba install pandas biopython pybiolib cialign orthofinder ete3 seaborn multiqc

# from repo apt
sudo apt install mafft -y
sudo apt install iqtree -y
sudo apt install hmmer -y
sudo apt install ncbi-blast+ -y
sudo apt install emboss -y
sudo apt install busco -y

# install cath-resolve-hits
curl -O -L https://github.com/UCLOrengoGroup/cath-tools/releases/download/v0.16.10/cath-resolve-hits.ubuntu-20.04
sudo mv cath-resolve-hits.ubuntu-20.04 /usr/bin/cath-resolve-hits
chmod 755 /usr/bin/cath-resolve-hits

# install signalp6
tar -xzvf signalp6.fast.tar.gz
pip install signalp6_fast/signalp-6-package/
SIGNALP_DIR=$(python3 -c "import signalp; import os; print(os.path.dirname(signalp.__file__))" )
sudo cp -r signalp6_fast/signalp-6-package/models/* $SIGNALP_DIR/model_weights/

# install deeploc2
tar -xzvf deeploc-2.0.All.tar.gz 
pip install deeploc2_package/

```
This setup will create a dedicated environment named 'genemining' and install the necessary Python packages and tools for the GeneMining project.

You can also verify the accessibility of the expected programs by executing the following step:

```bash
python dependencies.py
```