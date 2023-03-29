# bioanalysis

This repository contains two Python scripts: filo.py and config.py.

### filo.py
The filo.py script contains two functions:

###### pfam
This function receives a list of fasta files and uses the hmmscan program to generate a database for each species/input file of the domains present in each sequence contained in the files.

###### analysis
This function receives the database generated by the pfam function and filters the sequences that contain the domain specified by the user. Then, the filtered sequences are analyzed by the pepstats program from the emboss package, deeploc2, mafft, CiAlign and iqtree2. At the end, a database is generated containing the sequence data, such as species, genome assembly, sequence ID, sequence length, isoelectric point, subcellular localization, and also a phylogenetic tree separated in another file.

### config.py
This script contains the configuration variables used by the filo.py script, such as the paths to the input files, the domain to be filtered, and the paths to the necessary programs.

###### Pay attention

To use the filo.py script, make sure to have the necessary programs installed. Then, run the config.py script by passing the variables that are necessary.

Note that this script was developed for research purposes and may require modifications to adapt to different scenarios.
