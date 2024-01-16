# Domain-based gene mining for phylogenetic database construction

This repository contains two Python scripts: gene_mining.py, and dependencies.py:

## dependencies.py
This script checks the dependencies required for running a specific program. It verifies the presence of both Python libraries and external programs in the environment.

### Usage
To run the script, open a terminal or command prompt and navigate to the directory where the dependencies.py script is located. Then, use the following command:

```
python dependencies.py
```

## gene_mining.py
This script performs domain analysis by processing sequences of proteins, and filtering these sequences based on domain information. It also generates various outputs such as metadata, phylogenetic trees, and annotations.

The Class DomainAnalysis run two process:

- First, performs domain analysis by obtaining domains present in each sequence using the hmmscan program and cath-resolve hits.
- It then filters the specified domains and analyzes the filtered sequences using various programs such as pepstats, deeploc2, signalp6, DeepTMHMM, MAFFT, CIAlign, and iqtree2.

Finally, it generates a database containing the sequence data, phylogenetic trees, and annotation files for domains (functional, signal  peptides, and transmembrane) and subcellular localization, which can be used for tree annotation with the web software iTOL.

### Usage

The script can be executed using two strategies: either a domain-based selection of the target sequences or the other strategy, which involves pre-selection based first on similarity and then on domain. 

In the first strategy, it is necessary to enter a directory with protein files and a dictionary of these files. 

In the other strategy, it is necessary to specify the database of the target specie and a fasta file containing the proteins known as references for the gene family.

**Domain-Based Gene Mining**

First, save the protein FASTA files in the 'in.files.db' directory and the Pfam database in the 'in.pfam' directory. See below for instructions on how to prepare the Pfam database.
Then, create the script run.py with the following commands:

```python
from gene_mining import Parameters, DomainAnalysis

codes = {'Code_pattern_for_specie_one': ['Specie Abbreviation', 'Specie']}
domain = ['pfam_short_name']
outdir = 'outdir_of_the_full_analysis_based_on_domain'

parameters = Parameters(param_seq_dicio = codes,
                        param_seq_ext= '.faa',
                        param_seq_regex= '^\w*.\d',
                        param_seq_path = 'in.files.db',
                        param_pfam_in= 'in.pfam/Pfam-A.hmm',
                        param_pfam_out = 'out.pfam',
                        param_domain = domain,
                        param_domain_group = False,
                        param_outdir = outdir,
                        param_hmm_analysis = True,
                        param_full_analysis = True,
                        param_cpu = 4)

DomainAnalysis(parameters).run()
```
Note that the script run.py should be in the same directory as gene_mining.py.

After obtaining the domain analysis for the first time, setting `param_full_analysis=False`. 

Then you can run the analysis again, this time setting `param_hmm_analysis = False`, and obtain subsequent analyses for as many targets as desired. 

For example:

```python
codes = {'Code_pattern_for_specie_one': ['Specie Abbreviation', 'Specie1']}

metadata = Parameters(param_seq_dicio=codes,
                      param_full_analysis=False)

DomainAnalysis(metadata).run()

params = [['outdir_domain_one', ['Domain1']],
          ['outdir_domain_two', ['Domain2']]]

for outdir, domain in params:
   
    print(f'Output analysis: {outdir}')

    metadata = Parameters(param_seq_dicio=codes,
                          param_domain=domain,
                          param_outdir=outdir,
                          param_hmm_analysis=False)

DomainAnalysis(metadata).run()
```

**Similarity and Domain-Based Gene Mining**

First, save the protein FASTA files of the target specie in the 'in.files.blastp.db' directory, the FASTA file of the known proteins (references) in the 'in.pfiles.blastp.reference' directory, and the Pfam database in the 'in.pfam' directory. Refer below for instructions on how to prepare the Pfam database.
Next, create the 'run.py' script with the following commands:

```python
from gene_mining import Parameters, DomainAnalysis

domain = ['pfam_short_name']
outdir = 'outdir_of_the_full_analysis_based_on_domain'

metadata = Parameters(param_seq_ext = '.faa',
                      param_pfam_in = 'in.pfam/Pfam-A.hmm',
                      param_blastdb = 'in.files.blastp.db/blastdb.faa',
                      param_blast_reference = 'in.files.blastp.reference/reference.faa',
                      param_seq_path = 'in.files.db',
                      param_blast_out = 'out.blastp',
                      param_pfam_out = 'out.pfam',
                      param_outdir = outdir,
                      param_domain = domain,
                      param_domain_group = False,
                      param_blast_analysis=True,
                      param_cpu = 4)

DomainAnalysis(parameters).run()
```

**Here's an explanation for each of the parameters:**

- `param_seq_dicio`: A dictionary mapping assembly codes to their corresponding species names. It provides information about the relationship between assembly codes and species names.

- `param_seq_ext`: The file extension for the sequence files. It specifies the file format used for the input sequences.

- `param_seq_regex`: A regular expression pattern used to match and extract the assembly code from the sequence file names. It helps identify the assembly code within the file names.

- `param_seq_path`: The path to the directory containing the input sequence files. It indicates the location of the sequence files to be processed.

- `param_pfam_in`: The path to the Pfam HMM file used for domain analysis. It specifies the location of the Pfam HMM file to be used as input.

- `param_pfam_out`: The output directory for Pfam domain analysis results. It determines where the results of the Pfam domain analysis will be saved.

- `param_domain`: A list of domain names to filter the sequences during the analysis. It allows specifying specific domains of interest for filtering the sequences.

- `param_domain_group`: A boolean flag indicating whether the domain filtering should consider exact matches (`False`) or group matches (`True`). It controls whether domain matching should be strict or allow group matches.

- `param_outdir`: The output directory for the analysis results. It determines where the generated output files and directories will be saved.

- `param_blastdb`: It specifies the path of the database of the target specie used for performing blastp analysis.

- `param_blast_reference`: It indicates the path of the FASTA file containing known protein sequences that will be used as references in the blastp analysis.

- `param_blast_out`: It specifies the location where the table fmt6 format of the blastp analysis will be saved.

- `param_hmm_analysis`: A boolean flag indicating whether to perform HMM analysis. It determines whether the script will execute the HMM analysis step.

- `param_full_analysis`: A boolean flag indicating whether to perform the full analysis. It determines whether the script will execute the full analysis, including HMM analysis and subsequent steps.

- `param_blast_analysis`: A boolean flag indicating whether to perform the blastp analysis. It mandates that the `param_blastdb` and `param_blast_reference` must be included in the list of parameters.

- `param_cpu`: The number of CPUs to use for parallel processing. It specifies the desired level of parallelism for the analysis.

Adjust these parameters based on your specific needs and preferences before running the script run.py.

**To prepare the Pfam database, follow these steps:**

1. Download the [Pfam database](https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/) from the Pfam website.
2. Once downloaded, navigate to the terminal.
3. Uncompress the Pfam database file by running the following command:

   ```bash
   gunzip Pfam-A.hmm.gz
   ```
   This will extract the Pfam-A.hmm file.

4. Prepare the profile database for hmmscan using the following command:

   ```bash
   hmmpress Pfam-A.hmm
   ```
   This command will create additional files required for efficient searching with hmmscan.

The Pfam database is now ready to be used for domain analysis. Make sure to provide the correct path to the Pfam-A.hmm file as the value for the `param_pfam_in` parameter in your script.

Make sure to follow these steps carefully to properly prepare the Pfam database for your domain analysis.

# Requirements
To use the gene_mining.py script, you need to have the following programs installed:

- [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
- [HMMER (hmmscan)](http://hmmer.org/) 
- [CATH Tools (cath-resolve-hits)](https://cath-tools.readthedocs.io/en/latest/tools/cath-resolve-hits/)
- [EMBOSS (pepstats)](http://emboss.open-bio.org/rel/rel6/apps/pepstats.html)
- [deeploc2](https://services.healthtech.dtu.dk/services/DeepLoc-2.0/)
- [Signalp6](https://dtu.biolib.com/SignalP-6)
- [DeepTMHMM(to run locally, you need to install pybiolib, and Docker)](https://dtu.biolib.com/DeepTMHMM)
- [Docker](https://docs.docker.com/engine/install/ubuntu/#install-using-the-convenience-script)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [CiAlign](https://pypi.org/project/cialign/)
- [IQ-TREE2](http://www.iqtree.org/)
- [Orthofinder](https://davidemms.github.io/)
- [ete3](http://etetoolkit.org/docs/latest/index.html)

You also need to have Python 3.x installed, along with the following Python modules:

- [pybiolib](https://pypi.org/project/pybiolib/)
- [biopython](https://biopython.org/)
- [pandas](https://pandas.pydata.org/)
- [seaborn](https://seaborn.pydata.org/index.html)

## Attention

To use the gene_mining.py script, make sure to have the necessary programs installed.

Note that this script was developed for research purposes and may require modifications to adapt to different scenarios.

