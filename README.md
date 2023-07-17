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

- First, performs domain analysis by obtaining domains present in each sequence using the hmmscan program.
- It then filters the specified domains and analyzes the filtered sequences using various programs such as pepstats, deeploc2, signalp6, DeepTMHMM, MAFFT, CIAlign, and iqtree2.

Finally, it generates a database containing the sequence data, phylogenetic trees, and annotation files for domains and subcellular localization, which can be used for tree annotation with the web software iTOL.

### Usage

Create the script run.py with the following commands:

```python
from gene_mining import Parameters, DomainAnalysis

codes = {'Code_pattern_for_specie_one': ['Abbreviation', 'Specie1']}

parameters = Parameters(param_seq_dicio = codes,
                        param_seq_ext= '.faa',
                        param_seq_regex= '^\w*.\d',
                        param_seq_path = 'in.files',
                        param_pfam_in= 'in.pfam//Pfam-A.hmm',
                        param_pfam_out = 'out.pfam',
                        param_domain = None,
                        param_domain_group = False,
                        param_outdir = None,
                        param_cpu = 4,
                        param_hmm_analysis = True,
                        param_full_analysis = True)

DomainAnalysis(parameters).run()
```
Note that the script run.py should be in the same directory as gene_mining.py.

After obtaining the domain analysis for the first time, setting `param_full_analysis=False`. 
Then you can run the analysis again, this time setting `param_hmm_analysis = False`, and obtain subsequent analyses for as many targets as desired. 

For example:

```python
codes = {'Code_pattern_for_specie_one': ['Abbreviation', 'Specie1']}

params = [['domain_one', ['Domain1']],
          ['domain_two', ['Domain2']]

metadata = Parameters(param_seq_dicio=codes,
                      param_domain=domain,
                      param_outdir=outdir,
                      param_full_analysis=False)

DomainAnalysis(metadata).run()

for outdir, domain in params:
    print(f'Output analysis: {outdir}')

    metadata = Parameters(param_seq_dicio=codes,
                          param_domain=domain,
                          param_outdir=outdir,
                          param_hmm_analysis=False)

    DomainAnalysis(metadata).run()
```

Here's an explanation for each of the parameters:

- `param_seq_dicio`: A dictionary mapping assembly codes to their corresponding species names. It provides information about the relationship between assembly codes and species names.

- `param_seq_ext`: The file extension for the sequence files. It specifies the file format used for the input sequences.

- `param_seq_regex`: A regular expression pattern used to match and extract the assembly code from the sequence file names. It helps identify the assembly code within the file names.

- `param_seq_path`: The path to the directory containing the input sequence files. It indicates the location of the sequence files to be processed.

- `param_pfam_in`: The path to the Pfam HMM file used for domain analysis. It specifies the location of the Pfam HMM file to be used as input.

- `param_pfam_out`: The output directory for Pfam domain analysis results. It determines where the results of the Pfam domain analysis will be saved.

- `param_domain`: A list of domain names to filter the sequences during the analysis. It allows specifying specific domains of interest for filtering the sequences.

- `param_domain_group`: A boolean flag indicating whether the domain filtering should consider exact matches (`False`) or group matches (`True`). It controls whether domain matching should be strict or allow group matches.

- `param_outdir`: The output directory for the analysis results. It determines where the generated output files and directories will be saved.

- `param_cpu`: The number of CPUs to use for parallel processing. It specifies the desired level of parallelism for the analysis.

- `param_hmm_analysis`: A boolean flag indicating whether to perform HMM analysis. It determines whether the script will execute the HMM analysis step.

- `param_full_analysis`: A boolean flag indicating whether to perform the full analysis. It determines whether the script will execute the full analysis, including HMM analysis and subsequent steps.

Adjust these parameters based on your specific needs and preferences before running the script run.py.

# Requirements
To use the gene_mining.py script, you need to have the following programs installed:

- HMMER (hmmscan)
- EMBOSS (pepstats)
- deeploc2
- Signalp6
- DeepTMHMM (to run locally, you need to install Docker)
- MAFFT
- CiAlign
- IQ-TREE2

You also need to have Python 3.x installed, along with the following Python modules:

- biolib
- biopython
- pandas

## Attention

To use the gene_mining.py script, make sure to have the necessary programs installed.

Note that this script was developed for research purposes and may require modifications to adapt to different scenarios.
