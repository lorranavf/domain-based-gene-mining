import biolib
import os
import pandas as pd
import random
import re
import shutil
import textwrap
import time
import warnings
from Bio import SeqIO
from dataclasses import dataclass
from typing import List, Pattern

warnings.filterwarnings("ignore")


@dataclass
class Parameters:
    """
    """
    param_seq_dicio: dict
    param_seq_ext: str = '.faa'
    param_seq_regex: Pattern[str] = '^\w*.\d'
    param_seq_path: str = 'in.files'
    param_pfam_in: str = 'in.pfam//Pfam-A.hmm'
    param_pfam_out: str = 'out.pfam'
    param_domain: List[str] = None
    param_domain_group: bool = False
    param_outdir: str = None
    param_cpu: int = 4
    param_hmm_analysis: bool = True
    param_full_analysis: bool = True


class DomainAnalysis:
    """
    """

    def __init__(self, parameters: Parameters):
        self.parameters = parameters

    def run(self):

        if self.parameters.param_hmm_analysis:
            self.create_output_directory(analysis='hmm_analysis')
            self.process_sequences()

        if self.parameters.param_full_analysis:
            dirs = self.create_output_directory(analysis='full_analysis')
            self.filter_sequences_per_domain(dirs)

    def create_output_directory(self, analysis):

        if analysis == 'hmm_analysis':
            if os.path.exists(self.parameters.param_pfam_out):
                shutil.rmtree(self.parameters.param_pfam_out)
            os.mkdir(self.parameters.param_pfam_out)

        if analysis == 'full_analysis':
            outdir = f'out.{self.parameters.param_outdir}'
            dirs = [outdir,
                    f'{outdir}/metadados',
                    f'{outdir}/pfam',
                    f'{outdir}/fasta',
                    f'{outdir}/itol',
                    f'{outdir}/pepstats',
                    f'{outdir}/deeploc',
                    f'{outdir}/filogeny',
                    f'{outdir}/filogeny/mafft',
                    f'{outdir}/filogeny/cialign',
                    f'{outdir}/filogeny/iqtree',
                    f'{outdir}/log',
                    f'{outdir}/signalp',
                    f'{outdir}/deep_tmhmm']

            for directory in dirs:
                if os.path.exists(directory):
                    shutil.rmtree(directory)
                os.mkdir(directory)

            return dirs

    def process_sequences(self):
        regex = re.compile(self.parameters.param_seq_regex)
        ext = self.parameters.param_seq_ext
        files = os.listdir(self.parameters.param_seq_path)
        files = [filename for filename in files if os.path.splitext(filename)[1] == ext]

        for i in range(len(files)):
            assembly = regex.search(files[i])[0]
            seq = f'{self.parameters.param_seq_path}/{files[i]}'
            self.get_domains(assembly, seq)

    def get_domains(self, assembly, seq):

        print(f'Analyzing sample: {assembly}, {time.strftime("%H:%M %d/%m/%Y", time.localtime(time.time()))}')
        infile = self.parameters.param_pfam_in
        outfile = f'{self.parameters.param_pfam_out}/{assembly}'
        cpu = self.parameters.param_cpu
        os.system(f'hmmscan --cpu {cpu} -E 1e-5 -o {outfile}.out --domtblout {outfile}.pfam {infile} {seq}')

    def filter_sequences_per_domain(self, dirs):
        regex = re.compile(self.parameters.param_seq_regex)
        ext = self.parameters.param_seq_ext
        files = os.listdir(self.parameters.param_seq_path)
        files = [filename for filename in files if os.path.splitext(filename)[1] == ext]
        dfs_pfam = []
        dfs_itol = []
        dbs_faa = []

        for i in range(len(files)):
            domain = self.parameters.param_domain
            domain_group = self.parameters.param_domain_group
            dicio = self.parameters.param_seq_dicio
            assembly = regex.search(files[i])[0]
            seq = f'{self.parameters.param_seq_path}/{files[i]}'
            name_ab = next((v[0] for k,v in self.parameters.param_seq_dicio.items() if k == assembly), None)

            # pfamm hmmscan output
            table = f'{self.parameters.param_pfam_out}/{assembly}.pfam'
            table = pd.read_table(table, header=None, skiprows=3, skipfooter=10, engine='python')

            table_profile = pd.DataFrame({'Query': [i.split()[3] for i in table[0]],
                                          'Domain': [i.split()[0] for i in table[0]]})

            table_profile = table_profile.groupby('Query')['Domain'].apply(list).reset_index(name='Domains')
            table_profile['Assembly'] = assembly
            table_profile['Specie'] = next((v[1] for k, v in dicio.items() if k == assembly), None)
            table_profile = table_profile[['Specie', 'Assembly', 'Query', 'Domains']]

            print(f'Selecting domain in: {assembly}, {time.strftime("%H:%M %d/%m/%Y", time.localtime(time.time()))}')

            if domain_group:
                df_pfam_selected = table_profile[table_profile['Domains'] == str(domain)]
            else:
                df_pfam_selected = table_profile[table_profile['Domains'].apply(lambda x: any(d in x for d in domain))]

            dfs_pfam.append(df_pfam_selected)
            seq_to_select = list(df_pfam_selected['Query'])

            # iTOL domain anotation
            df_itol = pd.DataFrame({'query': [i.split()[3] for i in table[0]],
                                    'evalue': [i.split()[6] for i in table[0]],
                                    'from': [i.split()[17] for i in table[0]],
                                    'to': [i.split()[18] for i in table[0]],
                                    'domain': [i.split()[0] for i in table[0]]})
            df_itol = df_itol[df_itol['query'].isin(seq_to_select)]
            df_itol['query'] = df_itol['query'].apply(lambda x: f'{name_ab}_{x}')
            dfs_itol.append(df_itol)

            # faa database
            database = [record for record in SeqIO.parse(seq, 'fasta')]
            database_selected = [record for record in database if record.id in seq_to_select]
            for seq in database_selected:
                seq.id = f'{name_ab}_{seq.id}'
                seq.name = ''
                seq.description = ''
            dbs_faa.append(database_selected)

        df_pfam_selected_united = pd.concat(dfs_pfam)
        df_pfam_selected_united.to_csv(f'{dirs[2]}/pfam.csv', index=False)

        df_selected_itol = pd.concat(dfs_itol)
        df_selected_itol.to_csv(f'{dirs[2]}/pfam.itol', index=False)

        database_selected_united = [record for database in dbs_faa for record in database]
        SeqIO.write(database_selected_united, f'{dirs[3]}/database_selected_united.fasta', 'fasta')

        db_path = f'{dirs[3]}/database_selected_united.fasta'

        df_pepstats = self.get_pepstats(dirs, db_path)
        df_deeploc = self.get_deeploc(dirs, db_path)
        self.get_signalp(dirs, db_path)
        self.get_deep_tmhmm(dirs, db_path)
        self.get_metadados(dirs, df_pepstats, df_deeploc, df_pfam_selected_united)
        self.get_filogeny(dirs, db_path)
        self.get_domain_anot(dirs, db_path)
        self.get_local_anot(dirs, df_deeploc)

    def get_signalp(self, dirs, db_path):
        print('Calculating signal peptides')
        os.system(f'signalp6 --fastafile {db_path} --organism eukarya --output_dir {dirs[12]} --format txt --mode fast > {dirs[11]}/signalp.log 2>&1')

    def get_deep_tmhmm(self, dirs, db_path):
        print('Calculating transmembrane domains')
        deeptmhmm = biolib.load('DTU/DeepTMHMM')
        deeptmhmm_job = deeptmhmm.cli(args=f'--fasta {db_path}')
        deeptmhmm_job.save_files(dirs[13])

    def get_pepstats(self, dirs, db_path):
        print('Calculating protein stats')
        out_stats = f'{dirs[5]}/pepstats.out'
        os.system(f'pepstats -sequence {db_path} -outfile {out_stats} > {dirs[11]}/pepstats.log 2>&1')
        lines = open(out_stats, "r").readlines()
        pi = [float(line.replace("Isoelectric Point = ", "")) for line in lines if line.startswith("Isoelectric Point")]
        pm = [float(line.replace("Molecular weight = ", "").split()[0]) for line in lines if
              line.startswith("Molecular weight")]
        protein_id = [record.id for record in SeqIO.parse(db_path, 'fasta')]
        protein_length = [len(record.seq) for record in SeqIO.parse(db_path, 'fasta')]
        df_pepstats = pd.DataFrame({'Protein_ID': protein_id,
                                    'ProteinLength': protein_length,
                                    'IsoelectricPoint': pi,
                                    'MolecularWeight': pm})
        df_pepstats.to_csv(f'{dirs[5]}/pepstats.csv', index=False)
        return df_pepstats

    def get_deeploc(self, dirs, db_path):
        print('Calculating subcellular localization')
        os.system(f'deeploc2 -f {db_path} -o {dirs[6]} > {dirs[11]}/deeploc.log 2>&1')
        file_deeploc = f'{dirs[6]}/{os.listdir(dirs[6])[0]}'
        df_deeploc = pd.read_csv(file_deeploc, usecols=['Protein_ID', 'Localizations'])
        return df_deeploc

    def get_metadados(self, dirs, df_pepstats, df_deeploc, df_pfam_selected_united):
        print('Collecting metadata')
        df_merged = pd.merge(df_pepstats, df_deeploc, on='Protein_ID')
        # ([^_]+)
        inregex = re.compile("^...([A-Za-z0-9._]+)")
        df_merged['Query'] = df_merged['Protein_ID'].apply(lambda x: inregex.search(x).group(1))
        df_merged = df_merged[['Query', 'ProteinLength', 'IsoelectricPoint', 'MolecularWeight', 'Localizations']]
        df_merged_stats = pd.merge(df_merged, df_pfam_selected_united, on='Query')
        df_merged_stats = df_merged_stats[['Specie',
                                           'Assembly',
                                           'Query',
                                           'ProteinLength',
                                           'IsoelectricPoint',
                                           'MolecularWeight',
                                           'Localizations',
                                           'Domains']]
        df_merged_stats.to_csv(f'{dirs[1]}/metadados.csv', index=False)

    def get_filogeny(self, dirs, db_path):
        print('Computing phylogenetic tree')
        cpu = self.parameters.param_cpu
        in_mafft = db_path
        out_mafft = f'{dirs[8]}/out.fasta'
        out_cialign = f'{dirs[9]}/out'
        out_iqtree = f'{dirs[10]}/out_cleaned.fasta'

        # Inicia o processo de alinhamento das sequências usando o MAFFT e CIAlign
        os.system(f'mafft --thread {cpu} --quiet {in_mafft} > {out_mafft}')

        os.system(f'CIAlign --silent \
                            --infile {out_mafft} \
                            --outfile_stem {out_cialign} \
                            --remove_insertions \
                            --crop_ends \
                            --unalign_output \
                            --plot_input \
                            --plot_output ')

        # Copia os dados de sequências alinhadas limpas (cleaned) para a pasta de saída do IQ-TREE
        os.system(f'cp {out_cialign}_cleaned.fasta {dirs[10]}')

        # Usa o IQ-TREE para construir uma árvore filogenética a partir das sequências limpas (cleaned)
        os.system(f'iqtree2 -s {out_iqtree} -nt {cpu} -quiet -B 1000 -alrt 1000')

    def get_local_anot(self, dirs, df_deeploc):
        print('Generating sublocalization anotation')
        ref_colors = ['#ff0000', '#0000ff', '#00ff00', '#ffff00', '#00ffff', '#ff00ff', '#000000', '#ffffff', '#808080']
        labels = list(df_deeploc['Localizations'].unique())
        dict_colors = {label: ref_colors[i % len(ref_colors)] for i, label in enumerate(labels)}
        colors = [values for keys, values in dict_colors.items()]
        shapes = ['1' for i in labels]
        shape_scales = ['1' for i in labels]

        anot_data = {}
        for _, row in df_deeploc.iterrows():
            query = row['Protein_ID']
            label = row['Localizations']
            labels_anot = f"1,1,{dict_colors[label]},1,1,{label}"
            anot_data[query] = labels_anot

        anot = textwrap.dedent(f"""\
        DATASET_SYMBOL
        SEPARATOR COMMA
        DATASET_LABEL,Localization
        COLOR,#ffff00
        GRADIENT_FILL,1
        LEGEND_TITLE,Subcellular location
        LEGEND_SHAPES,{','.join(shapes)}
        LEGEND_COLORS,{','.join(colors)}
        LEGEND_LABELS,{','.join(labels)}
        LEGEND_SHAPE_SCALES,{','.join(shape_scales)}
        DATA
        """)

        # {Protein_ID, symbol, size, color, fill, position, label}
        anot += '\n'.join([f'{key},{values}' for key, values in anot_data.items()])

        log = open(f'{dirs[4]}/itol_localizations.txt', 'w')
        log.write(anot)
        log.close()

    def get_domain_anot(self, dirs, db_path):
        print('Generating domains anotations')

        database = SeqIO.to_dict(SeqIO.parse(db_path, 'fasta'))
        df_selected_itol = pd.read_csv(f'{dirs[2]}/pfam.itol')

        # Signalp6
        table_signalp = pd.read_table(f'{dirs[12]}/output.gff3', header=None, skiprows=1, engine='python')
        table_signalp = table_signalp[[0, 2, 3, 4]]
        table_signalp.columns = ['query', 'domain', 'from', 'to']
        table_signalp['symbol'] = 'RE'
        table_signalp['color'] = '#ff0000'

        # Deep TMHMM
        keywords = ['TMhelix', 'Beta sheet']
        lines = open(f'{dirs[13]}/TMRs.gff3', "r").readlines()
        lines_filtered = [line.strip().split('\t') for line in lines if any(word in line for word in keywords)]
        table_deeptmhmm = pd.DataFrame(lines_filtered, columns=['query', 'domain', 'from', 'to'])
        table_deeptmhmm['symbol'] = 'RE'
        table_deeptmhmm['color'] = '#0000ff'

        merged_df = pd.concat([table_signalp, table_deeptmhmm])

        dicio_itol = {
            'RE': ['#00ff00', '#ffff00', '#00ffff', '#ff00ff', '#000000', '#ffffff', '#808080'],
            'HH': ['#ff0000', '#0000ff', '#00ff00', '#ffff00', '#00ffff', '#ff00ff', '#000000', '#808080'],
            'HV': ['#ff0000', '#0000ff', '#00ff00', '#ffff00', '#00ffff', '#ff00ff', '#000000', '#808080'],
            'EL': ['#ff0000', '#0000ff', '#00ff00', '#ffff00', '#00ffff', '#ff00ff', '#000000', '#808080'],
            'DI': ['#ff0000', '#0000ff', '#00ff00', '#ffff00', '#00ffff', '#ff00ff', '#000000', '#808080'],
            'TR': ['#ff0000', '#0000ff', '#00ff00', '#ffff00', '#00ffff', '#ff00ff', '#000000', '#808080'],
            'TL': ['#ff0000', '#0000ff', '#00ff00', '#ffff00', '#00ffff', '#ff00ff', '#000000', '#808080'],
            'PL': ['#ff0000', '#0000ff', '#00ff00', '#ffff00', '#00ffff', '#ff00ff', '#000000', '#808080'],
            'PR': ['#ff0000', '#0000ff', '#00ff00', '#ffff00', '#00ffff', '#ff00ff', '#000000', '#808080'],
            'PU': ['#ff0000', '#0000ff', '#00ff00', '#ffff00', '#00ffff', '#ff00ff', '#000000', '#808080'],
            'PD': ['#ff0000', '#0000ff', '#00ff00', '#ffff00', '#00ffff', '#ff00ff', '#000000', '#808080'],
            'OC': ['#ff0000', '#0000ff', '#00ff00', '#ffff00', '#00ffff', '#ff00ff', '#000000', '#808080'],
            'GP': ['#ff0000', '#0000ff', '#00ff00', '#ffff00', '#00ffff', '#ff00ff', '#000000', '#808080']}

        list_leg = [[k, value] for k, values in dicio_itol.items() for value in values]
        list_dom = df_selected_itol['domain'].unique().tolist()
        list_num = random.sample(range(len(list_leg)), len(list_dom))
        list_sim = [list_leg[i][0] for i in list_num]
        list_col = [list_leg[i][1] for i in list_num]
        dici_dom = {dom: [sim, col] for dom, sim, col in zip(list_dom, list_sim, list_col)}

        df_selected_itol['symbol'] = df_selected_itol['domain'].apply(
            lambda x: next(v[0] for k, v in dici_dom.items() if x == k))
        df_selected_itol['color'] = df_selected_itol['domain'].apply(
            lambda x: next(v[1] for k, v in dici_dom.items() if x == k))

        # dfs (itol, signalp, deeptmhmm)
        df_selected_itol = df_selected_itol[['query', 'from', 'to', 'domain', 'symbol', 'color']]
        df_domains_architecture = pd.concat([df_selected_itol, merged_df])
        df_domains_architecture['qlen'] = df_domains_architecture['query'].apply(
            lambda x: len(database[x].seq) if x in database else None)

        list_dom.extend(['Signal Peptide', 'Transmembrane Domains'])
        list_sim.extend(['RE', 'RE'])
        list_col.extend(['#ff0000', '#0000ff'])

        anot_data = {}
        for _, row in df_domains_architecture.iterrows():
            query_anot = f"{row['query']},{row['qlen']}"
            domain_anot = f"{row['symbol']}|{row['from']}|{row['to']}|{row['color']}|{row['domain']}"

            if query_anot in anot_data:
                anot_data[query_anot].append(domain_anot)
            else:
                anot_data[query_anot] = [domain_anot]

        anot = textwrap.dedent(f"""\
        DATASET_DOMAINS
        SEPARATOR COMMA
        DATASET_LABEL,Domains
        COLOR,#0000ff
        BORDER_WIDTH,1
        GRADIENT_FILL,1
        SHOW_DOMAIN_LABELS,1
        LEGEND_TITLE,Domain architecture
        LEGEND_LABELS,{','.join(list_dom)}
        LEGEND_SHAPES,{','.join(list_sim)}
        LEGEND_COLORS,{','.join(list_col)}
        DATA
        """)

        anot += '\n'.join([f'{key},{",".join(values)}' for key, values in anot_data.items()])

        log = open(f'{dirs[4]}/itol_domains.txt', 'w')
        log.write(anot)
        log.close()
