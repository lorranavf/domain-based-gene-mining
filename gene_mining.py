import biolib
import contextlib
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
from tqdm import trange


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
        
        print(f'Obtaining domains through hmmscan')
        
        for i in trange(len(files), ncols = 100):
            assembly = regex.search(files[i])[0]
            seq = f'{self.parameters.param_seq_path}/{files[i]}'
            self.get_domains(assembly, seq)
            self.get_resolved_hits(assembly)

    def get_domains(self, assembly, seq):

        infile = self.parameters.param_pfam_in
        outfile = f'{self.parameters.param_pfam_out}/{assembly}'
        cpu = self.parameters.param_cpu
        os.system(f'hmmscan --cpu {cpu} -E 1e-5 -o {outfile}.out --domtblout {outfile}.pfam {infile} {seq}')

    def get_resolved_hits(self, assembly):
        file = f'{self.parameters.param_pfam_out}/{assembly}'
        comm = f'cath-resolve-hits {file}.out --input-format hmmscan_out --hits-text-to-file {file}.resolved  --input-hits-are-grouped --quiet > /dev/null 2>&1'
        os.system(comm)

    def filter_sequences_per_domain(self, dirs):
        regex = re.compile(self.parameters.param_seq_regex)
        ext = self.parameters.param_seq_ext
        files = os.listdir(self.parameters.param_seq_path)
        files = [filename for filename in files if os.path.splitext(filename)[1] == ext]
        tables_pfam = []
        tables_itol = []
        dbs_faa = []

        for i in range(len(files)):
            domain = self.parameters.param_domain
            domain_group = self.parameters.param_domain_group
            dicio = self.parameters.param_seq_dicio
            assembly = regex.search(files[i])[0]
            seq = f'{self.parameters.param_seq_path}/{files[i]}'
            name_ab = next((v[0] for k, v in self.parameters.param_seq_dicio.items() if k == assembly), None)

            # pfamm hmmscan output
            path = f'{self.parameters.param_pfam_out}/{assembly}.resolved'
            columns = ['Query', 'Domain', 'Score', 'Boundaries', 'Resolved', 'Cond-Evalue', 'Indp-Evalue']
            table = pd.read_table(path, sep=' ', header=None, skiprows=2, names=columns, engine='python')

            table_profile = table.groupby('Query')['Domain'].apply(list).reset_index(name='Domains')
            table_profile['Assembly'] = assembly
            table_profile['Specie'] = next((v[1] for k, v in dicio.items() if k == assembly), None)
            table_profile = table_profile[['Specie', 'Assembly', 'Query', 'Domains']]

            print(f'Selecting domain in: {assembly}, {time.strftime("%H:%M %d/%m/%Y", time.localtime(time.time()))}')

            if domain_group:
                table_pfam = table_profile[table_profile['Domains'] == str(domain)]
            else:
                table_pfam = table_profile[table_profile['Domains'].apply(lambda x: any(d in x for d in domain))]

            tables_pfam.append(table_pfam)

            seq_to_select = list(table_pfam['Query'])

            table['From'] = table['Resolved'].str.split('-').str[0]
            table['To'] = table['Resolved'].str.split('-').str[1]
            table_itol = table[['Query', 'Domain', 'From', 'To']]
            table_itol = table_itol[table_itol['Query'].isin(seq_to_select)]
            table_itol['Query'] = table_itol['Query'].apply(lambda x: f'{name_ab}_{x}')
            tables_itol.append(table_itol)

            # faa database
            database = [record for record in SeqIO.parse(seq, 'fasta')]
            database_selected = [record for record in database if record.id in seq_to_select]
            for seq in database_selected:
                seq.id = f'{name_ab}_{seq.id}'
                seq.name = ''
                seq.description = ''
            dbs_faa.append(database_selected)

        table_selected_united = pd.concat(tables_pfam)
        table_selected_united.to_csv(f'{dirs[2]}/pfam_profile.csv', index=False)

        tables_itol_united = pd.concat(tables_itol)
        tables_itol_united.to_csv(f'{dirs[2]}/pfam_coordinates.itol', index=False)

        database_selected_united = [record for database in dbs_faa for record in database]
        SeqIO.write(database_selected_united, f'{dirs[3]}/database_selected_united.fasta', 'fasta')

        db_path = f'{dirs[3]}/database_selected_united.fasta'

        df_pepstats = self.get_pepstats(dirs, db_path)
        df_deeploc = self.get_deeploc(dirs, db_path)
        self.get_signalp(dirs, db_path)
        self.get_deep_tmhmm(dirs, db_path)
        self.get_metadados(dirs, df_pepstats, df_deeploc, table_selected_united)
        self.get_filogeny(dirs, db_path)
        self.get_anotations(dirs, db_path)
        self.get_local_anot(dirs, df_deeploc)

    def get_signalp(self, dirs, db_path):
        print('Calculating signal peptides')
        os.system(f'signalp6 --fastafile {db_path} --organism eukarya --output_dir {dirs[12]} --format txt --mode fast > {dirs[11]}/signalp.log 2>&1')

    def get_deep_tmhmm(self, dirs, db_path):
        print('Calculating transmembrane domains')
        deeptmhmm = biolib.load('DTU/DeepTMHMM')

        deeptmhmm_log = f'{dirs[11]}/deeptmhmm.log'

        with open(deeptmhmm_log, 'w') as log:
            with contextlib.redirect_stdout(log):
                deeptmhmm_job = deeptmhmm.cli(args=f'--fasta {db_path}', machine='local')
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
        inregex = re.compile("_(.*)")
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

    def get_file_anot(self, df, outfile, label, title, list_dom, list_sim, list_col):
        anot = textwrap.dedent(f"""\
                DATASET_DOMAINS
                SEPARATOR COMMA
                DATASET_LABEL,{label}
                COLOR,#0000ff
                BORDER_WIDTH,1
                GRADIENT_FILL,1
                SHOW_DOMAIN_LABELS,1
                LEGEND_TITLE,{title}
                LEGEND_LABELS,{','.join(list_dom)}
                LEGEND_SHAPES,{','.join(list_sim)}
                LEGEND_COLORS,{','.join(list_col)}
                DATA
                """)

        anot_data = {}

        for _, row in df.iterrows():
            query_anot = f"{row['Query']},{row['Qlen']}"
            domain_anot = f"{row['Symbol']}|{row['From']}|{row['To']}|{row['Color']}|{row['Domain']}"

            if query_anot in anot_data:
                anot_data[query_anot].append(domain_anot)
            else:
                anot_data[query_anot] = [domain_anot]

        anot += '\n'.join([f'{key},{",".join(values)}' for key, values in anot_data.items()])

        log = open(outfile, 'w')
        log.write(anot)
        log.close()

    def get_anotations(self, dirs, db_path):
        print('Generating domains anotations')

        # Fasta
        database = SeqIO.to_dict(SeqIO.parse(db_path, 'fasta'))

        # Pfam
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

        df_selected_itol = pd.read_csv(f'{dirs[2]}/pfam_coordinates.itol')
        list_leg = [[k, value] for k, values in dicio_itol.items() for value in values]
        list_dom = df_selected_itol['Domain'].unique().tolist()
        list_num = random.sample(range(len(list_leg)), len(list_dom))
        list_sim = [list_leg[i][0] for i in list_num]
        list_col = [list_leg[i][1] for i in list_num]
        dici_dom = {dom: [sim, col] for dom, sim, col in zip(list_dom, list_sim, list_col)}
        df_selected_itol['Symbol'] = df_selected_itol['Domain'].apply(lambda x: next(v[0] for k, v in dici_dom.items() if x == k))
        df_selected_itol['Color'] = df_selected_itol['Domain'].apply(lambda x: next(v[1] for k, v in dici_dom.items() if x == k))
        df_selected_itol['Qlen'] = df_selected_itol['Query'].apply(lambda x: len(database[x].seq) if x in database else None)
        df_selected_itol = df_selected_itol[['Query', 'Qlen', 'From', 'To', 'Domain', 'Symbol', 'Color']]

        table_signalp = pd.DataFrame(columns=['Query', 'Domain', 'From', 'To', 'Symbol', 'Color'])
        table_deeptmhmm = pd.DataFrame(columns=['Query', 'Domain', 'From', 'To', 'Symbol', 'Color'])
        merged_df = pd.DataFrame(columns=['Query', 'Domain', 'From', 'To', 'Symbol', 'Color'])


        try:
            table_signalp_temp = pd.read_table(f'{dirs[12]}/output.gff3', header=None, skiprows=1, engine='python')

            if not table_signalp_temp.empty:
                table_signalp = table_signalp_temp[[0, 2, 3, 4]]
                table_signalp.columns = ['Query', 'Domain', 'From', 'To']
                table_signalp['Symbol'] = 'RE'
                table_signalp['Color'] = '#ff0000'

        except FileNotFoundError:
            print("Arquivo output.gff3 não encontrado.")
        except pd.errors.EmptyDataError:
            print("O arquivo output.gff3 está vazio ou não contém dados.")
        except Exception as e:
            print(f"Erro ao processar output.gff3: {e}")

        # DeepTMHMM
        try:
            keywords = ['TMhelix', 'Beta sheet']
            lines = open(f'{dirs[13]}/TMRs.gff3', "r").readlines()
            lines_filtered = [line.strip().split('\t') for line in lines if any(word in line for word in keywords)]
            table_deeptmhmm_temp = pd.DataFrame(lines_filtered, columns=['Query', 'Domain', 'From', 'To'])

            if not table_deeptmhmm_temp.empty:
                table_deeptmhmm = table_deeptmhmm_temp.copy()
                table_deeptmhmm['Symbol'] = 'RE'
                table_deeptmhmm['Color'] = '#0000ff'

        except FileNotFoundError:
            print("Arquivo TMRs.gff3 não encontrado.")
        except pd.errors.EmptyDataError:
            print("O arquivo TMRs.gff3 está vazio ou não contém dados.")
        except Exception as e:
            print(f"Erro ao processar TMRs.gff3: {e}")

        # Merge Signalp6 + DeepTMHMM
        if not table_signalp.empty:
            merged_df = pd.concat([table_signalp])

        if not table_deeptmhmm.empty:
            merged_df = pd.concat([table_deeptmhmm])

        if not merged_df.empty:
            merged_df['Qlen'] = merged_df['Query'].apply(lambda x: len(database[x].seq) if x in database else None)
            merged_df = merged_df[['Query', 'Qlen', 'From', 'To', 'Domain', 'Symbol', 'Color']]
            df_domains_architecture = pd.concat([df_selected_itol, merged_df])
        else:
            df_domains_architecture = df_selected_itol.copy()

        lista = [[df_selected_itol, 'Pfam', 'Pfam Architecture', f'{dirs[4]}/itol_pfam_domains.txt'],
                 [df_domains_architecture, 'Domains', 'Domains Arquitecture', f'{dirs[4]}/itol_domains.txt'],
                 [merged_df, 'SP & TM', 'Signal Peptides & Transmembrane Arquitecture', f'{dirs[4]}/itol_tm_domains.txt']]


        for df, label, title, outfile in lista:

            if label == 'Pfam':
                pass

            if label == 'Domains':
                list_dom.extend(['Signal Peptide', 'Transmembrane Domains'])
                list_sim.extend(['RE', 'RE'])
                list_col.extend(['#ff0000', '#0000ff'])

            if label == 'SP & TM':
                list_dom = ['Signal Peptide', 'Transmembrane Domains']
                list_sim = ['RE', 'RE']
                list_col = ['#ff0000', '#0000ff']

            self.get_file_anot(df, outfile, label, title, list_dom, list_sim, list_col)

