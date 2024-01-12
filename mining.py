import biolib
import contextlib
import colorsys
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
from ete3 import Tree, faces, TreeStyle, SeqMotifFace

warnings.filterwarnings("ignore")


@dataclass
class Parameters:
    """
    """
    param_seq_ext: str = '.faa'
    param_seq_regex: Pattern[str] = '^\w*.\d'
    param_pfam_in: str = 'in.pfam/Pfam-A.hmm'
    param_blastdb: str = 'in.files.blastp.db/blastdb.faa'
    param_blast_reference: str = 'in.files.blastp.reference/reference.faa'
    param_seq_path: str = 'in.files.db'
    param_pfam_out: str = 'out.pfam'
    param_blast_out: str = 'out.blastp'
    param_outdir: str = None
    param_domain: List[str] = None
    param_seq_dicio: dict = None
    param_domain_group: bool = False
    param_hmm_analysis: bool = True
    param_full_analysis: bool = True
    param_filogeny_analysis: bool = True
    param_blast_analysis: bool = False
    param_orthogroup_analysis: bool = False
    param_cpu: int = 8

class DomainAnalysis:
    """
    """

    def __init__(self, parameters: Parameters):
        self.parameters = parameters

    def run(self):

        if self.parameters.param_blast_analysis:
            self.parameters.param_seq_dicio = {'out_blastp.6': ['', 'Multiple species']}
            self.create_output_directory(analysis='blast_analysis')
            self.get_similar_genes()

        if self.parameters.param_hmm_analysis:
            self.create_output_directory(analysis='hmm_analysis')
            self.process_sequences()

        if self.parameters.param_full_analysis:
            dirs = self.create_output_directory(analysis='full_analysis')
            self.filter_sequences_per_domain(dirs)

    def create_output_directory(self, analysis):

        if analysis == 'blast_analysis':

            dirs = [self.parameters.param_blast_out,
                    self.parameters.param_seq_path]

            for directory in dirs:
                if os.path.exists(directory):
                    shutil.rmtree(directory)
                os.mkdir(directory)

        if analysis == 'hmm_analysis':
            if os.path.exists(self.parameters.param_pfam_out):
                shutil.rmtree(self.parameters.param_pfam_out)
            os.mkdir(self.parameters.param_pfam_out)

        if analysis == 'full_analysis':
            outdir = f'out.{self.parameters.param_outdir}'
            dirs = [outdir,
                    f'{outdir}/metadata',
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
                    f'{outdir}/deep_tmhmm',
                    f'{outdir}/orthofinder',
                    f'{outdir}/orthofinder/input',
                    f'{outdir}/figures']

            for directory in dirs:
                if os.path.exists(directory):
                    shutil.rmtree(directory)
                os.mkdir(directory)

            return dirs

    def get_similar_genes(self):

        now = time.strftime("%H:%M %d/%m/%Y", time.localtime(time.time()))
        print(f'Obtaining proteins similar to the reference, {now}')

        blastdb = self.parameters.param_blastdb
        reference = self.parameters.param_blast_reference
        outfmt6 = f'{self.parameters.param_blast_out}/blastp.out'
        outfasta = f'{self.parameters.param_seq_path}/out_blastp.6.faa'

        comm_make_blastdb = f'makeblastdb -in {blastdb} -dbtype prot -parse_seqids > /dev/null 2>&1'
        os.system(comm_make_blastdb)

        comm_blastp = f'blastp -query {reference} -db {blastdb} -out {outfmt6} -evalue 1e-5 -outfmt 6'
        os.system(comm_blastp)

        table = pd.read_table(outfmt6, header=None)
        target_list = table[1].unique().tolist()

        database = [record for record in SeqIO.parse(blastdb, 'fasta')]
        database_selected = [record for record in database if record.id in target_list]
        database_reference = [record for record in SeqIO.parse(reference, 'fasta')]
        databases = [database_reference, database_selected]
        database_selected_united = [record for database in databases for record in database]

        for seq in database_selected_united:
            seq.name = ''
            seq.description = ''
        SeqIO.write(database_selected_united, outfasta, 'fasta')

    def process_sequences(self):
        regex = re.compile(self.parameters.param_seq_regex)
        ext = self.parameters.param_seq_ext
        files = os.listdir(self.parameters.param_seq_path)
        files = [filename for filename in files if os.path.splitext(filename)[1] == ext]
        
        print(f'Obtaining domains through hmmscan, {time.strftime("%H:%M %d/%m/%Y", time.localtime(time.time()))}')
        
        for i in trange(len(files), ncols=100):
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
        comm = f'cath-resolve-hits {file}.out \
                --input-format hmmscan_out \
                --hits-text-to-file {file}.resolved  \
                --input-hits-are-grouped \
                --quiet > /dev/null 2>&1'
        os.system(comm)

    def filter_sequences_per_domain(self, dirs):
        
        now = time.strftime("%H:%M %d/%m/%Y", time.localtime(time.time()))
        print(f'Choosing sequences containing the target domain, {now}')

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

            path = f'{self.parameters.param_pfam_out}/{assembly}.resolved'
            columns = ['Query', 'Domain', 'Score', 'Boundaries', 'Resolved', 'Cond-Evalue', 'Indp-Evalue']
            table = pd.read_table(path, sep=' ', header=None, skiprows=2, names=columns, engine='python')

            table_profile = table.groupby('Query')['Domain'].apply(list).reset_index(name='Domains')
            table_profile['Assembly'] = assembly
            table_profile['Specie'] = next((v[1] for k, v in dicio.items() if k == assembly), None)
            table_profile = table_profile[['Specie', 'Assembly', 'Query', 'Domains']]

            if domain_group:
                table_pfam = table_profile[table_profile['Domains'].apply(lambda x: x == domain)]

            else:
                table_pfam = table_profile[table_profile['Domains'].apply(lambda x: any(d in x for d in domain))]

            seq_to_select = list(table_pfam['Query'])

            table['From'] = table['Resolved'].str.split('-').str[0]
            table['To'] = table['Resolved'].str.split('-').str[1]
            table_itol = table[['Query', 'Domain', 'From', 'To']]
            table_itol = table_itol[table_itol['Query'].isin(seq_to_select)]
            table_itol['Query'] = table_itol['Query'].apply(lambda x: f'{name_ab}_{x}')

            database = [record for record in SeqIO.parse(seq, 'fasta')]
            database_selected = [record for record in database if record.id in seq_to_select]
            for seq in database_selected:
                seq.id = f'{name_ab}_{seq.id}'
                seq.name = ''
                seq.description = ''
            
            SeqIO.write(database_selected, f'{dirs[15]}/{assembly}.fasta', 'fasta')

            tables_pfam.append(table_pfam)
            tables_itol.append(table_itol)
            dbs_faa.append(database_selected)

        table_selected_united = pd.concat(tables_pfam)
        if self.parameters.param_blast_analysis:
            table_selected_united = table_selected_united[['Query', 'Domains']]
        table_selected_united.to_csv(f'{dirs[2]}/pfam_profile.csv', index=False)

        tables_itol_united = pd.concat(tables_itol)
        tables_itol_united.to_csv(f'{dirs[2]}/pfam_coordinates.itol', index=False)

        database_selected_united = [record for database in dbs_faa for record in database]
        SeqIO.write(database_selected_united, f'{dirs[3]}/database_selected_united.fasta', 'fasta')

        db_path = f'{dirs[3]}/database_selected_united.fasta'

        df_pepstats = self.get_pepstats(dirs, db_path)
        df_deeploc = self.get_deeploc(dirs, db_path)
        self.get_metadados(dirs, df_pepstats, df_deeploc, table_selected_united)
        self.get_signalp(dirs, db_path)
        self.get_deep_tmhmm(dirs, db_path)
        self.get_annotations(dirs, db_path)
        self.get_local_annot(dirs, df_deeploc)

        if self.parameters.param_filogeny_analysis:
            self.get_filogeny(dirs, db_path)

        if self.parameters.param_orthogroup_analysis:

            self.get_orthogroups(dirs) 

            orthodir = os.listdir(f'{dirs[14]}/output')[0]

            tree_species = os.path.join(f'{dirs[14]}/output', orthodir, 'Species_Tree/SpeciesTree_rooted.txt')   
            
            self.get_species_tree(dirs, tree_species)    

            tree_orthogroups = os.path.join(f'{dirs[14]}/output', orthodir, 'Resolved_Gene_Trees')  

            self.get_orthogroups_trees(dirs, tree_orthogroups)

    @staticmethod
    def get_signalp(dirs, db_path):
        
        now = time.strftime("%H:%M %d/%m/%Y", time.localtime(time.time()))
        print(f'Calculating signal peptides, {now}')

        log = f'{dirs[11]}/signalp.log 2>&1'
        comm = f'signalp6 --fastafile {db_path} \
                          --organism eukarya \
                          --output_dir {dirs[12]} \
                          --format txt \
                          --mode fast > {log}'
        os.system(comm)

    @staticmethod
    def get_deep_tmhmm(dirs, db_path):
        
        now = time.strftime("%H:%M %d/%m/%Y", time.localtime(time.time()))
        print(f'Calculating transmembrane domains, {now}')
        
        deeptmhmm = biolib.load('DTU/DeepTMHMM')
        deeptmhmm_log = f'{dirs[11]}/deeptmhmm.log'

        with open(deeptmhmm_log, 'w') as log:
            with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                deeptmhmm_job = deeptmhmm.cli(args=f'--fasta {db_path}', machine='local')
                deeptmhmm_job.save_files(dirs[13])

    @staticmethod
    def get_pepstats(dirs, db_path):
        
        now = time.strftime("%H:%M %d/%m/%Y", time.localtime(time.time()))
        print(f'Calculating protein stats, {now}')
        
        out_stats = f'{dirs[5]}/pepstats.out'
        log = f'{dirs[11]}/pepstats.log 2>&1'
        comm = f'pepstats -sequence {db_path} -outfile {out_stats} > {log}'
        os.system(comm)

        lines = open(out_stats, "r").readlines()
        pi = [float(line.replace("Isoelectric Point = ", "")) for line in lines if
              line.startswith("Isoelectric Point")]
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

    @staticmethod
    def get_deeploc(dirs, db_path):
        now = time.strftime("%H:%M %d/%m/%Y", time.localtime(time.time()))
        print(f'Calculating subcellular localization, {now}')

        log = f'{dirs[11]}/deeploc.log 2>&1'
        comm = f'deeploc2 -f {db_path} -o {dirs[6]} > {log}'
        os.system(comm)

        file_deeploc = f'{dirs[6]}/{os.listdir(dirs[6])[0]}'
        df_deeploc = pd.read_csv(file_deeploc, usecols=['Protein_ID', 'Localizations'])
        return df_deeploc

    def get_metadados(self, dirs, df_pepstats, df_deeploc, df_pfam_selected_united):
        now = time.strftime("%H:%M %d/%m/%Y", time.localtime(time.time()))
        print(f'Collecting metadata, {now}')

        df_merged = pd.merge(df_pepstats, df_deeploc, on='Protein_ID')
        inregex = re.compile("_(.*)")
        df_merged['Query'] = df_merged['Protein_ID'].apply(lambda x: inregex.search(x).group(1))
        df_merged = df_merged[['Query', 'ProteinLength', 'IsoelectricPoint', 'MolecularWeight', 'Localizations']]
        df_merged_stats = pd.merge(df_merged, df_pfam_selected_united, on='Query')

        if self.parameters.param_blast_analysis:
            columns = ['Query', 'ProteinLength', 'IsoelectricPoint', 'MolecularWeight', 'Localizations', 'Domains']
        else:
            columns = ['Specie',
                       'Assembly',
                       'Query',
                       'ProteinLength',
                       'IsoelectricPoint',
                       'MolecularWeight',
                       'Localizations',
                       'Domains']

        df_merged_stats = df_merged_stats[columns]
        df_merged_stats.to_csv(f'{dirs[1]}/metadata.csv', index=False)

    def get_filogeny(self, dirs, db_path):
        now = time.strftime("%H:%M %d/%m/%Y", time.localtime(time.time()))
        print(f'Computing phylogenetic tree, {now}')

        cpu = self.parameters.param_cpu
        in_mafft = db_path
        out_mafft = f'{dirs[8]}/out.fasta'
        out_cialign = f'{dirs[9]}/out'
        out_iqtree = f'{dirs[10]}/out_cleaned.fasta'

        os.system(f'mafft --thread {cpu} --quiet {in_mafft} > {out_mafft}')
        os.system(f'CIAlign --silent \
                            --infile {out_mafft} \
                            --outfile_stem {out_cialign} \
                            --remove_insertions \
                            --crop_ends \
                            --unalign_output \
                            --plot_input \
                            --plot_output ')
        os.system(f'cp {out_cialign}_cleaned.fasta {dirs[10]}')
        os.system(f'iqtree2 -s {out_iqtree} -nt {cpu} -quiet -B 1000 -alrt 1000')

    def get_orthogroups(self, dirs):
        now = time.strftime("%H:%M %d/%m/%Y", time.localtime(time.time()))
        print(f'Obtaining orthogroups, {now}')

        infile = dirs[15]
        outfile = f'{dirs[14]}/output'
        cpu = self.parameters.param_cpu
        log = f'{dirs[11]}/orthofinder.log 2>&1'
        command = f'orthofinder -f {infile} -o {outfile} -a {cpu} > {log}' 
        os.system(command)
    
    def get_species_tree(self, dirs, filename):

        code2names = self.parameters.param_seq_dicio

        orders = set([value[2] for value in code2names.values()])
        colors = self.get_colors(orders)

        order2color = {order:color for order,color in zip(orders, colors)}

        def get_node_name(node):

            if node.is_leaf():

                name = code2names[node.name][0]
                order = code2names[node.name][2]
                color = order2color[order]

                # node name
                longNameFace = faces.TextFace(name, fsize=15, fgcolor='dark', bold=True)
                faces.add_face_to_node(longNameFace, node, column=0, position="aligned")

                # node description
                descFace = faces.TextFace(order, fsize=15)
                descFace.hz_align = 1
                descFace.background.color = color
                descFace.margin_left = 10
                faces.add_face_to_node(descFace, node, column=1, aligned=True)
                
                # node shape
                node.img_style['size'] = 0

            else:
                # node shape
                node.img_style["size"] = 0

        t = Tree(filename)
        ts = TreeStyle()
        ts.layout_fn = get_node_name
        ts.show_leaf_name = True
        ts.show_branch_support = True
        ts.show_branch_length = True
        
        outdir = f'{dirs[16]}/out.tree.species'
        os.mkdir(outdir)

        outfiles = [f'{outdir}/{os.path.basename(outfile)}{extension}' for extension in ['.pdf','.png','.svg'] for outfile in [os.path.splitext(filename)[0]]]
        
        for outfile in outfiles:
            t.render(outfile, tree_style=ts)

    def get_orthogroups_trees(self, dirs, indir):

        outdir = f'{dirs[16]}/out.tree.orthogroups'

        os.mkdir(outdir)
        
        database = f'{dirs[0]}/fasta/database_selected_united.fasta'
        pfam = f'{dirs[0]}/itol/pfam.txt'
        signals = f'{dirs[0]}/itol/signals.txt'
        domains = f'{dirs[0]}/itol/domains.txt'
        localizations = f'{dirs[0]}/itol/localizations.txt'

        def get_domains(domains):

            domains_list = [
                [
                    int(domain.split(sep='|')[1]),
                    int(domain.split(sep='|')[2]),
                    '[]',
                    None,
                    10,
                    domain.split(sep='|')[3],
                    domain.split(sep='|')[3],
                    f'arial|10|black|{domain.split(sep="|")[4]}'

                ] 
                for domain in domains if len(domain.split(sep='|')) == 5]

            return domains_list

        lines = lambda x: [line.strip() for line in open(x)][12:]

        name2seq = {record.name: str(record.seq) for record in SeqIO.parse(database, 'fasta')}
        name2motifs = {line.split(sep=',')[0]: get_domains(line.split(sep=',')[2:]) for line in lines(pfam)}
        name2signals = {line.split(sep=',')[0]: get_domains(line.split(sep=',')[2:]) for line in lines(signals)}
        name2domains = {line.split(sep=',')[0]: get_domains(line.split(sep=',')[2:]) for line in lines(domains)}
        name2loc = {line.split(',')[0]: [line.split(',')[3], line.split(',')[6]] for line in lines(localizations)}
        
        code2names = self.parameters.param_seq_dicio

        orders = set([value[2] for value in code2names.values()])
        colors = self.get_colors(orders)
        order2color = {order:color for order,color in zip(orders, colors)}

        def get_motif_tree(filename, outdir,
                           name2seq, name2motifs, name2signals, name2domains, name2loc, 
                           code2names, order2color,
                           domain_structure=True, signal_structure=False, pfam_structure=False):

        
            

            def get_node_faces(node):
        
                if node.is_leaf():

                    name = node.name[15:]
                    longNameFace = faces.TextFace(name, fsize=10, fgcolor='dark', bold=True)
                    longNameFace.margin_right = 40
                    faces.add_face_to_node(longNameFace, node, column=0, position="aligned")

                    precode = node.name[0:15]
                    subcode = precode.rsplit('_', 1)
                    code = f'{subcode[0]}.{subcode[1]}'
                    order = code2names[code][2]
                    color = order2color[order]
                    descFace = faces.TextFace(order, fsize=10)
                    descFace.margin_left = 10
                    descFace.margin_right = 10
                    descFace.hz_align = 1
                    descFace.inner_background.color = color
                    faces.add_face_to_node(descFace, node, column=1, aligned=True)

                    if name in name2loc:
                        loc = name2loc.get(name)[1]
                        color = name2loc.get(name)[0]

                    else:
                        loc = ''
                        color = None
                    loc_face = faces.TextFace(loc, fsize=10)
                    loc_face.inner_background.color = color
                    loc_face.hz_align = 1
                    loc_face.margin_left = 10
                    loc_face.margin_right = 10
                    faces.add_face_to_node(loc_face, node, column=2, aligned=True)

                    if domain_structure:
                        sequence = name2seq.get(name)
                        motifs = name2domains.get(name)
                        domain_face = SeqMotifFace(sequence, motifs=motifs, seq_format="-")
                        domain_face.margin_left = 10
                        faces.add_face_to_node(domain_face, node, column=3, aligned=True)

                    if signal_structure:
                        sequence = name2seq.get(name)
                        motifs = name2signals.get(name)
                        signal_face = SeqMotifFace(sequence, motifs=motifs, seq_format="-")
                        signal_face.margin_left = 10
                        faces.add_face_to_node(signal_face, node, column=3, aligned=True)

                    if pfam_structure:
                        sequence = name2seq.get(name)
                        motifs = name2motifs.get(name)
                        pfam_face = SeqMotifFace(sequence, motifs=motifs, seq_format="-")
                        pfam_face.margin_left = 10
                        faces.add_face_to_node(pfam_face, node, column=3, aligned=True)
                    
                    node.img_style['size'] = 0

                    
                else:
                    node.img_style["size"] = 0

            t = Tree(filename, format=1)
            ts = TreeStyle()
            ts.layout_fn = get_node_faces
            # ts.scale =  120
            ts.show_leaf_name = False
            ts.show_branch_support = True
            ts.show_branch_length = True

            outfiles = [f'{outdir}/{os.path.basename(outfile)}{extension}' for extension in ['.pdf','.png','.svg'] for outfile in [os.path.splitext(filename)[0]]]

            for outfile in outfiles:
                t.render(outfile, tree_style=ts)

        trees = [os.path.join(indir,tree) for tree in os.listdir(indir)]

        for tree in trees:
            get_motif_tree(filename, tree,
                           name2seq, name2motifs, name2signals, name2domains, name2loc, 
                           code2names, order2color,
                           domain_structure=True, signal_structure=False, pfam_structure=False)
        
    @staticmethod
    def get_colors(labels):

        def random_color():
            hue = random.random()
            saturation = 0.5
            luminosity = 0.5
            rgb = colorsys.hls_to_rgb(hue, luminosity, saturation)
            color = "#{:02x}{:02x}{:02x}".format(
            int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))
            return color

        colors = set()

        while len(colors) < len(labels):

            color = random_color()

            if color == '#FFFFFF':
                pass
            elif color == '#ff0000':
                pass
            elif color == '#0000ff':
                pass
            else:
                colors.add(color)

        return colors

    def get_local_annot(self, dirs, df_deeploc):
        print(f'Generating sublocalization annotation, {time.strftime("%H:%M %d/%m/%Y", time.localtime(time.time()))}')
        
        labels = list(df_deeploc['Localizations'].unique())
        shapes = ['1' for _ in labels]
        scales = ['1' for _ in labels]
        colors = self.get_colors(labels)

        label2color = {label:color for label,color in zip(labels, colors)}

        annot_data = {}
        for _, row in df_deeploc.iterrows():
            query = row['Protein_ID']
            label = row['Localizations']
            labels_annot = f"1,1,{label2color[label]},1,1,{label}"
            annot_data[query] = labels_annot

        annot = textwrap.dedent(f"""\
        DATASET_SYMBOL
        SEPARATOR COMMA
        DATASET_LABEL,Localization
        COLOR,#ffff00
        GRADIENT_FILL,1
        LEGEND_TITLE,Subcellular location
        LEGEND_SHAPES,{','.join(shapes)}
        LEGEND_COLORS,{','.join(colors)}
        LEGEND_LABELS,{','.join(labels)}
        LEGEND_SHAPE_SCALES,{','.join(scales)}
        DATA
        """)

        annot += '\n'.join([f'{key},{values}' for key, values in annot_data.items()])

        log = open(f'{dirs[4]}/localizations.txt', 'w')
        log.write(annot)
        log.close()

    @staticmethod
    def get_file_annot(df, outfile, label, title, domains, shapes, colors):
        annot = textwrap.dedent(f"""\
                DATASET_DOMAINS
                SEPARATOR COMMA
                DATASET_LABEL,{label}
                COLOR,#0000ff
                BORDER_WIDTH,1
                GRADIENT_FILL,1
                SHOW_DOMAIN_LABELS,1
                LEGEND_TITLE,{title}
                LEGEND_LABELS,{','.join(domains)}
                LEGEND_SHAPES,{','.join(shapes)}
                LEGEND_COLORS,{','.join(colors)}
                DATA
                """)

        annot_data = {}

        for _, row in df.iterrows():
            query_annot = f"{row['Query']},{row['Qlen']}"
            domain_annot = f"{row['Symbol']}|{row['From']}|{row['To']}|{row['Color']}|{row['Domain']}"

            if query_annot in annot_data:
                annot_data[query_annot].append(domain_annot)
            else:
                annot_data[query_annot] = [domain_annot]

        annot += '\n'.join([f'{key},{",".join(values)}' for key, values in annot_data.items()])

        log = open(outfile, 'w')
        log.write(annot)
        log.close()

    def get_annotations(self, dirs, db_path):
        print(f'Generating domain annotations, {time.strftime("%H:%M %d/%m/%Y", time.localtime(time.time()))}')

        database = SeqIO.to_dict(SeqIO.parse(db_path, 'fasta'))
        df_itol = pd.read_csv(f'{dirs[2]}/pfam_coordinates.itol')

        domains = list(df_itol['Domain'].unique())
        shapes = ['RE' for _ in range(len(domains))]
        colors = self.get_colors(domains)

        domain2color = {domain: [shape, color] for domain, shape, color in zip(domains, shapes, colors)}

        df_itol['Symbol'] = df_itol['Domain'].apply(lambda x: next(v[0] for k, v in domain2color.items() if x == k))
        df_itol['Color'] = df_itol['Domain'].apply(lambda x: next(v[1] for k, v in domain2color.items() if x == k))
        df_itol['Qlen'] = df_itol['Query'].apply(lambda x: len(database[x].seq) if x in database else None)
        df_itol = df_itol[['Query', 'Qlen', 'From', 'To', 'Domain', 'Symbol', 'Color']]

        table_signalp = pd.DataFrame(columns=['Query', 'Domain', 'From', 'To', 'Symbol', 'Color'])
        table_deeptmhmm = pd.DataFrame(columns=['Query', 'Domain', 'From', 'To', 'Symbol', 'Color'])
        merged_df = pd.DataFrame(columns=['Query', 'Domain', 'From', 'To', 'Symbol', 'Color'])

        try:
            table_signalp_temp = pd.read_table(f'{dirs[12]}/output.gff3', header=None, skiprows=1, engine='python')

            if not table_signalp_temp.empty:
                table_signalp = table_signalp_temp.copy()
                table_signalp = table_signalp[[0, 2, 3, 4]]
                table_signalp.columns = ['Query', 'Domain', 'From', 'To']
                table_signalp['Symbol'] = 'RE'
                table_signalp['Color'] = '#ff0000'

        except pd.errors.EmptyDataError:
            error_signalp = 'Empty table'

        keywords = ['TMhelix', 'Beta sheet']
        lines = open(f'{dirs[13]}/TMRs.gff3', "r").readlines()
        lines_filtered = [line.strip().split('\t') for line in lines if any(word in line for word in keywords)]
        table_deeptmhmm_temp = pd.DataFrame(lines_filtered, columns=['Query', 'Domain', 'From', 'To'])

        if not table_deeptmhmm_temp.empty:
            table_deeptmhmm = table_deeptmhmm_temp.copy()
            table_deeptmhmm['Symbol'] = 'RE'
            table_deeptmhmm['Color'] = '#0000ff'

        if not table_signalp.empty:
            merged_df = pd.concat([merged_df, table_signalp])

        if not table_deeptmhmm.empty:
            merged_df = pd.concat([merged_df, table_deeptmhmm])

        if not merged_df.empty:
            merged_df['Qlen'] = merged_df['Query'].apply(lambda x: len(database[x].seq) if x in database else None)
            merged_df = merged_df[['Query', 'Qlen', 'From', 'To', 'Domain', 'Symbol', 'Color']]
            df_domains_architecture = pd.concat([df_itol, merged_df])
        else:
            df_domains_architecture = df_itol.copy()

        objects = [[df_itol, 'Pfam', 'Pfam Architecture', f'{dirs[4]}/pfam.txt'],
                 [df_domains_architecture, 'Domains', 'Domains Arquitecture', f'{dirs[4]}/domains.txt'],
                 [merged_df, 'SP & TM', 'Signal Peptides & Transmembrane Arquitecture', f'{dirs[4]}/signals.txt']]

        for df, label, title, outfile in objects:

            if label == 'Pfam':
                pass

            if label == 'Domains':
                domains.extend(['Signal Peptide', 'Transmembrane Domains'])
                shapes.extend(['RE', 'RE'])
                colors.update(['#ff0000', '#0000ff'])

            if label == 'SP & TM':
                domains = ['Signal Peptide', 'Transmembrane Domains']
                shapes = ['RE', 'RE']
                colors = ['#ff0000', '#0000ff']

            self.get_file_annot(df, outfile, label, title, domains, shapes, colors)
