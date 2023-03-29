import re
import os
import shutil
import time
from tqdm import tqdm
import pandas as pd
from Bio import SeqIO
import warnings
warnings.filterwarnings("ignore")

def pfam(seq_list, seq_path, dicio, pfam_path, regex):
    """
    Parâmetros:
    
    - seq_list (list): Lista de nomes de arquivos de sequências genômicas em formato FASTA.
    - seq_path (str): Caminho para o diretório onde os arquivos de sequências genômicas estão armazenados.
    - dicio (dict): Dicionário com as informações de montagem do genoma para cada espécie.
    - pfam_path (str): Caminho para o diretório onde o arquivo do banco de dados PFAM está armazenado.
    - regex (re.Pattern): Expressão regular que identifica o código de montagem do genoma nos nomes dos arquivos de sequências genômicas.

    Retorno:
    
    - Não há retorno, mas a função salva os arquivos de saída no diretório já especificado.
    
    """
    # Registra o início da análise
    print(f'Início pfam(): {time.strftime("%d/%m/%Y %H:%M:%S %A", time.localtime(time.time()))}')
    
    # Nomeia os diretórios de saída
    dir_name = f'out.pfam'

    # Cria os diretórios de saída
    if os.path.exists(dir_name):
        shutil.rmtree(dir_name)
    os.mkdir(dir_name)

    # Extrai o código de montagem do genoma
    regex = re.compile(regex)
    
    for i in tqdm(range(len(seq_list))):
        
        # extrai o código do genoma        
        assembly = regex.search(seq_list[i])[0]

        # Path para cada arquivo a ser analisado 
        seq = f'{seq_path}/{seq_list[i]}'

        # Análise de domínios conservados
        outfile = f'{dir_name}/{assembly}'
        os.system(f'hmmscan --cpu 8 -o {outfile}.out --domtblout {outfile}.pfam {pfam_path} {seq}')

        # Lê a tabela gerada pelo hmmscan usando pandas
        table = f'{dir_name}/{assembly}.pfam'
        df = pd.read_table(table, header=None, skiprows=3, skipfooter=10, engine='python')

        # Cria um novo dataframe com colunas 'Query' e 'Domain'
        df = pd.DataFrame({'Query': [i.split()[3] for i in df[0]],
                           'Domain': [i.split()[0] for i in df[0]]})

        # Agrupa os domínios por consulta em uma lista, renomeia a coluna 'Domain' para 'Domains'
        df = df.groupby('Query')['Domain'].apply(list).reset_index(name='Domains')

        # Cria coluna Assembly
        df['Assembly'] = assembly

        # Cria a coluna de identificação da espécie
        df['Specie'] = next((k for k,v in dicio.items() if v[0] == assembly), None)

        # Reorganiza o dataframe
        df = df[['Specie','Assembly','Query','Domains']]

        # Salva os resultados
        df.to_csv(f'{dir_name}/{assembly}.csv', index=False)
        
    # Registra o fim da análise
    print(f'Fim pfam(): {time.strftime("%d/%m/%Y %H:%M:%S %A", time.localtime(time.time()))}')
        
def metadados_filogeny(outdir, domain, seq_list, seq_path, dicio,regex):
    """
    Parâmetros:
    
    - outdir: string contendo o nome do diretório de saída
    - domain: string contendo o nome do domínio a ser analisado
    - seq_list: lista de strings contendo o nome dos arquivos de sequência a serem analisados
    - seq_path: string contendo o caminho para o diretório contendo os arquivos de sequência a serem analisados
    - dicio: dicionário contendo informações adicionais sobre as sequências a serem analisadas
    - regex: objeto regex a ser utilizado na extração do código de montagem do genoma

    Retorno:
    
    - Não há retorno, mas a função salva os arquivos de saída no diretório especificado.

    """
    # Registra o início da análise 
    print(f'\n Início metadados_filogeny(): {time.strftime("%d/%m/%Y %H:%M:%S %A", time.localtime(time.time()))}')
    #####################################################################################################   
    # Nomeia os diretórios de saída
    metadados_dirs = [f'out.{directory}' for directory in [outdir,
                                                           f'{outdir}/metadados',
                                                           f'{outdir}/pfam',
                                                           f'{outdir}/fasta']]

    # Cria os diretórios de saída
    for dir_name in metadados_dirs:
        if os.path.exists(dir_name):
            shutil.rmtree(dir_name)
        os.mkdir(dir_name)

    # Extrai o código de montagem do genoma
    regex = re.compile(regex)
    
    # Cria database para análise filogenética
    dfs = []
    database_filogeny = []
    
    for i in tqdm(range(len(seq_list))):
        
        # extrai o código do genoma        
        assembly = regex.search(seq_list[i])[0]
        
        # Path para cada arquivo a ser analisado 
        seq = f'{seq_path}/{seq_list[i]}'
        
        df = pd.read_csv(f'out.pfam/{assembly}.csv')
        
        # Seleciona as sequências que contém ao menos um domínio alvo
        # df_selected = df.loc[df['Domains'].apply(lambda x: domain in x)]
        df_selected = df[df['Domains'].apply(lambda x: any(d in x for d in domain))]

        # Gera uma lista das sequências selecionadas
        seq_to_select = list(df_selected['Query'])

        # Obtém o nome abreviado da espécie
        name_ab = next((v[1] for v in dicio.values() if v[0] == assembly), None)

        # Lê as sequências de um arquivo FASTA usando BioPython
        database = [record for record in SeqIO.parse(seq, 'fasta')]

        # Filtra as sequências que pertencem às espécies identificados no passo anterior
        database_selected = [record for record in database if record.id in seq_to_select]

        # Renomeia as sequências adicionando a abreviação da espécie
        for seq in database_selected:
            seq.id  = f'{name_ab}_{seq.id}'
            seq.name = ''
            seq.description = ''
            
        # Adiciona database_selected à database_filogeny
        database_filogeny.append(database_selected)
        
        # Adiciona df_selected à lista dfs
        dfs.append(df_selected)

        # Salva os resultados
        # df_selected.to_csv(f'{metadados_dirs[1]}/metadados_seq_selecionadas_{assembly}.csv', index=False)
        # SeqIO.write(database_selected, f'{metadados_dirs[1]}/database_seq_selecionadas_{assembly}.fasta', 'fasta')
    
    #####################################################################################################
    # Une os dataframes em um único dataframe
    df_selected_united = pd.concat(dfs)
    df_selected_united.to_csv(f'{metadados_dirs[2]}/pfam.csv', index=False)
    
    
    # Une e salva as sequências selecionadas e salva em um único arquivo fasta
    database_selected_united = [record for database in database_filogeny for record in database]
    SeqIO.write(database_selected_united, f'{metadados_dirs[3]}/database_selected_united.fasta', 'fasta')
    
    # Defini caminho para a base de dados das próximas análises
    db_path = f'{metadados_dirs[3]}/database_selected_united.fasta'
    
    ####################################################################################################    
    protein_stats_dirs = [f'out.{directory}' for directory in [f'{outdir}/pepstats',
                                                               f'{outdir}/deeploc']]
    
    # Cria os diretórios de saída
    for dir_name in protein_stats_dirs:
        if os.path.exists(dir_name):
            shutil.rmtree(dir_name)
        os.mkdir(dir_name)
    
    out_stats = f'{protein_stats_dirs[0]}/metadados.pepstats'
    
    # Obtém as propriedades de uma sequência de proteína presente em um arquivo fasta, utilizando o programa pepstats
    os.system(f'pepstats -sequence {db_path} -outfile {out_stats}')    
    

    # Lê output do pepstats
    lines = open(out_stats, "r").readlines()
    
    # extraindo valores para colunas
    pi = [float(line.replace("Isoelectric Point = ", "")) for line in lines if line.startswith("Isoelectric Point")]
    pm = [float(line.replace("Molecular weight = ", "").split()[0]) for line in lines if line.startswith("Molecular weight")]
    protein_id = [record.id for record in SeqIO.parse(db_path, 'fasta')]
    protein_length = [len(record.seq) for record in SeqIO.parse(db_path, 'fasta')]

    # Cria o dataframe com os dados gerados pelo pepstats
    df_pepstats = pd.DataFrame({'Protein_ID': protein_id,
                                'ProteinLength': protein_length,
                                'IsoelectricPoint': pi,
                                'MolecularWeight': pm})
    
    # Salva o dataframe     
    df_pepstats.to_csv(f'{protein_stats_dirs[0]}/pepstats.csv', index=False)
    
    #Obtém a localização subcelular das proteínas utilizando o programa deeploc2
    os.system(f'deeploc2 -f {db_path} -o {protein_stats_dirs[1]}')

    # Lê um arquivo CSV contendo informações sobre localizações de proteínas e armazena-o no dataframe df_deeploc.
    file_deeploc = os.path.join(protein_stats_dirs[1], os.listdir(protein_stats_dirs[1])[0])
    df_deeploc = pd.read_csv(file_deeploc, usecols=['Protein_ID', 'Localizations'])

    # Padroniza a coluna "Query" nos dois dataframes e une-os
    df_merged = pd.merge(df_pepstats, df_deeploc, on='Protein_ID')
    
    # Defini a expressão regular que extrai o código das sequências
    inregex = re.compile("^...([A-Za-z0-9._]+)")
    
    # Extrai o código das sequências e salva na coluna Query
    df_merged['Query'] = df_merged['Protein_ID'].apply(lambda x: inregex.search(x).group(1))
    
    # Seleciona as colunas de interesse do dataframe     
    df_merged = df_merged[['Query','ProteinLength','IsoelectricPoint','MolecularWeight','Localizations']]
    
    # Une o dataframe merged (pepstats + deeploc) com o dataframe de metadados do pfam
    df_merged_stats = pd.merge(df_merged, df_selected_united, on='Query')
    
    # Seleciona as colunas de interesse    
    df_merged_stats = df_merged_stats[['Specie',
                                       'Assembly',
                                       'Query',
                                       'ProteinLength',
                                       'IsoelectricPoint',
                                       'MolecularWeight',
                                       'Localizations', 
                                       'Domains']]
    
    # Salva o dataframe final    
    df_merged_stats.to_csv(f'{metadados_dirs[1]}/metadados.csv', index=False)
    
    ####################################################################################################  
    # Nomeia os diretórios de saída da análise filogenética
    filogeny_dirs = [f'out.{directory}' for directory in [f'{outdir}/filogeny',
                                                          f'{outdir}/filogeny/mafft',
                                                          f'{outdir}/filogeny/cialign',
                                                          f'{outdir}/filogeny/iqtree']]

    # Cria os diretórios de saída
    for dir_name in filogeny_dirs:
        if os.path.exists(dir_name):
            shutil.rmtree(dir_name)
        os.mkdir(dir_name)
        
    # Define as variáveis que serão utilizadas ao longo da análise filogenética
    in_mafft = db_path
    out_mafft = f'{filogeny_dirs[1]}/out.fasta'
    out_cialign = f'{filogeny_dirs[2]}/out'
    out_iqtree = f'{filogeny_dirs[3]}/out_cleaned.fasta'

    # Inicia o processo de alinhamento das sequências usando o MAFFT e CIAlign
    os.system(f'mafft --thread 8 --quiet {in_mafft} > {out_mafft}')
    
  
    os.system(f'CIAlign --silent \
                        --infile {out_mafft} \
                        --outfile_stem {out_cialign} \
                        --remove_insertions \
                        --crop_ends \
                        --unalign_output \
                        --plot_input \
                        --plot_output')

    # Copia os dados de sequências alinhadas limpas (cleaned) para a pasta de saída do IQ-TREE
    os.system(f'cp {out_cialign}_cleaned.fasta {filogeny_dirs[3]}')

    # Usa o IQ-TREE para construir uma árvore filogenética a partir das sequências limpas (cleaned)
    os.system(f'iqtree2 -s {out_iqtree} -nt 8 -quiet -B 1000 -alrt 1000')
    
    # Registra o fim da análise
    print(f'Fim metadados_filogeny(): {time.strftime("%d/%m/%Y %H:%M:%S %A", time.localtime(time.time()))}')
    ###############################################################################################       