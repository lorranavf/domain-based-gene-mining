import sys

sys.path.append('../..')

from mining import DomainAnalysis, Parameters

codes = {
    'REF_000000000.0': ['Ata', 'Arabdopsis thaliana', 'Ata'],
    'REF_000000000.1': ['Bra', 'Brassica rapa', 'Bra'],
    'REF_000000000.2': ['Csi', 'Citrus sinensis', 'Csi'],
}

metadata = Parameters(
    seq_path='infiles.domain_based',
    seq_dicio=codes,
    seq_ext='.faa',
    # pfam_in='/home/lcbm/databases/in.pfam/Pfam-A.hmm',
    pfam_in='/media/lorrana/LorranaHD/Projetos/DB/in.pfam/Pfam-A.hmm',
    domain=['Peptidase_S8'],
    domain_group=False,
    outdir='peptidase.v1',
    hmm_analysis=True,
    busco_analysis=True,
    full_analysis=True,
    filogeny_analysis=True,
    orthogroup_analysis=True,
)

DomainAnalysis(metadata).run()
