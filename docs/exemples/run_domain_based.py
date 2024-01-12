import sys
sys.path.append("../..")

from mining import Parameters, DomainAnalysis


codes = {'REF_000000000.0': ['Ata', 'Arabdopsis thaliana', 'Ata'],
         'REF_000000000.1': ['Bra', 'Brassica rapa', 'Bra'],
         'REF_000000000.2': ['Csi', 'Citrus sinensis', 'Csi']}

metadata = Parameters(
            param_seq_path='infiles.domain_based',
            param_seq_dicio=codes,
            param_seq_ext='.faa',
            # param_pfam_in='/home/lcbm/databases/in.pfam/Pfam-A.hmm',
            param_pfam_in='/media/lorrana/LorranaHD/Projetos/DB/in.pfam/Pfam-A.hmm',
            param_domain=['Peptidase_S8'],
            param_domain_group=False,
            param_outdir='peptidase.v1',
            param_hmm_analysis=False,
            param_full_analysis=True,
            param_filogeny_analysis=True,
            param_orthogroup_analysis=True)

DomainAnalysis(metadata).run()
