from mining import Parameters, DomainAnalysis


codes = {'REF_000000000.0': ['Ata', 'Arabdopsis thaliana'],
         'REF_000000000.1': ['Bra', 'Brassica rapa'],
         'REF_000000000.2': ['Csi', 'Citrus sinensis'],
         'REF_000000000.3': ['Egr', 'Eucalyptus grandis'],
         'REF_000000000.4': ['Gma', 'Glycine max'],
         'REF_000000000.5': ['Mac', 'Musa acuminata'],
         'REF_000000000.6': ['Mtr', 'Medicago truncatula'],
         'REF_000000000.7': ['Sly', 'Solanum lycopersicum'],
         'REF_000000000.8': ['Vvi', 'Vitis vinifera'],
         'REF_000000000.9': ['Zma', 'Zea mays']}

metadata = Parameters(
            param_seq_path='exemple',
            param_seq_dicio=codes,
            param_seq_ext='.faa',
            param_pfam_in='in.pfam/Pfam-A.hmm',
            param_domain=['Inhibitor_I9', 'Peptidase_S8'],
            param_domain_group=True,
            param_outdir='peptidase.v1',
            param_hmm_analysis=True,
            param_full_analysis=True)

DomainAnalysis(metadata).run()
