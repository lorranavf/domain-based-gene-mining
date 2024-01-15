from mining import DomainAnalysis, Parameters

metadata = Parameters(
    param_outdir='wox',
    param_domain=['Homeodomain'],
    param_blast_analysis=True,
    param_pfam_in='in.pfam/Pfam-A.hmm',
)

DomainAnalysis(metadata).run()
