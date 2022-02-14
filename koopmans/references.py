'''
pseudo_dojo
GNU GPL
D. R. Hamann, Phys. Rev. B 88, 085117 (2013)

sg15
CC BY-SA 4.0
M. Schlipf and F. Gygi, Computer Physics Communications 196, 36 (2015).
http: // dx.doi.org/10.1016/j.cpc.2015.05.011
FR -> P. Scherpelz, M. Govoni, I. Hamada, G. Galli
J. Chem. Theory Comput. (2016)
http: // dx.doi.org/10.1021/acs.jctc.6b00114".

SSSP
Mix of licenses
SSSP: G. Prandini, A. Marrazzo, I. E. Castelli, N. Mounet and N. Marzari, npj Computational Materials 4, 72 (2018).
WEB: http: // materialscloud.org/sssp.
K. Lejaeghere et al., Science 351 (6280), 1415 (2016).
DOI: 10.1126/science.aad3000, WEB: http: // molmod.ugent.be/deltacodesdft
As well as a mix of pseudopotential attributions(see https: // www.materialscloud.org/discover/sssp/table/efficiency for details)
'''

from pathlib import Path
from pybtex.database.input import bibtex

parser = bibtex.Parser()
bib_file = Path(__file__).parents[1] / 'docs/refs.bib'
bib_data = parser.parse_file(bib_file)
