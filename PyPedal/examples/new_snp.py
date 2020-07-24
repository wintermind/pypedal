#!/usr/bin/python

###############################################################################
# NAME: new_graphics3.py
# VERSION: 2.0.0b15 (18SEPTEMBER2006)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_newclasses
from PyPedal import pyp_snp

options = {
	'pedname': 'Mrode, 2nd Edition, Example 11.1',
	'pedformat': 'asdP',
	'pedfile': 'mrode_genotypes.ped',
	'snpfile': 'mrode_genotypes.txt',
	'sepchar': "\t",
	'snp_sepchar': "\t",
	'messages': 'verbose',
	'debug': True,
}

if __name__ == '__main__':

	example = pyp_newclasses.loadPedigree(options=options)

	#print example.snp

	G = pyp_snp.form_grm_from_snp(example, scale_m=True, method=1, debug=True)

	g_inbr, g_rels = pyp_snp.compute_genomic_inbreeding_from_grm(example, g_matrix=G, rels=False, update_pedigree=True)

	print g_inbr['metadata']
