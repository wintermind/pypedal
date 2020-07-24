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
	'pedname': 'Made-up data for testing new SNP code',
	'pedformat': 'asdP',
	'pedfile': 'fake_genotypes.ped',
	'snpfile': 'fake_genotypes.txt',
	'sepchar': "\t",
	'snp_sepchar': "\t",
	'messages': 'verbose',
	'debug': True,
}

if __name__ == '__main__':

	example = pyp_newclasses.loadPedigree(options=options)

	g_homo = pyp_snp.compute_genomic_homozygosity_from_snp(example, debug=True)

	print g_homo['metadata']
