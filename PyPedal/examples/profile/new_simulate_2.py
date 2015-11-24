#!/usr/bin/python

###############################################################################
# NAME: new_simulate.py
# VERSION: 2.0.0b15 (08JUNE2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_newclasses
from PyPedal import pyp_graphics
from PyPedal import pyp_network
from PyPedal import pyp_nrm
from PyPedal import pyp_utils

if __name__ == '__main__':

	print 'Simulating pedigree at %s' % pyp_utils.pyp_nice_time()
        example = pyp_newclasses.loadPedigree(optionsfile='new_simulate.ini')
	#print 'Calculating inbreeding at %s' % pyp_utils.pyp_nice_time()
	#@profile
	#example_inbreeding = pyp_nrm.inbreeding(example, method='vanraden')
	print 'Finished at %s' % pyp_utils.pyp_nice_time()
