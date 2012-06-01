#!/usr/bin/python

###############################################################################
# NAME: new_inbreeding.py
# VERSION: 2.0.0b10 (10MAY2006)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################
# Test the new inbreeding routines and compare their performance to the old
# routine using simulated pedigrees of various sizes.
###############################################################################
from PyPedal import *

options = {}
options['simulate_pedigree'] = 1
options['simulate_mp'] = 0
options['pedigree_save'] = 1
options['pedigree_summary'] = 2
options['messages'] = 'quiet'
options['renumber'] = 1
options['pedformat'] = 'asdxg'

if __name__ == '__main__':

	sizes = [100, 1000, 10000, 100000]
	for size in sizes:

		options['simulate_n'] = size
		options['pedname'] = 'Simulated Pedigree ' + str(size)
		options['pedfile'] = 'simulated_pedigree_' +  str(size) + '.ped'
		print '\nStarted pedigree simulation with %s animals at %s' % (size, pyp_utils.pyp_nice_time())
	        test = pyp_newclasses.loadPedigree(options)
		print 'Finished pedigree simulation with %s animals at %s' % (size, pyp_utils.pyp_nice_time())

		print '\n\tStarted computing inbreeding using VanRaden\'s method  at %s' % (pyp_utils.pyp_nice_time())
		test_inbreeding_vr = pyp_nrm.inbreeding(test, method='vanraden')
		print '\tFinished computing inbreeding using VanRaden\'s method  at %s' % (pyp_utils.pyp_nice_time())

		print '\n\tStarted computing inbreeding using Meuwissen and Luo\'s method  at %s' % (pyp_utils.pyp_nice_time())
		test_inbreeding_ml = pyp_nrm.inbreeding(test, method='meu_luo')
		print '\tFinished computing inbreeding using Meuwissen and Luo\'s method  at %s' % (pyp_utils.pyp_nice_time())

		print '\n\tStarted computing inbreeding using Meuwissen and Luo\'s modified method  at %s' % (pyp_utils.pyp_nice_time())
		test_inbreeding_qu = pyp_nrm.inbreeding(test, method='mod_meu_luo')
		print '\tFinished computing inbreeding using Meuwissen and Luo\'s modified method  at %s' % (pyp_utils.pyp_nice_time())

		print '\n\tVanRaden                  : ', test_inbreeding_vr['metadata']
		print '\n\tMeuwissen and Luo         : ', test_inbreeding_ml['metadata']
		print '\n\tModified Meuwissen and Luo: ', test_inbreeding_qu['metadata']
