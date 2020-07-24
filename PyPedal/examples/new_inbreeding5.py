#!/usr/bin/python

###############################################################################
# NAME: new_inbreeding.py
# VERSION: 2.0.0b10 (10MAY2006)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

# Testind changed to sparse matrix routines on 04/23/2020.

from PyPedal import  pyp_newclasses
from PyPedal import  pyp_nrm
from PyPedal.pyp_utils import pyp_nice_time

options = {}
options['messages'] = 'verbose'
options['renumber'] = 0
options['pedigree_is_renumbered'] = 1
options['pedfile'] = 'mrode.ped'
options['pedformat'] = 'asd'
options['pedname'] = 'My Pedigree'
options['matrix_type'] = 'sparse'
options['debug_messages'] = 0
options['missing_parent'] = 0

if __name__ == '__main__':

    example = pyp_newclasses.loadPedigree(options)
    print example

    for nrm_method in ['nrm', 'frm']:
	print '='*120
        print '[DEBUG]: nrm_method = ', nrm_method
        for matrix_type in ['sparse', 'dense']:
	    print '-'*120
            print '[DEBUG]: matrix_type = ', matrix_type
	    example.kw['matrix_type'] = matrix_type
	    example.kw['nrm_method'] = nrm_method
            print 'Started computing inbreeding at %s' % (pyp_nice_time())
            example_inbreeding_vr = pyp_nrm.inbreeding(example, method='vanraden')
            print 'Finished computing inbreeding at %s' % (pyp_nice_time())
            print '\nVanRaden                  : ', example_inbreeding_vr['metadata']
