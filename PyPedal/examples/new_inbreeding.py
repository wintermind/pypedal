#!/usr/bin/python

###############################################################################
# NAME: new_inbreeding.py
# VERSION: 2.0.0b10 (10MAY2006)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import  pyp_newclasses
from PyPedal import  pyp_nrm
from PyPedal.pyp_utils import pyp_nice_time

options = {}
options['messages'] = 'verbose'
options['renumber'] = 0
#options['pedfile'] = 'new_renumbering.ped'
options['pedfile'] = 'mrode.ped'
#options['pedformat'] = 'asdgb'
options['pedformat'] = 'asd'
options['pedname'] = 'My Pedigree'
options['matrix_type'] = 'sparse'
options['renumber'] = 0
options['pedigree_is_renumbered'] = 1

if __name__ == '__main__':

    # example = pyp_newclasses.loadPedigree(optionsfile='new_inbreeding.ini')
    example = pyp_newclasses.loadPedigree(options)
    print example

    #print '[DEBUG]: matrix_type = ', example.kw['matrix_type']

    print 'Started computing inbreeding at %s' % (pyp_nice_time())
    example_inbreeding_vr = pyp_nrm.inbreeding(example,method='vanraden')
    example_inbreeding_ml = pyp_nrm.inbreeding(example,method='meu_luo')
    example_inbreeding_qu = pyp_nrm.inbreeding(example,method='mod_meu_luo')
    print 'Finished computing inbreeding at %s' % (pyp_nice_time())

    # print example_inbreeding
    # print example_inbreeding['fx'][28]
    print '\nVanRaden                  : ', example_inbreeding_vr['metadata']
    print '\nMeuwissen and Luo         : ', example_inbreeding_ml['metadata']
    print '\nModified Meuwissen and Luo: ', example_inbreeding_qu['metadata']

    # example_2 = pyp_newclasses.loadPedigree(optionsfile='new_inbreeding_2.ini')
    # example_inbreeding_2 = pyp_nrm.inbreeding(example_2)
    # print example_inbreeding_2
    # print 'f_x for 73543: ', example_inbreeding_2['fx'][example_2.idmap[int('73543')]]
