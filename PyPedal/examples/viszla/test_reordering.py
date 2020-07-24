#!/usr/bin/python

###############################################################################
# NAME: new_amatrix.py
# VERSION: 2.0.0rc7 (07MAY2008)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import  pyp_newclasses
from PyPedal import  pyp_nrm
from PyPedal.pyp_utils import pyp_nice_time

import numpy
import time

options = {}
options['messages'] = 'verbose'
options['pedfile'] = 'vizsla_ped3.csv'
options['pedformat'] = 'asd'
options['reorder'] = 1
options['reorder_max_rounds'] = 1000
options['has_header'] = 1
options['sepchar'] = ','
options['pedname'] = 'Viszla Pedigree'
options['debug_messages'] = 0

if __name__ == '__main__':

    print 'PyPedal started at %s' % (pyp_nice_time())
    viszla = pyp_newclasses.loadPedigree(options, debugLoad=True)
    print 'PyPedal stopped at %s' % (pyp_nice_time())
