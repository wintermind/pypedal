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

if __name__ == '__main__':

    print 'Starting pypedal.py at %s' % (pyp_nice_time())

    example = pyp_newclasses.loadPedigree(optionsfile='new_amatrix.ini')

    amatrix = pyp_newclasses.NewAMatrix(example.kw)
    amatrix.form_a_matrix(example.pedigree)

    # Here's how to save a matrix to a binary file.
#    amatrix.save('boichard2_pedigree.bin')

    # Here's how to load a matrix from a binary file.
#    amatrix2 = pyp_newclasses.NewAMatrix(example.kw)
#    amatrix2.load('boichard2_pedigree.bin')

    # Calculate coefficients of inbreeding on this pedigree.
    print '\tEntering pyp_nrm.inbreeding() at %s' % (pyp_nice_time())
    example_inbreeding = pyp_nrm.inbreeding(example,method='vanraden')
    print '\tReturning from pyp_nrm.inbreeding() at %s' % (pyp_nice_time())
    print example_inbreeding['metadata']
