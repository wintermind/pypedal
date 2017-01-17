#!/usr/bin/python
from __future__ import print_function
###############################################################################
# NAME: new_decompose.py
# VERSION: 2.0.0rc7 (07MAY2008)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import  pyp_newclasses
from PyPedal import  pyp_nrm
from PyPedal.pyp_utils import pyp_nice_time

import numpy
import numpy.linalg
import time

if __name__ == '__main__':

    print('Starting pypedal.py at %s' % (pyp_nice_time()))
    example = pyp_newclasses.loadPedigree(optionsfile='new_decompose.ini')
    amatrix = pyp_newclasses.NewAMatrix(example.kw)
    amatrix.form_a_matrix(example.pedigree)

    print('-'*80)
    print('Calling a_decompose()')
    print('=====================')
    D, T  = pyp_nrm.a_decompose(example))
    print('D: ', D)
    print('T: ', T)

    print('-'*80)
    print('Calculating Ainv from D and T')
    print('============================='
    l = example.metadata.num_records
    Tinv = numpy.linalg.inv(T)
    print('Tinv: ', Tinv)
    Tpinv = numpy.linalg.inv(T.T)
    print('Tpinv: ', Tpinv)
    Dinv = numpy.linalg.inv(D)
    print('Dinv: ', Dinv)
    Ainvhalf = numpy.dot(Tpinv,Dinv)
    Ainv = numpy.dot(Ainvhalf,Tinv)
    print('Ainv: ', Ainv)

    print('-'*80)
    print('Calling form_d_nof()')
    print('====================')
    D  = pyp_nrm.form_d_nof(example)
    print('D: ', D)

    print('-'*80)
    print('Calling a_inverse_dnf()')
    print('=======================')
    Ainv = pyp_nrm.a_inverse_dnf(example)
    print('Ainv: ', Ainv)

    print('-'*80)
    print('Calling a_inverse_df()')
    print('======================')
    Ainv = pyp_nrm.a_inverse_df(example)
    print('Ainv: ', Ainv)

    print('Stopping pypedal.py at %s' % (pyp_nice_time()))
