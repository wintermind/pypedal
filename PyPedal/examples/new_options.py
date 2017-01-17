#!/usr/bin/python
from __future__ import print_function
###############################################################################
# NAME: new_methods.py
# VERSION: 2.0.0a5 (12DECEMBER2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_newclasses
from PyPedal.pyp_utils import pyp_nice_time

if __name__=='__main__':

    print('Starting pypedal.py at %s' % (pyp_nice_time()))

    print('='*80)

    # This should fail because we are providing neither a dictionary
    # of options or a configuration file name.
    print('This load should fail')
    myped1 = pyp_newclasses.loadPedigree()
    print(myped1)

    print('-'*80)

    # This should fail because we are providing an empty dictionary
    # and there is no configuration file named pypedal.ini in the
    # examples directory.
    options = {}
    print('This load should fail')
    myped2 = pyp_newclasses.loadPedigree(options)
    print(myped2)

    print('-'*80)

    # This should work because we are providing a dictionary of
    # options and no configuration file name.
    options = {}
    options['pedfile'] = 'new_lacy.ped'
    options['pedformat'] = 'asd'
    options['pedname'] = 'Lacy Pedigree'
    myped3 = pyp_newclasses.loadPedigree(options)
    print(myped3)

    print('-'*80)

    # This should work because we are providing an empty dictionary
    # and the name of a valid configuration file.
    options = {}
    myped4 = pyp_newclasses.loadPedigree(optionsfile='new_options.ini')
    print(myped4)

    print('='*80)

    print('Stopping pypedal.py at %s' % (pyp_nice_time()))
