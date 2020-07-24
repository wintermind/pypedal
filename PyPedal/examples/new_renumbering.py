#!/usr/bin/python

###############################################################################
# NAME: new_graphics.py
# VERSION: 2.0.0b5 (19DECEMBER2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_graphics
from PyPedal import pyp_newclasses
from PyPedal import pyp_nrm

if __name__ == '__main__':

    example = pyp_newclasses.loadPedigree(optionsfile='new_renumbering.ini')

    print '-'*80
    print 'INFO: Pedigree ID map: %s' % (example.idmap)
    print '-'*80
    print 'INFO: Pedigree ID reverse map: %s' % (example.backmap)
    print '-'*80
    myinbr = pyp_nrm.inbreeding(example)
    print 'INFO: Coefficients of inbreeding: %s' % (myinbr)

    pyp_graphics.draw_pedigree(example,gfilename='new_renumbering',gtitle='My  Pedigree')