#!/usr/bin/python

###############################################################################
# NAME: new_jbc.py
# VERSION: 2.0.0b5 (14DECEMBER2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_graphics
from PyPedal import pyp_jbc
from PyPedal import pyp_newclasses

if __name__=='__main__':

    example = pyp_newclasses.loadPedigree(optionsfile='new_jbc.ini')

    pyp_jbc.color_pedigree(example,gfilename='BoichardPedigree_colored_dot', \
        ghatch='dna', metric='sons', gtitle='Nodes are colored by number of sons.')