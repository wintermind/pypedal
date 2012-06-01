#!/usr/bin/python

###############################################################################
# NAME: new_ids.py
# VERSION: 2.0.0b5 (14DECEMBER2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_newclasses

if __name__ == '__main__':

    example = pyp_newclasses.loadPedigree(optionsfile='new_ids.ini')

    for _p in example.pedigree:
        _p.printme()