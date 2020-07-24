#!/usr/bin/python

###############################################################################
# NAME: new_hartl.py
# VERSION: 2.0.0b12 (15MAY2060)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_graphics
from PyPedal import pyp_newclasses
from PyPedal import pyp_nrm

if __name__ == '__main__':

    example = pyp_newclasses.loadPedigree(optionsfile='new_hartl.ini')

    pyp_graphics.draw_pedigree(example, gfilename='hartlandclark',
        gtitle='Pedigree from van Noordwijck and Scharloo (1981)',
        gorient='p', gname=0, gdirec='', gfontsize=12, garrow=0)

    example_inbreeding = pyp_nrm.inbreeding(example)
    print example_inbreeding

	
