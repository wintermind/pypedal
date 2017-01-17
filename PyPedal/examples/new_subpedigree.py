#!/usr/bin/python
from __future__ import print_function
###############################################################################
# NAME: new_subpedigree.py
# VERSION: 2.0.0b14 (16MAY2006)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_graphics
from PyPedal import pyp_newclasses
from PyPedal import pyp_nrm
from PyPedal import pyp_utils
from PyPedal.pyp_utils import pyp_nice_time

if __name__ == '__main__':

    example = pyp_newclasses.loadPedigree(optionsfile='new_subpedigree.ini')
    print('='*70)
    keepid = []
    for k in ['Clover','Dewey','Wicket']:
        ped = pyp_nrm.recurse_pedigree(example,example.idmap[example.namemap[k]],[])
        for p in ped:
            if p.animalID not in keepid:
                keepid.append(p.animalID)
        if example.idmap[example.namemap[k]] not in keepid:
            keepid.append(example.idmap[example.namemap[k]])
    example2 = pyp_utils.subpedigree(example,keepid)
    example2.metadata.printme()
