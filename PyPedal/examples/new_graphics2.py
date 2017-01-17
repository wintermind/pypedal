#!/usr/bin/python
from __future__ import print_function
###############################################################################
# NAME: new_graphics2.py
# VERSION: 2.0.0b10 (27APRIL2006)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_demog
from PyPedal import pyp_graphics
from PyPedal import pyp_newclasses
from PyPedal import pyp_nrm
from PyPedal import pyp_metrics
from PyPedal.pyp_utils import pyp_nice_time
#import pyp_demog, pyp_graphics, pyp_newclasses, pyp_nrm, pyp_metrics
#from pyp_utils import pyp_nice_time

if __name__ == '__main__':

    print('Starting pypedal.py at %s' % (pyp_nice_time()))

    example = pyp_newclasses.loadPedigree(optionsfile='new_graphics.ini')
    if example.kw['messages'] == 'verbose':
        print('[INFO]: Calling pyp_graphics.draw_pedigree() at %s' % (pyp_nice_time()))
    pyp_graphics.draw_pedigree(example, gfilename='graphics2', gtitle='graphics2 pedigree', gorient='p')
