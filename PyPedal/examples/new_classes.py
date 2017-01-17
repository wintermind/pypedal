#!/usr/bin/python
from __future__ import print_function
###############################################################################
# NAME: new_classes.py
# VERSION: 2.0.0b5 (13DECEMBER2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_graphics
from PyPedal import pyp_newclasses
from PyPedal import pyp_nrm
from PyPedal.pyp_utils import pyp_nice_time

if __name__ == '__main__':

    print('Starting pypedal.py at %s' % ( pyp_nice_time() ))

    example = pyp_newclasses.loadPedigree(optionsfile='new_classes.ini')

    if example.kw['messages'] == 'verbose':
        print('[INFO]: Forming numerator relationship matrix at %s' % ( pyp_nice_time() ))

    my_a = pyp_nrm.fast_a_matrix_r(example.pedigree,example.kw)

    if example.kw['messages'] == 'verbose':
        print('[INFO]: Visualizing NRM sparsity at %s' % ( pyp_nice_time() ))

    pyp_graphics.rmuller_spy_matrix_pil(my_a,fname='boichard2_spy.png')

    if example.kw['messages'] == 'verbose':
        print('[INFO]: Visualizing NRM in pseudocolor at %s' % ( pyp_nice_time() ))

    pyp_graphics.rmuller_pcolor_matrix_pil(my_a,fname='boichard2_pcolor.png')

    pyp_graphics.draw_pedigree(example, gfilename='boichard2_pedigree', \
        gtitle='', gorient='p', gdirec='RL', gfontsize=12)

    print('Stopping pypedal.py at %s' % ( pyp_nice_time() ))
