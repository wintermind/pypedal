#!/usr/bin/python

###############################################################################
# NAME: new_graphics3.py
# VERSION: 2.0.0b15 (18SEPTEMBER2006)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_demog
from PyPedal import pyp_graphics
from PyPedal import pyp_jbc
from PyPedal import pyp_newclasses
from PyPedal import pyp_nrm
from PyPedal import pyp_metrics
from PyPedal.pyp_utils import pyp_nice_time

if __name__ == '__main__':

    print 'Starting pypedal.py at %s' % (pyp_nice_time())

    example = pyp_newclasses.loadPedigree(optionsfile='new_graphics3.ini')
    if example.kw['messages'] == 'verbose':
        print '[INFO]: Calling pyp_graphics.new_draw_pedigree() at %s' % (pyp_nice_time())
    pyp_graphics.new_draw_pedigree(example, gfilename='graphics3', gtitle='graphics3 pedigree', gorient='p')

    pyp_jbc.color_pedigree(example,gfilename='graphics3', ghatch='0', \
        metric='sons', gtitle='Nodes are colored by number of sons.', \
        gprog='dot', gname=1)