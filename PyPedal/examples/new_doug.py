#!/usr/bin/python

###############################################################################
# NAME: new_doug.py
# VERSION: 2.0.0b5 (13DECEMBER2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_graphics
from PyPedal import pyp_newclasses
from PyPedal.pyp_utils import pyp_nice_time

if __name__ == '__main__':

    example = pyp_newclasses.loadPedigree(optionsfile='new_doug.ini')

    if example.kw['messages'] == 'verbose':
        print '[INFO]: Calling pyp_graphics.draw_pedigree() at %s' % (pyp_nice_time())

    pyp_graphics.draw_pedigree(example, gfilename='doug_below',
        gtitle='Doug the German Shepherd (B)', gorient='p', gname=1, gdirec='',
        gfontsize=12, garrow=0, gtitloc='b')

    pyp_graphics.draw_pedigree(example, gfilename='doug_above',
        gtitle='Doug the German Shepherd (A)', gorient='p', gname=1, gdirec='',
        gfontsize=12, garrow=0, gtitloc='t')

    pyp_graphics.draw_pedigree(example, gfilename='doug_p_rl_notitle',
        gtitle='', gorient='p', gname=1, gdirec='RL', gfontsize=12)