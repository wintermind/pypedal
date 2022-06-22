#!/usr/bin/python

###############################################################################
# NAME: new_graphics3.py
# VERSION: 2.0.0b15 (18SEPTEMBER2006)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_graphics
from PyPedal import pyp_utils
from PyPedal import pyp_newclasses
from PyPedal import pyp_metrics

options = {
        'pedname': 'Test Duplicate Handling',
        'pedformat': 'asd',
        'pedfile': 'duplicates.ped',
        'sepchar': ",",
        'messages': 'verbose',
        'debug': True,
	'assign_sexes': True,
	'renumber': 1,
	'resolve_duplicates': 1,
}

if __name__ == '__main__':

	example = pyp_newclasses.loadPedigree(options=options, debugLoad=True)

        print
        print 'List of duplciate animals'
        print '-------------------------'
	for d in example.duplicates:
            print '\t', d

        print
        print 'Drawing pedigree with duplicate animals'
	gtitle = 'Pedigree of duplicate animals'
	pyp_graphics.new_draw_pedigree(example, gfilename='duplicates', gtitle=gtitle,
	    gformat='jpg', gorient='l', gdirec='RL', gtitloc='t', gpenwidth=2, \
	    gfontsize=14)
