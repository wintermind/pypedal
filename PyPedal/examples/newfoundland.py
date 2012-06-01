#!/usr/bin/python

###############################################################################
# NAME: new_reporting.py
# VERSION: 2.0.0b15 (18SE{SEPTEMBET2006)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_graphics
from PyPedal import pyp_jbc
from PyPedal import pyp_network
from PyPedal import pyp_newclasses
from PyPedal import pyp_nrm
from PyPedal import pyp_metrics

if __name__ == '__main__':

    example = pyp_newclasses.loadPedigree(optionsfile='newfoundland.ini')

    print example.idmap
    print example.backmap
    print example.namemap
    print example.namebackmap

    newf_f = pyp_nrm.inbreeding(example)
    print newf_f['fx'][example.idmap[example.namemap['Kaptn Kvols von Widdersdorf']]]

##    dussel_id = example.idmap[example.namemap['King von der Dussel']]
##    print 'dussel_id: ', dussel_id
##
##    print 'Empirical proof that ancestor loss coefficients are the same'
##    print 'as pedigree completeness.'
##
##    pedcomp3 = pyp_metrics.pedigree_completeness(example,gens=3)
##    print example.pedigree[dussel_id-1].pedcomp
##    pedcomp4 = pyp_metrics.pedigree_completeness(example,gens=4)
##    print example.pedigree[dussel_id-1].pedcomp
##    pedcomp5 = pyp_metrics.pedigree_completeness(example,gens=5)
##    print example.pedigree[dussel_id-1].pedcomp

#    pyp_jbc.color_pedigree(example,gfilename='newfoundland', ghatch='0', \
#        gtitle='Nodes are colored by number of descendantss.', \
#        gprog='dot', gname=1, gdirec='RL', gformat='eps')
