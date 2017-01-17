#!/usr/bin/python
from __future__ import print_function
###############################################################################
# NAME: new_methods.py
# VERSION: 2.0.0b5 (19DECEMBER2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_newclasses
from PyPedal import pyp_metrics
from PyPedal import pyp_nrm
from PyPedal import pyp_utils
from PyPedal.pyp_utils import pyp_nice_time

if __name__=='__main__':

    print('Starting pypedal.py at %s' % ( pyp_nice_time() ))

    example = pyp_newclasses.loadPedigree(optionsfile='new_format.ini')
    #pyp_utils.assign_offspring(example)
    #print('First ID: %s' % ( example.kw['id_first'] ))
    #print('Last ID: %s' % ( example.kw['id_last'] ))

    inbr, reln = inbr, reln = pyp_nrm.inbreeding(example,method='vanraden',rels=1)
    print('inbr: ', inbr)
    print('reln: ', reln)

    example.nrm = pyp_newclasses.NewAMatrix(example.kw)
    example.nrm.form_a_matrix(example.pedigree)
    example.nrm.save('Amatrix.txt')

    inbr2, reln2 = pyp_nrm.inbreeding(example,method='vanraden',rels=1)
    print('inbr2: ', inbr2)
    print('reln2: ', reln2)

    example.nrm2 = pyp_newclasses.NewAMatrix(example.kw)
    example.nrm2.load('Amatrix.txt')
    #example.nrm2.printme()

    print(example.nrm.nrm[1][4])
    print(example.nrm.nrm[4][1])

    thesame = (example.nrm == example.nrm2)
    print('thesame: ', thesame)

    print('\tCalling related_animals() at %s' % ( pyp_nice_time() ))
    list_a = pyp_metrics.related_animals(example.pedigree[4].animalID,example)
    print(list_a)

    print('\tCalling related_animals() at %s' % ( pyp_nice_time() ))
    list_b = pyp_metrics.related_animals(example.pedigree[13].animalID,example)
    print(list_b)

    print('\tCalling common_ancestors() at %s' % ( pyp_nice_time() ))
    list_r = pyp_metrics.common_ancestors(example.pedigree[4].animalID,example.pedigree[13].animalID,example)
    print(list_r)

    #print('\tTesting NewPedigree::addanimal() at %s' % ( pyp_nice_time() ))
    #_added = example.addanimal(15,10,11)
    #if _added:
    #    print('\t\tAdded animal 15 to %s' % ( example.kw['pedname'] ))

    #print('\tTesting NewPedigree::delanimal() at %s' % ( pyp_nice_time() ))
    #_deleted = example.delanimal(15)
    #if _deleted:
    #    print('\t\tDeleted animal 15 from %s' % ( example.kw['pedname'] ))

    print('\tCalling mating_coi() at %s' % ( pyp_nice_time() ))
    f = pyp_metrics.mating_coi(example.pedigree[0].animalID,example.pedigree[5].animalID,example)
    print(f)

    print('\tCalling mating_coi() at %s' % ( pyp_nice_time() ))
    f = pyp_metrics.mating_coi(example.pedigree[0].animalID,example.pedigree[13].animalID,example,0)
    print(f)

    print('\tCalling mating_coi() at %s' % ( pyp_nice_time() ))
    f = pyp_metrics.mating_coi(example.pedigree[4].animalID,example.pedigree[13].animalID,example,1)
    print(f)

    matings = []
    matings.append('%s_%s'%(example.pedigree[0].animalID, example.pedigree[4].animalID))
    matings.append('%s_%s'%(example.pedigree[0].animalID, example.pedigree[13].animalID))
    matings.append('%s_%s'%(example.pedigree[4].animalID, example.pedigree[13].animalID))
    fgrp = pyp_metrics.mating_coi_group(matings,example)
    print('fgrp: ', fgrp['matings'])

    desc1 = pyp_metrics.descendants(5,example,{})
    print('desc1: ', desc1)

    descf = pyp_metrics.founder_descendants(example)
    print('descf: ', descf)

    print('Stopping pypedal.py at %s' % ( pyp_nice_time() ))
