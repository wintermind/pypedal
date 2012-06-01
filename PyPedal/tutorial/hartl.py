#!/usr/bin/python

###############################################################################
# NAME: hartl.py
# VERSION: 2.0.0a4 (19APR2004)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from pyp_classes import *
from pyp_io import *
from pyp_metrics import *
from pyp_nrm import *
from pyp_utils import *

if __name__=='__main__':

	print 'Starting pypedal.py at %s' % asctime(localtime(time()))
	print '\tPreprocessing pedigree at %s' % asctime(localtime(time()))
	example = preprocess('hartl.ped',sepchar=' ')
	example = renumber(example,'example',io='yes')
	print '\tCalling set_ancestor_flag at %s' % asctime(localtime(time()))
	set_ancestor_flag(example,'example',io='yes')
	print '\tCollecting pedigree metadata at %s' % asctime(localtime(time()))
	example_meta = Pedigree(example,'example.ped','example_meta')
	example_meta.printme()
	print '\tCalling inbreeding() at %s' % asctime(localtime(time()))
	f_dict = inbreeding(example,'inbreeding','tabular')
	print f_dict
	print '\tCalling a_effective_founders_lacy() at %s' % asctime(localtime(time()))
	a_effective_founders_lacy(example,filetag='example')
	print '\tCalling a_effective_founders_boichard() at %s' % asctime(localtime(time()))
	a_effective_founders_boichard(example,filetag='example')
	print '\tCalling a_effective_ancestors_definite() at %s' % asctime(localtime(time()))
	a_effective_ancestors_definite(example,filetag='example')
	print '\tCalling a_effective_ancestors_indefinite() at %s' % asctime(localtime(time()))
	a_effective_ancestors_indefinite(example,filetag='example',n=10)
	print '\tCalling related_animals() at %s' % asctime(localtime(time()))
	list_a = related_animals(example[14].animalID,example)
	print list_a
	print '\tCalling related_animals() at %s' % asctime(localtime(time()))
	list_b = related_animals(example[9].animalID,example)
	print list_b
	print '\tCalling common_ancestors() at %s' % asctime(localtime(time()))
	list_r = common_ancestors(example[14].animalID,example[9].animalID,example)
	print list_r
	#print '\tPrinting pedigree at %s' % asctime(localtime(time()))
	#for e in example:
	#	e.printme()
	print '\tCalling relationship() at %s' % asctime(localtime(time()))
	_r = relationship(example[11].animalID,example[12].animalID,example,filetag='example')
	print _r
	print '\tCalling mating_coi() at %s' % asctime(localtime(time()))
	_f = mating_coi(example[11].animalID,example[12].animalID,example,filetag='example')
	print _f
	print 'Stopping pypedal.py at %s' % asctime(localtime(time()))
