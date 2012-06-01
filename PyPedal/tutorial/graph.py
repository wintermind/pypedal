#!/usr/bin/python

###############################################################################
# NAME: graph.py
# VERSION: 2.0.0a5 (20APR2004)
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
	print '\tDrawing pedigree at %s' % asctime(localtime(time()))
	draw_pedigree(example,gfilename='hartl',gtitle='Hartl_and_Clark_Great_Tit_Pedigree')
	print '='*80
	print '\tPreprocessing pedigree at %s' % asctime(localtime(time()))
	example = preprocess('../examples/boichard2.ped',sepchar=' ')
	example = renumber(example,'example',io='yes')
	print '\tDrawing pedigree at %s' % asctime(localtime(time()))
	draw_pedigree(example,gfilename='boichard',gtitle='Boichard_Pedigree')
	print 'Stopping pypedal.py at %s' % asctime(localtime(time()))
