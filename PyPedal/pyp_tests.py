#!/usr/bin/python

###############################################################################
# NAME: pyp_tests.py
# VERSION: 2.0.0a18 (19JUL2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################
#   Provide unit tests for PyPedal.
###############################################################################

##
# pyp_tests contains a series of unit tests for verifying that PyPedal works
# correctly.
##
from __future__ import print_function
import pyp_demog
import pyp_metrics
import pyp_newclasses
import pyp_nrm
import pyp_utils
import logging, numpy, os, string, sys, time
import unittest

class PyPedalMetricsTestCases(unittest.TestCase):

    def testMetricsMinMaxF(self):
        # Pedigree from van Noordwijck and Scharloo (1981) as presented
        # in Hartl and Clark (1989), p. 242.
        options = {}
        options['renumber'] = 0
        options['messages'] = 'quiet'
        options['pedfile'] = 'tests/hartlandclark.ped'
        options['pedformat'] = 'asd'
        options['sepchar'] = ' '
        options['pedigree_is_renumbered'] = '1'
        options['form_nrm'] = '1'
        example = pyp_newclasses.NewPedigree(options)
        example.load()
        #example.nrm.info()

        high_coi, low_coi = pyp_metrics.min_max_f(example,n=5)
        print(high_coi)
        print(low_coi)
        #self.assertEqual(2.91, round(fe,2))

    def testMetricsEffectiveFoundersLacy(self):
        # Example taken from Lacy (1989), Appendix A.
        # Used to test calculation of effective population size.
        options = {}
        options['renumber'] = 0
        options['messages'] = 'quiet'
        options['pedfile'] = 'tests/lacy/lacy.ped'
        options['pedformat'] = 'asd'
        options['sepchar'] = ' '
        options['pedigree_is_renumbered'] = '1'
        example = pyp_newclasses.NewPedigree(options)
        example.load()
        fe = pyp_metrics.effective_founders_lacy(example)
        self.assertEqual(2.91, round(fe,2))

#     def testMetricsEffectiveFounderGenomes(self):
#         # Example taken from Lacy (1989), Appendix A.
#         # Used to test calculation of founder genome equivalents.
#         options = {}
#         options['renumber'] = 0
#         options['messages'] = 'quiet'
#         options['pedfile'] = 'tests/lacy/lacy.ped'
#         options['pedformat'] = 'asd'
#         options['sepchar'] = ' '
#         options['pedigree_is_renumbered'] = '1'
#         options['set_alleles'] = 1
#         example = pyp_newclasses.NewPedigree(options)
#         example.load()
#         fg = pyp_metrics.effective_founder_genomes(example,rounds=250)
#         self.assertEqual(round(2.18,2), round(fg,2))

    def testMetricsEffectiveFoundersBoichardA(self):
        # Example taken from Boichard et al. (1997), Figure 2.
        # Used to test calculation of effective founder number.
        options = {}
        options['messages'] = 'quiet'
        options['pedfile'] = 'tests/boichard/boichard2a.ped'
        options['pedname'] = 'Boichard Pedigree (Family 1 only)'
        options['pedformat'] = 'asdg'
        options['pedigree_is_renumbered'] = '1'
        options['filetag'] = 'example2a'
        example2a = pyp_newclasses.NewPedigree(options)
        example2a.load()
        fe = pyp_metrics.a_effective_founders_boichard(example2a)
        self.assertEqual(round(4.0,1), round(fe,1))

    def testMetricsEffectiveFoundersBoichardB(self):
        # Example taken from Boichard et al. (1997), Figure 2.
        # Used to test calculation of effective founder number.
        options = {}
        options['messages'] = 'quiet'
        options['pedfile'] = 'tests/boichard/boichard2b.ped'
        options['pedname'] = 'Boichard Pedigree (Family 2 only)'
        options['pedformat'] = 'asdg'
        options['pedigree_is_renumbered'] = '1'
        options['filetag'] = 'example2b'
        example2b = pyp_newclasses.NewPedigree(options)
        example2b.load()
        fe = pyp_metrics.a_effective_founders_boichard(example2b)
        self.assertEqual(round(2.0,1), round(fe,1))

    def testMetricsEffectiveFoundersBoichardC(self):
        # Example taken from Boichard et al. (1997), Figure 2.
        # Used to test calculation of effective founder number.
        options = {}
        options['messages'] = 'quiet'
        options['pedfile'] = 'tests/boichard/boichard2.ped'
        options['pedname'] = 'Boichard Pedigree (Family 1 and 2)'
        options['pedformat'] = 'asdg'
        options['pedigree_is_renumbered'] = '1'
        options['filetag'] = 'example2'
        example2 = pyp_newclasses.NewPedigree(options)
        example2.load()
        fe = pyp_metrics.a_effective_founders_boichard(example2)
        self.assertEqual(round(5.6,1), round(fe,1))

    def testMetricsEffectiveAncestorsDefiniteBoichardA(self):
        # Example taken from Boichard et al. (1997), Figure 2.
        # Used to test calculation of effective ancestor number.
        options = {}
        options['messages'] = 'quiet'
        options['pedfile'] = 'tests/boichard/boichard2a.ped'
        options['pedname'] = 'Boichard Pedigree (Family 1 only)'
        options['pedformat'] = 'asdg'
        options['pedigree_is_renumbered'] = '1'
        options['filetag'] = 'example2a'
        example2a = pyp_newclasses.NewPedigree(options)
        example2a.load()
        fa = pyp_metrics.a_effective_ancestors_definite(example2a)
        self.assertEqual(round(2.0,2), round(fa,2))

    def testMetricsEffectiveAncestorsDefiniteBoichardB(self):
        # Example taken from Boichard et al. (1997), Figure 2.
        # Used to test calculation of effective ancestor number.
        options = {}
        options['messages'] = 'quiet'
        options['pedfile'] = 'tests/boichard/boichard2b.ped'
        options['pedname'] = 'Boichard Pedigree (Family 2 only)'
        options['pedformat'] = 'asdg'
        options['pedigree_is_renumbered'] = '1'
        options['filetag'] = 'example2a'
        example2b = pyp_newclasses.NewPedigree(options)
        example2b.load()
        fa = pyp_metrics.a_effective_ancestors_definite(example2b)
        self.assertEqual(round(2.0,2), round(fa,2))

    def testMetricsEffectiveAncestorsDefiniteBoichardC(self):
        # Example taken from Boichard et al. (1997), Figure 2.
        # Used to test calculation of effective ancestor number.
        options = {}
        options['messages'] = 'quiet'
        options['pedfile'] = 'tests/boichard/boichard2.ped'
        options['pedname'] = 'Boichard Pedigree (Family 1 and 2)'
        options['pedformat'] = 'asdg'
        options['pedigree_is_renumbered'] = '1'
        options['filetag'] = 'example2'
        example2 = pyp_newclasses.NewPedigree(options)
        example2.load()
        fa = pyp_metrics.a_effective_ancestors_definite(example2)
        self.assertEqual(round(2.94,2), round(fa,2))

    def testMetricsEffectiveAncestorsIndefiniteBoichardA(self):
        # This is just a stub because I do not yet have a good way to test this routine.
        pass

class PyPedalNrmTestCases(unittest.TestCase):
    pass

class PyPedalUtilsTestCases(unittest.TestCase):
    pass

if __name__ == '__main__':
    try:
        import testoob
        testoob.main()
    except ImportError:
        print('Could not import testoob module (https://opensvn.csie.org/traccgi/testoob/trac.cgi/wiki)!')
        sys.exit(0)
