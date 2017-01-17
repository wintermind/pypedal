#!/usr/bin/python
from __future__ import print_function
###############################################################################
# NAME: new_lacy.py
# VERSION: 2.0.0b5 (14DECEMBER2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_newclasses
from PyPedal import pyp_nrm
from PyPedal import pyp_metrics
from PyPedal.pyp_utils import pyp_nice_time

if __name__ == '__main__':

    print('Starting pypedal.py at %s' % (pyp_nice_time()))

# Example taken from Lacy (1989), Appendix A.
    example = pyp_newclasses.loadPedigree(optionsfile='new_lacy.ini',debugLoad=True)
    example.printoptions()
    if example.kw['messages'] == 'verbose':
        print('='*80)
        print('[INFO]: Calling pyp_metrics.effective_founders_lacy at %s' % (pyp_nice_time()))
    pyp_metrics.effective_founders_lacy(example)

# Example taken from Boichard et al. (1997), Figure 2 / Table II.
    example2a = pyp_newclasses.loadPedigree(optionsfile='new_format.ini')
    if example2a.kw['messages'] == 'verbose':
        print('='*80)
        print('[INFO]: Calling pyp_metrics.effective_founders_lacy at %s' % (pyp_nice_time()))
    pyp_metrics.effective_founders_lacy(example2a)
    if example2a.kw['messages'] == 'verbose':
        print('='*80)
        print('Calling a_effective_founders_boichard() at %s' % (pyp_nice_time()))
    pyp_metrics.a_effective_founders_boichard(example2a)
    if example2a.kw['messages'] == 'verbose':
        print('Calling a_effective_ancestors_definite() at %s' % (pyp_nice_time()))
        print('='*80)
    pyp_metrics.a_effective_ancestors_definite(example2a)

# # Example taken from Boichard et al. (1997), Figure 2 / Table II.
#     options['pedfile'] = 'boichard2b.ped'
#     options['pedname'] = 'Boichard Pedigree (Family B only)'
#     options['pedformat'] = 'asdg'
#     example2b = pyp_newclasses.NewPedigree(options)
#     example2b.load()
#
#     print('='*80)
#
#     if example2b.kw['messages'] == 'verbose':
#         print('[INFO]: Calling pyp_metrics.effective_founders_lacy at %s' % (pyp_nice_time()))
#     pyp_metrics.effective_founders_lacy(example2b)
#
#     print('='*80)
#
#     if example2b.kw['messages'] == 'verbose':
#         print('Calling a_effective_founders_boichard() at %s' % (pyp_nice_time()))
#     pyp_metrics.a_effective_founders_boichard(example2b)
#
#     print('='*80)
#
#     if example2b.kw['messages'] == 'verbose':
#         print('Calling a_effective_ancestors_definite() at %s' % (pyp_nice_time()))
#     pyp_metrics.a_effective_ancestors_definite(example2b)
#
#     print('='*80)
#
# # Example taken from Boichard et al. (1997), Figure 2 / Table II.
#     options['pedfile'] = 'boichard2.ped'
#     options['pedname'] = 'Boichard Pedigree'
#     options['pedformat'] = 'asdg'
#     example2 = pyp_newclasses.NewPedigree(options)
#     example2.load()
#
#     print('='*80)
#
#     if example2.kw['messages'] == 'verbose':
#         print('[INFO]: Calling pyp_metrics.effective_founders_lacy at %s' % (pyp_nice_time()))
#     pyp_metrics.effective_founders_lacy(example2)
#
#     print('='*80)
#
#     if example2.kw['messages'] == 'verbose':
#         print('Calling a_effective_founders_boichard() at %s' % (pyp_nice_time()))
#     pyp_metrics.a_effective_founders_boichard(example2)
#
#     print('='*80)
#
#     if example2.kw['messages'] == 'verbose':
#         print('Calling a_effective_ancestors_definite() at %s' % (pyp_nice_time()))
#     pyp_metrics.a_effective_ancestors_definite(example2)
#
#     print('='*80)
#
# # Example taken from Boichard et al. (1997), Figure 1 / Table I.
#     options['pedfile'] = 'boichard1.ped'
#     options['pedname'] = 'Boichard Pedigree'
#     options['pedformat'] = 'asdg'
#     example1 = pyp_newclasses.NewPedigree(options)
#     example1.load()
#
#     print('='*80)
#
#     if example1.kw['messages'] == 'verbose':
#         print('[INFO]: Calling pyp_metrics.effective_founders_lacy at %s' % (pyp_nice_time()))
#     pyp_metrics.effective_founders_lacy(example1)
#
#     print('='*80)
#
#     if example1.kw['messages'] == 'verbose':
#         print('Calling a_effective_founders_boichard() at %s' % (pyp_nice_time()))
#     pyp_metrics.a_effective_founders_boichard(example1)
#
#     print('='*80)
#
#     if example1.kw['messages'] == 'verbose':
#         print('Calling a_effective_ancestors_definite() at %s' % (pyp_nice_time()))
#     pyp_metrics.a_effective_ancestors_definite(example1)
#
#     print('='*80)

    print('Stopping pypedal.py at %s' % (pyp_nice_time()))
