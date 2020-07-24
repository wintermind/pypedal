#!/usr/bin/python

###############################################################################
# NAME: new_inbreeding2.py
# VERSION: 2.0.0b5 (15DECEMBER2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import  pyp_graphics
from PyPedal import  pyp_network
from PyPedal import  pyp_newclasses
from PyPedal import  pyp_nrm

if __name__ == '__main__':

#     example = pyp_newclasses.loadPedigree(optionsfile='new_inbreeding2.ini')
#
#     ng = pyp_network.ped_to_graph(example)
#     _1g = pyp_network.find_ancestors_g(ng, example.idmap['1'], [], 1)
#     _1g.append(example.idmap['1'])
#     print '1 gen pedigree: %s' % ( _1g )
#     _2g = pyp_network.find_ancestors_g(ng, example.idmap['1'], [], 2)
#     _2g.append(example.idmap['1'])
#     print '2 gen pedigree: %s' % ( _2g )
#     _3g = pyp_network.find_ancestors_g(ng, example.idmap['1'], [], 3)
#     _3g.append(example.idmap['1'])
#     print '3 gen pedigree: %s' % ( _3g )
#     _4g = pyp_network.find_ancestors_g(ng, example.idmap['1'], [], 4)
#     _4g.append(example.idmap['1'])
#     print '4 gen pedigree: %s' % ( _4g )
#
# # This example tests pyp_nrm.inbreeding_vanraden()
#     example_inbreeding = pyp_nrm.inbreeding(example, method='vanraden', gens=3)
#     print '\n3 gen inbreeding: %s' % ( example_inbreeding )
#     print '\nMean 3 gen CoI: %s' % ( example_inbreeding['metadata']['all']['f_avg'] )
#
# # This example tests pyp_nrm.inbreeding_tabular()
#     example_inbreeding = pyp_nrm.inbreeding(example, method='tabular', gens=3)
#     print '\n3 gen inbreeding: %s' % ( example_inbreeding )
#     print '\nMean 3 gen CoI: %s' % ( example_inbreeding['metadata']['all']['f_avg'] )

# This is an advanced example that uses dict4ini to read
# configurations for multiple pedigrees from the same file.
    from PyPedal import dict4ini
    import sys
    kw_all = dict4ini.DictIni('new_inbreeding2multiple.ini')
    if len(kw_all) == 0:
        sys.exit(0)
    kw_all = kw_all.dict()
    kw_noinbred = kw_all['noinbreeding']
    kw_horse = kw_all['horse']

#     noinbred = pyp_newclasses.loadPedigree(kw_noinbred)
    horse = pyp_newclasses.loadPedigree(kw_horse)

#     noinbred_inbreeding = pyp_nrm.inbreeding(noinbred, method='vanraden', gens=3)
#     print '\nNoinbred 3 gen inbreeding: %s' % ( noinbred_inbreeding )
#     print '\nNoinbred mean 3 gen CoI: %s' % ( noinbred_inbreeding['metadata']['all']['f_avg'] )
#
#     pyp_graphics.draw_pedigree(horse, gfilename='blue_joe',
#         gtitle='Blue Joe\'s Pedigree', gorient='p', gname=1, gdirec='RL',
#         gfontsize=12, garrow=0, gtitloc='b')

##    print '\nRen ID\tName\n------\t----'
##    for _k in horse.backmap.keys():
##        print '%s\t%s' % ( _k, horse.pedigree[int(_k)-1].name )

##    print '-'*80

##    horse_inbreeding = pyp_nrm.inbreeding(horse, method='vanraden', gens=2)
##    print '\nHorse 2 gen inbreeding: %s' % ( horse_inbreeding )
##    print '\nHorse mean 2 gen CoI: %s' % ( horse_inbreeding['metadata']['all']['f_avg'] )

    horse_inbreeding, horse_rels = pyp_nrm.inbreeding(horse, method='tabular', rels=1)
    print '\nHorse inbreeding: %s' % ( horse_inbreeding )
    print '\nHorses wth non-zero coefficients of inbreeding:'
    print '-----------------------------------------------'
    for _h in horse_inbreeding['fx'].keys():
        if horse_inbreeding['fx'][_h] > 0.:
            print '%s\t\t%s' % ( horse.pedigree[int(_h)-1].name, horse_inbreeding['fx'][_h] )
    print '\nHorse mean CoI:\t\t%s' % ( horse_inbreeding['metadata']['all']['f_avg'] )

##    print '-'*80
##
##    hg = pyp_network.ped_to_graph(horse)
##    _2g = pyp_network.find_ancestors_g(hg, len(horse.idmap), {}, 2)
##    _2g[len(horse.idmap)] = 1
##    _names = []
##    for _h in _2g.keys():
##        _names.append(horse.pedigree[int(_h)-1].name)
##    print '2 gen pedigree: %s' % ( _names )

##    print '-'*80
##
##    _allg = pyp_network.find_ancestors(hg, len(horse.idmap), [])
##    _allg.append(len(horse.idmap))
##    _names = []
##    for _h in _allg:
##        _names.append(horse.pedigree[int(_h)-1].name)
##    print 'All gen pedigree: %s' % ( _names )
##
##    print '-'*80