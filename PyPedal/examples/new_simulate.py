#!/usr/bin/python
from __future__ import print_function
###############################################################################
# NAME: new_simulate.py
# VERSION: 2.0.0b15 (08JUNE2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import  pyp_newclasses
from PyPedal import  pyp_graphics
from PyPedal import  pyp_network
import networkx

if __name__ == '__main__':

    for i in range(1):
        print('Simulating and drawing pedigree %d' % ( i ))
        example = pyp_newclasses.loadPedigree(optionsfile='new_simulate.ini')
        print(example)

        #_gfn = 'simulated_pedigree_%d' % ( i )
        _gfn = 'simulated_pedigree'
        _gt = 'Simulated Pedigree %d' % ( i )
        pyp_graphics.draw_pedigree(example, gfilename=_gfn, gtitle=_gt, gformat='ps', garrow=0, gclusters=0)

        #print('Animal\tOriginal\tSire\tDam\tSex\tGen')
        #for _p in example.pedigree:
            #print('%s\t%s\t\t%s\t%s\t%s\t%s' % ( _p.animalID, _p.originalID, _p.sireID, _p.damID, _p.sex, _p.gen ))

        #print(example.idmap)

        pg = pyp_network.ped_to_graph(example)

        pgd = networkx.is_directed_acyclic_graph(pg)
        print('is_directed_acyclic_graph: %s' % ( pgd ))

        pgu = pg.to_undirected()
        pgc = networkx.is_connected(pgu)
        print('is_connected: %s' % ( pgc ))

        pgncc = networkx.number_connected_components(pgu)
        print('number_connected_components: %s' % ( pgncc ))

        ped_deg = pyp_network.get_node_degrees(pg)
        ped_hist = pyp_network.get_node_degree_histograms(ped_deg)
        print('ped_hist:\t')
        for k,v in ped_hist.iteritems():
            print('\t\t', k, '\t', v)

        print('nodes:\t\t', pg.order())
        print('edges:\t\t', pg.size())
        density = pyp_network.graph_density(pg)
        print('density:\t', density)

        census = pyp_network.dyad_census(pg)
        print('max dyads:\t', ( ( pg.order()*(pg.order()-1) ) / 2 ))
        print('census:\t\t', census)

        geodesic = pyp_network.mean_geodesic(pg)
        print('geodesic:\t', geodesic)

        centrality = pyp_network.mean_degree_centrality(pg)
        print('centrality:\t', centrality)

        centrality = pyp_network.mean_degree_centrality(pg,normalize=1)
        print('normed central:\t', centrality)

        close = pyp_network.get_closeness_centrality(pg)
        mean_close = pyp_network.mean_value(close)
        print('mean closeness:\t', mean_close)

        clust = pyp_network.get_clustering_coefficient(pg)
        mean_clust = pyp_network.mean_value(clust)
        print('mean clustering:\t', mean_clust)

        between = pyp_network.get_betweenness_centrality(pg)
        mean_between = pyp_network.mean_value(between)
        print('mean betweenness:\t', mean_between)

        between = pyp_network.get_node_betweenness(pg)
        mean_between = pyp_network.mean_value(between)
        print('mean betweenness:\t', mean_between)
