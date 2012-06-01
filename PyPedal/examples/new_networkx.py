#!/usr/bin/python

###############################################################################
# NAME: generations.py
# VERSION: 2.0.0b14 (15MAY2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_newclasses
from PyPedal import pyp_network
from PyPedal.pyp_utils import pyp_nice_time

try:
    import networkx
except ImportError:
    logging.error('The networkx module could not be imported in pyp_network.  Routines using networkx functionality are not available.')

if __name__=='__main__':

    print 'Starting pypedal.py at %s' % ( pyp_nice_time() )
    print '\tLoading pedigree at %s' % ( pyp_nice_time() )

    example = pyp_newclasses.loadPedigree(optionsfile='new_networkx.ini')

    print 'Calling pyp_network.ped_to_graph()'
    ng = pyp_network.ped_to_graph(example)

#    print 'The graph has %d nodes' % len(ng.nodes())
#    print ng.nodes()
#    print ng.edges()
#
##     print 'Drawing graph madeup.ps'
##     networkx.drawing.draw(ng)
##     networkx.drawing.savefig("madeup.ps")
##
##     print 'Drawing graph madeup.dot'
##     networkx.drawing.write_dot(ng,path="madeup.dot")
##     networkx.drawing.draw_nxpydot(ng)
#
#    print 'Number of animals in pedigree: %s' % ( ng.order() )
#    print ng.nodes()
#    print 'Number of edges in pedigree: %s' % ( ng.size() )
#    print ng.edges()
#    print ng.in_edges(7)
#
#    print 'Ancestors of 7: %s' % ( ng.predecessors(7) )
#    print 'Descendants of 7: %s' % ( ng.successors(7) )
#
#    print 'Edges in graph: %s' % ( ng.edges() )
#
#    print 'Neighbors of animal 7: %s' % ( networkx.neighbors(ng,7) )
#
#    print 'All ancestors of animal 2: %s' % ( pyp_network.find_ancestors(ng,2) )
#    print 'All ancestors of animal 7: %s' % ( pyp_network.find_ancestors(ng,7) )
#    print 'All ancestors of animal 13: %s' % ( pyp_network.find_ancestors(ng,13) )
#
#    print 'All descendants of animal 2: %s' % ( pyp_network.find_descendants(ng,2,[]) )
#    print 'All descendants of animal 7: %s' % ( pyp_network.find_descendants(ng,7,[]) )
#    print 'All descendants of animal 13: %s' % ( pyp_network.find_descendants(ng,13,[]) )
#
#    print 'Immediate family of animal 2: %s' % ( pyp_network.immediate_family(ng,2) )
#    print 'Immediate family of animal 7: %s' % ( pyp_network.immediate_family(ng,7) )
#    print 'Immediate family of animal 13: %s' % ( pyp_network.immediate_family(ng,13) )
#
#    print 'Number of offspring of animal 2: %s' % ( pyp_network.count_offspring(ng,2) )
#    print 'Number of offspring of animal 7: %s' % ( pyp_network.count_offspring(ng,7) )
#    print 'Number of offspring of animal 13: %s' % ( pyp_network.count_offspring(ng,13) )
#
#    print 'Most influential offspring of animal 2: %s' % ( pyp_network.most_influential_offspring(ng,2) )
#    print 'Most influential offspring of animal 7: %s' % ( pyp_network.most_influential_offspring(ng,7) )
#    print 'Most influential offspring of animal 13: %s' % ( pyp_network.most_influential_offspring(ng,13) )
#
#    pyp_network.get_founder_descendants(ng)
#
#    print 'Degrees: %s' % ( networkx.degree(ng) )
#
#    print '='*70
#    print 'Creating a PyPedal pedigree from the graph \'ng\''
#    print '-'*70
#
#    print 'Printout of graph'
#    print ng.nodes(data=True)
#    print '-'*70
#
    options = {}
    options['pedfile'] = 'dummy'
    options['messages'] = 'verbose'
    options['renumber'] = 1
    options['pedname'] = 'Testing fromgraph()'
    options['pedformat'] = 'asd'
    options['set_offspring'] = 1
    options['set_ancestors'] = 1
    options['set_sexes'] = 1
    options['set_generations'] = 1
    example2 = pyp_newclasses.loadPedigree(options,pedsource='graph',pedgraph=ng, debugLoad=True)
    example2.metadata.printme()
#
#    print 'Stopping pypedal.py at %s' % ( pyp_nice_time() )
