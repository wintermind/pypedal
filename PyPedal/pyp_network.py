#!/usr/bin/python

###############################################################################
# NAME: pyp_network.py
# VERSION: 2.0.0 (29SEPTEMBER2010)
# AUTHOR: John B. Cole, PhD (john.cole@ars.usda.gov)
# LICENSE: LGPL
###############################################################################
# FUNCTIONS:
#   ped_to_graph()
#   find_ancestors()
#   find_descendants()
#   immediate_family()
#   count_offspring()
#   offspring_influence()
#   most_influential_offspring()
#   get_founder_descendants()
#   ---------------------------------------------------------------------------
#   get_node_degrees()
#   get_node_degree_histograms()
#   mean_geodesic()
#   graph_density()
#   dyad_census()
#   mean_degree_centrality()
#   mean_value()
#   get_closeness_centrality()
#   get_clustering_coefficient()
#   get_betweenness_centrality()
#   get_node_betweenness()
###############################################################################

## @package pyp_network
# pyp_network contains a set of procedures for working with pedigrees as directed
# graphs.
##

import copy, logging

try:
    import networkx
except ImportError:
    print '[WARNING]: The NetworkX module could not be imported in pyp_network.  Routines using networkx functionality are not available.'

##
# ped_to_graph() Takes a PyPedal pedigree object and returns a networkx DiGraph
# object.
# @param pedobj A PyPedal pedigree object.
# @param oid Flag indicating if original (1) or renumbered (0) IDs should be used.
# @retval DiGraph object
def ped_to_graph(pedobj, oid=0):
    """
    ped_to_graph() Takes a PyPedal pedigree object and returns a networkx DiGraph
    object.
    """
    l = len(pedobj.pedigree)
    G = networkx.DiGraph(name=pedobj.kw['pedname'], selfloops=False, multiedges=True)
    for i in range(l):
        # The order in which we pass arguments to add_edge() is important -- the
        # parent has to be the first argument and the offspring the second if the
        # graph is to be ordered in the correct direction.
        #print 'Adding node %s' % ( i )
        if oid:
            G.add_node(int(pedobj.pedigree[i].originalID),
               sire=str(pedobj.pedigree[int(pedobj.pedigree[i].sireID)].originalID),
               dam=str(pedobj.pedigree[int(pedobj.pedigree[i].damID)].originalID))
        else:
            G.add_node(int(pedobj.pedigree[i].animalID),
               sire=str(pedobj.pedigree[i].sireID),
               dam=str(pedobj.pedigree[i].damID))
    # Add parent-to-child edges.
    if str(pedobj.pedigree[i].sireID) != str(pedobj.kw['missing_parent']):
        if oid:
            G.add_edge(pedobj.pedigree[int(pedobj.pedigree[i].sireID)].originalID, int(pedobj.pedigree[i].originalID), sex='s')
        else:
            G.add_edge(int(pedobj.pedigree[i].sireID), int(pedobj.pedigree[i].animalID), sex='s')
    if str(pedobj.pedigree[i].damID) != str(pedobj.kw['missing_parent']):
        if oid:
            G.add_edge(pedobj.pedigree[int(pedobj.pedigree[i].damID)].originalID, int(pedobj.pedigree[i].originalID), sex='d')
        else:
            G.add_edge(int(pedobj.pedigree[i].damID), int(pedobj.pedigree[i].animalID), sex='d')
    return G

##
# find_ancestors() identifies the ancestors of an animal and returns them in a list.
# @param pedgraph An instance of a networkx DiGraph.
# @param anid The animal for whom ancestors are to be found.
# @param _ancestors The list of ancestors already found.
# @retval List of ancestors of anid.
def find_ancestors(pedgraph, anid,_ancestors=[]):
    """
    find_ancestors() identifies the ancestors of an animal and returns them in a list.
    """
#    anid = int(anid)
#    #print 'anid:\t\t%d' % ( anid )
#    try:
#        _pred = pedgraph.predecessors(anid)
#        for _p in _pred:
#            if _p not in _ancestors:
#                _ancestors.append(int(_p))
#            #print _ancestors
#            find_ancestors(pedgraph,_p,_ancestors)
#    except:
#        pass
    _ancestors = list(networkx.algorithms.dag.ancestors(pedgraph, anid))
    return _ancestors

##
# find_ancestors_g() identifies the ancestors of an animal going back a user-specified number of generations and returns them in a list.
# @param pedgraph An instance of a networkx DiGraph.
# @param anid The animal for whom ancestors are to be found.
# @param _ancestors The list of ancestors already found.
# @param gens Number of generations to go back in the pedigree.
# @retval List of ancestors of anid.
def find_ancestors_g(pedgraph, anid, _ancestors={}, gens=3):
    """
    find_ancestors_g() identifies the ancestors of an animal and returns them in a list.
    """
    #print 'anid:\t\t%d' % ( anid )
    #print 'gens:\t\t%d' % ( gens )
    if gens == 0:
        return _ancestors
    else:
        try:
            _pred = pedgraph.predecessors(anid)
            #print '_pred:\t\t%s' % ( _pred )
            for _p in _pred:
                if _p not in _ancestors:
                    _ancestors[_p] = gens
                    find_ancestors_g(pedgraph, _p, _ancestors, gens-1)
        except:
            pass
    #print 'ancestors:\t%s' % ( _ancestors )
    return _ancestors

##
# find_descendants() identifies the descendants of an animal and returns them in a list.
# @param pedgraph An instance of a networkx DiGraph.
# @param anid The animal for whom descendants are to be found.
# @param _descendants The list of descendants already found.
# @retval List of descendants of anid.
def find_descendants(pedgraph, anid, _descendants=[]):
    """
    find_descendants() identifies the descendants of an animal and returns them in a list.
    """
#    try:
#        #print '\t%s\t%s' % ( anid, _descendants )
#        _desc = pedgraph.successors(anid)
#        for _d in _desc:
#            if _d not in _descendants:
#                _descendants.append(_d)
#                find_descendants(pedgraph,_d,_descendants)
#    except:
#        pass
    _descendants = list(networkx.algorithms.dag.descendants(pedgraph, anid))
    return _descendants

##
# immediate_family() returns parents and offspring of an animal.
# @param pedgraph An instance of a networkx DiGraph.
# @param anid The animal for whom immediate family are to be found.
# @retval List of immediate family members of anid.
def immediate_family(pedgraph, anid):
    """
    immediate_family() returns parents and offspring of an animal.
    """
    try:
        _family = networkx.neighbors(pedgraph,anid)
    except:
        pass
    return _family

##
# count_offspring() returns the number of offspring of an animal.
# @param pedgraph An instance of a networkx DiGraph.
# @param anid The animal for whom offspring are to be counted.
# @retval Count of offspring.
def count_offspring(pedgraph, anid):
    """
    count_offspring() returns the number of offspring of an animal.
    """
    try:
      _count = len(networkx.neighbors(pedgraph,anid)) - len(pedgraph.predecessors(anid))
    except:
        _count = 0
    return _count

##
# offspring_influence() returns the number of grand-children by each child of a given animal.
# @param pedgraph An instance of a networkx DiGraph.
# @param anid The animal for whom grand-progeny are to be counted.
# @retval A dictionary of counts of progeny per child.
def offspring_influence(pedgraph, anid):
    """
    offspring_influence() returns the number of grand-children by each child of a given
    animal.
    """
    try:
      _offspring = pedgraph.successors(anid)
      _degrees = networkx.degree(pedgraph, nbunch=_offspring, with_labels=True)
      # Correct for edges contributed by parents
      for _d in _degrees:
          _degrees[_d] = _degrees[_d] - len(pedgraph.predecessors(_d))
    except:
        pass
    return _degrees

##
# most_influential_offspring() returns the most influential offspring of an animal as measured by their number of offspring.
# @param pedgraph An instance of a networkx DiGraph.
# @param anid The animal for whom the most influential offspring is to be found.
# @param resolve Indicates how ties should be handled ('first'|'last'|'all').
# @retval The most influential offspring of anid.
def most_influential_offspring(pedgraph, anid, resolve='all'):
    """
    most_influential_offspring() returns the most influential offspring of an animal as
    measured by their number of offspring.
    """
    try:
        _max_off = -999
        _offid = -999
        _offdict = {}
        if resolve == 'all':
            _tempoffdict = {}
        _offspring = offspring_influence(pedgraph,anid)
        #print _offspring
        for _o in _offspring:
            # If there are ties, return the ID for the first offspring
            # seen with that number of progeny.
            if resolve == 'first':
                if _offspring[_o] > _max_off:
                    _max_off = _offspring[_o]
                    _offid = _o
            # If there are ties, return the ID for the last offspring
            # seen with that number of progeny.
            elif resolve == 'last':
                if _offspring[_o] >= _max_off:
                    _max_off = _offspring[_o]
                    _offid = _o
            # If there are ties, return the IDs for all offspring
            # with the highest number of progeny.
            else:
                if  _offspring[_o] > _max_off:
                    _tempoffdict = {}
                    _tempoffdict[_o] = _offspring[_o]
                    _max_off = _offspring[_o]
                    _offid = _o
                elif _offspring[_o] == _max_off:
                    _tempoffdict[_o] = _offspring[_o]
                    _offid = _o
        if resolve == 'all':
            _offdict = _tempoffdict
        else:
            _offdict[_offid] = _max_off
    except:
        pass
    return _offdict

##
# get_founder_descendants() returns a dictionary containing a list of descendants of
# each founder in the pedigree.
# @param pedgraph An instance of a NetworkX DiGraph.
# @retval A dictionary containing a list of descendants for each founder in the graph.
def get_founder_descendants(pedgraph):
    """
    get_founder_descendants() returns a dictionary containing a list of descendants of
    each founder in the pedigree.
    """
    try: logging.info('Entered get_founder_descendants()')
    except: pass
    founder_desc = {}
    for _n in pedgraph.nodes_iter():
        if pedgraph.in_degree(_n) < 2:
            _desc = {}
            _desc = find_descendants(pedgraph,_n,[])
            founder_desc[_n] = _desc
            print '\t%s\t%s' % ( _n, founder_desc[_n] )
    try: logging.info('Exited get_founder_descendants()')
    except: pass
    return founder_desc

##
# get_node_degrees() returns a dictionary containing the number of
# vertices (nodes) in pg with a given number of incoming, outgoing,
# or total edges.
# @param pg An instance of a NetworkX DiGraph.
# @retval A dictionary of dictionaries containing in, out, and total node degrees.
def get_node_degrees(pg):
    """
    get_node_degrees() returns a dictionary containing the number of
    vertices (nodes) in pg with a given number of incoming, outgoing,
    or total edges.
    """
    try:
        node_degrees = {}
        node_degrees['in'] = {}
        node_degrees['out'] = {}
        node_degrees['all'] = {}
        for n in pg.nodes_iter():
            node_degrees['in'][n] = pg.in_degree(n)
            node_degrees['out'][n] = pg.out_degree(n)
            node_degrees['all'][n] = pg.degree(n)
        return node_degrees
    except:
        return 0

##
# get_node_degree_histograms() returns a dictionary containing histograms of
# the number of vertices (nodes) in pg with a given number of incoming,
# outgoing, or total edges.
# @param node_degrees A dictionary of dictionaries of in, out, and total node degrees.
# @retval A dictionary of dictionaries containing in, out, and total degree histograms.
def get_node_degree_histograms(node_degrees):
    """
    get_node_degree_histograms() returns a dictionary containing histograms of
    the number of vertices (nodes) in pg with a given number of incoming,
    outgoing, or total edges.
    """
    try:
        histogram = {}
        histogram['in'] = {}
        histogram['out'] = {}
        histogram['all'] = {}
        for i in range(0,max(node_degrees['in'].values())+1): histogram['in'][i] = 0
        for i in range(0,max(node_degrees['out'].values())+1): histogram['out'][i] = 0
        for i in range(0,max(node_degrees['all'].values())+1): histogram['all'][i] = 0

        for k in node_degrees['in'].keys():
            histogram['in'][node_degrees['in'][k]] = histogram['in'][node_degrees['in'][k]] + 1

        for k in node_degrees['out'].keys():
            histogram['out'][node_degrees['out'][k]] = histogram['out'][node_degrees['out'][k]] + 1

        for k in node_degrees['all'].keys():
            histogram['all'][node_degrees['all'][k]] =         histogram['all'][node_degrees['all'][k]] + 1
        return histogram
    except:
        return 0

##
# mean_geodesic() calculates the mean geodesic (shortest) distance
# between two vertices in a network.
# @param pg An instance of a NetworkX DiGraph.
# @param debug Flag to turn debugging messages on (1) or off (0).
# @retval The mean geodesic for the pedigree represented in pg.
def mean_geodesic(pg, debug=0):
    """
    mean_geodesic() calculates the mean geodesic (shortest) distance
    between two vertices in a network.
    """
    length_sum = 0
    if networkx.is_directed_acyclic_graph(pg):
        n_pairs_with_paths = 0
    else:
        n_pairs_with_paths = ( pg.order() * ( pg.order() + 1 ) ) / 2
    tg = networkx.subgraph(pg, pg.nodes())
    for u in pg.nodes_iter():
        tg.delete_node(u)
        for v in tg.nodes_iter():
            try:
                length = networkx.shortest_path_length(pg,u,v)
                if length > 0:
                    length_sum = length_sum + length
                    if networkx.is_directed_acyclic_graph(pg):
                        n_pairs_with_paths = n_pairs_with_paths + 1
            except networkx.exception.NetworkXError:
                pass
    try:
        geodesic = float(length_sum) / float(n_pairs_with_paths)
    except:
        geodesic = -999.
    if debug:
        print 'length_sum:\t', length_sum
        print 'n_pairs_with_paths:\t', n_pairs_with_paths
    return geodesic

##
# graph_density() calculates the density of a digraph, which is the
# ratio of edges in the gaph to the maximum possible number of edges.
# @param pg An instance of a NetworkX DiGraph.
# @retval The density of the pedigree represented in pg.
def graph_density(pg):
    """
    graph_density() calculates the density of a digraph, which is the
    ratio of edges in the gaph to the maximum possible number of edges.
    """
    pgn = pg.order()    # Number of nodes in pg
    pgne = pg.size()    # Number of edges in pg
    denom = pgn * ( pgn - 1 )

    try:
        density = float(pgne) / float(denom)
    except:
        return -999.
    return density

##
# dyad_census() calculates the number of null, asymmetric, and
# mutual edges between all pairs of nodes in a directed graph.
# @param pg An instance of a NetworkX DiGraph.
# @param debug Flag to turn debugging messages on (1) or off (0).
# @param debuglog Flag to turn debugging messages to the logfile  on (1) or off (0).
# @retval A dictionary of counts.
def dyad_census(pg, debug=0, debuglog=0):
    """
    dyad_census() calculates the number of null, asymmetric, and
    mutual edges between all pairs of nodes in a directed graph.
    """
    if not networkx.is_directed_acyclic_graph(pg):
        logging.error('pyp_network.dyad_census() requires a directed graph as input!')
        return 0
    else:
        census = {}
        census['null'] = 0
        census['asymmetric'] = 0
        census['mutual'] = 0
        tg = networkx.subgraph(pg, pg.nodes())
        for u in pg.nodes_iter():
            tg.delete_node(u)
            for v in tg.nodes_iter():
                if not pg.has_neighbor(u,v):
                    census['null'] = census['null'] + 1
                elif u in pg.predecessors(v) and v in pg.successors(u):
                    census['mutual'] = census['mutual'] + 1
                    if debug:
                        print 'Nodes %s and %s link to one another!' % ( u, v )
                    if debuglog:
                        logging.error('Nodes %s and %s link to one another!',u, v)
                elif u in pg.predecessors(v) and v not in pg.successors(u):
                    census['asymmetric'] = census['asymmetric'] + 1
                elif u not in pg.predecessors(v) and v in pg.successors(u):
                    census['asymmetric'] = census['asymmetric'] + 1
                else:
                    pass
        del(tg)
        return census

##
# mean_degree_centrality(pg) calculates mean in- and out-degree
# centralities for directed graphs and simple degree-centralities
# for undirected graphs. If the normalize flag is set, each node's
# centralities are weighted by the number of edges in the (di)graph.
# @param pg An instance of a NetworkX DiGraph.
# @param normalize Flag to turn normalization on (1) or off (0).
# @retval A dictionary of mean degree centralities.
def mean_degree_centrality(pg, normalize=0):
    """
    mean_degree_centrality(pg) calculates mean in- and out-degree
    centralities for directed graphs and simple degree-centralities
    for undirected graphs. If the normalize flag is set, each node's
    centralities are weighted by the number of edges in the (di)graph.
    """
    centrality = {}
    try:
        if networkx.is_directed_acyclic_graph(pg):
            cent_sum_in, cent_sum_out = 0, 0
            for n in pg.nodes():
                n_cent_in = pg.in_degree(n)
                n_cent_out = pg.out_degree(n)
                if normalize:
                    n_cent_in = float(n_cent_in) / float(pg.size()-1)
                    n_cent_out = float(n_cent_out) / float(pg.size()-1)
                cent_sum_in = cent_sum_in + n_cent_in
                cent_sum_out = cent_sum_out + n_cent_out
            centrality['in'] = cent_sum_in / float(pg.order())
            centrality['out'] = cent_sum_out / float(pg.order())
        else:
            cent_sum = 0
            for n in pg.nodes():
                if not normalize:
                    n_cent = pg.degree(n)
                else:
                    n_cent = networkx.degree_centrality(pg,n)
                cent_sum = cent_sum + n_cent
            centrality['all'] = cent_sum / float(pg.order())
    except:
        logging.error('pyp_network.mean_degree_centrality() failed!')
    return centrality

##
# mean_value() calculates the mean from all values in a dictionary.
# @param mydict A dictionary whose values are integers or reals.
# @retval The mean of values in mydict.
def mean_value(mydict):
    """
    mean_value() calculates the mean from all values in a
    dictionary.
    """
    try:
        mymean = float(sum(mydict.values())) / float(len(mydict.values()))
    except:
        logging.error('pyp_network.mean_values() failed!')
        mymean =  -999.
    return mymean

##
# get_closeness_centrality(pg) returns a dictionary of the closeness
# centrality (1/(average distance to all nodes from n)) for each
# node in the graph.
# @param pg An instance of a NetworkX DiGraph.
# @retval A dictionary of closeness centralities.
def get_closeness_centrality(pg):
    """
    get_closeness_centrality(pg) returns a dictionary of the close-
    ness centrality (1/(average distance to all nodes from n)) for
    each node in the graph.
    """
    centrality = {}
    try:
        for n in pg.nodes():
            centrality[n] = networkx.closeness_centrality(pg,n)
    except:
        logging.error('pyp_network.get_closeness_centrality() failed!')
    return centrality

##
# get_clustering_coefficient(pg) returns a dictionary of the clustering
# coefficient of each node in the graph based on the definition of Watts
# and Strogatz (1998).
# @param pg An instance of a NetworkX DiGraph.
# @retval A dictionary of clustering coefficients.
def get_clustering_coefficient(pg):
    """
    get_clustering_coefficient(pg) returns a dictionary of the clustering
    coefficient of each node in the graph based on the definition of Watts
    and Strogatz (1998).
    """
    clustering = {}
    try:
        for n in pg.nodes():
            clustering[n] = networkx.clustering(pg,n)
    except:
        logging.error('pyp_network.get_clustering() failed!')
    return clustering

##
# get_betweenness_centrality(pg) returns a dictionary of the
# betweenness centrality of each node in the graph.
# @param pg An instance of a NetworkX DiGraph.
# @retval A dictionary of closeness centralities.
def get_betweenness_centrality(pg):
    """
    get_betweenness_centrality(pg) returns a dictionary of the
    betweenness centrality of each node in the graph.
    """
    between = {}
    #try:
    between = networkx.betweenness_centrality(pg)
    #except:
        #logging.error('pyp_network.get_betweenness_centrality() failed!')
    return between

##
# get_node_betweenness(pg) returns a dictionary of the
# number of shortest paths that go through each node in
# the graph.
# @param pg An instance of a NetworkX DiGraph.
# @retval A dictionary of node betweennesses.
def get_node_betweenness(pg):
    """
    get_node_betweenness(pg) returns a dictionary of the
    number of shortest paths that go through each node in
    the graph.
    """
    between = {}
    #try:
    for n in pg.nodes():
        #between[n] = networkx.node_betweenness(pg,n)
        between[n] = networkx.betweenness_centrality(pg,n)
    #except:
        #logging.error('pyp_network.get_node_betweenness() failed!')
    return between

#eccentricity(G,n)             - maximum of shortest-path lengths from n to anywhere in G.
#triangles(G,n)                - number of triangles which include n.
