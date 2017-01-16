###############################################################################
# NAME: pyp_jbc.py
# VERSION: 2.0.0 (29SEPTEMBER2010)
# AUTHOR: John B. Cole, PhD (john.cole@ars.usda.gov)
# LICENSE: LGPL`
###############################################################################
# FUNCTIONS:
#     get_color_32()
#     color_pedigree()
#     draw_colored_pedigree()
#     new_draw_colored_pedigree()
###############################################################################

## @package pyp_jbc
# pyp_template provides a skeleton on which user-defined modules may be built.
##
from __future__ import print_function
import logging, math, numpy
from PyPedal import  pyp_graphics
from PyPedal import pyp_network
from PyPedal import pyp_utils

##
# get_color_32() Converts a float value to one of a continuous range of colors
# using recipe 9.10 from the Python Cookbook.  Returns 32-bit colors rather than
# 16-bit colors.
# @param a Float value to convert to a color.
# @param cmin Minimum value in array (?).
# @param cmax Maximum value in array (?).
# @retval An RGB triplet.
def get_color_32(a, cmin, cmax):
    """
    Convert a float value to one of a continuous range of colors.
    Rewritten to use recipe 9.10 from the O'Reilly Python Cookbook.
    Returns 32-bit colors rather than 16-bit colors.
    """
    try:
        a = float( a - cmin ) / ( cmax - cmin )
    except ZeroDivisionError:
        a = 0.5 # cmax == cmin
    blue = min((max((4*(0.75-a),0.)),1.))
    red = min((max((4*(a-0.25),0.)),1.))
    green = min((max((4*math.fabs(a-0.5)-1.,0)),1.))
    _r = '%2x' % int(255*red)
    if _r[0] == ' ':
        _r = '0%s' % _r[1]
    _g = '%2x' % int(255*green)
    if _g[0] == ' ':
        _g = '0%s' % _g[1]
    _b = '%2x' % int(255*blue)
    if _b[0] == ' ':
        _b = '0%s' % _b[1]
    _triple = '#%s%s%s' % (_r,_g,_b)
    return _triple

##
# color_pedigree() forms a graph object from a pedigree object and determines the
# proportion of animals in a pedigree that are descendants of each animal in the
# pedigree. The results will be used to feed draw_colored_pedigree().
# @param pedobj A PyPedal pedigree object.
# @param metric Count for coloring nodes in the pedigree (descendants|sons).
# @param places Number of decimal palces to carry when rounding proportions.
# @param drawer Routine to use for drawing the pedigree (new|old).
# @param **kw Optional keywork arguments passed through to other functions.
# @retval A 1 for success and a 0 for failure.
def color_pedigree(pedobj, metric='descendants', places=2, drawer='new', **kw):
    """
    color_pedigree() forms a graph object from a pedigree object and determines the
    proportion of animals in a pedigree that are descendants of each animal in the
    pedigree.  The results will be used to feed draw_colored_pedigree().
    """
    _dprop = {}
    if metric == 'descendants':
        _pedgraph = pyp_network.ped_to_graph(pedobj)
        # Walk the pedigree and compute proportion of animals in the pedigree that are
        # descended from each animal.
        for _p in pedobj.pedigree:
            _dcount = pyp_network.find_descendants(_pedgraph,_p.animalID,[])
            if len(_dcount) < 1:
                _dprop[_p.animalID] = 0.0
            else:
                _dprop[_p.animalID] = round(float(len(_dcount)) / float(pedobj.metadata.num_records), places)
        del(_pedgraph)
    elif metric == 'sons':
        for _p in pedobj.pedigree:
            _dprop[_p.animalID] = float(len(_p.sons))
        #print(_dprop)
    else:
      return 0
    if drawer == 'new':
        new_draw_colored_pedigree(pedobj, _dprop, **kw)
    else:
        draw_colored_pedigree(pedobj, _dprop, **kw)

##
# draw_colored_pedigree() uses the pydot bindings to the graphviz library to produce a
# directed graph of your pedigree with paths of inheritance as edges and animals as
# nodes.  If there is more than one generation in the pedigree as determind by the "gen"
# attributes of the animals in the pedigree, draw_pedigree() will use subgraphs to try
# and group animals in the same generation together in the drawing.  Nodes will be colored
# based on the number of outgoing connections (number of offspring).
# @param pedobj A PyPedal pedigree object.
# @param shading A dictionary mapping animal IDs to levels that will be used to color nodes.
# @param gfilename The name of the file to which the pedigree should be drawn.
# @param gtitle The title of the graph.
# @param gformat The format in which the output file should be written  (JPG|PNG|PS).
# @param gsize The size of the graph: 'f': full-size, 'l': letter-sized page.
# @param gdot Whether or not to write the dot code for the pedigree graph to a file (can produce large files).
# @param gorient The orientation of the graph: 'p': portrait, 'l': landscape.
# @param gdirec Direction of flow from parents to offspring: 'TB': top-bottom, 'LR': left-right, 'RL': right-left.
# @param gname Flag indicating whether ID numbers (0) or names (1) should be used to label nodes.
# @param gfontsize Integer indicating the typeface size to be used in labelling nodes.
# @param garrow Flag indicating whether or not arrowheads should be drawn.
# @param gtitloc Indicates if the title be drawn or above ('t') or below ('b') the graph.
# @param gtitjust Indicates if the title should be center- ('c'), left- ('l'), or right-justified ('r').
# @param ghatch Shade a node if the value of its userField attribute is equal to ghatch.
# @param gprog Program to use to layout graph ('dot'|'neato').
# @retval A 1 for success and a 0 for failure.
def draw_colored_pedigree(pedobj, shading, gfilename='pedigree', gtitle='My_Pedigree', gformat='jpg', gsize='f', gdot='1', gorient='l', gdirec='', gname=0, gfontsize=10, garrow=1, gtitloc='b', gtitjust='c', ghatch='hatch', gprog='dot'):
    """
    draw_colored_pedigree() uses the pydot bindings to the graphviz library to produce a
    directed graph of your pedigree with paths of inheritance as edges and animals as
    nodes.  If there is more than one generation in the pedigree as determind by the "gen"
    attributes of the animals in the pedigree, draw_pedigree() will use subgraphs to try
    and group animals in the same generation together in the drawing.  Nodes will be
    colored based on the number of outgoing connections (number of offspring).
    """
    from pyp_utils import string_to_table_name
    _gtitle = string_to_table_name(gtitle)

    if gtitloc not in ['t','b']:
        gtitloc = 'b'
    if gtitjust not in ['c','l','r']:
        gtitjust = 'c'

    print('[DEBUG]: Entered draw_colored_pedigree()')

#     try:
    import pydot

    # Build a list of generations -- if we have more than on, we can use the
    # "rank=same" option in dot to get nicer output.
    gens = pedobj.metadata.unique_gen_list
    # Set some properties for the graph.
    g = pydot.Dot(label=gtitle, labelloc=gtitloc, labeljust=gtitjust, graph_name=_gtitle, type='graph', strict=False, suppress_disconnected=True, simplify=True)

    # Make sure that gfontsize has a valid value.
    try:
        gfontsize = int(gfontsize)
    except:
        gfontsize = 10
    if gfontsize < 10:
        gfontsize = 10
    gfontsize = str(gfontsize)
#     print('gfontsize = %s' % (gfontsize))
    g.set_page("8.5,11")
    g.set_size("7.5,10")
    if gorient == 'l':
        g.set_orientation("landscape")
    else:
        g.set_orientation("portrait")
    if gsize != 'l':
        g.set_ratio("auto")
    if gdirec == 'RL':
        g.set_rankdir('RL')
    elif gdirec == 'LR':
        g.set_rankdir('LR')
    else:
        pass
    g.set_center('true')
    g.set_concentrate('true')
    g.set_ordering('out')
    if gformat not in g.formats:
        gformat = 'jpg'
    # If we do not have any generations, we have to draw a less-nice graph.
    colormin = min(shading.values())
    colormax = max(shading.values())
    color_map = {}
    if len(gens) <= 1:
        animalCounter = 0
        print('\t[DEBUG]: Only one generation')
        for _m in pedobj.pedigree:
            animalCounter = animalCounter + 1
            if numpy.fmod(animalCounter,pedobj.kw['counter']) == 0:
                print('\t[DEBUG]: Records read: %s ' % ( animalCounter ))
            # Add a node for the current animal and set some properties.
            if gname:
                _node_name = _m.name
            else:
                _node_name = _m.animalID
            _an_node = pydot.Node(_node_name)
            _an_node.set_fontname('Helvetica')
            # _an_node.set_fontsize('10')
            _an_node.set_fontsize(gfontsize)
            _an_node.set_height('0.35')
            if _m.sex == 'M' or _m.sex == 'm':
                _an_node.set_shape('box')
            elif _m.sex == 'F' or _m.sex == 'f':
                _an_node.set_shape('ellipse')
            else:
                pass
            #print(_m.userField, ghatch, ( _m.userField == ghatch ))
            if _m.userField == ghatch:
                _an_node.set_style('filled,peripheries=2')
            else:
                _an_node.set_style('filled')
            try:
                _color = color_map[shading[_m.animalID]]
            except KeyError:
                _color = get_color_32(shading[_m.animalID], colormin, colormax)
                color_map[shading[_m.animalID]] = _color
                print('\t[DEBUG]: %s added to cache' % ( _color ))
            _an_node.set_fillcolor(_color)
            g.add_node(_an_node)
            # Add the edges to the parent nodes, if any.
            if _m.sireID != pedobj.kw['missing_parent']:
                if gname:
                    if garrow:
                        g.add_edge(pydot.Edge(pedobj.pedigree[int(_m.sireID)-1].name, _m.name))
                    else:
                        g.add_edge(pydot.Edge(pedobj.pedigree[int(_m.sireID)-1].name, _m.name, dir='none'))
                else:
                    if garrow:
                        g.add_edge(pydot.Edge(pedobj.pedigree[int(_m.sireID)-1].originalID,_m.originalID))
                    else:
                        g.add_edge(pydot.Edge(pedobj.pedigree[int(_m.sireID)-1].originalID,_m.originalID, dir='none'))
            if _m.damID != pedobj.kw['missing_parent']:
                if gname:
                    if garrow:
                        g.add_edge(pydot.Edge(pedobj.pedigree[int(_m.damID)-1].name, _m.name))
                    else:
                        g.add_edge(pydot.Edge(pedobj.pedigree[int(_m.damID)-1].name, _m.name, dir='none'))
                else:
                    if garrow:
                        g.add_edge(pydot.Edge(pedobj.pedigree[int(_m.damID)-1].originalID,_m.originalID))
                    else:
                        g.add_edge(pydot.Edge(pedobj.pedigree[int(_m.damID)-1].originalID,_m.originalID, dir='none'))
    # Otherwise we can draw a nice graph.
    else:
        for _g in gens:
            print('\t[DEBUG]: Looping over generations')
            _sg_anims = []
            _sg_name = 'sg%s' % (_g)
            sg = pydot.Subgraph(graph_name=_sg_name, suppress_disconnected=True, simplify=True)
            sg.set_simplify(True)
            animalCounter = 0
            for _m in pedobj.pedigree:
                animalCounter = animalCounter + 1
                if numpy.fmod(animalCounter,pedobj.kw['counter']) == 0:
                    print('\t[DEBUG]: Records read: %s ' % ( animalCounter ))
                if int(_m.gen) == int(_g):
                    _sg_anims.append(_m.animalID)
                # Add a node for the current animal and set some properties.
                if gname:
                    _node_name = _m.name
                else:
                    _node_name = _m.animalID
                _an_node = pydot.Node(_node_name)
                _an_node.set_fontname('Helvetica')
                _an_node.set_fontsize(gfontsize)
                _an_node.set_height('0.35')
                if _m.sex == 'M' or _m.sex == 'm':
                    _an_node.set_shape('box')
                elif _m.sex == 'F' or _m.sex == 'f':
                    _an_node.set_shape('ellipse')
                else:
                    pass
                if _m.userField == ghatch:
                    _an_node.set_style('filled,peripheries=2')
                else:
                    _an_node.set_style('filled')
                _color = get_color_32(shading[_m.animalID], colormin, colormax)
                _an_node.set_fillcolor(_color)
                sg.add_node(_an_node)
                # Add the edges to the parent nodes, if any.
                if _m.sireID != pedobj.kw['missing_parent']:
                    if gname:
                        if garrow:
                            sg.add_edge(pydot.Edge(pedobj.pedigree[int(_m.sireID)-1].name,_m.name))
                        else:
                            sg.add_edge(pydot.Edge(pedobj.pedigree[int(_m.sireID)-1].name,_m.name, dir='none'))
                    else:
                        if garrow:
                            sg.add_edge(pydot.Edge(pedobj.pedigree[int(_m.sireID)-1].originalID,_m.originalID))
                        else:
                            sg.add_edge(pydot.Edge(pedobj.pedigree[int(_m.sireID)-1].originalID,_m.originalID, dir='none'))
                #if int(_m.damID) != 0:
                if _m.damID != pedobj.kw['missing_parent']:
                    if gname:
                        if garrow:
                            sg.add_edge(pydot.Edge(pedobj.pedigree[int(_m.damID)-1].name,_m.name))
                        else:
                            sg.add_edge(pydot.Edge(pedobj.pedigree[int(_m.damID)-1].name,_m.name, dir='none'))
                    else:
                        if garrow:
                            sg.add_edge(pydot.Edge(pedobj.pedigree[int(_m.damID)-1].originalID,_m.originalID))
                        else:
                            sg.add_edge(pydot.Edge(pedobj.pedigree[int(_m.damID)-1].originalID,_m.originalID, dir='none'))
            if len(_sg_anims) > 0:
                _sg_list = ''
                for _a in _sg_anims:
                    if len(_sg_list) == 0:
                        _sg_list = 'same,%s' % (_a)
                    else:
                        _sg_list = '%s,%s' % (_sg_list,_a)
            sg.set_rank(_sg_list)
            g.add_subgraph(sg)
    # For large graphs it is nice to write out the .dot file so that it does not have to be recreated
    # whenever draw_pedigree is called.  Especially when I am debugging.  :-)
    if gdot:
        dfn = '%s.dot' % (gfilename)
#             try:
        g.write(dfn)
#             except:
#                 pass
    # Write the graph to an output file.
    outfile = '%s.%s' % (gfilename,gformat)
    if gprog not in ['dot','neato','none']:
        gprog = 'dot'
    if gprog != 'none':
        g.write(outfile,prog=gprog,format=gformat)
    return 1
#     except:
#         return 0

##
# new_draw_colored_pedigree() uses the pyygraphviz to produce a directed graph of your
# pedigree with paths of inheritance as edges and animals as nodes.  If there
# is more than one generation in the pedigree as determind by the "gen"
# attributes of the animals in the pedigree, draw_pedigree() will use subgraphs
# to try and group animals in the same generation together in the drawing.
# @param pedobj A PyPedal pedigree object.
# @param shading A dictionary mapping animal IDs to levels that will be used to color nodes.
# @param gfilename The name of the file to which the pedigree should be drawn
# @param gtitle The title of the graph.
# @param gformat The format in which the output file should be written  (JPG|PNG|PS).
# @param gsize The size of the graph: 'f': full-size, 'l': letter-sized page.
# @param gdot Whether or not to write the dot code for the pedigree graph to a file (can produce large files).
# @param gorient The orientation of the graph: 'p': portrait, 'l': landscape.
# @param gdirec Direction of flow from parents to offspring: 'TB': top-bottom, 'LR': left-right, 'RL': right-left.
# @param gname Flag indicating whether ID numbers (0) or names (1) should be used to label nodes.
# @param garrow Flag indicating whether or not arrowheads should be drawn.
# @param gtitloc Indicates if the title be drawn or above ('t') or below ('b') the graph.
# @param gtitjust Indicates if the title should be center- ('c'), left- ('l'), or right-justified ('r').
# @param gshowall Draws animals with no links to other ancestors in the pedigree (1) or suppresses them (0).
# @param gprog Specify which program should be used to position and render the graph.
# @param ghatch Shade a node if the value of its userField attribute is equal to ghatch.
# @retval A 1 for success and a 0 for failure.
def new_draw_colored_pedigree(pedobj, shading, gfilename='pedigree', \
    gtitle='', gformat='jpg', gsize='f', gdot=1, gorient='p', gdirec='', \
    gname=0, garrow=1, gtitloc='b', gtitjust='c', gshowall=1, gprog='dot', \
    ghatch='hatch'):
    """
    draw_pedigree() uses the pydot bindings to the graphviz library -- if they
    are available on your system -- to produce a directed graph of your pedigree
    with paths of inheritance as edges and animals as nodes.  If there is more than
    one generation in the pedigree as determind by the "gen" attributes of the animals
    in the pedigree, draw_pedigree() will use subgraphs to try and group animals in the
    same generation together in the drawing.
    """

    try:
        import pygraphviz
    except ImportError:
        if pedobj.kw['messages'] == 'verbose':
            print('[ERROR]: pyp_graphics/new_draw_pedigree() was unable to import the pygraphviz module!')
        logging.error('pyp_graphics/new_draw_pedigree() was unable to import the pygraphviz module!')
        return 0

    # Maps the 0/1 flags taken by the function and maps them to Python
    # True/False for settine edge and node attributes.
    _tf = {0:False, 1:True}

    from pyp_utils import string_to_table_name
    _gtitle = string_to_table_name(gtitle)

    if gtitloc not in ['t','b']:
        gtitloc = 'b'
    if gtitjust not in ['c','l','r']:
        gtitjust = 'c'

    if not pedobj.kw['pedigree_is_renumbered']:
        if pedobj.kw['messages'] != 'quiet':
            print('[GRAPH]: The pedigree that you passed to pyp_graphics/draw_pedigree() is not renumbered. Because of this, there may be errors in the rendered pedigree. In order to insure that the pedigree drawing is accurate, you should renumber the pedigree before calling draw_pedigree().')
        logging.error('The pedigree that you passed to pyp_graphics/draw_pedigree() is not renumbered. Because of this, there may be errors in the rendered pedigree. In order to insure that the pedigree drawing is accurate, you should renumber the pedigree before calling draw_pedigree().')

    # Create an empty pygraphviz graph using the Agraph class.
    g = pygraphviz.AGraph(directed=True,strict=False)

    # I'm not sure if I need to have this here or not.
    g.graph_attr['type'] = 'graph'

    # Name the graph -- _gtitle has the characters removed which would cause
    # dotfile processing to break.
    g.graph_attr['name'] = _gtitle

    # Set some properties for the graph.  The label attribute is based on the gtitle.
    # In cases where an empty string, e.g. '', is provided as the gtitle dot engine
    # processing breaks.  In such cases, don't add a label.

    if gtitle != '':
        g.graph_attr['label'] = gtitle
        g.graph_attr['labelloc'] = gtitloc
        g.graph_attr['labeljust'] = gtitjust

    # Set the page paper size and writeable area.
    g.graph_attr['page'] = '8.5,11'
    g.graph_attr['size'] = '7.5,10'

    # Set the page orientation.
    if gorient == 'l':
        g.graph_attr['orientation'] = 'landscape'
    else:
        g.graph_attr['orientation'] = 'portrait'

    if gsize != 'l':
        g.graph_attr['ratio'] = 'auto'
    if gdirec == 'RL':
        g.graph_attr['rankdir'] = 'RL'
    elif gdirec == 'LR':
        g.graph_attr['rankdir'] = 'LR'
    else:
        pass

    # Set a few other graph properties.
    g.graph_attr['center'] = 'True'
    g.graph_attr['concentrate'] = 'True'
    g.graph_attr['fontsize'] = str(pedobj.kw['default_fontsize'])
    g.graph_attr['ordering'] = 'out'

    # We need this for coloring the pedigree
    colormin = min(shading.values())
    colormax = max(shading.values())
    color_map = {}

    for _m in pedobj.pedigree:
        # Add a node for the current animal and set some node properties.
        if gname:
            _node_name = _m.name
        else:
            _node_name = _m.animalID
        g.add_node(_node_name)
        n = g.get_node(_node_name)
        n.attr['shape'] = 'box'
        n.attr['fontname'] = 'Helvetica'
        n.attr['fontsize'] = str(pedobj.kw['default_fontsize'])
        n.attr['height'] = '0.35'
        #print('[DEBUG]: sex = ', _m.sex)
        if _m.sex == 'M' or _m.sex == 'm':
            n.attr['shape'] = 'box'
        elif _m.sex == 'F' or _m.sex == 'f':
            n.attr['shape'] = 'ellipse'
        else:
            n.attr['shape'] = 'octagon'
            #pass

        # Color the nodes
        if _m.userField == ghatch:
            n.attr['style'] = 'filled,peripheries=2'
        else:
            n.attr['style'] = 'filled'
            _color = get_color_32(shading[_m.animalID], colormin, colormax)
            n.attr['fillcolor'] = _color
            # Add values to the color map
            if not color_map.has_key(shading[_m.animalID]):
                color_map[shading[_m.animalID]] = _color

        # Add the edges to the parent nodes, if any.
        if _m.sireID != pedobj.kw['missing_parent']:
            if gname:
                _sire_edge = pedobj.pedigree[int(_m.sireID)-1].name
            else:
                # Check some outputs -- should I be using the animalID or the
                # originalID to assign edges? Nodes are based on the animalID,
                # so edges should also be in order to maintain consistency.
                #_sire_edge = pedobj.pedigree[int(_m.sireID)-1].originalID
                _sire_edge = pedobj.pedigree[int(_m.sireID)-1].animalID
            g.add_edge(_sire_edge,_node_name)
            if not _tf[garrow]:
                e = g.get_edge(_sire_edge,_anim_node)
                e.attr['dir'] = 'none'
        if _m.damID != pedobj.kw['missing_parent']:
            if gname:
                _dam_edge = pedobj.pedigree[int(_m.damID)-1].name
            else:
                _dam_edge = pedobj.pedigree[int(_m.damID)-1].animalID
            g.add_edge(_dam_edge,_node_name)
            if not _tf[garrow]:
                e = g.get_edge(_dam_edge,_anim_node)
                e.attr['dir'] = 'none'

    # For large graphs it is nice to write out the .dot file so that it does
    # not have to be recreated whenever new_draw_pedigree is called.
    # Especially when I am debugging.
    if gdot:
        dfn = '%s.dot' % (gfilename)
        try:
            g.write(dfn)
        except:
            if pedobj.kw['messages'] == 'verbose':
                print('[ERROR]: pyp_graphics/new_draw_pedigree() was unable to write the dotfile %s.' % (dfn))
            logging.error('pyp_graphics/new_draw_pedigree() was unable to draw the dotfile %s.', (dfn))


    # Write the color map to a file.
    #try:
    mapfile = '%s_color_map.txt' % (gfilename)
    mf = file(mapfile,'w')
    mf.write('# Color map data\n')
    mf.write('# Data are metric (number of sons/descendants/etc.) followed\n')
    mf.write('# by color in RGB.\n')
    for k,v in color_map.iteritems():
        line = '%s\t%s\n' % (k,v)
        mf.write(line)
    mf.close()
    #except:
        #outfile = '%s_color_map.txt' % (gfilename)
        #if pedobj.kw['messages'] == 'verbose':
            #print('[ERROR]: pyp_jbc/new_draw_colored_pedigree() was unable to write the color map %s.' % (outfile))
        #logging.error('pyp_jbc/new_draw_colored_pedigree() was unable to write the color map %s.', (outfile))
        #return 0

    # Write the graph to an output file.
    #try:
    outfile = '%s.%s' % (gfilename,gformat)
    g.draw(outfile,prog=gprog)
    return 1
    #except:
        #outfile = '%s.%s' % (gfilename,gformat)
        #if pedobj.kw['messages'] == 'verbose':
            #print('[ERROR]: pyp_jbc/new_draw_colored_pedigree() was unable to draw the pedigree %s.' % (outfile))
        #logging.error('pyp_jbc/new_draw_colored_pedigree() was unable to draw the pedigree %s.', (outfile))
        #return 0
