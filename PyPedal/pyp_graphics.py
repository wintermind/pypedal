#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# NAME: pyp_graphics.py
# VERSION: 2.0.0 (29SEPTEMBER2010)
# AUTHOR: John B. Cole, PhD (john.cole@ars.usda.gov)
# LICENSE: LGPL
###############################################################################
# FUNCTIONS:
#   rmuller_spy_matrix_pil() [1]
#   rmuller_pcolor_matrix_pil() [1]
#   rmuller_get_color() [1]
#   draw_pedigree()
#   plot_founders_by_year()
#   plot_founders_pct_by_year()
#   plot_line_xy()
#   pcolor_matrix_pylab()
#   spy_matrix_pylab()
#   new_draw_pedigree()
###############################################################################
# [1] These routines were taken from the ASPN Python Cookbook
#     (http://aspn.activestate.com/ASPN/Cookbook/Python/) and are used under
#     terms of the license, "Python Cookbook code is freely available for use
#     and review" as I understand it.  I did not write them; Rick Muller
#     properly deserves credit for that.  I THINK Rick's homepage's is:
#     http://www.cs.sandia.gov/~rmuller/.
# Python Cookbook notes:
#     Title: Matlab-like 'spy' and 'pcolor' functions
#     Submitter: Rick Muller (other recipes)
#     Last Updated: 2005/03/02
#     Version no: 1.0
#     Category: Image
#     Description:
#         I really like the 'spy' and 'pcolor' functions, which are useful in viewing
#         matrices. 'spy' prints colored blocks for values that are above a threshold,
#         and 'pcolor' prints out each element in a continuous range of colors.  The
#         attached is a little Python/PIL script that does these functions for Numpy
#         arrays.
#     URL: http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/390208
###############################################################################

## @package pyp_graphics
# pyp_graphics contains routines for working with graphics in PyPedal, such as
# creating directed graphs from pedigrees using PyDot and visualizing relationship
# matrices using Rick Muller's spy and pcolor routines
# (http://aspn.activestate.com/ASPN/Cookbook/Python/).
#
# The Python Imaging Library (http://www.pythonware.com/products/pil/),
# matplotlib (http://matplotlib.sourceforge.net/), Graphviz (http://www.graphviz.org/),
# and pydot (http://dkbza.org/pydot.html) are required by one or more routines in this
# module.  They ARE NOT distributed with PyPedal and must be installed by the end-user!
# Note that the matplotlib functionality in PyPedal requires only the Agg backend, which
# means that you do not have to install GTK/PyGTK or WxWidgets/PyWxWidgets just to use
# PyPedal. Please consult the sites above for licensing and installation information.

##
# rmuller_spy_matrix_pil() implements a matlab-like 'spy' function to display the
# sparsity of a matrix using the Python Imaging Library.
# @param A Input Numpy matrix (such as a numerator relationship matrix).
# @param fname Output filename to which to dump the graphics (default 'tmp.png')
# @param cutoff Threshold value for printing an element (default 0.1)
# @param do_outline Whether or not to print an outline around the block (default 0)
# @param height The height of the image (default 300)
# @param width The width of the image (default 300)
# @retval None
def rmuller_spy_matrix_pil(A,fname='tmp.png',cutoff=0.1,do_outline=0, height=300, width=300):
    """
    rmuller_spy_matrix_pil() implements a matlab-like 'spy' function to display the
    sparsity of a matrix using the Python Imaging Library.
    """
    try:
        import Image, ImageDraw
    except:
        return 0
    img = Image.new("RGB",(width,height),(255,255,255))
    draw = ImageDraw.Draw(img)
    n,m = A.shape
    if n > width or m > height:
        raise "Rectangle too big %d %d %d %d" % (n,m,width,height)
    for i in range(n):
        xmin = width*i/float(n)
        xmax = width*(i+1)/float(n)
        # If the matrix is square, we can save a few cycles by only
        # doing the calculations for the upper triangle and reflecting
        # around the diagonal.
        if n == m:
            for j in range(i,n-1):
#                 print 'i: %s, j: %s' % ( i, j )
                ymin = height*j/float(m)
                ymax = height*(j+1)/float(m)
                if abs(A[i,j]) > cutoff:
                    if do_outline:
                        draw.rectangle((xmin,ymin,xmax,ymax),fill=(0,0,255),
                            outline=(0,0,0))
                        draw.rectangle((ymin,xmin,ymax,xmax),fill=(0,0,255),
                            outline=(0,0,0))
                    else:
                        draw.rectangle((xmin,ymin,xmax,ymax),fill=(0,0,255))
                        draw.rectangle((ymin,xmin,ymax,xmax),fill=(0,0,255))
        # If the matrix is not square, we can't shave a bit
        else:
            for j in range(m):
                ymin = height*j/float(m)
                ymax = height*(j+1)/float(m)
                if abs(A[i,j]) > cutoff:
                    if do_outline:
                        draw.rectangle((xmin,ymin,xmax,ymax),fill=(0,0,255),
                            outline=(0,0,0))
                    else:
                        draw.rectangle((xmin,ymin,xmax,ymax),fill=(0,0,255))
    img.save(fname)
    return

##
# rmuller_pcolor_matrix_pil() implements a matlab-like 'pcolor' function to
# display the large elements of a matrix in pseudocolor using the Python Imaging
# Library.
# @param A Input Numpy matrix (such as a numerator relationship matrix).
# @param fname Output filename to which to dump the graphics (default 'tmp.png')
# @param do_outline Whether or not to print an outline around the block (default 0)
# @param height The height of the image (default 300)
# @param width The width of the image (default 300)
# @retval A list of Animal() objects; a pedigree metadata object.
def rmuller_pcolor_matrix_pil(A, fname='tmp.png', do_outline=0, height=300, width=300):
    """
    rmuller_pcolor_matrix_pil() implements a matlab-like 'pcolor' function to
    display the large elements of a matrix in pseudocolor using the Python Imaging
    Library.
    """
    try:
        import Image, ImageDraw
    except:
        return 0
    key_dict = {}
    color_cache = {}

    img = Image.new("RGB",(width,height),(255,255,255))
    draw = ImageDraw.Draw(img)

    # For Numeric
    #mina = min(min(A))
    #maxa = max(max(A))

    # For NumPy
    mina = A.min()
    maxa = A.max()

    n,m = A.shape
    if n > width or m > height:
        raise "Rectangle too big %d %d %d %d" % (n,m,width,height)
    for i in range(n):
        xmin = width*i/float(n)
        xmax = width*(i+1)/float(n)
        for j in range(m):
            ymin = height*j/float(m)
            ymax = height*(j+1)/float(m)
            # JBC added a dictionary to cache colors to reduce the number of calls
            # to rmuller_get_color().  This may lead to a dramatic improvement in the
            # performance of rmuller_pcolor_matrix_pil(), which currently makes a call
            # to rmuller_get_color() for each of the n**2 elements of A.  The cache will
            # reduce that to the number of unique values in A, which should be much
            # smaller.
            _cache_key = '%s_%s_%s' % (A[i,j],mina,maxa)
            try:
                color = color_cache[_cache_key]
            except KeyError:
                color = rmuller_get_color(A[i,j],mina,maxa)
                color_cache[_cache_key] = color
            # JBC added this to generate a color key.
            try:
                _v = key_dict[color]
            except KeyError:
                #_key = round( (A[i,j] * 1000), 0 )
                key_dict[A[i,j]] = color
            if do_outline:
                draw.rectangle((xmin,ymin,xmax,ymax),fill=color,outline=(0,0,0))
            else:
                draw.rectangle((xmin,ymin,xmax,ymax),fill=color)

    #print key_dict
    img.save(fname)
    return

##
# rmuller_get_color() Converts a float value to one of a continuous range of colors
# using recipe 9.10 from the Python Cookbook.
# @param a Float value to convert to a color.
# @param cmin Minimum value in array.
# @param cmax Maximum value in array.
# @retval An integer containins an RGB triplet.
def rmuller_get_color(a, cmin, cmax):
    """
    Convert a float value to one of a continuous range of colors.
    Rewritten to use recipe 9.10 from the O'Reilly Python Cookbook.
    """
    try:
        a = float(a-cmin)/(cmax-cmin)
    except ZeroDivisionError:
        a = 0.5 # cmax == cmin
    blue = min((max((4*(0.75-a),0.)),1.))
    red = min((max((4*(a-0.25),0.)),1.))
    green = min((max((4*math.fabs(a-0.5)-1.,0)),1.))
    return '#%1x%1x%1x' % (int(15*red),int(15*green),int(15*blue))

##
# draw_pedigree() uses the pydot bindings to the graphviz library -- if they
# are available on your system -- to produce a directed graph of your pedigree
# with paths of inheritance as edges and animals as nodes.  If there is more than
# one generation in the pedigree as determind by the "gen" attributes of the animals
# in the pedigree, draw_pedigree() will use subgraphs to try and group animals in the
# same generation together in the drawing.
# @param pedobj A PyPedal pedigree object.
# @param gfilename The name of the file to which the pedigree should be drawn
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
# @param gshowall Draws animals with no links to other ancestors in the pedigree (1) or suppresses them (0).
# @param gclusters Indicates if subgraph clusters should be used when drawing the pedigree, which may improve layout in some pedigrees.
# @param gwidth Width of the canvas to plot on, in inches
# @param gheight Height of the canvas to plot on, in inches
# @param gdpi Resolution of image in DPI
# @retval A 1 for success and a 0 for failure.
def draw_pedigree(pedobj, gfilename='pedigree', gtitle='', gformat='jpg', gsize='f', gdot='1', gorient='p', gdirec='', gname=0, gfontsize=10, garrow=1, gtitloc='b', gtitjust='c', gshowall=1, gclusters=0, gwidth=8.5, gheight=11, gdpi=150):
    """
    draw_pedigree() uses the pydot bindings to the graphviz library -- if they
    are available on your system -- to produce a directed graph of your pedigree
    with paths of inheritance as edges and animals as nodes.  If there is more than
    one generation in the pedigree as determind by the "gen" attributes of the animals
    in the pedigree, draw_pedigree() will use subgraphs to try and group animals in the
    same generation together in the drawing.
    """
    if gtitle == '':
        gtitle = pedobj.kw['pedname']
    from pyp_utils import string_to_table_name
    _gtitle = string_to_table_name(gtitle)
#     if pedobj.kw['messages'] == 'verbose':
#         print 'gtitle: %s' % ( gtitle )
#         print '_gtitle: %s' % ( _gtitle )

    if gtitloc not in ['t','b']:
        gtitloc = 'b'
    if gtitjust not in ['c','l','r']:
        gtitjust = 'c'

    #print pedobj.metadata.unique_gen_list

    if not pedobj.kw['pedigree_is_renumbered']:
        if pedobj.kw['messages'] != 'quiet':
            print '[GRAPH]: The pedigree that you passed to pyp_graphics/draw_pedigree() is not renumbered. Because of this, there may be errors in the rendered pedigree. In order to insure that the pedigree drawing is accurate, you should renumber the pedigree before calling draw_pedigree().'
        logging.error('The pedigree that you passed to pyp_graphics/draw_pedigree() is not renumbered. Because of this, there may be errors in the rendered pedigree. In order to insure that the pedigree drawing is accurate, you should renumber the pedigree before calling draw_pedigree().')

    #try:
    import pydot
    # Set some properties for the graph.  The label attribute is based on the gtitle.
    # In cases where an empty string, e.g. '', is provided as the gtitle dot engine
    # processing breaks.  In such cases, don't add a label.
    if gtitle == '':
        if gshowall:
            g = pydot.Dot(graph_name=str(_gtitle), graph_type='digraph', strict=False, suppress_disconnected=False, simplify=True)
        else:
            g = pydot.Dot(graph_name=str(_gtitle), graph_type='digraph', strict=False, suppress_disconnected=True, simplify=True)
    else:
        if gshowall:
            g = pydot.Dot(label=str(_gtitle), labelloc=str(gtitloc), labeljust=str(gtitjust), graph_name=str(_gtitle), graph_type='digraph', strict=False, suppress_disconnected=False, simplify=True)
        else:
            g = pydot.Dot(label=str(_gtitle), labelloc=str(gtitloc), labeljust=str(gtitjust), graph_name=str(_gtitle), graph_type='digraph', strict=False, suppress_disconnected=True, simplify=True)
    # Make sure that gfontsize has a valid value.
    try:
        gfontsize = int(gfontsize)
    except:
        gfontsize = 10
    if gfontsize < 10:
        gfontsize = 10
    gfontsize = str(gfontsize)
#     print 'gfontsize = %s' % (gfontsize)
    # Check the resolution and set if it's a positive integer, otherwise set it to the default.
    dpi = int(gdpi)
    if dpi <= 0: dpi = 150
    g.set_dpi(dpi)
    # If the user specifies a height and width for the canvas use it if it's a plausible
    # value. Anything smaller than 8.5 x 11 is set to that size. Values larger than 8.5 x 11
    # are left alone.
    if gheight and gwidth:
        gheight = float(gheight)
	gwidth = float(gwidth)
        if gwidth <= 8.5: gwidth = 8.5
        if gheight <= 11.: gwidth = 11.
        try:
            page = "%.1f,%.1f" % ( gwidth, gheight )
            size = "%.1f,%.1f" % ( gwidth-1, gheight-1 )
        except:
            page = "8.5,11"
            size = "7.5,10"
    else:
        page = "8.5,11"
        size = "7.5,10"
    g.set_page(page)
    g.set_size(size)
    # Handle the orientation of the document.
    if gorient == 'l':
        g.set_orientation("landscape")
    else:
        g.set_orientation("portrait")
    if gsize != 'l':
        g.set_ratio("auto")
    if gdirec not in ['RL', 'LR', 'TB', 'BT']:
        gdirec = 'TB'
    g.set_rankdir(gdirec)
    g.set_center('true')
    g.set_concentrate('true')
    g.set_ordering('out')
    if gformat not in g.formats:
        gformat = 'jpg'
    # If we do not have any generations, we have to draw a less-nice graph.
    if len(pedobj.metadata.unique_gen_list) <= 1:
        for _m in pedobj.pedigree:
            # Add a node for the current animal and set some properties.
            if str(_m.animalID) != str(pedobj.kw['missing_parent']):
                if gname:
                    _node_name = _m.name
                else:
                    _node_name = _m.animalID
                _an_node = pydot.Node(str(_node_name))
                _an_node.set_fontname('Helvetica')
                # _an_node.set_fontsize('10')
                _an_node.set_fontsize(str(gfontsize))
                _an_node.set_height('0.35')
                if _m.sex == 'M' or _m.sex == 'm':
                    _an_node.set_shape('box')
                elif _m.sex == 'F' or _m.sex == 'f':
                    _an_node.set_shape('ellipse')
                else:
                    pass
                g.add_node(_an_node)
                # Add the edges to the parent nodes, if any.
                if str(_m.sireID) != str(pedobj.kw['missing_parent']):
                    if gname:
                        if garrow:
                            g.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.sireID)-1].name), str(_m.name)))
                        else:
                            g.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.sireID)-1].name), str(_m.name), dir='none'))
                    else:
                        if garrow:
                            g.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.sireID)-1].originalID), str(_m.originalID)))
                        else:
                            g.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.sireID)-1].originalID), str(_m.originalID), dir='none'))
                if str(_m.damID) != str(pedobj.kw['missing_parent']):
                    if gname:
                        if garrow:
                            g.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.damID)-1].name), str(_m.name)))
                        else:
                            g.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.damID)-1].name), str(_m.name), dir='none'))
                    else:
                        if garrow:
                            g.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.damID)-1].originalID), str(_m.originalID)))
                        else:
                            g.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.damID)-1].originalID), str(_m.originalID), dir='none'))
    # Or test the new subgraph clusters
    elif gclusters:
        for _g in pedobj.metadata.unique_gen_list:
            _sg_anims = []
            _sg_name = 'sg%s' % (_g)
            if gshowall:
                sg = pydot.Cluster(graph_name=str(_sg_name), suppress_disconnected=False, simplify=True)
            else:
                sg = pydot.Cluster(graph_name=str(_sg_name), suppress_disconnected=True, simplify=True)
            for _m in pedobj.pedigree:
                if int(_m.gen) == int(_g):
                    _sg_anims.append(_m.animalID)
                    # Add a node for the current animal and set some properties.
                    if str(_m.animalID) != str(pedobj.kw['missing_parent']):
                        if gname:
                            _node_name = _m.name
                        else:
                            _node_name = str(_m.animalID)
                        _an_node = pydot.Node(str(_node_name))
                        _an_node.set_fontname('Helvetica')
                        _an_node.set_fontsize(gfontsize)
                        _an_node.set_height('0.35')
                        if _m.sex == 'M' or _m.sex == 'm':
                            _an_node.set_shape('box')
                        if _m.sex == 'F' or _m.sex == 'f':
                            _an_node.set_shape('ellipse')
                        sg.add_node(_an_node)
            g.add_subgraph(sg)
        # Now that we've added the clusters to the graph, define the
        # edges.
        for _m in pedobj.pedigree:
            if str(_m.sireID) != str(pedobj.kw['missing_parent']):
                if gname:
                    if garrow:
                        g.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.sireID)-1].name), str(_m.name)))
                    else:
                        g.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.sireID)-1].name), str(_m.name), dir='none'))
                else:
                    if garrow:
                        g.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.sireID)-1].animalID), str(_m.animalID)))
                    else:
                        g.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.sireID)-1].animalID), str(_m.animalID), dir='none'))
            if str(_m.damID) != str(pedobj.kw['missing_parent']):
                if gname:
                    if garrow:
                        g.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.damID)-1].name), str(_m.name)))
                    else:
                        g.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.damID)-1].name), str(_m.name), dir='none'))
                else:
                    if garrow:
                        g.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.damID)-1].animalID), str(_m.animalID)))
                    else:
                        g.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.damID)-1].animalID), str(_m.animalID), dir='none'))
    # Otherwise we can draw a nice graph.
    else:
        for _g in pedobj.metadata.unique_gen_list:
            _sg_anims = []
            _sg_name = 'sg%s' % (_g)
            if gshowall:
                sg = pydot.Subgraph(graph_name=str(_sg_name), suppress_disconnected=False, simplify=True, rank='same')
            else:
                sg = pydot.Subgraph(graph_name=str(_sg_name), suppress_disconnected=True, simplify=True, rank='same')
            for _m in pedobj.pedigree:
                if str(_m.animalID) != str(pedobj.kw['missing_parent']):
                    if int(_m.gen) == int(_g):
                        _sg_anims.append(_m.animalID)
                        ## ->
                        # Add a node for the current animal and set some properties.
                        if gname:
                            _node_name = _m.name
                        else:
                            _node_name = str(_m.animalID)
                        _an_node = pydot.Node(str(_node_name))
                        _an_node.set_fontname('Helvetica')
                        _an_node.set_fontsize(str(gfontsize))
                        _an_node.set_height('0.35')
                        if _m.sex == 'M' or _m.sex == 'm':
                            _an_node.set_shape('box')
                        if _m.sex == 'F' or _m.sex == 'f':
                            _an_node.set_shape('ellipse')
                        sg.add_node(_an_node)
                        ## <-
                        # Add the edges to the parent nodes, if any.
                        if str(_m.sireID) != str(pedobj.kw['missing_parent']):
                            if gname:
                                if garrow:
                                    sg.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.sireID)-1].name), str(_m.name)))
                                else:
                                    sg.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.sireID)-1].name), str(_m.name), dir='none'))
                            else:
                                if garrow:
                                    sg.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.sireID)-1].animalID), str(_m.animalID)))
                                else:
                                    sg.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.sireID)-1].animalID), str(_m.animalID), dir='none'))
                        if str(_m.damID) != str(pedobj.kw['missing_parent']):
                            if gname:
                                if garrow:
                                    sg.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.damID)-1].name), str(_m.name)))
                                else:
                                    sg.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.damID)-1].name), str(_m.name), dir='none'))
                            else:
                                if garrow:
                                    sg.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.damID)-1].animalID), str(_m.animalID)))
                                else:
                                    sg.add_edge(pydot.Edge(str(pedobj.pedigree[int(_m.damID)-1].animalID), str(_m.animalID), dir='none'))
                        ## <-
            if len(_sg_anims) > 0:
                _sg_list = ''
                for _a in _sg_anims:
                    if len(_sg_list) == 0:
                        _sg_list = '%s' % (_a)
                    else:
                        _sg_list = '%s,%s' % (_sg_list,_a)
            sg.set_rank(_sg_list)
            g.add_subgraph(sg)
            #try:
                #print sg.get_node_list()
            #except:
                #print 'Could not get node list for subgraph!'
            #try:
                #print sg.get_edge_list()
            #except:
                #print 'Could not get edge list for subgraph!'

    # For large graphs it is nice to write out the .dot file so that it does
    # not have to be recreated whenever new_draw_pedigree is called.
    # Especially when I am debugging.
    if gdot:
        dfn = '%s.dot' % (gfilename)
        #try:
        g.write(dfn)
        #except:
        #    if pedobj.kw['messages'] == 'verbose':
        #        print '[ERROR]: pyp_graphics/draw_pedigree() was unable to write the dotfile %s.' % (dfn)
        #    logging.error('pyp_graphics/draw_pedigree() was unable to draw the dotfile %s.', (dfn))

    # Write the graph to an output file.
    #try:
    outfile = '%s.%s' % (gfilename,gformat)
    g.write(outfile, prog='dot', format=gformat)
    #    return 1
    #except:
    #    outfile = '%s.%s' % (gfilename,gformat)
    #    if pedobj.kw['messages'] == 'verbose':
    #        print '[ERROR]: pyp_graphics/draw_pedigree() was unable to draw the pedigree %s.' % (outfile)
    #    logging.error('pyp_graphics/draw_pedigree() was unable to draw the pedigree %s.', (outfile))
    #    return 0

##
# founders_by_year() uses matplotlib -- if available on your system -- to produce a
# bar graph of the number (count) of founders in each birthyear.
# @param pedobj A PyPedal pedigree object.
# @param gfilename The name of the file to which the pedigree should be drawn
# @param gtitle The title of the graph.
# @retval A 1 for success and a 0 for failure.
def plot_founders_by_year(pedobj, gfilename='founders_by_year', gtitle='Founders by Birthyear'):
    """
    founders_by_year() uses matplotlib -- if available on your system -- to produce a
    bar graph of the number (count) of founders in each birthyear.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import pylab
    except ImportError:
        if pedobj.kw['messages'] == 'verbose':
            print '[ERROR]: pyp_graphics/plot_founders_by_year() was unable to import the matplotlib module!'
        logging.error('pyp_graphics/plot_founders_by_year() was unable to import the matplotlib module!')
        return 0

    fby = pyp_demog.founders_by_year(pedobj)
    #print fby
    try:
        pylab.clf()
        pylab.bar(fby.keys(),fby.values())
        pylab.title(gtitle)
        pylab.xlabel('Year')
        pylab.ylabel('Number of founders')
        plotfile = '%s.png' % (gfilename)
        #print 'plotfile: ', plotfile
        myplotfile = open(plotfile,'w')
        #print 'myplotfile: ', myplotfile
        pylab.savefig(myplotfile)
        myplotfile.close()
        return 1
    except:
        if pedobj.kw['messages'] == 'verbose':
            print '[ERROR]: pyp_graphics/plot_founders_by_year() was unable to create the plot \'%s\' (%s.png).' % ( gtitle, gfilename )
        #print gtitle, gfilename
        logging.error('pyp_graphics/plot_founders_by_year() was unable to create the plot \'%s\' (%s.png).', gtitle, gfilename)
        return 0

##
# founders_pct_by_year() uses matplotlib -- if available on your system -- to produce a
# line graph of the frequency (percentage) of founders in each birthyear.
# @param pedobj A PyPedal pedigree object.
# @param gfilename The name of the file to which the pedigree should be drawn
# @param gtitle The title of the graph.
# @retval A 1 for success and a 0 for failure.
def plot_founders_pct_by_year(pedobj, gfilename='founders_pct_by_year', gtitle='Founders by Birthyear'):
    """
    founders_pct_by_year() uses matplotlib -- if available on your system -- to produce a
    line graph of the frequency (percentage) of founders in each birthyear.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import pylab
    except ImportError:
        if pedobj.kw['messages'] == 'verbose':
            print '[ERROR]: pyp_graphics/plot_founders_pct_by_year() was unable to import the matplotlib module!'
        logging.error('pyp_graphics/plot_founders_pct_by_year() was unable to import the matplotlib module!')
        return 0

    fby = pyp_demog.founders_by_year(pedobj)
    _freqdict = {}
    for _k in fby.keys():
        _freqdict[_k] = float(fby[_k]) / float(pedobj.metadata.num_unique_founders)
    try:
        import matplotlib
        matplotlib.use('Agg')
        import pylab
        pylab.clf()
        pylab.plot(fby.keys(),_freqdict.values())
        pylab.title(gtitle)
        pylab.xlabel('Year')
        pylab.ylabel('% founders')
        plotfile = '%s.png' % (gfilename)
        myplotfile = open(plotfile,'w')
        pylab.savefig(myplotfile)
        myplotfile.close()
        return 1
    except:
        if pedobj.kw['messages'] == 'verbose':
            print '[ERROR]: pyp_graphics/plot_pct_founders_by_year() was unable to create the plot \'%s\' (%s.png).' % (gtitle,gfilename)
        logging.error('pyp_graphics/plot_pct_founders_by_year() was unable to create the plot \'%s\' (%s.png).' % (gtitle,gfilename))
        return 0

##
# pcolor_matrix_pylab() implements a matlab-like 'pcolor' function to
# display the large elements of a matrix in pseudocolor using the Python Imaging
# Library.
# @param A Input Numpy matrix (such as a numerator relationship matrix).
# @param fname Output filename to which to dump the graphics (default 'tmp.png')
# @retval A list of NewAnimal() objects; a pedigree metadata object.
def pcolor_matrix_pylab(A, fname='pcolor_matrix_matplotlib'):
    """
    pcolor_matrix_pylab() implements a matlab-like 'pcolor' function to
    display the large elements of a matrix in pseudocolor using the Python Imaging
    Library.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import pylab
    except ImportError:
        if pedobj.kw['messages'] == 'verbose':
            print '[ERROR]: pyp_graphics/pcolor_matrix_pylab() was unable to import the matplotlib module!'
        logging.error('pyp_graphics/pcolor_matrix_pylab() was unable to import the matplotlib module!')
        return 0

    try:
        import numpy
        pylab.clf()
        x = pylab.arange(A.shape[0])
        X, Y = pylab.meshgrid(x,x)

        xmin = min(pylab.ravel(X))
        xmax = max(pylab.ravel(X))
        pylab.xlim(xmin, xmax)
        ymin = min(pylab.ravel(Y))
        ymax = max(pylab.ravel(Y))
        pylab.ylim(ymin, ymax)
        pylab.axis('off')

        pylab.pcolor(X, Y, pylab.transpose(A))#, shading='flat')
        pylab.clim(0.0, 1.0)
        plotfile = '%s.png' % (fname)
        myplotfile = open(plotfile,'w')
        pylab.savefig(myplotfile)
        myplotfile.close()
        return 1
    except:
        if pedobj.kw['messages'] == 'verbose':
            print '[ERROR]: pyp_graphics/pcolor_matrix_pylab() was unable to create the plot %s.' % (plotfile)
        logging.error('pyp_graphics/pcolor_matrix_pylab() was unable to create the plot %s.', (plotfile))
        return 0

##
# spy_matrix_pylab() implements a matlab-like 'pcolor' function to
# display the large elements of a matrix in pseudocolor using the Python Imaging
# Library.
# @param A Input Numpy matrix (such as a numerator relationship matrix).
# @param fname Output filename to which to dump the graphics (default 'tmp.png')
# @retval A list of NewAnimal() objects; a pedigree metadata object.
def spy_matrix_pylab(A, fname='spy_matrix_matplotlib'):
    """
    spy_matrix_pylab() implements a matlab-like 'pcolor' function to
    display the large elements of a matrix in pseudocolor using the Python Imaging
    Library.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import pylab
    except ImportError:
        if pedobj.kw['messages'] == 'verbose':
            print '[ERROR]: pyp_graphics/spy_matrix_pylab() was unable to import the matplotlib module!'
        logging.error('pyp_graphics/spy_matrix_pylab() was unable to import the matplotlib module!')
        return 0

    try:
        import numpy
        pylab.clf()
        pylab.spy2(A)
        plotfile = '%s.png' % (fname)
        myplotfile = open(plotfile,'w')
        pylab.savefig(myplotfile)
        myplotfile.close()
        return 1
    except:
        if pedobj.kw['messages'] == 'verbose':
            print '[ERROR]: pyp_graphics/spy_matrix_pylab() was unable to create the plot %s.' % (plotfile)
        logging.error('pyp_graphics/spy_matrix_pylab() was unable to create the plot %s.', (plotfile))
        return 0

##
# plot_line_xy() uses matplotlib -- if available on your system -- to produce a
# line graph of the values in a dictionary for each level of key.
# @param xydict A Python dictionary of x- (keys) and y-values (values) to be plotted. 
# @param gfilename The name of the file to which the figure should be written.
# @param gtitle The title of the graph.
# @param gxlabel The label for the x-axis.
# @param gylabel The label for the y-axis.
# @param gformat The format in which the output file should be written  (JPG|PNG|PS).
# @retval A 1 for success and a 0 for failure.
def plot_line_xy(xydict, gfilename='plot_line_xy', gtitle='Value by key', gxlabel='X', gylabel='Y', gformat='png'):
    """
    plot_line_xy() uses matplotlib -- if available on your system -- to produce a
    line graph of the values in a dictionary for each level of key.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import pylab
    except ImportError:
        if pedobj.kw['messages'] == 'verbose':
            print '[ERROR]: pyp_graphics/plot_line_xy() was unable to import the matplotlib module!'
        logging.error('pyp_graphics/plot_line_xy() was unable to import the matplotlib module!')
        return 0

    if gformat not in ['png']:
        gformat = 'png'

    try:
        pylab.clf()
        # Dictionary keys are not guaranteed to be listed
        # in ascending order, so we have to do it.
        keylist = xydict.keys()
        keylist.sort()
        valuelist = [xydict[key] for key in keylist]
        # Done with the sorting
        pylab.plot(keylist,valuelist)
        pylab.title(gtitle)
        pylab.xlabel(gxlabel)
        pylab.ylabel(gylabel)
        plotfile = '%s.%s' % (gfilename, gformat)
        myplotfile = open(plotfile,'w')
        pylab.savefig(myplotfile)
        myplotfile.close()
        _status = 1
    except:
        if pedobj.kw['messages'] == 'verbose':
            print '[ERROR]: pyp_graphics/plot_line_xy() was unable to create the plot \'%s\' (%s.%s).' % (gtitle,gfilename,gformat)
        logging.error('pyp_graphics/plot_line_xy() was unable to create the plot \'%s\' (%s.%s).' % (gtitle,gfilename,gformat))
        _status = 0
    return _status

##
# new_draw_pedigree() uses the pygraphviz to produce a directed graph of your
# pedigree with paths of inheritance as edges and animals as nodes.  If there
# is more than one generation in the pedigree as determind by the "gen"
# attributes of the animals in the pedigree, draw_pedigree() will use subgraphs
# to try and group animals in the same generation together in the drawing.
# @param pedobj A PyPedal pedigree object.
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
# @retval A 1 for success and a 0 for failure.
def new_draw_pedigree(pedobj, gfilename='pedigree', gtitle='', gformat='jpg', \
    gsize='f', gdot=1, gorient='p', gdirec='', gname=0, garrow=1, \
    gtitloc='b', gtitjust='c', gshowall=1, gprog='dot'):
    """
    new_draw_pedigree() uses the pygraphviz to produce a directed graph of your
    pedigree with paths of inheritance as edges and animals as nodes.  If there
    is more than one generation in the pedigree as determined by the "gen"
    attributes of the animals in the pedigree, draw_pedigree() will use subgraphs
    to try and group animals in the same generation together in the drawing.

    Note that there is not (AFAIK) a nice way to get pygraphviz installed on vanilla
    Windows, so you'll have to make due with draw_pedigree() and plain old dot.
    """

    try:
        import pygraphviz
    except ImportError:
        if pedobj.kw['messages'] == 'verbose':
            print '[ERROR]: pyp_graphics/new_draw_pedigree() was unable to import the pygraphviz module!'
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
            print '[GRAPH]: The pedigree that you passed to pyp_graphics/new_draw_pedigree() is not renumbered. Because of this, there may be errors in the rendered pedigree. In order to insure that the pedigree drawing is accurate, you should renumber the pedigree before calling new_draw_pedigree().'
        logging.error('The pedigree that you passed to pyp_graphics/new_draw_pedigree() is not renumbered. Because of this, there may be errors in the rendered pedigree. In order to insure that the pedigree drawing is accurate, you should renumber the pedigree before calling new_draw_pedigree().')

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
        #print '[DEBUG]: sex = ', _m.sex
        if _m.sex == 'M' or _m.sex == 'm':
            n.attr['shape'] = 'box'
        elif _m.sex == 'F' or _m.sex == 'f':
            n.attr['shape'] = 'ellipse'
        else:
            n.attr['shape'] = 'octagon'
            #pass

        # Add the edges to the parent nodes, if any.
        if int(_m.sireID) != pedobj.kw['missing_parent']:
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
                #e = g.get_edge(_sire_edge,_anim_node)
		e = g.get_edge(_sire_edge,n)
                e.attr['dir'] = 'none'
        if int(_m.damID) != pedobj.kw['missing_parent']:
            if gname:
                _dam_edge = pedobj.pedigree[int(_m.damID)-1].name
            else:
                _dam_edge = pedobj.pedigree[int(_m.damID)-1].animalID
            g.add_edge(_dam_edge,_node_name)
            if not _tf[garrow]:
                #e = g.get_edge(_dam_edge,_anim_node)
		e = g.get_edge(_dam_edge,n)
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
                print '[ERROR]: pyp_graphics/new_draw_pedigree() was unable to write the dotfile %s.' % (dfn)
            logging.error('pyp_graphics/new_draw_pedigree() was unable to draw the dotfile %s.', (dfn))

    # Write the graph to an output file.
    try:
        outfile = '%s.%s' % (gfilename,gformat)
        g.draw(outfile,prog=gprog)
        return 1
    except:
        outfile = '%s.%s' % (gfilename,gformat)
        if pedobj.kw['messages'] == 'verbose':
            print '[ERROR]: pyp_graphics/new_draw_pedigree() was unable to draw the pedigree %s.' % (outfile)
        logging.error('pyp_graphics/new_draw_pedigree() was unable to draw the pedigree %s.', (outfile))
        return 0
