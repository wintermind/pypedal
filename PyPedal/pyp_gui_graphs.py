#!/usr/bin/python

###############################################################################
# NAME: pyp_gui_graphs.py
# VERSION: 2.0.0rc9 (31AUGUST2010)
# AUTHOR: John B. Cole, PhD (john.cole@ars.usda.gov)
# LICENSE: LGPL
###############################################################################
# PyPedalGraphDialogInbreeding()
###############################################################################

##
# I am not smart enough to properly subclass Wax's DIalog class in the manner that
# I would like.  As a result, pyp_gui_graphs contains a subclass of Dialog for each
# type of graph that pyp_gui can draw.
##

import wx
#from PyPedal import pyp_gui
#from PyPedal import pyp_gui_metrics
#from PyPedal import pyp_gui_utils
from PyPedal import pyp_demog
from PyPedal import pyp_graphics
from PyPedal import pyp_io
from PyPedal import pyp_newclasses
from PyPedal import pyp_nrm
from PyPedal import pyp_metrics
from PyPedal import pyp_utils
from PyPedal import pyp_db
from PyPedal import pyp_reports

##
# The PyPedalGraphDialogInbreeding() class provides the dialogue box used
# to display the "inbreeding by birthyear" graph.
class PyPedalGraphDialogInbreeding(wx.Dialog):
    def Body(self):
        print 'In PyPedalGraphDialogInbreeding.Body()'
        print '\tCalling wx.InitAllImageHandlers()'
        wx.InitAllImageHandlers()
        print '\tCalling wx.Image()'
        _i = wx.Image('_coi_by_year.png', wx.BITMAP_TYPE_PNG)
        print '\tCalling wx.BitmapFromImage()'
        _b = wx.BitmapFromImage(_i)
        print '\tCalling wx.StaticBitmap()'
        _sb = wx.StaticBitmap(self, -1, _b, (10, 10), (_b.GetWidth(), _b.GetHeight()))
        #wx.StaticBitmap(panel, -1, bmp, (10, pos), (bmp.GetWidth(), bmp.GetHeight()))
        self.AddComponent(_sb)
        self.Pack()

