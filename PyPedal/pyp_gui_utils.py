#!/usr/bin/python

###############################################################################
# NAME: pyp_gui_utils.py
# VERSION: 2.0.0rc9 (31AUGUST2010)
# AUTHOR: John B. Cole, PhD (john.cole@ars.usda.gov)
# LICENSE: LGPL
###############################################################################
# PyPedalShowErrorDialog()
# PyPedalOptionsDialog() (I know, dialogue is spelled incorrectly.)
# UtilsViewLog()
###############################################################################

##
# pyp_gui_utils various bits and pieces for the PyPedal GUI, such as the class
# definition for the options dialogue.
##

import PIL.Image
import os, string
import wx
#from PyPedal import pyp_gui
#from PyPedal import pyp_gui_graphs
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
# The PyPedalOptionsDialog() class provides the dialogue box used for viewing and
# setting options.
class PyPedalOptionsDialog(wx.Dialog):
    def Body(self):
        self.AddComponent(wx.Bitmap(self,_b))
        self.Pack()

##
# PyPedalShowErrorDialog()
# @param sedTitle String used to construct the title of the dialogue box.
# @param sedMessage Message to be displayed in the dialogue box.
# @return None
# @defreturn None
def PyPedalShowErrorDialog(self,sedTitle="PyPedal - Unknown error", sedMessage="An unknown error occurred!"):
    sedTitle = 'PyPedal - %s' % (sedTitle)
    dlg = wx.MessageDialog(self,
        sedMessage,
        sedTitle,
        wx.OK |
        wx.ICON_INFORMATION)
    dlg.ShowModal()
    dlg.Destroy()

##
# UtilsViewLog() produces the file dialog used to select a log file to be
# opened.
# @param None
# @return None
# @defreturn None
def UtilsViewLog(self):
    if hasattr(self,'_pedigree'):
        self.SetFilename(self._pedigree.kw['logfile'])
        f = open(self._pedigree.kw['logfile'], 'r')
        data = f.read()
        f.close()
        shortfilename = string.split(self._pedigree.kw['logfile'],'/')[-1]
        _header = 'Viewing logfile %s\n' % (shortfilename)
        self.textbox.Clear()
        self.textbox.AppendText(_header)
        self.textbox.AppendText(pyp_io.LINE1)
        self.textbox.AppendText('\n')
        self.textbox.AppendText(data)
    else:
        pyp_gui_utils.PyPedalShowErrorDialog(self,sedTitle='View log', sedMessage='You cannot view a log file until you load a pedigree!')
