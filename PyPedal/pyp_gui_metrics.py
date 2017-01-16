#!/usr/bin/python

###############################################################################
# NAME: pyp_gui_metrics.py
# VERSION: 2.0.0rc9 (31AUGUST2010)
# AUTHOR: John B. Cole, PhD (john.cole@ars.usda.gov)
# LICENSE: LGPL
###############################################################################
# MetricsInbreeding()
# MetricsEffectiveFounders()
###############################################################################

##
# pyp_gui_metrics contains convenience functions for entries in the Metrics
# menu to reduce repetitive code in pyp_gui.
##
from __future__ import print_function
import PIL.Image
import os
import wx
from PyPedal import pyp_demog
from PyPedal import pyp_graphics
from PyPedal import pyp_io
from PyPedal import pyp_newclasses
from PyPedal import pyp_nrm
from PyPedal import pyp_metrics
from PyPedal import pyp_utils
from PyPedal import pyp_db
from PyPedal import pyp_reports
#from PyPedal import pyp_gui
#from PyPedal import pyp_gui_graphs
#from PyPedal import pyp_gui_utils

def MetricsInbreeding(self):
    if hasattr(self,'_pedigree'):
        self._inbreeding = pyp_nrm.inbreeding(self._pedigree)
        self.textbox.Clear()
        #print(pyp_io.summary_inbreeding(self._inbreeding['metadata']))
        self.textbox.AppendText(pyp_io.summary_inbreeding(self._inbreeding['metadata']))
    else:
        pyp_gui_graphs.PyPedalShowErrorDialog(self,sedTitle='Calculate inbreeding', sedMessage='This calculation cannot be displayed because you have not yet loaded a pedigree!')

def MetricsEffectiveFounders(self):
    if hasattr(self,'_pedigree'):
        self._inbreeding = pyp_nrm.inbreeding(self._pedigree)
        self.textbox.Clear()
        print(pyp_io.summary_inbreeding(self._inbreeding['metadata']))
        self.textbox.AppendText(pyp_io.summary_inbreeding(self._inbreeding['metadata']))
    else:
        pyp_gui_graphs.PyPedalShowErrorDialog(self,sedTitle='Calculate effective founder', sedMessage='This calculation cannot be displayed because you have not yet loaded a pedigree!')
