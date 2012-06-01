#!/usr/bin/python

###############################################################################
# NAME: pyp_gui.py
# VERSION: 2.0.0rc9 (31AUGUST2010)
# AUTHOR: John B. Cole, PhD (john.cole@ars.usda.gov)
# LICENSE: LGPL
###############################################################################
# class MainFrame(Frame)
#     def Body()
#     def CreateMenu()
#     def CreateToolbar()
#     def OnAbout()
#     def OnAppExit()
#     def OnGraphInbreed()
#     def OnMetricsEffectiveFounders()
#     def OnMetricsInbreeding()
#     def OnPedList()
#     def OnPedMeta()
#     def OnPedView()
#     def OnSettingsViewLog()
#     def Open()
#     def _OpenFile()
#     def Save()
#     def SetFilename()
#     def ToDo()
###############################################################################

##
# pyp_gui provides a simple graphical user interface for PyPedal using the wxPython
# (http://www.wxpython.org/) and Wax (http://sourceforge.net/projects/waxgui) extensions.
##

import PIL.Image
import logging, os
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
from PyPedal import pyp_gui_graphs
from PyPedal import pyp_gui_metrics
from PyPedal import pyp_gui_utils

# File menu items
ID_FILE_OPEN = 101
ID_FILE_SAVE = 102
ID_FILE_SAVEAS = 103
ID_FILE_CLOSE = 104
ID_FILE_EXIT = 110
# Pedigree menu items
ID_PED_META = 201
ID_PED_LIST = 202
ID_PED_VIEW = 203
# Metrics menu items
ID_METRIC_INBREEDING = 301
ID_METRIC_CALC_FA = 302
# Graphics menu items
ID_GRAPH_INBREEDING = 401
# Settings menu items
ID_SETTING_OPTIONS = 501
ID_SETTING_LOG = 502
# Help menu items
ID_HELP_HELP = 601
ID_HELP_ABOUT = 602
# Text box items
ID_MAIN_TB = 901

# Window sizes
MAIN_WINDOW_X = 640
MAIN_WINDOW_Y = 480

FIXED_FONT = ('Courier New', 10)

options = {}
options['pedformat'] = 'asdxb'
options['renumber'] = 1
options['messages'] = 'quiet'

##
# The MainFrame() class provides the main pyp_gui interface.
class MainWindow(wx.Frame):
    def __init__(self,parent,id,title):
        wx.Frame.__init__(self,parent,wx.ID_ANY, title,
            size = ( MAIN_WINDOW_X, MAIN_WINDOW_Y ),
            style = wx.DEFAULT_FRAME_STYLE |
                wx.NO_FULL_REPAINT_ON_RESIZE )
        self.textbox = wx.TextCtrl(self, 1,
            style = wx.TE_MULTILINE)
        self.CreateStatusBar()  # A Statusbar in the bottom of the window

        # This is the "File" menu
        filemenu = self.filemenu = wx.Menu()
        filemenu.Append(ID_FILE_OPEN, "&Open", " Open a pedigree file")
        filemenu.Append(ID_FILE_SAVE, "&Save", " Save a pedigree file")
        filemenu.Append(ID_FILE_SAVEAS, "Save &As", " Save a pedigree file as...")
        filemenu.Append(ID_FILE_CLOSE, "&Close", " Close a pedigree file")
        filemenu.AppendSeparator()
        filemenu.Append(ID_FILE_EXIT,"E&xit", " Exit PyPedal")

        # This is the "Pedigree" menu
        pedmenu= wx.Menu()
        pedmenu.Append(ID_PED_META, "&Metadata", " View pedigree metadata")
        pedmenu.Append(ID_PED_LIST, "&List Animals", " View a list of animal records")
        pedmenu.Append(ID_PED_VIEW, "&View", " View a diagram of the pedigree")

        # This is the "Metrics" menu
        metricmenu= wx.Menu()
        metricmenu.Append(ID_METRIC_INBREEDING, "&Inbreeding", " Calculate coefficients of inbreeding")
        metricmenu.Append(ID_METRIC_CALC_FA, "Effective &Founders", " Calculate the effective founder number")

        # This is the "Graphics" menu
        graphmenu= wx.Menu()
        graphmenu.Append(ID_GRAPH_INBREEDING, "&Inbreeding", " View inbreeding by birth year")

        # This is the "Settings" menu
        settingmenu= wx.Menu()
        settingmenu.Append(ID_SETTING_OPTIONS, "Pedigree &Options", " Pedigree options")
        settingmenu.AppendSeparator()
        settingmenu.Append(ID_SETTING_LOG, "View &Log", " View logfile")

        # This is the "Help" menu
        helpmenu= wx.Menu()
        helpmenu.Append(ID_HELP_HELP, "&Help", " Help with PyPedal")
        helpmenu.AppendSeparator()
        helpmenu.Append(ID_HELP_ABOUT,"&About", " About PyPedal")

        ###
        ### Options -- will change later to store in a configuration file
        ###
        self.pedigreedirty = False
        self.filetag='_test_gui_'
        self.sepchar=' '
        self.debug=0
        self.io='no'
        self.renum=1
        self.outformat='0'
        self.name='GUI Test Pedigree'
        self.alleles=0

        ###
        ### Creating the menubar.
        ###
        menuBar = self.menuBar = wx.MenuBar()
        menuBar.Append(filemenu,"&File") # Adding the "filemenu" to the MenuBar
        menuBar.Append(pedmenu,"&Pedigree") # Adding the "pedmenu" to the MenuBar
        menuBar.Append(metricmenu,"&Metrics")
        menuBar.Append(graphmenu,"&Graphics")
        menuBar.Append(settingmenu,"&Settings")
        menuBar.Append(helpmenu,"&Help") # Adding the "helpmenu" to the MenuBar
        self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
        # Attach handlers to events
        ## "File" menu events
        wx.EVT_MENU(self, ID_FILE_OPEN, self.OnFileOpenDialog)
        wx.EVT_MENU(self, ID_FILE_SAVE, self.ToDo)
        wx.EVT_MENU(self, ID_FILE_SAVEAS, self.ToDo)
        wx.EVT_MENU(self, ID_FILE_CLOSE, self.ToDo)
        wx.EVT_MENU(self, ID_FILE_EXIT, self.OnAppExit)
        ## "Pedigree" menu events
        wx.EVT_MENU(self, ID_PED_META, self.OnPedMeta)
        wx.EVT_MENU(self, ID_PED_LIST, self.OnPedList)
        wx.EVT_MENU(self, ID_PED_VIEW, self.OnPedView)
        ## "Metrics" menu events
        wx.EVT_MENU(self, ID_METRIC_INBREEDING, self.OnMetricsInbreeding)
        wx.EVT_MENU(self, ID_METRIC_CALC_FA, self.OnMetricsEffectiveFounders)
        ## "Graphics" menu events
        wx.EVT_MENU(self, ID_GRAPH_INBREEDING, self.OnGraphInbreed)
        ## "Settings" menu events
        wx.EVT_MENU(self, ID_SETTING_OPTIONS, self.OnSettingsOptions)
        wx.EVT_MENU(self, ID_SETTING_LOG, self.OnSettingsViewLog)
        ## "Help" menu events
        wx.EVT_MENU(self, ID_HELP_HELP, self.ToDo)
        wx.EVT_MENU(self, ID_HELP_ABOUT, self.OnAbout)
        # Show the window
        self.Show(True)

        ###
        ### Add a file history
        ###
        self.filehistory = wx.FileHistory()
        self.filehistory.UseMenu(filemenu)

        ###
        ### Add some event bindings for the file history.
        ###
        self.Bind(wx.EVT_RIGHT_UP, self.OnRightClick)
        self.Bind(wx.EVT_MENU, self.OnFileOpenDialog, id=ID_FILE_OPEN)
        self.Bind( wx.EVT_MENU_RANGE,
            self.OnFileHistory,
            id = wx.ID_FILE1,
            id2 = wx.ID_FILE9 )
        self.Bind(wx.EVT_WINDOW_DESTROY, self.Cleanup)

    def Cleanup(self, *args):
        # A little extra cleanup is required for the FileHistory control
        del self.filehistory
        #self.filemenu.Destroy()

    def OnRightClick(self, evt):
        self.PopupMenu(self.menu, evt.GetPosition())

    ##
    # OnFileOpenDialog() produces the file dialog used to select a pedigree file to be opened.
    # @param event The type of event causing invocation of this method.
    # @return None
    # @defreturn None
    def OnFileOpenDialog(self, event):
        dlg = wx.FileDialog(self,
            'PyPedal - Select a pedigree to open',
            '.',
            '',
            '*.ped',
            wx.OPEN)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                self._OpenFile(dlg.GetPath())
                logging.info('OnFileOpenDialog(): Pedigree file %s selected.\n', path)
                self.filehistory.AddFileToHistory(path)
        except:
            dlg.Destroy()


    def OnFileHistory(self, evt):
        # get the file based on the menu ID
        fileNum = evt.GetId() - wx.ID_FILE1
        path = self.filehistory.GetHistoryFile(fileNum)
        logging.info('OnFileHistory(): Pedigree file %s selected.\n', path)
        # add it back to the history so it will be moved up the list
        self.filehistory.AddFileToHistory(path)

    ##
    # _OpenFile() actually loads a pedigree from a file.  Once the pedigree is
    # loaded the pedigree metadata is written to the screen.
    # @param filename The pedigree file name.
    # @return None
    # @defreturn None
    def _OpenFile(self, filename):
        self.SetFilename(filename)
        f = open(filename, 'r')
        data = f.read()
        f.close()
        self.textbox.Clear()
        #self.textbox.AppendText(data)
        try:
            options['pedfile'] = filename
            self._pedigree = pyp_newclasses.NewPedigree(options)
            self._pedigree.load()
            self.textbox.AppendText(self._pedigree.metadata.stringme())
        except:
            pass

    ##
    # SetFilename() writes the name of the opened pedigree file to the statusbar.
    # @param filename The pedigree file name.
    # @return None
    # @defreturn None
    def SetFilename(self, filename):
        self.filename = filename
        self.SetStatusText(self.filename)
        #self.statusbar[1] = 'Filename: ' + str(self.filename)

    ##
    # Save() produces the file dialog used when a pedigree file is to be saved.
    # @param event The type of event causing invocation of this method.
    # @return None
    # @defreturn None
    def Save(self, event):
        pass

    ##
    # OnPedMeta() produces the file dialog used when a pedigree file is to be saved.
    # @param event The type of event causing invocation of this method.
    # @return None
    # @defreturn None
    def OnPedMeta(self, event):
        if hasattr(self,'_pedigree'):
            self.textbox.Clear()
            self.textbox.AppendText(self._pedigree.metadata.stringme())
        else:
            pyp_gui_utils.PyPedalShowErrorDialog(self,sedTitle='Display pedigree metadata', sedMessage='These data cannot be displayed because you have not yet loaded a pedigree!')

    ##
    # OnPedList() presents a list of the animals in the pedigree.
    # @param event The type of event causing invocation of this method.
    # @return None
    # @defreturn None
    def OnPedList(self, event):
        pass

    ##
    # OnPedView() uses pydot and graphviz to produce a drawing of the pedigree,
    # if they are installed.
    # @param event The type of event causing invocation of this method.
    # @return None
    # @defreturn None
    def OnPedView(self, event):
        pass

    ##
    # OnAbout() produces an "about this application" dialogue box.
    # @param event The type of event causing invocation of this method.
    # @return None
    # @defreturn None
    def OnAbout(self, event):
        dlg = wx.MessageDialog(self, ' PyPedal 2.0.0a19\n'
            ' A software package for pedigree analysis.\n'
            ' (c) 2002-2005 John B. Cole\n'
            ' http://pypedal.sourceforge.net/\n'
            ' jcole@aipl.arsusda.gov',
            'About PyPedal',
            wx.OK |
            wx.ICON_INFORMATION)
        # Create a message dialog box
        dlg.ShowModal() # Shows it
        dlg.Destroy() # finally destroy it when finished.

    ##
    # OnAppExit() produces a "do you really want to exit this application" dialogue box.
    # @param event The type of event causing invocation of this method.
    # @return None
    # @defreturn None
    def OnAppExit(self, event):
        changed = 0
        dlg = wx.MessageDialog(self,
            'Are you sure you want to exit PyPedal?',
            'Exit PyPedal?',
            wx.YES_NO |
            wx.ICON_QUESTION)
        if dlg.ShowModal() == wx.ID_YES:
            changed = 1
        dlg.Destroy()
        if changed:
            self.Destroy()

    ##
    # OnSettingsViewLog() dumps the logfile to a window.
    # @param event The type of event causing invocation of this method.
    # @return None
    # @defreturn None
    def OnSettingsViewLog(self, event):
        pyp_gui_utils.UtilsViewLog(self)

    ##
    # ToDo() produces a "this feature has not yet been implemented" dialogue box.
    # @param event The type of event causing invocation of this method.
    # @return None
    # @defreturn None
    def ToDo(self,event):
        dlg = wx.MessageDialog(self,
            'This feature has not yet been implemented!',
            'PyPedal - Unknown Feature',
            wx.OK |
            wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()

    ##
    # OnGraphInbreed() produces a "do you really want to exit this application" dialogue box.
    # @param event The type of event causing invocation of this method.
    # @return None
    # @defreturn None
    def OnGraphInbreed(self,event):
#         try:
        if hasattr(self,'_pedigree'):
            if not self._pedigree.kw['f_computed']:
                _ped_inbr = pyp_nrm.inbreeding(self._pedigree)
            pyp_db.loadPedigreeTable(self._pedigree)
            _coi_by_year = pyp_reports.meanMetricBy(self._pedigree,metric='fa',byvar='by')
            del(_coi_by_year[1900])
            _graph = pyp_graphics.plot_line_xy(_coi_by_year,gfilename='_coi_by_year',gtitle='Average coefficients of inbreeding',gxlabel='Birth year',gylabel='f')
            if _graph:
#                 wx.InitAllImageHandlers()
#                 _i = wx.Image('_coi_by_year.png', wx.BITMAP_TYPE_PNG)
#                 _b = wx.BitmapFromImage(_i)
#                 _sb = wx.StaticBitmap(self, -1, _b, (10, 10), (_b.GetWidth(), _b.GetHeight()))
                dlg = pyp_gui_graphs.PyPedalGraphDialogInbreeding(self,ID_MAIN_TB,'PyPedal - View inbreeding by birth year')
                try:
                    dlg.ShowModal()
                finally:
                    dlg.Destroy()
#             else:
#                 pyp_gui_utils.PyPedalShowErrorDialog(self,sedTitle='View inbreeding by birth year', sedMessage='This graph cannot be displayed because you have not yet loaded a pedigree!')
#         except:
#                 pyp_gui_utils.PyPedalShowErrorDialog(self,sedTitle='View inbreeding by birth year', sedMessage='An unknown error occurred while trying to display this graph!')

    ##
    # OnMetricsInbreeding() writes summary statistics for inbreeding to the main textbox.
    # @param event The type of event causing invocation of this method.
    # @return None
    # @defreturn None
    def OnMetricsInbreeding(self,event):
        pyp_gui_metrics.MetricsInbreeding(self)

    ##
    # OnMetricsEffectiveFounders() writes effective founders information to the
    # main textbox.
    # @param event The type of event causing invocation of this method.
    # @return None
    # @defreturn None
    def OnMetricsEffectiveFounders(self,event):
        pyp_gui_metrics.MetricsEffectiveFounders(self)

    def OnSettingsOptions(self,event):
        pass

class MyApp(wx.App):
    def OnInit(self):
        frame = MainWindow(None, -1, "PyPedal")
        frame.Show(True)
        self.SetTopWindow(frame)
        return True

if __name__ == '__main__':
    app = MyApp(0)
    app.MainLoop()