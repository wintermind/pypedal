###############################################################################
# NAME: pyp_reports.py
# VERSION: 2.0.0b7 (29SEPTEMBER2010)
# AUTHOR: John B. Cole, PhD (john.cole@ars.usda.gov)
# LICENSE: LGPL
###############################################################################
# FUNCTIONS:
#   meanMetricBy()
#   pdfMeanMetricBy()
#   pdfPedigreeMetadata()
#   pdf3GenPed()
#   _pdfInitialize()
#   _pdfDrawPageFrame()
#   _pdfCreateTitlePage()
###############################################################################

## @package pyp_reports
# pyp_reports contains a set of procedures for ...
##
from __future__ import print_function
import logging, os, string, sys
import pyp_db
import pyp_io
import pyp_nrm
import pyp_utils
import pyp_network
import networkx

from reportlab.lib.pagesizes import letter, A4
from reportlab.lib.units import inch, cm
from reportlab.pdfgen import canvas

global metric_to_column
global byvar_to_column

metric_to_column = {'fa':'coi'}
byvar_to_column = {'by':'birthyear',
                  'gen':'generation'}

##
# meanMetricBy() returns a dictionary of means keyed by levels of the 'byvar' that
# can be used to draw graphs or prepare reports of summary statistics.
# @param pedobj A PyPedal pedigree object.
# @param metric The variable to summarize on a BY variable.
# @param byvar The variable on which to group the metric.
# @param createpdf Flag indicating whether or not a PDF version of the report should be created.
# @param conn Handle to a database connection, or False.
# @retval A dictionary containing means for the metric variable keyed to levels of the byvar.
def meanMetricBy(pedobj, metric='fa', byvar='by', createpdf=False, conn=False):
    """
    meanMetricBy() returns a dictionary of means keyed by levels of the 'byvar' that
    can be used to draw graphs or prepare reports of summary statistics.
    """
    # If the user doesn't pass us a db_dict then try and connect to the database
    conn_created = False
    result_dict = {}
    if conn == False:
        conn = pyp_db.connectToDatabase(pedobj)
        conn_created = True
    try:
        if metric not in ['fa']:
            logging.warning('You passed an unrecognized variable, %s, to pyp_reports/animalMetricBy() in the METRIC field.  It has been changed to \'fa\' (coefficient of inbreeding).', metric)
            metric = 'fa'
        if byvar not in ['gen','sex','birthyear','by']:
            logging.warning('You passed an unrecognized variable, %s, to pyp_reports/animalMetricBy() in the BYVAR field.  It has been changed to \'by\' (birth year).', byvar)
            byvar = 'by'

        if pyp_db.doesTableExist(pedobj, conn=conn):
            sql = 'select %s, avg(%s) from %s group by %s order by %s' % ( byvar_to_column[byvar], \
                metric_to_column[metric], pedobj.kw['database_table'], byvar_to_column[byvar], \
                byvar_to_column[byvar] )
            #print(sql)
            cursor = conn.Execute(sql)
            while not cursor.EOF:
                result_dict[cursor.fields[0]] = cursor.fields[1]
                cursor.MoveNext()
            cursor.Close()
            logging.info('pyp_reports/meanMetricBy() report completed.')
        else:
            logging.error('pyp_reports/meanMetricBy() report failed!')
        if conn_created:
            conn.Close()

        if createpdf:
            try:
                mmbPdfTitle = '%s_mean_metric_%s_%s' % \
                    (pedobj.kw['default_report'], metric, byvar)
                _mmbPdf = pdfMeanMetricBy(pedobj, result_dict, 1, mmbPdfTitle)
                if _mmbPdf:
                    logging.info('pyp_reports/pdfMeanMetricBy() succeeded.')
                else:
                    logging.error('pyp_reports/pdfMeanMetricBy() failed.')
            except:
                pass
    except:
        pass
    return result_dict

##
# pdfMeanMetricBy() returns a dictionary of means keyed by levels of the 'byvar' that
# can be used to draw graphs or prepare reports of summary statistics.
# @param pedobj A PyPedal pedigree object.
# @param results A dictionary containing means for the metric variable keyed to levels of the byvar.
# @param titlepage Show (1) or hide (0) the title page.
# @param reporttitle Title of report; if '', _pdfTitle is used.
# @param reportauthor Author/preparer of report.
# @param reportfile Optional name of file to which the report should be written.
# @retval 1 on success, 0 on failure
def pdfMeanMetricBy(pedobj, results, titlepage=0, reporttitle='', reportauthor='', reportfile=''):
    """
    pdfMeanMetricBy() returns a dictionary of means keyed by levels of the 'byvar' that
    can be used to draw graphs or prepare reports of summary statistics.
    """
    try:
        import reportlab
    except ImportError:
        logging.error('Unable to import ReportLab in pyp_reports/pdfMeanMetricBy().')
        return 0

    try:
        if reportfile == '':
            _pdfOutfile = '%s_mean_metric_by.pdf' % ( pedobj.kw['default_report'] )
        else:
            _pdfOutfile = reportfile
        if pedobj.kw['messages'] == 'verbose':
            print('Writing meanMetricBy report to %s' % ( _pdfOutfile ))
        logging.info('Writing meanMetricBy report to %s', _pdfOutfile )

        _pdfSettings = _pdfInitialize(pedobj)
        canv = canvas.Canvas(_pdfOutfile, invariant=1)
        canv.setPageCompression(1)

        if titlepage:
            if reporttitle == '':
                reporttitle = 'meanMetricBy Report for Pedigree\n%s' % (pedobj.kw['pedname'])
            _pdfCreateTitlePage(canv, _pdfSettings, reporttitle, reportauthor)
        _pdfDrawPageFrame(canv, _pdfSettings)

        canv.setFont("Times-Bold", 12)
        tx = canv.beginText( _pdfSettings['_pdfCalcs']['_left_margin'],
            _pdfSettings['_pdfCalcs']['_top_margin'] - 0.5 * _pdfSettings['_pdfCalcs']['_unit'] )
        # This is where the actual content is written to a text object that will be displayed on
        # a canvas.
        for _k, _v in results.iteritems():
            if len(str(_k)) <= 14:
                _line = '\t%s:\t\t%s' % (_k, _v)
            else:
                _line = '\t%s:\t%s' % (_k, _v)
            tx.textLine(_line)
            if tx.getY() < _pdfSettings['_pdfCalcs']['_bottom_margin'] + \
                0.5 * _pdfSettings['_pdfCalcs']['_unit']:
                canv.drawText(tx)
                canv.showPage()
                _pdfDrawPageFrame(canv, _pdfSettings)
                canv.setFont('Times-Roman', 12)
                tx = canv.beginText( _pdfSettings['_pdfCalcs']['_left_margin'],
                    _pdfSettings['_pdfCalcs']['_top_margin'] -
                    0.5 * _pdfSettings['_pdfCalcs']['_unit'] )
        if tx:
            canv.drawText(tx)
            canv.showPage()
        canv.save()
        return 1
    except:
        return 0

##
# pdfPedigreeMetadata() produces a report, in PDF format, of the metadata from
# the input pedigree.  It is intended for use as a template for custom printed
# reports.
# @param pedobj A PyPedal pedigree object.
# @param titlepage Show (1) or hide (0) the title page.
# @param reporttitle Title of report; if '', _pdfTitle is used.
# @param reportauthor Author/preparer of report.
# @param reportfile Optional name of file to which the report should be written.
# @retval A 1 on success, 0 otherwise.
def pdfPedigreeMetadata(pedobj, titlepage=0, reporttitle='', reportauthor='', reportfile=''):
    """
    pdfPedigreeMetadata() produces a report, in PDF format, of the metadata from
    the input pedigree.  It is intended for use as a template for custom printed
    reports.
    """
    try:
        import reportlab
    except ImportError:
        logging.error('Unable to import ReportLab in pyp_reports/pdfPedigreeMetadata().')
        return 0

    if reportfile == '':
        _pdfOutfile = '%s_metadata.pdf' % ( pedobj.kw['default_report'] )
    else:
        _pdfOutfile = reportfile
    if pedobj.kw['messages'] == 'verbose':
        print('Writing metadata report to %s' % ( _pdfOutfile ))
    logging.info('Writing metadata report to %s', _pdfOutfile )

    # The _pdfSettings dictionary contains several settings, such as page size,
    # page height and width, and margins, that are used several times.
    _pdfSettings = _pdfInitialize(pedobj)

    # The actual report is written to an instance of a canvas object, which is
    # stored in a file whose name is _pdfOutfile.  We have to create this canvas
    # before we can start assembling our report.
    canv = canvas.Canvas(_pdfOutfile, invariant=1)
    canv.setPageCompression(1)

    # Add a title page to the report if the user wants one.
    if titlepage:
        if reporttitle == '':
            reporttitle = 'Metadata for Pedigree\n%s' % (pedobj.kw['pedname'])
        _pdfCreateTitlePage(canv, _pdfSettings, reporttitle, reportauthor)

    # Start a new page of output.  Split the metadata output returned by the
    # stringme() method on linebreak characters and write each of the resulting
    # tokens to the canvas.
    _pdfDrawPageFrame(canv, _pdfSettings)
    canv.setFont("Times-Bold", 12)
    tx = canv.beginText( _pdfSettings['_pdfCalcs']['_left_margin'],
        _pdfSettings['_pdfCalcs']['_top_margin'] - 0.5 * _pdfSettings['_pdfCalcs']['_unit'] )
    _metadata_string = string.split(pedobj.metadata.stringme(), '\n')
    for _m in _metadata_string:
        tx.textLine(_m)
        # Page breaking has to be done manually.  Oh, well, I was able to steal, er,
        # borrow this but from odyssey.py.
        if tx.getY() < _pdfSettings['_pdfCalcs']['_bottom_margin'] + \
            0.5 * _pdfSettings['_pdfCalcs']['_unit']:
            canv.drawText(tx)
            canv.showPage()
            _pdfDrawPageFrame(canv, _pdfSettings)
            canv.setFont('Times-Roman', 12)
            tx = canv.beginText( _pdfSettings['_pdfCalcs']['_left_margin'],
                _pdfSettings['_pdfCalcs']['_top_margin'] -
                0.5 * _pdfSettings['_pdfCalcs']['_unit'] )

            pg = canv.getPageNumber()
            if pedobj.kw['messages'] == 'verbose' and pg % 10 == 0:
                print('Printed page %d' % pg)

    # Finish the document.  You need to draw the text object to the canvas and
    # call the showPage() method on the canvas so that your changes are visible.
    if tx:
        canv.drawText(tx)
        canv.showPage()
        # The following line inserts a blank page.  It was in the odyssey.py
        # example distributed with ReportLab, but I'm not sure why.  If the
        # last page of your document is getting chopped off try uncommenting
        # this line.
        #_pdfDrawPageFrame(canv, _pdfSettings)

    # Once your report has been assembled it must be written to disc.  If you
    # omit the 'canv.save()' line then your PDF will never be written to a file.
    canv.save()

##
# pdf3GenPed() draws a three-generation pedigree for animal 'animalID'.
# @param animalID An animal ID or list of animal IDs.
# @param pedobj A PyPedal pedigree object.
# @param titlepage Show (1) or hide (0) the title page.
# @param reporttitle Title of report; if '', _pdfTitle is used.
# @param reportauthor Author/preparer of report.
# @param reportfile Optional name of file to which the report should be written.
# @retval 1 on success, 0 on failure
def pdf3GenPed(animalID, pedobj, titlepage=0, reporttitle='', reportauthor='', reportfile=''):
    """
    pdf3GenPed() draws a three-generation pedigree for animal 'animalID'.
    """

    try:
        import reportlab
    except ImportError:
        logging.error('Unable to import ReportLab in pyp_reports/pdf3GenPed().')
        return 0

    # If the user specified only a single ID as an integer
    # we need to put it in a list.
    if type(animalID) != list:
        animalID = [animalID]

    if reportfile == '':
        _pdfOutfile = 'three_generation_pedigrees.pdf'
    else:
        _pdfOutfile = reportfile
    if pedobj.kw['messages'] == 'verbose':
        print('Writing 3GenPed to %s' % ( _pdfOutfile ))
    logging.info('Writing 3GenPed to %s', _pdfOutfile )

    _pdfSettings = _pdfInitialize(pedobj)
    #print(_pdfSettings)
    canv = canvas.Canvas(_pdfOutfile, invariant=1)
    canv.setPageCompression(1)

    if titlepage:
        if reporttitle == '':
            reporttitle = 'Three-generation Pedigrees'
        _pdfCreateTitlePage(canv, _pdfSettings, reporttitle, reportauthor)
        _pdfDrawPageFrame(canv, _pdfSettings)

    # This is where the actual content is written to a text object
    # that will be displayed on a canvas.
    try:
        for _anid in animalID:
            # Each pedigree is on a separate page, and page settings do
            # not persist across pages.
            canv.setFont("Times-Bold", 12)
            tx = canv.beginText( _pdfSettings['_pdfCalcs']['_left_margin'],
            _pdfSettings['_pdfCalcs']['_top_margin'] - 0.5 * \
                _pdfSettings['_pdfCalcs']['_unit'] )
            _pdfDrawPageFrame(canv, _pdfSettings)
            # We need to know the animal's renumbered ID. If the pedigree
            # we are using was based on animal names rather than integral
            # IDs we need to map the name back to the correct renumbered
            # ID.
            if 'A' in pedobj.kw['pedformat']:
                _anidx = pedobj.idmap[pedobj.namemap[_anid]] - 1
                _anid = pedobj.backmap[pedobj.idmap[pedobj.namemap[_anid]]]
            else:
                _anidx = pedobj.idmap[_anid] - 1
            # Oh, sure, it's ugly and inelegant and not easily extensible,
            # but it gets the job done. There are 15 possible animals
            # in a three-generation pedigree, and with the data structures
            # I'm using it's easiest to just populate the dictionary by
            # hand.
            _places = {}
            _places['a'] = pedobj.idmap[_anid]
            # Go down the sire side of the pedigree
            _places['s'] = pedobj.pedigree[_places['a']-1].sireID
            if int(_places['s']) != pedobj.kw['missing_parent']:
                _places['ss'] = pedobj.pedigree[_places['s']-1].sireID
                _places['sd'] = pedobj.pedigree[_places['s']-1].damID
            else:
                _places['ss'] = pedobj.kw['missing_parent']
                _places['sd'] = pedobj.kw['missing_parent']
            if int(_places['ss']) != pedobj.kw['missing_parent']:
                _places['sss'] = pedobj.pedigree[_places['ss']-1].sireID
                _places['ssd'] = pedobj.pedigree[_places['ss']-1].damID
            else:
                _places['sss'] = pedobj.kw['missing_parent']
                _places['ssd'] = pedobj.kw['missing_parent']
            if int(_places['sd']) != pedobj.kw['missing_parent']:
                _places['sds'] = pedobj.pedigree[_places['sd']-1].sireID
                _places['sdd'] = pedobj.pedigree[_places['sd']-1].damID
            else:
                _places['sds'] = pedobj.kw['missing_parent']
                _places['sdd'] = pedobj.kw['missing_parent']
            # Go down the dam side of the pedigree
            _places['d'] = pedobj.pedigree[_places['a']-1].damID
            if _places['d'] != pedobj.kw['missing_parent']:
                _places['ds'] = pedobj.pedigree[_places['d']-1].sireID
                _places['dd'] = pedobj.pedigree[_places['d']-1].damID
            else:
                _places['ds'] = pedobj.kw['missing_parent']
                _places['dd'] = pedobj.kw['missing_parent']
            if _places['ds'] != pedobj.kw['missing_parent']:
                _places['dss'] = pedobj.pedigree[_places['ds']-1].sireID
                _places['dsd'] = pedobj.pedigree[_places['ds']-1].damID
            else:
                _places['dss'] = pedobj.kw['missing_parent']
                _places['dsd'] = pedobj.kw['missing_parent']
            if _places['dd'] != pedobj.kw['missing_parent']:
                _places['dds'] = pedobj.pedigree[_places['dd']-1].sireID
                _places['ddd'] = pedobj.pedigree[_places['dd']-1].damID
            else:
                _places['dds'] = pedobj.kw['missing_parent']
                _places['ddd'] = pedobj.kw['missing_parent']
            #for k,v in _places.iteritems():
            #        print(k,v)

            _sill_width = _pdfSettings['_pdfCalcs']['_frame_width'] * 0.25
            _16 = _pdfSettings['_pdfCalcs']['_frame_height'] / 16.
            _64 = _pdfSettings['_pdfCalcs']['_frame_height'] / 64.
            _x = _pdfSettings['_pdfCalcs']['_left_margin']
            _y = _pdfSettings['_pdfCalcs']['_bottom_margin']

            # We're going to use a lookup table to store the x and y
            # coordinates for drawing sills. This will allow us to use
            # a single loop to populate the pedigree.
            os = {}
            for k in _places.keys(): os[k] = {}
            os['a']['x'] = _x; os['a']['y'] = _y+(8*_16)
            os['d']['x'] = _x+_sill_width; os['d']['y'] = _y+(4*_16)
            os['s']['x'] = _x+_sill_width; os['s']['y'] = _y+(12*_16)
            os['dd']['x'] = _x+(2*_sill_width); os['dd']['y'] = _y+(2*_16)
            os['ds']['x'] = _x+(2*_sill_width); os['ds']['y'] = _y+(6*_16)
            os['sd']['x'] = _x+(2*_sill_width); os['sd']['y'] = _y+(10*_16)
            os['ss']['x'] = _x+(2*_sill_width); os['ss']['y'] = _y+(14*_16)
            os['ddd']['x'] = _x+(3*_sill_width); os['ddd']['y'] = _y+(1*_16)
            os['dds']['x'] = _x+(3*_sill_width); os['dds']['y'] = _y+(3*_16)
            os['dsd']['x'] = _x+(3*_sill_width); os['dsd']['y'] = _y+(5*_16)
            os['dss']['x'] = _x+(3*_sill_width); os['dss']['y'] = _y+(7*_16)
            os['sdd']['x'] = _x+(3*_sill_width); os['sdd']['y'] = _y+(9*_16)
            os['sds']['x'] = _x+(3*_sill_width); os['sds']['y'] = _y+(11*_16)
            os['ssd']['x'] = _x+(3*_sill_width); os['ssd']['y'] = _y+(13*_16)
            os['sss']['x'] = _x+(3*_sill_width); os['sss']['y'] = _y+(15*_16)

            _line = 'Pedigree for %s (%s)' % ( pedobj.pedigree[_anidx].name, \
                pedobj.pedigree[_anidx].originalID )
            canv.setFont("Times-Bold", 12)
            canv.drawString( _x, _pdfSettings['_pdfCalcs']['_top_margin']-0.25*_16, _line )
            canv.setLineWidth(1)
            # Remember that the sills (and labels) are drawn from
            # bottom-to-top. The "+2" fudge factor in the drawString()
            # calls for _sill_text_1 is there so that the text sits a
            # hair above the bracket lines; it looks better that way.
            for k in _places.keys():
                canv.line( os[k]['x'], os[k]['y'], os[k]['x']+_sill_width, os[k]['y'] )
                if _places[k] == pedobj.kw['missing_parent']:
                    _sill_text_1 = '(%s)' % ( 'Unknown Parent' )
                    _sill_text_2 = ''
                else:
                    _sill_text_1 = '%s' % ( pedobj.pedigree[_places[k]-1].name )
                    _sill_text_2 = '(%s)' % ( pedobj.pedigree[_places[k]-1].originalID )
                canv.setFont("Times-Bold", 12)
                canv.drawString( os[k]['x'], os[k]['y']+2, _sill_text_1 )
                canv.setFont("Times-Roman", 12)
                canv.drawString( os[k]['x'], os[k]['y']-_64, _sill_text_2 )
                if k == 'a':
                    canv.setFont("Times-Roman", 12)
                    if pedobj.pedigree[_places[k]-1].herd == 'u':
                        _herd = '%s' % ( pedobj.kw['missing_herd'] )
                    else:
                        _herd = '%s' % ( pedobj.pedigree[_places[k]-1].herd )
                    _breed = '%s' % ( pedobj.pedigree[_places[k]-1].breed )
                    _inbreed = '%s' % ( pedobj.pedigree[_places[k]-1].fa )
                    _pedcomp = '%5.3f' % ( pedobj.pedigree[_places[k]-1].pedcomp )
                    canv.drawString( _x, _y+0.25*_16, 'Pedigree completeness:' )
                    canv.drawString( _x+_sill_width, _y+0.25*_16, _pedcomp )
                    canv.drawString( _x, _y+0.5*_16, 'Inbreeding:' )
                    canv.drawString( _x+_sill_width, _y+0.5*_16, _inbreed )
                    canv.drawString( _x, _y+0.75*_16, 'Breed:' )
                    canv.drawString( _x+_sill_width, _y+0.75*_16, _breed )
                    canv.drawString( _x, _y+_16, 'Herd:' )
                    canv.drawString( _x+_sill_width, _y+_16, _herd )
            # Draw the vertical lines that join the sills.
            ## Parents
            canv.line( os['s']['x'],os['d']['y'],os['s']['x'],os['s']['y'] )
            ## Maternal grandparents
            canv.line( os['ds']['x'],os['dd']['y'],os['ds']['x'],os['ds']['y'] )
            ## Paternal grandparents
            canv.line( os['ss']['x'],os['sd']['y'],os['ss']['x'],os['ss']['y'] )
            ## Mat-mat grandparents
            canv.line( os['dds']['x'],os['ddd']['y'],os['dds']['x'],os['dds']['y'] )
            ## Mat-pat grandparents
            canv.line( os['dss']['x'],os['dsd']['y'],os['dss']['x'],os['dss']['y'] )
            ## Pat-mat grandparents
            canv.line( os['sds']['x'],os['sdd']['y'],os['sds']['x'],os['sds']['y'] )
            ## Pat-pat grandparents
            canv.line( os['sss']['x'],os['ssd']['y'],os['sss']['x'],os['sss']['y'] )
            # Make sure we draw the page to the canvas when we've
            # finished with it.
            canv.drawText(tx)
            canv.showPage()
        canv.save()
        return 1
    except:
        logging.error('Unable to fetch the pedigree for animal %s in pyp_reports/pdf3GenPed()!',animalID)
        return 0

###############################################################################
# _pdfDrawPageFrame() was taken from the procedure drawPageFrame() included in
# the demo program odyssey.py in the ReportLab distribution.  _pdfInitialize()
# is rolled together using some of the code in odyssey.py as an example, as
# well as some of my own work.
###############################################################################

##
# _pdfInitialize() returns a dictionary of metadata that is used for report
# generation.
# @param pedobj A PyPedal pedigree object.
# @retval A dictionary of metadata that is used for report generation.
def _pdfInitialize(pedobj):
    """
    _pdfInitialize() returns a dictionary of metadata that is used for report
    generation.
    """
    _pdfSettings = {}
    _pdfSettings['_pdfTitle'] = pedobj.kw['pedname']
    _pdfSettings['_pdfPageinfo'] = pedobj.kw['filetag']
    # Calculate margins, etc.
    _pdfCalcs = {}
    if pedobj.kw['default_unit'] == 'inch':
        _pdfCalcs['_unit'] = inch
    else:
        _pdfCalcs['_unit'] = cm
    if pedobj.kw['paper_size'] == 'letter':
        _pdfCalcs['_page'] = letter
        _pdfCalcs['_top_margin'] = letter[1] - inch
        _pdfCalcs['_bottom_margin'] = inch
        _pdfCalcs['_left_margin'] = inch
        _pdfCalcs['_right_margin'] = letter[0] - inch
        _pdfCalcs['_frame_width'] = _pdfCalcs['_right_margin'] - _pdfCalcs['_left_margin']
        _pdfCalcs['_frame_height'] = _pdfCalcs['_top_margin'] - _pdfCalcs['_bottom_margin']
        _pdfCalcs['_page_width'] = letter[0]
        _pdfCalcs['_page_height'] = letter[1]
    else:
        _pdfCalcs['_page'] = A4
        _pdfCalcs['_top_margin'] = A4[1] - inch
        _pdfCalcs['_bottom_margin'] = inch
        _pdfCalcs['_left_margin'] = inch
        _pdfCalcs['_right_margin'] = A4[0] - inch
        _pdfCalcs['_frame_width'] = _pdfCalcs['_right_margin'] - _pdfCalcs['_left_margin']
        _pdfCalcs['_frame_height'] = _pdfCalcs['_top_margin'] - _pdfCalcs['_bottom_margin']
        _pdfCalcs['_page_width'] = A4[0]
        _pdfCalcs['_page_height'] = A4[1]
    _pdfSettings['_pdfCalcs'] = _pdfCalcs
    return _pdfSettings

##
# _pdfDrawPageFrame() nicely frames page contents and includes the
# document title in a header and the page number in a footer.
# @param canv An instance of a ReportLab Canvas object.
# @param _pdfSettings An options dictionary created by _pdfInitialize().
# @retval None
def _pdfDrawPageFrame(canv, _pdfSettings):
    """
    _pdfDrawPageFrame() nicely frames page contents and includes the
    document title in a header and the page number in a footer.
    """
    # Write the report title in the top, left-hand corner of the page.
    canv.line(_pdfSettings['_pdfCalcs']['_left_margin'],
        _pdfSettings['_pdfCalcs']['_top_margin'],
        _pdfSettings['_pdfCalcs']['_right_margin'],
        _pdfSettings['_pdfCalcs']['_top_margin'])
    canv.setFont('Times-Italic', 12)
    canv.drawString(_pdfSettings['_pdfCalcs']['_left_margin'],
        _pdfSettings['_pdfCalcs']['_top_margin'] + 2,
        _pdfSettings['_pdfTitle'])

    # Write the date/time in the top, right-hand corner of the page.
    canv.line(_pdfSettings['_pdfCalcs']['_left_margin'],
        _pdfSettings['_pdfCalcs']['_top_margin'],
        _pdfSettings['_pdfCalcs']['_right_margin'],
        _pdfSettings['_pdfCalcs']['_top_margin'])
    canv.setFont('Times-Italic', 12)
    canv.drawString(_pdfSettings['_pdfCalcs']['_right_margin'] - \
        1.85 * _pdfSettings['_pdfCalcs']['_unit'],
        _pdfSettings['_pdfCalcs']['_top_margin'] + 2,
        pyp_utils.pyp_nice_time())

    # Write the page number bottom center.
    canv.line(_pdfSettings['_pdfCalcs']['_left_margin'],
        _pdfSettings['_pdfCalcs']['_top_margin'],
        _pdfSettings['_pdfCalcs']['_right_margin'],
        _pdfSettings['_pdfCalcs']['_top_margin'])
    canv.line(_pdfSettings['_pdfCalcs']['_left_margin'],
        _pdfSettings['_pdfCalcs']['_bottom_margin'],
        _pdfSettings['_pdfCalcs']['_right_margin'],
        _pdfSettings['_pdfCalcs']['_bottom_margin'])
    canv.drawCentredString(0.5 * _pdfSettings['_pdfCalcs']['_page'][0],
        0.5 * _pdfSettings['_pdfCalcs']['_unit'],
        "Page %d" % canv.getPageNumber())

##
# _pdfCreateTitlePage() adds a title page to a ReportLab canvas object.
# @param canv An instance of a ReportLab Canvas object.
# @param _pdfSettings An options dictionary created by _pdfInitialize().
# @param reporttitle Title of report; if '', _pdfTitle is used.
# @param reportauthor Author/preparer of report.
# @retval None
def _pdfCreateTitlePage(canv, _pdfSettings, reporttitle = '', reportauthor = ''):
    """
    _pdfCreateTitlePage() adds a title page to a ReportLab canvas object.
    """
    import textwrap
    _pdfDrawPageFrame(canv, _pdfSettings)
    # _title_y is the y-coordinate at which a given line of the title should be
    # printed.  It is defined here because both the title and the author renderers
    # need to see it.
    _title_y = 7 * _pdfSettings['_pdfCalcs']['_unit']
    canv.setFont("Times-Bold", 36)
    if reporttitle == '':
        canv.drawCentredString(0.5 * _pdfSettings['_pdfCalcs']['_page'][0],
            7 * _pdfSettings['_pdfCalcs']['_unit'], _pdfSettings['_pdfTitle'])
    else:
        # This is potentially confusing, so here goes: we're going to split
        # the title on '\n's and write each resulting bit of text on its own
        # line.  drawCentredString() does not automatically do this so I had
        # to hack it on.
        _bits = string.split(reporttitle, '\n')
        for _b in _bits:
            # Based strictly on observation we need to wrap titles longer than 26
            # characters when using the default 36-point typeface.
            if len(_b) > 26:
                _b_wrapped = textwrap.wrap(_b, 26, break_long_words=True)
#                 print(_b_wrapped)
                for _bw in _b_wrapped:
                    canv.drawCentredString( 0.5 * _pdfSettings['_pdfCalcs']['_page'][0], \
                        _title_y, _bw)
                    _title_y = _title_y - 1 * _pdfSettings['_pdfCalcs']['_unit']
            else:
                canv.drawCentredString( 0.5 * _pdfSettings['_pdfCalcs']['_page'][0], _title_y, _b)
                _title_y = _title_y - 1 * _pdfSettings['_pdfCalcs']['_unit']
    # Only print the author's name if there is one.  A report is not required to have an
    # author.
    if reportauthor != '':
        canv.setFont("Times-Bold", 18)
        canv.drawCentredString(0.5 * _pdfSettings['_pdfCalcs']['_page'][0],
            _title_y, reportauthor)
    canv.showPage()
