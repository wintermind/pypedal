###############################################################################
# NAME: pyp_reports_templates.py
# VERSION: 2.0.0 (29SEPTEMBER2010)
# AUTHOR: John B. Cole, PhD (john.cole@ars.usda.gov)
# LICENSE: LGPL
###############################################################################
# FUNCTIONS:
###############################################################################

# Consider the following code sequence which provides a very simple "hello world" example
# for Platypus.

# First we import some constructors, some paragraph styles and other conveniences from
# other modules.
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.rl_config import defaultPageSize
from reportlab.lib.units import inch

styles = getSampleStyleSheet()

##
# We define the fixed features of the first page of the document with this function.
def myFirstPage(_pdfSettings, canvas, doc):
    """
    We define the fixed features of the first page of the document with this function.
    """
    canvas.saveState()
    canvas.setFont('Times-Bold',16)
    canvas.drawCentredString(_pdfSettings['_pdfCalcs']['_page_width']/2.0, _pdfSettings['_pdfCalcs']['_page_height']-108, _pdfTitle)
    canvas.setFont('Times-Roman',9)
    canvas.drawString(inch, 0.75 * inch, "First Page / %s" % _pdfSettings['_pdfPageinfo'])
    canvas.restoreState()

##
# Since we want pages after the first to look different from the first we define an
# alternate layout for the fixed features of the other pages. Note that the two functions
# above use the pdfgen level canvas operations to paint the annotations for the pages.
def myLaterPages(_pdfSettings, canvas, doc):
    """
    Since we want pages after the first to look different from the first we define an
    alternate layout for the fixed features of the other pages. Note that the two functions
    above use the pdfgen level canvas operations to paint the annotations for the pages.
    """
    canvas.saveState()
    canvas.setFont('Times-Roman', 9)
    canvas.drawString(inch, 0.75 * inch, "Page %d %s" % (doc.page, _pdfSettings['_pdfPageinfo']))
    canvas.restoreState()

def go(_pdfSettings):
      doc = SimpleDocTemplate("phello.pdf")
      print 'Writing PDF to %s' % ( "phello.pdf" )
      Story = [Spacer(1,2*inch)]
      style = styles["Normal"]
      for i in range(100):
          bogustext = ("Paragraph number %s. " % i) *20
          p = Paragraph(bogustext, style)
          Story.append(p)
          Story.append(Spacer(1,0.2*inch))
      doc.build(Story, onFirstPage=myFirstPage(_pdfSettings),
                    onLaterPages=myLaterPages(_pdfSettings))
