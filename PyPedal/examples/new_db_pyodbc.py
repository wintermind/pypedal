#!/usr/bin/python

###############################################################################
# NAME: new_db.py
# VERSION: 2.0.0b5 (29DECEMBER2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_newclasses
from PyPedal import pyp_nrm
from PyPedal import pyp_db
from PyPedal import pyp_reports
from PyPedal import pyp_graphics

options = {}
options['messages'] = 'verbose'
options['pedfile'] = 'hartlandclark.ped'
options['pedname'] = 'Pedigree from van Noordwijck and Scharloo (1981)'
options['pedformat'] = 'asdb'
options['pedigree_is_renumbered'] = 1
options['database_name'] = 'new_db_sqlite3_test.db'
options['dbtable_name'] = 'test'

if __name__ == '__main__':

    example = pyp_newclasses.loadPedigree(options)
    pyp_nrm.inbreeding(example)

    print 'Dropping existing table...'
    pyp_db.deleteTable(example)
    print 'Checking to see if the table is gone...'
    pyp_db.doesTableExist(example)
    print 'Creating the table...'
    pyp_db.createPedigreeTable(example)
    print 'Populating the table...'
    pyp_db.populatePedigreeTable(example)

    mean_inbreeding = pyp_reports.meanMetricBy(example,metric='fa',byvar='by')
    print mean_inbreeding

    pyp_graphics.plot_line_xy(mean_inbreeding, gfilename='great_tit_coi_by_year', \
        gtitle='', gxlabel='Birth year', gylabel='Coefficient of inbreeding')
