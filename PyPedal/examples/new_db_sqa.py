#!/usr/bin/python
from __future__ import print_function
###############################################################################
# NAME: new_db_sqa.py
# VERSION: 2.0.0b15 (23MAY2006)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_newclasses
from PyPedal import pyp_nrm
from PyPedal import pyp_db
from PyPedal import pyp_reports

if __name__ == '__main__':

    example = pyp_newclasses.loadPedigree(optionsfile='new_db.ini')

    db = pyp_db.getCursorSQA(dbdebug=example.kw['database_debug'])
    if db:
        print('Successfully connected to the database %s.' % ( example.kw['database_name'] ))
    if pyp_db.tableExistsSQA(db,example.kw['dbtable_name']):
        print('The table %s exists in %s' % ( example.kw['dbtable_name'], example.kw['database_name'] ))
    if not pyp_db.tableExistsSQA(db,'cheesegrommit'):
        print('The table cheesegrommit does not exist in %s' % ( example.kw['database_name'] ))

#    pyp_db.tableDropTable(example.kw['database_name'], example.kw['dbtable_name'])
#    pyp_db.loadPedigreeTable(example)
#    mean_inbreeding = pyp_reports.meanMetricBy(example,metric='fa',byvar='by')
#    print(mean_inbreeding)
