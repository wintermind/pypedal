#!/usr/bin/python
from __future__ import print_function
###############################################################################
# NAME: new_graphics.py
# VERSION: 2.0.0b5 (13DECEMBER2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_db
from PyPedal import pyp_graphics
from PyPedal import pyp_newclasses
from PyPedal import pyp_nrm
from PyPedal import pyp_reports
from PyPedal.pyp_utils import pyp_nice_time

if __name__ == '__main__':

    example = pyp_newclasses.loadPedigree(optionsfile='new_graphics.ini')
#     example_inbreeding = pyp_nrm.inbreeding(example)

    pyp_graphics.rmuller_spy_matrix_pil(example.nrm.nrm, fname='dog_sparsity.png', cutoff=0.01, do_outline=0, height=example.nrm.nrm.shape[0], width=example.nrm.nrm.shape[0])

    pyp_graphics.rmuller_pcolor_matrix_pil(example.nrm.nrm, fname='dog_long_pcolored.png', do_outline=0, height=example.nrm.nrm.shape[0], width=example.nrm.nrm.shape[0])

    #pyp_db.loadPedigreeTable(example)
    pyp_db.deleteTable(example)
    pyp_db.populatePedigreeTable(example)
    coi_by_year = pyp_reports.meanMetricBy(example, metric='fa', byvar='by')
    print('coi_by_year: ', coi_by_year)
    cby = coi_by_year
    del(cby[1900])
    pyp_graphics.plot_line_xy(coi_by_year, gfilename='dog_coi_by_year',
        gtitle='', gxlabel='Birth year', gylabel='Coefficient of inbreeding')

    pyp_graphics.draw_pedigree(example, gfilename='dog_pedigree', gtitle='Dog pedigree',gformat='jpg',gsize='f')

    pyp_graphics.plot_founders_by_year(example, gfilename='dog_founders_by_year', gtitle='Number of founders within each birthyear')

    pyp_graphics.plot_founders_pct_by_year(example, gfilename='dog_founder_pct', gtitle='Percentage of founders within each birthyear')
