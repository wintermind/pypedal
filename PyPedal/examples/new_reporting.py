#!/usr/bin/python
from __future__ import print_function
###############################################################################
# NAME: new_reporting.py
# VERSION: 2.0.0b7 (11APRIL2006)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

from PyPedal import pyp_graphics
from PyPedal import pyp_network
from PyPedal import pyp_newclasses
from PyPedal import pyp_nrm
#from PyPedal import pyp_reports
from PyPedal import pyp_metrics

if __name__ == '__main__':

    example = pyp_newclasses.loadPedigree(optionsfile='new_reporting.ini')
    #example.printoptions()

    #pyp_reports.pdfPedigreeMetadata(example, titlepage = 1, reporttitle = 'Metadada for My Pedigree', reportauthor = 'John B. Cole, PhD')

    #pyp_reports.pdfPedigreeMetadata(example, titlepage = 1, reportauthor = 'John B. Cole, PhD', reportfile = 'metadata_report.pdf')

    ib = pyp_nrm.inbreeding(example)
    #coi_by_year = pyp_reports.meanMetricBy(example,metric='fa',byvar='by',createpdf=1)
    #coi_by_year

    #print(ib)
    #for _e in example.pedigree:
        #print(_e.name, _e.fa)

    print(example.kw['paper_size'])

    #print(example.backmap)
    ## Use with new_renumbering.ped.
    #pyp_reports.pdf3GenPed([56,72], example)
    ## Use with horse.ped.
    #pyp_reports.pdf3GenPed(["Pie's Joseph","Green's Dingo"], example)
    #pyp_reports.pdf3GenPed("Green's Dingo", example,reportfile='greens_dingo_pedigree.pdf')
    #pyp_reports.pdf3GenPed(example.namemap.keys(), example)

    pyp_graphics.draw_pedigree(example, gfilename='greens_dingo_pedigree', \
		gtitle="Green's Dingo pedigree", gname=1, gformat='ps', garrow=1)

    #matings = {}
    #for s in example.metadata.unique_sire_list:
    #    for d in example.metadata.unique_dam_list:
    #        matings[example.pedigree[s-1].name] = example.pedigree[d-1].name
    #pyp_metrics.mating_coi_group(matings,example,names=1)

    pg = pyp_network.ped_to_graph(example)
    #print(pg, pg.degree(), pg.nodes())
    census = pyp_network.dyad_census(pg,debug=1)
    print('max dyads:\t', ( ( pg.order()*(pg.order()-1) ) / 2 ))
    print('census:\t\t', census)
