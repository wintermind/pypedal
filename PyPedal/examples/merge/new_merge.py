#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      jcole
#
# Created:     05/06/2012
# Copyright:   (c) jcole 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

from PyPedal import *

if __name__ == '__main__':

    options = {}
    options['pedname'] = 'Fake Pedigree 1'
    options['messages'] = 'verbose'
    options['renumber'] = 1
    options['pedfile'] = 'merge1.ped'
    options['pedformat'] = 'asd'
    merge1 = pyp_newclasses.loadPedigree(options)

    options2 = {}
    options2['pedname'] = 'Fake Pedigree 2'
    options2['messages'] = 'verbose'
    options2['renumber'] = 1
    options2['pedfile'] = 'merge2.ped'
    options2['pedformat'] = 'asd'
    merge2 = pyp_newclasses.loadPedigree(options2)

    print merge1.kw['pedname']
    print merge2.kw['pedname']

    merge3 = merge1 + merge2