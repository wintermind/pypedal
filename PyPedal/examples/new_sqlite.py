from PyPedal import pyp_newclasses
from PyPedal import pyp_metrics
from PyPedal import pyp_nrm
from PyPedal import pyp_io
from PyPedal import pyp_graphics

options = {}
# This is the name of the input tile.
options['pedfile'] = ''
# This is a descriptor used in some output.
options['pedname'] = 'Text Stream'
# Don't mess with these options.
options['messages'] = 'verbose'
options['renumber'] = 1
options['pedformat'] = 'ASD'
options['assign_sexes'] = 1

if __name__ == "__main__":

    pedstream = 'a1,s1,d1\na2,s2,d2\na3,a1,a2\n'
    test = pyp_newclasses.loadPedigree(options,pedsource='textstream',pedstream=pedstream)

    test.kw['database_name'] = 'test_pypedal_save'
    test.kw['dbtable_name'] = 'test'
    test.savedb(drop=True)

    test2 = pyp_newclasses.loadPedigree(options,pedsource='db')
    test2.metadata.printme()

    test2.savegraph(pedoutfile='test.adj')
