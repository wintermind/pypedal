from PyPedal import *

options = {}
options['pedfile'] = 'illhaveanother.ped'
options['pedformat'] = 'ASDyx'
options['missing_parent'] = '0'
options['renumber'] = 1

if __name__ == '__main__':

	iha = pyp_newclasses.loadPedigree(options)
