from PyPedal import *
import copy

def main():
	options1 = {
		'pedname': 'Fake Pedigree 1',
		'messages': 'verbose',
		'renumber': 1,
		'pedfile': 'set1.ped',
		'pedformat': 'asd',
		'debug_messages': True,
	}

	options2 = copy.copy(options1)
	options2['pedname'] = 'Fake Pedigree 2'
	options2['pedfile'] = 'set2.ped'

	set1 = pyp_newclasses.loadPedigree(options1, debugLoad=True)
	print 'Animals in set1.ped:'
	print set1.idmap.keys()

	set2 = pyp_newclasses.loadPedigree(options2, debugLoad=True)
	print 'Animals in set2.ped:'
 	print set2.idmap.keys()

	print 'Testing the "+" operator...'
	added = set1 + set2
	print added.idmap.keys()

        print '='*80

	options3 = copy.copy(options1)
	options3['pedname'] = 'Fake Pedigree 3'
	options3['pedfile'] = 'set3.ped'

        set3 = pyp_newclasses.loadPedigree(options3, debugLoad=True)
	print 'Animals in set3.ped:'
	print set3.idmap.keys()

        print 'Testing the "+" operator...'
	added2 = set1 + set3
	print added2.idmap.keys()

if __name__ == '__main__':
	main()
