#!/usr/bin/python
from __future__ import print_function
###############################################################################
# NAME: pyp_classes.py
# VERSION: 2.0.0a10 (14APRIL2005)
# AUTHOR: John B. Cole, PhD (jcole@aipl.arsusda.gov)
# LICENSE: LGPL
###############################################################################

##
# pyp_classes contains two base classes that are used by PyPedal, the Animal() class
# and the Pedigree() class.  What most PyPedal routines recognize as a pedigree is
# actually just a Python list of Animal() objects.  An instance of a Pedigree() object
# is a collection of METADATA about a list of Animals().  I know that this is confusing,
# and it is going to change by the time that PyPedal 2.0.0 final is released.
##

import string
from time import *
from numpy import *

# This is my shot at supporting translations.
# The following three lines work but I need to do some work on paths to make sure
# that they will always work correctly.
#import gettext
#lang = gettext.translation('pypedal','/work1/jcole/pypedal/PyPedal/locale',languages=['en'])
#lang.install()

# This will move somewhere else later, but for new set the test munging fprmat here.
PYPEDAL_OUTPUT_TYPE = 'text' # 'html' or 'html'

##
# The Animal() class is holds animals records read from a pedigree file.
class Animal:
    """A simple class to hold animals records read from a pedigree file."""
    ##
    # __init__() initializes an Animal() object.
    # @param self Reference to the current Animal() object
    # @param animalID Animal ID number
    # @param sireID Sire ID number
    # @param damID Dam ID number
    # @param gen Generation to which the animal belongs
    # @param by Birthyear of the animal
    # @param sex Sex of the animal (m|f|u)
    # @param fa Coefficient of inbreeding of the animal
    # @param name Name of animal
    # @param alleles A two-element array of strings, which represent allelotypes.
    # @param breed Breed of animal
    # @param age Age of animal
    # @param alive Status of animal (alive or dead)
    # @return An instance of an Animal() object populated with data
    # @defreturn object
    def __init__(self,animalID,sireID,damID,gen='0',by=1900,sex='u',fa=0.,name='u',alleles=['',''],breed='u',age=-999,alive=-999):
        """Initialize an animal record."""
        self.animalID = string.strip(animalID)
        self.renumberedID = -999
        self.sireID = string.strip(sireID)
        self.damID = string.strip(damID)
        self.gen = gen
        self.igen = -999
        self.sex = sex
        self.by = int(by)
        self.fa = fa
        if self.sireID == '0' and self.damID == '0':
            self.founder = 'y'
        else:
            self.founder = 'n'
        # Take them out, put them back. Take them out, put them back...
        self.sons = {}    
        self.daus = {}    
        self.unks = {}    
        self.paddedID = self.pad_id()
        self.ancestor = 0
        self.name = name
        self.breed = breed
        # Assign alleles for use in gene-dropping runs.  Automatically assign two distinct alleles
        # to founders if no genotypes information is provided.  Otherwise, use the alleles passed
        # to Animal.__init__().
        if alleles == ['',''] and self.founder == 'y':
            _allele_1 = '%s%s' % (self.paddedID,'__1')
            _allele_2 = '%s%s' % (self.paddedID,'__2')
            self.alleles = [_allele_1,_allele_2]
        else:
            self.alleles = alleles
        self.pedcomp = -999.9
        self.age = age
        self.alive = alive
    ##
    # printme() prints a summary of the data stored in the Animal() object.
    # @param self Reference to the current Animal() object
    def printme(self):
        """Print the contents of an animal record - used for debugging."""
        self.animalID = int(self.animalID)
        self.sireID = int(self.sireID)
        self.damID = int(self.damID)
        self.by = int(self.by)
        print('ANIMAL %s RECORD' % (self.animalID))
        print('\tAnimal ID:\t%s' % (self.animalID))
        print('\tAnimal name:\t%s' % (self.name))
        print('\tSire ID:\t%s' % (self.sireID))
        print('\tDam ID:\t\t%s' % (self.damID))
        print('\tGeneration:\t%s' % (self.gen))
        print('\tInferred gen.:\t%s' % (self.igen))
        print('\tBirth Year:\t%s' % (self.by))
        print('\tSex:\t\t%s' % (self.sex))
        print('\tCoI (f_a):\t%s' % (self.fa))
        print('\tFounder:\t%s' % (self.founder))
        print('\tAncestor:\t%s' % (self.ancestor))
        print('\tAlleles:\t%s' % (self.alleles))
        print('\tRenumbered ID:\t%s' % (self.renumberedID))
        print('\tPedigree Comp.:\t%s' % (self.pedcomp))
        print('\tBreed:\t%s' % (self.breed))
        print('\tAge:\t%s' % (self.age))
        print('\tAlive:\t%s' % (self.alive))
    ##
    # stringme() returns a summary of the data stored in the Animal() object
    # as a string.
    # @param self Reference to the current Animal() object
    def stringme(self):
        """Return the contents of an animal record as a string."""
        self.animalID = int(self.animalID)
        self.sireID = int(self.sireID)
        self.damID = int(self.damID)
        self.by = int(self.by)
        _me = ''
#         _me = '%s%s' % ( _me, _('ANIMAL %s RECORD\n') % (self.animalID) )
#         _me = '%s%s' % ( _me, _('\tAnimal ID:\t%s\n') % (self.animalID) )
#         _me = '%s%s' % ( _me, _('\tAnimal name:\t%s\n') % (self.name) )
#         _me = '%s%s' % ( _me, _('\tSire ID:\t%s\n') % (self.sireID) )
#         _me = '%s%s' % ( _me, _('\tDam ID:\t\t%s\n') % (self.damID) )
#         _me = '%s%s' % ( _me, _('\tGeneration:\t%s\n') % (self.gen) )
#         _me = '%s%s' % ( _me, _('\tInferred gen.:\t%s\n') % (self.igen) )
#         _me = '%s%s' % ( _me, _('\tBirth Year:\t%s\n') % (self.by) )
#         _me = '%s%s' % ( _me, _('\tSex:\t\t%s\n') % (self.sex) )
#         _me = '%s%s' % ( _me, _('\tCoI (f_a):\t%s\n') % (self.fa) )
#         _me = '%s%s' % ( _me, _('\tFounder:\t%s\n') % (self.founder) )
#         _me = '%s%s' % ( _me, _('\tAncestor:\t%s\n') % (self.ancestor) )
#         _me = '%s%s' % ( _me, _('\tAlleles:\t%s\n') % (self.alleles) )
#         _me = '%s%s' % ( _me, _('\tRenumbered ID:\t%s\n') % (self.renumberedID) )
#         _me = '%s%s' % ( _me, _('\tPedigree Comp.:\t%s\n') % (self.pedcomp) )
#         _me = '%s%s' % ( _me, _('\tBreed:\t%s') % (self.breed) )
#         _me = '%s%s' % ( _me, _('\tAge:\t%s') % (self.age) )
#         _me = '%s%s' % ( _me, _('\tAlive:\t%s') % (self.alive) )
        _me = '%s%s' % ( _me, 'ANIMAL %s RECORD\n' % (self.animalID) )
        _me = '%s%s' % ( _me, '\tAnimal ID:\t%s\n' % (self.animalID) )
        _me = '%s%s' % ( _me, '\tAnimal name:\t%s\n' % (self.name) )
        _me = '%s%s' % ( _me, '\tSire ID:\t%s\n' % (self.sireID) )
        _me = '%s%s' % ( _me, '\tDam ID:\t\t%s\n' % (self.damID) )
        _me = '%s%s' % ( _me, '\tGeneration:\t%s\n' % (self.gen) )
        _me = '%s%s' % ( _me, '\tInferred gen.:\t%s\n' % (self.igen) )
        _me = '%s%s' % ( _me, '\tBirth Year:\t%s\n' % (self.by) )
        _me = '%s%s' % ( _me, '\tSex:\t\t%s\n' % (self.sex) )
        _me = '%s%s' % ( _me, '\tCoI (f_a):\t%s\n' % (self.fa) )
        _me = '%s%s' % ( _me, '\tFounder:\t%s\n' % (self.founder) )
        _me = '%s%s' % ( _me, '\tAncestor:\t%s\n' % (self.ancestor) )
        _me = '%s%s' % ( _me, '\tAlleles:\t%s\n' % (self.alleles) )
        _me = '%s%s' % ( _me, '\tRenumbered ID:\t%s\n' % (self.renumberedID) )
        _me = '%s%s' % ( _me, '\tPedigree Comp.:\t%s\n' % (self.pedcomp) )
        _me = '%s%s' % ( _me, '\tBreed:\t%s' % (self.breed) )
        _me = '%s%s' % ( _me, '\tAge:\t%s' % (self.age) )
        _me = '%s%s' % ( _me, '\tAlive:\t%s' % (self.alive) )

        return _me
    ##
    # trap() checks for common errors in Animal() objects
    # @param self Reference to the current Animal() object
    def trap(self):
        """Trap common errors in pedigree file entries."""
        if int(self.animalID) == int(self.sireID):
            print('[ERROR]: Animal %s has an ID number equal to its sire\'s ID (sire ID %s).\n' % (self.animalID,self.sireID))
        if int(self.animalID) == int(self.damID):
            print('[ERROR]: Animal %s has an ID number equal to its dam\'s ID (dam ID %s).\n' % (self.animalID,self.damID))
        if int(self.animalID) < int(self.sireID):
            print('[ERROR]: Animal %s is older than its sire (sire ID %s).\n' % (self.animalID,self.sireID))
        if int(self.animalID) < int(self.damID):
            print('[ERROR]: Animal %s is older than its dam (dam ID %s).\n' % (self.animalID,self.damID))

    ##
    # pad_id() takes an Animal ID, pads it to fifteen digits, and prepends the birthyear
    # (or 1950 if the birth year is unknown).  The order of elements is: birthyear, animalID,
    # count of zeros, zeros.
    # @param self Reference to the current Animal() object
    # @return A padded ID number that is supposed to be unique across animals
    # @defreturn integer
    def pad_id(self):
        """Take an Animal ID, pad it to fifteen digits, and prepend the birthyear (or 1950 if the birth year is unknown)"""
        l = len(self.animalID)
        pl = 15 - l - 1
        if pl > 0:
            zs = '0'*pl
            pid = '%s%s%s%s' % (self.by,zs,self.animalID,l)
        else:
            pid = '%s%s%s' % (self.by,self.animalID,l)
        return pid

##
# The Pedigree() class stores metadata about pedigrees.  Hopefully this will help improve performance in some procedures,
# as well as provide some useful summary data.
class Pedigree:
    """A class to hold metadata about pedigrees.  Hopefully this will help improve performance in some procedures, as well as
    provide some useful summary data."""
    ##
    # __init__() initializes a Pedigree metata object.
    # @param self Reference to the current Pedigree() object
    # @param myped A PyPedal pedigree
    # @param inputfile The name of the file from which the pedigree was loaded
    # @param name The name assigned to the PyPedal pedigree
    # @param pedcode The format code for the PyPedal pedigree
    # @param reord Flag indicating whether or not the pedigree is reordered (0|1)
    # @param renum Flag indicating whether or not the pedigree is renumbered (0|1)
    # @return An instance of a Pedigree() object populated with data
    # @defreturn object
    def __init__(self,myped,inputfile,name,pedcode='asd',reord=0,renum=0,debug=0):
        """Initialize a pedigree record."""
       	if debug == 1:
            print('\t\t[DEBUG]:  Instantiating a new Pedigree() object...')
        if debug == 1:
            print('\t\t[DEBUG]:  Naming the Pedigree()...')
            self.name = name
        if debug == 1:
            print('\t\t[DEBUG]:  Assigning a filename...')
            self.filename = inputfile
        if debug == 1:
            print('\t\t[DEBUG]:  Attaching a pedigree...')
            self.myped = myped
        if debug == 1:
            print('\t\t[DEBUG]:  Setting the pedcode...')
            self.pedcode = pedcode
        if debug == 1:
            print('\t\t[DEBUG]:  Counting the number of animals in the pedigree...')
            self.num_records = len(self.myped)
        if debug == 1:
            print('\t\t[DEBUG]:  Counting and finding unique sires...')
            self.num_unique_sires, self.unique_sire_list = self.nus()
        if debug == 1:
            print('\t\t[DEBUG]:  Counting and finding unique dams...')
            self.num_unique_dams, self.unique_dam_list = self.nud()
        if debug == 1:
            print('\t\t[DEBUG]:  Setting reord flag...')
            self.reordered = reord
        if debug == 1:
            print('\t\t[DEBUG]:  Setting renum flag...')
            self.renumbered = renum
        if debug == 1:
            print('\t\t[DEBUG]:  Counting and finding unique generations...')
            self.num_unique_gens, self.unique_gen_list = self.nug()
        if debug == 1:
            print('\t\t[DEBUG]:  Counting and finding unique birthyears...')
            self.num_unique_years, self.unique_year_list = self.nuy()
        if debug == 1:
            print('\t\t[DEBUG]:  Counting and finding unique founders...')
            self.num_unique_founders, self.unique_founder_list = self.nuf()
        if debug == 1:
            print('\t\t[DEBUG]:  Detaching pedigree...')
            self.myped = []
    ##
    # printme() prints a summary of the metadata stored in the Pedigree() object.
    # @param self Reference to the current Pedigree() object
    def printme(self):
        """Print the pedigree metadata."""
        print('PEDIGREE %s (%s)' % (self.name,self.filename))
        print('\tRecords:\t\t%s' % (self.num_records))
        print('\tUnique Sires:\t\t%s' % (self.num_unique_sires))
        print('\tUnique Dams:\t\t%s' % (self.num_unique_dams))
        print('\tUnique Gens:\t\t%s' % (self.num_unique_gens))
        print('\tUnique Years:\t\t%s' % (self.num_unique_years))
        print('\tUnique Founders:\t%s' % (self.num_unique_founders))
        print('\tPedigree Code:\t\t%s' % (self.pedcode))
    ##
    # stringme() returns a summary of the metadata stored in the pedigree as
    # a string.
    # @param self Reference to the current Pedigree() object
    def stringme(self):
        """Build a string from the pedigree metadata."""
        _me = ''
        _me = '%s%s' % ( _me, 'PEDIGREE %s (%s)\n' % (self.name,self.filename) )
        _me = '%s%s' % ( _me, '\tRecords:\t\t\t%s\n' % (self.num_records) )
        _me = '%s%s' % ( _me, '\tUnique Sires:\t\t%s\n' % (self.num_unique_sires) )
        _me = '%s%s' % ( _me, '\tUnique Dams:\t\t%s\n' % (self.num_unique_dams) )
        _me = '%s%s' % ( _me, '\tUnique Gens:\t\t%s\n' % (self.num_unique_gens) )
        _me = '%s%s' % ( _me, '\tUnique Years:\t\t%s\n' % (self.num_unique_years) )
        _me = '%s%s' % ( _me, '\tUnique Founders:\t%s\n' % (self.num_unique_founders) )
        _me = '%s%s' % ( _me, '\tPedigree Code:\t\t%s\n' % (self.pedcode) )
        return _me
    ##
    # fileme() writes the metada stored in the Pedigree() object to disc.
    # @param self Reference to the current Pedigree() object
    def fileme(self):
        """Save the pedigree metadata to a file."""
        outputfile = '%s%s%s' % (self.name,'_ped_metadata_','.dat')
        aout = open(outputfile,'w')
        line1 = 'PEDIGREE %s (%s)\n' % (self.name,self.filename)
        line2 = '\tRecords:\t%s\n' % (self.num_records)
        line3 = '\tUnique Sires:\t%s\n' % (self.num_unique_sires)
        line4 = '\tUnique Dams:\t%s\n' % (self.num_unique_dams)
        line5 = '\tPedigree Code:\t%s\n' % (self.pedcode)
        line6 = '\tUnique Founders:\t%s\n' % (self.num_unique_founders)
        line7 =  '\tUnique Gens:\t%s\n' % (self.num_unique_gens)
        line8 = '\tUnique Years:\t%s\n' % (self.num_unique_years)
        line9 = '='*80
        line10 = '\tUnique Sire List:\t%s\n' % (self.unique_sire_list)
        line11 = '\tUnique Dam List:\t%s\n' % (self.unique_dam_list)
        line12 = '\tUnique Gen List:\t%s\n' % (self.unique_gen_list)
        line13 = '\tUnique Year List:\t%s\n' % (self.unique_year_list)
        line14 = '\tUnique Founder List:\t%s\n' % (self.num_founder_list)
        aout.write(line1)
        aout.write(line2)
        aout.write(line3)
        aout.write(line4)
        aout.write(line5)
        aout.write(line6)
        aout.write(line7)
        aout.write(line8)
        aout.write(line9)
        aout.write(line10)
        aout.write(line11)
        aout.write(line12)
        aout.write(line13)
        aout.write(line14)
        aout.close()

    ##
    # nus() returns the number of unique sires in the pedigree along with a list of the sires
    # @param self Reference to the current Pedigree() object
    # @return The number of unique sires in the pedigree and a list of those sires
    # @defreturn integer-and-list
    def nus(self):
        """Count the number of unique sire IDs in the pedigree.  Returns an integer count and a Python list of the
        unique sire IDs."""
        siredict = {}
        for l in range(self.num_records):
            if int(self.myped[l].sireID) != 0:
                try:
                    _s = siredict[self.myped[l].sireID]
                except KeyError:
                    siredict[self.myped[l].sireID] = self.myped[l].sireID
        n = len(siredict.keys())
        return n, siredict.keys()
    ##
    # nud() returns the number of unique dams in the pedigree along with a list of the dams
    # @param self Reference to the current Pedigree() object
    # @return The number of unique dams in the pedigree and a list of those dams
    # @defreturn integer-and-list
    def nud(self):
        """Count the number of unique dam IDs in the pedigree.  Returns an integer count and a Python list of the
        unique dam IDs."""
        damdict = {}
        for l in range(self.num_records):
            if int(self.myped[l].damID) != 0:
                try:
                    _d = damdict[self.myped[l].damID]
                except KeyError:
                    damdict[self.myped[l].damID] = self.myped[l].damID
        n = len(damdict.keys())
        return n, damdict.keys()
    ##
    # nug() returns the number of unique generations in the pedigree along with a list of the generations
    # @param self Reference to the current Pedigree() object
    # @return The number of unique generations in the pedigree and a list of those generations
    # @defreturn integer-and-list
    def nug(self):
        """Count the number of unique generations in the pedigree.  Returns an integer count and a Python list of the unique generations."""
        gendict = {}
        for l in range(self.num_records):
            try:
                _g = gendict[self.myped[l].gen]
            except KeyError:
                gendict[self.myped[l].gen] = self.myped[l].gen
        n = len(gendict.keys())
        return n, gendict.keys()
    ##
    # nuy() returns the number of unique birthyears in the pedigree along with a list of the birthyears
    # @param self Reference to the current Pedigree() object
    # @return The number of unique birthyears in the pedigree and a list of those birthyears
    # @defreturn integer-and-list
    def nuy(self):
        """Count the number of unique birth years in the pedigree.  Returns an integer count and a Python list of the
        unique birth years."""
        yeardict = {}
        for l in range(self.num_records):
            try:
                _y = yeardict[self.myped[l].by]
            except KeyError:
                yeardict[self.myped[l].by] = self.myped[l].by
        n = len(yeardict.keys())
        return n, yeardict.keys()
    ##
    # nuf() returns the number of unique founders in the pedigree along with a list of the founders
    # @param self Reference to the current Pedigree() object
    # @return The number of unique founders in the pedigree and a list of those founders
    # @defreturn integer-and-list
    def nuf(self):
        """Count the number of unique founders in the pedigree."""
        founderdict = {}
        for l in range(self.num_records):
            if self.myped[l].founder == 'y':
                try:
                    _f = founderdict[self.myped[l].animalID]
                except KeyError:
                    founderdict[self.myped[l].animalID] = self.myped[l].animalID
        n = len(founderdict.keys())
        return n, founderdict.keys()
