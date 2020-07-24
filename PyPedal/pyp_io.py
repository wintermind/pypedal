#!/usr/bin/python

###############################################################################
# NAME: pyp_io.py
# VERSION: 2.0.0 (29SEPTEMBER2010)
# AUTHOR: John B. Cole, PhD (john.cole@ars.usda.gov)
# LICENSE: LGPL
###############################################################################
# FUNCTIONS:
#   a_inverse_from_file()
#   a_inverse_to_file()
#   dissertation_pedigree_to_file()
#   dissertation_pedigree_to_pedig_format()
#   dissertation_pedigree_to_pedig_interest_format()
#   dissertation_pedigree_to_pedig_format_mask()
#   pyp_file_header()
#   pyp_file_footer()
#   renderTitle()
#   renderBodyText()
#   pickle_pedigree()
#   unpickle_pedigree()
#   summary_inbreeding()
#   save_ijk()
#   load_from_gedcom()
#   save_from_gedcom()
#   save_to_gedcom()
#   load_from_genes()
#   save_from_genes()
#   save_to_genes()
#   save_newanimals_to_file()
###############################################################################

## @package pyp_io
# pyp_io contains several procedures for writing structures to and reading them from
# disc (e.g. using pickle() to store and retrieve A and A-inverse).  It also includes a set
# of functions used to render strings as HTML or plaintext for use in generating output
# files.
##

import logging, numpy, pickle, string, time
import pyp_utils

global LINE1
global LINE2
LINE1 = '%s' % ('='*80)
LINE2 = '%s' % ('-'*80)

##
# a_inverse_to_file() uses the Python pickle system for persistent objects to write the
# inverse of a relationship matrix to a file.
# @param pedobj A PyPedal pedigree object.
# @param ainv The inverse of a numerator relationship matrix, A, or an empty string if A is to be calculated.
# @retval True (1) on success, false (0) on failure
def a_inverse_to_file(pedobj, ainv=''):
    """
    Use the Python pickle system for persistent objects to write the inverse of a relationship matrix to a file.
    """
    try: logging.info('Entered a_inverse_to_file()')
    except: pass
    try:
        from pickle import Pickler
        if not ainv:
            ainv = a_inverse_df(pedobj.pedigree,pedobj.kw['filetag'])
        a_outputfile = '%s%s%s' % (filetag,'_a_inverse_pickled_','.pkl')
        aout = open(a_outputfile,'w')
        ap = pickle.Pickler(aout)
        ap.dump(a)
        aout.close()
        _r = 1
    except:
        _r = 0

    try: logging.info('Exited a_inverse_to_file()')
    except: pass
    return _r

##
# a_inverse_from_file() uses the Python pickle system for persistent objects to read the inverse of
# a relationship matrix from a file.
# @param inputfile The name of the input file.
# @retval The inverse of a numerator relationship matrix.
def a_inverse_from_file(inputfile):
    """
    Use the Python pickle system for persistent objects to read the inverse of a relationship matrix from a file.
    """
    try: logging.info('Entered a_inverse_from_file()')
    except: pass
    try:
        from pickle import Pickler
        ain = open(inputfile,'r')
        au = pickle.Unpickler(ain)
        a_inv = au.load()
    except:
        a_inv = numpy.zeros([1,1],Float)
    try: logging.info('Exited a_inverse_from_file()')
    except: pass
    return a_inv

##
# dissertation_pedigree_to_file() takes a pedigree in 'asdxfg' format and writes is to a file.
# @param pedobj A PyPedal pedigree object.
# @retval True (1) on success, false (0) on failure
def dissertation_pedigree_to_file(pedobj):
    """
    dissertation_pedigree_to_file() takes a pedigree in 'asdxfg' format and writes is to
    a file.
    """
    # This procedure assumes that the pedigree passed to it is in 'asdxfg' format.
    try: logging.info('Entered dissertation_pedigree_to_file()')
    except: pass
    try:
        length = len(pedobj.pedigree)
        #print 'DEBUG: length of pedigree is %s' % (length)
        outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_diss','.ped')
        print '\t\tWriting dissertation pedigree to %s' % (outputfile)
        aout = open(outputfile,'w')
        aout.write('# DISSERTATION pedigree produced by PyPedal.\n')
        aout.write('% asdbxfg\n')
        for l in range(length):
            aout.write('%s,%s,%s,%s,%s,%s,%s\n' % pedobj.pedigree[l].animalID,pedobj.pedigree[l].sireID,pedobj.pedigree[l].damID,pedobj.pedigree[l].by, pedobj.pedigree[l].sex,pedobj.pedigree[l].fa,pedobj.pedigree[l].gen)
        aout.close()
        _r = 1
    except:
        _r = 0
    try: logging.info('Exited dissertation_pedigree_to_file()')
    except: pass
    return _r

##
# dissertation_pedigree_to_pedig_format() takes a pedigree in 'asdbxfg' format, formats it into
# the form used by Didier Boichard's 'pedig' suite of programs, and writes it to a file.
# @param pedobj A PyPedal pedigree object.
# @retval True (1) on success, false (0) on failure
def dissertation_pedigree_to_pedig_format(pedobj):
    """
    dissertation_pedigree_to_pedig_format() takes a pedigree in 'asdbxfg' format, formats
    it into the form used by Didier Boichard's 'pedig' suite of programs, and writes it
    to a file.
    """
    try: logging.info('Entered dissertation_pedigree_to_pedig_format()')
    except: pass
    try:
        length = len(pedobj.pedigree)
        outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_pedig','.ped')
        aout = open(outputfile,'w')
        for l in range(length):
            if pedobj.pedigree[l].sex == 'm' or pedobj.pedigree[l].sex == 'M':
                    sex = 1
            else:
                sex = 2
            aout.write('%s %s %s %s %s %s %s\n' % pedobj.pedigree[l].animalID,pedobj.pedigree[l].sireID,pedobj.pedigree[l].damID,pedobj.pedigree[l].by,sex,'1','1')
        aout.close()
        _r = 1
    except:
        _r = 0
    try: logging.info('Exited dissertation_pedigree_to_pedig_format()')
    except: pass
    return _r

##
# dissertation_pedigree_to_pedig_interest_format() takes a pedigree in 'asdbxfg' format,
# formats it into the form used by Didier Boichard's parente program for the studied
# individuals file.
# @param pedobj A PyPedal pedigree object.
# @retval True (1) on success, false (0) on failure
def dissertation_pedigree_to_pedig_interest_format(pedobj):
    """
    dissertation_pedigree_to_pedig_interest_format() takes a pedigree in 'asdbxfg' format,
    formats it into the form used by Didier Boichard's parente program for the studied
    individuals file.
    """
    try: logging.info('Entered dissertation_pedigree_to_pedig_interest_format()')
    except: pass
    try:
        length = len(pedobj.pedigree)
        outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_parente','.ped')
        aout = open(outputfile,'w')
        for l in range(length):
            aout.write('%s %s\n' % pedobj.pedigree[l].animalID,'1')
        aout.close()
        _r = 1
    except:
        _r = 0
    try: logging.info('Exited dissertation_pedigree_to_pedig_interest_format()')
    except: pass
    return _r

##
# dissertation_pedigree_to_pedig_format_mask() Takes a pedigree in 'asdbxfg' format,
# formats it into the form used by Didier Boichard's 'pedig' suite of programs, and
# writes it to a file. THIS FUNCTION MASKS THE GENERATION ID WITH A FAKE BIRTH YEAR
# AND WRITES THE FAKE BIRTH YEAR TO THE FILE INSTEAD OF THE TRUE BIRTH YEAR. THIS IS
# AN ATTEMPT TO FOOL PEDIG TO GET f_e, f_a et al. BY GENERATION.
# @param pedobj A PyPedal pedigree object.
# @retval True (1) on success, false (0) on failure
def dissertation_pedigree_to_pedig_format_mask(pedobj):
    """
    dissertation_pedigree_to_pedig_format_mask() Takes a pedigree in 'asdbxfg' format,
    formats it into the form used by Didier Boichard's 'pedig' suite of programs, and
    writes it to a file. THIS FUNCTION MASKS THE GENERATION ID WITH A FAKE BIRTH YEAR
    AND WRITES THE FAKE BIRTH YEAR TO THE FILE INSTEAD OF THE TRUE BIRTH YEAR. THIS IS
    AN ATTEMPT TO FOOL PEDIG TO GET f_e, f_a et al. BY GENERATION.
    """
    try: logging.info('Entered dissertation_pedigree_to_pedig_format_mask()')
    except: pass
    try:
        length = len(pedobj.pedigree)
        outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_pedig_mask','.ped')
        aout = open(outputfile,'w')
        for l in range(length):
            ## mask generations (yes, this could be shorter - but this is easy to debug
            mygen = float(pedobj.pedigree[l].gen)
            if ( mygen > 0 and mygen <= 1.25 ):         _gen = 10
            elif ( mygen > 1.25 and mygen <= 1.75 ):    _gen = 15
            elif ( mygen > 1.75 and mygen <= 2.25 ):    _gen = 20
            elif ( mygen > 2.25 and mygen <= 2.75 ):    _gen = 25
            elif ( mygen > 2.75 and mygen <= 3.25 ):    _gen = 30
            elif ( mygen > 3.25 and mygen <= 3.75 ):    _gen = 35
            elif ( mygen > 3.75 and mygen <= 4.25 ):    _gen = 40
            elif ( mygen > 4.25 and mygen <= 4.75 ):    _gen = 45
            elif ( mygen > 4.75 and mygen <= 5.25 ):    _gen = 50
            elif ( mygen > 5.25 and mygen <= 5.75 ):    _gen = 55
            elif ( mygen > 5.75 and mygen <= 6.25 ):    _gen = 60
            elif ( mygen > 6.25 and mygen <= 6.75 ):    _gen = 65
            elif ( mygen > 6.75 and mygen <= 7.25 ):    _gen = 70
            elif ( mygen > 7.25 and mygen <= 7.75 ):    _gen = 75
            else:                                       _gen = 0
            _maskgen = 1950 + _gen
            ## convert sexes
            if pedobj.pedigree[l].sex == 'm' or pedobj.pedigree[l].sex == 'M':
                sex = 1
            else:
                sex = 2
            aout.write('%s %s %s %s %s %s %s\n' % pedobj.pedigree[l].animalID,pedobj.pedigree[l].sireID,pedobj.pedigree[l].damID,_maskgen,sex,'1','1')
        aout.close()
        _r = 1
    except:
        _r = 0
    try: logging.info('Exited dissertation_pedigree_to_pedig_format_mask()')
    except: pass
    return _r

##
# pyp_file_header() writes a header to a page of PyPedal output.
# @param ofhandle A Python file handle.
# @param caller A string indicating the name of the calling routine.
# @retval None
def pyp_file_header(ofhandle, caller="Unknown PyPedal routine"):
    """
    pyp_file_header() writes a header to a page of PyPedal output.
    """
    try:
        ofhandle.write('%s\n' % ('-'*80))
        ofhandle.write('Created by %s at %s\n' % (caller,pyp_utils.pyp_nice_time()))
        ofhandle.write('%s\n' % ('-'*80))
    except:
        pass

##
# pyp_file_footer() writes a footer to a page of PyPedal output.
# @param ofhandle A Python file handle.
# @param caller A string indicating the name of the calling routine.
# @retval None
def pyp_file_footer(ofhandle, caller="Unknown PyPedal routine"):
    """
    pyp_file_footer() writes a footer to a page of PyPedal output.
    """
    try:
        ofhandle.write('%s\n' % ('-'*80))
    except:
        pass

##
# renderTitle() renders page titles (produces HTML output by default).
# @param title_string String to be enclosed in HTML "H" tags
# @param title_level Size to be attached to "H" tags, e.g., "H1"
# @retval None
def renderTitle(title_string, title_level="1"):
    """
    renderTitle() renders page titles (produces HTML output by default).
    """
    if not PYPEDAL_OUTPUT_TYPE:
        PYPEDAL_OUTPUT_TYPE = 'html'
    # We are trying to keep it simple here.
    if ( not title_level ) or ( title_level < 1 ) or ( title_level > 3 ):
        title_level = 1
    if PYPEDAL_OUTPUT_TYPE == 'h':
        renderedTitle = '<H%s>%s</H%s>\n' % (title_level,title_string,title_level)
    else:
        _underline = '='*len(title_string)
        renderedTitle = '%s\n' % (title_string,_underline)
    return renderedTitle

##
# renderBodyText() renders page contents (produces HTML output by default).
# @param text_string String to be rendered with either a trailing newline or enclosed in HTNL "P" tags
# @retval None
def renderBodyText(text_string):
    """
    renderBodyText() renders page contents (produces HTML output by default).
    """
    if not PYPEDAL_OUTPUT_TYPE:
        PYPEDAL_OUTPUT_TYPE = 'html'
    if PYPEDAL_OUTPUT_TYPE == 'h':
        renderedBodyText = '<p>%s</p>' % (text_string)
    else:
        renderedBodyText = '%s\n' % (text_string)
    return renderedBodyText

##
# pickle_pedigree() pickles a pedigree.
# @param pedobj An instance of a PyPedal pedigree object.
# @param filename The name of the file to which the pedigree object should be pickled (optional).
# @retval A 1 on success, a 0 otherwise.
def pickle_pedigree(pedobj, filename=''):
    """
    pickle_pedigree() pickles a pedigree.
    """
    try: logging.info('Entered pickle_pedigree()')
    except: pass
    try:
        _r = 1
        if not filename:
          _pfn = '%s.pkl' % ( pedobj.kw['filetag'] )
        else:
            _pfn = '%s.pkl' % ( filename )
        print _pfn
        pickle.dump(pedobj,open(_pfn,'w'))
        logging.info('Pickled pedigree %s to file %s', pedobj.kw['pedname'], _pfn )
        if pedobj.kw['messages'] == 'verbose':
            print 'Pickled pedigree %s to file %s' % ( pedobj.kw['pedname'], _pfn )
    except:
        logging.error('Unable to pickle pedigree %s to file %s', self.kw['pedname'], _pfn )
        _r = 0
    try: logging.info('Exited pickle_pedigree()')
    except: pass
    return _r

##
# unpickle_pedigree() reads a pickled pedigree in from a file and returns the unpacked
# pedigree object.
# @param filename The name of the pickle file.
# @retval An instance of a NewPedigree object on success, a 0 otherwise.
def unpickle_pedigree(filename=''):
    """
    unpickle_pedigree() reads a pickled pedigree in from a file and returns the unpacked
    pedigree object.
    """
    try: logging.info('Entered unpickle_pedigree()')
    except: pass
    try:
        if not filename:
            logging.error('No filename provided for pedigree unpickling!' )
            _r = 0
        else:
            _ck_pfn = string.split(filename,'.')
            if len(_ck_pfn) == 2:
                _pfn = filename
            else:
                _pfn = '%s.pkl' % ( filename )
                logging.info('No file extension provided for %s.  An extension (.pkl) was added.', filename )
            _myped = pickle.load(open(_pfn))
            logging.info('Unpickled pedigree %s from file %s', _myped.kw['pedname'], _pfn )
            _r = _myped
    except:
        logging.error('Unable to unpickle pedigree from file!' )
        _r = 0
    try: logging.info('Exited unpickle_pedigree()')
    except: pass
    return _r

##
# summary_inbreeding() returns a string representation of the data contained in
# the 'metadata' dictionary contained in the output dictionary returned by
# pyp_nrm/pyp_inbreeding().
# @param f_metadata Dictionary of inbreeding metadata.
# @retval A string on success, a 0 otherwise.
def summary_inbreeding(f_metadata):
    """
    summary_inbreeding() returns a string representation of the data contained in
    the 'metadata' dictionary contained in the output dictionary returned by
    pyp_nrm/pyp_inbreeding().
    """
    try:
        _summary = ''
        _summary = '%s' % (LINE1)
        _summary = '%s\n%s' % (_summary, 'Inbreeding Statistics')
        _summary = '\n%s\n%s' % (_summary, LINE1)
        _summary = '%s\n%s' % (_summary, 'All animals:')
        _summary = '\n%s\n%s' % (_summary, LINE2)
        for k,v in f_metadata['all'].iteritems():
            _line = '\t%s\t%s' % (k,v)
            _summary = '%s\n%s' % (_summary, _line)
        _summary = '\n%s\n%s' % (_summary, LINE1)
        _summary = '%s\n%s' % (_summary, 'Animals with non-zero CoI:')
        _summary = '\n%s\n%s' % (_summary, LINE2)
        for k,v in f_metadata['nonzero'].iteritems():
            _line = '\t%s\t%s' % (k,v)
            _summary = '%s\n%s' % (_summary, _line)
        _summary = '\n%s\n%s' % (_summary, LINE1)
        return _summary
    except:
        return '0'

##
# save_ijk() saves an NRM to a file in the form "animal A" "animal B" "rAB".
# @param pedobj: The pedigree to which the NRM is attached
# @param nrm_filename: The file to which the matrix should be written.
# @retval A save status indicator (0: failed, 1: success).
def save_ijk(pedobj, nrm_filename):
    """
    save_ijk() saves an NRM to a file in the form "animal A" "animal B" "rAB".
    """
    if pedobj.kw['messages'] == 'verbose':
        print '[INFO]: Saving A-matrix to file %s at %s.' % ( nrm_filename, pyp_utils.pyp_nice_time() )
    logging.info('Saving A-matrix to file %s', nrm_filename)
    #try:
    of = file(nrm_filename,"w")
    for i in range(pedobj.metadata.num_records):
        for j in range(i,pedobj.metadata.num_records):
            line = '%s %s %s\n' % ( pedobj.backmap[i+1], pedobj.backmap[j+1], pedobj.nrm.nrm[i,j])
            of.write(line)
    of.close()
    if pedobj.kw['messages'] == 'verbose':
        print '[INFO]: A-matrix successfully saved to file %s at %s.' % ( nrm_filename, pyp_utils.pyp_nice_time() )
    logging.info('A-matrix successfully saved to file %s', nrm_filename)
    return 1

##
# load_from_gedcom() reads and parses pedigree data that conform to
# a subset of the GEDCOM 5.5 specification. Not all valid GEDCOM
# are supported; unsupported tags are ignored.
# @param infilename The file to which the matrix should be written.
# @param messages Controls output to the screen
# @param standalone Uses logging if called by a NewPedigree method
# @param missing_sex Value assigned to an animal with unknown sex
# @param missing_parent Value assigned to unknown parents
# @param missing_name Name assigned by default
# @param missing_byear VAlue assigned to unknown birth years
# @param debug Flag turning debugging messages on (1) and off (0)
# @retval A save status indicator (0: failed, 1: success).
def load_from_gedcom(infilename, messages='verbose', standalone=1, missing_sex='u', \
    missing_parent=0, missing_name='Unknown Name', missing_byear='0001', debug=0):
    """
    load_from_gedcom() reads and parses pedigree data that conform to
    a subset of the GEDCOM 5.5 specification. Not all valid GEDCOM
    are supported; unsupported tags are ignored.
    """

    # NOTE: There is a lot of error-checking that could be done here but isn't. That
    # behavior is a deliberate design decision. load_from_gedcom() is NOT intended to
    # be a general-purpose GEDCOM parser, but rather is supposed to convert a GEDCOM
    # file into a pedigree format that PyPedal can load into a NewPedigree object with
    # no trouble. The error-checking and related operations, such as automatically adding
    # pedigree entries for parents with no records themselves, are left to the preprocess()
    # and load() methods of the NewPedigree and NewAnimal classes.
    #
    # If someone wants to take this function and extend it to be more general that's fine,
    # but I ask that they create a new function that extends this function rather than
    # altering the no-doubt imperfect behavior of load_from_gedcom().
    known_tags = ['BIRT','CHIL','DATE','FAM','FAMC','FAMS','HUSB', \
        'INDI','NAME','SEX','WIFE']
    current_level, any_names, any_sexes, any_birth, get_birth = 0, 0, 0, 0, 0
    indi, fam = {}, {}
    # Structures for mapping data
    fam2husb, fam2wife, fam2chil, indi2name, indi2sex, indi2famc, indi2fams, indi2birth = {}, {}, {}, {}, {}, {}, {}, {}
    # Read the file
    if standalone == 0: logging.info('[load_from_gedcom]: Opening GEDCOM pedigree file %s.',infilename)
    infile = open(infilename,'r')
    inlines = infile.readlines()
    infile.close()
    all_done = 0  # Have we processed all of the input lines?
    next_zero = 0 # Have we seen the next zero?
    zero_mark = 0 # If so, where?
    last_round = 0

    try:
        while not all_done:
            next_zero = 0
            _curr_rec = []
            for l in xrange(zero_mark,len(inlines)+1):
                # I don't think that all GEDCOM 5.5 0-level records contain three tokens, but
                # the INDI and FAM records do, and those are the records we want.
                try:
                    linelist = inlines[l].strip().split(' ')
                    if len(linelist) > 1:
                        if linelist[0] == '0' and l > zero_mark:
                            next_zero = 1
                            zero_mark = l
                            break
                        else:
                            _curr_rec.append(inlines[l])
                except IndexError:
                    all_done = 1
            if messages == 'debug':
                print '-'*80
                print _curr_rec
                print '-'*80
            firstlinelist = _curr_rec[0].strip().rstrip('\r').rstrip('\n').split(' ')
            if firstlinelist[-1] == 'INDI':
                _tag, _name, _sex, _famc, _fams, _byr = 0, 0, 0, 0, [], 0
                any_sexes, any_birth, get_birth = 0, 0, 0
                for c in xrange(1,len(_curr_rec)):
                    linelist = _curr_rec[c].strip().rstrip('\r').rstrip('\n').split(' ')
                    if c == 1:
                        _indi = _curr_rec[0].strip().split(' ')[1][1:-1]
                    if len(linelist) > 1:
                        _tag = linelist[1].upper()
                    else: break
                    if _tag in known_tags:
                        if messages == 'debug':
                            print 'Processing tag %s' % ( _tag )
                        if _tag.upper() == 'NAME':
                            _name =  ' '.join(linelist[2:])
                        elif _tag.upper() == 'SEX':
                            _sex =  linelist[2]
                            any_sexes = 1
                        elif _tag.upper() == 'FAMC':
                            _famc = linelist[2][1:-1]
                        elif _tag.upper() == 'FAMS':
                            _fams.append(linelist[2][1:-1])
                        # The only dates handled right now are birth dates
                        elif _tag.upper() == 'BIRT':
                            get_birth = 1
                        elif _tag.upper() == 'DATE' and get_birth:
                            _byr = ' '.join(linelist[2:])[-4:]
                            get_birth = 0
                            any_birth = 1
                        else: pass
                    else:
                        if messages == 'debug':
                            print 'Skipping unknown tag %s' % ( _tag )
                    # Once we've read the entire record we can process the data.
                    # Here's the basic idea: put the details from the INDI and
                    # FAM records into dictionaries (lookup tables) as we sweep
                    # the input data. Once we've done that we'll assemble the
                    # complete animal records from the pieces in the relevant
                    # lookup tables.
                    if c == len(_curr_rec)-1:
                        if not any_sexes:
                            _sex = missing_sex.upper()
                        indi2sex[_indi] = _sex
                        if _byr:
                            indi2birth[_indi] = _byr
                        else:
                            indi2birth[_indi] = missing_byear
                        if _famc:
                            indi2famc[_indi] = _famc
                        if _fams:
                            indi2fams[_indi] = _fams
                        if '_name':
                            if not any_names: any_names = 1
                            indi2name[_indi] = _name
                        else:
                            id2name[_id] = missing_name
            elif firstlinelist[-1] == 'FAM':
                _husb, _wife,_child = 0, 0, []
                end_family = 0
                for c in xrange(1,len(_curr_rec)):
                    linelist = _curr_rec[c].strip().rstrip('\r').rstrip('\n').split(' ')
                    if c == 1:
                        _fam = _curr_rec[0].strip().split(' ')[1][1:-1]
                    if len(linelist) > 1:
                        _tag = linelist[1].upper()
                    else: break
                    if _tag in known_tags:
                        if messages == 'debug':
                            print 'Processing tag %s' % ( _tag )
                        if _tag.upper() == 'HUSB':
                            _husb =  linelist[2][1:-1]
                        elif _tag.upper() == 'WIFE':
                            _wife =  linelist[2][1:-1]
                        elif _tag.upper() == 'CHIL':
                            _child.append(linelist[2][1:-1])
                        else:
                            if messages == 'debug':
                                print 'Skipping unknown tag %s' % ( _tag )
                    if c == len(_curr_rec)-1:
                        if _husb:
                            fam2husb[_fam] = _husb
                        if _wife:
                            fam2wife[_fam] = _wife
                        if _child:
                            for _ch in _child:
                                fam2chil[_fam] = {}
                                fam2chil[_fam][_ch] = _ch

        if debug:
            print 'indi2sex: ', indi2sex
            print 'indi2birth: ', indi2birth
            print 'indi2name: ', indi2name
            print 'indi2famc: ', indi2famc
            print 'indi2fams: ', indi2fams
            print 'fam2husb: ', fam2husb
            print 'fam2wife: ', fam2wife

        # Now we walk through the INDI records and assemble everything. Note
        # that there is no error-checking here. I don't want to have error-
        # checking code in two different places, namely here and in the
        # NewPedigree::preprocess() method. Therefore, all error-checking is
        # deferred to the class method.
        assembled = {}
        for i in indi2sex.keys():
            assembled[i] = {}
            assembled[i]['indi'] = i
            assembled[i]['sex'] = indi2sex[i]
            assembled[i]['birth'] = indi2birth[i]
            assembled[i]['name'] = indi2name[i]
            try:
                assembled[i]['sire'] = fam2husb[indi2famc[i]]
            except KeyError:
                assembled[i]['sire'] = missing_parent
            try:
                assembled[i]['dam'] = fam2wife[indi2famc[i]]
            except KeyError:
                assembled[i]['dam'] = missing_parent
        if debug:
            print 'assembled: ', assembled
        if messages == 'verbose':
            print '[INFO]: Successfully imported pedigree from the GEDCOM file %s!' \
                % ( infilename )
        logging.info('Successfully imported pedigree from the GEDCOM file %s!',infilename)
    except:
        if messages == 'verbose':
            print '[ERROR]: Unable to import pedigree from the GEDCOM file %s!' \
                % ( infilename )
        logging.error('Unable to import pedigree from the GEDCOM file %s!',infilename)

    # Save the GEDCOM file in ASD format and update the
    # the value of pedfile.
    try:
        outfilename = '%s.tmp' % ( infilename )
        pedformat = save_from_gedcom(outfilename,assembled)
    except:
        pedformat = 'xxxx'
    return pedformat

##
# save_from_gedcom() takes pedigree data parsed by load_from_gedcom() and
# writes it to a text file in an ASD format that PyPedal can easily read.
# @param outfilename The file to which the records should be written.
# @param assembled A list of records read from a GEDCOM input file
# @retval A string containing the pedigree format code. 'xxxx' if there was a problem.
def save_from_gedcom(outfilename, assembled):
    """
    save_from_gedcom() takes pedigree data parsed by load_from_gedcom() and
    writes it to a text file in an ASD format that PyPedal can easily read.
    """
    pedformat = 'xxxx'
    try:
        ofh = file(outfilename,'w')
        for _i in assembled.keys():
            pedformat = 'ASDxbu'
            outstring = '%s,%s,%s,%s,%s,%s\n' % ( \
                assembled[_i]['indi'], \
                assembled[_i]['sire'], \
                assembled[_i]['dam'], \
                assembled[_i]['sex'], \
                assembled[_i]['birth'], \
                assembled[_i]['name'], \
            )
            ofh.write(outstring)
        ofh.close()
        logging.info('Saved GEDCOM pedigree to the file %s!',outfilename)
    except:
        logging.error('Unable to save GEDCOM pedigree to the file %s!',outfilename)
    return pedformat

##
# save_to_gedcom() writes a PyPedal NewPedigree object to a file in
# GEDCOM 5.5 format.
# @param pedobj An instance of a PyPedal NewPedigree object
# @param outfilename The file to which the matrix should be written
# @retval A save status indicator (0: failed, 1: success).
def save_to_gedcom(pedobj, outfilename):
    """
    save_to_gedcom() writes a PyPedal NewPedigree object to a file in
    GEDCOM 5.5 format.
    """
    try:
        ofh = file(outfilename,'w')
        # Write file header
        ofh.write('0 HEAD\n')
        ofh.write('1 SOUR PYPEDAL\n')
        ofh.write('2 VERS V2.0\n')
        ofh.write('2 CORP USDA-ARS-BA-ANRI-AIPL\n')
        ofh.write('1 DEST PYPEDAL\n')
        ofh.write('1 DATE %s\n' % (time.strftime('%m %d %Y', \
            (time.localtime(time.time())))) )
        ofh.write('1 FILE %s\n' % (pedobj.kw['pedfile']) )
        ofh.write('1 GEDC\n')
        ofh.write('2 VERS 5.5\n')
        ofh.write('2 FORM Lineage-Linked\n')
        ofh.write('1 CHAR ASCII\n')
        # Fill the file
        indi = {}
        fam = {}
        par2spouses = {}
        # Sweep the pedigree once to find allnon-founders and map the founders
        # to families so that we can correctly assign FAMS tags to founder records.
        for p in pedobj.pedigree:
            if p.sireID != pedobj.kw['missing_parent'] or p.damID != pedobj.kw['missing_parent']:
                if not par2spouses.has_key(p.sireID):
                    par2spouses[p.sireID] = []
                if not par2spouses.has_key(p.damID):
                    par2spouses[p.damID] = []
                if p.sireName == pedobj.kw['missing_name'] and p.damName == pedobj.kw['missing_name']:
                    _spouses = 'F0_0'
                    if p.sireID != pedobj.kw['missing_parent']:
                        if _spouses not in par2spouses[p.sireID]:
                            par2spouses[p.sireID].append(_spouses)
                    if p.damID != pedobj.kw['missing_parent']:
                        if _spouses not in par2spouses[p.damID]:
                            par2spouses[p.damID].append(_spouses)
                elif p.sireName == pedobj.kw['missing_name']:
                    _spouses = 'F%s' % ( p.damName )
                    if p.damID != pedobj.kw['missing_parent']:
                        if _spouses not in par2spouses[p.damID]:
                            par2spouses[p.damID].append(_spouses)
                elif p.damName == pedobj.kw['missing_name']:
                    _spouses = 'F%s' % ( p.sireName )
                    if p.sireID != pedobj.kw['missing_parent']:
                        if _spouses not in par2spouses[p.sireID]:
                            par2spouses[p.sireID].append(_spouses)
                else:
                    _spouses = 'F%s_%s' % ( p.sireName, p.damName )
                    if p.sireID != pedobj.kw['missing_parent']:
                        if _spouses not in par2spouses[p.sireID]:
                            par2spouses[p.sireID].append(_spouses)
                    if p.damID != pedobj.kw['missing_parent']:
                        if _spouses not in par2spouses[p.damID]:
                            par2spouses[p.damID].append(_spouses)
        # Create INDI and FAM records for each animal in the pedigree.
        for p in pedobj.pedigree:
            if p.sireName == pedobj.kw['missing_name'] and p.damName == pedobj.kw['missing_name']:
                _fam = 'F0_0'
            elif p.sireName == pedobj.kw['missing_name']:
                _fam = 'F%s' % ( p.damName )
            elif p.damName == pedobj.kw['missing_name']:
                _fam = 'F%s' % ( p.sireName )
            else:
                _fam = 'F%s_%s' % ( p.sireName, p.damName )
            if not indi.has_key(p.animalID):
                indi[p.animalID] = '0 @%s@ INDI\n' % \
                    ( pedobj.namebackmap[pedobj.backmap[p.animalID]] )
                indi[p.animalID] = '%s1 SEX %s\n' % ( indi[p.animalID], p.sex.upper() )
                if 'n' in pedobj.kw['pedformat']:
                    indi[p.animalID] = '%s1 NAME %s\n' % ( indi[p.animalID], p.name )
                elif 'u' in pedobj.kw['pedformat']:
                    indi[p.animalID] = '%s1 NAME %s\n' % ( indi[p.animalID], p.userField )
                else:
                    pass
                if 'y' in pedobj.kw['pedformat'] and p.by != pedobj.kw['missing_byear']:
                    indi[p.animalID] = '%s1 BIRT\n' % ( indi[p.animalID] )
                    indi[p.animalID] = '%s2 DATE %s\n' % ( indi[p.animalID], p.by )
                if 'b' in pedobj.kw['pedformat'] and ( p.bd != pedobj.kw['missing_bdate'] and p.bd != str(pedobj.kw['missing_byear']) ):
                    indi[p.animalID] = '%s1 BIRT\n' % ( indi[p.animalID] )
                    indi[p.animalID] = '%s2 DATE %s\n' % ( indi[p.animalID], p.bd )
                if _fam != 'F0_0':
                    indi[p.animalID] = '%s1 FAMC @%s@\n' % ( indi[p.animalID], _fam )
                if par2spouses.has_key(p.animalID):
                    for _p2s in par2spouses[p.animalID]:
                        if _p2s != 'F0_0':
                            indi[p.animalID] = '%s1 FAMS @F%s@\n' % ( indi[p.animalID], _p2s )
            # Create the family if it does not yet exist
            if not fam.has_key(_fam) and _fam[3:7] != 'F0_0':
                fam[_fam] = '0 @%s@ FAM\n' % ( _fam )
            if 'HUSB' not in fam[_fam] and p.sireName != pedobj.kw['missing_name']:
                fam[_fam] = '%s1 HUSB @%s@\n' % ( fam[_fam],p.sireName )
            if 'WIFE' not in fam[_fam] and p.damName != pedobj.kw['missing_name']:
                fam[_fam] = '%s1 WIFE @%s@\n' % ( fam[_fam], p.damName )
            fam[_fam] = '%s1 CHIL @%s@\n' % ( fam[_fam], \
                pedobj.namebackmap[pedobj.backmap[p.animalID]] )
            #print fam[_fam]
        # Now loop and write the contents of indi and fam to the file.
        for i in indi.values():
            ofh.write(i)
        for f in fam.values():
            ofh.write(f)
        # Write footer and close file
        ofh.write('0 TRLR\n')
        ofh.close()
        if pedobj.kw['messages'] == 'verbose':
            print '[INFO]: Successfully exported pedigree to the GEDCOM file %s!' \
                % ( outfilename )
        logging.info('Successfully exported pedigree to the GEDCOM file %s!',outfilename)
        return 1
    except:
        if pedobj.kw['messages'] == 'verbose':
            print '[ERROR]: Unable to save pedigree to the GEDCOM file %s!' \
                % ( outfilename )
        logging.error('Unable to export pedigree to the GEDCOM file %s!',outfilename)
        return 0

##
# load_from_genes() reads and parses pedigree data that conforms to
# the DBF format used by GENES software for pedigree management v1.2
# (R. Lacey, http://www.vortex9.org/genes.html). When possible data
# are mapped into similar PyPedal fields.
# @param infilename The file from which the pedigree should be read.
# @param messages Controls output to the screen
# @param standalone Uses logging if called by a NewPedigree method
# @param missing_sex Value assigned to an animal with unknown sex
# @param missing_parent Value assigned to unknown parents
# @param missing_name Name assigned by default
# @param missing_bdate Value assigned to unknown birthdates
# @param debug Flag turning debugging messages on (1) and off (0)
# @retval A save status indicator (0: failed, 1: success).
def load_from_genes(infilename, messages='verbose', standalone=1, missing_sex='u', \
    missing_parent=0, missing_name='Unknown Name', missing_bdate='01011900', debug=0):
    """
    load_from_genes() reads and parses pedigree data from the dBase III
    files used by GENES 1.20 (http://www.vortex9.org/genes.html).
    """

    import struct, datetime, decimal, itertools

    record_list = []
    f = open(infilename, 'rb')

    ### Begin code of Raymond Hettinger's taken from http://code.activestate.com/recipes/362715/
    ### and modified slightly for PyPedal.
    # fields is list of field names and field specs: (type, size, decimal places).
    # record_list contains data records.
    # If a record is marked as deleted, it is skipped.
    #
    # File should be opened for binary reads.
    numrec, lenheader = struct.unpack('<xxxxLH22x', f.read(32))
    numfields = (lenheader - 33) // 32
    fields = []
    for fieldno in xrange(numfields):
        name, typ, size, deci = struct.unpack('<11sc4xBB14x', f.read(32))
        name = name.replace('\0', '')       # eliminate NULs from string
        fields.append((name, typ, size, deci))
    #yield [field[0] for field in fields]
    #print fields
    #if debug:
    #    print [field[0] for field in fields]
    #yield [tuple(field[1:]) for field in fields]
    #if debug:
    #    print [tuple(field[1:]) for field in fields]
    terminator = f.read(1)
    assert terminator == '\r'
    fields.insert(0, ('DeletionFlag', 'C', 1, 0))
    fmt = ''.join(['%ds' % fieldinfo[2] for fieldinfo in fields])
    fmtsiz = struct.calcsize(fmt)
    for i in xrange(numrec):
        record = struct.unpack(fmt, f.read(fmtsiz))
        if record[0] != ' ':
            continue                        # deleted record
        result = []
        for (name, typ, size, deci), value in itertools.izip(fields, record):
            if name == 'DeletionFlag':
                continue
            if typ == "N":
                value = value.replace('\0', '').lstrip()
                if value == '':
                    value = 0
                elif deci:
                    value = value.replace('\x1a', '').lstrip() # String ctrl-Z line ending
                    value = decimal.Decimal(value)
                    value = float(value)
                else:
                    value = int(value)
            elif typ == 'D':
                y, m, d = value[:4], value[4:6], value[6:8]
                value=d+"/"+m+"/"+y
            elif typ == 'L':
                value = (value in 'YyTt' and 'T') or (value in 'NnFf' and 'F') or '?'
            elif typ == 'F':
                value = float(value)
            result.append(value)
        record_list.append(result)
    ### End code of Raymond Hettinger's taken from http://code.activestate.com/recipes/362715/
    f.close()
    # Check inputs visually
    #if debug:
    #    for r in record_list:
    #        print r
    #    print [field[0] for field in fields]
    # Here we go...
    field_list = [field[0] for field in fields]
    assembled = {}
    # As is also the case for GEDCOM record processing, no consistency checks are performed in this
    # subroutine. The only "edits" made are to insure that values, e.g. sex codes, are consistent with
    # the values expected by PyPedal.
    for r in record_list:
        anidx = field_list.index('STUD_ID') - 1         # The subtraction deals with the DeletionFlag
        assembled[r[anidx]] = {}
        assembled[r[anidx]]['indi'] = r[anidx].strip()
        # Unlike GENES, PyPedal does not distinguish among animals that are known founders and those that
        # have unknown parents but are not founders.
        assembled[r[anidx]]['sire'] = r[field_list.index('SIRE_ID') - 1].strip()     # "   UNK" = unknown ("  WILD" = founder)
        assembled[r[anidx]]['dam'] = r[field_list.index('DAM_ID') - 1].strip()      # "   UNK" = unknown ("  WILD" = founder)
        if assembled[r[anidx]]['sire'].strip() == 'UNK': assembled[r[anidx]]['sire'] = missing_parent
        if assembled[r[anidx]]['sire'].strip() == 'WILD': assembled[r[anidx]]['sire'] = missing_parent
        if assembled[r[anidx]]['dam'].strip() == 'UNK': assembled[r[anidx]]['dam'] = missing_parent
        if assembled[r[anidx]]['dam'].strip() == 'WILD': assembled[r[anidx]]['dam'] = missing_parent
        assembled[r[anidx]]['name'] = r[anidx].strip()
        # Sex codes are different
        assembled[r[anidx]]['sex'] = r[field_list.index('SEX') - 1]          # 0 = female, 1 = male
        if str(assembled[r[anidx]]['sex']) == '0': assembled[r[anidx]]['sex'] = 'f'
        elif str(assembled[r[anidx]]['sex']) == '1': assembled[r[anidx]]['sex'] = 'm'
        else: assembled[r[anidx]]['sex'] = missing_sex
        # Flip alive/dead status
        assembled[r[anidx]]['alive'] = r[field_list.index('DEAD') - 1]       # T/F indicating whether the animal is dead (T = dead)
        if assembled[r[anidx]]['alive'] == 'T': assembled[r[anidx]]['alive'] = '1'
        else: assembled[r[anidx]]['alive'] = '0'
        # Inbreeding
        assembled[r[anidx]]['fa'] = r[field_list.index('INBREED') - 1]
        # Birthdates can be missing
        assembled[r[anidx]]['bd'] = r[field_list.index('BDATE') - 1]
        pieces = assembled[r[anidx]]['bd'].split("/")
        assembled[r[anidx]]['bd'] = '%s%s%s' % ( pieces[0], pieces[1], pieces[2] )
        if assembled[r[anidx]]['bd'] == '': assembled[r[anidx]]['bd'] = missing_bdate
        # Use herd as a proxy for location
        assembled[r[anidx]]['herd'] = r[field_list.index('LOCATION') - 1].strip()
        # Age should be okay as-is
        assembled[r[anidx]]['age'] = r[field_list.index('AGE') - 1]
        #if debug:
        #    print assembled[r[anidx]]
    if messages == 'verbose':
        print '[INFO]: Successfully imported pedigree from the GENES 1.20 file %s!' \
            % ( infilename )
    logging.info('Successfully imported pedigree from the GENES 1.20 file %s!',infilename)

    # Save the GENES file in ASD format and update the
    # the value of pedfile.
    try:
        outfilename = '%s.tmp' % ( infilename )
        pedformat = save_from_genes(outfilename,assembled)
    except:
        pedformat = 'xxxx'
    return pedformat

##
# save_from_genes() takes pedigree data parsed by load_from_genes() and
# writes it to a text file in an ASD format that PyPedal can easily read.
# @param outfilename The file to which the records should be written.
# @param assembled A list of records read from a GENES 1.20 input file
# @retval A string containing the pedigree format code. 'xxxx' if there was a problem.
def save_from_genes(outfilename, assembled):
    """
    save_from_genes() takes pedigree data parsed by load_from_genes() and
    writes it to a text file in an ASD format that PyPedal can easily read.
    """
    pedformat = 'xxxx'
    try:
        ofh = file(outfilename,'w')
        for _i in assembled.keys():
            pedformat = 'ASDefHlnxy'
            outstring = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % ( \
                assembled[_i]['indi'], \
                assembled[_i]['sire'], \
                assembled[_i]['dam'], \
                assembled[_i]['age'], \
                assembled[_i]['fa'], \
                assembled[_i]['herd'], \
                assembled[_i]['alive'], \
                assembled[_i]['name'], \
                assembled[_i]['sex'], \
                assembled[_i]['bd'], \
            )
            ofh.write(outstring)
        ofh.close()
        logging.info('Saved GENES 1.20 pedigree to the file %s!',outfilename)
    except:
        logging.error('Unable to save GENES 1.20 pedigree to the file %s!',outfilename)
    return pedformat

##
# save_to_genes() writes a PyPedal NewPedigree object to a file in
# GENES 1.20 (dBase III) format.
# @param pedobj An instance of a PyPedal NewPedigree object
# @param outfilename The file to which the matrix should be written
# @retval A save status indicator (0: failed, 1: success).
def save_to_genes(pedobj, outfilename):
    """
    save_to_genes() writes a PyPedal NewPedigree object to a file in
    GENES 1.20 (dBase III) format.
    """
    import struct, datetime, decimal, itertools

    # These are the field names and data-types per the specification in Lacy's GENES.DOC.
    fields = [ ('STUD_ID', 'C', 6, 0), ('DAM_ID', 'C', 6, 0), ('SIRE_ID', 'C', 6, 0), ('BDATE', 'D', 8, 0),
               ('SEX', 'N', 1, 0), ('DATEOUT', 'D', 8, 0), ('DEATHDATE', 'D', 8, 0), ('LOCATION', 'C', 9, 0),
               ('SELECTED', 'L', 1, 0), ('DEAD', 'L', 1, 0), ('INBREED', 'N', 8, 4), ('AGE', 'N', 3, 0),
               ('KNOWN', 'N', 8, 4), ('INBREED_KN', 'N', 8, 4), ('MK', 'N', 8, 4), ('MK_KN', 'N', 8, 4),
               ('KV', 'N', 8, 4), ('KV_KN', 'N', 8, 4), ('VX', 'N', 8, 4), ('GU_ALL', 'N', 8, 4),
               ('GU_DESC', 'N', 8, 4), ('PR_LOST', 'N', 8, 4) ]

    # Dictionary to map GENES fields to PyPedal NewAnimal attributes
    pyp2genes = { 'STUD_ID': 'originalID', 'SIRE_ID': 'sireName', 'DAM_ID': 'damName', 'STUD_ID': 'animalID', 'BDATE': 'bd',
          'SEX': 'sex', 'LOCATION': 'originalHerd', 'DEAD': 'alive', 'INBREED': 'fa', 'AGE': 'age'
    }

    # Go ahead and try exporting the pedigree.
    try:
        ### Begin code of Raymond Hettinger's taken from http://code.activestate.com/recipes/362715/
        ### and modified slightly for PyPedal.
        f = open(outfilename, 'wb')
        # Write header info
        ver = 3
        now = datetime.datetime.now()
        yr, mon, day = now.year-1900, now.month, now.day
        numrec = len(pedobj.pedigree)
        numfields = len(fields)
        lenheader = numfields * 32 + 33
        lenrecord = sum(field[2] for field in fields) + 1
        hdr = struct.pack('<BBBBLHH20x', ver, yr, mon, day, numrec, lenheader, lenrecord)
        f.write(hdr)
        # Write field specs
        for fld in fields:
            (name, typ, size, deci) = fld
            name = name.ljust(11, '\x00')
            fld = struct.pack('<11sc4xBB14x', name, typ, size, deci)
            f.write(fld)
        # terminator
        f.write('\r')
        # Write individual records for record in records:
        for p in pedobj.pedigree:
            f.write(' ')                        # deletion flag
            for fld in fields:
                (name, typ, size, deci) = fld
                # Assign values or blanks, depending on whether or not PyPedal and GENES information
                # are correspondent.
                if pyp2genes.has_key(name):
                    value = getattr(p, pyp2genes[name])
                else:
                    value = ' '
                # We have to do a little re-mapping to deal with, e.g., sex codes.
                if value == pedobj.kw['missing_name']: value = '   UNK'
                if name == 'SEX':
                    if value == 'f': value = 0
                    elif value == 'm': value = 1
                    else: value = 0                 # Assume that unknown sex animals are females
                if name == 'DEAD':
                    if value == 0: value = 'T'
                    else: value = 'F'
                # Now cast everything and pad out fields to the correct length.
                if typ == 'N':
                    value = str(value).rjust(size, ' ')
                elif typ == 'D':
                    if value == ' ':
                        value = ' '*size
                    else:
                        value = '%s%s' % ( value[4:8], value[0:4] )
                elif typ == 'L':
                    value = str(value)[0].upper()
                else:
                    # If we have values that exceed the width of the field truncate them and warn the user.
                    if len(str(value)) > size:
                        value = str(value[0:size+1])
                        if messages == 'verbose':
                            print '[WARNING]: Truncated field %s while exporting to GENES 1.20 file %s!' % ( name, outfilename )
                        logging.warn('Truncated field %s while exporting to GENES 1.20 file %s!', name, outfilename)
                    value = str(value)[:size].ljust(size, ' ')
                # Double-check that lengths are correct before writing the record
                assert len(value) == size
                f.write(value)
        # End-of-file marker
        f.write('\x1A')
        ### End code of Raymond Hettinger's taken from http://code.activestate.com/recipes/362715/
        ### and modified slightly for PyPedal.
        f.close()
        if pedobj.kw['messages'] == 'verbose':
            print '[INFO]: Successfully exported pedigree to the GENES file %s!' \
                % ( outfilename )
        logging.info('Successfully exported pedigree to the GENES file %s!',outfilename)
        return 1
    # If the pedigree could not be exported then tell the user that something went wrong.
    except:
        if pedobj.kw['messages'] == 'verbose':
            print '[ERROR]: Unable to export pedigree to the GENES file %s!' \
                % ( outfilename )
        logging.error('Unable to export pedigree to the GENES file %s!',outfilename)
        return 0

##
# save_newanimals_to_file() take a list of PyPedal NewAnimal objects as input and writes them to a pedigree file.
# @param animal_list A list of PyPedal NewAnimal bjects
# @param filename The nae of the file t which the animals should be written
# @param pedformat Pedigree format code for the output file
def save_newanimals_to_file(animal_list, filename, pedformat, sepchar):
	if len(animal_list) == 0:
		pass
	else:
		# First, save the unique animals from the union of pedigrees a and
		# b based on the match rule. Note that the pedformat from the first
		# pedigree passed to __add__() will be used for both pedigrees. This
		# makes sense because you cannot have two different pedformats in
		# the same file.
		try:
			f = open(filename, 'w')
			for animal in animal_list:
				_outstring = ''
				for pf in pedformat:
					if originalID == False:
						value = getattr(_a, self.new_animal_attr[pf])
					else:
						if pf in['a','A']:
							value = _a.originalID
						# This cascade may break if the pedigree is not
						# renumbered...
						elif pf in['s','S']:
							if _a.sireID != self.kw['missing_parent']:
								value = self.pedigree[_a.sireID-1].originalID
							else:
								value = 0
						elif pf in['d','D']:
							if _a.damID != self.kw['missing_parent']:
								value = self.pedigree[_a.damID-1].originalID
							else:
								value = 0
						else:
							value = getattr(_a, self.new_animal_attr[pf])
				# If we don't catch the special case of the first entry
				# in an output line the a sepchar always will be the
				# first character in the line.
				if len(_outstring) > 0:
					_outstring = '%s%s%s' % ( _outstring, sepchar, value )
				else:
					_outstring = '%s' % ( value )
				ofh.write( '%s\n' % (_outstring) )
			f.close()
		except:
			pass
