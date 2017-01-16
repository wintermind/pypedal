#!/usr/bin/python

###############################################################################
# NAME: pyp_nrm.py
# VERSION: 2.0.0 (29SEPTEMBER2010)
# AUTHOR: John B. Cole, PhD (john.cole@ars.usda.gov)
# LICENSE: LGPL
################################################################################ FUNCTIONS:
#   a_matrix()
#   fast_a_matrix()
#   fast_a_matrix_r()
#   inbreeding()
#   inbreeding_vanraden()
#   recurse_pedigree()
#   recurse_pedigree_n()
#   recurse_pedigree_onesided()
#   recurse_pedigree_idonly()
#   inbreeding_tabular()
#   inbreeding_meuwissen_luo()
#   inbreeding_modified_meuwissen_luo()
#   a_decompose()
#   form_d_nof()
#   a_inverse_dnf()
#   a_inverse_df()
#   partial_inbreeding()
#   fast_partial_a_matrix()
###############################################################################

## @package pyp_nrm
# pyp_nrm contains several procedures for computing numerator relationship matrices and for
# performing operations on those matrices.  It also contains routines for computing CoI on
# large pedigrees using the recursive method of VanRaden (1992).
##
from __future__ import print_function
try:
    from pysparse import spmatrix
except ImportError:
    #logging.info('Could not import the spmatrix module from PySparse! Using NumPY dense matrices instead.')
    print('[INFO]: Could not import the spmatrix module from PySparse in pyp_nrm! NumPY dense matrices will be used instead.')
import pyp_utils

##
# a_matrix() is used to form a numerator relationship matrix from a pedigree.  DEPRECATED.
# use fast_a_matrix() instead.
# @param pedobj A PyPedal pedigree object.
# @param save Flag to indicate whether or not the relationship matrix is written to a file.
# @retval The NRM as a NumPy matrix.
def a_matrix(pedobj, save=0):
    """
    Form a numerator relationship matrix from a pedigree.  DEPRECATED.
    """
    try: logging.info('Entered a_matrix()')
    except: pass
    l = pedobj.medata.num_records
    # Grab some array tools
    try:
        a = numpy.zeros([l,l],'d')  # initialize a matrix of zeros of appropriate dimension
        for row in xrange(l):
            for col in xrange(row,l):
                # cast these b/c items are read from the pedigree file as characters, not integers
                pedobj.pedigree[col].animalID = int(pedobj.pedigree[col].animalID)
                pedobj.pedigree[col].sireID = int(pedobj.pedigree[col].sireID)
                pedobj.pedigree[col].damID = int(pedobj.pedigree[col].damID)
                #if pedobj.pedigree[col].sireID == 0 and pedobj.pedigree[col].damID == 0:
                if str(pedobj.pedigree[col].sireID) == str(pedobj.kw['missing_parent']) and str(pedobj.pedigree[col].damID) == str(pedobj.kw['missing_parent']):
                    if row == col:
                        # both parents unknown and assumed unrelated
                        a[row,col] = 1.
                    else:
                        a[row,col] = 0.
                        a[col,row] = a[row,col]
                elif str(pedobj.pedigree[col].sireID) == str(pedobj.kw['missing_parent']):
                    # sire unknown, dam known
                    if row == col:
                        a[row,col] = 1.
                    else:
                        a[row,col] = 0.5 * a[row,pedobj.pedigree[col].damID-1]
                        a[col,row] = a[row,col]
                elif str(pedobj.pedigree[col].damID) == str(pedobj.kw['missing_parent']):
                    # sire known, dam unknown
                    if row == col:
                        a[row,col] = 1.
                    else:
                        a[row,col] = 0.5 * a[row,pedobj.pedigree[col].sireID-1]
                        a[col,row] = a[row,col]
                elif str(pedobj.pedigree[col].sireID) != str(pedobj.kw['missing_parent']) and str(pedobj.pedigree[col].damID) != str(pedobj.kw['missing_parent']):
                    # both parents known
                    if row == col:
                        a[row,col] = 1. + ( 0.5 * a[pedobj.pedigree[col].sireID-1,pedobj.pedigree[col].damID-1] )
                    else:
                        intermediate = a[row,pedobj.pedigree[col].sireID-1] + a[row,pedobj.pedigree[col].damID-1]
                        finprod = 0.5 * intermediate
                        a[row,col] = 0.5 * intermediate
                        a[col,row] = a[row,col]
                else:
                    print('[ERROR]: There is a problem with the sire (ID %s) and/or dam (ID %s) of animal %s' % (pedobj.pedigree[col].sireID,pedobj.pedigree[col].damID,pedobj.pedigree[col].animalID))
                    break
    except:
        a = numpy.zeros([1,1],'d')

    if save:
        a_outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_a_matrix_','.dat')
        aout = open(a_outputfile,'w')
        for row in xrange(l):
            line = ''
            for col in xrange(l):
                if col == 0:
                    line = '%7.5f' % (a[row,col])
                else:
                    line = '%s%s%s' % (line,',',a[row,col])
            line = '%s%s' % (line,'\n')
            aout.write(line)
        aout.close()

    try: logging.info('Exited a_matrix()')
    except: pass
    return a

##
# Form a numerator relationship matrix from a pedigree.  fast_a_matrix() is a hacked version of a_matrix()
# modified to try and improve performance.  Lists of animal, sire, and dam IDs are formed and accessed rather
# than myped as it is much faster to access a member of a simple list rather than an attribute of an object in a
# list.  Further note that only the diagonal and upper off-diagonal of A are populated.  This is done to save
# n(n+1) / 2 matrix writes.  For a 1000-element array, this saves 500,500 writes.
# @param pedigree A PyPedal pedigree.
# @param pedopts PyPedal options.
# @param save Flag to indicate whether or not the relationship matrix is written to a file.
# @param method Use dense or sparse matrix storage.
# @param debug Print NRM for debugging
# @param fill Fill both the upper and lower off-diagonals when 1, only the upper otherwise.
# @retval The NRM as Numarray matrix.
#@profile
def fast_a_matrix(pedigree, pedopts, save=0, method='dense', debug=0, fill=1):
    """
    Form a numerator relationship matrix from a pedigree.  fast_a_matrix() is a hacked
    version of a_matrix() modified to try and improve performance.  Lists of animal,
    sire, and dam IDs are formed and accessed rather than myped as it is much faster to
    access a member of a simple list rather than an attribute of an object in a
    list.  Further note that only the diagonal and upper off-diagonal of A are
    populated.  This is done to save n(n+1) / 2 matrix writes.  For a 1000-element
    array, this saves 500,500 writes.
    """
#     try: logging.info('Entered fast_a_matrix()')
#     except: pass
#     print('\t\t\tEntered pyp_nrm.fast_a_matrix()')
#    print('Entered pyp_nrm.fast_a_matrix()')
#    print(pedopts)
#    print('[DEBUG]: method = ', method)

    # Check the foundercoi flag and do what it says.
    foundercoi = int(pedopts['foundercoi'])
    # Check to make sure that foundercoi is either 0 or 1.
    if foundercoi not in [0, 1]: foundercoi = 0

    _animals = {}
    _sires = {}
    _dams = {}
    _str_sires = {}
    _str_dams = {}
    l = len(pedigree)
    #method = 'dense'
    if method not in ['dense','sparse']:
        method = 'dense'
    #try:
    # Use PySparse to provide sparse matrix storage for large
    # relationship matrices.
    if method == 'sparse':
        try:
#            from pysparse import spmatrix
            a = spmatrix.ll_mat_sym(l*l)
            for i in xrange(l):
                a[i,i] = 1.
        except:
#            logging.error('Could not import the spmatrix module from PySparse! Using NumPY dense matrices instead.')
#            print('[ERROR]: Could not import the spmatrix module from PySparse! Using NumPY dense matrices instead.')
            #a = numpy.zeros([l,l],'d')
	    try:
                a = numpy.zeros([l,l],'d')  # initialize a matrix of zeros of appropriate size
	    except MemoryError:
	        a = numpy.memmap('fast_a_matrix_mmap.bin', dtype='float32', mode='w+', shape=(l,l))
	    except:
	        print('[ERROR]: Unable to allocate a matrix of rank %s in pyp_nrm.fast_a_matrix()!' % ( l ))
	        logging.error('Unable to allocate a matrix of rank %s in pyp_nrm.fast_a_matrix()!', l)
	        return False
    # Otherwise, use Numpy and its dense matrices
    else:
	# First, try and allocate the vectors in RAM. If that does not work, try and allocate them using
        # memory-mapped files. If that does not work, well, give up.
	try:
            a = numpy.zeros([l,l],'d')  # initialize a matrix of zeros of appropriate size
	except MemoryError:
	    a = numpy.memmap('fast_a_matrix_mmap.bin', dtype='float32', mode='w+', shape=(l,l))
	except:
	    print('[ERROR]: Unable to allocate a matrix of rank %s in pyp_nrm.fast_a_matrix()!' % ( l ))
	    logging.error('Unable to allocate a matrix of rank %s in pyp_nrm.fast_a_matrix()!', l)
	    return False
        # print(a)
    if pedopts['debug_messages'] and pedopts['messages'] != 'quiet':
        print('\t\t[pyp_nrm/fast_a_matrix()] Started forming animal, sire, and dam lists at %s' %  pyp_utils.pyp_nice_time())
    #print('\tIdx\tID\tName')
    for i in xrange(l):
	a[i,i] = 1.0
	if foundercoi == 1:
            if str(pedigree[i].sireID) == str(pedopts['missing_parent']) and str(pedigree[i].damID) == str(pedopts['missing_parent']):
                a[i,i] = 1.0 + pedigree[i].fa
        try:
            _a = _animals[i]
        except KeyError:
            _animals[i] = int(pedigree[i].animalID)
            #print('A:\t', i, '\t', pedigree[i].animalID, '\t', pedigree[i].name)
        try:
            _s = _sires[i]
        except KeyError:
            _sires[i] = int(pedigree[i].sireID)
	    _str_sires[i] = str(_sires[i])
            #print('S:\t', i, '\t', pedigree[i].sireID, '\t', pedigree[i].sireName)
        try:
            _d = _dams[i]
        except KeyError:
            _dams[i] = int(pedigree[i].damID)
            _str_dams[i] = str(_dams[i])
            #print('D:\t', i, '\t', pedigree[i].damID, '\t', pedigree[i].damName)
    if pedopts['debug_messages'] and pedopts['messages'] != 'quiet':
        print('\t\t[pyp_nrm/fast_a_matrix()] Finished forming animal, sire, and dam lists at %s' %  pyp_utils.pyp_nice_time())
        print('\t\t[pyp_nrm/fast_a_matrix()] Started computing A at %s' %  pyp_utils.pyp_nice_time())
    #print('_animals: ', _animals)
    #print('_sires: ', _sires)
    #print('_dams: ', _dams)
    #print('missing_parent: ', pedopts['missing_parent'])
    #print(l)
    this_msg = str(pedopts['missing_parent'])
    for row in xrange(l):
        for col in xrange(row,l):
            if _str_sires[col] != this_msg and _str_dams[col] != this_msg:
                # both parents known
                if row == col:
                    a[row,col] = a[row,col] + ( 0.5 * a[_sires[col]-1,_dams[col]-1] )
                else:
                    a[row,col] = 0.5 * ( a[row,_sires[col]-1] + a[row,_dams[col]-1] )
                    a[col,row] = a[row,col]
            else:
                if _str_sires[col] == this_msg and _str_dams[col] == this_msg:
                    # sire and dam unknown
                    pass
                elif _str_sires[col] == this_msg and _str_dams[col] != this_msg:
                    # sire unknown, dam known
                    if row != col:
                        a[row,col] = 0.5 * a[row,_dams[col]-1]
                        a[col,row] = a[row,col]
		else:
                    # sire known, dam unknown
                    if row != col:
                        a[row,col] = 0.5 * a[row,_sires[col]-1]
                        a[col,row] = a[row,col]
    if pedopts['debug_messages'] and pedopts['messages'] != 'quiet':
        print('\t\t[pyp_nrm/fast_a_matrix()] Finished computing A at %s' %  pyp_utils.pyp_nice_time())
    #except:
    #    a = numpy.zeros([l,l],'d')  # initialize a matrix of zeros of appropriate order

    if save == 1:
        a_outputfile = '%s%s%s' % (pedopts['filetag'],'_new_a_matrix_','.dat')
        aout = open(a_outputfile,'w')
        label = 'Produced by pyp_nrm/fast_a_matrix()\n'
        aout.write(label)
        for row in xrange(l):
            line = ''
            for col in xrange(l):
                if col == 0:
                    line = '%7.5f' % (a[row,col])
                else:
                    line = '%s%s%s' % (line,',',a[row,col])
            line = '%s%s' % (line,'\n')
            aout.write(line)
        aout.close()

#     try: logging.info('Exited fast_a_matrix()')
#     except: pass
#     print('\t\t\tExited pyp_nrm.fast_a_matrix()')
    if debug:
        print(a)
    return a

##
# Form a relationship matrix from a pedigree.  fast_a_matrix_r() differs from fast_a_matrix() in that the
# coefficients of relationship are corrected for the inbreeding of the parents.
# @param pedigree A PyPedal pedigree.
# @param pedopts PyPedal options.
# @param save Flag to indicate whether or not the relationship matrix is written to a file.
# @param method Use dense or sparse matrix storage.
# @retval A relationship as Numarray matrix.
def fast_a_matrix_r(pedigree, pedopts, save=0, method='dense'):
    """
    Form a relationship matrix from a pedigree.  fast_a_matrix_r() differs from
    fast_a_matrix() in that the coefficients of relationship are corrected for the
    inbreeding of the parents.
    """
    #try: logging.info('Entered fast_a_matrix_r()')
    #except: pass
    import math   # We need it for sqrt()
    animals = []
    sires = []
    dams = []
    #print(pedigree)
    l = len(pedigree)
    #method = 'dense'
    if method not in ['dense','sparse']:
        method = 'dense'
    try:
        # Use PySparse to provide sparse matrix storage for large
        # relationship matrices.
        if method == 'sparse':
            try:
                a = spmatrix.ll_mat_sym(l*l)
                for i in xrange(l):
                    a[i,i] = 1.
            except:
                #a = numpy.zeros([l,l],'d')
	        try:
                    a = numpy.zeros([l,l],'d')  # initialize a matrix of zeros of appropriate size
	        except MemoryError:
	            a = numpy.memmap('fast_a_matrix_r_mmap.bin', dtype='float32', mode='w+', shape=(l,l))
	        except:
	            print('[ERROR]: Unable to allocate a matrix of rank %s in pyp_nrm.fast_a_matrix_r()!' % ( l ))
	            logging.error('Unable to allocate a matrix of rank %s in pyp_nrm.fast_a_matrix_r()!', l)
                    return False

        # Otherwise, use NumPy and its dense matrices
        else:
            #a = numpy.zeros([l,l],'d')
	    try:
                a = numpy.zeros([l,l],'d')  # initialize a matrix of zeros of appropriate size
	    except MemoryError:
	        a = numpy.memmap('fast_a_matrix_r_mmap.bin', dtype='float32', mode='w+', shape=(l,l))
	    except:
	        print('[ERROR]: Unable to allocate a matrix of rank %s in pyp_nrm.fast_a_matrix_r()!' % ( l ))
	        logging.error('Unable to allocate a matrix of rank %s in pyp_nrm.fast_a_matrix_r()!', l)
                return False

        for i in xrange(l):
            animals.append(int(pedigree[i].animalID))
            sires.append(int(pedigree[i].sireID))
            dams.append(int(pedigree[i].damID))
        # Poorly-written code -- loops twice.  Once to compute CoI and a second time to
        # correct CoR for parental inbreeding.
        for row in xrange(l):
            for col in xrange(row,l):
                #print('[DEBUG]: row = %s\tcol = %s' % (row,col))
                if str(sires[col]) == str(pedopts['missing_parent']) and str(dams[col]) == str(pedopts['missing_parent']):
                    if row == col:
                        # both parents unknown and assumed unrelated
                        a[row,col] = 1.
                elif str(sires[col]) == str(pedopts['missing_parent']):
                    # sire unknown, dam known
                    if row == col:
                        a[row,col] = 1.
                    else:
                        a[row,col] = 0.5 * a[row,dams[col]-1]
                        a[col,row] = a[row,col]
                elif str(dams[col]) == str(pedopts['missing_parent']):
                    # sire known, dam unknown
                    if row == col:
                        a[row,col] = 1.
                    else:
                        a[row,col] = 0.5 * a[row,sires[col]-1]
                        a[col,row] = a[row,col]
                elif sires[col] != pedopts['missing_parent'] and dams[col] != pedopts['missing_parent']:
                    # both parents known
                    if row == col:
                        a[row,col] = 1. + ( 0.5 * a[sires[col]-1,dams[col]-1] )
                    else:
                        intermediate = a[row,sires[col]-1] + a[row,dams[col]-1]
                        a[row,col] = 0.5 * intermediate
                        a[col,row] = a[row,col]
                else:
                    if pedopts['debug_messages'] and pedopts['messages'] != 'quiet':
                        print('[ERROR]: There is a problem with the sire (ID %s) and/or dam (ID %s) of animal %s' % (pedigree[col].sireID, pedigree[col].damID, pedigree[col].animalID))
                    break
        for row in xrange(l):
            for col in xrange(row,l):
                if str(sires[col]) == str(pedopts['missing_parent']) and str(dams[col]) == str(pedopts['missing_parent']):
                    pass
                elif str(sires[col]) == str(pedopts['missing_parent']):
                    # sire unknown, dam known
                    if row != col and a[row,col] > 0.:
                        numerator = 0.5 * a[row,dams[col]-1]
                        denominator = sqrt ( a[dams[col]-1,dams[col]-1] )
                        try:
                            coefficient = numerator / denominator
                        except:
                            coefficient = 0.
                        a[row,col] = coefficient
                        a[col,row] = a[row,col]
                elif str(dams[col]) == str(pedopts['missing_parent']):
                    # sire known, dam unknown
                    if row != col and a[row,col] > 0.:
                        numerator = 0.5 * a[row,sires[col]-1]
                        denominator = sqrt ( a[sires[col]-1,sires[col]-1] )
                        try:
                            coefficient = numerator / denominator
                        except:
                            coefficient = 0.
                        a[row,col] = coefficient
                        a[col,row] = a[row,col]
                elif sires[col] != pedopts['missing_parent'] and dams[col] != pedopts['missing_parent']:
                    # both parents known
                    if row != col and a[row,col] > 0.:
                        numerator = 0.5 * ( a[row,sires[col]-1] + a[row,dams[col]-1] )
                        denominator = math.sqrt ( ( a[sires[col]-1,sires[col]-1] ) * ( a[dams[col]-1,dams[col]-1] ) )
                        try:
                            coefficient = numerator / denominator
                        except:
                            coefficient = 0.
                        a[row,col] = coefficient
                        a[col,row] = a[row,col]
                else:
                    pass
    except:
        a = numpy.zeros([1,1],'d')
    # print(a)
    if save == 1:
        a_outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_a_matrix_r_','.dat')
        aout = open(a_outputfile,'w')
        label = 'Produced by pyp_nrm/fast_a_matrix_r()\n'
        aout.write(label)
        for row in xrange(l):
            line = ''
            for col in xrange(l):
                if col == 0:
                    line = '%7.5f' % (a[row,col])
                else:
                    line = '%s%s%s' % (line,',',a[row,col])
            line = '%s%s' % (line,'\n')
            aout.write(line)
        aout.close()

    #try: logging.info('Exited fast_a_matrix_r()')
    #except: pass
    return a

##
# inbreeding() is a proxy function used to dispatch pedigrees to the appropriate
# function for computing CoI.  By default, small pedigrees < 10,000 animals) are
# processed with the tabular method directly.  For larger pedigrees, or if requested,
# the recursive method of VanRaden (1992) is used.
# @param pedobj A PyPedal pedigree object.
# @param method Keyword indicating which method of computing CoI should be used (tabular|vanraden).
# @param gens The number of generations from the pedigree to be used for calculating CoI.  By default, gens=0, which uses the complete pedigree.
# @param rels Flag indicating whether or not summary statistics should be computed for coefficients of relationship.
# @param output Flag indicating whether or not output files should be written.
# @param amethod The method parameter used by Aguilar's INBUPGF90 program.
# @param force Flag to override use of NRM attached to pedigree for finding COI (0: use NRM, 1: ignore NRM)
# @retval A dictionary of CoI keyed to renumbered animal IDs.
#@profile
def inbreeding(pedobj, method='tabular', gens=0, rels=0, output=1, force=0, amethod=3):
    """
    inbreeding() is a proxy function used to dispatch pedigrees to the appropriate
    function for computing CoI.  By default, small pedigrees < 10,000 animals) are
    processed with the tabular method directly.  For larger pedigrees, or if requested,
    the recursive method of VanRaden (1992) is used.
    """
    try: logging.info('Entered inbreeding()')
    except: pass
    fx = {}
    metadata = {}
    if method not in ['vanraden','tabular','meu_luo', 'mod_meu_luo', 'aguilar']:
        try: logging.warning('You passed an unrecognized method, %s, to pyp_nrm/inbreeding(); the method was changed to the default of \'tabular\'.', method)
        except: pass
        method = 'tabular'
    if int(gens) < 0:
        try: logging.warning('You passed an invalid value of gens, %s, to pyp_nrm/inbreeding(); gens was changed to the default of 0.', gens)
        except: pass
        gens = 0

    if rels:
        rel_dict = {}
        #rel_dict['r_count'] = 0
        rel_dict['r_count'] = (pedobj.metadata.num_records*(pedobj.metadata.num_records+1))/2
        rel_dict['r_nonzero_count'] = 0
        rel_dict['r_min'] = 0.
        rel_dict['r_max'] = 0.
        rel_dict['r_rng'] = 0.
        rel_dict['r_avg'] = 0.
        rel_dict['r_nonzero_avg'] = 0.
        rel_dict['r_sum'] = 0.
        rel_dict['r_nonzero_sum'] = 0.

        reldict = {}
        reldict['r_count'] = 0
        reldict['r_nonzero_count'] = 0
        reldict['r_nonzero_sum'] = 0.
        reldict['r_max'] = 0.
        reldict['r_min'] = 1.
        reldict['r_sum'] = 0.

    # If the user has already computed the NRM for this pedigree we will use it
    # to get the COI rather than recompute it all over again.
    if pedobj.kw['form_nrm'] and pedobj.nrm.nrm.shape[0] == pedobj.metadata.num_records and force == 0:
        try: logging.info('pyp_nrm/inbreeding() is using the NRM attached to the pedigree to get COI.')
        except: pass
        for _i in xrange(pedobj.metadata.num_records):
            fx[pedobj.pedigree[_i].animalID] = pedobj.nrm.nrm[_i][_i] - 1.
        if rels == 1:
            n = pedobj.nrm.nrm.shape[0]
            reldict['r_count'] = ( n * ( n + 1 ) ) / 2
            for i in xrange(n):
                for j in xrange (i, n):
                    if i != j:
                        if pedobj.nrm.nrm[i][j] > 0.:
                            reldict['r_nonzero_count'] = \
                                reldict['r_nonzero_count'] + 1
                            reldict['r_nonzero_sum'] = reldict['r_nonzero_sum'] + \
                                pedobj.nrm.nrm[i][j]
                            if pedobj.nrm.nrm[i][j] > reldict['r_max']:
                                reldict['r_max'] = pedobj.nrm.nrm[i][j]
                            if pedobj.nrm.nrm[i][j] < reldict['r_min']:
                                reldict['r_min'] = pedobj.nrm.nrm[i][j]
                        reldict['r_sum'] = reldict['r_sum'] + pedobj.nrm.nrm[i][j]
            #print('[DEBUG]: reldict: ', reldict)
    else:
        if method == 'vanraden':# or pedobj.metadata.num_records > 1000:
            #if pedobj.metadata.num_records > 1000:
            #    try: logging.warning('pyp_nrm.inbreeding() dispatched the pedigree %s to pyp_nrm/inbreeding_vanraden() because it contains more than 1000 records.', pedobj.kw['pedname'])
            #    except: pass
            if rels:
                fx, reldict = inbreeding_vanraden(pedobj, gens=gens, rels=rels)
                #print('[DEBUG]: reldict: ', reldict)
            else:
                fx = inbreeding_vanraden(pedobj, gens=gens)
        elif method == 'meu_luo':
            if rels != 0:
                logging.warning('You asked pyp_nrm.inbreeding() to compute relationships as well as inbreeding, but requested method %s, which does not provide coefficients of relationship. Only coefficients of inbreeding will be returned.', method)
            fx = inbreeding_meuwissen_luo(pedobj, gens=gens, rels=rels)
        elif method == 'mod_meu_luo':
            if rels != 0:
                logging.warning('You asked pyp_nrm.inbreeding() to compute relationships as well as inbreeding, but requested method %s, which does not provide coefficients of relationship. Only coefficients of inbreeding will be returned.', method)
            fx = inbreeding_modified_meuwissen_luo(pedobj, gens=gens, rels=rels)
        elif method == 'aguilar':
            logging.info('Using the INBUPGF90 program to compute inbreeding with method ', amethod)
            fx = inbreeding_aguilar(pedobj, amethod)
        else:
            if rels:
                fx, reldict = inbreeding_tabular(pedobj, gens=gens, rels=rels)
                #print('[DEBUG]: reldict: ', reldict)
            else:
                fx = inbreeding_tabular(pedobj, gens=gens)
    # Write summary stats to a file.
    if output:
        a_outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_inbreeding','.dat')
        aout = open(a_outputfile,'w')
        aout.write('# Inbreeding coefficients\n')
        # If the pedigree uses names do the same in the output file because
        # the original IDs won't mean anything when compared to the input
        # pedugree.
        if 'ASD' in pedobj.kw['pedformat']:
            aout.write('# Name\tRenum ID\tf_x\n')
        else:
            aout.write('# Orig ID\tRenum ID\tf_x\n')
    f_sum = 0.0
    f_min = 999.0
    f_max = -999.0
    f_nonzero_sum = 0.0
    f_nonzero_min = 999.0
    f_nonzero_max = -999.0
    f_nonzero_count = 0
    for k, v in fx.iteritems():
        if output:
            if 'ASD' in pedobj.kw['pedformat']:
                aout.write('%s\t%s\t%s\n'%(pedobj.pedigree[int(k)-1].name,k,v))
            else:
                aout.write('%s\t%s\t%s\n'%(pedobj.pedigree[int(k)-1].originalID,k,v))
        # Update self.fa for each Animal object in the pedigree.
        pedobj.pedigree[int(k)-1].fa = v
        f_sum = f_sum + v
        if v > f_max:
            f_max = v
        if v < f_min:
            f_min = v
        if v > 0.:
            f_nonzero_count = f_nonzero_count + 1
            f_nonzero_sum = f_nonzero_sum + v
            if v > f_nonzero_max:
                f_nonzero_max = v
            if v < f_nonzero_min:
                f_nonzero_min = v
    if len(fx.keys()) == 0:
        f_rng = 0.
        f_avg = 0.
        f_min = 0.
        f_max = 0.
    else:
        f_avg = f_sum / len(fx.keys())
        f_rng = f_max - f_min
    # If there are no inbred animals in the pedigree we need
    # to make sure that meaningful values are returned, rather
    # that, e.g., -999.
    if f_nonzero_count == 0:
        f_nonzero_rng = 0.
        f_nonzero_avg = 0.
        f_nonzero_min = 0.
        f_nonzero_max = 0.
    else:
        f_nonzero_rng = f_nonzero_max - f_nonzero_min
        f_nonzero_avg = f_nonzero_sum / f_nonzero_count
    # Summary statistics including all CoI
    metadata['all'] = {}
    metadata['all']['f_count'] = len(fx.keys())
    metadata['all']['f_sum'] = f_sum
    metadata['all']['f_min'] = f_min
    metadata['all']['f_max'] = f_max
    metadata['all']['f_rng'] = f_rng
    metadata['all']['f_avg'] = f_avg
    # Summary statistics including only nonzero CoI
    metadata['nonzero'] = {}
    metadata['nonzero']['f_count'] = f_nonzero_count
    metadata['nonzero']['f_sum'] = f_nonzero_sum
    metadata['nonzero']['f_min'] = f_nonzero_min
    metadata['nonzero']['f_max'] = f_nonzero_max
    metadata['nonzero']['f_rng'] = f_nonzero_rng
    metadata['nonzero']['f_avg'] = f_nonzero_avg

    if rels:
	if pedobj.kw['debug_messages'] == 1: print('[DEBUG]: reldict: ', reldict)
        if reldict['r_count'] > 0:
            if reldict['r_min'] < rel_dict['r_min']:
                rel_dict['r_min'] = reldict['r_min']
            if reldict['r_max'] > rel_dict['r_max']:
                rel_dict['r_max'] = reldict['r_max']
            rel_dict['r_rng'] = rel_dict['r_max'] - rel_dict['r_min']
            rel_dict['r_sum'] = reldict['r_sum'] + reldict['r_sum']
            rel_dict['r_avg'] = rel_dict['r_sum'] / rel_dict['r_count']
        if reldict['r_nonzero_count'] > 0:
            rel_dict['r_nonzero_count'] = reldict['r_nonzero_count'] + rel_dict['r_nonzero_count']
            rel_dict['r_nonzero_sum'] = rel_dict['r_nonzero_sum'] + reldict['r_nonzero_sum']
            rel_dict['r_nonzero_avg'] = rel_dict['r_nonzero_sum'] / rel_dict['r_nonzero_count']
        #print('\n\nrel_dct: %s\n\n' % ( rel_dict ))

    if output:
        line = '='*80
        aout.write('%s\n' % line)
        aout.write('Inbreeding Statistics\n')
        line = '-'*80
        aout.write('All animals:\n')
        aout.write('%s\n' % line)
        aout.write('\tCount:\t%s\n'%len(fx.keys()))
        aout.write('\tMean:\t%s\n'%f_avg)
        aout.write('\tMin:\t%s\n'%f_min)
        aout.write('\tMax:\t%s\n'%f_max)
        line = '-'*80
        aout.write('Animals with non-zero CoI:\n')
        aout.write('%s\n' % line)
        aout.write('\tCount:\t%s\n'%f_nonzero_count)
        aout.write('\tMean:\t%s\n'%f_nonzero_avg)
        aout.write('\tMin:\t%s\n'%f_nonzero_min)
        aout.write('\tMax:\t%s\n'%f_nonzero_max)
        aout.close()

    try: logging.info('Exited inbreeding()')
    except: pass

    pedobj.kw['f_computed'] = 1
    out_dict = {}
    out_dict['metadata'] = metadata
    out_dict['fx'] = fx
    if rels:
        return out_dict, rel_dict
    else:
        return out_dict

##
# inbreeding_vanraden() uses VanRaden's (1992) method for computing coefficients of
# inbreeding in a large pedigree.  The method works as follows:
#   1.  Take a large pedigree and order it from youngest animal to oldest (n, n-1, ..., 1);
#   2.  Recurse through the pedigree to find all of the ancestors of that animal n;
#   3.  Reorder and renumber that "subpedigree";
#   4.  Compute coefficients of inbreeding for that "subpedigree" using the tabular
#       method (Emik and Terrill, 1949);
#   5.  Put the coefficients of inbreeding in a dictionary;
#   6.  Repeat 2 - 5 for animals n-1 through 1; the process is slowest for the early
#       pedigrees and fastest for the later pedigrees.
# @param pedobj A PyPedal pedigree object.
# @param cleanmaps Flag to denote whether or not subpedigree ID maps should be deleted after they are used (0|1).
# @param gens The number of generations from the pedigree to be used for calculating CoI.  By default, gens=0, which uses the complete pedigree.
# @param rels Flag indicating whether or not summary statistics should be computed for coefficients of relationship.
# @retval A dictionary of CoI keyed to renumbered animal IDs
# @profile
def inbreeding_vanraden(pedobj, cleanmaps=1, gens=0, rels=0):
    """
    inbreeding_vanraden() uses VanRaden's (1992) method for computing coefficients of
    inbreeding in a large pedigree.  The method works as follows:
    1.  Take a large pedigree and order it from youngest animal to oldest (n, n-1, ...,
        1);
    2.  Recurse through the pedigree to find all of the ancestors of that animal n;
    3.  Reorder and renumber that "subpedigree";
    4.  Compute coefficients of inbreeding for that "subpedigree" using the tabular
        method (Emik and Terrill, 1949);
    5.  Put the coefficients of inbreeding in a dictionary;
    6.  Repeat 2 - 5 for animals n-1 through 1; the process is slowest for the early
        pedigrees and fastest for the later pedigrees.
    """
    try: logging.info('Entered inbreeding_vanraden()')
    except: pass

    #print('Entered inbreeding_vanraden()')
    #print('\tConverting pedigree to graph.')
    from PyPedal import pyp_network
    # Using pyp_network.ped_to_graph() instead of pyp_nrm.recurse_pedigree()
    # provides a gain in performance of at least an order of magnitude.
    ng = pyp_network.ped_to_graph(pedobj)

    _ped = []       # This is a temporary pedigree
    top_ped = []

    if int(gens) > 0:
#         print('gens: %s' % ( gens ))
        top_peddict = pyp_network.find_ancestors_g(ng, len(pedobj.idmap), {}, gens)
        top_peddict[len(pedobj.idmap)] = 1
        top_ped = top_peddict.keys()
        #print('top_ped: ', top_ped)
        top_r = []
        _anids = []
        for _j in top_ped:
            # Make sure that we only include animals from the user-
            # specified generations.
            if top_peddict[_j] <= gens:
                top_r.append(copy.copy(pedobj.pedigree[int(_j)-1]))
                # Animals in the earliest generation need to have
                # their sire and dam IDs set to unknown.
                if top_peddict[_j] == 1:
                    top_r[-1].sireID = 0
                    top_r[-1].damID = 0
                _anids.append(top_r[-1].animalID)
    else:
        _anids = pedobj.backmap.keys()     # Distinct animal IDs in the pedigree

    fx = {}         # This will hold our coefficients of inbreeding
    _parents = {}   # Stores a list of sire-dam pairs along with the youngest offspring
                    # of that pair.  Used as a lookup table to avoid lots of redundant
                    # calculations for full-sibs.
    _anids.sort()       # sort from oldest to youngest
    _anids.reverse()    # reverse the list to put the youngest animals first
    _counter = 0
    _cum_f_counter = 0
    _vanraden_round = 0
    _cum_pct_proc = 0.
    _related = {}   # Dictionary for looking-up animals with non-zero
                    # relationships.

    # If the user wants summary stats on coefficients of relationship,
    # prepare the dictionary, counters, and accumulators.
    if rels:
        reldict = {}
        reldict['r_count'] = 0
        reldict['r_nonzero_count'] = 0
        reldict['r_nonzero_sum'] = 0.
        reldict['r_max'] = 0.
        reldict['r_min'] = 1.
        reldict['r_sum'] = 0.
        #print('[DEBUG]: inbreeding_vanraden(): reldict: ', reldict)
    #print('_anids: ', _anids)

    for i in _anids:
        #print('\t_anid: ', i)
        if int(gens) == 0:
            _parent_key = '%s_%s' % ( pedobj.pedigree[int(i)-1].sireID,
                pedobj.pedigree[int(i)-1].damID )
        else:
            _parent_key = '%s_%s' % ( top_r[top_peddict[int(i)]].sireID, top_r[top_peddict[int(i)]].damID )
        try:
            _k = fx[i]  # If an exception is thrown, an animal is not in the
                        # dictionary yet.
        except KeyError:
            try:
                # If the parental combination has been seen already then we
                # already know the COI of the mating and do not need to do
                # the calculations again.
                #print('_parent_key: ', _parent_key)
                #print('fx[_parents[_parent_key]]: ', fx[_parents[_parent_key]])
                fx[i] = fx[_parents[_parent_key]]
            except:
                #print('\tAnimal %s not in inbreeding dict' % ( i ))
                _f_counter = 0
                _vanraden_round = _vanraden_round + 1
                if _vanraden_round == 1:
                    try: logging.info('Starting round %s of pyp_nrm/inbreeding_vanraden().', _vanraden_round)
                    except: pass
                _tag = '%s_%s' % (pedobj.kw['filetag'],i)
                if int(gens) > 0:
                    _ped = top_peddict
                else:
                    _ped = pyp_network.find_ancestors(ng, i, [])
                    _ped.append(i)
                #print('ped len: %d' % ( len(_ped) ))
                #print('ped: %s' % ( _ped ))
                if int(gens) > 0:
                    _r = top_r
                else:
                    _r = []     # This list will hold a copy of the objects in _ped
                                # so that we can renumber animal i's pedigree without
                                # changing the data in pedobj.pedigree.
                    _map = {}
                    for j in _ped:
                        # This is VERY important -- rather than append a reference
                        # to _ped[j-1] to _r we need to append a COPY of _ped[j-1]
                        # to _r.  If you change this code and get rid of the call to
                        # copy.copy() then things will not work correctly.  You will
                        # realize what you have done when your renumberings seem to
                        # be spammed.
                        _r.append(copy.copy(pedobj.pedigree[int(j)-1]))
                # We also need to honor the slow_reorder option.
                if pedobj.kw['slow_reorder']:
                    _r = pyp_utils.reorder(_r,_tag)      # Reorder the pedigree
                else:
                    _r = pyp_utils.fast_reorder(_r,_tag)      # Reorder the pedigree
                _s, _map = pyp_utils.renumber(_r,_tag, returnmap=1, debug=pedobj.kw['debug_messages'],animaltype=pedobj.kw['animal_type'])  # Renumber the pedigree
                # _map maps IDs from original IDs to renumbered IDs.
                # _backmap allows renumbered ID => original ID reverse lookups.
                _backmap = {}
                for _mk, _mv in _map.iteritems():
                    _backmap[_mv] = _mk
                # There is a potential error lurking here!  The filetag passed to
                # fast_a_matrix as "_tag" is expected to be a pedoptions dictionary.
                # Hm...I think that passing a copy of the kw dictionary from pedobj
                # with the filetag changed as appropriate will do the trick.
                _opts = copy.copy(pedobj.kw)
                _opts['filetag'] = _tag
                # We need to accomodate the 'nrm_method' option, too.  The need for
                # this is clearly demonstrated by horse.ped in the examples/ subdirectory -
                # the inbreeding is so intense in that pedigree that four of the animals
                # have r_xy >= 1. if we do not adjust the elements of A for parental in-
                # breeding.
                if pedobj.kw['nrm_method'] == 'nrm':
                    _a = fast_a_matrix(_s,_opts,method=pedobj.kw['matrix_type'])     # Form the NRM w/the tabular method
                else:
                    _a = fast_a_matrix_r(_s,_opts,method=pedobj.kw['matrix_type'])
                #print('len(_ped): ', len(_ped))
                for j in xrange(len(_ped)):
                    _orig_id = _backmap[_s[j].animalID]
                    # The same animal can appear in many different pedigrees, but
                    # it should always have the same CoI.  We don't want to waste
                    # a lot of time writing the same value to a dictionary many
                    # times, so we are going to check and see if the animal already
                    # has an entry in the CoI dictionary.  It it does, do not write
                    # to it again.
                    try:
                        _check = fx[_orig_id]
                    except KeyError:
                        #print(_a)
                        fx[_orig_id] = _a[j][j] - 1.
                        _f_counter = _f_counter + 1
                    if rels:
                        for k in xrange(j, len(_ped)):
                            if j != k:
                                _rxykey = '%s_%s' % (_backmap[_s[j].animalID], \
                                    _backmap[_s[k].animalID])
                                if _a[j][k] > 0.:
                                    try:
                                        _rxy = _related[_rxykey]
                                    except KeyError:
                                        _related[_rxykey] = _a[j][k]
                                    reldict['r_nonzero_count'] = \
                                        reldict['r_nonzero_count'] + 1
                                    reldict['r_nonzero_sum'] = \
                                        reldict['r_nonzero_sum'] + _a[j][k]
                                    if _a[j][k] > reldict['r_max']:
                                        reldict['r_max'] = _a[j][k]
                                    if _a[j][k] < reldict['r_min']:
                                        reldict['r_min'] = _a[j][k]
                                reldict['r_count'] = reldict['r_count'] + 1
                                reldict['r_sum'] = reldict['r_sum'] + _a[j][k]
                    try:
                        _ptest = _parents[_parent_key]
                    except KeyError:
                        #_parents[_parent_key] = _orig_id
                        _parents[_parent_key] = _parent_key
                    # We only got into this loop because this combination of  parents
                    # did not have an entry in fx, so put one there.
                    fx[_parent_key] = fx[_orig_id]
                    #print('Parent combination: ', _parent_key, '\t\tf: ', fx[_parent_key])
                #print('\n\n', reldict, '\n\n')
                #print('_related (%d): %s' % ( len(_related), _related ))

                if cleanmaps:               # Clean up the subpedigree ID maps that we are
                    pyp_utils.delete_id_map(_tag)     # not going to use again.
                _map = {}               # Empty our working dictionary and lists
                _a = []
                _s = []
                _r = []
                _ped = []
                _pct_proc = float(_f_counter) / float(pedobj.metadata.num_records)
                _cum_pct_proc = _cum_pct_proc + _pct_proc
                _cum_f_counter = _cum_f_counter + _f_counter
    #             if _pct_proc > 0.01:
    #            if pedobj.kw['messages'] == 'verbose':
                    #print('%s of animals processed in round %s of #pyp_nrm/inbreeding_vanraden().' % (_pct_proc,_vanraden_round))
                    #print('%s of all animals have been processed in #pyp_nrm/inbreeding_vanraden().' % (_cum_pct_proc))
                #try: logging.info('%s of animals processed in round %s of #pyp_nrm/inbreeding_vanraden().', _pct_proc, _vanraden_round)
                #except: pass
                logging.info('%s pct (%s) of all animals have been processed in pyp_nrm/inbreeding_vanraden().', _cum_pct_proc, _cum_f_counter)
            _counter = _counter + 1

    # Clean-up the parent-combination entries from fx
    for k in _parents.keys():
        del fx[k]

    try: logging.info('Exited inbreeding_vanraden()')
    except: pass
#     print(fx)
    if rels:
        return fx, reldict
    else:
        return fx

##
# inbreeding_aguilar() uses use Ignacio Aguilar's INBUPGF90 programto compute coefficients of
# inbreeding in large pedigrees.
# @param pedobj A PyPedal pedigree object.
# @retval A dictionary of CoI keyed to renumbered animal IDs
# @profile
def inbreeding_aguilar(pedobj, amethod=3):
    """
    inbreeding_aguilar() uses use Ignacio Aguilar's INBUPGF90 program to compute coefficients of
    inbreeding in large pedigrees.
    :param pedobj:
    :param amethod:
    :return:
    """
    # Before we do anything let's see if the INBUPGF90 executable is visible to PyPedal.
    if not pyp_utils.which('inbupgf90') and not pyp_utils.which('inbupgf90.exe'):
        if pedobj.kw['messages'] == 'verbose':
            print('\t[inbreeding_aguilar]: Cannot find the INBUPGF90 executable tat %s!' % \
                  datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
            logging.error('[inbreeding_aguilar]: cannot find the INBUPGF90 executable at %s!',
                          datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
        return None

    # Per an e-mail from Ignacio Aguilar on 06/25/2014, INBUPGF90 does NOT emit a proper
    # status return code when it exits, which makes it tricky to know for sure when the
    # job is done. I've observed a number of cases where the simulation appears to stall
    # because subprocess.call() does not recognize that INBUPGF90 has finished a job. So,
    # I've cobbled-together a solution using ideas from Ezequiel Nicolazzi
    # (https://github.com/nicolazzie/AffyPipe/blob/master/AffyPipe.py) and a post on
    # Stack Overflow (http://stackoverflow.com/questions/12057794/
    # python-using-popen-poll-on-background-process). I'm not 100% sure that this works
    # as intended, but I'm out of ideas.
    logfile = '%s_aguilar.log' % pedobj.kw['pedname']
    # Several methods can be used:
    # 1 - recursive as in Aguilar & Misztal, 2008 (default)
    # 2 - recursive but with coefficients store in memory, faster with large number of
    #     generations but more memory requirements
    # 3 - method as in Meuwissen & Luo 1992
    if amethod not in [1, 2, 3]: amethod = 3
    if pedobj.kw['messages'] == 'verbose':
        print('\t[inbreeding_aguilar]: Started inbupgf90 to calculate COI at %s' % datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    ###
    # 1. Put some code here to check for the inbupgf90 binary in the user's path.
    # ...
    # 2. Put in code to make the pedigree flatfile that INBUPGF90 needs
    pedfile = 'aguilar_pedigree_%s.txt' % pedobj.kw['pedname']
    callinbupgf90 = ['inbupgf90', '--pedfile', pedfile, '--method', '3', '--yob', '>', logfile, '2>&1&']
    time_waited = 0
    p = subprocess.Popen(callinbupgf90, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while p.poll() is None:
        # Wait 1 second between pokes with a sharp stick.
        time.sleep(10)
        time_waited += 10
        p.poll()
        if time_waited % 60 == 0 and pedobj.kw['messages'] == 'verbose':
            print('\t\t[inbreeding_aguilar]: Waiting for INBUPGF90 to finish -- %s minutes so far...' % int(time_waited/60))
            logging.info('[inbreeding_aguilar]: Waiting for INBUPGF90 to finish -- %s minutes so far...', int(time_waited/60))

    # Pick-up the output from INBUPGF90
    (results, errors) = p.communicate()
    if errors == '':
        if pedobj.kw['messages'] == 'verbose':
            print('\t\t[inbreeding_aguilar]: INBUPGF90 finished without problems at %s!' % \
                  datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
            logging.info('[inbreeding_aguilar]: INBUPGF90 finished without problems at %s!',
                datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    else:
        if pedobj.kw['messages'] == 'verbose':
            print('\t\t[inbreeding_aguilar]: INBUPGF90 finished with errors: %s' % errors)
        logging.error('[inbreeding_aguilar]: INBUPGF90 finished with errors: ', errors)
    if pedobj.kw['messages'] == 'verbose':
        print('\t[aguilar_inbreeding]: Finished using INBUPGF90 to calculate COI at %s' % \
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    logging.info('[inbreeding_aguilar]: Finished using INBUPGF90 to calculate COI at %s',
        datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))

    # Load the COI into a dictionary keyed by original animal ID
    coifile = '%s.solinb' % pedfile
    if pedobj.kw['messages'] == 'verbose':
        print('\t[inbreeding_aguilar]: Putting coefficients of inbreeding from %s.solinb in a dictionary at %s' \
            % (pedfile, datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
    logging.info('[inbreeding_aguilar]: Putting coefficients of inbreeding from %s.solinb in a dictionary at %s' % \
            (pedfile, datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
    inbr = {}
    ifh = open(coifile, 'r')
    for line in ifh:
        pieces = line.split()
        inbr[pieces[0]] = float(pieces[1])
    ifh.close()

    # Now, assign the coefficients of inbreeding to the animal records
#    for c in cows: c[10] = inbr[c[0]]
#    for dc in dead_cows: dc[10] = inbr[dc[0]]
#    for b in bulls: b[10] = inbr[b[0]]
#    for db in dead_bulls: db[10] = inbr[db[0]]

    # Clean-up
    os.remove(pedfile)
    os.remove('%.solinb') % (pedfile)
    os.remove('%s.errors') % (pedfile)
    os.remove('%.inbavgs') % (pedfile)

    # Send back the coefficients of inbreeding to the caller.
    return inbr

##
# recurse_pedigree() performs the recursion needed to build the subpedigrees used by
# inbreeding_vanraden().  For the animal with animalID anid recurse_pedigree() will
# recurse through the pedigree myped and add references to the relatives of anid to
# the temporary pedigree, _ped.
# @param pedobj A PyPedal pedigree.
# @param anid The ID of the animal whose relatives are being located.
# @param _ped A temporary PyPedal pedigree that stores references to relatives of anid.
# @retval A list of references to the relatives of anid contained in myped.
def recurse_pedigree(pedobj, anid, _ped):
    """
    recurse_pedigree() performs the recursion needed to build the subpedigrees used by
    inbreeding_vanraden().  For the animal with animalID anid recurse_pedigree() will
    recurse through the pedigree myped and add references to the relatives of anid to
    the temporary pedigree, _ped.
    """
    try:
        anid = int(anid)
        if anid != 0:
            if pedobj.pedigree[anid-1] not in _ped:
                _ped.append(pedobj.pedigree[anid-1])
        _sire = pedobj.pedigree[anid-1].sireID
        _dam = pedobj.pedigree[anid-1].damID
        if _sire != pedobj.kw['missing_parent']:
            recurse_pedigree(pedobj,_sire,_ped)
        if _dam != pedobj.kw['missing_parent']:
            recurse_pedigree(pedobj,_dam,_ped)
    except:
        pass
    return _ped

##
# recurse_pedigree_n() recurses to build a pedigree of depth n.  A depth less than 1 returns
# the animal whose relatives were to be identified.
# @param pedobj A PyPedal pedigree.
# @param anid The ID of the animal whose relatives are being located.
# @param _ped A temporary PyPedal pedigree that stores references to relatives of anid.
# @param depth The depth of the pedigree to return.
# @retval A list of references to the relatives of anid contained in myped.
def recurse_pedigree_n(pedobj, anid, _ped, depth=3):
    """
    recurse_pedigree_n() recurses to build a pedigree of depth n.  A depth
    less than 1 returns the animal whose relatives were to be identified.
    """
    try:
        anid = int(anid)
        if anid != pedobj.kw['missing_parent']:
            if pedobj.pedigree[anid-1] not in _ped:
                _ped.append(pedobj.pedigree[anid-1])
        if depth > 0:
            _sire = pedobj.pedigree[anid-1].sireID
            _dam = pedobj.pedigree[anid-1].damID
            if _sire != pedobj.kw['missing_parent']:
                recurse_pedigree_n(pedobj,_sire,_ped,depth-1)
            if _dam != pedobj.kw['missing_parent']:
                recurse_pedigree_n(pedobj,_dam,_ped,depth-1)
    except:
        pass
    return _ped

##
# recurse_pedigree_onesided() recurses to build a subpedigree from either the sire
# or dam side of a pedigree.
# @param pedobj A PyPedal pedigree.
# @param side The side to build: 's' for sire and 'd' for dam.
# @param anid The ID of the animal whose relatives are being located.
# @param _ped A temporary PyPedal pedigree that stores references to relatives of anid.
# @retval A list of references to the relatives of anid contained in myped.
def recurse_pedigree_onesided(pedobj, anid, _ped, side):
    """
    recurse_pedigree_onesided() recurses to build a subpedigree from either the sire
    or dam side of a pedigree.
    """
    try:
        anid = int(anid)
        if anid != 0:
            if pedobj.pedigree[anid-1] not in _ped:
                _ped.append(pedobj.pedigree[anid-1])
        if side == 's':
            _sire = pedobj.pedigree[anid-1].sireID
            if _sire != pedobj.kw['missing_parent']:
                recurse_pedigree(pedobj,_sire,_ped)
        else:
            _dam = pedobj.pedigree[anid-1].damID
            if _dam != pedobj.kw['missing_parent']:
                recurse_pedigree(pedobj,_dam,_ped)
    except:
        pass
    return _ped

##
# recurse_pedigree_idonly() performs the recursion needed to build subpedigrees.
# @param pedobj A PyPedal pedigree.
# @param anid The ID of the animal whose relatives are being located.
# @param _ped A PyPedal list that stores the animalIDs of relatives of anid.
# @retval A list of animalIDs of the relatives of anid contained in myped.
def recurse_pedigree_idonly(pedobj, anid, _ped):
    """
    recurse_pedigree_idonly() performs the recursion needed to build subpedigrees.
    """
    try:
        anid = int(anid)
        if anid != 0:
            #print('\t\tanid %s\tsireid %s\tdamid %s' % ( anid, pedobj.pedigree[anid-1].sireID, pedobj.pedigree[anid-1].damID ))
            if pedobj.pedigree[anid-1].animalID not in _ped:
                _ped.append(pedobj.pedigree[anid-1].animalID)
        _sire = pedobj.pedigree[anid-1].sireID
        _dam = pedobj.pedigree[anid-1].damID
        if _sire != pedobj.kw['missing_parent']:
            recurse_pedigree_idonly(pedobj,_sire,_ped)
        if _dam != pedobj.kw['missing_parent']:
            recurse_pedigree_idonly(pedobj,_dam,_ped)
    except:
        pass
    return _ped

##
# recurse_pedigree_idonly_side() performs the recursion needed to build
# a subpedigree containing only animal IDs for either all sires or all
# dams.  That is, a pedigree would go sire-paternal grandsire-paternal
# great-grandsire, etc.
# @param pedobj A PyPedal pedigree.
# @param anid The ID of the animal whose relatives are being located.
# @param _ped A PyPedal list that stores the animalIDs of relatives of anid.
# @param side The side of the pedigree to follow ('s'|'d').
# @retval A list of animalIDs of the relatives of anid contained in myped.
def recurse_pedigree_idonly_side(pedobj, anid, _ped, side='s'):
    """
    recurse_pedigree_idonly_side() performs the recursion needed to build
    a subpedigree containing only animal IDs for either all sires or all
    dams.  That is, a pedigree would go sire-paternal grandsire-paternal
    great-grandsire, etc.
    """
    if side not in ['s','d']:
        side = 's'
    try:
        anid = int(anid)
        if anid != 0:
            if pedobj.pedigree[anid-1].animalID not in _ped:
                _ped.append(pedobj.pedigree[anid-1].animalID)
        _sire = pedobj.pedigree[anid-1].sireID
        _dam = pedobj.pedigree[anid-1].damID
        if side == 's':
            if _sire != pedobj.kw['missing_parent']:
                recurse_pedigree_idonly_side(pedobj,_sire,_ped,side='s')
        if side == 'd':
            if _dam != pedobj.kw['missing_parent']:
                recurse_pedigree_idonly_side(pedobj,_dam,_ped,side='d')
    except:
        pass
    return _ped

##
# inbreeding_tabular() computes CoI using the tabular method by calling
# fast_a_matrix() to form the NRM directly.  In order for this routine
# to return successfully requires that you are able to allocate a matrix
# of floats of dimension len(myped)**2.
# @param pedobj A PyPedal pedigree object.
# @param gens The number of generations from the pedigree to be used for calculating CoI.  By default, gens=0, which uses the complete pedigree.
# @param rels Flag indicating whether or not summary statistics should be computed for coefficients of relationship.
# @retval A dictionary of CoI keyed to renumbered animal IDs
def inbreeding_tabular(pedobj, gens=0, rels=0):
    """
    inbreeding_tabular() computes CoI using the tabular method by calling
    fast_a_matrix() to form the NRM directly.  In order for this routine
    to return successfully requires that you are able to allocate a matrix
    of floats of dimension len(myped)**2.
    """
    try: logging.info('Entered inbreeding_tabular()')
    except: pass

    # If the user wants summary stats on coefficients of relationship,
    # prepare the dictionary, counters, and accumulators.
    if rels:
        reldict = {}
        reldict['r_count'] = 0
        reldict['r_nonzero_count'] = 0
        reldict['r_nonzero_sum'] = 0.
        reldict['r_max'] = 0.
        reldict['r_min'] = 1.
        reldict['r_sum'] = 0.

    # See pyp_nrm.inbreeding_vanraden() for detailed notes on what
    # the code in this loop does.
    if int(gens) > 0:
        _ped = pyp_network.find_ancestors_g(ng, i, [], gens)
        _ped.append(i)
        _a, _s, _r = [], [], []
        _map = {}
        for j in _ped:
            _r.append(copy.copy(pedobj.pedigree[int(j)-1]))
        if pedobj.kw['slow_reorder']:
            _r = pyp_utils.reorder(_r,_tag)      # Reorder the pedigree
        else:
            _r = pyp_utils.fast_reorder(_r,_tag)      # Reorder the pedigree
        _s, _map = pyp_utils.renumber(_r,_tag, returnmap=1, debug=pedobj.kw['debug_messages'],animaltype=pedobj.kw['animal_type'])
        _backmap = {}
        for _mk, _mv in _map.iteritems():
            _backmap[_mv] = _mk
        _opts = copy.copy(pedobj.kw)
        _opts['filetag'] = _tag
        if pedobj.kw['nrm_method'] == 'nrm':
            _a = fast_a_matrix(pedobj.pedigree,pedobj.kw,method=pedobj.kw['matrix_type'])
        else:
            _a = fast_a_matrix_r(pedobj.pedigree,pedobj.kw,method=pedobj.kw['matrix_type'])
        fx = {}
        for i in xrange(len(_ped)):
            fx[pedobj.pedigree[i].animalID] = _a[i][i] - 1.
        del(_a)
    else:
        try:
            if pedobj.kw['nrm_method'] == 'nrm':
                _a = fast_a_matrix(pedobj.pedigree, pedobj.kw,method=pedobj. kw['matrix_type'])
            else:
                _a = fast_a_matrix_r(pedobj.pedigree, pedobj.kw, method=pedobj.kw['matrix_type'])
            fx = {}
            for i in xrange(pedobj.metadata.num_records):
                fx[pedobj.pedigree[i].animalID] = _a[i][i] - 1.
                if rels:
                    for j in xrange(i, pedobj.metadata.num_records):
                        if i != j:
                            if _a[i][j] > 0.:
                                reldict['r_nonzero_count'] = \
                                    reldict['r_nonzero_count'] + 1
                                if _a[i][j] > reldict['r_max']:
                                    reldict['r_max'] = _a[i][j]
                                if _a[i][j] < reldict['r_min']:
                                    reldict['r_min'] = _a[i][j]
                            reldict['r_count'] = reldict['r_count'] + 1
                            reldict['r_sum'] = reldict['r_sum'] + _a[i][j]
            del(_a)
        except:
            pass
    try: logging.info('Exited inbreeding_tabular()')
    except: pass
    if rels:
        return fx, reldict
    else:
        return fx

##
# inbreeding_meuwissen_luo() computes CoI using the method of Meuwissen and
# Luo (1992). It calculates only inbreeding coefficients, not relationships.
# @param pedobj A PyPedal pedigree object.
# @param gens The number of generations from the pedigree to be used for calculating CoI.  By default, gens=0, which uses the complete pedigree.
# @param rels Flag indicating whether or not summary statistics should be computed for coefficients of relationship.
# @retval A dictionary of CoI keyed to renumbered animal IDs
#@profile
def inbreeding_meuwissen_luo(pedobj, gens=0,**kw):
    """
    inbreeding_meuwissen_luo() computes CoI using the method of Meuwissen and
    Luo (1992). It calculates only inbreeding coefficients, not relationships.
    This code is a pretty direct implementation of the algorithm presented on
    pp. 311-312 of Meuiwissen, T.H.E., and Z. Luo. 1992. Computing inbreeding
    coefficients in large populations. Genet. Sel. Evol. 24:305-313.
    """
    try: logging.info('Entered inbreeding_meuwissen_luo()')
    except: pass

    # Setup dictionary to accumulate coefficients of inbreeding
    fx = {}
    for p in pedobj.pedigree:
        fx[p.animalID] = 0.0

    # First, try and allocate the vectors in RAM. If that does not work, try and allocate them using
    # memory-mapped files. If that does not work, well, give up.
    try:
        logging.info('Allocating vectors in pyp_nrm.inbreeding_meuwissen_luo().')
        lvec = numpy.zeros((len(pedobj.pedigree)),'d')
        #print('lvec:\t', lvec)
        avec = numpy.zeros((len(pedobj.pedigree)),'d')
        dvec = numpy.zeros((len(pedobj.pedigree)),'d')
    except MemoryError:
        logging.info('Unable to allocate a matrix of rank %s in RAM, trying to allocate a memory-mapped file, in pyp_nrm.inbreeding_meuwissen_luo()!', l)
        lvec = numpy.memmap('lvec_memmap.bin', dtype='float32', mode='w+', shape=(len(pedobj.pedigree)))
        lvec = 0.0
        avec = numpy.memmap('avec_memmap.bin', dtype='float32', mode='w+', shape=(len(pedobj.pedigree)))
        avec = 0.0
        dvec = numpy.memmap('dvec_memmap.bin', dtype='float32', mode='w+', shape=(len(pedobj.pedigree)))
        dvec = 0.0
    except:
        print('[ERROR]: Unable to allocate a matrix of rank %s in pyp_nrm.inbreeding_meuwissen_luo()!' % ( l ))
        logging.error('[ERROR]: Unable to allocate a matrix of rank %s in pyp_nrm.inbreeding_meuwissen_luo()!', l)
	return False

    if pedobj.kw['debug_messages']:
        print('[DEBUG]: Starting loop over pedigree with ', len(pedobj.pedigree), ' animals')
    for i in xrange(len(pedobj.pedigree)):
        if pedobj.kw['debug_messages']:
            print('\t[DEBUG]: Initializing local data structures for animal %s (idx: %s)' % ( pedobj.pedigree[i].animalID, i))
            print('\t\t[DEBUG]: a[%s] = %s' % ( i, avec[i] ))
        anc = []
	lvec[:] = 0.0
	lvec[i] = 1.0
        if pedobj.kw['debug_messages']:
	    print('\t\t[DEBUG]: l[%s] = %s' % ( i, lvec[i] ))
	# We're using 0-indexing, so the little F0 = -1 trick that M&L use does not help us here.
	# Check for the most common case first -- both parents unknown
	if pedobj.pedigree[i].sireID != pedobj.kw['missing_parent'] and pedobj.pedigree[i].damID != pedobj.kw['missing_parent']:
	    dvec[i] = 0.5 - ( 0.25 * ( fx[pedobj.pedigree[i].sireID] + fx[pedobj.pedigree[i].damID] ) )
	# Then check for both parents unknown
	elif pedobj.pedigree[i].sireID == pedobj.kw['missing_parent'] and pedobj.pedigree[i].damID == pedobj.kw['missing_parent']:
	    dvec[i] = 1.0
	# Finally, deal with either a known sire, unknown dam or known dam, unknown sire.
	else:
	    # The "-1.0" in the expressions below comes from M&L setting the coefficient of inbreeding for unknown parents to -1.0
	    # in order to ensure that the within-family variance is correct.
	    # Knwon sire, unknown dam
	    if pedobj.pedigree[i].sireID != pedobj.kw['missing_parent'] and pedobj.pedigree[i].damID == pedobj.kw['missing_parent']:
                dvec[i] = 0.5 - ( 0.25 * ( fx[pedobj.pedigree[i].sireID] - 1.0 ) )
	    else:
		dvec[i] = 0.5 - ( 0.25 * ( fx[pedobj.pedigree[i].damID] - 1.0 ) )
        if pedobj.kw['debug_messages']:
	    print('\t\t[DEBUG]: dvec[%s] = %s' % ( i, dvec[i] ))
        anc.append(i)
        if pedobj.kw['debug_messages']:
	    print('\t\t[DEBUG]: Starting loop over ancestor list for animal %s (idx: %s)' % ( pedobj.pedigree[i].animalID, pedobj.pedigree[i].animalID-1 ))
	while len(anc) > 0:
            if pedobj.kw['debug_messages']:
	        print('\t\t\t[DEBUG]: Starting ancestor list: %s' % ( anc ))
	    jidx = max(anc)				# This is the index of the ID in the pedigree
	    j = pedobj.pedigree[jidx].animalID		# This is an animal ID
            if pedobj.kw['debug_messages']:
	        print('\t\t\t[DEBUG]: Processing Ancestor %s (idx: %s)\tSire: %s\tDam: %s' % \
		( j, jidx, pedobj.pedigree[jidx].sireID, pedobj.pedigree[jidx].damID ))
	    if pedobj.pedigree[jidx].sireID != pedobj.kw['missing_parent'] and pedobj.pedigree[jidx].sireID-1 not in anc:
	        anc.append(pedobj.pedigree[jidx].sireID-1)
            lvec[pedobj.pedigree[jidx].sireID-1] = lvec[pedobj.pedigree[jidx].sireID-1] + ( 0.5 * lvec[jidx] )
	    if pedobj.pedigree[jidx].damID != pedobj.kw['missing_parent'] and pedobj.pedigree[jidx].damID-1 not in anc:
                anc.append(pedobj.pedigree[jidx].damID-1)
            lvec[pedobj.pedigree[jidx].damID-1] = lvec[pedobj.pedigree[jidx].damID-1] + ( 0.5 * lvec[jidx] )
	    if pedobj.kw['debug_messages']:
	        print('\t\t\t\t[DEBUG]: anc      = %s' % ( anc ))
                print('\t\t\t\t[DEBUG]: l[%s]    = %s' % ( jidx, lvec[jidx] ))
	        print('\t\t\t\t[DEBUG]: d[%s]    = %s' % ( jidx, dvec[jidx] ))
	        print('\t\t\t\t[DEBUG]: RHS      = %s' % ( lvec[jidx] * lvec[jidx] * dvec[jidx] ))
	    avec[i] = avec[i] + lvec[jidx] * lvec[jidx] * dvec[jidx]
            if pedobj.kw['debug_messages']:
	        print('\t\t\t\t[DEBUG]: avec[%s] = %s' % ( i, avec[i] ))
	    anc.remove(jidx)
        if pedobj.kw['debug_messages']:
            print('\t\t[DEBUG]: avec[%s] - 1 = %s' % ( i, avec[i]-1 ))
	fx[pedobj.pedigree[i].animalID] = avec[i] - 1.0
        if pedobj.kw['debug_messages']:
	    print('\t\t[DEBUG]: fx[%s] = %s' % ( pedobj.pedigree[i].animalID, fx[pedobj.pedigree[i].animalID] ))
	    print('\t[DEBUG]: Current coefficients of inbeeding: ', fx)
    # We need to clean-up so that we don't have things like memory-ammped files laying around.
    del lvec; del avec; del dvec
    try: logging.info('Exited inbreeding_meuwissen_luo()')
    except: pass
    return fx

##
# inbreeding_modified_meuwissen_luo() computes CoI using the method of Meuwissen and
# Luo (1992). as modified by Quaas (1995). It calculates only inbreeding coefficients,
# not relationships.
# @param pedobj A PyPedal pedigree object.
# @param gens The number of generations from the pedigree to be used for calculating CoI.  By default, gens=0, which uses the complete pedigree.
# @param rels Flag indicating whether or not summary statistics should be computed for coefficients of relationship.
# @retval A dictionary of CoI keyed to renumbered animal IDs
#@profile
def inbreeding_modified_meuwissen_luo(pedobj, gens=0,**kw):
    """
    inbreeding_modified_meuwissen_luo() computes CoI using the method of Meuwissen
    and Luo (1992) as modified by Quaas (1995). It calculates only inbreeding coefficients,
    not relationships. This code is a pretty direct implementation of the algorithm presented
    in Appendix B.2 of Mrode (2005). Mrode cites Quaas's method as: Quaas, R. L. 1995. Fx
    algorithms. An unpublished note.
    """
    try: logging.info('Entered inbreeding_modified_meuwissen_luo()')
    except: pass

    # Setup dictionary to accumulate coefficients of inbreeding
    fx = {}
    for p in pedobj.pedigree:
        fx[p.animalID] = 0.0

    # First, try and allocate the vectors in RAM. If that does not work, try and allocate them using
    # memory-mapped files. If that does not work, well, give up.
    try:
        logging.info('Allocating vectors in pyp_nrm.inbreeding_modified_meuwissen_luo().')
        lvecs = numpy.zeros((len(pedobj.pedigree)),'d')
        #print('lvecs:\t', lvecs)
        lvecd = numpy.zeros((len(pedobj.pedigree)),'d')
        avec = numpy.zeros((len(pedobj.pedigree)),'d')
        dvec = numpy.zeros((len(pedobj.pedigree)),'d')
    except MemoryError:
	logging.info('Unable to allocate a matrix of rank %s in RAM, trying to allocate a memory-mapped file, in pyp_nrm.inbreeding_modified_meuwissen_luo()!', l)
        lvecs = numpy.memmap('lvecs_memmap.bin', dtype='float32', mode='w+', shape=(len(pedobj.pedigree)))
	lvecs = 0.0
        lvecd = numpy.memmap('lvecd_memmap.bin', dtype='float32', mode='w+', shape=(len(pedobj.pedigree)))
	lvecd = 0.0
	avec = numpy.memmap('avec_memmap.bin', dtype='float32', mode='w+', shape=(len(pedobj.pedigree)))
	avec = 0.0
	dvec = numpy.memmap('dvec_memmap.bin', dtype='float32', mode='w+', shape=(len(pedobj.pedigree)))
	dvec = 0.0
    except:
        print('[ERROR]: Unable to allocate a matrix of rank %s in pyp_nrm.inbreeding_modified_meuwissen_luo()!' % ( l ))
        logging.error('[ERROR]: Unable to allocate a matrix of rank %s in pyp_nrm.inbreeding_modified_meuwissen_luo()!', l)
	return False

    if pedobj.kw['debug_messages']: print('[DEBUG]: Starting loop over pedigree with ', len(pedobj.pedigree), ' animals')
    for i in xrange(len(pedobj.pedigree)):
        if pedobj.kw['debug_messages']: print('\t[DEBUG]: Initializing local data structures for animal %s (idx: %s)' % \
	    ( pedobj.pedigree[i].animalID, i))
        ancs = []
	ancd = []
        lvecs[:] = 0.0
	lvecd[:] = 0

        # We're using 0-indexing, so the little F0 = -1 trick that M&L use does not help us here.
        # Check for the most common case first -- both parents unknown
        if pedobj.pedigree[i].sireID != pedobj.kw['missing_parent'] and pedobj.pedigree[i].damID != pedobj.kw['missing_parent']:
            dvec[i] = 0.5 - ( 0.25 * ( fx[pedobj.pedigree[i].sireID] + fx[pedobj.pedigree[i].damID] ) )
        # Then check for both parents unknown
        elif pedobj.pedigree[i].sireID == pedobj.kw['missing_parent'] and pedobj.pedigree[i].damID == pedobj.kw['missing_parent']:
	    if pedobj.kw['debug_messages']: print('\t\t[DEBUG]: Animal %s (idx = %s) has unknown parents' % ( pedobj.pedigree[i].animalID, i ))
            dvec[i] = 1.0
        # Finally, deal with either a known sire, unknown dam or known dam, unknown sire.
        else:
            # Knwon sire, unknown dam
            if pedobj.pedigree[i].sireID != pedobj.kw['missing_parent'] and pedobj.pedigree[i].damID == pedobj.kw['missing_parent']:
                dvec[i] = 0.5 - ( 0.25 * ( fx[pedobj.pedigree[i].sireID] - 1.0 ) )
	    # Knowm dam, unknown sire
            else:
                dvec[i] = 0.5 - ( 0.25 * ( fx[pedobj.pedigree[i].damID] - 1.0 ) )
        if pedobj.kw['debug_messages']: print('\t\t[DEBUG]: dvec[%s] = %s' % ( i, dvec[i] ))
        
        if pedobj.pedigree[i].sireID != pedobj.kw['missing_parent'] and pedobj.pedigree[i].sireID-1 not in ancs:
	    ancs.append(pedobj.pedigree[i].sireID-1)
	    if pedobj.kw['debug_messages']: print('\t\t[DEBUG]: Adding animal %s (idx: %s) to ancs' % ( pedobj.pedigree[i].sireID, pedobj.pedigree[i].sireID-1 ))
	    lvecs[pedobj.pedigree[i].sireID-1] = 1.0
	    if pedobj.kw['debug_messages']: print('\t\t[DEBUG]: lvecs[%s]: %s' % ( pedobj.pedigree[i].sireID-1, lvecs[pedobj.pedigree[i].sireID-1] ))

	if pedobj.pedigree[i].damID != pedobj.kw['missing_parent'] and pedobj.pedigree[i].damID-1 not in ancd:
            ancd.append(pedobj.pedigree[i].damID-1)
	    if pedobj.kw['debug_messages']: print('\t\t[DEBUG]: Adding animal %s (idx: %s) to ancd' % ( pedobj.pedigree[i].damID, pedobj.pedigree[i].damID-1 ))
            lvecd[pedobj.pedigree[i].damID-1] = 1.0
	    if pedobj.kw['debug_messages']: print('\t\t[DEBUG]: lvecd[%s]: %s' % ( pedobj.pedigree[i].damID-1, lvecd[pedobj.pedigree[i].damID-1] ))

	# This loop was miserable to code due in large part to the publisher's decision to use two-point italic typefaces
	# for setting subscripts. Thanks, CABI, that was awesome. It would have been much easier to read the text if 1) it
        # has been larger, and 2) it has been typeset as an algorithm using proper indentation and notation. Lesson learned:
	# use extreme magnification.
        while len(ancs) > 0 and len(ancd) > 0:
            j = max(ancs)
            k = max(ancd)
            if pedobj.kw['debug_messages']: print('\t\t[DEBUG]: j = %s (idx: %s)\tk = %s (idx: %s)' % ( pedobj.pedigree[j].animalID, \
		j, pedobj.pedigree[k].animalID, k ))
	    if j > k:
	        if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: Animal: %s\tSire: %s\tDam: %s' % ( pedobj.pedigree[j].animalID, \
		    pedobj.pedigree[j].sireID, pedobj.pedigree[j].damID ))
		if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: Starting ANCS: %s' % ( ancs ))
                if pedobj.pedigree[j].sireID != pedobj.kw['missing_parent'] and pedobj.pedigree[j].sireID-1 not in ancs:
                    ancs.append(pedobj.pedigree[j].sireID-1)
		    if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: Adding animal %s (idx: %s) to ancs' % ( pedobj.pedigree[j].sireID, \
			pedobj.pedigree[j].sireID-1 ))
                lvecs[pedobj.pedigree[j].sireID-1] = lvecs[pedobj.pedigree[j].sireID-1] + ( 0.5 * lvecs[j] )
		if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: lvecs[%s]: %s' % ( pedobj.pedigree[j].sireID-1, lvecs[pedobj.pedigree[j].sireID-1] ))
		if pedobj.pedigree[j].damID != pedobj.kw['missing_parent'] and pedobj.pedigree[j].damID-1 not in ancs:
                    ancs.append(pedobj.pedigree[j].damID-1)
		    if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: Adding animal %s (idx: %s) to ancs' % ( pedobj.pedigree[j].damID, \
			pedobj.pedigree[j].damID-1 ))
                lvecs[pedobj.pedigree[j].damID-1] = lvecs[pedobj.pedigree[j].damID-1] + ( 0.5 * lvecs[j] )
		if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: lvecs[%s]: %s' % ( pedobj.pedigree[j].damID-1, lvecs[pedobj.pedigree[j].damID-1] ))
	        ancs.remove(j)
		if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: Ending ANCS: %s' % ( ancs ))

            elif k > j:
		if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: Animal: %s\tSire: %s\tDam: %s' % ( pedobj.pedigree[k].animalID, \
		    pedobj.pedigree[k].sireID, pedobj.pedigree[k].damID ))
		if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: Starting ANCD: %s' % ( ancd ))
                if pedobj.pedigree[k].sireID != pedobj.kw['missing_parent'] and pedobj.pedigree[k].sireID-1 not in ancd:
                    ancd.append(pedobj.pedigree[k].sireID-1)
		    if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: Adding animal %s (idx: %s) to ancd' % ( pedobj.pedigree[k].sireID, \
			pedobj.pedigree[k].sireID-1 ))
                lvecd[pedobj.pedigree[k].sireID-1] = lvecd[pedobj.pedigree[k].sireID-1] + ( 0.5 * lvecd[k] )
		if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: lvecd[%s]: %s' % ( pedobj.pedigree[k].sireID-1, lvecd[pedobj.pedigree[k].sireID-1] ))
                if pedobj.pedigree[k].damID != pedobj.kw['missing_parent'] and pedobj.pedigree[k].damID-1 not in ancd:
                    ancd.append(pedobj.pedigree[k].damID-1)
		    if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: Adding animal %s (idx: %s) to ancd' % ( pedobj.pedigree[k].damID, \
			pedobj.pedigree[k].damID-1 ))
                lvecd[pedobj.pedigree[k].damID-1] = lvecd[pedobj.pedigree[k].damID-1] + ( 0.5 * lvecd[k] )
		if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: lvecd[%s]: %s' % ( pedobj.pedigree[k].damID-1, lvecd[pedobj.pedigree[k].damID-1] ))
                ancd.remove(k)
		if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: Ending ANCD: %s' % ( ancd ))

	    else:
		if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: The next youngest ancestor, %s (idx: %s), is a common ancestor' % \
		    ( pedobj.pedigree[j].animalID, j ))
		if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: Animal: %s\tSire: %s\tDam: %s' % ( pedobj.pedigree[k].animalID, \
		    pedobj.pedigree[k].sireID, pedobj.pedigree[k].damID ))
		if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: Ending ANCS: %s' % ( ancs ))
		if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: Ending ANCD: %s' % ( ancd ))
		# If the sire is known
		if pedobj.pedigree[j].sireID != pedobj.kw['missing_parent']:
		    if pedobj.pedigree[j].sireID-1 not in ancs:
                        ancs.append(pedobj.pedigree[j].sireID-1)
		        if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: Adding animal %s (idx: %s) to ancs' % \
			    ( pedobj.pedigree[j].sireID, pedobj.pedigree[j].sireID-1 ))
                    lvecs[pedobj.pedigree[j].sireID-1] = lvecs[pedobj.pedigree[j].sireID-1] + ( 0.5 * lvecs[j] )
		    if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: lvec[%s]: %s' % ( pedobj.pedigree[i].sireID-1, \
			lvecs[pedobj.pedigree[i].sireID-1] ))
		    if pedobj.pedigree[j].sireID-1 not in ancd:
  		        ancd.append(pedobj.pedigree[j].sireID-1)
		        if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: Adding animal %s (idx: %s) to ancd' % \
			    ( pedobj.pedigree[j].sireID, pedobj.pedigree[j].sireID-1 ))
                    lvecd[pedobj.pedigree[j].sireID-1] = lvecd[pedobj.pedigree[j].sireID-1] + ( 0.5 * lvecd[j] )
		    if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: lvec[%s]: %s' % ( pedobj.pedigree[j].sireID-1, \
			lvecd[pedobj.pedigree[j].sireID-1] ))
		# If the dam is known
		if pedobj.pedigree[j].damID != pedobj.kw['missing_parent']:
		    if pedobj.pedigree[j].damID-1 not in ancs:
                        ancs.append(pedobj.pedigree[j].damID-1)
		        if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: Adding animal %s (idx: %s) to ancs' % \
			    ( pedobj.pedigree[j].damID, pedobj.pedigree[j].damID-1 ))
                    lvecs[pedobj.pedigree[j].damID-1] = lvecs[pedobj.pedigree[j].damID-1] + ( 0.5 * lvecs[j] )
                    if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: lvec[%s]: %s' % ( pedobj.pedigree[j].damID-1, \
			lvecs[pedobj.pedigree[j].damID-1] ))
		    if pedobj.pedigree[j].damID-1 not in ancd:
                        ancd.append(pedobj.pedigree[j].damID-1)
		        if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: Adding animal %s (idx: %s) to ancd' % \
			    ( pedobj.pedigree[j].damID, pedobj.pedigree[j].damID-1 ))
                    lvecd[pedobj.pedigree[j].damID-1] = lvecd[pedobj.pedigree[j].damID-1] + ( 0.5 * lvecd[j] )
		    if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: lvec[%s]: %s' % ( pedobj.pedigree[j].damID-1, \
			lvecd[pedobj.pedigree[j].damID-1] ))
 	        fx[pedobj.pedigree[i].animalID] = fx[pedobj.pedigree[i].animalID] + ( lvecs[j] * lvecd[j] * 0.5 * dvec[j] )
		if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: fx[%s] = %s * %s * 0.5 * %s = %s' % ( pedobj.pedigree[i].animalID, \
		    lvecs[j], lvecd[j], dvec[j], fx[pedobj.pedigree[i].animalID] ))
		if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: dvec   = %s' % ( dvec ))
		if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: lvecs  = %s' % ( lvecs ))
		if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: lvecd  = %s' % ( lvecd ))
	        ancs.remove(j)
		if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: Ending ANCS: %s' % ( ancs ))
 	        ancd.remove(j)
		if pedobj.kw['debug_messages']: print('\t\t\t[DEBUG]: Ending ANCD: %s' % ( ancd ))

    # We need to clean-up so that we don't have things like memory-ammped files laying around.
    del lvecs; del lvecd; del avec; del dvec
    try: logging.info('Exited inbreeding_modified_meuwissen_luo()')
    except: pass
    return fx

##
# Form the decomposed form of A, TDT', directly from a pedigree (after
# Henderson, 1976; Thompson, 1977; Mrode, 1996).  Return D, a diagonal
# matrix, and T, a lower triagular matrix such that A = TDT'.
# @param pedobj A PyPedal pedigree object.
# @retval A diagonal matrix, D, and a lower triangular matrix, T.
def a_decompose(pedobj):
    """
    Form the decomposed form of A, TDT', directly from a pedigree (after
    Henderson, 1976; Thompson, 1977; Mrode, 1996).  Return D, a diagonal
    matrix, and T, a lower triagular matrix such that A = TDT'.
    """
    try: logging.info('Entered a_decompose()')
    except: pass
    l = pedobj.metadata.num_records

    if not ( pedobj.kw['form_nrm'] and pedobj.nrm.nrm.shape[0] == pedobj.metadata.num_records ):
        if pedobj.kw['nrm_method'] == 'nrm':
            a = fast_a_matrix(pedobj.pedigree, pedobj.kw, method=pedobj.kw['matrix_type'])
        else:
            a = fast_a_matrix_r(pedobj.pedigree, pedobj.kw, method=pedobj.kw['matrix_type'])
    else:
        a = pedobj.nrm

    try:
        T = numpy.identity(l, dtype=numpy.float)
        D = numpy.identity(l, dtype=numpy.float)
        for row in xrange(l):
            for col in xrange(row+1):
                # cast these b/c items are read from the pedigree file as characters, not  integers
                pedobj.pedigree[col].animalID = int(pedobj.pedigree[col].animalID)
                pedobj.pedigree[col].sireID = int(pedobj.pedigree[col].sireID)
                pedobj.pedigree[col].damID = int(pedobj.pedigree[col].damID)
                if pedobj.pedigree[row].sireID == pedobj.kw['missing_parent'] and pedobj.pedigree[row].damID == pedobj.kw['missing_parent']:
                    if row == col:
                        # both parents unknown and assumed unrelated
                        T[row,col] = 1.
                        D[row,col] = 1.
                    else:
                        T[row,col] = 0.
                elif pedobj.pedigree[row].sireID == pedobj.kw['missing_parent']:
                    # sire unknown, dam known
                    if row == col:
                        T[row,col] = 1.
                        fd = a[pedobj.pedigree[row].damID-1,pedobj.pedigree[row].damID-1] - 1.
                        D[row,col] = 0.75 - ( 0.5 * fd )
                    else:
                        T[row,col] = 0.5 * T[pedobj.pedigree[row].damID-1,col]
                elif pedobj.pedigree[row].damID == pedobj.kw['missing_parent']:
                    # sire known, dam unknown
                    if row == col:
                        T[row,col] = 1.
                        fs = a[pedobj.pedigree[row].sireID-1,pedobj.pedigree[row].sireID-1] - 1.
                        D[row,col] = 0.75 - ( 0.5 * fs )
                    else:
                        T[row,col] = 0.5 * T[pedobj.pedigree[row].sireID-1,col]
                elif pedobj.pedigree[row].sireID != pedobj.kw['missing_parent'] and pedobj.pedigree[row].damID != pedobj.kw['missing_parent']:
                    # both parents known
                    if row == col:
                        T[row,col] = 1.
                        fs = a[pedobj.pedigree[row].sireID-1,pedobj.pedigree[row].sireID-1] - 1.
                        fd = a[pedobj.pedigree[row].damID-1,pedobj.pedigree[row].damID-1] - 1.
                        D[row,col] = 0.5 - ( 0.25 * ( fs + fd ) )
                    else:
                        T[row,col] = 0.5 * ( T[int(pedobj.pedigree[row].sireID)-1,col] + T[int(pedobj.pedigree[row].damID)-1,col] )
                else:
                    print('[ERROR]: There is a problem with the sire (ID %s) and/or dam (ID %s) of animal %s' % (pedobj.pedigree[col].sireID,pedobj.pedigree[col].damID,pedobj.pedigree[col].animalID))
                    break
    except:
        D = numpy.identity(1, dtype=numpy.float)
        T = numpy.identity(1, dtype=numpy.float)

    outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_a_decompose_d_','.dat')
    aout = open(outputfile,'w')
    for row in xrange(l):
        line = ''
        for col in xrange(l):
            if col == 0:
                line = '%7.5f' % (D[row,col])
            else:
                line = '%s%s%s' % (line,',',D[row,col])
        line = '%s%s' % (line,'\n')
        aout.write(line)
    aout.close()

    outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_a_decompose_t_','.dat')
    aout = open(outputfile,'w')
    for row in xrange(l):
        line = ''
        for col in xrange(l):
            if col == 0:
                line = '%7.5f' % (T[row,col])
            else:
                line = '%s%s%s' % (line,',',T[row,col])
        line = '%s%s' % (line,'\n')
        aout.write(line)
    aout.close()

    try: logging.info('Exited a_decompose()')
    except: pass
    return D,T

##
# Form the diagonal matrix, D, used in decomposing A and forming the direct
# inverse of A.  This function does not write output to a file - if you need D in
# a file, use the a_decompose()  function.  form_d() is a convenience function
# used by other functions.  Note that inbreeding is not considered in the
# formation of D.
# @param pedobj A PyPedal pedigree object.
# @retval A diagonal matrix, D.
def form_d_nof(pedobj):
    """
    Form the diagonal matrix, D, used in decomposing A and forming the direct
    inverse of A.  This function does not write output to a file - if you need D in
    a file, use the a_decompose()  function.  form_d() is a convenience function
    used by other functions.  Note that inbreeding is not considered in the
    formation of D.
    """
    try: logging.info('Entered form_d_nof()')
    except: pass
    try:
        l = pedobj.metadata.num_records
        D = numpy.identity(l, dtype=numpy.float)
        for row in xrange(l):
            for col in xrange(row+1):
                # cast these b/c items are read from the pedigree file as characters, not integers
                pedobj.pedigree[col].animalID = int(pedobj.pedigree[col].animalID)
                pedobj.pedigree[col].sireID = int(pedobj.pedigree[col].sireID)
                pedobj.pedigree[col].damID = int(pedobj.pedigree[col].damID)
                if pedobj.pedigree[row].sireID == pedobj.kw['missing_parent'] and pedobj.pedigree[row].damID == pedobj.kw['missing_parent']:
                    if row == col:
                        # both parents unknown and assumed unrelated
                        D[row,col] = 1.
                    else:
                        pass
                elif pedobj.pedigree[row].sireID == pedobj.kw['missing_parent']:
                    # sire unknown, dam known
                    if row == col:
                        D[row,col] = 0.75
                    else:
                        pass
                elif pedobj.pedigree[row].damID == pedobj.kw['missing_parent']:
                    # sire known, dam unknown
                    if row == col:
                        D[row,col] = 0.75
                    else:
                        pass
                elif pedobj.pedigree[row].sireID != pedobj.kw['missing_parent'] and pedobj.pedigree[row].damID != pedobj.kw['missing_parent']:
                    # both parents known
                    if row == col:
                        D[row,col] = 0.5
                    else:
                        pass
                else:
                    print('[ERROR]: There is a problem with the sire (ID %s) and/or dam (ID %s) of animal %s' % (pedobj.pedigree[col].sireID,pedobj.pedigree[col].damID,pedobj.pedigree[col].animalID))
                    break
    except:
        D = numpy.identity(1,dtype=numpy.float)
    try: logging.info('Exited form_d_nof()')
    except: pass
    return D

##
# Form the inverse of A directly using the method of Henderson (1976) which
# does not account for inbreeding.
# @param pedobj A PyPedal pedigree object.
# @param filetag Prefix added to output file names.
# @retval The inverse of the NRM, A, not accounting for inbreeding.
def a_inverse_dnf(pedobj,filetag='_a_inverse_dnf_'):
    """
    Form the inverse of A directly using the method of Henderson (1976) which
    does not account for inbreeding.
    """
    try: logging.info('Entered a_inverse_dnf()')
    except: pass
    l = pedobj.metadata.num_records
    try:
        # grab the diagonal matrix, d, and form its inverse
        d_inv = form_d_nof(pedobj)
        for i in xrange(l):
            d_inv[i,i] = 1. / d_inv[i,i]
        a_inv = numpy.zeros([l,l], dtype=numpy.float)
        for i in xrange(l):
            # cast these b/c items are read from the pedigree file as characters, not integers
            pedobj.pedigree[i].animalID = int(pedobj.pedigree[i].animalID)
            pedobj.pedigree[i].sireID = int(pedobj.pedigree[i].sireID)
            pedobj.pedigree[i].damID = int(pedobj.pedigree[i].damID)
            s = pedobj.pedigree[i].sireID-1
            d = pedobj.pedigree[i].damID-1
            if pedobj.pedigree[i].sireID == pedobj.kw['missing_parent'] and pedobj.pedigree[i].damID == pedobj.kw['missing_parent']:
                # both parents unknown and assumed unrelated
                a_inv[i,i] = a_inv[i,i] + d_inv[i,i]
            elif pedobj.pedigree[i].sireID == pedobj.kw['missing_parent']:
                # sire unknown, dam known
                a_inv[i,i] = a_inv[i,i] + d_inv[i,i]
                a_inv[d,i] = a_inv[d,i] + ( (-0.5) * d_inv[i,i] )
                a_inv[i,d] = a_inv[i,d] + ( (-0.5) * d_inv[i,i] )
                a_inv[d,d] = a_inv[d,d] + ( 0.25 * d_inv[i,i] )
            elif pedobj.pedigree[i].damID == pedobj.kw['missing_parent']:
                # sire known, dam unknown
                a_inv[i,i] = a_inv[i,i] + d_inv[i,i]
                a_inv[s,i] = a_inv[s,i] + ( (-0.5) * d_inv[i,i] )
                a_inv[i,s] = a_inv[i,s] + ( (-0.5) * d_inv[i,i] )
                a_inv[s,s] = a_inv[s,s] + ( 0.25 * d_inv[i,i] )
            elif pedobj.pedigree[i].sireID != pedobj.kw['missing_parent'] and pedobj.pedigree[i].damID != pedobj.kw['missing_parent']:
                # both parents known
                a_inv[i,i] = a_inv[i,i] + d_inv[i,i]
                a_inv[s,i] = a_inv[s,i] + ( (-0.5) * d_inv[i,i] )
                a_inv[i,s] = a_inv[i,s] + ( (-0.5) * d_inv[i,i] )
                a_inv[d,i] = a_inv[d,i] + ( (-0.5) * d_inv[i,i] )
                a_inv[i,d] = a_inv[i,d] + ( (-0.5) * d_inv[i,i] )
                a_inv[s,s] = a_inv[s,s] + ( 0.25 * d_inv[i,i] )
                a_inv[s,d] = a_inv[s,d] + ( 0.25 * d_inv[i,i] )
                a_inv[d,s] = a_inv[d,s] + ( 0.25 * d_inv[i,i] )
                a_inv[d,d] = a_inv[d,d] + ( 0.25 * d_inv[i,i] )
            else:
                print('[ERROR]: There is a problem with the sire (ID %s) and/or dam (ID %s) of animal %s' % (pedobj.pedigree[col].sireID,pedobj.pedigree[col].damID,pedobj.pedigree[col].animalID))
                break
    except:
        a_inv = numpy.zeros([1,1],dtype=numpy.float)

    outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_a_inverse_dnf_a_inv','.dat')
    aout = open(outputfile,'w')
    for row in xrange(l):
        line = ''
        for col in xrange(l):
            if col == 0:
                line = '%7.5f' % (a_inv[row,col])
            else:
                line = '%s%s%s' % (line,',',a_inv[row,col])
        line = '%s%s' % (line,'\n')
        aout.write(line)
    aout.close()

    outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_a_inverse_dnf_d_inv','.dat')
    aout = open(outputfile,'w')
    for row in xrange(l):
        line = ''
        for col in xrange(l):
            if col == 0:
                line = '%7.5f' % (d_inv[row,col])
            else:
                line = '%s%s%s' % (line,',',d_inv[row,col])
        line = '%s%s' % (line,'\n')
        aout.write(line)
    aout.close()

    try: logging.info('Exited a_inverse_dnf()')
    except: pass
    return a_inv

##
# Directly form the inverse of A from the pedigree file - accounts for
# inbreeding - using the method of Quaas (1976).
# @param pedobj A PyPedal pedigree object.
# @retval The inverse of the NRM, A, accounting for inbreeding.
def a_inverse_df(pedobj):
    """
    Directly form the inverse of A from the pedigree file - accounts for
    inbreeding - using the method of Quaas (1976).
    """
    try: logging.info('Entered a_inverse_df()')
    except: pass
    l = pedobj.metadata.num_records
    try:
        from math import sqrt
        # Grab some array tools
        d_inv = numpy.zeros([l,l], dtype=numpy.float)
        a_inv = numpy.zeros([l,l], dtype=numpy.float)
        LL = numpy.zeros([l,l], dtype=numpy.float)
        # Form L and D-inverse
        for row in xrange(l):
            for col in xrange(row+1):
                # cast these b/c items are read from the pedigree file as characters, not integers
                pedobj.pedigree[col].animalID = int(pedobj.pedigree[col].animalID)
                pedobj.pedigree[col].sireID = int(pedobj.pedigree[col].sireID)
                pedobj.pedigree[col].damID = int(pedobj.pedigree[col].damID)
                s = pedobj.pedigree[row].sireID-1
                d = pedobj.pedigree[row].damID-1
                s_sq = d_sq = 0.
                if row == col:
                    for m in xrange(s+1):
                        s_sq = s_sq + ( LL[s,m] * LL[s,m] )
                    s_sq = 0.25 * s_sq
                    for m in xrange(d+1):
                        d_sq = d_sq + ( LL[d,m] * LL[d,m] )
                    d_sq = 0.25 * d_sq
                    LL[row,col] = sqrt(1. - s_sq - d_sq)
                    d_inv[row,col] = 1. / ( LL[row,col] * LL[row,col] )
                else:
                    LL[row,col] = 0.5 * ( LL[s,col] + LL[d,col] )
        # use D-inverse to compute A-inverse
        for i in xrange(l):
            s = pedobj.pedigree[i].sireID-1
            d = pedobj.pedigree[i].damID-1
            if pedobj.pedigree[i].sireID == pedobj.kw['missing_parent'] and pedobj.pedigree[i].damID == pedobj.kw['missing_parent']:
                # both parents unknown and assumed unrelated
                a_inv[i,i] = a_inv[i,i] + d_inv[i,i]
            elif pedobj.pedigree[i].sireID == pedobj.kw['missing_parent']:
                # sire unknown, dam known
                a_inv[i,i] = a_inv[i,i] + d_inv[i,i]
                a_inv[d,i] = a_inv[d,i] + ( (-0.5) * d_inv[i,i] )
                a_inv[i,d] = a_inv[i,d] + ( (-0.5) * d_inv[i,i] )
                a_inv[d,d] = a_inv[d,d] + ( 0.25 * d_inv[i,i] )
            elif pedobj.pedigree[i].damID == pedobj.kw['missing_parent']:
                # sire known, dam unknown
                a_inv[i,i] = a_inv[i,i] + d_inv[i,i]
                a_inv[s,i] = a_inv[s,i] + ( (-0.5) * d_inv[i,i] )
                a_inv[i,s] = a_inv[i,s] + ( (-0.5) * d_inv[i,i] )
                a_inv[s,s] = a_inv[s,s] + ( 0.25 * d_inv[i,i] )
            elif pedobj.pedigree[i].sireID != pedobj.kw['missing_parent'] and pedobj.pedigree[i].damID != pedobj.kw['missing_parent']:
                # both parents known
                a_inv[i,i] = a_inv[i,i] + d_inv[i,i]
                a_inv[s,i] = a_inv[s,i] + ( (-0.5) * d_inv[i,i] )
                a_inv[i,s] = a_inv[i,s] + ( (-0.5) * d_inv[i,i] )
                a_inv[d,i] = a_inv[d,i] + ( (-0.5) * d_inv[i,i] )
                a_inv[i,d] = a_inv[i,d] + ( (-0.5) * d_inv[i,i] )
                a_inv[s,s] = a_inv[s,s] + ( 0.25 * d_inv[i,i] )
                a_inv[s,d] = a_inv[s,d] + ( 0.25 * d_inv[i,i] )
                a_inv[d,s] = a_inv[d,s] + ( 0.25 * d_inv[i,i] )
                a_inv[d,d] = a_inv[d,d] + ( 0.25 * d_inv[i,i] )
    except:
        a_inv = numpy.zeros([1,1], dtype=numpy.float)

    outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_a_inverse_df_a_inv','.dat')
    aout = open(outputfile,'w')
    for row in xrange(l):
        line = ''
        for col in xrange(l):
            if col == 0:
                line = '%7.5f' % (a_inv[row,col])
            else:
                line = '%s%s%s' % (line,',',a_inv[row,col])
        line = '%s%s' % (line,'\n')
        aout.write(line)
    aout.close()

    outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_a_inverse_df_l','.dat')
    aout = open(outputfile,'w')
    for row in xrange(l):
        line = ''
        for col in xrange(l):
            if col == 0:
                line = '%7.5f' % (LL[row,col])
            else:
                line = '%s%s%s' % (line,',',LL[row,col])
        line = '%s%s' % (line,'\n')
        aout.write(line)
    aout.close()

    outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_a_inverse_df_d_inv','.dat')
    aout = open(outputfile,'w')
    for row in xrange(l):
        line = ''
        for col in xrange(l):
            if col == 0:
                line = '%7.5f' % (d_inv[row,col])
            else:
                line = '%s%s%s' % (line,',',d_inv[row,col])
        line = '%s%s' % (line,'\n')
        aout.write(line)
    aout.close()

    try: logging.info('Exited a_inverse_df()')
    except: pass
    return a_inv

##
# partial_inbreeding() computes coefficients of partial inbreeding,
# which is the probability that an individual, i, is IDB at a locus
# and that the alleles were derived from ancestor j.
# @param pedobj A PyPedal pedigree object.
# @param animals An empty list of renumbered animal IDs to process (do all when list empty).
# @param gens The number of generations from the pedigree to be used for calculating CoI.  By default, gens=0, which uses the complete pedigree.
# @param rels Flag indicating whether or not summary statistics should be computed for coefficients of relationship.
# @param cleanmaps Flag to denote whether or not subpedigree ID maps should be deleted after they are used (0|1).
# @retval A dictionary of partial CoI keyed to renumbered animal IDs.
def partial_inbreeding(pedobj, animals=[], gens=0, rels=1, cleanmaps=1):
    """
    partial_inbreeding() computes coefficients of partial inbreeding,
    which is the probability that an individual, i, is IDB at a locus
    and that the alleles were derived from ancestor j.
    Input: A PyPedal pedigree object.
    Output: A dictionary of partial CoI keyed to renumbered animal IDs.
    """

    try: logging.info('Entered partial_inbreeding()')
    except: pass

    from PyPedal import pyp_network
    ng = pyp_network.ped_to_graph(pedobj)

    _ped = []       # This is a temporary pedigree
    top_ped = []

    if int(gens) > 0:
        top_peddict = pyp_network.find_ancestors_g(ng, len(pedobj.idmap), {}, gens)
        top_peddict[len(pedobj.idmap)] = 1
        top_ped = top_peddict.keys()
        top_r = []
        _anids = []
        for _j in top_ped:
            # Make sure that we only include animals from the user-
            # specified generations.
            if top_peddict[_j] <= gens:
                top_r.append(copy.copy(pedobj.pedigree[int(_j)-1]))
                # Animals in the earliest generation need to have
                # their sire and dam IDs set to unknown.
                if top_peddict[_j] == 1:
                    top_r[-1].sireID = 0
                    top_r[-1].damID = 0
                _anids.append(top_r[-1].animalID)
    else:
        _anids = pedobj.backmap.keys()     # Distinct animal IDs in the pedigree

    fx = {}         # This will hold our coefficients of inbreeding
    _parents = {}   # Stores a list of sire-dam pairs along with the youngest offspring
                    # of that pair.  Used as a lookup table to avoid lots of redundant
                    # calculations for full-sibs.
    _anids.sort()       # sort from oldest to youngest
    _anids.reverse()    # reverse the list to put the youngest animals first
    _counter = 0
    _cum_f_counter = 0
    _vanraden_round = 0
    _cum_pct_proc = 0.
    _related = {}   # Dictionary for looking-up animals with non-zero
                    # relationships.

    # If the user wants summary stats on coefficients of relationship,
    # prepare the dictionary, counters, and accumulators.
    if rels:
        reldict = {}
        reldict['r_count'], reldict['r_nonzero_count'] = 0, 0
        reldict['r_sum'], reldict['r_max'] = 0., 0.
        reldict['r_min'] = 1.

    partial_inbreeding_dict = {}

    for i in _anids:
        if int(gens) == 0:
            _parent_key = '%s_%s' % ( pedobj.pedigree[int(i)-1].sireID,
                pedobj.pedigree[int(i)-1].damID )
        else:
            _parent_key = '%s_%s' % ( top_r[top_peddict[int(i)]].sireID, top_r[top_peddict[int(i)]].damID )
        try:
            _k = fx[i]  # If an exception is thrown, an animal is not in the
                        # dictionary yet.
        except KeyError:
            try:
                _pk = _parents[_parent_key]
                fx[i] = fx[_pk]
            except:
                _f_counter = 0
                _vanraden_round = _vanraden_round + 1
                if _vanraden_round == 1:
                    try: logging.info('Starting round %s of pyp_nrm/partial_inbreeding().', _vanraden_round)
                    except: pass
                _tag = '%s_%s' % (pedobj.kw['filetag'],i)

                if int(gens) > 0:
                    _ped = top_peddict
                else:
                    _ped = pyp_network.find_ancestors(ng, i, [])
                    _ped.append(i)
                if int(gens) > 0:
                    _r = top_r
                else:
                    _r = []     # This list will hold a copy of the objects in _ped
                                # so that we can renumber animal i's pedigree without
                                # changing the data in pedobj.pedigree.
                    _map = {}
                    for j in _ped:
                        _r.append(copy.copy(pedobj.pedigree[int(j)-1]))
                _r = pyp_utils.reorder(_r,_tag)      # Reorder the pedigree
                _s, _map = pyp_utils.renumber(_r,_tag, returnmap=1, \
                    debug=pedobj.kw['debug_messages'], \
                    animaltype=pedobj.kw['animal_type'])
                #print('_map: ', _map)
                # We need to get a new founder list here so that we can loop over them.
                _flist = pyp_utils.founders_from_list(_r,pedobj.kw['missing_parent'])
                #print('founders: ', _flist)
                # _map maps IDs from original IDs to renumbered IDs.
                # _backmap allows renumbered ID => original ID reverse lookups.
                _backmap = {}
                for _mk, _mv in _map.iteritems():
                    _backmap[_mv] = _mk
                _opts = copy.copy(pedobj.kw)
                for _f in _flist:
                    _opts['filetag'] = '%s_%s' % ( _tag, _f )
                    _f = fast_partial_a_matrix(_s, _f, _flist, _opts,method=pedobj.kw['matrix_type'])
                    # _f is a dictionary that contains a dictionary, The
                    # key if the founder (_f) that was passed to
                    # fast_partial_a_matrix and the dictionary keyed to
                    # that founder contains animal -> coefficients of partial
                    # inbreeding between an animal and _f.
                    for k,v in _f.iteritems():
                        try: _hask = partial_inbreeding_dict[k]
                        except KeyError: partial_inbreeding_dict[k] = {}
                        for k2,v2 in _f[k].iteritems():
                            partial_inbreeding_dict[k][k2] = v2
                if cleanmaps:               # Clean up the subpedigree ID maps that we are
                    pyp_utils.delete_id_map(_tag)     # not going to use again.
                _map = {}               # Empty our working dictionary and lists
                _a = []
                _s = []
                _r = []
                _ped = []
                logging.info('%s pct (%s) of all animals have been processed in pyp_nrm/partial_inbreeding().', _cum_pct_proc, _cum_f_counter)
            _counter = _counter + 1

        #print('partial_inbreeding_dict')
        #print(partial_inbreeding_dict)

    try:logging.info('Exited partial_inbreeding()')
    except: pass
    #if rels:
        #return fx
    #else:
    return partial_inbreeding_dict

##
# fast_partial_a_matrix() calculates a partial kinship matrix for a given
# founder in a pedigree, and returns a dictionary of partial inbreeding
# coefficients between that founder and and descendants (non-founders) in
# the pedigree.
# @param pedigree A PyPedal pedigree.
# @param founder Founder of interest.
# @param founderlist List of founders in the pedigree.
# @param pedopts PyPedal options.
# @param method Use dense or sparse matrix storage.
# @param debug Print NRM for debugging
# @retval The NRM as Numarray matrix.
def fast_partial_a_matrix(pedigree, founder, founderlist, pedopts, method='dense', debug=0):
    """
    fast_partial_a_matrix() calculates a partial kinship matrix for a given
    founder in a pedigree, and returns a dictionary of partial inbreeding
    coefficients between that founder and and descendants (non-founders) in
    the pedigree.
    """

    _animals = {}
    _sires = {}
    _dams = {}
    l = len(pedigree)
    if method not in ['dense','sparse']:
        method = 'dense'
    # Use PySparse to provide sparse matrix storage for large
    # relationship matrices.
    if method == 'sparse':
        try:
            from pysparse import spmatrix
            a = spmatrix.ll_mat_sym(l*l)
            a = 0.0
        except ImportError:
            logging.error('Could not import spmatrix from PySparse; using Numpy instead!')
            a = numpy.zeros([l,l],'d')  # initialize a matrix of zeros
    # Otherwise, use Numpy and its dense matrices
    else:
        a = numpy.zeros([l,l],'d')  # initialize a matrix of zeros
    if pedopts['debug_messages'] and pedopts['messages'] != 'quiet':
        print('\t\t[pyp_nrm/fast_partial_a_matrix()] Started forming animal, sire, and dam lists at %s' %  pyp_utils.pyp_nice_time())
    for i in xrange(l):
        try:
            _a = _animals[i]
        except KeyError:
            _animals[i] = int(pedigree[i].animalID)
        try:
            _s = _sires[i]
        except KeyError:
            _sires[i] = int(pedigree[i].sireID)
        try:
            _d = _dams[i]
        except KeyError:
            _dams[i] = int(pedigree[i].damID)
    if pedopts['debug_messages'] and pedopts['messages'] != 'quiet':
        print('\t\t[pyp_nrm/fast_partial_a_matrix()] Finished forming animal, sire, and dam lists at %s' %  pyp_utils.pyp_nice_time())
        print('\t\t[pyp_nrm/fast_partial_a_matrix()] Started computing A at %s' %  pyp_utils.pyp_nice_time())

    partial_f = {}
    if debug:
        print('n_founders: ', len(founderlist))
        print('founder   : ', founder)
    # Step 1: Set the len(founderlist)-square block of a to 0. and the element
    #         a[founder,founder] to 0.5. The first part is already done.
    fidx = founderlist.index(founder)
    a[fidx,fidx] = 0.5

    # Step 2:   Intermediate ancestors (non-founders):
    #           2a: Founder row
    row = fidx
    #print('fidx      : ', fidx)
    for col in xrange(len(founderlist),l):
        #print(row, col)
        if str(_sires[col]) == str(pedopts['missing_parent']) and str(_dams[col] == pedopts['missing_parent']):
            pass
        elif str(_sires[col]) == str(pedopts['missing_parent']):
            if row != col:
                a[row,col] = 0.5 * a[row,_dams[col]-1]
                a[col,row] = a[row,col]
        elif str(_dams[col]) == str(pedopts['missing_parent']):
            if row != col:
                a[row,col] = 0.5 * a[row,_sires[col]-1]
                a[col,row] = a[row,col]
        elif str(_sires[col]) != str(pedopts['missing_parent']) and str(_dams[col]) != str(pedopts['missing_parent']):
            # both parents known
            #print(_sires[col])
            #print(_dams[col])
            if row == col:
                a[row,row] = a[row,fidx] + ( 0.5 * a[_sires[col]-1,_dams[col]-1] )
            else:
                a[row,col] = 0.5 * ( a[row,_sires[col]-1] + a[row,_dams[col]-1] )
                a[col,row] = a[row,col]
        else:
            print('[ERROR]: There is a problem with the sire (ID %s) and/or dam (ID %s) of animal %s' % (pedigree[col].sireID, pedigree[col].damID, pedigree[col].animalID))
            break
    #           2a: Non-founder rows
    for row in xrange(len(founderlist),l):
        for col in xrange(row,l):
            if str(_sires[col]) == str(pedopts['missing_parent']) and str(_dams[col]) == str(pedopts['missing_parent']):
                pass
            elif str(_sires[col]) == str(pedopts['missing_parent']):
                if row != col:
                    a[row,col] = 0.5 * a[row,_dams[col]-1]
                    a[col,row] = a[row,col]
            elif str(_dams[col]) == str(pedopts['missing_parent']):
                if row != col:
                    a[row,col] = 0.5 * a[row,_sires[col]-1]
                    a[col,row] = a[row,col]
            elif str(_sires[col]) != str(pedopts['missing_parent']) and str(_dams[col]) != str(pedopts['missing_parent']):
                # both parents known
                if row == col:
                    a[row,row] = a[row,fidx] + ( 0.5 * a[_sires[col]-1,_dams[col]-1] )
                else:
                    a[row,col] = 0.5 * ( a[row,_sires[col]-1] + a[row,_dams[col]-1] )
                    a[col,row] = a[row,col]
            else:
                print('[ERROR]: There is a problem with the sire (ID %s) and/or dam (ID %s) of animal %s' % (pedigree[col].sireID, pedigree[col].damID, pedigree[col].animalID))
                break

    # Step 3:   Partial inbreeding coefficients
    #print('-'*80)
    partial_f[founder] = {}
    for i in xrange(len(founderlist),l):
        partial_f[founder][_animals[i]] = a[_sires[i]-1,_dams[i]-1]
        #print(_animals[i], _sires[i], _dams[i])

    if debug:
        numpy.set_printoptions(precision=4,linewidth=100)
        print(a)
        print(partial_f)

    if pedopts['debug_messages'] and pedopts['messages'] != 'quiet':
        print('\t\t[pyp_nrm/fast_partial_a_matrix()] Finished computing A at %s' %  pyp_utils.pyp_nice_time())

    return partial_f
