#!/usr/bin/python

###############################################################################
# NAME: pyp_metrics.py
# VERSION: 2.0.4 (28JUNE2022)
# AUTHOR: John B. Cole (john.b.cole@gmail.com)
# LICENSE: LGPL
###############################################################################
# FUNCTIONS:
#   min_max_f()
#   a_effective_founders_lacy()
#   effective_founders_lacy()
#   a_effective_founders_boichard()
#   a_effective_ancestors_definite()
#   a_effective_ancestors_indefinite()
#   a_coefficients()
#   fast_a_coefficients()
#   theoretical_ne_from_metadata()
#   pedigree_completeness()
#   related_animals()
#   common_ancestors()
#   relationship()
#   mating_coi()
#   mating_coi_group()
#   effective_founder_genomes()
#   generation_lengths()
#   generation_lengths_all()
#   founder_descendants()
#   descendants()
#   dropped_ancestral_inbreeding()
#   ballou_ancestral_inbreeding()
###############################################################################

## @package pyp_metrics
# pyp_metrics contains a set of procedures for calculating metrics on PyPedal
# pedigree objects.  These metrics include coefficients of inbreeding and
# relationship as well as effective founder number, effective population size,
# and effective ancestor number.
##

import copy
import logging
import numpy
# import os
# import pickle
import random
# import string
# import sys
import pyp_io
import pyp_network
import pyp_nrm
import pyp_utils
from numpy import random

##
# min_max_f() takes a pedigree and returns a list of the individuals with the n
# largest and n smallest coefficients of inbreeding.  Individuals with CoI of
# zero are not included.
# @param pedobj A PyPedal pedigree object.
# @param a A numerator relationship matrix (optional).
# @param n An integer (optional, default is 10) number of coefficients to return (e.g., 10 smallest/largest).
# @param forma If A must be formed should dense or sparse matrices be used?
# @retval Lists of the n largest and n smallest CoI in the pedigree, or False on failure.
def min_max_f(pedobj, a='', n=10, forma='dense'):
    """
    Given a pedigree or relationship matrix, return a list of the
    individuals with the n largest and n smallest coefficients of
    inbreeding; individuals with CoI of zero are not included.
    :param pedobj: A PyPedal pedigree object.
    :param a: A numerator relationship matrix (optional).
    :param n: An integer (optional, default is 10) number of coefficients to return (e.g., 10 smallest/largest).
    :param forma: If A must be formed should dense or sparse matrices be used?
    :return: Lists of the n largest and the  n smallest CoI in the pedigree, or False on failure.
    """
    logging.info('Entered min_max_f()')
    if forma not in ['dense', 'sparse']:
        forma = 'dense'
    if not pedobj.kw['form_nrm'] and not a:
        a = pyp_nrm.fast_a_matrix(pedobj.pedigree, pedobj.kw, method=forma)
        if pedobj.kw['debug_messages']:
            print 'Value of matrix a in pyp_metrics/min_max_f():'
            print a
        individual_coi = fast_a_coefficients(pedobj, a=a)
        if pedobj.kw['debug_messages']:
            print 'Value of individual_coi in pyp_metrics/min_max_f():'
            print individual_coi
    else:
        individual_coi = fast_a_coefficients(pedobj)

    _mycoi = individual_coi.items()
    _mycoi.sort(pyp_utils.sort_dict_by_values)

    # If the user asks for more high-low COI than there are animals in
    # the pedigree with non-zero COI set the new n to approximately half of the
    # animals returned by fast_a_coefficients().
    if n > len(_mycoi):
        old_n = n
        if len(_mycoi) % 2 == 0:
            n = int(round(float(len(_mycoi))/2., 0))
        else:
            n = int(round(float(len(_mycoi))/2. - 1.0, 0))
        print '[INFO]: You asked for more high and low COI, %s, than there are animals in the pedigree with ' \
              'non-zero COI, %s. n was adjusted from %s to %s.' % (old_n, len(_mycoi), old_n, n)
        logging.info('[INFO]: You asked for more high and low COI, %s, than there are animals in the pedigree '
                     'with non-zero COI, %s. n was adjusted from %s to %s.',old_n, len(_mycoi), old_n, n)

    high_coi = []
    low_coi = []

    for _i in range(n):
        high_coi.append(_mycoi[len(_mycoi)-_i-1])
        low_coi.append(_mycoi[_i])

    if pedobj.kw['debug_messages']:
        logging.info('Exited min_max_f()')

    return high_coi, low_coi


##
# a_effective_founders_lacy() calculates the number of effective founders in a
# pedigree using the exact method of Lacy.
# @param pedobj A PyPedal pedigree object.
# @param a A numerator relationship matrix (optional).
# @param half Use half founders or not, default is to use only full founders (Boolean, optional).
# @retval A dictionary of results, including the effective founder number.
def a_effective_founders_lacy(pedobj, a='', half=False):
    """
    Calculate the number of effective founders in a pedigree using the exact method of Lacy.
    :param pedobj: A PyPedal pedigree object.
    :param a: A numerator relationship matrix (optional).
    :param half: Use half founders or not, default is to use only full founders (Boolean, optional).
    :return: A dictionary of results, including the effective founder number.
    """
    if pedobj.kw['debug_messages']:
        logging.info('Entered a_effective_founders_lacy()')
    if not a:
        try:
            a = pyp_nrm.fast_a_matrix(pedobj.pedigree, pedobj.kw)
        except:
            return -999.9
    lenped = len(pedobj.pedigree)
    # form lists of founders and descendants
    out_dict = {}
    n_f = 0
    n_d = 0
    fs = []
    ds = []
    for i in range(lenped):
        animalid = int(pedobj.pedigree[i].animalID)
        sireid = int(pedobj.pedigree[i].sireID)
        damid = int(pedobj.pedigree[i].damID)
        #
        if pedobj.pedigree[i].founder == 'y' and animalid != pedobj.kw['missing_parent']:
            n_f = n_f + 1
            fs.append(animalid)
        # An animal with a missing parent is a half-founder. Note that this is not a part of Lacy's original
        # derivation, and I haven't proven it works, either. This is a strictly heuristic approach at the
        # moment. An alternative approach might be to create metafounders in place of half founders.
        elif half and pedobj.pedigree[i].founder == 'n' and (sireid == pedobj.kw['missing_parent'] or
                                                             damid == pedobj.kw['missing_parent']):
            if pedobj.kw['debug_messages']:
                print '[pyp_metrics]: Animal % is a half-founder...' % animalid
            # n_f = n_f + 0.5
            # Can't be 0.5 b/c of table below listing relationship between founders and descendants
            n_f += 1
            fs.append(animalid)
        elif animalid != pedobj.kw['missing_parent']:
            n_d = n_d + 1
            ds.append(animalid)
        else:
            pass
    # print 'fs : %s' % fs
    # print 'ds : %s' % ds
    p = numpy.zeros([n_d, n_f], 'd')
    # create a table listing relationship between founders and descendants
    for row in range(n_d):
        for col in range(n_f):
            p[row, col] = a[fs[col]-1, ds[row]-1]
    # print 'p : %s' % p
    # sum each column
    p_sums = []
    for col in range(n_f):
        p_sum = 0.
        for row in range(n_d):
            if p[row, col] != 0:
                p_sum = p_sum + p[row, col]
        p_sums.append(p_sum)
    # weight sums by counts to get relative contributions
    rel_p = []
    rel_p_sq = []
    for i in range(len(p_sums)):
        rel_p.append(p_sums[i] / n_d)
        rel_p_sq.append(rel_p[i] * rel_p[i])
    # sum  the squared relative contributions and take the reciprocal to get f_e
    sum_rel_p_sq = 0.
    for i in range(len(rel_p_sq)):
        sum_rel_p_sq = sum_rel_p_sq + rel_p_sq[i]
    print '='*60
    # print 'p_sums:\t%s' % p_sums
    # print 'rel_ps:\t%s' % rel_p
    if sum_rel_p_sq == 0.:
        f_e = 0.
    else:
        f_e = 1. / sum_rel_p_sq
    if pedobj.kw['messages'] == 'verbose':
        print 'animals:\t%s' % len(fs)+len(ds)
        print 'founders:\t%s' % n_f
        print 'descendants:\t%s' % n_d
        print 'f_e:\t\t%5.3f' % f_e
        print '='*60

    out_dict['fa_animal_count'] = len(fs) + len(ds)
    out_dict['fa_founder_count'] = n_f
    out_dict['fa_descendant_count'] = n_d
    out_dict['fa_effective_founders'] = f_e

    # write some output to a file for later use
    outputfile = '%s%s%s' % (pedobj.kw['filetag'], '_fe_lacy_', '.dat')
    aout = open(outputfile, 'w')
    line = '='*60+'\n'
    aout.write('%s\n' % line)
    aout.write('%s animals\n' % lenped)
    aout.write('%s founders: %s\n' % (n_f, fs))
    aout.write('%s descendants: %s\n' % (n_d, ds))
    aout.write('effective number of founders: %s\n' % f_e)
    aout.write('%s\n' % line)
    aout.close()

    if pedobj.kw['debug_messages']:
        logging.info('Exited a_effective_founders_lacy()')

    return out_dict


##
# effective_founders_lacy() calculates the number of effective founders in a pedigree
# using the exact method of Lacy.  This version of the routine a_effective_founders_lacy()
# is designed to work with larger pedigrees as it forms "family-wise" relationship matrices
# rather than a "population-wise" relationship matrix.
# @param pedobj A PyPedal pedigree object.
# @retval A dictionary of results, including the effective founder number.
def effective_founders_lacy(pedobj):
    """
    Calculate the number of effective founders in a pedigree using the exact method of Lacy.
    :param pedobj: A PyPedal pedigree object.
    :return: A dictionary of results, including the effective founder number.
    """
    if pedobj.kw['debug_messages']:
        logging.info('Entered effective_founders_lacy()')

    caller = 'pyp_metrics.effective_founders_lacy'
    out_dict = {}

    # founder_descendants is expecting a renumbered pedigree.
    # print 'pedobj.kw[\'pedigree_is_renumbered\'] = %s' % ( pedobj.kw['pedigree_is_renumbered'] )
    if not pedobj.kw['pedigree_is_renumbered']:
        if pedobj.kw['debug_messages']:
            print '[NOTE]: The pedigree passed to pyp_metrics/effective_founders_lacy() is not renumbered! Fixing...'
        pedobj.renumber()
        pedobj.kw['renumber'] = True        # Set this last in case renumbering fails to avoid inconsistent state.
    else:
        if pedobj.kw['debug_messages']:
            print '[NOTE]: The pedigree passed to pyp_metrics/effective_founders_lacy() is renumbered!'
        pass

    _f_peds = founder_descendants(pedobj)
    _f_contribs = {}
    _f_contribs_sum = 0.0
    _f_contribs_weighted = {}
    _f_contribs_weighted_sum_sq = 0.0
    _f_e = 0.0
    for f in _f_peds.keys():
        # Note that the "pedigrees" returned by founder_descendants are really
        # just dictionaries of dictionaries.  We need to form real pedigrees from
        # them before we can continue.
        # print 'Working on animal %s' % ( pedobj.idmap[f] )
        # lenped = len(_f_peds[f])
        _r = []
        for _k, _v in _f_peds[f].iteritems():
            _r.append(copy.deepcopy(pedobj.pedigree[int(_k)-1]))
        # Don't forget to add the founder to the pedigree!
        _r.append(copy.deepcopy(pedobj.pedigree[int(pedobj.idmap[f])-1]))
        _tag = '%s_%s' % (pedobj.kw['filetag'], pedobj.idmap[f])
        _r = pyp_utils.fast_reorder(_r, _tag)       # Reorder the pedigree
        _s = pyp_utils.renumber(_r, _tag)           # Renumber the pedigree
        _opts = copy.deepcopy(pedobj.kw)
        _opts['filetag'] = _tag
        _a = pyp_nrm.fast_a_matrix(_s, _opts)       # Form the NRM w/the tabular method
        _map = pyp_utils.load_id_map(_tag)          # Load the ID map from the renumbering
                                                    # procedure so that the CoIs are assigned
                                                    # to the correct animals.
        # We are going to sum over each slice to get the (unweighted) founder
        # contribution associated with founder f.
        _f_contribs[f] = _a[0, 1:].sum()
        _f_contribs_sum = _f_contribs_sum + _f_contribs[f]
        # Clean up the subpedigree ID maps that we are not going to use again.
        pyp_utils.delete_id_map(_tag)
        # Empty our working dictionary and lists
        _a = []
        _s = []
        _r = []
        _map = {}
    for _k,_v in _f_contribs.iteritems():
        try:
            _f_contribs_weighted[_k] = _v / _f_contribs_sum
        except ZeroDivisionError:
            _f_contribs_weighted[_k] = 0.0
        _f_contribs_weighted_sum_sq = _f_contribs_weighted_sum_sq + (_f_contribs_weighted[_k] ** 2)
    try:
        _f_e = 1. / _f_contribs_weighted_sum_sq
    except ZeroDivisionError:
        _f_e = 0.0

    if pedobj.kw['messages'] == 'verbose':
        print '\tFounder contributions: %s' % _f_contribs
        print '\tProportional founder contributions: %s' % _f_contribs_weighted
        print '\tSum of founder contributions: %s' % _f_contribs_sum
        print '\tSum of squared proportional founder contributions: %s' % _f_contribs_weighted_sum_sq
        print 'Effective founder number (f_e): %s' % _f_e

    # write some output to a file for later use
    outputfile = '%s%s%s' % (pedobj.kw['filetag'], '_fe_lacy', '.dat')
    aout = open(outputfile, 'w')
    pyp_io.pyp_file_header(aout, caller)
    aout.write('%s animals in pedigree\n' % (len(pedobj.pedigree)))
    aout.write('%s founders: %s\n' % (pedobj.metadata.num_unique_founders,pedobj.metadata.unique_founder_list))
    aout.write('Founder contributions: %s\n' % _f_contribs)
    aout.write('Proportional founder contributions: %s\n' % _f_contribs_weighted)
    aout.write('Sum of founder contributions: %s\n' % _f_contribs_sum)
    aout.write('Sum of squared proportional founder contributions: %s\n' % _f_contribs_weighted_sum_sq)
    aout.write('Effective founder number: %s\n' % _f_e)
    pyp_io.pyp_file_footer(aout, caller)
    aout.close()

    out_dict['fa_animal_count'] = len(pedobj.pedigree)
    out_dict['fa_founder_count'] = pedobj.metadata.num_unique_founders
    out_dict['fa_descendant_count'] = len(pedobj.pedigree) - pedobj.metadata.num_unique_founders
    out_dict['fa_effective_founders'] = _f_e

    if pedobj.kw['debug_messages']:
        logging.info('Exited effective_founders_lacy()')

    return out_dict


##
# a_effective_founders_boichard() uses the algorithm in Appendix A of Boichard et al.
# (1996) to compute the effective founder number for myped.  Note that results from
# this function will not necessarily match those from a_effective_founders_lacy().
# Boichard's algorithm requires information about the GENERATION of animals.  If you
# do not provide an input pedigree with generations things may not work.  By default
# the most recent generation -- the generation with the largest generation ID -- will
# be used as the reference population.
# @param pedobj A PyPedal pedigree object.
# @param a A numerator relationship matrix (optional).
# @param gen Generation of interest.
# @retval The effective founder number.
def a_effective_founders_boichard(pedobj, a='', gen=''):
    """
    The algorithm in Appendix A of Boichard et al. (1996) is not very well written.
    a_effective_founders_boichard() implements that algorithm (successfully, I hope).
    Note that answers from this function will not necessarily match those from
    a_effective_founders_lacy().
    :param pedobj: A PyPedal pedigree object.
    :param a: A numerator relationship matrix (optional).
    :param gen: Generation of interest.
    :return: A dictionary of results, including the effective founder number.
    """

    if pedobj.kw['debug_messages']:
        logging.info('Entered a_effective_founders_boichard()')

    if not a:
        try:
            a = pyp_nrm.fast_a_matrix(pedobj.pedigree,pedobj.kw)
        except:
            return -999.9
    lenped = len(pedobj.pedigree)
    # count founders and descendants
    n_f = 0
    n_d = 0
    fs = []
    ds = []
    gens = []
    ngen = 0    # number of individuals in the most recent generation
    # loop through the pedigree quickly and count founders and descendants
    # also, make a list of generations present in the pedigree
    for i in range(lenped):
        animalid = int(pedobj.pedigree[i].animalID)
        if pedobj.pedigree[i].founder == 'y' and str(animalid) != str(pedobj.kw['missing_parent']):
            n_f = n_f + 1
            fs.append(animalid)
        elif animalid != str(pedobj.kw['missing_parent']):
            n_d = n_d + 1
            ds.append(animalid)
        else:
            pass
        # g = int(pedobj.pedigree[i].gen)
        g = pedobj.pedigree[i].gen
        if g in gens:
            pass
        else:
            gens.append(g)
    # OK - now we have a list of generations sorted in reverse (descending) order
    gens.sort()
    gens.reverse()
    # print 'gens : %s' % (gens)
    # print 'fs : %s' % (fs)
    # print 'ds : %s' % (ds)
    # make a copy of pedobj.pedigree
    tempped = pedobj.pedigree[:]
    # reverse the elements of tempped in place
    # now animals are ordered from oldest to youngest in tempped
    tempped.reverse()
    # form q, a vector of that will contain the probabilities of gene origin when we are done
    # We are going to initialize a vector of numpy.zeros to form q, and then we will add ones to q that
    # correspond to members of the youngest generation.
    q = numpy.zeros([lenped], 'd')
    for i in range(lenped):
        # If the user did not explicitly ask for an analysis of a particular generation then use
        # the most recent generation.
        if not gen:
            gen = gens[0]
        if pedobj.pedigree[i].gen == gens[0]:
            # be careful messing with this or the elements of q will end up in the wrong
            # columns
            q[i] = 1.
            ngen = ngen + 1
        else:
            pass
    # print 'DEBUG (e_f_b): q : %s' % (q)
    # loop through the pedigree and form the final version of q (the vector of
    # individual contributions)
    for i in range(lenped):
        animalid = int(tempped[i].animalID)
        sireid = int(tempped[i].sireID)
        damid = int(tempped[i].damID)
        if sireid == str(pedobj.kw['missing_parent']) and damid == str(pedobj.kw['missing_parent']):
            # both parents unknown
            pass
        elif str(sireid) == str(pedobj.kw['missing_parent']):
            # sire unknown, dam known
            q[i] = q[animalid-1] * 0.5
            q[damid-1] = q[damid-1] + (0.5 * q[animalid-1])
        elif str(damid) == str(pedobj.kw['missing_parent']):
            # sire known, dam unknown
            q[i] = q[animalid-1] * 0.5
            q[sireid-1] = q[sireid-1] + (0.5 * q[animalid-1])
        else:
            # both parents known
            q[sireid-1] = q[sireid-1] + (0.5 * q[animalid-1])
            q[damid-1] = q[damid-1] + (0.5 * q[animalid-1])
    # print q
    # divide the elements of q by the number of individuals in the pedigree.  this should
    # ensure that the founder contributions sum to 1.
    q = q / ngen
    # print 'DEBUG (e_f_b): q : %s' % (q)
    # accumulate the sum of squared founder contributions
    sum_sq = 0.
    sum_fn = 0.
    # print fs
    for i in fs:
        sum_sq = sum_sq + (q[i-1] * q[i-1])
        sum_fn = sum_fn + q[i-1]
    # print '='*60
    # print 'sum_fn:\t%s' % (sum_fn)
    # print 'sum_sq:\t%s' % (sum_sq)
    if sum_sq == 0.:
        f_e = 0.
    else:
        f_e = 1. / sum_sq
    if pedobj.kw['messages'] == 'verbose':
        print '='*60
        print 'animals:\t%s' % lenped
        print 'founders:\t%s' % n_f
        print 'descendants:\t%s' % n_d
        print 'f_e:\t\t%5.3f' % f_e
        print '='*60

    # write some output to a file for later use
    outputfile = '%s%s%s' % (pedobj.kw['filetag'], '_fe_boichard_', '.dat')
    aout = open(outputfile, 'w')
    line = '='*60+'\n'
    aout.write('%s\n' % line)
    aout.write('q: %s\n' % q)
    aout.write('%s founders: %s\n' % (n_f, fs))
    aout.write('%s descendants: %s\n' % (n_d, ds))
    aout.write('generations: %s\n' % gens)
    aout.write('%s animals in generation %s\n' % (ngen, gens[0]))
    aout.write('effective number of founders: %s\n' % f_e)
    aout.write('%s\n' % line)
    aout.close()

    if pedobj.kw['debug_messages']:
        logging.info('Exited a_effective_founders_boichard()')

    return f_e


##
# a_effective_ancestors_definite() uses the algorithm in Appendix B of Boichard et al. (1996) to compute the effective
# ancestor number for a myped pedigree. NOTE: One problem here is that if you pass a pedigree WITHOUT generations an
# error is not thrown. You simply end up wth a list of generations that contains the default value for Animal() objects,
# 0. Boichard's algorithm requires information about the GENERATION of animals.  If you do not provide an input pedigree
# with generations things may not work. By default the most recent generation -- the generation with the largest
# generation ID -- will be used as the reference population.
# @param pedobj A PyPedal pedigree object.
# @param a A numerator relationship matrix (optional).
# @param gen Generation of interest.
# @retval The effective ancestor number, False on failure.
def a_effective_ancestors_definite(pedobj, a='', gen=''):
    """
    The algorithm in Appendix B of Boichard et al. (1996) is not very well written.
    a_effective_ancestors_definite() implements that algorithm (successfully, I hope).

    NOTE: One problem here is that if you pass a pedigree WITHOUT generations an error
    is not thrown.  You simply end up wth a list of generations that contains the
    default value for Animal() objects, 0.
    :param pedobj: A PyPedal pedigree object.
    :param a: A numerator relationship matrix (optional).
    :param gen: Generation of interest.
    :return: The effective ancestor number, False on failure.
    """
    if pedobj.kw['debug_messages']:
        logging.info('Entered a_effective_ancestors_definite()')

    if not a:
        try:
            a = pyp_nrm.fast_a_matrix(pedobj.pedigree, pedobj.kw)
        except:
            return False
    lenped = len(pedobj.pedigree)  # number of animals in the pedigree file
    # count founders and descendants
    n_f = 0     # number of founders
    n_d = 0     # number of descendants
    fs = []     # list of founders
    ds = []     # list of descendants
    gens = []       # list of generation IDs in the pedigree
    ancestors = []  # list of ancestors already processed
    contribs = {}   # ancestor contributions
    ngen = 0        # number of individuals in the most recent generation
    # Loop through the pedigree quickly and count founders and descendants
    # also, make a list of generations present in the pedigree
    for i in range(lenped):
        g = pedobj.pedigree[i].gen
        if g in gens:
            pass
        else:
            gens.append(g)
    gens.sort()
    # print '[DEBUG]: gens: %s' % gens
    if int(gens[-1]) == -999:
        logging.warning('pyp_metrics/a_effective_ancestors_definite() assumes that generations are defined in the '
                        'pedigree. That is not the case with %s. Solutions from this routine may be inaccurate or '
                        'nonsensical.', pedobj.kw['pedname'])
    for i in range(lenped):
        # print 'DEBUG: Animal: %s\tGen: %s' % (pedobj.pedigree[i].animalID, pedobj.pedigree[i].gen)
        if pedobj.pedigree[i].gen != gens[len(gens)-1]:
            n_f = n_f + 1
            fs.append(int(pedobj.pedigree[i].animalID))
        else:
            n_d = n_d + 1
            ds.append(int(pedobj.pedigree[i].animalID))
    # OK - now we have a list of generations sorted in reverse (descending) order
    ngen = len(gens)
    gens.sort()
    # print 'DEBUG: gens: %s' % gens
    # print 'DEBUG: ancestors: %s' % fs
    gens.reverse()
    # make a copy of pedobj.pedigree - note that tempped = pedobj.pedigree would only have created a reference to
    # pedobj.pedigree, not an actual separate copy of pedobj.pedigree.
    tempped = pedobj.pedigree[:]
    # now animals are ordered from oldest to youngest in tempped
    tempped.reverse()
    # form q, a vector of that will contain the probabilities of gene origin when we are done
    # We are going to initialize a vector of numpy.zeros to form q, and then we will add ones to q that
    # correspond to members of the youngest generation.
    younglist = []
    q = numpy.zeros([lenped], 'd')
    for i in range(lenped):
        # If the user did not explicitly ask for an analysis of a particular generation then use
        # the most recent generation.
        if not gen:
            gen = gens[0]
        # print 'DEBUG: Most recent generation = %s' % gen
        if pedobj.pedigree[i].gen == gen:
            q[i] = 1.
            younglist.append(pedobj.pedigree[i].animalID)
    ngen = len(younglist)
    # print 'DEBUG: Young animals (n=%s) = %s' % (ngen, younglist)
    #
    # Algorithm B, Step 1
    #
    # print 'DEBUG: q : %s' % q
    # loop through the pedigree and form the final version of q (the vector of
    # individual contributions)
    for i in range(lenped):
        animalid = int(tempped[i].animalID)
        sireid = int(tempped[i].sireID)
        damid = int(tempped[i].damID)
        if str(sireid) == str(pedobj.kw['missing_parent']) and str(damid) == str(pedobj.kw['missing_parent']):
            # both parents unknown
            pass
        elif str(sireid) == str(pedobj.kw['missing_parent']):
            # sire unknown, dam known
            q[i] = q[animalid-1] * 0.5
            q[damid-1] = q[damid-1] + (0.5 * q[animalid-1])
        elif str(damid) == str(pedobj.kw['missing_parent']):
            # sire known, dam unknown
            q[i] = q[animalid-1] * 0.5
            q[sireid-1] = q[sireid-1] + (0.5 * q[animalid-1])
        else:
            # both parents known
            q[sireid-1] = q[sireid-1] + (0.5 * q[animalid-1])
            q[damid-1] = q[damid-1] + (0.5 * q[animalid-1])
    # Divide the elements of q by the number of individuals in the pedigree.  this should
    # ensure that the founder contributions sum to 1.
    for y in younglist:
        q[int(y)-1] = 0.
    # print 'DEBUG: Uncorrected q: %s' % q
    q = q / ngen
    # print 'DEBUG: q: %s' % q

    # Find largest value of q
    max_p_index = numpy.argmax(q)
    max_p = q[max_p_index]
    # print 'DEBUG: Animal %s had the largest marginal contribution (%s) (index: %s) this round.' % \
    #       (tempped[l-max_p_index-1].animalID, max_p,max_p_index)
    contribs[pedobj.pedigree[max_p_index].animalID] = max_p
    picked = []
    picked.append(lenped-max_p_index-1)
    # print '\t\tWas sire: %s, dam %s' % (tempped[lenped-max_p_index-1].sireID, tempped[lenped-max_p_index-1].damID)
    tempped[lenped-max_p_index-1].sireID = pedobj.kw['missing_parent']      # delete sire in pedobj.pedigree
                                                                            # (forward order)
    tempped[lenped-max_p_index-1].damID = pedobj.kw['missing_parent']       # delete dam in pedobj.pedigree
                                                                            # (forward order)
    # print '\t\tNow sire: %s, dam %s' % (tempped[lenped-max_p_index-1].sireID,tempped[lenped-max_p_index-1].damID)
    ancestors.append(pedobj.pedigree[max_p_index].animalID)         # add the animal with largest q to the
                                                                    # list of ancestors

    for j in range(n_f-1):
        # form q, the vector of contributions we are going to use
        q = numpy.zeros([lenped], 'd')
        a = numpy.zeros([lenped], 'd')
        for i in range(lenped):
            if pedobj.pedigree[i].gen == gens[0]:
                q[int(pedobj.pedigree[i].animalID)-1] = 1.
        for j in ancestors:
            a[int(j)-1] = 1.
        # print 'DEBUG: a: %s' % a

        # Loop through pedigree to process q
        # q must be processed from YOUNGEST to OLDEST
        for i in range(len(tempped)):
            if str(tempped[i].sireID) == str(pedobj.kw['missing_parent']) and \
                    str(tempped[i].damID) == str(pedobj.kw['missing_parent']):
                # both parents unknown
                pass
            elif str(tempped[i].sireID) == str(pedobj.kw['missing_parent']):
                # sire unknown, dam known
                q[int(tempped[i].damID)-1] = q[int(tempped[i].damID)-1] + (0.5 * q[int(tempped[i].animalID)-1])
            elif str(tempped[i].damID) == str(pedobj.kw['missing_parent']):
                # sire known, dam unknown
                q[int(tempped[i].sireID)-1] = q[int(tempped[i].sireID)-1] + (0.5 * q[int(tempped[i].animalID)-1])
            else:
                # both parents known
                q[int(tempped[i].sireID)-1] = q[int(tempped[i].sireID)-1] + (0.5 * q[int(tempped[i].animalID)-1])
                q[int(tempped[i].damID)-1] = q[int(tempped[i].damID)-1] + (0.5 * q[int(tempped[i].animalID)-1])
        # a must be processed from OLDEST to YOUNGEST
        tempped.reverse()
        for i in range(len(tempped)):
            # print '[DEBUG]: animal: %s, sire: %s, dam %s' % (tempped[i].animalID, tempped[i].sireID, tempped[i].damID)
            if str(tempped[i].sireID) == str(pedobj.kw['missing_parent']) and \
                    str(tempped[i].damID) == str(pedobj.kw['missing_parent']):
                # both parents unknown
                pass
            elif str(tempped[i].sireID) == str(pedobj.kw['missing_parent']):
                # sire unknown, dam known
                a[i] = a[i] + (0.5 * a[int(pedobj.pedigree[i].damID)-1])
            elif str(tempped[i].damID) == str(pedobj.kw['missing_parent']):
                # sire known, dam unknown
                a[i] = a[i] + (0.5 * a[int(tempped[i].sireID)-1])
            else:
                # both parents known
                a[i] = a[i] + (0.5 * a[int(tempped[i].sireID)-1])
                a[i] = a[i] + (0.5 * a[int(tempped[i].damID)-1])
        tempped.reverse()
        for y in younglist:
            q[int(y)-1] = 0.

        # print 'DEBUG: post q: %s' % q
        # print 'DEBUG: post a: %s' % a
        # Loop through the pedigree to process p
        p = numpy.zeros([lenped], 'd')
        for i in range(lenped):
            p[i] = q[i] * (1. - a[i])
            # print 'DEBUG: p[%s] = q[%s]*(1.-a[%s]) = %s*(1-%s) = %s' % (i, i, i, q[i], a[i], p[i]/ngen)
        p = p / ngen

        # Find largest p
        p_temp = p[:]
        # print 'DEBUG: p_temp: %s' % p_temp
        for y in younglist:
            p_temp[int(y)-1] = -1.
        p_temp = p_temp[::-1]
        # print 'DEBUG: picked: %s' % picked
        for c in picked:
            p_temp[int(c)] = -1.
            # print 'DEBUG: p_temp: %s' % p_temp
            max_p_index = numpy.argmax(p_temp)
            max_p = p_temp[max_p_index]
        # print 'DEBUG: Animal %s had the largest marginal contribution (%s) this round.' % \
        #       (tempped[max_p_index].animalID, max_p)
        contribs[tempped[max_p_index].animalID] = max_p
        picked.append(max_p_index)

        # Delete the pedigree info for the animal with largest q
        # print 'DEBUG: Deleting parent information for animal %s' % (tempped[max_p_index].animalID)
        # print '\t\tWas sire: %s, dam %s' % (tempped[max_p_index].sireID,tempped[max_p_index].damID)
        tempped[max_p_index].sireID = pedobj.kw['missing_parent']       # delete sire in pedobj.pedigree
                                                                        # (forward order)
        tempped[max_p_index].damID = pedobj.kw['missing_parent']        # delete dam in pedobj.pedigree
                                                                        # (forward order)
        # print '\t\tNow sire: %s, dam %s' % (tempped[max_p_index].sireID,tempped[max_p_index].damID)
        ancestors.append(tempped[max_p_index].animalID)     # add the animal with largest q to the
                                                            # list of ancestors
        # print '[DEBUG]: ancestors: %s' % (ancestors)
    # print '[DEBUG]: contribs: %s' % (contribs)
    sum_p_sq = 0.
    for i in contribs.values():
        sum_p_sq = sum_p_sq + (i * i)
    try:
        f_a = 1. / sum_p_sq
    except:
        f_a = 0.0
    if pedobj.kw['messages'] == 'verbose':
        print '='*60
        print 'animals:\t%s' % lenped
        print 'ancestors:\t%s' % n_f
        print 'descendants:\t%s' % n_d
        print 'f_a:\t\t%5.3f' % f_a
        print '='*60
    if pedobj.kw['messages'] == 'verbose':
        print 'generations: %s' % gens
        print '%s animals in generation %s' % (ngen, gen)
        print 'ancestors: %s' % ancestors
        print 'ancestor contributions: %s' % contribs
        print 'DEBUG: f_a: %s' % f_a
    # write some output to a file for later use
    outputfile = '%s%s%s' % (pedobj.kw['filetag'], '_fa_boichard_definite_', '.dat')
    aout = open(outputfile, 'w')
    line = '='*60+'\n'
    aout.write('%s\n' % line)
    aout.write('%s ancestors: %s\n' % (n_f, fs))
    aout.write('%s descendants: %s\n' % (n_d, ds))
    aout.write('generations: %s\n' % gens)
    aout.write('%s animals in generation %s\n' % (ngen, gens[0]))
    aout.write('effective number of ancestors: %s\n' % f_a)
    aout.write('ancestors: %s\n' % ancestors)
    aout.write('ancestor contributions: %s\n' % contribs)
    aout.write('%s\n' % line)
    aout.close()

    if pedobj.kw['debug_messages']:
        logging.info('Exited a_effective_ancestors_definite()')

    return f_a


##
# a_effective_ancestors_indefinite() uses the approach outlined on pages 9 and 10 of
# Boichard et al. (1996) to compute approximate upper and lower bounds for f_a.  This
# is much more tractable for large pedigrees than the exact computation provided in
# a_effective_ancestors_definite().
# NOTE: One problem here is that if you pass a pedigree WITHOUT generations an error
# is not thrown.  You simply end up wth a list of generations that contains the default
# value for NewAnimal() objects, 0.
# NOTE: If you pass a value of n that is greater than the actual number of ancestors in
# the pedigree then strange things happen.  As a stop-gap, a_effective_ancestors_indefinite()
# will detect that case and replace n with the number of founders - 1.
# Boichard's algorithm requires information about the GENERATION of animals.  If you
# do not provide an input pedigree with generations things may not work.  By default
# the most recent generation -- the generation with the largest generation ID -- will
# be used as the reference population.
# @param pedobj A PyPedal pedigree object.
# @param a A numerator relationship matrix (optional).
# @param gen Generation of interest.
# @param n The number of rounds to iterate for solutions.
# @retval The effective ancestor number, or False on failure.
def a_effective_ancestors_indefinite(pedobj, a='', gen='', n=25):
    """
    a_effective_ancestors_indefinite() uses the approach outlined on pages 9 and 10 of
    Boichard et al. (1996) to compute approximate upper and lower bounds for f_a.  This
    is much more tractable for large pedigrees than the exact computation provided in
    a_effective_ancestors_definite().
    :param pedobj: A PyPedal pedigree object.
    :param a: A numerator relationship matrix (optional).
    :param gen: Generation of interest.
    :param n: The number of rounds to iterate for solutions.
    :return: The effective ancestor number, or False on failure.
    """

    if pedobj.kw['debug_messages']:
        logging.info('Entered a_effective_ancestors_indefinite()')

    if not a:
        try:
            a = pyp_nrm.fast_a_matrix(pedobj.pedigree, pedobj.kw)
        except:
            return False
    lenped = len(pedobj.pedigree)  # number of animals in the pedigree file
    # count founders and descendants
    n_f = 0     # number of founders
    n_d = 0     # number of descendants
    fs = []     # list of founders
    ds = []     # list of descendants
    gens = []       # list of generation IDs in the pedigree
    ancestors = []  # list of ancestors already processed
    contribs = {}   # ancestor contributions
    ngen = 0        # number of individuals in the most recent generation
    # Loop through the pedigree quickly and count founders and descendants
    # also, make a list of generations present in the pedigree
    for i in range(lenped):
        g = pedobj.pedigree[i].gen
        if g in gens:
            pass
        else:
            gens.append(g)
    gens.sort()
    # 05/15/2018: This code is supposed to quickly count founders and descendants, but it's not actually doing that.
    #             What it's really doing is counting the number of animals in the earliest generation and calling
    #             them founders, while all other animals are labeled as descendants.
    for i in range(lenped):
        if pedobj.pedigree[i].gen != gens[len(gens)-1]:
            n_f = n_f + 1
            fs.append(int(pedobj.pedigree[i].animalID))
        else:
            n_d = n_d + 1
            ds.append(int(pedobj.pedigree[i].animalID))
    # OK - now we have a list of generations sorted in reverse (descending) order
    ngen = len(gens)
    gens.sort()
    gens.reverse()
    # make a copy of pedobj.pedigree - note that tempped = pedobj.pedigree would only have created a reference to
    # pedobj.pedigree, not an actual separate copy of pedobj.pedigree.
    tempped = pedobj.pedigree[:]
    # now animals are ordered from oldest to youngest in tempped
    tempped.reverse()
    # form q, a vector of that will contain the probabilities of gene origin when we are done
    # We are going to initialize a vector of numpy.zeros to form q, and then we will add ones to q that
    # correspond to members of the youngest generation.
    younglist = []
    q = numpy.zeros([lenped], 'd')
    for i in range(lenped):
        # If the user did not explicitly ask for an analysis of a particular generation then use
        # the most recent generation.
        if not gen:
            gen = gens[0]
        if pedobj.pedigree[i].gen == gen:
            q[i] = 1.
            younglist.append(pedobj.pedigree[i].animalID)
    ngen = len(younglist)
    # loop through the pedigree and form the final version of q (the vector of
    # individual contributions)
    for i in range(lenped):
        if str(tempped[i].sireID) == str(pedobj.kw['missing_parent']) and \
                str(tempped[i].damID) == str(pedobj.kw['missing_parent']):
            # both parents unknown
            pass
        elif str(tempped[i].sireID) == str(pedobj.kw['missing_parent']):
            # sire unknown, dam known
            q[i] = q[tempped[i].animalID-1] * 0.5
            q[tempped[i].damID-1] = q[tempped[i].damID-1] + (0.5 * q[tempped[i].animalID-1])
        elif str(tempped[i].damID) == str(pedobj.kw['missing_parent']):
            # sire known, dam unknown
            q[i] = q[tempped[i].animalID-1] * 0.5
            q[tempped[i].sireID-1] = q[tempped[i].sireID-1] + (0.5 * q[tempped[i].animalID-1])
        else:
            # both parents known
            q[int(tempped[i].sireID)-1] = q[int(tempped[i].sireID)-1] + (0.5 * q[int(tempped[i].animalID)-1])
            q[int(tempped[i].damID)-1] = q[int(tempped[i].damID)-1] + (0.5 * q[int(tempped[i].animalID)-1])
    # divide the elements of q by the number of individuals in the pedigree.  this should
    # ensure that the founder contributions sum to 1.
    for y in younglist:
        q[int(y)-1] = 0.
    q = q / ngen

    # Find largest value of q
    max_p_index = numpy.argmax(q)
    max_p = q[max_p_index]
    contribs[pedobj.pedigree[max_p_index].animalID] = max_p
    picked = []
    picked.append(lenped-max_p_index-1)
    tempped[lenped-max_p_index-1].sireID = pedobj.kw['missing_parent']          # delete sire in pedobj.pedigree
                                                                                # (forward order)
    tempped[lenped-max_p_index-1].damID = pedobj.kw['missing_parent']           # delete dam in pedobj.pedigree
                                                                                # (forward order)
    ancestors.append(pedobj.pedigree[max_p_index].animalID)                     # add the animal with largest q to the
                                                                                # list of ancestors

    # Tricky here... now we have to deal with the fact that we may not want as many contributions are there
    # are ancestors.
    print 'Number of founders in the pedigree file:    ', n
    print 'Number of animals in the pedigree file:     ', lenped
    print 'Number of descendants in the pedigree file: ', n_d
    if n >= (lenped - n_d):
        print '-'*60
        print 'WARNING: (pyp_metrics/a_effective_ancestors_indefinite()): Setting n (%s) to be equal to the actual ' \
              'number of founders (%s) in the pedigree!' % (n, n_f)
        n = n - 1
    if n_f > n:
        _this_many = n
    else:
        _this_many = n_f
    for _t in range(_this_many-1):
        for j in range(n_f-1):
            # form q, the vector of contributions we are going to use
            q = numpy.zeros([lenped], 'd')
            a = numpy.zeros([lenped], 'd')
            for i in range(lenped):
                if pedobj.pedigree[i].gen == gens[0]:
                    q[int(pedobj.pedigree[i].animalID)-1] = 1.
            for j in ancestors:
                a[int(j)-1] = 1.

            # Loop through pedigree to process q
            # q must be processed from YOUNGEST to OLDEST
            for i in range(len(tempped)):
                if str(tempped[i].sireID) == str(pedobj.kw['missing_parent']) and \
                        str(tempped[i].damID) == str(pedobj.kw['missing_parent']):
                    # both parents unknown
                    pass
                elif str(tempped[i].sireID) == str(pedobj.kw['missing_parent']):
                    # sire unknown, dam known
                    q[int(tempped[i].damID)-1] = q[int(tempped[i].damID)-1] + (0.5 * q[int(tempped[i].animalID)-1])
                elif str(tempped[i].damID) == str(pedobj.kw['missing_parent']):
                    # sire known, dam unknown
                    q[int(tempped[i].sireID)-1] = q[int(tempped[i].sireID)-1] + (0.5 * q[int(tempped[i].animalID)-1])
                else:
                    # both parents known
                    q[int(tempped[i].sireID)-1] = q[int(tempped[i].sireID)-1] + (0.5 * q[int(tempped[i].animalID)-1])
                    q[int(tempped[i].damID)-1] = q[int(tempped[i].damID)-1] + (0.5 * q[int(tempped[i].animalID)-1])
            # a must be processed from OLDEST to YOUNGEST
            tempped.reverse()
            for i in range(len(tempped)):
                if str(tempped[i].sireID) == str(pedobj.kw['missing_parent']) and \
                        str(tempped[i].damID) == str(pedobj.kw['missing_parent']):
                    # both parents unknown
                    pass
                elif str(tempped[i].sireID) == str(pedobj.kw['missing_parent']):
                    # sire unknown, dam known
                    a[i] = a[i] + (0.5 * a[int(pedobj.pedigree[i].damID)-1])
                elif str(tempped[i].damID) == str(pedobj.kw['missing_parent']):
                    # sire known, dam unknown
                    a[i] = a[i] + (0.5 * a[int(tempped[i].sireID)-1])
                else:
                    # both parents known
                    a[i] = a[i] + (0.5 * a[int(tempped[i].sireID)-1])
                    a[i] = a[i] + (0.5 * a[int(tempped[i].damID)-1])
            tempped.reverse()
            for y in younglist:
                q[int(y)-1] = 0.

            # Loop through the pedigree to process p
            p = numpy.zeros([lenped], 'd')
            for i in range(lenped):
                p[i] = q[i] * (1. - a[i])
            p = p / ngen

            # Find largest p
            p_temp = p[:]
            for y in younglist:
                p_temp[int(y)-1] = -1.
            p_temp = p_temp[::-1]
            for c in picked:
                p_temp[int(c)] = -1.

            max_p_index = numpy.argmax(p_temp)
            max_p = p_temp[max_p_index]
            contribs[tempped[max_p_index].animalID] = max_p
            picked.append(max_p_index)

            # Delete the pedigree info for the animal with largest q
            tempped[max_p_index].sireID = pedobj.kw['missing_parent']       # delete sire in pedobj.pedigree
                                                                            # (forward order)
            tempped[max_p_index].damID = pedobj.kw['missing_parent']        # delete dam in pedobj.pedigree
                                                                            # (forward order)
            ancestors.append(tempped[max_p_index].animalID)     # add the animal with largest q to the
                                                                # list of ancestors

    # Now compute the upper and lower bounds for f_a based on the founder contributions.
    # print 'DEBUG: n: %s, n_f: %s' % (n, n_f)
    # print 'DEBUG: contribs.values(): %s' % contribs.values()
    n_contrib = len(contribs.values())-1
    _c = 0.
    sum_p_sq = 0.
    # Compute f_u, the upper bound of f_a
    for i in contribs.values():
        if i >= 0.0:
            _c = _c + i
        sum_p_sq = sum_p_sq + (i * i)
    # print 'DEBUG: _c: %s' % _c
    try:
        if n_f <= n_f:
            _df = 1
        else:
            _df = n_f - n
        f_u =  1. / (sum_p_sq + (((1. - _c) ** 2) / _df))
    except:
        f_u = 0.
    # Compute f_l, the lower bound of f_a
    try:
        _denom = 0.
        for i in contribs.values():
            _p_sq = i ** 2
            try:
                _m = ((1. - _c) / i)
            except:
                _m = 0.
            _denom = _denom + (_p_sq + (_m * (i ** 2)))
        f_l = 1. / _denom
    except:
        f_l = 0.
    if pedobj.kw['messages'] == 'verbose':
        print '='*60
        print 'animals:\t%s' % lenped
        print 'founders:\t%s' % n_f
        print 'descendants:\t%s' % n_d
        print 'f_l:\t\t%5.3f' % f_l
        print 'f_u:\t\t%5.3f' % f_u
        print '='*60

    # write some output to a file for later use
    outputfile = '%s%s%s' % (pedobj.kw['filetag'], '_fa_boichard_indefinite_', '.dat')
    aout = open(outputfile, 'w')
    line = '='*60+'\n'
    aout.write('%s\n' % line)
    aout.write('%s founders: %s\n' % (n_f, fs))
    aout.write('%s descendants: %s\n' % (n_d, ds))
    aout.write('generations: %s\n' % gens)
    aout.write('%s animals in generation %s\n' % (ngen, gens[0]))
    aout.write('f_l:\t\t%5.3f' % f_l)
    aout.write('f_u:\t\t%5.3f' % f_u)
    aout.write('ancestors: %s\n' % ancestors)
    aout.write('ancestor contributions: %s\n' % contribs)
    aout.write('%s\n' % line)
    aout.close()

    if pedobj.kw['debug_messages']:
        logging.info('Exited a_effective_ancestors_indefinite()')
    return f_l, f_u


##
# a_coefficients() writes population average coefficients of inbreeding and
# relationship to a file, as well as individual animal IDs and coefficients of
# inbreeding.  Some pedigrees are too large for fast_a_matrix() or fast_a_matrix_r()
# -- an array that large cannot be allocated due to memory restrictions -- and will
# result in a value of -999.9 for all outputs.
# @param pedobj A PyPedal pedigree object.
# @param a A numerator relationship matrix (optional).
# @param method If no relationship matrix is passed, determines which procedure should be called to build one (nrm|frm).
# @retval A dictionary of non-zero individual inbreeding coefficients.
def a_coefficients(pedobj, a='', method='nrm'):
    """
    Write population average coefficients of inbreeding and relationship to a
    file, as well as individual animal IDs and coefficients of inbreeding.  Some
    pedigrees are too large for fast_a_matrix() or fast_a_matrix_r()
    -- an array that large cannot be allocated due to memory restrictions -- and will
    result in a value of -999.9 for all outputs.
    """

    if pedobj.kw['debug_messages']:
        logging.info('Entered a_coefficients()')

    if method not in ['nrm', 'frm']:
        method = 'nrm'
    if pedobj.kw['form_nrm']:
        a = pedobj.nrm.nrm
    if not pedobj.kw['form_nrm'] and not a:
        try:
            if method == 'nrm':
                a = pyp_nrm.fast_a_matrix(pedobj.pedigree, pedobj.kw)
            else:
                a = pyp_nrm.fast_a_matrix_r(pedobj.pedigree, pedobj.kw)
        except:
            return False
    # print a
    lenped = len(pedobj.pedigree)
    f_avg = f_sum = f_n = 0.
    fnz_avg = fnz_sum = fnz_n = 0.
    r_avg = r_sum = r_n = 0.
    rnz_avg = rnz_sum = rnz_n = 0.

    # Populate a dictionary with individual non-zero COI
    individual_coi = {}

    # calculate average coefficients of inbreeding
    for row in range(lenped):
        f_sum = f_sum + a[row, row] - 1.
        f_n = f_n + 1.
    f_avg = f_sum / f_n

    # calculate average non-zero coefficients of inbreeding
    for row in range(lenped):
        if a[row, row] > 1.:
            fnz_sum = fnz_sum + a[row, row] - 1.
            fnz_n = fnz_n + 1.
            individual_coi[pedobj.pedigree[row].animalID] = a[row, row]-1.
    if fnz_sum > 0.:
        fnz_avg = fnz_sum / fnz_n
    else:
        fnz_avg = 0.

    # calculate average coefficients of relationship
    for row in range(lenped):
        for col in range(row):
            r_sum = r_sum + a[row, col]
            r_n = r_n + 1.
    r_avg = r_sum / r_n

    # calculate average non-zero coefficients of relationship
    for row in range(lenped):
        for col in range(row):
            if a[row, col] > 0.:
                rnz_sum = rnz_sum + a[row, col]
                rnz_n = rnz_n + 1.
    if rnz_sum > 0.:
        rnz_avg = rnz_sum / rnz_n
    else:
        rnz_avg = 0.

    # calculate the average relationship between each individual in the population
    # and all other animals in the population and write it to a file
    outputfile2 = '%s%s%s' % (pedobj.kw['filetag'], '_rel_to_pop_', '.dat')
    aout2 = open(outputfile2, 'w')
    line1_2 = '# Average relationship to population (renumbered ID, r)\n'
    for row in range(lenped):
        r_pop_avg = 0.
        for col in range(lenped):
            if row == col:
                pass
            else:
                r_pop_avg = r_pop_avg + a[row, col]
        r_pop_avg = r_pop_avg / lenped
        line = '%s %s\n' % (pedobj.pedigree[row].animalID, r_pop_avg)
        aout2.write(line)
    aout2.close()

    # output population average coefficients
    outputfile = '%s%s%s' % (pedobj.kw['filetag'], '_population_coefficients_', '.dat')
    aout = open(outputfile, 'w')
    line = '='*60+'\n'
    aout.write('%s\n' % line)
    aout.write('# Population average coefficients of inbreeding and relationship\n')
    aout.write('#   f_avg [fnz_avg] = average [nonzero] coefficient of inbreeding\n')
    aout.write('#   f_n [fnz_n] = number of diagonal elements (animals) [>1] in the relationship matrix\n')
    aout.write('#   f_sum [fnz_sum] = sum of [nonzero] coefficients of inbreeding\n')
    aout.write('#   r_avg [rnz_avg] = average [nonzero] coefficient of relationship\n')
    aout.write('#   r_n [rnz_n] = number of [non-zero] elements in the upper off-diagonal of A\n')
    aout.write('#   r_sum [rnz_sum] = sum of [non-zero] elements in the upper off-diagonal of A\n')
    aout.write('f_n: %s\nf_sum: %s\nf_avg: %5.3f\n' % (f_n, f_sum, f_avg))
    aout.write('fnz_n: %s\nfnz_sum: %s\nfnz_avg: %5.3f\n' % (fnz_n, fnz_sum, fnz_avg))
    aout.write('r_n: %s\nr_sum: %s\nr_avg: %5.3f\n' % (r_n, r_sum, r_avg))
    aout.write('rnz_n: %s\nrnz_sum: %s\nrnz_avg: %5.3f\n' % (rnz_n, rnz_sum, rnz_avg))
    aout.write('%s\n' % line)
    aout.close()
    # output individual coefficients of inbreeding
    outputfile = '%s%s%s' % (pedobj.kw['filetag'], '_individual_coefficients_', '.dat')
    aout = open(outputfile, 'w')
    aout.write('# individual coefficients of inbreeding\n')
    aout.write('# animalID f_a\n')
    for row in range(lenped):
        aout.write('%s\t%6.4f\n' % pedobj.pedigree[row].animalID, a[row, row]-1.)
    aout.close()

    if pedobj.kw['debug_messages']:
        logging.info('Exited a_coefficients()')

    return individual_coi


##
# a_fast_coefficients() writes population average coefficients of inbreeding and
# relationship to a file, as well as individual animal IDs and coefficients of
# inbreeding.  It returns a list of non-zero individual CoI.
# @param pedobj A PyPedal pedigree object.
# @param a A numerator relationship matrix (optional).
# @param method If no relationship matrix is passed, determines which procedure should be called to build one (nrm|frm).
# @param debug Print deubgging messages if 1, don't print otherwise.
# @param storage Use dense or sparse matrix storage.
# @retval A dictionary of non-zero individual inbreeding coefficients.
def fast_a_coefficients(pedobj, a='', method='nrm', debug=0, storage='dense'):
    """
    a_fast_coefficients() writes population average coefficients of inbreeding and
    relationship to a file, as well as individual animal IDs and coefficients of
    inbreeding.  It returns a list of non-zero individual CoI.
    """
    try: logging.info('Entered fast_a_coefficients()')
    except: pass
    if method not in ['nrm','frm']:
        method = 'nrm'
    if storage not in ['dense','sparse']:
        storage = 'dense'
    if pedobj.kw['form_nrm']:
        a = pedobj.nrm.nrm
    if not pedobj.kw['form_nrm'] and not a:
        try:
            #if method == 'nrm':
            #    a = pyp_nrm.fast_a_matrix(pedobj.pedigree,pedobj.kw,method=storage)
            #else:
            #    a = pyp_nrm.fast_a_matrix_r(pedobj.pedigree,pedobj.kw,method=storage)
            a = pyp_nrm.fast_a_matrix(pedobj.pedigree,pedobj.kw,method=storage)
        except:
            a = []
            a.append(-999.9)
            return a

    l = len(pedobj.pedigree)
    f_avg = f_sum = f_n = 0.
    fnz_avg = fnz_sum = fnz_n = 0.
    r_avg = r_sum = r_n = 0.
    rnz_avg = rnz_sum = rnz_n = 0.

    # Populate a dictionary with individual non-zero COI
    individual_coi = {}

    for row in range(l):

        # Do inbreeding things here in the outer loop
        f_sum = f_sum + ( a[row,row] - 1. )
        f_n = f_n + 1
        if ( a[row,row] > 1. ) :
                fnz_sum = fnz_sum + ( a[row,row] - 1. )
                fnz_n = fnz_n + 1
                individual_coi[pedobj.pedigree[row].animalID] = a[row,row]-1.

        for col in range(row):
        # Do relationship things here in the inner loop
            r_sum = r_sum + a[row,col]
            r_n = r_n + 1
            if ( a[row,col] > 0. ):
                rnz_sum = rnz_sum + a[row,col]
                rnz_n = rnz_n + 1

    if f_sum > 0.:
        f_avg = f_sum / f_n
    else:
        f_avg = 0.
    if fnz_sum > 0.:
        fnz_avg = fnz_sum / fnz_n
    else:
        fnz_avg = 0.
    if r_sum > 0.:
        r_avg = r_sum / r_n
    else:
        r_avg = 0.
    if rnz_sum > 0.:
        rnz_avg = rnz_sum / rnz_n
    else:
        rnz_avg = 0.

    if pedobj.kw['file_io']:

        # calculate the average relationship between each individual in the population
        # and all other animals in the population and write it to a file
        outputfile2 = '%s%s%s' % (pedobj.kw['filetag'],'_rel_to_pop_','.dat')
        aout2 = open(outputfile2,'w')
        line1_2 = '# Average relationship to population (renumbered ID, r)\n'
        for row in range(l):
            r_pop_avg = 0.
            for col in range(l):
                if ( row == col ):
                    pass
                else:
                    r_pop_avg = r_pop_avg + a[row,col]
            r_pop_avg = r_pop_avg / l
            line = '%s %s\n' % (pedobj.pedigree[row].animalID,r_pop_avg)
            aout2.write(line)
        aout2.close()
        # output population average coefficients
        outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_population_coefficients_','.dat')
        aout = open(outputfile,'w')
        line = '='*60+'\n'
        aout.write('%s\n' % line)
        aout.write('# Population average coefficients of inbreeding and relationship (fast_a_coefficients)\n')
        aout.write('#   f_avg [fnz_avg] = average [nonzero] coefficient of inbreeding\n')
        aout.write('#   f_n [fnz_n] = number of diagonal elements (animals) [>1] in the relationship matrix\n')
        aout.write('#   f_sum [fnz_sum] = sum of [nonzero] coefficients of inbreeding\n')
        aout.write('#   r_avg [rnz_avg] = average [nonzero] coefficient of relationship\n')
        aout.write('#   r_n [rnz_n] = number of [non-zero] elements in the upper off-diagonal of A\n')
        aout.write('#   r_sum [rnz_sum] = sum of [non-zero] elements in the upper off-diagonal of A\n')
        aout.write('f_n: %s\nf_sum: %s\nf_avg: %5.3f\n' % (f_n,f_sum,f_avg))
        aout.write('fnz_n: %s\nfnz_sum: %s\nfnz_avg: %5.3f\n' % (fnz_n,fnz_sum,fnz_avg))
        aout.write('r_n: %s\nr_sum: %s\nr_avg: %5.3f\n' % (r_n,r_sum,r_avg))
        aout.write('rnz_n: %s\nrnz_sum: %s\nrnz_avg: %5.3f\n' % (rnz_n,rnz_sum,rnz_avg))
        aout.write('%s\n' % line)
        aout.close()
        # output individual coefficients of inbreeding
        outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_individual_coefficients_','.dat')
        aout = open(outputfile,'w')
        aout.write('# individual coefficients of inbreeding\n')
        aout.write('# animalID f_a\n')
        for row in range(l):
            aout.write('%s\t%6.4f\n' % (pedobj.pedigree[row].animalID,a[row,row]-1.))
        aout.close()

    try: logging.info('Exited fast_a_coefficients()')
    except: pass
    return individual_coi

##
# theoretical_ne_from_metadata() computes the theoretical effective population
# size based on the number of sires and dams contained in a pedigree metadata
# object.  Writes results to an output file.
# @param pedobj A PyPedal pedigree object.
# @retval True (1) on success, false (0) on failure
def theoretical_ne_from_metadata(pedobj):
    """
    theoretical_ne_from_metadata() computes the theoretical effective population
    size based on the number of sires and dams contained in a pedigree metadata
    object.  Writes results to an output file.
    """
    try: logging.info('Entered theoretical_ne_from_metadata()')
    except: pass
    try:
        ns = float(pedobj.metadata.num_unique_sires)
        nd = float(pedobj.metadata.num_unique_dams)
        ne = 1. / ( (1./(4.*ns)) + (1./(4.*nd)) )
        outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_ne_from_metadata_','.dat')
        aout = open(outputfile,'w')
        line = '='*60+'\n'
        aout.write('%s' % line)
        aout.write('# Theoretical effective population size (N_e)\n')
        aout.write('#   n_sires = number of sires\n')
        aout.write('#   n_dams = number of dams\n')
        aout.write('#   n_e = effective population size\n')
        aout.write('n_sires: %s\n' % ns)
        aout.write('n_dams: %s\n' % nd)
        aout.write('n_e: %s\n' % ne)
        aout.write('%s\n' % line)
        aout.close()
        _return = 1
    except:
        _return = 0

    try: logging.info('Exited theoretical_ne_from_metadata()')
    except: pass
    return _return

##
# pedigree_completeness() computes the proportion of known ancestors in the pedigree of
# each animal in the population for a user-determined number of generations.    Also,
# the mean pedcomps for all animals and for all animals that are not founders are
# computed as summary statistics.
# @param pedobj A PyPedal pedigree object.
# @param gens The number of generations the pedigree should be traced for completeness.
# @retval Dictionary of summary statistics
def pedigree_completeness(pedobj, gens=4):
    """
    pedigree_completeness() computes the proportion of known ancestors in the
    pedigree of each animal in the population for a user-determined number of
    generations. Also, the mean pedcomps for all animals and for all animals
    that are not founders are computed as summary statistics.
    """
    try: logging.info('Entered pedigree_completeness()')
    except: pass
    l = len(pedobj.pedigree)
    #print l
    c_summary = {}
    mp = str(pedobj.kw['missing_parent'])

    # All-animal summary stats
    c_max = 0.0
    c_min = 1.0
    c_cnt = 0
    c_sum = 0
    # Non-founder summary stats
    nf_c_max = 0.0
    nf_c_min = 1.0
    nf_c_cnt = 0
    nf_c_sum = 0

    for i in range(l):
        animalid = int(pedobj.pedigree[i].animalID)
        sireid = int(pedobj.pedigree[i].sireID)
        damid = int(pedobj.pedigree[i].damID)
        animalname = pedobj.pedigree[i].name
        sirename = pedobj.pedigree[pedobj.pedigree[i].sireID-1].name
        damname = pedobj.pedigree[pedobj.pedigree[i].damID-1].name
        if pedobj.kw['debug_messages']:
            print 'Animal: %s (%s)' % (animalid, animalname)
            if str(sireid) != mp:
                print '\tSire: %s (%s)' % (sireid, sirename)
            if str(damid) != mp:
                print '\tDam: %s (%s)' % (damid, damname)
        # Founders are easy to deal with!
        if pedobj.pedigree[i].founder == 'y':
            pedobj.pedigree[i].pedcomp = 0.0
        if str(sireid) == mp and str(damid) == mp:
            _compl = 0.0
        else:
            # Use pyp_nrm/recurse_pedigree_n()...
            _sire_ped = []
            n_max_ancestors = 2 * ( ( 2 ** gens ) - 1 )
            _sire_ped = pyp_nrm.recurse_pedigree_n(pedobj,sireid,_sire_ped,gens-1)
            _l_sire_ped = len(_sire_ped)
            #print '\tLen sire ped: ', len(_sire_ped)#, _sire_ped
            _dam_ped = []
            _dam_ped = pyp_nrm.recurse_pedigree_n(pedobj,damid,_dam_ped,gens-1)
            _l_dam_ped = len(_dam_ped)
            #print '\tLen dam ped: ', len(_dam_ped)#, _dam_ped
            #print '\tMax ancestors: ', n_max_ancestors
            _compl = float( _l_sire_ped + _l_dam_ped ) / float( n_max_ancestors )
        #if pedobj.kw['debug_messages']:
        #print '\tPedigree completeness: %s' % (_compl)
        pedobj.pedigree[i].pedcomp = _compl

        # Summary statistics
        c_sum = c_sum + _compl
        c_cnt = c_cnt + 1
        if _compl > c_max:
            c_max = _compl
        if _compl < c_min:
            c_min = _compl

        c_summary['sum'] = c_sum
        c_summary['n'] = c_cnt
        c_summary['min'] = c_min
        c_summary['max'] = c_max
        c_rng = c_max - c_min
        c_summary['range'] = c_rng
        c_avg = float(c_sum) / float(c_cnt)
        c_summary['average'] = c_avg

        # Non-founder summary stats
        if pedobj.pedigree[i].founder != 'y':
            nf_c_sum = nf_c_sum + _compl
            nf_c_cnt = nf_c_cnt + 1
            if _compl > nf_c_max:
                    nf_c_max = _compl
            if _compl < nf_c_min:
                    nf_c_min = _compl

    c_summary['nonfounder_sum'] = nf_c_sum
    c_summary['nonfounder_n'] = nf_c_cnt
    c_summary['nonfounder_min'] = nf_c_min
    c_summary['nonfounder_max'] = nf_c_max
    nf_c_rng = nf_c_max - nf_c_min
    c_summary['nonfounder_range'] = nf_c_rng
    nf_c_avg = float(nf_c_sum) / float(nf_c_cnt)
    c_summary['nonfounder_average'] = nf_c_avg

    if pedobj.kw['messages'] == 'verbose':
        print '-'*80
        print 'Pedigree completeness summary statistics (all animals):'
        print '\tN\t%s' % (c_cnt)
        print '\tSum\t%s' % (c_sum)
        print '\tMean\t%s' % (c_avg)
        print '\tMin\t%s' % (c_min)
        print '\tMax\t%s' % (c_max)
        print '\tRange\t%s' % (c_rng)
        print '-'*80
        print 'Pedigree completeness summary statistics (non-founders):'
        print '\tN\t%s' % (nf_c_cnt)
        print '\tSum\t%s' % (nf_c_sum)
        print '\tMean\t%s' % (nf_c_avg)
        print '\tMin\t%s' % (nf_c_min)
        print '\tMax\t%s' % (nf_c_max)
        print '\tRange\t%s' % (nf_c_rng)
        print '-'*80

    try: logging.info('Exited pedigree_completeness()')
    except: pass
    return c_summary

##
# common_ancestors() returns a list of the ancestors that two animals share in common.
# @param anim_a The renumbered ID of the first animal, a.
# @param anim_b The renumbered ID of the second animal, b.
# @param pedobj A PyPedal pedigree object.
# @retval A list of animals related to anim_a AND anim_b
def common_ancestors(anim_a, anim_b, pedobj):
    """
    common_ancestors() returns a list of the ancestors that two animals share in common.
    """
    #try: logging.info('Entered common_ancestors()')
    #except: pass
    #print anim_a
    #print anim_b
    ped_a = related_animals(anim_a,pedobj)
    #print 'ped_a %s: %s' % ( anim_a, ped_a )
    ped_b = related_animals(anim_b,pedobj)
    #print 'ped_b %s: %s' % ( anim_b, ped_b )
    shared = []
    try:
        ped_a.sort()
        ped_b.sort()
        # I know that this is a bad way to do this!
        for _a in ped_a:
            if _a in ped_b:
                shared.append(_a)
    except:
        pass

    #try: logging.info('Exited common_ancestors()')
    #except: pass
    return shared

##
# related_animals() returns a list of the ancestors of an animal.
# @param anim The renumbered ID of an animal, a.
# @param pedobj A PyPedal pedigree object.
# @retval A list of animals related to anim_a
def related_animals(anim, pedobj):
    """
    related_animals() returns a list of the ancestors of an animal.
    """
    # Use dictionaries for the lookups!
    #try: logging.info('Entered related_animals()')
    #except: pass
    _ped = []
    try:
        _ped = pyp_nrm.recurse_pedigree_idonly(pedobj, anim,_ped)
    except:
        pass
    #try: logging.info('Exited related_animals()')
    #except: pass
    #print '_ped %s: %s' % ( anim, _ped )
    return _ped

##
# relationship() returns the coefficient of relationship for two
# animals, anim_a and anim_b.
# @param anim_a The renumbered ID of an animal, a.
# @param anim_b The renumbered ID of an animal, b.
# @param pedobj A PyPedal pedigree object.
# @param renumber Renumber the pedigree if it has not yet been renumbered.
# @retval The coefficient of relationship of anim_a and anim_b
def relationship(anim_a, anim_b, pedobj, renumber=False):
    """
    relationship() returns the coefficient of relationship for two
    animals, anim_a and anim_b.
    """
    if pedobj.kw['pedigree_is_renumbered'] == 0 and renumber == False:
        if pedobj.kw['messages'] != 'quiet':
            print '[WARNING]: The pedigree you passed to pyp_metrics/relationship() is not renumbered; this may result in incorrect calculations!'
        try: logging.warning('The pedigree you passed to pyp_metrics/relationship() is not renumbered; this may result in incorrect calculations!')
        except: pass
    elif pedobj.kw['pedigree_is_renumbered'] == 0 and renumber == True:
        if pedobj.kw['messages'] != 'quiet':
            print '[INFO]: Renumbering the pedigree in pyp_metrics/relationship().'
        try: logging.info('Renumbering the pedigree in pyp_metrics/relationship().')
        except: pass
        pedobj.kw['renumber'] = 1
        pedobj.renumber()
    else:
        pass
    _r = 0.0                                                # This is the default
    if anim_a == anim_b:
        _r = 1.0
    else:
        # If there is an NRM attached to the pedigree then use it.
        try:
            if pedobj.kw['form_nrm'] and pedobj.nrm.nrm.shape[0] == pedobj.metadata.num_records:
                _r = pedobj.nrm.nrm[int(anim_a)-1][int(anim_b)-1]
        # If the NRM lookup fails go ahead and calculate the relationship the long way.
        except:
            # Extract the pedigree for each animal, merge the two, and
            # compute the coefficient of relationship.
            try:
                # This is not very efficient, but let's get the bugs fixed
                # before we tty anything fancy.
                _ped_a, _ped_b, _ped, _seen = [], [], [], {}
                _ped_a = pyp_nrm.recurse_pedigree(pedobj,anim_a,_ped_a)
                _ped_b = pyp_nrm.recurse_pedigree(pedobj,anim_b,_ped_b)
                # The dictionary _ped_a tracks "seen" animals
                for _a in _ped_a:
                    try:
                        _seen[_a.animalID]
                    except KeyError:
                        _ped.append(_a)
                        _seen[_a.animalID] = _a.animalID
                for _b in _ped_b:
                    try:
                        _seen[_b.animalID]
                    except KeyError:
                        _ped.append(_b)
                        _seen[_b.animalID] = _b.animalID
                del(_ped_a)
                del(_ped_b)
                _tag = '%s' % (pedobj.kw['filetag'])
                _reord = []
                for j in range(len(_ped)):
                    _reord.append(copy.deepcopy(_ped[j-1]))
                if pedobj.kw['slow_reorder']:
                    _reord = pyp_utils.reorder(_reord,_tag,debug=pedobj.kw['debug_messages'])
                else:
                    _reord = pyp_utils.fast_reorder(_reord,_tag)
                _s, _map = pyp_utils.renumber(_reord,_tag, returnmap=1, \
                    debug=pedobj.kw['debug_messages'])
                _backmap = {}
                for _mk, _mv in _map.iteritems():
                    _backmap[_mv] = _mk
                _opts = copy.deepcopy(pedobj.kw)
                _opts['filetag'] = _tag
                if pedobj.kw['nrm_method'] == 'nrm':
                    _a = pyp_nrm.fast_a_matrix(_s,_opts)
                else:
                    _a = pyp_nrm.fast_a_matrix_r(_s,_opts)
                _r = _a[_map[anim_a]-1][_map[anim_b]-1]
            # If we can't calculate the relationship for some reason return a 0.
            except:
                try: logging.warn('Could not compute the relationship between animals %s and %s; defaulting to 0.0', anim_a, anim_b)
                except: pass
    return _r

##
# mating_coi() returns the coefficient of inbreeding of offspring of a
# mating between two animals, anim_a and anim_b.
# @param anim_a The renumbered ID of an animal, a.
# @param anim_b The renumbered ID of an animal, b.
# @param pedobj A PyPedal pedigree object.
# @param gens The number of generations from the pedigree to be used for calculating CoI.  By default, gens=0.
# @retval The coefficient of inbreeding of the offpsring of anim_a and anim_b
def mating_coi(anim_a, anim_b, pedobj, gens=0):
    """
    mating_coi() returns the coefficient of inbreeding of offspring of a
    mating between two animals, anim_a and anim_b.
    """
    try: logging.info('Entered mating_coi()')
    except: pass
    gens = int(gens)
    _f = -999.9
    #print 'gens: %d' % ( gens )
    if anim_a == anim_b:
        _f =  1.0
    else:
        # This is the original method that does not require than an animal
        # be added to, and then deleted from, the pedigree.
        if gens == -1:
            try:
                # This is simple -- the CoI of an animal is one-half of the
                # relationship between sire and dam.
                _r = relationship(anim_a,anim_b,pedobj)
                _f = 0.5 * _r
            except:
                _f = 0.0
        # This is the new (as of 04/10/2006) method that uses a new dummy
        # animal in the pedigree.
        elif gens >= 0:
            try: logging.warning('Using the new algorithm in pyp_metrics/mating_coi().')
            except: pass
            #print 'Using the new algorithm in pyp_metrics/mating_coi().'
            # Add the hypothetical animal to the pedigree.
            _newid = max(pedobj.idmap.keys())+1
            #print '_newid: ', _newid
            #print 'ID: ', pedobj.pedigree[-1].animalID, '\tName: ', \
                #pedobj.pedigree[-1].name, '\tSire:', pedobj.pedigree[-1].sireName, \
                #'\tDam:', pedobj.pedigree[-1].damName
            _added = pedobj.addanimal(_newid,anim_a,anim_b)
            #print '_added: ', _added
            #print pedobj.pedigree[-1].animalID, pedobj.pedigree[-1].name
            if _added:
                # Most of this code was lifted verbatim from
                # pyp_nrm/inbreeding_vanraden(). That's where you
                # should look for insight and commentary.
                #print '\tConverting pedigree to graph.'
                ng = pyp_network.ped_to_graph(pedobj)
                _ped, top_ped = [], []
                # If gens > 0 then we need to limit the pedigree to
                # that number of generations.
                if int(gens) > 0:
                    top_peddict = pyp_network.find_ancestors_g(ng, len(pedobj.idmap), {}, gens)
                    top_peddict[len(pedobj.idmap)] = 1
                    top_ped = top_peddict.keys()
                    top_r, _anids = [], []
                    for _j in top_ped:
                        if top_peddict[_j] <= gens:
                            top_r.append(copy.deepcopy(pedobj.pedigree[int(_j)-1]))
                            if top_peddict[_j] == 1:
                                top_r[-1].sireID = pedobj.kw['missing_parent']
                                top_r[-1].damID = pedobj.kw['missing_parent']
                            _anids.append(top_r[-1].animalID)
                else:
                    _anids = pedobj.backmap.keys()
                #print _anids
                # Now we need to actually build the subpedigree for
                # the hypothetical animal.
                if int(gens) > 0:
                    _ped = top_peddict
                else:
                    #pedobj.pedigree[-1].printme()
                    #print '\t\tFinding ancestors of %s.' % \
                        #( pedobj.pedigree[-1].animalID )
                    _ped = pyp_network.find_ancestors(ng, \
                        pedobj.pedigree[-1].animalID, [])
                    _ped.append(int(pedobj.pedigree[-1].animalID))
                    #print _ped
                    #print '\t\tThe pedigree for animal %s (orig ID: %s) has %s entries.' % ( pedobj.pedigree[-1].animalID, \
                        #pedobj.pedigree[-1].originalID, len(_ped))
                if int(gens) > 0:
                    _r = top_r
                else:
                    _r = []
                    _map = {}
                    #print '\t\tCopying pedigree for %s.' % ( pedobj.pedigree[-1].originalID )
                    for j in _ped:
                        _r.append(copy.deepcopy(pedobj.pedigree[int(j)-1]))
                #print '\t\tReordering pedigree for %s.' % ( pedobj.pedigree[-1].animalID )
                _tag = '%s_%s' % (pedobj.kw['filetag'],pedobj.pedigree[-1].animalID)
                if pedobj.kw['slow_reorder']:
                    _r = pyp_utils.reorder(_r,_tag)
                else:
                    _r = pyp_utils.fast_reorder(_r,_tag)
                #print '\t\tRenumbering pedigree for %s.' % ( pedobj.pedigree[-1].animalID )
                _s, _map = pyp_utils.renumber(_r,_tag, returnmap=1, debug=pedobj.kw['debug_messages'])
                _backmap = {}
                for _mk, _mv in _map.iteritems():
                    _backmap[_mv] = _mk
                _opts = copy.deepcopy(pedobj.kw)
                _opts['filetag'] = _tag
                #print '\t\tForming the A matrix for %s\'s pedigree.' % ( pedobj.pedigree[-1].animalID )
                if pedobj.kw['nrm_method'] == 'nrm':
                    _a = pyp_nrm.fast_a_matrix(_s,_opts)
                else:
                    _a = pyp_nrm.fast_a_matrix_r(_s,_opts)
                #print _map
                #print _a
                _f = _a[_map[pedobj.pedigree[-1].animalID]-1][_map[pedobj.pedigree[-1].animalID]-1] - 1.
                #print '_r: ', _r
                #print '_f: ', _f
                # Cleanup
                del(_r); del(_map); del(_backmap); del(_a)
                del(_s); del(ng); del(_ped); del(top_ped)
                try: del(top_peddict)
                except: pass
                try: del(top_r)
                except: pass
                try: del(anids)
                except: pass
                # Now that we're done, delete the hypothetical animal
                # from the pedigree.
                pedobj.delanimal(_newid)
        else:
            pass
    try: logging.info('Exited mating_coi()')
    except: pass
    return _f

##
# mating_coi_group() returns the coefficients of inbreeding of offspring
# of a series of matings, as well as a list of minimum-inbreeding matings.
# @param matings A list of proposed matings
# @param pedobj A PyPedal pedigree object.
# @param names Indicates if the identifiers in 'matings' are names or animalIDs
# @param gens The number of generations from the pedigree to be used for calculating CoI.  By default, gens=0.
# @retval The coefficient of inbreeding of the offpsring of anim_a and anim_b
def mating_coi_group(matings, pedobj, names=0, gens=0):
    _results = {}
    fx = {}
    metadata = {}
    matingf = {}
    try:
        if names == 1:
            _matings = []
            for m in matings:
                msplit = m.split('_')
                # Note that we need to map from names to original IDs
                # and from original IDs to renumbered IDs
                _matings[pedobj.namemap[k]] = pedobj.namemap[v]
                _matings.append('%s_%s'%(pedobj.namemap[msplit[0]], pedobj.namemap[msplit[1]]))
        else:
            _matings = matings
        #if pedobj.kw['messages'] == 'verbose':
            #print 'S\tD\tf'

        f_min=1.;f_max=-1.;f_mean=0.;f_range=0.;f_sum=0.;f_n = 0.
        f_min_nz=1.;f_max_nz=-1.;f_mean_nz=0.;f_range_nz=0.;f_sum_nz=0.;f_n_nz=0.

        #for k,v in _matings.iteritems():
        for m in matings:
            msplit = m.split('_')
            k = msplit[0]
            v = msplit[1]
            if not fx.has_key(k): fx[k] = {}
            coi = mating_coi(k,v,pedobj,gens)
            fx[k][v] = coi
            # Accumulate summary statistics
            f_sum = f_sum + coi
            f_n = f_n + 1
            if coi < f_min: f_min = coi
            if coi > f_max: f_max = coi
            # Stats for non-zero CoI
            if coi > 0.:
                f_sum_nz = f_sum_nz + coi
                f_n_nz = f_n_nz + 1
                if coi < f_min_nz: f_min_nz = coi
                if coi > f_max_nz: f_max_nz = coi
            #if pedobj.kw['messages'] == 'verbose':
                #if names == 1:
                    #print pedobj.namebackmap[k], '\t', pedobj.namebackmap[v], \
                        #'\t', coi
                #else:
                    #print k, '\t', v, '\t', coi
            matingf[m] = coi
        f_range = f_max - f_min
        f_range_nz = f_max_nz - f_min_nz
        try: f_mean = f_sum / f_n
        except: f_mean = 0.
        try: f_mean_nz = f_sum_nz / f_n_nz
        except: f_mean_nz = 0.
        # Summary statistics including all CoI
        metadata['all'] = {}
        metadata['all']['f_n'] = f_n
        metadata['all']['f_sum'] = f_sum
        metadata['all']['f_min'] = f_min
        metadata['all']['f_max'] = f_max
        metadata['all']['f_range'] = f_range
        metadata['all']['f_mean'] = f_mean
        # Summary statistics including only nonzero CoI
        metadata['nonzero'] = {}
        metadata['nonzero']['f_n'] = f_n_nz
        metadata['nonzero']['f_sum'] = f_sum_nz
        metadata['nonzero']['f_min'] = f_min_nz
        metadata['nonzero']['f_max'] = f_max_nz
        metadata['nonzero']['f_range'] = f_range_nz
        metadata['nonzero']['f_mean'] = f_mean_nz
        _results['metadata'] = metadata
        #_results['fx'] = fx
        _results['matings'] = matingf
        #print '_results: ', _results
    except:
        pass
    return _results

##
# effective_founder_genomes() simulates the random segregation of founder alleles through a pedigree.
# At present only two alleles are simulated for each founder.  Summary statistics are
# computed on the most recent generation.
# @param pedobj A PyPedal pedigree object.
# @param rounds The number of times to simulate segregation through the entire pedigree.
# @param chrometype The type of chromosome to simulate ('autosome'|'sex').
# @param heterogametic The heterogametic sex in your species ('f'|'m').
# @param quiet If False output is shown; if True output is not shown.
# @retval The effective number of founder genomes based on 'rounds' of gene-drop simulations.
def effective_founder_genomes(pedobj, rounds=10, chrometype='autosome', heterogametic='m', quiet=False):
    """
    effective_founder_genomes() simulates the random segregation of founder alleles
    through a pedigree.  At present only two alleles are simulated for each founder.
    Summary statistics are computed on the most recent generation.
    """
    logging.info('Entered effective_founder_genomes()')
    #except: pass
    #print '[DEBUG]: rounds: %s' % ( rounds )
    if rounds < 1:
        if pedobj.kw['debug_messages']:
            print '[ERROR]: The rounds parameter in pyp_metrics/gene_drop() must be positive (>0)!'
        rounds = 1
    # Check chrometype values for validity; default to 'autosome'
    if chrometype not in ['autosome','sex']:
        try: logging.warning('You provided an unrecognized value of the chrometype parameter, %s, in effective_founder_genomes(); defaulting to \'autosome\'',chrometype)
        except: pass
        chrometype = 'autosome'
    # Check heterogametic sex; default to 'm'ale  as do mammals
    if heterogametic not in ['m','f']:
        try: logging.warning('You provided an unrecognized value of the heterogametic parameter, %s, in effective_founder_genomes(); defaulting to \'m\'',heterogametic)
        except: pass
        heterogametic = 'm'
    # Build a list of generations
    l = len(pedobj.pedigree)
    # Build a list of unique generations so that we can find the most recent generation
    gens = []       # list of generation IDs in the pedigree
    ngen = 0    # Number of animals in the latest generation
    nfounders = 0
    n_g = 0.0   # Effective number of founder genomes
    for j in range(l):
        g = pedobj.pedigree[j].gen
        if g in gens:
            pass
        else:
            gens.append(g)
        if pedobj.pedigree[j].founder == 'y':
            nfounders = nfounders + 1
    gens.sort()
    allele_freqs = {}
    summary_freqs = {}
    summary_freqs['distinct_alleles'] = {}
    summary_stats = {}
    summary_stats['distinct_alleles'] = {}
    outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_gene_drop','.out')
    myline = '='*80
    myline2 = '*'*80
    for r in xrange(rounds):
        #print '[DEBUG]: Round %s in effextive_founder_genomes()' % ( r )
        random.seed()
        allele_freqs = {}
        for i in range(l):
            if pedobj.pedigree[i].founder == 'y':
                if pedobj.kw['debug_messages']:
                    print '[DEBUG]: Alleles for founder %s:\t%s' % (pedobj.pedigree[i].animalID,pedobj.pedigree[i].alleles)
            else:
                _error = 0
                if pedobj.pedigree[int(pedobj.pedigree[i].sireID)-1].alleles == ['','']:
                    if pedobj.kw['debug_messages']:
                        print '[ERROR]: The sire (%s) of animal %s does not have usable alleles: %s!' % (pedobj.pedigree[i].animalID,pedobj.pedigree[i].sireID,pedobj.pedigree[int(pedobj.pedigree[i].sireID)-1].alleles)
                    _error = 1
                if pedobj.pedigree[int(pedobj.pedigree[i].damID)-1].alleles == ['','']:
                    if pedobj.kw['debug_messages']:
                        print '[ERROR]: The dam (%s) of animal %s does not have usable alleles: %s!' % (pedobj.pedigree[i].animalID,pedobj.pedigree[i].damID,pedobj.pedigree[int(pedobj.pedigree[i].damID)-1].alleles)
                    _error = 1
                if _error == 1:
                    return 0
                # Pick a sire allele at random
                #random.seed(random.randint(1,1000000000))
                random.seed()
                _rs = random.random()
                if _rs < 0.5:
                    _as = pedobj.pedigree[int(pedobj.pedigree[i].sireID)-1].alleles[0]
                else:
                    _as = pedobj.pedigree[int(pedobj.pedigree[i].sireID)-1].alleles[1]
                # Pick a dam allele at random
                #random.seed(random.randint(1,1000000000))
                random.seed()
                _rd = random.random()
                if _rd < 0.5:
                    _ad = pedobj.pedigree[int(pedobj.pedigree[i].damID)-1].alleles[0]
                else:
                    _ad = pedobj.pedigree[int(pedobj.pedigree[i].damID)-1].alleles[1]
                if pedobj.kw['debug_messages']:
                    print '[DEBUG]: Animal %s (g: %s) (s:%s,d:%s) got sire allele %s and dam allele %s' % (pedobj.pedigree[i].animalID,pedobj.pedigree[i].gen,pedobj.pedigree[i].sireID,pedobj.pedigree[i].damID,_as,_ad)
                pedobj.pedigree[i].alleles = [_as,_ad]
                if pedobj.pedigree[i].gen == gens[-1:][0]:
                    try:
                        allele_freqs[pedobj.pedigree[i].alleles[0]] = allele_freqs[pedobj.pedigree[i].alleles[0]] + 1
                    except KeyError:
                        allele_freqs[pedobj.pedigree[i].alleles[0]] = 1
                    try:
                        allele_freqs[pedobj.pedigree[i].alleles[1]] = allele_freqs[pedobj.pedigree[i].alleles[1]] + 1
                    except KeyError:
                        allele_freqs[pedobj.pedigree[i].alleles[1]] = 1
                    #print '[DEBUG]: Animal %s (g: %s) (s:%s,d:%s) got sire allele %s and dam allele %s' % (pedobj.pedigree[i].animalID,pedobj.pedigree[i].gen,pedobj.pedigree[i].sireID,pedobj.pedigree[i].damID,_as,_ad)
        # Sumamarize allele data
        nalleles = len(allele_freqs.keys()) # Number of distinct alleles in latest generation
        allelecount = 0
        for k in allele_freqs.keys():
            allelecount = allelecount + int(allele_freqs[k])
        _ngs = 0.0
        for k in allele_freqs.keys():
            _freq = float(allele_freqs[k]) / float(allelecount)
            _ngs = _ngs + _freq**2
        _ng = 1. / (2.*_ngs)

        # Update summary dictionary
        try:
            summary_freqs['allele_count'] = summary_freqs['allele_count'] + allelecount
        except KeyError:
            summary_freqs['allele_count'] = allelecount
        try:
            summary_freqs['distinct_allele_count'] = summary_freqs['distinct_allele_count'] + nalleles
        except KeyError:
            summary_freqs['distinct_allele_count'] = nalleles
        for k in allele_freqs.keys():
            try:
                summary_freqs['distinct_alleles'][k] = summary_freqs['distinct_alleles'][k] + allele_freqs[k]
            except KeyError:
                summary_freqs['distinct_alleles'][k] = allele_freqs[k]
        try:
            summary_freqs['n_g'] = ( _ng + summary_freqs['n_g'] ) / 2.
        except KeyError:
            summary_freqs['n_g'] = _ng

        # Handle writing the output file
        if r == 0:
            dout = open(outputfile,'w')
            mname = '# FILE: %s\n' % (outputfile)
            dout.write(mname)
            dout.write('# Results from %s-round PyPedal gene-drop simulation.\n'%(rounds))
        else:
            dout = open(outputfile,'a')
        dout.write('%s\n'%(myline))
        dout.write('Allele frequency data from gene drop simulation, round %s\n' % (r+1))
        dout.write('\tNumber of distinct alleles: %s\n' % (nalleles))
        allelecount = 0
        for k in allele_freqs.keys():
            allelecount = allelecount + int(allele_freqs[k])
        dout.write('\tNumber of alleles in latest generation: %s\n' % (allelecount))
        _ngs = 0.0
        for k in allele_freqs.keys():
            _freq = float(allele_freqs[k]) / float(allelecount)
            dout.write('\t\tAllele %s:\t%s (%s)\n' % (k,_freq,_freq**2))
            _ngs = _ngs + _freq**2
        _ng = 1. / (2.*_ngs)
        dout.write('\tEffective number of founder genomes: %s\n' % (_ng))
        dout.close()

    # Compute summary stats from summary dictionary
    summary_stats['allele_count'] = float(summary_freqs['allele_count']) / float(rounds)
    summary_stats['distinct_allele_count'] = float(summary_freqs['distinct_allele_count']) / float(rounds)
    summary_stats['freq_sum'] = 0.0
    for k in summary_freqs['distinct_alleles'].keys():
        summary_stats['freq_sum'] = summary_stats['freq_sum'] + summary_freqs['distinct_alleles'][k]
    for k in summary_freqs['distinct_alleles'].keys():
        summary_stats['distinct_alleles'][k] = summary_freqs['distinct_alleles'][k] / float(summary_stats['freq_sum'])
    summary_stats['n_g'] = summary_freqs['n_g']

    # Print summary statistics to screen at the end of the last round
    if pedobj.kw['messages'] == 'verbose'  and quiet == False:
        print '*'*100
        print 'Summary statistics from %s-round gene-drop simulation' % (rounds)
        print '\tNumber of distinct founder alleles: %s' % (2*nfounders)
        print '\tMean allele count in latest generation: %s' % (summary_stats['allele_count'])
        print '\tMean number of distinct alleles in latest generation: %s' % (summary_stats['distinct_allele_count'])
        print '\tFrequency of distinct alleles sampled:'
        for k in summary_stats['distinct_alleles'].keys():
            print '\t\tAllele %s:\t%s (%s)' % (k,summary_stats['distinct_alleles'][k],summary_stats['distinct_alleles'][k]**2)
        print '\tMean effective number of founder genomes: %s' % (summary_stats['n_g'])
        print '*'*100

    # Write summary statistics to the output file at the end of the last round.
    dout = open(outputfile,'a')
    dout.write('%s\n'%(myline2))
    dout.write('Summary statistics from %s-round gene-drop simulation\n' % (rounds))
    dout.write('\tNumber of distinct founder alleles: %s\n' % (2*nfounders))
    dout.write('\tMean allele count in latest generation: %s\n' % (summary_stats['allele_count']))
    dout.write('\tMean number of distinct alleles in latest generation: %s\n' % (summary_stats['distinct_allele_count']))
    dout.write('\tFrequency of distinct alleles sampled:\n')
    for k in summary_stats['distinct_alleles'].keys():
        dout.write('\t\tAllele %s:\t%s (%s)\n' % (k,summary_stats['distinct_alleles'][k],summary_stats['distinct_alleles'][k]**2))
    dout.write('\tMean effective number of founder genomes: %s\n'%(summary_stats['n_g']))
    dout.close()

    try: logging.info('Exited effective_founder_genomes()')
    except: pass
    return summary_stats['n_g']

##
# generation_intervals() computes the average age of parents at the time of
# birth of their first (oldest) offspring.  This is implies that selection
# decisions are made at the time of birth of the first offspring.  Average
# ages are computed for each of four paths: sire-son, sire-daughter, dam-son,
# and dam-daughter.  An overall mean is computed, as well.  IT IS IMPORTANT
# to note that if you DO NOT provide birthyears in your pedigree file that the
# returned dictionary will contain only zeroes!  This is because when no birthyear
# is provided a default value (1900) is assigned to all animals in the pedigree.
# @param pedobj A PyPedal pedigree object.
# @param units A character indicating the units in which the generation lengths should be returned.
# @retval A dictionary containing the five average ages.
def generation_intervals(pedobj, units='y'):
    """
    generation_intervals() computes the average age of parents at the time of
    birth of their first (oldest) offspring.  This is implies that selection
    decisions are made at the time of birth of the first offspring.  Average
    ages are computed for each of four paths: sire-son, sire-daughter, dam-son,
    and dam-daughter.  An overall mean is computed, as well.  IT IS IMPORTANT
    to note that if you DO NOT provide birthyears in your pedigree file that the
    returned dictionary will contain only zeroes!  This is because when no birthyear
    is provided a default value (1900) is assigned to all animals in the pedigree.
    """
    try: logging.info('Entered generation_intervals()')
    except: pass
    _sire_son = {}
    _sire_dau = {}
    _dam_son = {}
    _dam_dau = {}
    _n_unks = 0

    if not pedobj.kw['set_offspring']:
        pyp_utils.assign_offspring(pedobj)
        pedobj.kw['set_offspring'] = 1

    for m in pedobj.pedigree:
        #print m.sex
        if pedobj.kw['debug_messages']:
            if len(m.sons) > 0 or len(m.daus) > 0:
                print 'Animal %s has sex %s' % (m.animalID,m.sex)
        if m.sex == 'u' or m.sex == 'U':
            _n_unks = _n_unks + 1
        if pedobj.kw['debug_messages']:
            print '\tAnimal: %s (%s)' % (m.animalID, m.originalID)
            print '\t\tsons: %s' % (m.sons)
            print '\t\tdaus: %s' % (m.daus)
        _by = m.by
        _oldestson = -999
        _oldestdau = -999
        #
        # Walk through the sons list for this animal and use the birthyear attribute
        # to assign the oldest son to _oldestson.  If more than one son has the same
        # age then keep the first-seen son as the oldest.  It does not matter which we
        # keep for the computation; keeping the first is simply the easiest thing to
        # do.
        #
        if len(m.sons) > 0:
            if pedobj.kw['debug_messages']:
                print '\tAnimal %s sons: %s' % (m.animalID,m.sons)
            for s in m.sons:
                s = int(s)
                if _oldestson == -999:
                    _oldestson = int(s)
                elif int(pedobj.pedigree[s-1].by)  < int(pedobj.pedigree[_oldestson-1].by):
                    _oldestson = int(s)
                else:
                    pass
        if pedobj.kw['debug_messages']:
            print '\t\t_oldestson: %s' % (_oldestson)
        #
        # Walk through the daus list for this animal and use the birthyear attribute
        # to assign the oldest dau to _oldestdau.  If more than one dau has the same
        # age then keep the first-seen dau as the oldest.  It does not matter which we
        # keep for the computation; keeping the first is simply the easiest thing to
        # do.
        #
        if len(m.daus) > 0:
            if pedobj.kw['debug_messages']:
                print '\tAnimal %s daus: %s' % (m.animalID,m.daus)
            for d in m.daus:
                d = int(d)
                if _oldestdau == -999:
                    _oldestdau = int(d)
                elif pedobj.pedigree[d-1].by  < pedobj.pedigree[_oldestdau-1].by:
                    _oldestdau = int(d)
                else:
                    pass
        if pedobj.kw['debug_messages']:
            print '\t\t_oldestdau: %s' % (_oldestdau)
        #
        # This is where we assign sons and daus to one of the four selection paths:
        # sire-son, dam-son, sire-daughter, and dam-daughter.
        #
        # If the animal is a sire and has offspring, he contributes paths to the
        # sire-son and sire-daughter dictionaries.
        if m.sex == 'm':
            if _oldestson != -999:
                if pedobj.kw['debug_messages']:
                    print '\tAdding sire-son pair %s-%s to _sire_son[]' % (m.animalID,_oldestson)
                _sire_son[m.animalID] = _oldestson
            if _oldestdau != -999:
                if pedobj.kw['debug_messages']:
                    print '\tAdding sire-dau pair %s-%s to _sire_dau[]' % (m.animalID,_oldestdau)
                _sire_dau[m.animalID] = _oldestdau
        # If the animal is a dam and has offspring, she contributes paths to the
        # dam-son and dam-daughter dictionaries.
        elif m.sex == 'f':
            if _oldestson != -999:
                if pedobj.kw['debug_messages']:
                    print '\tAdding dam-son pair %s-%s to _dam_son[]' % (m.animalID,_oldestson)
                _dam_son[m.animalID] = _oldestson
            if _oldestdau != -999:
                if pedobj.kw['debug_messages']:
                    print '\tAdding dam-dau pair %s-%s to _dam_dau[]' % (m.animalID,_oldestdau)
                _dam_dau[m.animalID] = _oldestdau
        else:
            pass

    if pedobj.kw['messages'] == 'verbose':
        if _n_unks > 0:
            print '\t[MESSAGE]: %s of %s animals in the pedigree were of unknown sex and were excluded from calculations.' \
                % (_n_unks,len(pedobj.pedigree))
        print '\tPaths:'
        print '\t\tSire-Son: %s' % (_sire_son)
        print '\t\tSire-Dau: %s' % (_sire_dau)
        print '\t\tDam-Son: %s' % (_dam_son)
        print '\t\tDam-Dau: %s' % (_dam_dau)

    #
    # Now that we have four dictionaries with the parent-offspring pathways
    # corresponding to the births of oldest offspring we need to compute the
    # actual generation lengths.  For now, we are going to compute them in years.
    # If you are a mouse or fly person, you may not get the answers that you are
    # expecting here.
    #
    _ssy = 0.
    _sdy = 0.
    _dsy = 0.
    _ddy = 0.
    _overall = 0.
    try:
        # For each path in the sire-son dictionary compute the difference between
        # the son's year of birth and the sire's year of birth.  Add that difference
        # to an accumulator.
        for k,v in _sire_son.iteritems():
            #print 'Son %s\'s birthyear: %s' % ( v, int(pedobj.pedigree[int(v)-1].by) )
            #print 'Sire %s\'s birthyear: %s' % ( k, int(pedobj.pedigree[int(k)-1].by) )
            _ssy = _ssy + ( int(pedobj.pedigree[int(v)-1].by) - int(pedobj.pedigree[int(k)-1].by) )
        _ssym = _ssy / len(_sire_son)
    except:
        _ssym = 0.
    try:
        for k,v in _sire_dau.iteritems():
            _sdy = _sdy + ( int(pedobj.pedigree[int(v)-1].by) - int(pedobj.pedigree[int(k)-1].by) )
        _sdym = _sdy / len(_sire_dau)
    except:
        _sdym = 0.
    try:
        for k,v in _dam_son.iteritems():
            _dsy = _dsy + ( int(pedobj.pedigree[int(v)-1].by) - int(pedobj.pedigree[int(k)-1].by) )
        _dsym = _dsy / len(_dam_son)
    except:
        _dsym = 0.
    try:
        for k,v in _dam_dau.iteritems():
            _ddy = _ddy + ( int(pedobj.pedigree[int(v)-1].by) - int(pedobj.pedigree[int(k)-1].by) )
        _ddym = _ddy / len(_dam_dau)
    except:
        _ddym = 0.
    try:
        _overall = ( _ssym + _sdym + _dsym + _ddym ) / 4.
    except:
        _overall = 0.

    _genlens = {}
    _genlens['ss'] = _ssym
    _genlens['sd'] = _sdym
    _genlens['ds'] = _dsym
    _genlens['dd'] = _ddym
    _genlens['mean'] = _overall

    if pedobj.kw['messages'] == 'verbose':
        print '\tMeans:'
        print '\t\tSire-Son: %s' % (_ssym)
        print '\t\tSire-Dau: %s' % (_sdym)
        print '\t\tDam-Son: %s' % (_dsym)
        print '\t\tDam-Dau: %s' % (_ddym)
        print '\t\tOverall: %s' % (_overall)

    try: logging.info('Exited generation_intervals()')
    except: pass
    return _genlens

##
# generation_intervals_all() computes the average age of parents at the time of
# birth of their offspring.  The computation is made using birth years for all
# known offspring of sires and dams, which implies discrete generations.  Average
# ages are computed for each of four paths: sire-son, sire-daughter, dam-son, and
# dam-daughter.  An overall mean is computed, as well. IT IS IMPORTANT to note that
# if you DO NOT provide birthyears in your pedigree file that the returned dictionary
# will contain only zeroes!  This is because when no birthyear is provided a default
# value (1900) is assigned to all animals in the pedigree.
# @param pedobj A PyPedal pedigree object.
# @param units A character indicating the units in which the generation lengths should be returned.
# @retval A dictionary containing the five average ages.
def generation_intervals_all(pedobj, units='y'):
    """
    generation_intervals_all() computes the average age of parents at the time of
    birth of their offspring.  The computation is made using birth years for all
    known offspring of sires and dams, which implies discrete generations.  Average
    ages are computed for each of four paths: sire-son, sire-daughter, dam-son, and
    dam-daughter.  An overall mean is computed, as well. IT IS IMPORTANT to note that
    if you DO NOT provide birthyears in your pedigree file that the returned dictionary
    will contain only zeroes!  This is because when no birthyear is provided a default
    value (1900) is assigned to all animals in the pedigree.
    """
    try: logging.info('Entered generation_intervals_all()')
    except: pass
    _sire_son = {}
    _sire_dau = {}
    _dam_son = {}
    _dam_dau = {}
    _n_unks = 0

    if not pedobj.kw['set_offspring']:
        pyp_utils.assign_offspring(pedobj)
        pedobj.kw['set_offspring'] = 1

    for m in pedobj.pedigree:
        if pedobj.kw['debug_messages']:
            m.printme()
            if len(m.sons) > 0 or len(m.daus) > 0:
                if pedobj.kw['debug_messages']:
                    print 'Animal %s has sex %s' % (m.animalID,m.sex)
        if m.sex == 'u' or m.sex == 'U':
            _n_unks = _n_unks + 1
        if pedobj.kw['debug_messages']:
            print '\tAnimal: %s' % (m.animalID)
            print '\t\tsons: %s' % (m.sons)
            print '\t\tdaus: %s' % (m.daus)
        #
        # Add dictionary entries for all dam-offspring pairs.
        #
        if m.sex == 'm':
            if len(m.sons) > 0:
                for s in m.sons:
                    s = int(s)
                    #if pedobj.kw['debug_messages']:
                    print '\tAdding sire-son pair %s-%s to _sire_son' % (m.animalID,s)
                    _sire_son[m.animalID] = s
            if len(m.daus) > 0:
                for d in m.daus:
                    d = int(d)
                    #if pedobj.kw['debug_messages']:
                    print '\tAdding sire-dau pair %s-%s to _sire_dau' % (m.animalID,d)
                    _sire_dau[m.animalID] = d
        #
        # Add dictionary entries for all dam-offspring pairs.
        #
        if m.sex == 'f':
            if len(m.sons) > 0:
                for s in m.sons:
                    s = int(s)
                    #if pedobj.kw['debug_messages']:
                    print '\tAdding sire-son pair %s-%s to _dam_son' % (m.animalID,s)
                    _dam_son[m.animalID] = s
            if len(m.daus) > 0:
                for d in m.daus:
                    d = int(d)
                    #if pedobj.kw['debug_messages']:
                    print '\tAdding sire-dau pair %s-%s to _dam_dau' % (m.animalID,d)
                    _dam_dau[m.animalID] = d

    if pedobj.kw['messages'] == 'verbose':
        if _n_unks > 0:
            print '\t[MESSAGE]: %s of %s animals in the pedigree were of unknown sex and were excluded from calculations.' \
                % (_n_unks,len(pedobj.pedigree))
        print '\tPaths:'
        print '\t\tSire-Son: %s' % (_sire_son)
        print '\t\tSire-Dau: %s' % (_sire_dau)
        print '\t\tDam-Son: %s' % (_dam_son)
        print '\t\tDam-Dau: %s' % (_dam_dau)

    #
    # Now that we have four dictionaries with the parent-offspring pathways corresponding to the
    # births of offspring we need to compute the actual generation lengths.  For now, we are going
    # to compute them in years.  If you are a mouse or fly person, you may not get he answers that
    # you are expecting.
    #
    _ssy = 0.
    _sdy = 0.
    _dsy = 0.
    _ddy = 0.
    _overall = 0.
    try:
        for k,v in _sire_son.iteritems():
            _ssy = _ssy + ( int(pedobj.pedigree[int(v)-1].by) - int(pedobj.pedigree[int(k)-1].by) )
        _ssym = _ssy / len(_sire_son)
    except:
        _ssym = 0.
    try:
        for k,v in _sire_dau.iteritems():
            _sdy = _sdy + ( int(pedobj.pedigree[int(v)-1].by) - int(pedobj.pedigree[int(k)-1].by) )
        _sdym = _sdy / len(_sire_dau)
    except:
        _sdym = 0.
    try:
        for k,v in _dam_son.iteritems():
            _dsy = _dsy + ( int(pedobj.pedigree[int(v)-1].by) - int(pedobj.pedigree[int(k)-1].by) )
        _dsym = _dsy / len(_dam_son)
    except:
        _dsym = 0.
    try:
        for k,v in _dam_dau.iteritems():
            _ddy = _ddy + ( int(pedobj.pedigree[int(v)-1].by) - int(pedobj.pedigree[int(k)-1].by) )
        _ddym = _ddy / len(_dam_dau)
    except:
        _ddym = 0.
    try:
        _overall = ( _ssym + _sdym + _dsym + _ddym ) / 4.
    except:
        _overall = 0.

    _genlens = {}
    _genlens['ss'] = _ssym
    _genlens['sd'] = _sdym
    _genlens['ds'] = _dsym
    _genlens['dd'] = _ddym
    _genlens['mean'] = _overall

    if pedobj.kw['messages'] == 'verbose':
        print '\tMeans:'
        print '\t\tSire-Son: %s' % (_ssym)
        print '\t\tSire-Dau: %s' % (_sdym)
        print '\t\tDam-Son: %s' % (_dsym)
        print '\t\tDam-Dau: %s' % (_ddym)
        print '\t\tOverall: %s' % (_overall)

    try: logging.info('Exited generation_intervals_all()')
    except: pass
    return _genlens

##
# founder_descendants() returns a dictionary containing a list of descendants of
# each founder in the pedigree.
# @param pedobj An instance of a PyPedal NewPedigree object.
# @retval A dictionary containing a list of descendants of each founder in the pedigree.
def founder_descendants(pedobj):
    """
    founder_descendants() returns a dictionary containing a list of descendants of
    each founder in the pedigree.
    """
    try: logging.info('Entered founder_descendants()')
    except: pass
    founder_peds = {}
    for f in pedobj.metadata.unique_founder_list:
        _desc = descendants(pedobj.idmap[f],pedobj,{})
        founder_peds[f] = _desc
    try: logging.info('Exited founder_descendants()')
    except: pass
    return founder_peds

##
# descendants() uses pedigree metadata to walk a pedigree and return a list of all
# of the descendants of a given animal.
# @param anid An animal ID
# @param pedobj A Python list of PyPedal Animal() objects.
# @param _desc A Python dictionary of descendants of animal anid.
# @retval A list of descendants of anid.
def descendants(anid, pedobj, _desc):
    """
    descendants() uses pedigree metadata to walk a pedigree and return a list of all
    of the descendants of a given animal.
    """
    try: logging.info('Entered descendants()')
    except: pass
    #pedobj.pedigree[int(anid)-1].printme()
    #print 'unks: ', pedobj.pedigree[int(anid)-1].unks
    #print 'sons: ', pedobj.pedigree[int(anid)-1].sons
    #print 'daus: ', pedobj.pedigree[int(anid)-1].daus

    if len(pedobj.pedigree[int(anid)-1].unks) > 0:
        for _u,_v in pedobj.pedigree[int(anid)-1].unks.iteritems():
            try:
                _utest = _desc[_u]
            except KeyError:
                _desc[_u] = _u
            _desc = descendants(_u,pedobj,_desc)
    if len(pedobj.pedigree[int(anid)-1].sons) > 0:
        for _u,_v in pedobj.pedigree[int(anid)-1].sons.iteritems():
            try:
                _utest = _desc[_u]
            except KeyError:
                _desc[_u] = _u
            _desc = descendants(_u,pedobj,_desc)
    if len(pedobj.pedigree[int(anid)-1].daus) > 0:
        for _u,_v in pedobj.pedigree[int(anid)-1].daus.iteritems():
            try:
                _utest = _desc[_u]
            except KeyError:
                _desc[_u] = _u
            _desc = descendants(_u,pedobj,_desc)
    try: logging.info('Exited descendants()')
    except: pass
    return _desc

##
# dropped_ancestral_inbreeding() uses a gene dropping approach to calculate
# ancestral inbreeding, the probability of an individual inheriting an allele
# that has undergone inbreeding in the past at least once.
# @param pedobj A PyPedal pedigree object.
# @param rounds The number of times to simulate segregation through the entire pedigree.
# @param loci The number biallelic, unlinked loci to simulate.
# @param frequency The minor allele frequency.
# @param seed The seed for the RNG.
# @retval A dictionary of ancestral inbreeding coefficients keyed to animal IDs.
def dropped_ancestral_inbreeding(pedobj, rounds=100, loci=100, frequency=0.05, seed=5048665):
    """
    dropped_ancestral_inbreeding() uses a gene dropping approach to calculate
    ancestral inbreeding, the probability of an individual inheriting an allele
    that has undergone inbreeding in the past at least once.
    """
    try: logging.info('Entered dropped_ancestral_inbreeding()')
    except: pass
    if rounds < 1:
        if pedobj.kw['debug_messages']:
            print '[ERROR]: The rounds parameter in pyp_metrics/dropped_ancestral_inbreeding() must be positive (>0)! Defaulting to 100.'
        rounds = 100
    if loci < 1:
        if pedobj.kw['debug_messages']:
            print '[ERROR]: The loci parameter in pyp_metrics/dropped_ancestral_inbreeding() must be positive (>0)! Defaulting to 100.'
        rounds = 100
    if frequency < 0.:
        if pedobj.kw['debug_messages']:
            print '[ERROR]: The frequency parameter in pyp_metrics/dropped_ancestral_inbreeding() must be positive (>0)! Defaulting to 0.01.'
        frequency = 0.01
    if frequency > 1.:
        if pedobj.kw['debug_messages']:
            print '[ERROR]: The frequency parameter in pyp_metrics/dropped_ancestral_inbreeding() cannot be >1! Defaulting to 0.01.'
        frequency = 0.01
    try:
        _seed = int(seed)
    except:
        if pedobj.kw['debug_messages']:
            print '[ERROR]: The seed parameter in pyp_metrics/dropped_ancestral_inbreeding() must be an integer! Defaulting to 5048665.'
        seed = 5048665
    # Build a list of generations
    l = len(pedobj.pedigree)
    # Build a list of unique generations so that we can find the most recent generation
    gens = []       # list of generation IDs in the pedigree
    ngen = 0    # Number of animals in the latest generation
    nfounders = 0
    n_g = 0.0   # Effective number of founder genomes
    for j in range(l):
        g = pedobj.pedigree[j].gen
        if g in gens:
            pass
        else:
            gens.append(g)
        if pedobj.pedigree[j].founder == 'y':
            nfounders = nfounders + 1
    gens.sort()
    allele_dict = {}
    summary_freqs = {}
    id2aic = {}     # Map animal IDs to ancestral ID coefficients
    summary_freqs['distinct_alleles'] = {}
    summary_stats = {}
    summary_stats['distinct_alleles'] = {}
    outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_gene_drop','.out')
    myline = '='*80
    myline2 = '*'*80
    numpy.random.seed(seed)
    for r in range(rounds):
        allele_inbred = {}
        # The dictionary ancestor_alleles contains alleles for
        # eah locus in the simulaton.
        for p in pedobj.pedigree:
            p.ancestor_alleles = {}
            # Create a new set of alleles to drop for each round.
            for lo in xrange(loci):
                p.ancestor_alleles[lo] = []
                # Founders contribute two novel alleles
                if p.founder == 'y':
                    _allele_1 = '%s_%s_1' % (p.paddedID,lo+1)
                    _allele_2 = '%s_%s_2' % (p.paddedID,lo+1)
                # Half-founders contribute one novel allele
                elif str(p.sireID) == str(pedobj.kw['missing_parent']):
                    _allele_1 = '%s_%s_1' % (p.paddedID,lo+1)
                    _allele_2 = ''
                elif str(p.damID) == str(pedobj.kw['missing_parent']):
                    _allele_1 = ''
                    _allele_2 = '%s_%s_1' % (p.paddedID,lo+1)
                else:
                    _allele_1 = ''
                    _allele_2 = ''
                p.ancestor_alleles[lo].append(_allele_1)
                p.ancestor_alleles[lo].append(_allele_2)
                allele_dict[lo] = {}
            allele_inbred[p.animalID] = numpy.zeros([loci],'d')
            id2aic[p.animalID] = 0.
        # Now that we have the ancestor loci populated we can start
        # the gene dropping
        for lo in xrange(loci):
            for i in range(l):
                if pedobj.pedigree[i].founder == 'y':
                    pass
                else:
                    # Pick a sire allele at random
                    if numpy.random.ranf() < frequency:
                        _as = pedobj.pedigree[int(pedobj.pedigree[i].sireID)-1].ancestor_alleles[lo][0]
                    else:
                        _as = pedobj.pedigree[int(pedobj.pedigree[i].sireID)-1].ancestor_alleles[lo][1]
                    # Pick a dam allele at random
                    if numpy.random.ranf() < frequency:
                        _ad = pedobj.pedigree[int(pedobj.pedigree[i].damID)-1].ancestor_alleles[lo][0]
                    else:
                        _ad = pedobj.pedigree[int(pedobj.pedigree[i].damID)-1].ancestor_alleles[lo][1]
                    # This looks funny, right? The whole idea here is that we're looking for
                    # inbreeding which occured BEFORE a pair of alleles arrived in this animal,
                    # so we've got to look back and see if the alleles were IBD in the parent.
                    pedobj.pedigree[i].ancestor_alleles[lo] = [_as,_ad]
                    if pedobj.pedigree[int(pedobj.pedigree[i].sireID)-1].ancestor_alleles[lo][0] == \
                        pedobj.pedigree[int(pedobj.pedigree[i].sireID)-1].ancestor_alleles[lo][1]:
                        allele_inbred[pedobj.pedigree[i].animalID][lo] = 1.
                    if pedobj.pedigree[int(pedobj.pedigree[i].damID)-1].ancestor_alleles[lo][0] == \
                        pedobj.pedigree[int(pedobj.pedigree[i].damID)-1].ancestor_alleles[lo][1]:
                        allele_inbred[pedobj.pedigree[i].animalID][lo] = 1.
        # Sumamarize allele data
        for p in pedobj.pedigree:
            id2aic[p.animalID] = id2aic[p.animalID] + allele_inbred[p.animalID].mean()
    try: logging.info('Exited dropped_ancestral_inbreeding()')
    except: pass
    return id2aic

##
# ballou_ancestral_inbreeding() calculates ancestral inbreeding,
# the probability of an individual inheriting an allele that has
# undergone inbreeding in the past at least once, using the method
# of Ballou (1997).
# @param pedobj A PyPedal pedigree object.
# @retval A dictionary of ancestral inbreeding coefficients keyed to animal IDs.
def ballou_ancestral_inbreeding(pedobj):
    """
    ballou_ancestral_inbreeding() calculates ancestral inbreeding,
    the probability of an individual inheriting an allele that has
    undergone inbreeding in the past at least once, using the method
    of Ballou (1997).
    """
    try: logging.info('Entered ballou_ancestral_inbreeding()')
    except: pass
    # Initialize the dictionary mapping animal IDs to ancestor inbreeding.
    id2aic = {}
    for p in pedobj.pedigree:
        id2aic[p.animalID] = 0.
    # Calculate coefficients of inbreeding if they're not already in the pedigree.
    if not pedobj.kw['f_computed']:
        pyp_nrm.inbreeding(pedobj)
    # Calculate ancestral inbreeding
    for p in pedobj.pedigree:
        if str(p.sireID) == str(pedobj.kw['missing_parent']):
            f_s = 0.
            f_as = 0.
        else:
            f_s = pedobj.pedigree[p.sireID-1].fa
            f_as = id2aic[p.sireID]
        if str(p.damID) == str(pedobj.kw['missing_parent']):
            f_d = 0.
            f_ad = 0.
        else:
            f_d = pedobj.pedigree[p.damID-1].fa
            f_ad = id2aic[p.damID]
        id2aic[p.animalID] = ( f_as + (1.-f_as)*f_s + \
            f_ad + (1.-f_ad) * f_d ) / 2.
    try: logging.info('Exited ballou_ancestral_inbreeding()')
    except: pass
    return id2aic
