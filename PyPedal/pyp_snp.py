#!/usr/bin/python

###############################################################################
# NAME: pyp_snp.py
# VERSION: 2.0.0 (11MAY2020)
# AUTHOR: John B. Cole, PhD (john.cole@usda.gov)
# LICENSE: LGPL
###############################################################################
# FUNCTIONS:
#   read_agil_chromosome_data()
#   read_agil_pedigree_file()
#   read_agil_genotypes_txt()
#   read_agil_true_frequency()
#   form_p_matrix_from_snp()
#   form_m_matrix_from_snp()
#   form_grm_from_snp()
#   compute_genomic_inbreeding_from_grm()
#   compute_genomic_homozygosity_from_snp()
#   renumber_snp_ids()
#   generate_random_genotype()
###############################################################################

## @package pyp_snp
# pyp_snp contains several procedures for working with single nucleotide poly-
# morphism (SNP) genotype data.
##

import logging
import numpy as np
import pandas as pd
import random

##
# read_agil_chromosome_data() loads SNP marker information from the chromosome.data file
# used by AGIL and CDCB. Note that this ONLY reads the first 5 columns (SNP name, chromosome
# number, within-chromosome marker number, overall marker number, and location in base pairs).
# @param filename The name of the input file.
# @retval A Pandas dataframe, or False if the file could not be read.
def read_agil_chromosome_data(filename='chromosome.data'):
    """
    read_agil_chromosome_data() loads SNP marker information from the chromosome.data file
    used by AGIL and CDCB. Note that this ONLY reads the first 5 columns (SNP name, chromosome
    number, within-chromosome marker number, overall marker number, and location in base pairs).
    """
    try:
        df = pd.read_csv(filename, header=0, skiprows=0, delim_whitespace=True,
                         usecols=[0, 1, 2, 3, 4], colnames=['snp_name', 'chrome', 'within', 'overall',
                                                            'location'])
        print '[INFO]: pyp_snp/read_agil_chromosome_data() read information about %s SNP from file %s.' % \
              (len(df), filename)
        logging.info('pyp_snp/read_agil_chromosome_data() read information about %s SNP from file %s.' % \
                     (len(df), filename))
        return df
    except:
        logging.error('pyp_snp/read_agil_chromosome_data() could not access the file %s.' % \
                      filename)
        print '[ERROR]: pyp_snp/read_agil_chromosome_data() could not access the file %s.' % \
              filename
        return False

##
# read_agil_pedigree_file() loads pedigree information from the file format used by AGIL and CDCB.
# @param filename The name of the input file.
# @retval A Pandas dataframe, or False if the file could not be read.
def read_agil_pedigree_file():
    pass

##
# read_agil_genotypes_txt() loads SNP genotypes from the genotypes.txt file used by AGIL and CDCB.
# @param filename The name of the input file.
# @retval A Pandas dataframe, or False if the file could not be read.
def read_agil_genotypes_txt(filename='genotypes.txt'):
    try:
        df = pd.read_csv(filename, delim_whitespace=True, dtype={'genotype': str},
                         names=['animalID', 'chip_type', 'n_snps', 'genotype'] )
        print '[INFO]: pyp_snp/read_agil_genotypes_txt() read information about %s SNP genotypes ' \
              'from file %s.' % (len(df), filename)
        logging.info('pyp_snp/read_agil_genotypes_txt() read information about %s SNP genotypes '
                     'from file %s.' % (len(df), filename))
        return df
    except:
        logging.error('pyp_snp/read_agil_genotypes_txt() could not access the file %s.' % \
                      filename)
        print '[ERROR]: pyp_snp/read_agil_genotypes_txt() could not access the file %s.' % \
              filename
        return False

##
# read_agil_true_frequency() loads SNP frequency data from the true.frequency file used by AGIL and CDCB.
# @param filename The name of the input file.
# @retval A Pandas dataframe, or False if the file could not be read.
def read_agil_true_frequency(filename='true.frequency'):
    """
    (base) ARSNYPLM53200LF:Example_Output john.cole$ head true.frequency
    Marker000000001 1   0.88990
    Marker000000002 2   0.32821
    """
    try:
        df = pd.read_csv('true.frequency', delim_whitespace=True,
                         names=['snp_name', 'overall', 'frequency'])
        print '[INFO]: pyp_snp/read_agil_true_frequency() read information about %s SNP ' \
              'from file %s.' % (len(df), filename)
        logging.info('pyp_snp/read_agil_true_frequency() read information about %s SNP '
                     'from file %s.' % (len(df), filename))
        return df
    except:
        logging.error('pyp_snp/read_agil_true_frequency() could not access the file %s.' % \
                      filename)
        print '[ERROR]: pyp_snp/read_agil_true_frequency() could not access the file %s.' % \
              filename
        return False

##
# form_p_matrix_from_snp() calculates individual SNP frequencies from the genotypes provided.
# @param pedobj A PyPedal pedigree object.
# @param debug Turns debug messages on (True) and off (False).
# @retval A NumPy array, or False if the file could not be read.
def form_p_matrix_from_snp(pedobj, debug=False):
    """
    :param pedobj:
    :param debug:
    :return:
    """

    if pedobj.snp.empty:

        logging.error('pyp_snp/form_p_matrix_from_snp() could not form P because the pedigree has no SNP data.')
        print '[ERROR]: pyp_snp/form_p_matrix_from_snp() could not form P because the pedigree has no SNP data.'
        return False

    else:

        # Form P, an n-by-m array of allele frequencies, as in VanRaden (2008; p. 4416).
        # We can estimate the frequencies from the data provided in the genotypes file.
        # f_i = \sum_j^n m_ij / 2*n where n is the number of genotyped animals
        P = np.zeros( [ len(pedobj.snp.iloc[0,3]) ] )
        for a in xrange(len( pedobj.snp )):
            for s in xrange( len(pedobj.snp.iloc[0,3]) ):
                P[s] = P[s] + float(pedobj.snp.iloc[a, 3][s])
        P = P / ( 2. * len( pedobj.snp ) )

        if debug:
            print
            print 'P:\t', P
            print

        return P

##
# form_m_matrix_from_snp() form the matrix, M, that specifies which marker alleles each individual inherited.
# Dimensions of M are the number of individuals (n) by the number of loci (m).
# @param pedobj A PyPedal pedigree object.
# @param scale_m Scale M to {-1, 0, +1}.
# @param debug Turns debug messages on (True) and off (False).
# @retval A NumPy array, or False if the file could not be read.
def form_m_matrix_from_snp(pedobj, scale_m = True, debug=False):
    """
    :param pedobj:
    :param scale_m:
    :param debug:
    :return:
    """

    if pedobj.snp.empty:

        logging.error('pyp_snp/form_m_matrix_from_snp() could not form M because the pedigree has no SNP data.')
        print '[ERROR]: pyp_snp/form_m_matrix_from_snp() could not form M because the pedigree has no SNP data.'
        return False

    else:

        P = form_p_matrix_from_snp(pedobj, debug=debug)

        # Form M, an n-by-m matrix of the m alleles inherited by each of the n individuals
        # in the population, as in VanRaden (2008; p. 4416).
        M = np.zeros( [ len(pedobj.snp), len(pedobj.snp.iloc[0,3]) ] )
        for a in range( len( pedobj.snp ) ):
            for s in range( len( pedobj.snp.iloc[0,3] ) ):
                M[a, s] = pedobj.snp.iloc[a, 3][s]
                if scale_m:
                    # Scale the elements of M to -1, 0, and 1 for the homozygote, heterozygote,
                    # and other homozygote, respectively
                    M[a, s] = M[a, s] - P[s]

        if debug:
            print
            print 'M:\t', M
            print

        return M

##
# form_grm_from_snp() forms the genomic relationship matrix, G, from the SNP information provided. G is a square
# matrix with an order of the number of individuals (n).
# @param pedobj A PyPedal pedigree object.
# @param scale_m Scale M to {-1, 0, +1}.
# @param method Which of VanRaden's (2008) methods to use when forming G (1|2|3).
# @param debug Turns debug messages on (True) and off (False).
# @retval A NumPy array, or False if the file could not be read.
def form_grm_from_snp(pedobj, scale_m=True, method=1, debug=False):
    """
    :param pedobj:
    :param scale_m:
    :param method:
    :param debug:
    :return:
    """

    if pedobj.snp.empty:

        logging.error('pyp_snp/form_grm_from_snp() could not form a GRM because the pedigree has no SNP data.')
        print '[ERROR]: pyp_snp/form_grm_from_snp() could not form a GRM because the pedigree has no SNP data.'
        return False

    else:

        P = form_p_matrix_from_snp(pedobj, debug=debug)

        M = form_m_matrix_from_snp(pedobj, scale_m=scale_m, debug=debug)

        # Initialize the genomic relationship matrix
        G = np.zeros( [ len(pedobj.snp), len(pedobj.snp) ] )

        # Form G using VanRaden's Method 1
        if method == 1:

            # Compute denominator: 2 * \sump_i(1-p_1)
            sum_freq = 0.
            for i in range( len( P ) ):
                sum_freq += P[i] * ( 1. - P[i] )

            if debug:
                print
                print 'sum_freq:\t', sum_freq
                print '2*sum_freq:\t', 2*sum_freq
                print

            Z = M - P

            if debug:
                print
                print 'Z:\t', Z
                print

            G = Z.dot(Z.T) / ( 2. * sum_freq )

            if debug:
                print
                print 'G:\t', G
                print

        else:

            logging.error('pyp_snp/form_grm_from_snp() could not form a GRM using VanRaden method %s because that'
                          'has not yet been implemented.', method)
            print '[ERROR]: pyp_snp/form_grm_from_snp() could not form a GRM using VanRaden method %s because that' \
                  'has not yet been implemented.' % method
            return False

##
# compute_genomic_inbreeding_from_grm() calculates genomic inbreeding from the diagonals of G, and computed summary
# statistics for genomic coefficients of inbreeding and coefficients of relationship.
# @param pedobj A PyPedal pedigree object.
# @param g_matrix A G matrix precomputed from the SNP genotypes in the attached pedigree object.
# @param scale_m Scale M to {-1, 0, +1}.
# @param rels Calculate summary statistics for coefficients of genomic inbreeding (True|False).
# @param method Which of VanRaden's (2008) methods to use when forming G (1|2|3).
# @param debug Turns debug messages on (True) and off (False).
# @retval A NumPy array, or False if the file could not be read.
def compute_genomic_inbreeding_from_grm(pedobj, g_matrix=False, scale_m=True, rels=False, update_pedigree=True,
                                        output=True, debug=False):
    """
    :param pedobj:
    :param g_matrix:
    :param rels:
    :param update_pedigree:
    :param output:
    :param debug:
    :return:
    """

    if not g_matrix:

        if pedobj.snp.empty:
            logging.error('pyp_snp/compute_genomic_inbreeding_from_grm() could not form a GRM because the pedigree has no SNP data.')
            print '[ERROR]: pyp_snp/compute_genomic_inbreeding_from_grm() could not form a GRM because the pedigree has no SNP data.'
            return False

        else:
            G = form_grm_from_snp(pedobj, scale_m=scale_m, method=1, debug=False)

    if len(pedobj.snp.index) != len(pedobj.pedigree):
        logging.warning('pyp_snp/compute_genomic_inbreeding_from_grm(): There are different numbers of SNP '
                        'genotypes and animals in the pedigree file, which can lead to errors in matching '
                        'genomic relationships and coefficients of inbreeding!')
        print '[WARNING]: pyp_snp/compute_genomic_inbreeding_from_grm(): There are different numbers of SNP ' \
              'genotypes and animals in the pedigree file, which can lead to errors in matching genomic ' \
              'relationships and coefficients of inbreeding!'

    # Setup data structures
    fx = {}
    metadata = {}

    if rels:
        rel_dict = {}
        rel_dict['r_count'] = (pedobj.snp.index * (pedobj.snp.index + 1)) / 2
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

    # Pull inbreeding coefficients out of the genomic relationship matrix.
    for _i in xrange(pedobj.snp.index):
        fx[pedobj.pedigree[_i].animalID] = G[_i, _i] - 1.
        if update_pedigree:
            pedobj.pedigree[_i].fg = G[_i, _i] - 1.

    # Pull coefficients of relationship  out of the genomic relationship matrix.
    if rels == 1:
        n = len(pedobj.snp.index)
        reldict['r_count'] = ( n * ( n + 1 ) ) / 2
        for i in xrange(n):
            for j in xrange (i, n):
                if i != j:
                    if G[i, j] > 0.:
                        reldict['r_nonzero_count'] = \
                            reldict['r_nonzero_count'] + 1
                        reldict['r_nonzero_sum'] = reldict['r_nonzero_sum'] + G[i, j]
                        if pedobj.nrm.nrm[i][j] > reldict['r_max']:
                            reldict['r_max'] = G[i, j]
                        if pedobj.nrm.nrm[i][j] < reldict['r_min']:
                            reldict['r_min'] = G[i, j]
                    reldict['r_sum'] = reldict['r_sum'] + G[i, j]

    # Write summary statistics to a file.
    if output:
        a_outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_genomic_inbreeding','.dat')
        aout = open(a_outputfile, 'w')
        aout.write('# Genomic inbreeding coefficients\n')
        # If the pedigree uses names do the same in the output file because the original IDs won't mean
        # anything when compared to the input pedigree.
        if 'ASD' in pedobj.kw['pedformat']:
            aout.write('# Name\tRenum ID\tf_g\n')
        else:
            aout.write('# Orig ID\tRenum ID\tf_g\n')

    # Compute summary statistics for genomic inbreeding
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
        # Update self.fg for each Animal object in the pedigree.
        if update_pedigree:
            pedobj.pedigree[int(k)-1].fg = v
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
    # If there are no inbred animals in the pedigree we need to make sure that meaningful values are returned, rather
    # than, e.g., -999.
    if f_nonzero_count == 0:
        f_nonzero_rng = 0.
        f_nonzero_avg = 0.
        f_nonzero_min = 0.
        f_nonzero_max = 0.
    else:
        f_nonzero_rng = f_nonzero_max - f_nonzero_min
        f_nonzero_avg = f_nonzero_sum / f_nonzero_count
    # Summary statistics including all genomic CoI
    metadata['all'] = {}
    metadata['all']['f_count'] = len(fx.keys())
    metadata['all']['f_sum'] = f_sum
    metadata['all']['f_min'] = f_min
    metadata['all']['f_max'] = f_max
    metadata['all']['f_rng'] = f_rng
    metadata['all']['f_avg'] = f_avg
    # Summary statistics including only nonzero genomic CoI
    metadata['nonzero'] = {}
    metadata['nonzero']['f_count'] = f_nonzero_count
    metadata['nonzero']['f_sum'] = f_nonzero_sum
    metadata['nonzero']['f_min'] = f_nonzero_min
    metadata['nonzero']['f_max'] = f_nonzero_max
    metadata['nonzero']['f_rng'] = f_nonzero_rng
    metadata['nonzero']['f_avg'] = f_nonzero_avg

    if rels:
        if pedobj.kw['debug_messages'] == 1:
            print '[DEBUG]: reldict: ', reldict
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

    if output:
        line = '='*80
        aout.write('%s\n' % line)
        aout.write('Genomic Inbreeding Statistics\n')
        line = '-'*80
        aout.write('All animals:\n')
        aout.write('%s\n' % line)
        aout.write('\tCount:\t%s\n'%len(fx.keys()))
        aout.write('\tMean:\t%s\n'%f_avg)
        aout.write('\tMin:\t%s\n'%f_min)
        aout.write('\tMax:\t%s\n'%f_max)
        line = '-'*80
        aout.write('Animals with non-zero genomic CoI:\n')
        aout.write('%s\n' % line)
        aout.write('\tCount:\t%s\n'%f_nonzero_count)
        aout.write('\tMean:\t%s\n'%f_nonzero_avg)
        aout.write('\tMin:\t%s\n'%f_nonzero_min)
        aout.write('\tMax:\t%s\n'%f_nonzero_max)
        aout.close()

        pedobj.kw['g_computed'] = 1
        out_dict = {}
        out_dict['metadata'] = metadata
        out_dict['fx'] = fx
        if rels:
            return out_dict, rel_dict
        else:
            return out_dict

##
# compute_genomic_homozygosity_from_snp() calculates genomic homozygosity for each SNP genotype as the proportion
# of homozygous loci.
# @param pedobj A PyPedal pedigree object.
# @param debug Turns debug messages on (True) and off (False).
# @retval A PyPedal pedigree object with genomic inbreeding coefficients assigned.
def compute_genomic_homozygosity_from_snp(pedobj, update_pedigree=True, output=True, debug=False):
    """
    :param pedobj:
    :param update_pedigree:
    :param debug:
    :return:
    """

    # Setup data structures
    fx = {}
    metadata = {}

    if pedobj.snp.empty:
        if pedobj.kw['debug_messages'] == 1:
            print '[ERROR]: There are no SNP data associated with this pedigree so no calculations can be performed!'
        logging.error('There are no SNP data associated with this pedigree so no calculations can be performed!')
        return False
    else:
        if len(pedobj.snp.index) != len(pedobj.pedigree):
            logging.warning('pyp_snp/compute_genomic_homozygosity_from_snp(): There are different numbers of SNP '
                            'genotypes and animals in the pedigree file, which can lead to errors in matching '
                            'genomic relationships and coefficients of inbreeding!')
            print '[WARNING]: pyp_snp/compute_genomic_homozygosity_from_snp(): There are different numbers of SNP ' \
                  'genotypes and animals in the pedigree file, which can lead to errors in matching genomic ' \
                  'relationships and coefficients of inbreeding!'

        if pedobj.kw['debug_messages'] == 1:
            print '[INFO]: pyp_snp/compute_genomic_homozygosity_from_snp(): Renumbering animal IDs in the SNP dataframe.'
        logging.info('pyp_snp/compute_genomic_homozygosity_from_snp(): Renumbering animal IDs in the SNP dataframe')
    # We count the 1s, which are the heterozygotes, so we need to subtract that from 1. to get the frequency of the
    # homozygous loci.
    for p in pedobj.pedigree:
        fx[p.animalID] = 1. - float(pedobj.snp[pedobj.snp['animalID'] == p.animalID]['genotype'].values[0].count('1')) / \
                         float(len(pedobj.snp[pedobj.snp['animalID'] == p.animalID]['genotype'].values[0]))

    # Write summary statistics to a file.
    if output:
        a_outputfile = '%s%s%s' % (pedobj.kw['filetag'],'_genomic_homozygosity','.dat')
        aout = open(a_outputfile, 'w')
        aout.write('# Genomic homozygosity\n')
        # If the pedigree uses names do the same in the output file because the original IDs won't mean
        # anything when compared to the input pedigree.
        if 'ASD' in pedobj.kw['pedformat']:
            aout.write('# Name\tRenum ID\thomozygosity\n')
        else:
            aout.write('# Orig ID\tRenum ID\thomozygosity\n')

    # Compute summary statistics for genomic inbreeding
    f_sum = 0.0
    f_min = 999.0
    f_max = -999.0
    for k, v in fx.iteritems():
        if output:
            if 'ASD' in pedobj.kw['pedformat']:
                aout.write('%s\t%s\t%s\n'%(pedobj.pedigree[int(k)-1].name, k, v))
            else:
                aout.write('%s\t%s\t%s\n'%(pedobj.pedigree[int(k)-1].originalID, k, v))
        # Update self.homozygosity for each Animal object in the pedigree.
        if update_pedigree:
            pedobj.pedigree[int(k)-1].homozygosity = v
        f_sum = f_sum + v
        if v > f_max:
            f_max = v
        if v < f_min:
            f_min = v
    if len(fx.keys()) == 0:
        f_rng = 0.
        f_avg = 0.
        f_min = 0.
        f_max = 0.
    else:
        f_avg = f_sum / len(fx.keys())
        f_rng = f_max - f_min
    # Summary statistics including all genomic CoI
    metadata = {}
    metadata['f_count'] = len(fx.keys())
    metadata['f_sum'] = f_sum
    metadata['f_min'] = f_min
    metadata['f_max'] = f_max
    metadata['f_rng'] = f_rng
    metadata['f_avg'] = f_avg
    # Prepare the output
    if output:
        line = '='*80
        aout.write('%s\n' % line)
        aout.write('Genomic Homozygosity Statistics\n')
        line = '-'*80
        aout.write('All animals:\n')
        aout.write('%s\n' % line)
        aout.write('\tCount:\t%s\n'%len(fx.keys()))
        aout.write('\tMean:\t%s\n'%f_avg)
        aout.write('\tMin:\t%s\n'%f_min)
        aout.write('\tMax:\t%s\n'%f_max)
        aout.close()
        out_dict = {}
        out_dict['metadata'] = metadata
        out_dict['fx'] = fx
        return out_dict

##
# compute_genomic_homozygosity_from_snp() calculates genomic homozygosity for each SNP genotype as the proportion
# of homozygous loci.
# @param pedobj A PyPedal pedigree object.
# @param debug Turns d ebug messages on (True) and off (False).
# @retval A PyPedal pedigree object with genomic inbreeding coefficients assigned.
def renumber_snp_ids(pedobj):
    """
    :param pedobj:
    :param debug:
    :return:
    """

    #if pedobj.snp.empty:
    if pedobj.snp:
        if pedobj.kw['debug_messages'] == 1:
            print '[INFO]: pyp_snp/renumber_snp_ids(): Renumbering animal IDs in the SNP dataframe.'
        logging.info('pyp_snp/renumber_snp_ids(): Renumbering animal IDs in the SNP dataframe')
        for p in pedobj.pedigree:
            pedobj.snp.loc[:, ['animalID']].replace(to_replace=p.originalID, value=p.animalID, inplace=True)
    else:
        if pedobj.kw['debug_messages'] == 1:
            print '[ERROR]: pyp_snp/renumber_snp_ids(): There are no SNP data associated with this pedigree so no ' \
                  'IDs need to be renumbered.'
        logging.error('pyp_snp/renumber_snp_ids(): There are no SNP data associated with this pedigree so no IDs '
                      'need to be renumbered.')
    return

##
# generate_random_genotype() generates a random string of 0s, 1, and 2s. The resulting sting has no actual
# genealogical significance -- it's just a random string for use as simple test data!
# @param string_length The length of the strength to generate.
# @retval A string of length <string_length> made up of 0s, 1s, and 2s.
def generate_random_genotype(string_length):
    random_genotype = ''.join([random.choice(['0', '1', '2']) for i in xrange(string_length)])
    return random_genotype