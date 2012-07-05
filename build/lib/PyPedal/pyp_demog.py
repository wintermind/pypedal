#!/usr/bin/python

###############################################################################
# NAME: pyp_demog.py
# VERSION: 2.0.0 (29SEPTEMBER2010)
# AUTHOR: John B. Cole, PhD (john.cole@ars.usda.gov)
# LICENSE: LGPL
###############################################################################
# FUNCTIONS:
#     set_base_year()
#     set_age_units()
#     age_distribution()
#     sex_ratio()
#     founders_by_year()
###############################################################################

## @package pyp_demog
# pyp_demog contains a set of procedures for demographic calculations on the
# population describe in a pedigree.

import pyp_utils

# Define some globals in case the user forgets to call set_base_year() and
# set_age_units().
global BASE_DEMOGRAPHIC_YEAR
global BASE_DEMOGRAPHIC_UNIT
global SEX_CODE_MAP
BASE_DEMOGRAPHIC_YEAR = 1900
BASE_DEMOGRAPHIC_UNIT = 'year'
SEX_CODE_MAP = {'m':'Male','f':'Female','u':'Unk'}

##
# set_base_year() defines a global variable, BASE_DEMOGRAPHIC_YEAR.
# @param year The year to be used as a base for computing ages.
# @retval None
def set_base_year(year=1900):
    """
    set_base_year() defines a global variable, BASE_DEMOGRAPHIC_YEAR.
    """
    global BASE_DEMOGRAPHIC_YEAR
    BASE_DEMOGRAPHIC_YEAR = year

##
# set_age_units() defines a global variable, BASE_DEMOGRAPHIC_UNIT.
# @param units The base unit for age computations ('year'|'month'|'day').
# @retval None
def set_age_units(units='year'):
    """
    set_age_units() defines a global variable, BASE_DEMOGRAPHIC_UNIT.
    """
    _units = ['year','month','day']
    global BASE_DEMOGRAPHIC_UNIT
    if units in _units:
        BASE_DEMOGRAPHIC_UNIT = units
    else:
        BASE_DEMOGRAPHIC_UNIT = 'year'

##
# age_distribution() computes histograms of the age distribution of
# males and females in the population.  You can also stratify by
# sex to get individual histograms.
# @param pedobj An instance of a PyPedal NewPedigree object.
# @param sex A flag which determines whether or not to stratify by sex.
# @retval None
def age_distribution(pedobj,sex=1):
    """
    age_distribution() computes histograms of the age distribution of
    males and females in the population.  You can also stratify by
    sex to get individual histograms.
    """
    age_dict = {}
    age_freq_total = 0.0
    if pedobj.pedigree[0].age == -999:
        if pedobj.pedigree[0].igen == -999:
            pyp_utils.set_generation(pedobj.pedigree)
        pyp_utils.set_age(pedobj.pedigree)
    if not sex:
        for i in range(len(pedobj.pedigree)):
            try:
                age_dict[pedobj.pedigree[i].age] = age_dict[pedobj.pedigree[i].age] + 1
            except KeyError:
                age_dict[pedobj.pedigree[i].age] = 1
        age_hist = pyp_utils.simple_histogram_dictionary(age_dict)
        if pedobj.kw['debug_messages']:
            print '-'*80
            print 'Population Age Distribution'
            print '-'*80
            print '\tAge\tCount\tFrequency\tHistogram'
            for key in age_dict.keys():
                age_freq_total = age_freq_total + float(age_dict[key])/float(len(pedobj.pedigree))
                print '\t%s\t%s\t%s\t%s' % (key,age_dict[key],float(age_dict[key])/float(len(pedobj.pedigree)),age_hist[key])
            print '\tTOTAL\t%s\t%s' % (len(pedobj.pedigree),age_freq_total)
            print '-'*80
    else:
        males = []
        females = []
        unknowns = []
        male_dict = {}
        female_dict = {}
        unknown_dict = {}
        for i in range(len(pedobj.pedigree)):
            if pedobj.pedigree[i].sex == 'm':
                males.append(pedobj.pedigree[i])
            elif pedobj.pedigree[i].sex == 'f':
                females.append(pedobj.pedigree[i])
            else:
                unknowns.append(pedobj.pedigree[i])
        for m in range(len(males)):
            try:
                male_dict[males[m].age] = male_dict[males[m].age] + 1
            except KeyError:
                male_dict[males[m].age] = 1
        for f in range(len(females)):
            try:
                female_dict[females[f].age] = female_dict[females[f].age] + 1
            except KeyError:
                female_dict[females[f].age] = 1
        for u in range(len(unknowns)):
            try:
                unknown_dict[unknowns[u].age] = unknown_dict[unknowns[u].age] + 1
            except KeyError:
                unknown_dict[unknowns[u].age] = 1
        male_hist = pyp_utils.simple_histogram_dictionary(male_dict)
        female_hist = pyp_utils.simple_histogram_dictionary(female_dict)
        unknown_hist = pyp_utils.simple_histogram_dictionary(unknown_dict)
        if pedobj.kw['messages'] == 'verbose':
            print '-'*80
            print 'Population Age Distribution by Sex'
            print '-'*80
            age_freq_total = 0.0
            print 'Males'
            print '\tAge\tCount\tFrequency\tHistogram'
            for key in male_dict.keys():
                age_freq_total = age_freq_total + float(male_dict[key])/float(len(males))
                print '\t%s\t%s\t%s\t%s' % (key,male_dict[key],float(male_dict[key])/float(len(males)),male_hist[key])
            print '\tTOTAL\t%s\t%s' % (len(males),age_freq_total)
            print '-'*80
            age_freq_total = 0.0
            print 'Females'
            print '\tAge\tCount\tFrequency\tHistogram'
            for key in female_dict.keys():
                age_freq_total = age_freq_total + float(female_dict[key])/float(len(females))
                print '\t%s\t%s\t%s\t%s' % (key,female_dict[key],float(female_dict[key])/float(len(females)),female_hist[key])
            print '\tTOTAL\t%s\t%s' % (len(females),age_freq_total)
            print '-'*80
            age_freq_total = 0.0
            print 'Unknowns'
            print '\tAge\tCount\tFrequency\tHistogram'
            for key in unknown_dict.keys():
                age_freq_total = age_freq_total + float(unknown_dict[key])/float(len(unknowns))
                print '\t%s\t%s\t%s\t%s' % (key,unknown_dict[key],float(unknown_dict[key])/float(len(unknowns)),unknown_hist[key])
            print '\tTOTAL\t%s\t%s' % (len(unknowns),age_freq_total)
            print '-'*80

##
# sex_ratio() returns a dictionary containing the proportion of males and females in the population.
# @param pedobj An instance of a PyPedal NewPedigree object.
# @retval dict A dictionary containing entries for each sex/gender code defined in the global SEX_CODE_MAP.
def sex_ratio(pedobj):
    """
    sex_ratio() returns a dictionary containing the proportion of males and females in
    the population.
    """
    sexratiodict = {}
    for s in SEX_CODE_MAP.keys():
        sexratiodict[s] = 0
    for i in range(len(pedobj.pedigree)):
        if sexratiodict.has_key(pedobj.pedigree[i].sex):
            sexratiodict[pedobj.pedigree[i].sex] = sexratiodict[pedobj.pedigree[i].sex] + 1
        else:
            sexratiodict[pedobj.pedigree[i].sex] = 1
    if pedobj.kw['messages'] == 'verbose':
        print '-'*80
        print 'Overall Sex Ratio'
        print '-'*80
        print '(n = %s)' % (len(pedobj.pedigree))
        print 'Sex\tCount\tFrequency'
        for s in sexratiodict.keys():
            print '%s:\t%s\t%s' % (SEX_CODE_MAP[s],sexratiodict[s],float(sexratiodict[s])/float(len(pedobj.pedigree)))
        print '-'*80
        if int(sexratiodict['u']) > 0:
            marginal = sexratiodict['m'] + sexratiodict['f']
            print 'Conditional Sex Ratio'
            print '-'*80
            print '(n = %s)' % (marginal)
            print 'Sex\tCount\tFrequency'
            print '%s:\t%s\t%s' % (SEX_CODE_MAP['m'],sexratiodict['m'],float(sexratiodict['m'])/marginal)
            print '%s:\t%s\t%s' % (SEX_CODE_MAP['f'],sexratiodict['f'],float(sexratiodict['f'])/marginal)
    return sexratiodict

##
# founders_by_year() returns a dictionary containing the number of founders in each
# birthyear.
# @param pedobj A PyPedal pedigree object.
# @retval dict A dictionary containing entries for each sex/gender code defined in the global SEX_CODE_MAP.
def founders_by_year(pedobj):
    """
    founders_by_year() returns a dictionary containing the number of founders in each
    birthyear.
    """
    founderbyyeardict = {}
    if 'b' not in pedobj.kw['pedformat'] and 'y' not in pedobj.kw['pedformat']:
        # Birthyears were not provided with the pedigree.  There only birthdate in
        # the pedigree is the proxy year of 1900.
        founderbyyeardict[BASE_DEMOGRAPHIC_YEAR] = pedobj.metadata.num_unique_founders
    else:
        for _f in pedobj.metadata.unique_founder_list:
            _by = pedobj.pedigree[int(_f)-1].by
            #print _by
            try:
                founderbyyeardict[_by] = founderbyyeardict[_by] + 1
            except KeyError:
                founderbyyeardict[_by] = 1
    # If the dictionary has more than one birthyear in it we should iterate through the list
    # to fill in any gaps in years.  This will make downstream graphing much easier.
    #print founderbyyeardict.keys()
    _years = founderbyyeardict.keys()
    _years.sort()
    #print _years
    for _f in range(_years[0],_years[-1]):
        try:
            _c = founderbyyeardict[_f]
        except KeyError:
            founderbyyeardict[_f] = 0
    #print founderbyyeardict
    return founderbyyeardict
