###########################################################################
#    Copyright (C) 2010 by John B. Cole, PhD
#    <john.cole@ars.usda.gov>
#
# Copyright: See COPYING file that comes with this distribution
###########################################################################
import ez_setup
ez_setup.use_setuptools()
from setuptools import setup, find_packages

setup(
    name = "PyPedal",
    version = "2.0.4",
    packages = find_packages(),

    classifiers = [
	'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Mathematics',
      ],
    author = "John B. Cole, PhD",
    author_email = "john.cole@ars.usda.gov",
    description = "Tools for pedigree analysis",
    long_description="PyPedal provides an API that may be used to maniuplate pedigrees in a number of ways.  Key metrics include coefficients of inbreeding and relationship; effective founder and ancestor numbers; and expected inbreeding from a given mating.  If you have Graphviz and pydot installed, PyPedal can be used to produce a graph from your pedigree.",
    license = "GNU LGPL",
    keywords = "Python pedigree genetic analysis diversity",
    url = "http://pypedal.sourceforge.net/",
)
