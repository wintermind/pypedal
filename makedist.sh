#!/bin/bash
# Quickie to brew up PyPedal distributions
python setup.py bdist_egg
python setup.py sdist --formats=gztar,zip
python setup.py register
