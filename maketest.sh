#!/bin/bash
# Quick shell script to run all of the examples in PyPedal/examples
# to check for breakage before a release is packaged. I know, it's
# not a real replacement for proper unit testing, but at least it's
# something.

cd PyPedal/examples
# python new_amatrix.py
python new_classes.py
python new_db.py
python new_doug.py
python new_format.py
python new_graphics.py
python new_graphics2.py
python new_hartl.py
python new_ids.py
python new_inbreeding.py
python new_inbreeding2.py
python new_jbc.py
python new_lacy.py
python new_methods.py
python new_networkx.py
python new_options.py
python new_renumbering.py
python new_reporting.py
python new_simulate.py
