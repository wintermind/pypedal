#!/bin/bash
# Make sure the Doxygen tree is up-to-date
doxygen doxygen.conf

cd PyPedal/doc

sh make_pypedal_manual.sh

cd -
