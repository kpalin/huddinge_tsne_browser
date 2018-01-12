#!/bin/bash

#  Fri Jan 12 15:33:36 2018
#  kpalin@merit.ltdk.helsinki.fi

# For debugging
#set -o verbose 

# Die on unset variables
set -o nounset
# Die on errors
set -o errexit
# Die if any part of a pipe fails
set -o pipefail


python setup.py install --single-version-externally-managed --record=record.txt

cd MODER
make
mv all_pairs_huddinge $PREFIX/bin/
