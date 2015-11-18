#! /bin/bash

# file: next.sh

# This bash script analyzes the record named in its command-line argument ($1),
# saving the results as an annotation file for the record, with annotator 'qrs'.
# This script is run once for each record in the Challenge test set.

# For example, if invoked as
#    next.sh 100
# it analyzes record 100 using 'gqrs', and saves the results as an annotation
# file named '100.qrs'.
#cat ~/.octaverc
OCTAVE='octave --quiet --eval '
#OCTAVE='octave --eval '
RECORD=$1
RPATH=`pwd`
STR="${OCTAVE} \"eval_record('$RECORD'); quit;\" "
echo "$STR"
eval ${STR}
