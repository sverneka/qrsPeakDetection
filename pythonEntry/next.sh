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
RESAMPLE_FILES="xform -i $RECORD -S fileHeaderFormat.txt"

DAT_TO_TEXT="rdsamp -r file > $RECORD'.txt'"

EVAL_RECORD="python sources/eval_record_python.py ${RECORD} "
WRITE_ANNOTATIONS="${OCTAVE} \"writeAnnotations('$RECORD'); quit;\" "

eval ls
eval pwd


echo ${RESAMPLE_FILES}
eval ${RESAMPLE_FILES}

echo ${DAT_TO_TEXT}
eval ${DAT_TO_TEXT}
#eval ${SAVE_FILE_DESCRIPTION}

echo "$EVAL_RECORD"
eval ${EVAL_RECORD}

echo "$WRITE_ANNOTATIONS"
eval ${WRITE_ANNOTATIONS}

