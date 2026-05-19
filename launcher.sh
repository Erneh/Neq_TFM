#!/bin/bash
MODIFIER_ID=$1
N_POT=$2
E=$3
TEMP=$4
MU=$5
GAMMA=$6
M=$7
NR=$8
NT=$9
NMEAST=${10}
STPERT=${11}
HAMTYPE=${12}
if [ $HAMTYPE = "basic" ]; then
    HAMTYPE=""
fi
LOG_LOC="Out/${MODIFIER_ID}${HAMTYPE}/G=${GAMMA}_E=${E}_Temp=${TEMP}_mu=${MU}/N=${N_POT}_M=${M}_R=${NR}_nT=${NT}_measT=${NMEAST}_stT=${STPERT}"
mkdir -p $LOG_LOC
echo "log can be found in $LOG_LOC"
python3 -u neq.py "$@" > "$LOG_LOC/log.txt" 2>&1
