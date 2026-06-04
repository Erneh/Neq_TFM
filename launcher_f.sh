#!/bin/bash
MODIFIER_ID=$1
N_POT=$2
E=$3
TEMP=$4
MU=$5
GAMMA=$6
M=$7
NR=$8
PATH_TYPE=$9
NK=${10}
NT=${11}
NMEAST=${12}
STPERT=${13}
HAMTYPE=${14}
MASS=${15}
if [ $HAMTYPE = "basic" ]; then
    HAMTYPE=""
    LOG_LOC="ARPES/Out/${MODIFIER_ID}${HAMTYPE}/${PATH_TYPE}/G=${GAMMA}_E=${E}_Temp=${TEMP}_mu=${MU}/N=${N_POT}_M=${M}_nT=${NT}_measT=${NMEAST}_stT=${STPERT}_nk=${NK}/logs"
elif [ $HAMTYPE = "hbn" ]; then
    LOG_LOC="ARPES/Out/${MODIFIER_ID}${HAMTYPE}/${PATH_TYPE}/G=${GAMMA}_E=${E}_Temp=${TEMP}_mu=${MU}_m=${MASS}/N=${N_POT}_M=${M}_nT=${NT}_measT=${NMEAST}_stT=${STPERT}_nk=${NK}/logs"  
elif [ $HAMTYPE = "jclhbn" ]; then
    LOG_LOC="ARPES/Out/${MODIFIER_ID}${HAMTYPE}/${PATH_TYPE}/G=${GAMMA}_E=${E}_Temp=${TEMP}_mu=${MU}_m=${MASS}/N=${N_POT}_M=${M}_nT=${NT}_measT=${NMEAST}_stT=${STPERT}_nk=${NK}/logs"  

else
    LOG_LOC="ARPES/Out/${MODIFIER_ID}${HAMTYPE}/${PATH_TYPE}/G=${GAMMA}_E=${E}_Temp=${TEMP}_mu=${MU}/N=${N_POT}_M=${M}_nT=${NT}_measT=${NMEAST}_stT=${STPERT}_nk=${NK}/logs"
fi

mkdir -p $LOG_LOC
echo "log can be found in $LOG_LOC"
N_LOGS=$(ls -1 ${LOG_LOC} | wc -l)
python3 -u ARPES/neq_f.py "$@" > "$LOG_LOC/log$N_LOGS.txt" 2>&1

# Eliminate the log in the case the code finds calculations already done
EXIT_CODE=${PIPESTATUS[0]}
if [ $EXIT_CODE = 42 ]; then
    rm "$LOG_LOC/log$N_LOGS.txt"
    echo "Calculation made, log erased"
fi