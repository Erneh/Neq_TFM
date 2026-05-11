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
LOG_LOC="Out/$MODIFIER_ID/G=${GAMMA}_E=${E}_Temp=${TEMP}_mu=${MU}/N=${N_POT}_M=${M}_R=${NR}_nT=${NT}_measT=${NMEAST}_stT=${STPERT}"
mkdir -p $LOG_LOC
echo "Launching calculation with:"
echo "- Number of atoms = 2^$N_POT"
echo "- Energy = ${E} eV"
echo "- Temperature = ${TEMP} K"
echo "- Chemical potential = $MU eV"
echo "- Gamma = $GAMMA"
echo "- # of moments = $M"
echo "- # of random vectors = $NR"
echo "- # of periods simulated = $NT"
echo "- # of measures per period = $NMEAST"
echo "- # of time steps per period = $STPERT"
echo "log can be found in $LOG_LOC"
python3 -u neq.py "$@" > "$LOG_LOC/log.txt" 2>&1
echo "Finished calculation with:"
echo "- Number of atoms = 2^$N_POT"
echo "- Energy = ${E} eV"
echo "- Temperature = ${T} K"
echo "- Chemical potential = $MU eV"
echo "- Gamma = $GAMMA"
echo "- # of moments = $M"
echo "- # of random vectors = $NR"
echo "- # of periods simulated = $NT"
echo "- # of measures per period = $NMEAST"
echo "- # of time steps per period = $STPERT"
