#!/bin/bash
#PBS -P m19
#PBS -q express
#PBS -N make_maps
#PBS -l walltime=01:00:00
#PBS -l mem=40000MB
#PBS -l cput=01:30:00
#PBS -l jobfs=10MB
#PBS -l wd
#PBS -l ncpus=4
#PBS -l software=idl
#PBS -j oe

if [ -z ${part+x} ]; then 
    echo "part is unset"
    echo "EG: qsub -v part=1 $0"
fi

# Make a dummy display  buffer
Xvfb :99 &
export DISPLAY=:99

idl -e "get_all_maps, parts=$part"

