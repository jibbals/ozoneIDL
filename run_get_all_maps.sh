#!/bin/bash

for i in `seq 1 6`; do
    qsub -v part=$i qsub_get_all_maps.sh
done
   

