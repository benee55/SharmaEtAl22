#!/bin/bash

niter=5
ens=1007
for args in `seq 3 7`;
do
    if [ "${args}" -eq "1" ]; then
        two=$(qsub MPI_cycle.PBS -v "args=$args $ens $niter")
    else
        two=$(qsub -W depend=afterany:$one MPI_cycle.PBS -v "args=$args $ens $niter")
    fi
    echo $two
    echo $args $niter $ens 
    one=$two
done
