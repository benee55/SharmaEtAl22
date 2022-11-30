#!/bin/bash

niter=6
ens=2015
for args in `seq 4 10`;
do
    if [ "${args}" -eq "4" ]; then
        two=$(qsub -v "args=$args $ens $niter" MPI_cycle.PBS )
    else
        two=$(qsub -W depend=afterany:$one -v "args=$args $ens $niter" MPI_cycle.PBS )
    fi
    echo $two
    echo $args $niter $ens 
    one=$two
done
