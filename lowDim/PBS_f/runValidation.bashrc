#!/bin/bash

for args in `seq 2 58`;
do
  qsub  -v args=$args validation.PBS
     echo $args
done
