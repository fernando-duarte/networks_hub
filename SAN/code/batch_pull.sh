#!/bin/bash
#$ -cwd
#$ -m abe
#$ -M collin.jones@ny.frb.org
echo "Starting job at `date`"
python /home/frb-ny/ruelaf/WRDS_query.py
#python /home/frb-ny/fduarte/Fire_Sale_Pull/WRDS_query.py
echo "Finished job at `date`"