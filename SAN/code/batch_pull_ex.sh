#!/bin/bash
#$ -cwd
#$ -m abe
#$ -M collin.jones@ny.frb.org
echo "Starting job at `date`"
# Create (or uncomment) line below that matches your WRDS cloud file structure.
# python /home/frb-ny/ruelaf/WRDS_query.py
# python /home/berkeley/cdj345/WRDS_query.py
#python /home/frb-ny/fduarte/Fire_Sale_Pull/WRDS_query.py
echo "Finished job at `date`"
