#!/bin/bash
# first argument $1 is the name of the julia file 
# second argument $2 is the address of the s3 container to transfer file
# third argument $3 is the number of nodes

# setting up access (only needed on first use)
aws s3 cp s3://ran-s3-install-pkgs/config/RanPocKP.pem .
chmod 400 RanPocKP.pem

# copy julia script to s3
aws s3 cp "$1" "$2"

sudo su - ec2-user <<ENDSUDO
ssh -o StrictHostKeyChecking=no -i RanPocKP.pem ec2-user@172.31.100.6 <<ENDSSH
aws s3 cp "$2""$1" ~/"$1"
#nohup time julia ~/"$1">output_"$1".txt
#nohup time env JULIA_NUM_THREADS=4 julia ~/"$1">output_"$1".txt #4 threads
#nohup time julia -p 4 ~/"$1">output_"$1".txt  # 4 workers
nohup time julia ~/"$1">output_"$3"_"$1".txt $3 #4 threads and the right amount of workers given the threads
aws s3 cp /home/ec2-user/output_"$3"_"$1".txt "$2"output_"$3"_"$1".txt
aws s3 cp /home/ec2-user/quadrature_results_"$3".csv "$2"quadrature_results_"$3".csv

ENDSSH
ENDSUDO
