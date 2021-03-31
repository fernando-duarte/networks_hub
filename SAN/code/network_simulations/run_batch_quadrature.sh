#!/bin/bash
# first argument $1 is the name of the julia file 
# second argument $2 is the address of the s3 container to transfer file
# third argument $3 is the number of nodes

# folder to save results
OUT_DIR="./quadrature_timing/"
mkdir -p $OUT_DIR

# setting up access (only needed on first use)
aws s3 cp s3://ran-s3-install-pkgs/config/RanPocKP.pem .
chmod 400 RanPocKP.pem

# copy package environment
aws s3 cp /home/ec2-user/SageMaker/networks_hub/SAN/code/network_simulations/quadrature_timing/Project.toml "$2"quadrature_timing/Project.toml
aws s3 cp /home/ec2-user/SageMaker/networks_hub/SAN/code/network_simulations/quadrature_timing/Manifest.toml "$2"quadrature_timing/Manifest.toml

# copy julia script to s3
aws s3 cp "$1" "$2"

sudo su - ec2-user <<ENDSUDO
ssh -o StrictHostKeyChecking=no -i RanPocKP.pem ec2-user@172.31.100.6 <<ENDSSH
aws s3 cp "$2""$1" ~/"$1"
# copy package environment
mkdir -p quadrature_timing
aws s3 cp "$2"quadrature_timing quadrature_timing/ --recursive 
#nohup time ~/julia-1.6.0/bin/julia  ~/"$1">output_"$1".txt
#nohup time env JULIA_NUM_THREADS=4 julia ~/"$1">output_"$1".txt # 4 threads
#nohup time ~/julia-1.6.0/bin/julia  -p 4 ~/"$1">output_"$1".txt  # 4 workers
nohup time julia ~/"$1">output_"$3"_"$1".txt $3 # calling julia
aws s3 cp /home/ec2-user/output_"$3"_"$1".txt "$2"output_"$3"_"$1".txt
aws s3 cp /home/ec2-user/meta_data"$3".csv "$2"meta_data_"$3".csv

ENDSSH
ENDSUDO

# copy output from s3 to local
aws s3 cp "$2"output_"$3"_"$1".txt $OUT_DIR
aws s3 cp "$2"meta_data_"$3".csv $OUT_DIR

# remove julia script and output files from s3
# aws s3 rm "$2""$1"
# aws s3 rm "$2"output_"$1".txt
# aws s3 rm "$2"/* --recursive 