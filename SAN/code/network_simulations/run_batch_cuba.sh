#!/bin/bash
# first argument $1 is the name of the julia file 
# second argument $2 is the address of the s3 container to transfer file
# third argument $3 is the number of nodes

# folder to save results
OUT_DIR="./cuba_results/"
mkdir -p $OUT_DIR

# setting up access (only needed on first use)
aws s3 cp s3://ran-s3-install-pkgs/config/RanPocKP.pem .
chmod 400 RanPocKP.pem

# copy package environment
aws s3 cp  /home/ec2-user/SageMaker/networks_hub/SAN/code/network_simulations/joint_timing/Project.toml "$2"Project.toml
aws s3 cp  /home/ec2-user/SageMaker/networks_hub/SAN/code/network_simulations/joint_timing/Manifest.toml "$2"Manifest.toml
aws s3 cp  /home/ec2-user/SageMaker/networks_hub/SAN/code/network_simulations/joint_timing/src/joint_timing.jl "$2"src/joint_timing.jl


mkdir -p /home/ec2-user/joint_timing
mkdir -p /home/ec2-user/joint_timing/src


cp /home/ec2-user/SageMaker/networks_hub/SAN/code/network_simulations/joint_timing/Project.toml /home/ec2-user/joint_timing/Project.toml
cp /home/ec2-user/SageMaker/networks_hub/SAN/code/network_simulations/joint_timing/Manifest.toml /home/ec2-user/joint_timing/Manifest.toml
cp /home/ec2-user/SageMaker/networks_hub/SAN/code/network_simulations/joint_timing/src/joint_timing.jl /home/ec2-user/joint_timing/src/joint_timing.jl

# copy julia script to s3
aws s3 cp "$1" "$2"

sudo su - ec2-user <<ENDSUDO
ssh -o StrictHostKeyChecking=no -i RanPocKP.pem ec2-user@172.31.100.6 <<ENDSSH
aws s3 cp "$2""$1" ~/"$1"
# copy package environment
mkdir -p quadrature_timing
#aws s3 cp "$2"quadrature_timing quadrature_timing/ --recursive 
#nohup time ~/julia-1.6.0/bin/julia  ~/"$1">output_"$1".txt
#nohup time env JULIA_NUM_THREADS=4 julia ~/"$1">output_"$1".txt # 4 threads
#nohup time ~/julia-1.6.0/bin/julia  -p 4 ~/"$1">output_"$1".txt  # 4 workers
nohup time  env JULIA_NUM_THREADS=4 julia ~/"$1">output_"$3"_"$1".txt $3 # calling julia
aws s3 cp /home/ec2-user/output_"$3"_"$1".txt "$2"cuba_results/output_"$3"_"$1".txt
aws s3 cp /home/ec2-user/data_out_"$3"_cuhre.csv "$2"cuba_results/data_out_"$3"_cuhre.csv
aws s3 cp /home/ec2-user/data_out_"$3"_divonne.csv "$2"cuba_results/data_out_"$3"_divonne.csv
aws s3 cp /home/ec2-user/data_out_"$3"_suave.csv "$2"cuba_results/data_out_"$3"_suave.csv

ENDSSH
ENDSUDO

# copy output from s3 to local
aws s3 cp "$2"cuba_results/output_"$3"_"$1".txt $OUT_DIR
aws s3 cp "$2"cuba_results/data_out_"$3"_cuhre.csv $OUT_DIR
aws s3 cp "$2"cuba_results/data_out_"$3"_suave.csv $OUT_DIR
aws s3 cp "$2"cuba_results/data_out_"$3"_divonne.csv $OUT_DIR

# remove julia script and output files from s3
# aws s3 rm "$2""$1"
# aws s3 rm "$2"output_"$1".txt
# aws s3 rm "$2"/* --recursive 

