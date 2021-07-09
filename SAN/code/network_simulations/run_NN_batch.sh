#!/bin/bash
# first argument $1 is the name of the julia file 
# second argument $2 is the address of the s3 container to transfer file
# third argument $3 is the number of nodes
# Have to run this file from /home/ec2-user/SageMaker

# folder to save results
OUT_DIR="./benchmark_results"
mkdir -p $OUT_DIR

# setting up access (only needed on first use)
aws s3 cp s3://ran-s3-install-pkgs/config/RanPocKP.pem .
chmod 400 RanPocKP.pem

# copy package environment
aws s3 cp /home/ec2-user/SageMaker/networks_hub/SAN/code/network_simulations/benchmark_timing/Project.toml "$2"Project.toml
aws s3 cp /home/ec2-user/SageMaker/networks_hub/SAN/code/network_simulations/benchmark_timing/Manifest.toml "$2"Manifest.toml

mkdir -p /home/ec2-user/joint_timing
mkdir -p /home/ec2-user/joint_timing/src

cp /home/ec2-user/SageMaker/networks_hub/SAN/code/network_simulations/joint_timing/Project.toml /home/ec2-user/joint_timing/Project.toml
cp /home/ec2-user/SageMaker/networks_hub/SAN/code/network_simulations/joint_timing/Manifest.toml /home/ec2-user/joint_timing/Manifest.toml

# copy julia script to s3
aws s3 cp "$1" "$2"

sudo su - ec2-user <<ENDSUDO
ssh -o StrictHostKeyChecking=no -i RanPocKP.pem ec2-user@172.31.100.6 <<ENDSSH
aws s3 cp "$2""$1" ~/"$1"
nohup time julia "$1">output_"$3"_"$1".txt $3 # calling julia and creating output file
aws s3 cp /home/ec2-user/output_"$3"_"$1".txt "$2"benchmark_results/output_"$3"_"$1".txt # copy output file over
aws s3 cp /home/ec2-user/clearing_p_NN_gpu_N"$3".bson "$2"benchmark_results/clearing_p_NN_gpu_N"$3".bson #copy NN file which stores variable solutions
ENDSSH
ENDSUDO

# copy output from s3 to local folder
aws s3 cp "$2"benchmark_results/output_"$3"_"$1".txt $OUT_DIR
aws s3 cp "$2"benchmark_results/clearing_p_NN_gpu_N"$3".bson $OUT_DIR

# remove julia script and output files from s3
# aws s3 rm "$2""$1"
# aws s3 rm "$2"output_"$1".txt
# aws s3 rm "$2"/* --recursive 