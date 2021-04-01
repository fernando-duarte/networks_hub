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

#mkdir -p benchmark_timing
#aws s3 cp "$2"benchmark_timing benchmark_timing/ --recursive 
#nohup time ~/julia-1.6.0/bin/julia  ~/"$1">output_"$1".txt
#nohup time env JULIA_NUM_THREADS=4 julia ~/"$1">output_"$1".txt # 4 threads
#nohup time ~/julia-1.6.0/bin/julia -p 4 ~/"$1">output_"$1".txt  # 4 workers
nohup time julia ./"$1">output_"$3"_"$1".txt $3 # calling julia and creating output file
aws s3 cp /home/ec2-user/output_"$3"_"$1".txt "$2"benchmark_results/output_"$3"_"$1".txt # copy file over
aws s3 cp /home/ec2-user/data_out_"$3".csv "$2"benchmark_results/data_out_"$3".csv #copy data file which stores variable solutions
aws s3 cp /home/ec2-user/variables_"$3".jld2 "$2"benchmark_results/variables_"$3".jld2 #copy file which stores variable solutions as .jld2 file
aws s3 cp /home/ec2-user/profile_"$3".txt "$2"benchmark_results/profile_"$3".txt #copies over profiling results
aws s3 cp /home/ec2-user/meta_data_"$3".csv "$2"benchmark_results/meta_data_"$3".csv #copies over data on the run (timing, mem, etc) 

ENDSSH
ENDSUDO

# copy output from s3 to local folder
aws s3 cp "$2"benchmark_results/output_"$3"_"$1".txt $OUT_DIR
aws s3 cp "$2"benchmark_results/data_out_"$3".csv $OUT_DIR
aws s3 cp "$2"benchmark_results/variables_"$3".jld2 $OUT_DIR
aws s3 cp "$2"benchmark_results/profile_"$3".txt $OUT_DIR
aws s3 cp "$2"benchmark_results/meta_data_"$3".csv $OUT_DIR

# remove julia script and output files from s3
# aws s3 rm "$2""$1"
# aws s3 rm "$2"output_"$1".txt
# aws s3 rm "$2"/* --recursive 
