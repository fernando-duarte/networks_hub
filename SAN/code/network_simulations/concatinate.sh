cd /home/ec2-user/SageMaker/benchmark_results

awk '(NR == 1) || (FNR > 1)' meta_data_*.csv > final_meta_data.csv #concatinating meta data
aws s3 cp final_meta_data.csv s3://ran-s3-shared/benchmark/final_meta_data.csv