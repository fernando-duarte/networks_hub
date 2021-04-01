cd /home/ec2-user/SageMaker/quadrature_results

awk '(NR == 1) || (FNR > 1)' meta_data_*.csv > final_meta_data.csv #concatinating meta data
aws s3 cp final_meta_data.csv s3://ran-s3-shared/benchmark/quadrature_meta_data.csv