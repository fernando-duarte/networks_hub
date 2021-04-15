for N in {5..40..5}
do
    echo "Running for N = $N"
    sh run_batch_cuba.sh cuba_beta.jl s3://ran-s3-shared/benchmark/ $N 
done
