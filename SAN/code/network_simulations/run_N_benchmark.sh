for i in {2..8}
do
    echo "Running for N = $i"
    sh run_batch_benchmark.sh benchmark_timing.jl s3://ran-s3-shared/benchmark/ $i
done 