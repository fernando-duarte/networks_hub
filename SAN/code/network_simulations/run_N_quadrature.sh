for i in {3..15}
do
    echo "Running for N = $i"
    sh run_batch_quadrature.sh quadrature_timing.jl s3://ran-s3-shared/benchmark/ $i
done

