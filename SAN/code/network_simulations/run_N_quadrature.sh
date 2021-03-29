for i in {2..100}
do
    echo "Running for N = $i"
    sh run_batch_quadrature.sh quadrature_test.jl s3://ran-s3-shared/temp/ $i
done