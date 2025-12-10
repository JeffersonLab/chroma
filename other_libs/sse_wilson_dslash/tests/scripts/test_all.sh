# Save original testvol files
mv ../sse_wilson_dslash/tests/testvol.h ../sse_wilson_dslash/tests/testvol.h.bak

cat > ../sse_wilson_dslash/tests/testvol.h << EOF
const int nrow_in[] = { 2,6,6,2 };
EOF
make ; make check
mpirun -np 72 -mca btl sm ./tests/t_dslash -geom 1 6 6 2
mpirun -np 36 --mca btl sm ./tests/t_dslash -geom 1 6 6 1
mpirun -np 9 --mca btl sm ./tests/t_dslash  -geom 1 3 3 1
mpirun -np 6 --mca btl sm ./tests/t_dslash  -geom 1 3 2 1
mpirun -np 6 --mca btl sm ./tests/t_dslash  -geom 1 2 3 1
mpirun -np 4 --mca btl sm ./tests/t_dslash  -geom 1 2 2 1
mpirun -np 2 --mca btl sm ./tests/t_dslash  -geom 1 2 1 1
mpirun -np 2 --mca btl sm ./tests/t_dslash  -geom 1 1 2 1
mpirun -np 2 --mca btl sm ./tests/t_dslash  -geom 1 1 1 2
mpirun -np 1 --mca btl sm ./tests/t_dslash  -geom 1 1 1 1


cat > ../sse_wilson_dslash/tests/testvol.h << EOF
const int nrow_in[] = { 2,6,6,4 };
EOF
make ; make check

mpirun -np 72 --mca btl sm ./tests/t_dslash -geom 1 6 6 2
mpirun -np 36 --mca btl sm ./tests/t_dslash -geom 1 6 6 1
mpirun -np 36 --mca btl sm ./tests/t_dslash -geom 1 3 3 4
mpirun -np 18 --mca btl sm ./tests/t_dslash -geom 1 3 3 2
mpirun -np 9 --mca btl sm ./tests/t_dslash  -geom 1 3 3 1
mpirun -np 6 --mca btl sm ./tests/t_dslash  -geom 1 3 2 1
mpirun -np 6 --mca btl sm ./tests/t_dslash  -geom 1 2 3 1
mpirun -np 4 --mca btl sm ./tests/t_dslash  -geom 1 2 2 1
mpirun -np 2 --mca btl sm ./tests/t_dslash  -geom 1 2 1 1
mpirun -np 2 --mca btl sm ./tests/t_dslash  -geom 1 1 2 1
mpirun -np 2 --mca btl sm ./tests/t_dslash  -geom 1 1 1 2
mpirun -np 1 --mca btl sm ./tests/t_dslash  -geom 1 1 1 1

cat > ../sse_wilson_dslash/tests/testvol.h << EOF
const int nrow_in[] = { 2,2,2,2 };
EOF
make ; make check

mpirun -np 8 --mca btl sm ./tests/t_dslash  -geom 1 2 2 2
mpirun -np 4 --mca btl sm ./tests/t_dslash  -geom 1 2 2 1
mpirun -np 2 --mca btl sm ./tests/t_dslash  -geom 1 2 1 1
mpirun -np 2 --mca btl sm ./tests/t_dslash  -geom 1 1 2 1
mpirun -np 2 --mca btl sm ./tests/t_dslash  -geom 1 1 1 2
mpirun -np 1 --mca btl sm ./tests/t_dslash  -geom 1 1 1 1

cat > ../sse_wilson_dslash/tests/testvol.h << EOF
const int nrow_in[] = { 2,2,2,4 };
EOF
make ; make  check

mpirun -np 16 --mca btl sm ./tests/t_dslash -geom 1 2 2 4
mpirun -np 8 --mca btl sm ./tests/t_dslash  -geom 1 2 2 2
mpirun -np 4 --mca btl sm ./tests/t_dslash  -geom 1 2 2 1
mpirun -np 2 --mca btl sm ./tests/t_dslash  -geom 1 2 1 1
mpirun -np 2 --mca btl sm ./tests/t_dslash  -geom 1 1 2 1
mpirun -np 2 --mca btl sm ./tests/t_dslash  -geom 1 1 1 2
mpirun -np 1 --mca btl sm ./tests/t_dslash  -geom 1 1 1 1

cat > ../sse_wilson_dslash/tests/testvol.h << EOF
const int nrow_in[] = { 2,2,2,6 };
EOF
make ; make  check

mpirun -np 24 --mca btl sm ./tests/t_dslash -geom 1 2 2 6
mpirun -np 12 --mca btl sm ./tests/t_dslash -geom 1 2 2 3
mpirun -np 8 --mca btl sm ./tests/t_dslash  -geom 1 2 2 2
mpirun -np 4 --mca btl sm ./tests/t_dslash  -geom 1 2 2 1
mpirun -np 2 --mca btl sm ./tests/t_dslash  -geom 1 2 1 1
mpirun -np 2 --mca btl sm ./tests/t_dslash  -geom 1 1 2 1
mpirun -np 2 --mca btl sm ./tests/t_dslash  -geom 1 1 1 2
mpirun -np 1 --mca btl sm ./tests/t_dslash  -geom 1 1 1 1

