# Save original testvol files
mv ../sse_wilson_dslash/tests/testvol.h ../sse_wilson_dslash/tests/testvol.h.bak


cat > ../sse_wilson_dslash/tests/testvol.h << EOF
const int nrow_in[] = { 2,6,6,1 };
EOF
make ; make check

mpirun -np 36 --mca btl sm ./tests/t_dslash_3d -geom 1 6 6 1 
mpirun -np 9 --mca btl sm ./tests/t_dslash_3d  -geom 1 3 3 1
mpirun -np 6 --mca btl sm ./tests/t_dslash_3d  -geom 1 3 2 1
mpirun -np 6 --mca btl sm ./tests/t_dslash_3d  -geom 1 2 3 1
mpirun -np 4 --mca btl sm ./tests/t_dslash_3d  -geom 1 2 2 1
mpirun -np 2 --mca btl sm ./tests/t_dslash_3d  -geom 1 2 1 1
mpirun -np 2 --mca btl sm ./tests/t_dslash_3d  -geom 1 1 2 1
mpirun -np 1 --mca btl sm ./tests/t_dslash_3d  -geom 1 1 1 1

cat > ../sse_wilson_dslash/tests/testvol.h << EOF
const int nrow_in[] = { 2,6,6,2 };
EOF
make ; make check

mpirun -np 36 --mca btl sm ./tests/t_dslash_3d -geom 1 6 6 1
mpirun -np 9 --mca btl sm ./tests/t_dslash_3d  -geom 1 3 3 1
mpirun -np 6 --mca btl sm ./tests/t_dslash_3d  -geom 1 3 2 1
mpirun -np 6 --mca btl sm ./tests/t_dslash_3d  -geom 1 2 3 1
mpirun -np 4 --mca btl sm ./tests/t_dslash_3d  -geom 1 2 2 1
mpirun -np 2 --mca btl sm ./tests/t_dslash_3d  -geom 1 2 1 1
mpirun -np 2 --mca btl sm ./tests/t_dslash_3d  -geom 1 1 2 1
mpirun -np 1 --mca btl sm ./tests/t_dslash_3d  -geom 1 1 1 1

cat > ../sse_wilson_dslash/tests/testvol.h << EOF
const int nrow_in[] = { 2,6,6,3 };
EOF
make ; make check

mpirun -np 36 --mca btl sm ./tests/t_dslash_3d -geom 1 6 6 1
mpirun -np 9 --mca btl sm ./tests/t_dslash_3d  -geom 1 3 3 1
mpirun -np 6 --mca btl sm ./tests/t_dslash_3d  -geom 1 3 2 1
mpirun -np 6 --mca btl sm ./tests/t_dslash_3d  -geom 1 2 3 1
mpirun -np 4 --mca btl sm ./tests/t_dslash_3d  -geom 1 2 2 1
mpirun -np 2 --mca btl sm ./tests/t_dslash_3d  -geom 1 2 1 1
mpirun -np 2 --mca btl sm ./tests/t_dslash_3d  -geom 1 1 2 1
mpirun -np 1 --mca btl sm ./tests/t_dslash_3d  -geom 1 1 1 1

cat > ../sse_wilson_dslash/tests/testvol.h << EOF
const int nrow_in[] = { 2,6,6,4 };
EOF
make ; make check

mpirun -np 36 --mca btl sm ./tests/t_dslash_3d -geom 1 6 6 1
mpirun -np 9 --mca btl sm ./tests/t_dslash_3d  -geom 1 3 3 1
mpirun -np 6 --mca btl sm ./tests/t_dslash_3d  -geom 1 3 2 1
mpirun -np 6 --mca btl sm ./tests/t_dslash_3d  -geom 1 2 3 1
mpirun -np 4 --mca btl sm ./tests/t_dslash_3d  -geom 1 2 2 1
mpirun -np 2 --mca btl sm ./tests/t_dslash_3d  -geom 1 2 1 1
mpirun -np 2 --mca btl sm ./tests/t_dslash_3d  -geom 1 1 2 1
mpirun -np 1 --mca btl sm ./tests/t_dslash_3d  -geom 1 1 1 1

cat > ../sse_wilson_dslash/tests/testvol.h << EOF
const int nrow_in[] = { 2,2,2,1 };
EOF
make ; make check

mpirun -np 4 --mca btl sm ./tests/t_dslash_3d  -geom 1 2 2 1
mpirun -np 2 --mca btl sm ./tests/t_dslash_3d  -geom 1 2 1 1
mpirun -np 2 --mca btl sm ./tests/t_dslash_3d  -geom 1 1 2 1
mpirun -np 1 --mca btl sm ./tests/t_dslash_3d  -geom 1 1 1 1

cat > ../sse_wilson_dslash/tests/testvol.h << EOF
const int nrow_in[] = { 2,2,2,2 };
EOF
make ; make check

mpirun -np 4 --mca btl sm ./tests/t_dslash_3d  -geom 1 2 2 1
mpirun -np 2 --mca btl sm ./tests/t_dslash_3d  -geom 1 2 1 1
mpirun -np 2 --mca btl sm ./tests/t_dslash_3d  -geom 1 1 2 1
mpirun -np 1 --mca btl sm ./tests/t_dslash_3d  -geom 1 1 1 1

cat > ../sse_wilson_dslash/tests/testvol.h << EOF
const int nrow_in[] = { 2,2,2,3 };
EOF
make ; make  check

mpirun -np 4 --mca btl sm ./tests/t_dslash_3d  -geom 1 2 2 1
mpirun -np 2 --mca btl sm ./tests/t_dslash_3d  -geom 1 2 1 1
mpirun -np 2 --mca btl sm ./tests/t_dslash_3d  -geom 1 1 2 1
mpirun -np 1 --mca btl sm ./tests/t_dslash_3d  -geom 1 1 1 1

cat > ../sse_wilson_dslash/tests/testvol.h << EOF
const int nrow_in[] = { 2,2,2,4 };
EOF
make ; make check

mpirun -np 4 --mca btl sm ./tests/t_dslash_3d  -geom 1 2 2 1
mpirun -np 2 --mca btl sm ./tests/t_dslash_3d  -geom 1 2 1 1
mpirun -np 2 --mca btl sm ./tests/t_dslash_3d  -geom 1 1 2 1
mpirun -np 1 --mca btl sm ./tests/t_dslash_3d  -geom 1 1 1 1

mv ../sse_wilson_dslash/tests/testvol.h.bak ../sse_wilson_dslash/tests/testvol.h
