#!/bin/bash
#
mkdir temp
cd temp
#
~/binc/f77split ../../blas0/blas0.f
#
~/binc/f77split ../../blas1_c/blas1_c.f
~/binc/f77split ../../blas1_d/blas1_d.f
~/binc/f77split ../../blas1_s/blas1_s.f
~/binc/f77split ../../blas1_z/blas1_z.f
#
~/binc/f77split ../../blas2_c/blas2_c.f
~/binc/f77split ../../blas2_d/blas2_d.f
~/binc/f77split ../../blas2_s/blas2_s.f
~/binc/f77split ../../blas2_z/blas2_z.f
#
~/binc/f77split ../../blas3_c/blas3_c.f
~/binc/f77split ../../blas3_d/blas3_d.f
~/binc/f77split ../../blas3_s/blas3_s.f
~/binc/f77split ../../blas3_z/blas3_z.f
#
for FILE in `ls -1 *.f`;
do
  gfortran -c $FILE
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
done
rm *.f
#
#  Create the library.
#
ar qc libblas.a *.o
rm *.o
#
mv libblas.a ~/libf77
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/libblas.a."
