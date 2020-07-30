#!/bin/bash
#
mkdir temp
cd temp
rm *
~/binc/f77split ../blas1_z.f
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
ar r ~/libf77/libblas.a *.o
rm *.o
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/libblas.a."
