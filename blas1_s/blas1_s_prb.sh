#!/bin/bash
#
gfortran -c blas1_s_prb.f
if [ $? -ne 0 ]; then
  echo "Errors compiling blas1_s_prb.f"
  exit
fi
#
gfortran blas1_s_prb.o -L$HOME/libf77 -lblas
if [ $? -ne 0 ]; then
  echo "Errors linking and loading blas1_s_prb.o"
  exit
fi
rm blas1_s_prb.o
#
mv a.out blas1_s_prb
./blas1_s_prb > blas1_s_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running blas1_s_prb"
  exit
fi
rm blas1_s_prb
#
echo "Test results written to blas1_s_prb_output.txt."
