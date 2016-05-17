#!/bin/bash

#if [ ! -e lib/libknapsack.a ]; then
#  echo 'A lib/libknapsack.a file should be on the current folder. Aborting.'
#  exit
#fi

if [ -e 'mtu1.out' ]; then
  echo 'A mtu1.out file already exists. It will be overwritten.'
fi

gfortran -O3 -c 'hbm_mtu1.f'
if [ $? -ne 0 ]; then
  echo "Errors compiling hbm_mtu1.f"
  exit
fi

gfortran -O3 'hbm_mtu1.o' # -Llib -lknapsack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hbm_mtu1.o"
  exit
fi
rm hbm_mtu1.o

mv a.out 'mtu1.out'
./mtu1.out > mtu1.txt
if [ $? -ne 0 ]; then
  echo "Errors running mtu1"
  exit
fi

echo "Test results written to mtu1.txt."
