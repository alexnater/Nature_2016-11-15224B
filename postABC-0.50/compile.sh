#!/bin/bash

cd ms_lib
gcc -O3 -c ms.c streec.c rand1.c
ar -rcs libms.a ms.o streec.o rand1.o
rm ms.o streec.o rand1.o

cd ../randlib1.3
gcc -O3 -c randlib.c linpack.c com.c
ar -rcs librand.a randlib.o linpack.o com.o
rm randlib.o linpack.o com.o

cd ..
g++ -O3 *.cpp -L ./randlib1.3 -L ./ms_lib -o postABC -lrand -lms
rm ./ms_lib/libms.a ./randlib1.3/librand.a

