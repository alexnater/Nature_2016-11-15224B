#!/bin/bash

gcc -O3 -c randlib.c linpack.c com.c
ar -rcs librand.a randlib.o linpack.o com.o
rm randlib.o linpack.o com.o

