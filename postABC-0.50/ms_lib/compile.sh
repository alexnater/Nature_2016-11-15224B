#!/bin/bash

gcc -O3 -c ms.c streec.c rand1.c
ar -rcs libms.a ms.o streec.o rand1.o
rm ms.o streec.o rand1.o

