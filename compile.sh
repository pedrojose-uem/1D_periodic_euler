#!/bin/bash
set -e

SOURCE_FILES="\
main.cpp \
DataStructs.cpp \
rk4.cpp \
Euler1D.cpp"

g++ -g -O3 -D_DOUBLE_ $SOURCE_FILES -Iincludes -o euler_double.p

g++ -g -O3 $SOURCE_FILES -Iincludes -o euler_single.p
