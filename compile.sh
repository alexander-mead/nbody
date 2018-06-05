#!/bin/bash

#Set the compiler and code
compiler=gfortran
code=nbody

#Normal compiler options
normal='-Warray-bounds -ffree-line-length-none -fmax-errors=4 -ffpe-trap=invalid,zero,overflow -fimplicit-none -O3'

#Set the precision (I think you need REAL*16)
#precision=''
#precision='-fdefault-real-8'
precision='-fdefault-real-16'

$compiler src/constants.f90 src/string_operations.f90 src/file_info.f90 src/vectors.f90 src/$code.f90 $normal $precision -o $code.e
#$compiler src/*.f90 $normal $precision -o $code.e #Does not work because order matters
