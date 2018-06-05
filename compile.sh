#!/bin/bash

#Set the compiler and code
compiler=gfortran
code=nbody

#Normal compiler options
normal='-Warray-bounds -ffree-line-length-none -fmax-errors=4 -ffpe-trap=invalid,zero,overflow -fimplicit-none -O3'

#Mead library location
mead='/Users/Mead/Physics/library'

#Set the precision (I think you need REAL*16)
#precision=''
#precision='-fdefault-real-8'
precision='-fdefault-real-16'

$compiler $mead/constants.f90 $mead/string_operations.f90 $mead/file_info.f90 $mead/vectors.f90 $code.f90 $normal $precision -o $code.e
#$compiler src/*.f90 $normal $precision -o $code.e #Does not work because order matters
