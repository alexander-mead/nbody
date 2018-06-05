#!/bin/bash

#Libraries used
libs+=('constants')
libs+=('file_info')
libs+=('string_operations')
libs+=('vectors')

#Copy the libraries
for code in "${libs[@]}"; do
    cp /Users/Mead/Physics/library/$code.f90 .
done    
