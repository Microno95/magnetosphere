#!/bin/sh

module load anaconda3/personal

f2py -c Streamtracer.f90 --f90flags='-fopenmp' -lgomp -m streamtracer
f2py -c null_finder.f90 --f90flags='-fopenmp' -lgomp -m null_finder


#f2py -c Streamtracer.f90 -m streamtracer
#f2py -c null_finder.f90 -m null_finder

