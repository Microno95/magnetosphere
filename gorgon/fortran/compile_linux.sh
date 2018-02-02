#!/bin/sh

#module load anaconda3/personal

f2py -c connectivity.f90 --f90flags='-fopenmp' -lgomp -m connectivity_tracer
f2py -c Streamtracer.f90 --f90flags='-fopenmp' -lgomp -m streamtracer
f2py -c null_finder.f90 --f90flags='-fopenmp' -lgomp -m null_finder