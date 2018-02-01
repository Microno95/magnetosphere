#!/bin/bash

#PBS -q pqplasma2
#PBS -l walltime=60:00:00
#PBS -l select=1:ncpus=12:mem=8192mb
#PBS -j oe

cd $PBS_O_WORKDIR

module load anaconda3/personal

python cx1_python.py
