#!/bin/bash
#$ -S /bin/bash
#$ -R yes
#$ -cwd
#$ -V
#$ -N Min
#$ -j y
#$ -pe openmpi 12
#$ -P kenprj
#$ -q cpu_short

pwd

module load amber/2020

mpirun -np 12 /data01/software/amber20/bin/pmemd.MPI  -O  -i min.in  -p fx.top  -c fx.crd  -r minx.rst  -o minx.out

