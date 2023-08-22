#! /bin/bash
#$ -S /bin/bash
#$ -R yes
#$ -V
#$ -cwd
#$ -N Clustering
#$ -j y
#$ -q cpu_long
#$ -pe openmpi_octa 8
#$ -P kenprj

module load amber/2020
module load openmpi/gcc/64/1.10.7
export OMP_NUM_THREADS=8

for r in 6
do

# create clustering directory
rm -r Cluster_$r
mkdir Cluster_$r
cd Cluster_$r
/data01/software/amber20/bin/cpptraj ../TEMPLATES/model.top << EOF >&ptraj_err
trajin  ../trajectory.00.dcd 1 1 
go
EOF

Frames=`grep reading ptraj_err|awk '{print $NF}'|sed 's/)//'`
#Let's use half of the trajectory for equilibration... read the last half of the trajectory
start=$(echo "scale=0;$Frames/2"|bc -l)

# cluster
cat << EOF > cpptraj.in
trajin ../trajectory.00.dcd $start $Frames
rms RMSD first @CA,CB out trajrmsd.dat
cluster hieragglo epsilon $r linkage rms @CA,CB sieve 10 summary summary singlerepout representative repout unique repfmt pdb clusterout clusttraj avgout avg avgfmt pdb
go
EOF

for j in 0 1 2 3 4 5 6 7 8 9 
do
   cat << EOF >> cpptraj.in
   clear trajin
   trajin  clusttraj.c$j
   trajout clusttraj.c$j.pdb model
   average average.c$j.pdb pdb
   reference unique.c$j.pdb [uniq$j]
   rms centr$j ref [uniq$j] @CA 
   atomicfluct out back.$j.apf @C,CA,N,O byres
   go
EOF
done

/data01/software/amber20/bin/cpptraj.MPI -p ../TEMPLATES/model.top -i cpptraj.in
cd ..

done

