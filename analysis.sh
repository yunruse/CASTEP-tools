#!/bin/bash -l

# Wallclock (h:mm:ss)
#$ -l h_rt=8:00:00

# Budget
#$ -P Free
#$ -A UKCP_ED_P

# Request nodes with RAM and TMPDIR
#$ -pe mpi 1
#$ -l mem=1G
#$ -l tmpfs=15G

#$ -N analysis_res
#$ -cwd

mkdir graphs_live

date +%s
for f in cells/*.md; do
  fname=$(echo $f | cut -d/ -f2 | cut -d- -f1,2)
  fcellstat="../lab/aimd_structures/$fname.hydrogens.txt"
  echo $f
  python ~/lab/analysis.py $f \
    'graphs_live/$name-$method.png' \
    --every 20 --hydropath $fcellstat -s \
    $@
done
date +%s
