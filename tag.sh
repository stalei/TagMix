#!/bin/bash

#SBATCH -n 1

#SBATCH --mem-per-cpu 32000

#SBATCH -o  /bighome/stalei/analysis/TagMix/tag.log

#SBATCH -p main

#SBATCH --qos main
#SBATCH --mail-user=stalei@crimson.ua.edu
#SBATCH --mail-type=ALL

cd /bighome/stalei/analysis/TagMix #/h32251z5/

python TagMix.py /bighome/stalei/analysis/h32251z5/z0.18/snap_c32251z5_z0.18 /bighome/stalei/analysis/h32251z5/z0.18/halos_c32251z5_z0.18.ascii /bighome/stalei/analysis/TagMix/gals-z0.18.csv 10
