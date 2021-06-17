#!/bin/bash

#SBATCH -n 1

#SBATCH --mem-per-cpu 25000

#SBATCH -o  /bighome/stalei/analysis/h32251z5/shape0.62C.log

#SBATCH -p main

#SBATCH --qos main
#SBATCH --mail-user=stalei@crimson.ua.edu
#SBATCH --mail-type=ALL

cd /bighome/stalei/analysis #/h32251z5/

#how to run: python ShapeAnalysis.py snapshot_file halo_catalog particles_list check_contamination extract_shape bin_number iteration
#example: $python ShapeAnalysisV2.py snap_264 halos_0.0.ascii halos_0.0.particles 1.4e12 1.1e12   1 1 5 3

#python ShapeAnalysisV2.py snap_g5z32251_264 halos_g32251z5.ascii halos_g32251z5.particles 1.4e12 9.0e11   0 1 8 3
#python ShapeAnalysisV2.py z1.5/snap_g5z32251_184 z1.5/halos_g5z32251_184.ascii z1.5/halos_g5z32251_184.particles 1.4e12 9.0e10   0 1 15 3

#python ShapeAnalysisV2.py z0.62/C32251snap_027 z0.62/halos_c32251z5.ascii z0.62/halos_c32251z5.particles 1.4e12 5.0e11   0 1 30 3

#python ShapeAnalysisV2.py m12f/z0/snap_m12f_Gz0 m12f/z0/halos_m12fG.ascii m12f/z0/halos_m12fG.particles 1.4e12 5.0e11   0 1 30 3

python ShapeAnalysisV2.py h32251z5/z0.18/snap_g32251z5_z0.18 h32251z5/z0.18/halos_g32251z5_z0.18.ascii h32251z5/z0.18/halos_g32251z5_z0.18.particles 1.4e12 5.0e11   0 1 15 3

python ShapeAnalysisV2.py h32251z5/z0/snap_g5z32251_264 h32251z5/z0/halos_g32251z5_z0.ascii h32251z5/z0/halos_g32251z5_z0.particles 1.4e12 5.0e11   0 1 30 5
