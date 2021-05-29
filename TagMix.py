#  © Shahram Talei @ 2021 The University of Alabama - All rights reserved.
#you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.


from __future__ import division
import yt
import numpy as np
from yt.analysis_modules.halo_finding.api import *
from yt.analysis_modules.halo_analysis.api import *
from os import environ
environ['CFLAGS'] = "-I"+np.get_include()

import pyximport; pyximport.install()
#import particle_ops
import argparse


import tempfile
import shutil
import os
import sys

from scipy.spatial.transform import Rotation as R
from numpy import linalg as LA
from operator import mul
from functools import reduce
import matplotlib.pyplot as plt

import csv

plt.rcParams["font.size"] =12

    #how to run: python ShapeAnalysis.py snapshot_file halo_catalog particles_list check_contamination extract_shape bin_number iteration
    #example: $python ShapeAnalysisV2.py snap_264 halos_0.0.ascii halos_0.0.particles 1.4e12 1.1e12   1 1 5 3
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("snap", type=str)
    parser.add_argument("halo",type=str)
    parser.add_argument("GalFile",type=str)
    parser.add_argument("fraction",type=float)
    args = parser.parse_args()
    Gals=np.genfromtxt(args.GalFile, delimiter = ',')
    Gx0=np.array(Gals[:,0])
    Gy0=np.array(Gals[:,1])
    Gz0=np.array(Gals[:,2])
    GMv0=np.array(Gals[:,3])
    GRv0=np.array(Gals[:,4])
    GS0=np.array(Gals[:,5])
    GSM0=np.array(Gals[:,6])
    GIndex0=np.array(Gals[:,7])
    HIndex0=np.array(Gals[:,8])
    CentralG0=np.array(Gals[:,9])
    MetalStellarG0=np.array(Gals[:,10])
    #
    snap = yt.load(args.snap)
    ad = snap.all_data()
    coordinatesDM = ad[("Halo","Coordinates")]
    velocitiesDM = ad[("Halo","Velocities")]
    p=np.array(coordinatesDM)
    v=np.array(velocitiesDM)
    #print(p[:,1].shape)
    #print(p[1:])
    px=p[:,0]
    py=p[:,1]
    pz=p[:,2]
    pVx=v[:,0]
    pVy=v[:,1]
    pVz=v[:,2]
