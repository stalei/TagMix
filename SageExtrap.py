#  © Shahram Talei @ 2021 The University of Alabama - All rights reserved.
#you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division
import numpy as np
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
    #how to run: python SageExtrap.py gals1 gals2
    #example: $python SageExtrap.py gals158.csv gals159.csv
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("GalFile1",type=str)
    parser.add_argument("GalFile2",type=str)
    args = parser.parse_args()
    #Load first galaxy file
    Gals1=np.genfromtxt(args.GalFile1, delimiter = ',')
    Gx1=np.array(Gals1[:,0])
    Gy1=np.array(Gals1[:,1])
    Gz1=np.array(Gals1[:,2])
    GMv1=np.array(Gals1[:,3])
    GRv1=np.array(Gals1[:,4])
    GS1=np.array(Gals1[:,5])
    GSM1=np.array(Gals1[:,6])
    GIndex1=np.array(Gals1[:,7])
    HIndex1=np.array(Gals1[:,8])
    CentralG1=np.array(Gals1[:,9])
    MetalStellarG1=np.array(Gals1[:,10])
    #2nd galaxy file
    Gals2=np.genfromtxt(args.GalFile2, delimiter = ',')
    Gx2=np.array(Gals2[:,0])
    Gy2=np.array(Gals2[:,1])
    Gz2=np.array(Gals2[:,2])
    GMv2=np.array(Gals2[:,3])
    GRv2=np.array(Gals2[:,4])
    GS2=np.array(Gals2[:,5])
    GSM2=np.array(Gals2[:,6])
    GIndex2=np.array(Gals2[:,7])
    HIndex2=np.array(Gals2[:,8])
    CentralG2=np.array(Gals2[:,9])
    MetalStellarG2=np.array(Gals2[:,10])
    for id in GIndex1:
        if len(Gx2[GIndex2==id])>0:
            dx=Gx2[GIndex2==id]-Gx1[GIndex1==id]
            #print("dx:%g"%dx)
            Gx3=2.*(Gx2[GIndex2==id])-Gx1[GIndex1==id]
            Gy3=2.*(Gy2[GIndex2==id])-Gy1[GIndex1==id]
            Gz3=2.*(Gz2[GIndex2==id])-Gz1[GIndex1==id]
            GMv3=np.abs(2.*(GMv2[GIndex2==id])-GMv1[GIndex1==id])
            GRv3=2.*(GRv2[GIndex2==id])-GRv1[GIndex1==id]
            if GRv3 < GRv2[GIndex2==id]/5:
                GRv3=GRv2[GIndex2==id]
            GS3=GS2[GIndex2==id]+1
            GSM3=np.abs(2.*(GSM2[GIndex2==id])-GSM1[GIndex1==id])
            GIndex3=id
            HIndex3=HIndex2[GIndex2==id]
            CentralG3=CentralG2[GIndex2==id]
            MetalStellarG3=2.*(MetalStellarG2[GIndex2==id])-MetalStellarG1[GIndex1==id]
            #print("id:%d"%id)
            #print("G2:")
            #print(GRv2[GIndex2==id])
            #print("G1:")
            #print(GRv1[GIndex1==id])
            #print(GRv2[GIndex2==id]-GRv1[GIndex1==id])
            print("%g,%g,%g,%g,%g,%d,%g,%d,%d,%d,%g"%(Gx3,Gy3,Gz3,GMv3,GRv3,GS3,GSM3,GIndex3,HIndex3,CentralG3,MetalStellarG3))
