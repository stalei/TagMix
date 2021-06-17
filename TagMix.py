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
import h5py as h5

plt.rcParams["font.size"] =12

def partition(array, begin, end):
    pivot = begin
    for i in range(begin+1, end+1):
        if array[i] <= array[begin]:
            pivot += 1
            array[i], array[pivot] = array[pivot], array[i]
    array[pivot], array[begin] = array[begin], array[pivot]
    return pivot



def quicksort(array, begin=0, end=None):
    if end is None:
        end = len(array) - 1
    def _quicksort(array, begin, end):
        if begin >= end:
            return
        pivot = partition(array, begin, end)
        _quicksort(array, begin, pivot-1)
        _quicksort(array, pivot+1, end)
    return _quicksort(array, begin, end)

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
    IDDM = ad[("Halo","ParticleIDs")]
    #BEDM = ad[("Halo","BindingEnergy")]
    print(sorted(snap.field_list))
    p=np.array(coordinatesDM)
    v=np.array(velocitiesDM)
    Id=np.array(IDDM)
    #print(p[:,1].shape)
    #print(p[1:])
    px=p[:,0]
    py=p[:,1]
    pz=p[:,2]
    pVx=v[:,0]
    pVy=v[:,1]
    pVz=v[:,2]
    print("Total number of galaxies:%d"%len(Gx0))
    for i in range(0,len(Gx0)):
        if True:
            print(i)
        #if GSM0[i] !=0:
            dx=px-Gx0[i]
            dy=py-Gy0[i]
            dz=pz-Gz0[i]
            r2=dx*dx+dy*dy+dz*dz
            r=np.sqrt(r2)
            Pcount=len(r[r<GRv0[i]])
            tagLimit=int((args.fraction/100.)*Pcount)
            rLim=GRv0[i]
            pxh=px[r<rLim]
            pyh=py[r<rLim]
            pzh=pz[r<rLim]
            pVxh=pVx[r<rLim]
            pVyh=pVy[r<rLim]
            pVzh=pVz[r<rLim]
            Idh=Id[r<rLim]
            size=len(pxh)
            PotE=[0.0]*size
            KinE=[0.0]*size
            BE=[0.0]*size
            rh=r[r<rLim]
            c=0
            for j in Idh:
                dxp=pxh[Idh==j]-pxh[Idh !=j]
                dyp=pyh[Idh==j]-pyh[Idh !=j]
                dzp=pzh[Idh==j]-pzh[Idh !=j]
                rp2=dxp*dxp+dyp*dyp+dzp*dzp
                rp=np.sqrt(rp2)
                PotE[c]=np.sum(1./rp)
                KinE[c]=0.5*(pVxh*pVxh+pVyh*pVyh+pVzh*pVzh)
                BE[c]=PotE[c]+KinE[c]
                c+=1
            BE2=np.array(np.sort(BE))
            #quicksort(BE2)
            BErev=BE2[::-1] #reverse it
            BELimit=BErev[tagLimit] #what is there are amny particles at the same BE?
            pxtag=pxh[BE<BELimit]
            pytag=pyh[BE<BELimit]
            pztag=pzh[BE<BELimit]
            pSM=[0.0]*len(pxtag)
            pZZ=[0.0]*len(pxtag)
            pAge=[0.0]*len(pxtag)
            for k in range(0,len(pxtag)):
                pSM[k]=GSM0[i]/tagLimit
                pZZ[k]= MetalStellarG0[i]/GSM0[i]
                pAge[k]=1
            #AllStars[id].ZZ=SageOutput[galaxy].MetalsStellarMass/SageOutput[galaxy].StellarMass
            hf = h5.File('%s.h5' %str(Id[i]), 'w')
            hf.create_dataset('X', data=pxtag)
            hf.create_dataset('Y', data=pytag)
            hf.create_dataset('Z', data=pztag)
            hf.close()
