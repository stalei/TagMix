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
from math import *

plt.rcParams["font.size"] =12
def Snap2Z(n):
    csvfile='../SnapZ.csv'
    d=np.genfromtxt(csvfile, delimiter = ',',comments='#')
    snapshot=np.array(d[:,0])
    redshift=np.array(d[:,1])
    zzz=redshift[snapshot==n]
    return zzz
def GetAge(s): # print '''Cosmology calculator ala Ned Wright (www.astro.ucla.edu/~wright)
    z=Snap2Z(s) #float(sys.argv[1+verbose])            # redshift
    H0 = 70.2                         # Hubble constant
    WM = 0.272                        # Omega(matter)
    WV = 1.0 - WM - 0.4165/(H0*H0)
    # initialize constants
    WR = 0.        # Omega(radiation)
    WK = 0.        # Omega curvaturve = 1-Omega(total)
    c = 299792.458 # velocity of light in km/sec
    Tyr = 977.8    # coefficent for converting 1/H into Gyr
    DTT = 0.5      # time from z to now in units of 1/H0
    DTT_Gyr = 0.0  # value of DTT in Gyr
    age = 0.5      # age of Universe in units of 1/H0
    age_Gyr = 0.0  # value of age in Gyr
    zage = 0.1     # age of Universe at redshift z in units of 1/H0
    zage_Gyr = 0.0 # value of zage in Gyr
    DCMR = 0.0     # comoving radial distance in units of c/H0
    DCMR_Mpc = 0.0 
    DCMR_Gyr = 0.0
    DA = 0.0       # angular size distance
    DA_Mpc = 0.0
    DA_Gyr = 0.0
    kpc_DA = 0.0
    DL = 0.0       # luminosity distance
    DL_Mpc = 0.0
    DL_Gyr = 0.0   # DL in units of billions of light years
    V_Gpc = 0.0
    a = 1.0        # 1/(1+z), the scale factor of the Universe
    az = 0.5       # 1/(1+z(object))
    h = H0/100.
    WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
    WK = 1-WM-WR-WV
    az = 1.0/(1+1.0*z)
    age = 0.
    n=1000         # number of points in integrals
    for i in range(n):
        a = az*(i+0.5)/n
        adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        age = age + 1./adot
    zage = az*age/n
    zage_Gyr = (Tyr/H0)*zage
    DTT = 0.0
    DCMR = 0.0
    # do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
    for i in range(n):
        a = az+(1-az)*(i+0.5)/n
        adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        DTT = DTT + 1./adot
        DCMR = DCMR + 1./(a*adot)
    DTT = (1.-az)*DTT/n
    DCMR = (1.-az)*DCMR/n
    age = DTT+zage
    age_Gyr = age*(Tyr/H0)
    DTT_Gyr = (Tyr/H0)*DTT
    DCMR_Gyr = (Tyr/H0)*DCMR
    DCMR_Mpc = (c/H0)*DCMR
    # tangential comoving distance
    ratio = 1.00
    x = sqrt(abs(WK))*DCMR
    if x > 0.1:
        if WK > 0:
            ratio =  0.5*(exp(x)-exp(-x))/x 
        else:
            ratio = sin(x)/x
    else:
        y = x*x
    if WK < 0: y = -y
    ratio = 1. + y/6. + y*y/120.
    DCMT = ratio*DCMR
    DA = az*DCMT
    DA_Mpc = (c/H0)*DA
    kpc_DA = DA_Mpc/206.264806
    DA_Gyr = (Tyr/H0)*DA
    DL = DA/(az*az)
    DL_Mpc = (c/H0)*DL
    DL_Gyr = (Tyr/H0)*DL
# comoving volume computation
    ratio = 1.00
    x = sqrt(abs(WK))*DCMR
    if x > 0.1:
        if WK > 0:
            ratio = (0.125*(exp(2.*x)-exp(-2.*x))-x/2.)/(x*x*x/3.)
        else:
            ratio = (x/2. - sin(2.*x)/4.)/(x*x*x/3.)
    else:
        y = x*x
    if WK < 0: y = -y
    ratio = 1. + y/5. + (2./105.)*y*y
    VCM = ratio*DCMR*DCMR*DCMR/3.
    V_Gpc = 4.*pi*((0.001*c/H0)**3)*VCM
    #if verbose == 1:
    print('For H_o = ' + '%1.1f' % H0 + ', Omega_M = ' + '%1.2f' % WM + ', Omega_vac = ',)
    print('%1.2f' % WV + ', z = ' + '%1.3f' % z)
    print('It is now ' + '%1.1f' % age_Gyr + ' Gyr since the Big Bang.')
    print('The age at redshift z was ' + '%1.1f' % zage_Gyr + ' Gyr.')
    print('The light travel time was ' + '%1.1f' % DTT_Gyr + ' Gyr.')
    print('The comoving radial distance, which goes into Hubbles law, is',)
    print('%1.1f' % DCMR_Mpc + ' Mpc or ' + '%1.1f' % DCMR_Gyr + ' Gly.')
    print('The comoving volume within redshift z is ' + '%1.1f' % V_Gpc + ' Gpc^3.')
    print('The angular size distance D_A is ' + '%1.1f' % DA_Mpc + ' Mpc or',)
    print('%1.1f' % DA_Gyr + ' Gly.')
    print('This gives a scale of ' + '%.2f' % kpc_DA + ' kpc/".')
    print('The luminosity distance D_L is ' + '%1.1f' % DL_Mpc + ' Mpc or ' + '%1.1f' % DL_Gyr + ' Gly.')
    print('The distance modulus, m-M, is '+'%1.2f' % (5*log10(DL_Mpc*1e6)-5))
    #else:
    #    print('%1.2f' % zage_Gyr,)
    #    print('%1.2f' % DCMR_Mpc,)
    #    print('%1.2f' % kpc_DA,)
    #    print('%1.2f' % (5*log10(DL_Mpc*1e6)-5))
    return DTT_Gyr

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

    #how to run: python ShapeAnalysis.py snapshot_num f_mb
    #example: $python ShapeAnalysisV2.py 264 10
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("snapNum", type=int)
    #parser.add_argument("halo",type=str)
    #parser.add_argument("GalFile",type=str)
    parser.add_argument("fraction",type=float)
    #parser.add_argument("Age",type=float)
    args = parser.parse_args()
    s=args.snapNum
    GalFile='../'+str(s)+'.csv'
    Gals=np.genfromtxt(GalFile, delimiter = ',')
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
    AgeGyr=GetAge(s)
    if s>100:
        snap='../snap_'+str(s)
    else:
        snap='../snap_0'+str(s)
    snap = yt.load(snap)
    ad = snap.all_data()
    coordinatesDM = ad[("Halo","Coordinates")]
    velocitiesDM = ad[("Halo","Velocities")]
    IDDM = ad[("Halo","ParticleIDs")]
    #BEDM = ad[("Halo","BindingEnergy")]
    #print(sorted(snap.field_list))
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
    #hf = h5.File('s1.h5', 'w')
    for i in range(1,len(Gx0)):
        if True:
            print(i)
        #if GSM0[i] !=0:
            dx=px-Gx0[i]
            dy=py-Gy0[i]
            dz=pz-Gz0[i]
            r2=dx*dx+dy*dy+dz*dz
            r=np.sqrt(r2)
            Pcount=len(r[r<=GRv0[i]])
            tagLimit=int((args.fraction/100.)*Pcount)
            if tagLimit>Pcount:
                tagLimit=Pcount
            if Pcount==0:
                continue
            #    tagLimit=1
            rLim=GRv0[i]
            pxh=px[r<=rLim]
            pyh=py[r<=rLim]
            pzh=pz[r<=rLim]
            pVxh=pVx[r<=rLim]
            pVyh=pVy[r<=rLim]
            pVzh=pVz[r<=rLim]
            Idh=Id[r<=rLim]
            size=len(pxh)
            print("Ps in Rv:%d"%size)
            #PotE=[0.0]*size
            #KinE=[0.0]*size
            BE=[0.0]*size
            print(np.array(BE).shape)
            rh=r[r<=rLim]
            c=0
            for j in Idh:
                dxp=pxh[Idh==j]-pxh[Idh !=j]
                dyp=pyh[Idh==j]-pyh[Idh !=j]
                dzp=pzh[Idh==j]-pzh[Idh !=j]
                vx=float(pVxh[Idh==j])
                vy=float(pVyh[Idh==j])
                vz=float(pVzh[Idh==j])
                #print("vx:%g"%vx)
                #rp2=dxp*dxp+dyp*dyp+dzp*dzp
                rp=np.sqrt(dxp*dxp+dyp*dyp+dzp*dzp)
                #PotE[c]=np.sum(1./rp)
                #KinE[c]=0.5*(pVxh*pVxh+pVyh*pVyh+pVzh*pVzh)
                BE[c]=float(np.sum(1./rp)+0.5*(vx*vx+vy*vy+vz*vz))#PotE[c]+KinE[c]
                c+=1
            print("counted:%d"%c)
            BE2=BE#np.array(np.sort(BE))
            #BE.sort(key=lambda x: x[0],reverse=True)
            BE2.sort(reverse=True)
            #print(BE.shape)
            #print("before sort:")
            print(np.array(BE))
            print("after sort:")
            print(np.array(BE2))
            #quicksort(BE2)
            #BErev=BE2[::-1] #reverse it
            #print(BE)
            #BELimit=BE[0][tagLimit] #what is there are amny particles at the same BE?
            BELimit=BE2[tagLimit]
            print("BELimit:%g"%BELimit)
            #print(BE[0][:])
            BE=np.array(BE)
            pxtag=pxh[BE>=BELimit]
            pytag=pyh[BE>=BELimit]
            pztag=pzh[BE>=BELimit]
            pIDtag=Idh[BE>=BELimit]
            #
            pVxTag=pVxh[BE>=BELimit]
            pVyTag=pVyh[BE>=BELimit]
            pVzTag=pVzh[BE>=BELimit]
            pGID=[GIndex0[i]]*len(Idh)
            pHID=[HIndex0[i]]*len(Idh)
            #
            print(" # of most bound Ps:%d"%len(pxtag))
            pSM=[0.0]*len(pxtag)
            pZZ=[0.0]*len(pxtag)
            pAge=[0.0]*len(pxtag)
            for k in range(0,len(pxtag)):
                pSM[k]=GSM0[i]/tagLimit
                pZZ[k]= MetalStellarG0[i]/GSM0[i]
                pAge[k]=AgeGyr
            #AllStars[id].ZZ=SageOutput[galaxy].MetalsStellarMass/SageOutput[galaxy].StellarMass
            hf = h5.File('%s.h5' %str(Id[i]), 'w')
            hf.create_dataset('ID', data=pIDtag)
            hf.create_dataset('X', data=pxtag)
            hf.create_dataset('Y', data=pytag)
            hf.create_dataset('Z', data=pztag)
            #
            hf.create_dataset('Vx', data=pVxTag)
            hf.create_dataset('Vy', data=pVyTag)
            hf.create_dataset('Vz', data=pVzTag)
            hf.create_dataset('GID', data=pGID)
            hf.create_dataset('HID', data=pHID)
            #
            hf.create_dataset('StellarMass', data=pSM)
            hf.create_dataset('Metallicity', data=pZZ)
            hf.create_dataset('Age', data=pAge)
            hf.close()
