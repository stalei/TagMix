#  © Shahram Talei @ 2021 The University of Alabama - All rights reserved.
#you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import division
import h5py as h5
import numpy as np
from os import environ
import os
import yt
environ['CFLAGS'] = "-I"+np.get_include()
import argparse
import glob

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("TagFile", type=str)
    parser.add_argument("SnapFile", type=str)
    #parser.add_argument("Hi", type=int)
    #parser.add_argument("Hf", type=int)
    args = parser.parse_args()
    h5name=args.TagFile
    with h5.File(h5name, "r") as f:
        # List all groups
        print("Keys: %s" % f.keys())
        a_group_key = list(f.keys())[0]
        age = list(f[a_group_key])
        a_group_key = list(f.keys())[1]
        GID = list(f[a_group_key])
        a_group_key = list(f.keys())[2]
        HID = list(f[a_group_key])
        a_group_key = list(f.keys())[3]
        ID = list(f[a_group_key])
        a_group_key = list(f.keys())[4]
        Metallicity= list(f[a_group_key])
        a_group_key = list(f.keys())[5]
        StellarMass = list(f[a_group_key])
        a_group_key = list(f.keys())[6]
        Vx = list(f[a_group_key])
        a_group_key = list(f.keys())[7]
        Vy = list(f[a_group_key])
        a_group_key = list(f.keys())[8]
        Vz= list(f[a_group_key])
        a_group_key = list(f.keys())[9]
        x = list(f[a_group_key])
        a_group_key = list(f.keys())[10]
        y = list(f[a_group_key])
        a_group_key = list(f.keys())[11]
        z = list(f[a_group_key])
        snap = yt.load(args.SnapFile)
        ad = snap.all_data()
        coordinatesDM = ad[("Halo","Coordinates")]
        px0=coordinatesDM[:,0]
        py0=coordinatesDM[:,1]
        pz0=coordinatesDM[:,2]
        pID0 = ad[("Halo","ParticleIDs")]
        #m12i: 29.3575,31.0276,32.4926,6.37801e+11,0.139977,264,5.1374e+10,0,0,0,3.99904e+08
        gx=29.3575
        gy=31.0276
        gz=32.4926
        Rv=0.139977
        dpx=px0-gx
        dpy=py0-gy
        dpz=pz0-gz
        rp=(dpx*dpx+dpy*dpy+dpz*dpz)**0.5
        pID = PID0[rp<Rv]
        px=px0[rp<Rv]
        py=py0[rp<Rv]
        pz=pz0[rp<Rv]
        S=len(pID)
        newX=[0.0]*S
        newY=[0.0]*S
        newZ=[0.0]*S
        newID[i]=[0]*S
        newVx[i]=[0.0]*S
        newVy[i]=[0.0]*S
        newVz[i]=[0.0]*S
        newGID[i]=[0.0]*S
        newHID[i]=[0.0]*S
        newStellarMass[i]=[0.0]*S
        newMetallicity[i]=[0.0]*S
        newAge[i]=[0.0]*S
        #i=0
        for i in range(0,S):
            id =pID[i]
            xx=px[pID==id]
            yy=py[pID==id]
            zz=pz[pID==id]
            newX[i]=xx
            newY[i]=yy
            newZ[i]=zz
            newID[i]=id
            newVx[i]=Vx[ID==id]
            newVy[i]=Vy[ID==id]
            newVz[i]=Vz[ID==id]
            newGID[i]=GID[ID==id]
            newHID[i]=HID[ID==id]
            newStellarMass[i]=StellarMass[ID==id]
            newMetallicity[i]=Metallicity[ID==id]
            newAge[i]=age[ID==id]
            #i+=1
        hf = h5.File('AllTagsPosFixed2.h5', 'w')
        hf.create_dataset('ID', data=newID)
        hf.create_dataset('X', data=newX)
        hf.create_dataset('Y', data=newY)
        hf.create_dataset('Z', data=newZ)
        #
        hf.create_dataset('Vx', data=newVx)
        hf.create_dataset('Vy', data=newVy)
        hf.create_dataset('Vz', data=newVz)
        hf.create_dataset('GID', data=newGID)
        hf.create_dataset('HID', data=newHID)
        #
        hf.create_dataset('StellarMass', data=newStellarMass)
        hf.create_dataset('Metallicity', data=newMetallicity)
        hf.create_dataset('Age', data=newAge)
        hf.close()
