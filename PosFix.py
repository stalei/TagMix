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
        snap = yt.load(args.snap)
        ad = snap.all_data()
        coordinatesDM = ad[("Halo","Coordinates")]
        px=p[:,0]
        py=p[:,1]
        pz=p[:,2]
        pID = ad[("Halo","ParticleIDs")]
        S=len(ID)
        newX=[0.0]*S
        newY=[0.0]*S
        newZ=[0.0]*S
        #i=0
        for i in range(0,S):
            id =ID[i]
            newX[i]=px[pID==id]
            newY[i]=py[pID==id]
            newZ[i]=pz[pID==id]
            #i+=1
        hf = h5.File('AllTagsPosFixed.h5', 'w')
        hf.create_dataset('ID', data=ID)
        hf.create_dataset('X', data=newX)
        hf.create_dataset('Y', data=newY)
        hf.create_dataset('Z', data=newZ)
        #
        hf.create_dataset('Vx', data=Vx)
        hf.create_dataset('Vy', data=Vy)
        hf.create_dataset('Vz', data=Vz)
        hf.create_dataset('GID', data=GID)
        hf.create_dataset('HID', data=HID)
        #
        hf.create_dataset('StellarMass', data=StellarMass)
        hf.create_dataset('Metallicity', data=Metallicity)
        hf.create_dataset('Age', data=age)
        hf.close()
