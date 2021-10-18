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
        #a_group_key = list(f.keys())[0]
        #age =np.array(list(f[a_group_key]))
        #a_group_key = list(f.keys())[1]
        #GID = np.array(list(f[a_group_key]))
        #a_group_key = list(f.keys())[2]
        #HID = np.array(list(f[a_group_key]))
        #a_group_key = list(f.keys())[3]
        #ID =np.array(list(f[a_group_key]))
        #a_group_key = list(f.keys())[4]
        #Metallicity=np.array(list(f[a_group_key]))
        #a_group_key = list(f.keys())[5]
        #StellarMass =np.array(list(f[a_group_key]))
        #a_group_key = list(f.keys())[6]
        #Vx =np.array(list(f[a_group_key]))
        #a_group_key = list(f.keys())[7]
        #Vy =np.array(list(f[a_group_key]))
        #a_group_key = list(f.keys())[8]
        #Vz=np.array(list(f[a_group_key]))
        #a_group_key = list(f.keys())[9]
        #x =np.array(list(f[a_group_key]))
        #a_group_key = list(f.keys())[10]
        #y =np.array(list(f[a_group_key]))
        #a_group_key = list(f.keys())[11]
        #z =np.array(list(f[a_group_key]))
        #['Age', 'GID', 'HID', 'ID', 'Metallicity', 'StellarMass', 'Vx', 'Vy', 'Vz', 'X', 'Y', 'Z']>
        print("Read keys")
        f_key=list(f.keys())
        a_group_key = f_key[0]
        age =np.array(f[a_group_key])
        print("0")
        a_group_key = f_key[1]
        GID = np.array(f[a_group_key])
        print("1")
        a_group_key = f_key[2]
        HID = np.array(f[a_group_key])
        print("2")
        a_group_key = f_key[3]
        ID =np.array(f[a_group_key])
        print("3")
        a_group_key = f_key[4]
        Metallicity=np.array(f[a_group_key])
        print("4")
        a_group_key = f_key[5]
        StellarMass =np.array(f[a_group_key])
        print("5")
        a_group_key = f_key[6]
        Vx =np.array(f[a_group_key])
        print("6")
        a_group_key = f_key[7]
        Vy =np.array(f[a_group_key])
        print("7")
        a_group_key = f_key[8]
        Vz=np.array(f[a_group_key])
        print("8")
        a_group_key = f_key[9]
        x =np.array(f[a_group_key])
        print("9")
        a_group_key = f_key[10]
        y =np.array(f[a_group_key])
        print("10")
        a_group_key = f_key[11]
        z =np.array(f[a_group_key])
        print("11")
        snap = yt.load(args.SnapFile)
        ad = snap.all_data()
        coordinatesDM0 = ad[("Halo","Coordinates")]
        coordinatesDM=np.array(coordinatesDM0.v)
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
        pID = pID0[rp<Rv]
        px=px0[rp<Rv]
        py=py0[rp<Rv]
        pz=pz0[rp<Rv]
        S=len(pID)
        print("Found %d paticles within Rv"%S)
        newX=[0.0]*S
        newY=[0.0]*S
        newZ=[0.0]*S
        newID=[0]*S
        newVx=[0.0]*S
        newVy=[0.0]*S
        newVz=[0.0]*S
        newGID=[0.0]*S
        newHID=[0.0]*S
        newStellarMass=[0.0]*S
        newMetallicity=[0.0]*S
        newAge=[0.0]*S
        #i=0
        for i in range(0,S):
            id =pID[i]
            if(len(Vx[ID==id])>0):
                xx=float(px[pID==id])
                yy=float(py[pID==id])
                zz=float(pz[pID==id])
                newX[i]=xx
                newY[i]=yy
                newZ[i]=zz
                newID[i]=id
                #print("len ID=%d len vx[id=ID]=%d"%(len(ID),len(ID==id)))
                #print("id =%d and ID:"%id)
                #print(Vx[ID==id])
                newVx[i]=float((Vx[ID==id])[0])
                newVy[i]=float((Vy[ID==id])[0])
                newVz[i]=float((Vz[ID==id])[0])
                newGID[i]=int((GID[ID==id])[0])
                newHID[i]=int((HID[ID==id])[0])
                newStellarMass[i]=float((StellarMass[ID==id])[0])
                newMetallicity[i]=float((Metallicity[ID==id])[0])
                newAge[i]=float((age[ID==id])[0])
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
