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
environ['CFLAGS'] = "-I"+np.get_include()
import argparse
import glob

if __name__ == "__main__":
    #parser = argparse.ArgumentParser()
    #parser.add_argument("Folder", type=str)
    #parser.add_argument("Hi", type=int)
    #parser.add_argument("Hf", type=int)
    #args = parser.parse_args()
    #
    # with h5.File('tag_merge.h5',mode='w') as h5fw:
    #     row1 = 0
    #     for h5name in glob.glob('testTag/*.h5'):
    #         print(h5name)
    #         h5fr = h5.File(h5name,'r')
    #         dset1 = list(h5fr.keys())[0]
    #         print(dset1)
    #         arr_data = h5fr[dset1][:]
    #         print(arr_data.shape[0])
    #         dslen = arr_data.shape[0]
    #         cols = arr_data.shape[1]
    #         if row1 == 0:
    #             h5fw.create_dataset('alldata', dtype="f",  shape=(dslen,cols), maxshape=(None, cols) )
    #         if row1+dslen <= len(h5fw['alldata']) :
    #             h5fw['alldata'][row1:row1+dslen,:] = arr_data[:]
    #         else :
    #             h5fw['alldata'].resize( (row1+dslen, cols) )
    #             h5fw['alldata'][row1:row1+dslen,:] = arr_data[:]
    #         row1 += dslen
    # d_names = os.listdir("testTag")#  (os.getcwd())
    # print(d_names)
    #glob   'testTag/*.h5'
    address='*.h5'
    AllID=[]
    c=0
    for h5name in glob.glob(address):
        print(h5name)
        with h5.File(h5name, "r") as f:
            # List all groups
            print("Keys: %s" % f.keys())
            a_group_key = list(f.keys())[2]
            # Get the data
            data = list(f[a_group_key])
            c+=len(data)
            #AllID.append(np.array(data))
            #print(data)
            f.close()
    #print(np.array(AllID).shape)
    print("Count:%d"%c)
    AllAge=[0.0]*c
    AllID=[0]*c
    AllHID=[0]*c
    AllGID=[0]*c
    AllMetallicity=[0.0]*c
    AllStellarMass=[0.0]*c
    AllVx=[0.0]*c
    AllVy=[0.0]*c
    AllVz=[0.0]*c
    AllX=[0.0]*c
    AllY=[0.0]*c
    AllZ=[0.0]*c
    i=0
    for h5name in glob.glob(address):
        print(h5name)
        with h5.File(h5name, "r") as f:
            # List all groups
            print("Keys: %s" % f.keys())
            #a_group_key = list(f.keys())[0]
            #age = list(f[a_group_key])
            #a_group_key = list(f.keys())[1]
            #GID = list(f[a_group_key])
            #a_group_key = list(f.keys())[2]
            #HID = list(f[a_group_key])
            #a_group_key = list(f.keys())[3]
            #ID = list(f[a_group_key])
            #a_group_key = list(f.keys())[4]
            #Metallicity= list(f[a_group_key])
            #a_group_key = list(f.keys())[5]
            #StellarMass = list(f[a_group_key])
            #a_group_key = list(f.keys())[6]
            #Vx = list(f[a_group_key])
            #a_group_key = list(f.keys())[7]
            #Vy = list(f[a_group_key])
            #a_group_key = list(f.keys())[8]
            #Vz= list(f[a_group_key])
            #a_group_key = list(f.keys())[9]
            #x = list(f[a_group_key])
            #a_group_key = list(f.keys())[10]
            #y = list(f[a_group_key])
            #a_group_key = list(f.keys())[11]
            #z = list(f[a_group_key])
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
            for j in range(0,len(age)):
                AllAge[i]=age[j]
                AllID[i]=ID[j]
                AllGID[i]=GID[j]
                AllHID[i]=HID[j]
                AllMetallicity[i]=Metallicity[j]
                AllStellarMass[i]=StellarMass[j]
                AllVx[i]=Vx[j]
                AllVy[i]=Vy[j]
                AllVz[i]=Vz[j]
                AllX[i]=x[j]
                AllY[i]=y[j]
                AllZ[i]=z[j]
                i+=1
    print(np.array(AllAge).shape)
    #hf = h5.File('%s/AllTags.h5' %args.Folder, 'w')
    hf = h5.File('AllTags.h5', 'w')
    hf.create_dataset('ID', data=AllID)
    hf.create_dataset('X', data=AllX)
    hf.create_dataset('Y', data=AllY)
    hf.create_dataset('Z', data=AllZ)
    #
    hf.create_dataset('Vx', data=AllVx)
    hf.create_dataset('Vy', data=AllVy)
    hf.create_dataset('Vz', data=AllVz)
    hf.create_dataset('GID', data=AllGID)
    hf.create_dataset('HID', data=AllHID)
    #
    hf.create_dataset('StellarMass', data=AllStellarMass)
    hf.create_dataset('Metallicity', data=AllMetallicity)
    hf.create_dataset('Age', data=AllAge)
    hf.close()
    # d_struct = {} #Here we will store the database structure
    # for i in d_names:
    #     f = h5.File(i,'r+')
    #     d_struct[i] = f.keys()
    # f.close()
    # for i in d_names:
    #     for j  in d_struct[i]:
    #         print(j)
    #         #h5copy -i '{i}' -o 'output.h5' -s {j} -d {j}
    #h5fr.close()
