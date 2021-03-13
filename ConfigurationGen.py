# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 08:51:07 2020

@author: mhdho
"""

import numpy as np
import math as MT
import random
import matplotlib.pyplot as plt
import os.path
from mpl_toolkits.mplot3d import Axes3D
import time

FileName="InitCordLammps.dat"

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

a=MT.pi/4
b=0
c=0

def bondLength(ATOM,PLID,SimulationLength,Polen):
    
    for i in range (0,len(ATOM)-1):
        if (not PLID[i] == PLID[i+1]):
            print(i)
            continue
        else:
            x = ATOM[i][0]-ATOM[i+1][0]
            y = ATOM[i][1]-ATOM[i+1][1]
            z = ATOM[i][2]-ATOM[i+1][2]
            
            dist = (x*x+y*y+z*z)**0.5
            if (dist > SimulationLength/2):
                print(str(dist) +" for this value of i " + str (i))


def Rotation (a,b,c):

    S1= [round(np.cos(a)*np.cos(b),4), round(np.cos(a)*np.sin(b)*np.sin(c)-np.sin(a)*np.cos(c),4),
       round( np.cos(a)*np.sin(b)*np.cos(c)+np.sin(a)*np.sin(c),4)]
    
    S2=[round(np.sin(a)*np.cos(b),4), round(np.sin(a)*np.sin(b)*np.sin(c),4) + round(np.cos(a)*np.cos(c),4),
        round(np.sin(a)*np.sin(b)*np.cos(c),4) - round(np.cos(a)*np.sin(c),4)]
    
    S3=[round(-1*np.sin(b),4), round(np.cos(b)*np.sin(c),4), round(np.cos(b)*np.cos(c),4)]
    
    R= [S1,S2,S3]
    return R


def CorCord(X,SimSizeX): # Correct of coordinate according to PBC
    
    #Nx=X
    
    
    if (X   <  0):
#        print("need lower")
#        print(X)
        X = X + SimSizeX
#        print("New Value",str(X))
    if (X > SimSizeX):
#        print("need Higher")
#        print(X)
        X = X - SimSizeX
#        print("New Value",str(X))
    
#    if (X > SimSizeX ):
#        Nx = X - SimSizeX
#    if (X < 0 ):
#        Nx = X + SimSizeX

#    
    return X
        
    
    

Rot = Rotation(a,b,c)
P1=[1,1,1]

ATOM=[]

PolID=[]

NumPoly=[]

BOND= True
ANGLE= True
Dehid = True


TotalNumberOfPolymer=200
PolymerLength = 200
MolarMass= 1 
SimSize= 50

UB = SimSize - SimSize/4
LB =  SimSize/4

FileNameSim = "P-"+str(TotalNumberOfPolymer)+"L-"+str(PolymerLength)+str(".txt")

P2= np.matmul(P1, Rot)
ChainCenter = [1,1,1]

BondLength = 0.8

for i in range (0,TotalNumberOfPolymer):
    ATOM.append([random.uniform(LB, UB),random.uniform(LB, UB),random.uniform(LB, UB)])
    PolID.append(i)
    
    NumPoly.append(i+1)
    j=0
    
    while (j<PolymerLength-1):
    #for j in range (0,PolymerLength-1):
        
        Rot = Rotation(2*MT.pi*random.random(),2*MT.pi*random.random(),2*MT.pi*random.random())
        #Rot = Rotation(1*MT.pi/2,1*MT.pi/2,0)
        
        if j%3 == 2:
            Rot=np.transpose(Rot)
        
        P2= BondLength*np.matmul(ChainCenter, Rot)
        
        Temp=[ATOM[-1][0]+P2[0],ATOM[-1][1]+P2[1],
              ATOM[-1][2]+P2[2]]
        
#        Temp=[CorCord(ATOM[-1][0]+P2[0],SimSize),CorCord(ATOM[-1][1]+P2[1],SimSize),
#              CorCord(ATOM[-1][2]+P2[2],SimSize)]
        if (0 > ATOM[-1][0]+P2[0] or ATOM[-1][0]+P2[0] > SimSize or 0 > ATOM[-1][2]+P2[2] 
        or ATOM[-1][2]+P2[2] > SimSize or 0 > ATOM[-1][1]+P2[1] or ATOM[-1][1]+P2[1] > SimSize):
            print("outside the system maybe unstable ")
            
        
        else :
            print(j)
            j=j+1
            ATOM.append(Temp)
            PolID.append(i)


print(len(ATOM))
f= open(FileName,"w")


f.write("# Model for POlYMER" +2*"\n\r")
#f.write("\r\n")
f.write(5*" "+str(TotalNumberOfPolymer*PolymerLength)+5*" "+"atoms"+"\n\r")
if (BOND):
    
    f.write(5*" "+str(TotalNumberOfPolymer*(PolymerLength-1))+5*" "+"bonds"+"\n\r")
    
if (ANGLE):

    f.write(5*" "+str(TotalNumberOfPolymer*(PolymerLength-2))+5*" "+"angles"+  "\n\r")
    
if (Dehid):
    
    f.write(5*" "+str(TotalNumberOfPolymer*(PolymerLength-3))+5*" "+"dihedrals"+  "\n\r")

f.write(5*" "+str(1)+5*" "+"atom types"+  2*"\n\r")
if (BOND):
    f.write(5*" "+str(1)+5*" "+"bond types"+  2*"\n\r")
if (ANGLE):
    f.write(5*" "+str(1)+5*" "+"angle types"+ 2*"\n\r")
if (Dehid):
    f.write(5*" "+str(1)+5*" "+"dihedral types"+  2*"\n\r")

f.write(5*" "+str(0)+5*" "+str(SimSize)+5*" "+"xlo "+"xhi "+  "\n\r")
f.write(5*" "+str(0)+5*" "+str(SimSize)+5*" "+"ylo "+"yhi "+  "\n\r")
f.write(5*" "+str(0)+5*" "+str(SimSize)+5*" "+"zlo "+"zhi "+  "\n\r")

f.write("\n\r")
f.write("Masses "+  2*"\n\r")


f.write(5*" "+str(1)+5*" "+str(MolarMass)+ 2*"\n\r")

f.write("Atoms"+2*"\n\r")

for i in range (0,len(NumPoly)):
    
    for j in range (0,PolymerLength):
        
        f.write(5*" "+5*" "+str(i*PolymerLength+j+1)+ 5*" " +str(i+1)+ 5*" "+str(1)+5*" ")
        f.write(str(ATOM[i*PolymerLength+j][0])+5*" "+str(ATOM[i*PolymerLength+j][1])+5*" "+str(ATOM[i*PolymerLength+j][2]))
    
        f.write( "\n\r")



if (BOND):
    f.write("Bonds"+2*"\n\r")
    BondNum=1
    
    for i in range (0,len(NumPoly)):
        
        for j in range (0,PolymerLength-1):
            
            
            f.write(5*" "+5*" "+str(BondNum)+ 5*" "+str(1)+5*" ")
            f.write(str(i*PolymerLength+j+1)+ 5*" "+str(i*PolymerLength+j+2))
            
            f.write( "\n\r")
            BondNum=BondNum+1
    

if (ANGLE):
    f.write("Angles"+ 2*"\n\r")
    AngNum=1
    
    
    for i in range (0,len(NumPoly)):
        
        for j in range (0,PolymerLength-2):
    
            
            f.write(5*" "+5*" "+str(AngNum)+ 5*" "+str(1)+5*" ")
            f.write(str(i*PolymerLength+j+1)+ 5*" "+str(i*PolymerLength+j+2)+5*" "+str(i*PolymerLength+j+3))
            
            f.write( "\n\r")
            AngNum=AngNum+1



if (Dehid):
    f.write("Dihedrals"+2*"\n\r")
    DahNum=1
    
    
    for i in range (0,len(NumPoly)):
        
        for j in range (0,PolymerLength-3):
    
            
            f.write(5*" "+5*" "+str(DahNum)+ 5*" "+str(1)+5*" ")
            f.write(str(i*PolymerLength+j+1)+ 5*" "+str(i*PolymerLength+j+2)+
                    5*" "+str(i*PolymerLength+j+3)+5*" "+str(i*PolymerLength+j+4))
            
            f.write( "\n\r")
            DahNum=DahNum+1

f.close()

#f = open(FileName).readlines()
#for s in f:
#    s.rstrip('\n')


#s.rstrip('\n') for s in flist
# replacement strings
WINDOWS_LINE_ENDING = b'\n\r'
UNIX_LINE_ENDING = b'\n'

# relative or absolute file path, e.g.:
file_path = FileName

with open(file_path, 'rb') as open_file:
    content = open_file.read()

content = content.replace(WINDOWS_LINE_ENDING, UNIX_LINE_ENDING)

with open(file_path, 'wb') as open_file:
    open_file.write(content)

bondLength(ATOM,PolID,SimSize,PolymerLength)

X=[]
Y=[]
Z=[]
for i in range (0,len(ATOM)):
    X.append(ATOM[i][0])
    Y.append(ATOM[i][1])
    Z.append(ATOM[i][2])
print("after")

ax.plot(X, Y, Z, marker=".")
#ax.set_xlim(int(min(X))/2,3*int(max(X))/4)
#ax.set_ylim(int(min(Y))/2,3*int(max(Y))/4)
#ax.set_zlim(int(min(Z))/2,3*int(max(Z))/4)
#


f= open(FileNameSim,"w")

f.write("ITEM: TIMESTEP")

f.write("\n")
f.write("0")
f.write("\n")
f.write("ITEM: NUMBER OF ATOMS")
f.write("\n")
f.write(str(TotalNumberOfPolymer*PolymerLength))
f.write("\n")
f.write("ITEM: BOX BOUNDS pp pp pp")
f.write("\n")
for i in range (0,3):
     f.write("0 "+str(SimSize))
     f.write("\n")
  
f.write("ITEM: ATOMS id mol type x y z " )
f.write("\n")

for i in range (0,len(NumPoly)):
    
    for j in range (0,PolymerLength):
        
        f.write(str(i*PolymerLength+j+1)+ " " +str(i+1)+ " "+str(1)+" ")
        f.write(str(ATOM[i*PolymerLength+j][0])+" "+str(ATOM[i*PolymerLength+j][1])+" "+str(ATOM[i*PolymerLength+j][2]))
        f.write( "\n")
        
f.close()
