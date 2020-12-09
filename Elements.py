# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 21:29:10 2020

@author: mhdho
"""

import matplotlib.pyplot as plt
import os.path
from mpl_toolkits.mplot3d import Axes3D
import time
import random
from scipy import signal
#from Energy import*
import numpy as np 
import math as MT
import constants as const


NumberOfPolymer = 20 # need to define 
class Monomer :
     def __init__(self,X,Y,Z,ID):
          
          self.Coordinate = [X,Y,Z]
          self.X=X
          self.Y=Y
          self.Z=Z
          self.ID=ID # the atome id for  plotting the data
          self.E=0 # monomer energy
          self.r=MT.sqrt(X*X+Y*Y+Z*Z) 
          self.Neighbor = []
          self.Bond =np.zeros(2)
          self.Bend =np.zeros(2)
          self.Tors =np.zeros(3)
          
     
     def SetCoordinate(self,X,Y,Z):
          self.Coordinate = [X,Y,Z]
          self.X=X
          self.Y=Y
          self.Z=Z
          self.r=MT.sqrt(X*X+Y*Y+Z*Z)
          
     
     def GetCoordinate(self):
          
          return np.asarray([self.X,self.Y,self.Z]) 
          

class Polymer:
     
     def __init__(self,PlyLen): # PlyLen the polymer length
          self.Chain=[]
          for i in range (0,PlyLen):
               
               MTemp=Monomer(0,0,0,0)
               
               self.Chain.append(MTemp)
          
          self.CoM=np.asarray([0,0,0]) # Center of Mass
          
          self.RoG=0 # Radius of gyration
          
          self.LinkBinding()
     
     def UpdatePolymer (self,Chainlocation):
          self.Chain=Chainlocation
          
          self.CenterOfMass()
          self.RadiousOfGyration()
          
     
     def Updatelocation(self,MonomerIndex,X,Y,Z):
          
          self.Chain[MonomerIndex].SetCoordinate(X,Y,Z)
          
     def BLC(self,Vol):
        for i in range (0, len(self.Chain)-1):
            #print(i)
            X2= (self.Chain[i+1].X -self.Chain[i].X)**2
            if (X2**0.5> (Vol)/4):
                #print("X BBC need correction")
                #print("I+1 = "+str(self.Chain[i+1].X)+" , I = "+str(self.Chain[i].X))
                self.Chain[i+1].X =abs( self.Chain[i+1].X-Vol)
                #print("New I+1 = "+str(self.Chain[i+1].X)+" , I = "+str(self.Chain[i].X))
                
                
                
                
            Y2= (self.Chain[i+1].Y -self.Chain[i].Y)**2
            if (Y2**0.5> (Vol)/4):
                #print("Y BBC need correction")
                #print("I+1 = "+str(self.Chain[i+1].Y)+" , I = "+str(self.Chain[i].Y))
                self.Chain[i+1].Y =abs( self.Chain[i+1].Y-Vol)
                #print("New I+1 = "+str(self.Chain[i+1].Y)+" , I = "+str(self.Chain[i].Y))
                
            Z2= (self.Chain[i+1].Z -self.Chain[i].Z)**2
            if (Z2**0.5> (Vol)/4):
                #print("Z BBC need correction")
                #print("I+1 = "+str(self.Chain[i+1].Z)+" , I = "+str(self.Chain[i].Z))
                self.Chain[i+1].Z =abs( self.Chain[i+1].Z-Vol)
                #print("New I+1 = "+str(self.Chain[i+1].Z)+" , I = "+str(self.Chain[i].Z))
        
        #for i in range (0, len(self.Chain)):
            
            #print("New X = "+str(self.Chain[i].X)+" , Y = "+str(self.Chain[i].Y)+" , Z = "+str(self.Chain[i].Z))
            
            
     def RadiousOfGyration(self):
        
        
        SumX = 0
        SumY = 0
        SumZ = 0
        
        N = len(self.Chain)

        for i in range (0,len(self.Chain)):
            SumX = SumX + self.Chain[i].X
            SumY = SumY + self.Chain[i].Y
            SumZ = SumZ + self.Chain[i].Z
        
        CoM = np.asarray([SumX / N, SumY / N, SumZ / N])
        
        
        SumX = 0
        SumY = 0
        SumZ = 0
        Sum = 0
        for i in range(0, len(self.Chain)):
            SumX = SumX + (self.Chain[i].X - CoM[0]) * (self.Chain[i].X - CoM[0])
            SumY = SumY + (self.Chain[i].Y - CoM[1]) * (self.Chain[i].Y - CoM[1])
            SumZ = SumZ + (self.Chain[i].Z - CoM[2]) * (self.Chain[i].Z - CoM[2])

        #   Sum = Sum + MT.pow(self.Chain[i].X - self.CoM[0], 2) + pow(self.Chain[i].Y - self.CoM[1], 2) + pow(
        #     self.Chain[i].Z - self.CoM[2], 2)
        SumX = SumX / len(self.Chain)
        SumY = SumY / len(self.Chain)
        SumZ = SumZ / len(self.Chain)
        Sum = SumX + SumY + SumZ
        RoG = MT.pow(Sum , 0.5)

        return RoG    
    
    
     def LinkBinding (self):
         
         for i in range (0,len(self.Chain)): # to  create the bonds
             
             if (i == 0):
                 
                 self.Chain[i].Bond = [i,i+1]
                 
             elif (i == len(self.Chain)-1):
                 
                 self.Chain[i].Bond = [i-1,i]
            
             else:
                 self.Chain[i].Bond = [i-1,i+1]
                 
         for i in range (0,len(self.Chain)):#
             
             if (i == 0):
                 
                 self.Chain[i].Bend = [i,i]
                 
             elif (i == len(self.Chain)-2):
                 
                 self.Chain[i].Bend = [i,i]
            
             else:
                 self.Chain[i].Bend = [i-1,i+1]
                 
         for i in range (0,len(self.Chain)):
             
             if (i < len(self.Chain)-3):
                 
                 self.Chain[i].Tors = [i+1,i+2,i+3]
                 
             else :
                 
                
                 self.Chain[i].Tors = [i-1,i-2,i-3]
                
                 
             
            
#PMMA=Polymer(200)    
#print(PMMA.Chain[0].Bond)       
             
"""
PMMA=[]

PMMA = [Polymer(1) for i in range (0,NumberOfPolymer)]
Directory="result."+str(int(i))
#Directory="result.202000"
f=open(Directory, "r")
        
f1=f.readlines()
Data=[]
l0= []
Density =[]
Correction = Data=f1[7].split(" ")
CX= float(Data[0])
Vol =  float(Data[1])-float(Data[0])
step= int( Vol/3)
#step= 3
#step=3

Subvol  =[[[0 for i in range (0,step)] for i in range (0,step)]for i in range (0,step)]


#Subvol = [[[0 for k in range(0,step)] for j in range(0,step)] i for range(0, step)]
for i in range (9,len(f1)):
    
    Data=f1[i].split(" ")
    MonomerTemp = Monomer(float(Data[3])-CX,float(Data[4])-CX,float(Data[5])-CX,float(Data[0]))
    PMMA[int (Data[1])-1].Chain.append(MonomerTemp)
    
    #l0.append(float(Data[0])/3.632093)
    #Density.append(float(Data[1]))
    
    
f.close()
for i in range (0, len(PMMA)):
    PMMA[i].Chain.pop(0)
                        
 """                       