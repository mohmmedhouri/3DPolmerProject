# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 08:50:03 2019

@author: alhouri
"""

# Generic python libs


import numpy 
import pylab 
import random 
import math as MT



# Project defined methods 
#from Global import * 


Free="f" # to determine status of the lattice cell
Occ="o" 

PLYlen=200 #polymerlength 
n = 1500 

# movment direction

#test=[North, North, North]



  #polymerlength 

# defining the number of steps 


# Class lattice 2d

class Lattice :
    
    def __init__(self,status,Tempr):
        self.status=status
        self.Tempr=Tempr



class SAWI :   #this will create the indicator for the polymer
    
    def __init__(self):
        self.Xc=0 # Current Xc location
        self.Yc=0 # Current Yc location
        self.Zc=0 # Current Zc location
        self.PointDensity=0


# a function to define the neighbor



class monomer:
    
    def __init__(self,INDEX,X,Y,Z):
        self.Index=INDEX  # index of the monomer in the chain 
        self.x=X # the location of the monomer in the lattice
        self.y=Y # the location of the monomer in the lattice
        self.z=Z # the location of the monomer in the lattice
    def setlocation(self,NewIndex,X,Y,Z):
        self.Index=NewIndex
        self.x=X
        self.y=Y
        self.z=Z

class polymer:
    
    def __init__(self):
        self.chain=[monomer(INDEX=i,X=0,Y=0,Z=0) for i in range(PLYlen)]
        
        
        self.E2E=0 # end to end length
        self.RoG=0 # radius of gyeriation
        self.PolymerRadius=0
        self.RadiusFirst=0
        self.RadiusLast=0
        self.CenterX=0
        self.CenterY=0
        self.CenterZ=0
        
    
    def FillPolymer (self,InputData):
        
        for i in range (0,len(self.chain)-1):
            
            self.chain[i].setlocation(i,InputData[i].Xc,InputData[i].Yc,InputData[i].Zc)
            
    
    def PolymerCircle (self):
        FPoint=0
        SPoint=0
        Diameter=0
        for i in range (0,len(self.chain)):
            for j in range (0,len(self.chain)):
                if j>i :
                    continue
                else:
                    ax=self.chain[i].x
                    ay=self.chain[i].y
                    az=self.chain[i].z
                    
                    bx=self.chain[j].x
                    by=self.chain[j].y
                    bz=self.chain[j].z
                    
                    Dia=MT.sqrt(MT.pow(ax-bx,2)+MT.pow(ay-by,2)+MT.pow(az-bz,2))
                    
                    if(Dia>Diameter):
                       Diameter = Dia
                       FPoint=i
                       SPoint=j
                    
                       
        self.PolymerRadius=Diameter/2
        self.RadiusFirst = FPoint
        self.RadiusLast  = SPoint
        self.CenterX= (self.chain[FPoint].x+self.chain[SPoint].x)//2 
        self.CenterY= (self.chain[FPoint].y+self.chain[SPoint].y)//2 
        self.CenterZ= (self.chain[FPoint].z+self.chain[SPoint].z)//2 
        
    
    def RadiusOfGyeriation(self):
        return 1
          
        
#def allocatePolymer (Xs,Ys,Grid,PLYlen) :
#    
#    WedInd = MT.ceil(MT.sqrt(PLYlen))
#    TempELe=SAWI()
#    PolymerCoodinates=[]
#    for i in range (0,WedInd):
#        
#        for j in range (0,WedInd-1):
#            TempELe=SAWI()
#            Grid[Xs][Ys].status=Occ
#            TempELe.Xc=Xs
#            TempELe.Yc=Ys
#            PolymerCoodinates.append(TempELe)
#            if (i%2 == 0 ):
#                Ys=Ys-1
#            else:
#                Ys=Ys+1
#            
#            #print(Xs,Ys)
#
#            
#            
#            
#        Xs=Xs+1
#        Ys=TempELe.Yc
#    return Grid,PolymerCoodinates
    

def ReadPolymer (Xs,Ys,Zs,Grid,Filename) :
    
   
    TempELe=SAWI()
    PolymerCoodinates=[]
    
    f=open(Filename,"r")
    f1=f.readlines()
    for J in range (0,len(f1)):
        Poly=f1[J].split(',')
        TempELe.Xc=int(Poly[0])+Xs
        TempELe.Yc=int(Poly[1])+Ys
        TempELe.Zc=int(Poly[2])+Zs
        PolymerCoodinates.append(TempELe)
        Grid[TempELe.Xc][TempELe.Yc][TempELe.Zc].status=Occ
        TempELe=SAWI()
    return Grid,PolymerCoodinates
            #print(Xs,Ys)




      
#    
#    
#            
#    
