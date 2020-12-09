# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 06:37:33 2020

@author: mhdho
"""

import math
import numpy as np
import os
import sys
import time
from os import path

import constants as const
import geomcalc
import Energycalc
from Elements import Polymer, Monomer

#from MC import *

#NumOfPoly=200

PMMA = [Polymer(1) for i in range (0,const.n)]


class Polypool :
    
    
    def __init__(self,Inputfile):
        
        self.PMMA=[]
        
        if (path.exists(Inputfile)):
            print("reading dataFile")
            self.PMMA = [Polymer(1) for i in range (0,const.n)]

            len(self.PMMA)
            f=open(Inputfile, "r")
                
            f1=f.readlines()

            for i in range (9,len(f1)):
                
                Data=f1[i].split(" ")
                #print(Data[3])
                MonomerTemp = Monomer(float(Data[3]),float(Data[4]),float(Data[5]),float(Data[0]))
                self.PMMA[int (Data[1])-1].Chain.append(MonomerTemp)
                
               
                
                
            f.close()
            for i in range (0, len(PMMA)):
                self.PMMA[i].Chain.pop(0)
                print()
                self.PMMA[i].LinkBinding()
            
            
        else:
            
            print("No DataFile was found or the file is unreadable")
            
            for i in range (0,const.N):
               
               self.PMMA.append(Polymer(const.n))
            
         

               
    def CalEnrMonv(self,PI,MI,dP):
        """Calculate the energy of the selected monomer after each diplecment.
        
         Args:
             PI (int): The index of the polymer
             MI (int): The index of the monomer
             dP (float*) 3 demplecement of the monomer location"""
        
        Ebo=0
        Ebe=0
        Eto=0
        Enb=0
        Ewa=0
        
        if (const.BOND):
            PointA =np.add( self.PMMA[PI].Chain[MI].Coordinate,dP)
            for i in (self.PMMA[PI].Chain[MI].Bond):
                PointB = self.PMMA[PI].Chain[i].Coordinate
                r_ij = geomcalc.GetRij(PointB,PointA)
                
                if (r_ij >const.BONDTHRESHOLD*1.5 ):
                    return math.inf
                 
                Ebo+=Energycalc.GetEBond(r_ij, const.r_eq, const.k_b)
        
        
        if (const.BEND and  (MI > 0 and  MI < len(self.PMMA[PI].Chain)-1)):#
            PointA = np.add( self.PMMA[PI].Chain[MI].Coordinate,dP)
            PointB = self.PMMA[PI].Chain[MI-1].Coordinate
            PointC = self.PMMA[PI].Chain[MI+1].Coordinate
            a_ijk =geomcalc.GetAijk (PointA,PointB,PointC)
            Ebe = Energycalc.GetEAngle(a_ijk,const.a_eq,const.k_a)
            
            
           
        if (const.Torsion):
            Points = []
            PointsA=np.add( self.PMMA[PI].Chain[MI].Coordinate,dP)
            Points.append(PointsA)
            #print ("Indicies", PI, self.PMMA[PI].Chain[MI].Tors)
            #print (len(self.PMMA[PI].Chain), PI, MI)
            for i in (self.PMMA[PI].Chain[MI].Tors):
               
                #print("Coordinate ", self.PMMA[PI].Chain[i].Coordinate)
                Points.append(self.PMMA[PI].Chain[i].Coordinate)
            t_ijkl=geomcalc.GetTijkl(Points[0],Points[1],Points[2],Points[3])
            
            Eto=Energycalc.GetETorsion(t_ijkl, const.v_n,const.gamma, const.nfold, const.paths)
            
        PointA =np.add( self.PMMA[PI].Chain[MI].Coordinate,dP)    
        for i in (self.PMMA[PI].Chain[MI].Neighbor):

            PointB =  self.PMMA[i[0]].Chain[i[1]].Coordinate #  check neighbors format 
            r_ij =geomcalc. GetRij(PointB,PointA)
            Enb+=Energycalc.GetE_LJ(r_ij, const.eps_ij,const.ro_ij)
            
        Ewa=Energycalc.WallEnergy(const.SimDimention,'cube',self.PMMA[PI].Chain[MI].Coordinate,const.Pres,const.Temp,const.k_box)    
            
            
        return Ebo+Ebe+Eto+Enb+Ewa    
            
            
             
    
    def UpdateNeighbors(self,Rc):
        """Calculate the neighboors of all monomers if the are located in a distance smaller than rc.
        
         Args:
             Rc (float*): the cutof raduis of LJ """     
        print("Updating Neighbors")
        for i in range (0,len(self.PMMA)):
            
            for j in range (0,len(self.PMMA[i].Chain)):
                self.PMMA[i].Chain[j].Neighbor = []
                for c in range (0,len(self.PMMA)):
                    
                    for k in range (0,len(self.PMMA[i].Chain)):
                    
                       if (i == c and j == k):
                           continue
                       else:
                           
                           Dist2 =geomcalc.GetR2ij (self.PMMA[i].Chain[j].Coordinate,self.PMMA[c].Chain[k].Coordinate)
                           
                           if (Dist2<Rc**2 or  not Dist2 == 0):
                               
                               self.PMMA[i].Chain[j].Neighbor.append([c,k,Dist2**0.5])
                           
                           