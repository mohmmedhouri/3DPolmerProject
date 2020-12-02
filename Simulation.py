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


from geomcalc import *
from Elements import Polymer, Monomer
from MC import *

NumOfPoly=200

class Polypool :
    
    
    def __init__(self,Inputfile):
        
        self.PMMA=[]
        
        if (path.exists(Inputfile)):
            print("reading dataFile")
        else:
            
             for i in range (0,NumOfPoly):
               
               self.PMMA.append(Polymer(200))
               
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
        
        if (BOND):
            PointA =np.add( self.PMMA[PI].Chain[MI].Coordinate,dP)
            for i in (self.PMMA[PI].Chain[MI].Bond):
                PointB = self.PMMA[PI].Chain[i].Coordinate
                r_ij = GetRij(PointB,PointA)
                 
                Ebo+=GetEBond(r_ij, r_eq, k_b)
        
        if (BEND and  (MI > 0 or  MI < len(self.PMMA[PI].Chain)-1)):#
            PointA = np.add( self.PMMA[PI].Chain[MI].Coordinate,dP)
            PointB = self.PMMA[PI].Chain[MI-1].Coordinate
            PointC = self.PMMA[PI].Chain[MI+1].Coordinate
            a_ijk = GetAijk (PointA,PointB,PointC)
            Ebe = GetEAngle(a_ijk, a_eq, k_a)
            
            #print("")
            
        if (Torsion):
            Points=np.add( self.PMMA[PI].Chain[MI].Coordinate,dP)
            for i in (self.PMMA[PI].Chain[MI].Tors):
                Points.append(self.PMMA[PI].Chain[i].Coordinate)
            t_ijkl=GetTijkl(Points[0],Points[1],Points[2],Points[3])
            
            Eto=GetETorsion(t_ijkl, v_n, gamma, nfold, paths)
            
        PointA =np.add( self.PMMA[PI].Chain[MI].Coordinate,dP)    
        for i in (self.PMMA[PI].Chain[MI].Neighbor):
            PointB =  PMMA[i[0]].Chain[i[1]].Coordinate
            r_ij = GetRij(PointB,PointA)
            Enb+=GetE_LJ(r_ij, eps_ij, ro_ij)
            
        Ewa=WallEnergy(SimDimention,'cube',Coordinates,Pres,Temp,k_box)    
            
            
        return Ebo+Ebe+Eto+Enb+Ewa    
            
            
             
    
    def UpdateNeighbors(self,Rc):
        """Calculate the neighboors of all monomers if the are located in a distance smaller than rc.
        
         Args:
             Rc (float*): the cutof raduis of LJ """
        for i in range (0,len(self.PMMA)):
            
            for j in range (0,len(self.PMMA[i].Chain)):
                self.PMMA[i].Chain[j].Neighbor = []
                for c in range (0,len(self.PMMA)):
                    
                    for k in range (0,len(self.PMMA[i].Chain)):
                    
                        
                       Dist2 = GetR2ij (self.PMMA[i].Chain[j].Coordinate,self.PMMA[c].Chain[k].Coordinate)
                       
                       if (Dist2<rc**2 or  not Dist2 == 0):
                           
                           self.PMMA[i].Chain[j].Neighbor.append([c,k,Dist2**0.5])
                           
                           