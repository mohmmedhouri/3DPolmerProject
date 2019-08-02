# -*- coding: utf-8 -*-
"""
Created on Sat May  4 06:15:22 2019

@author: mhdho
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 17:14:13 2019

@author: alhouri
"""

from Allocation import *
import math as MT
from helper import *  
import numpy as np

from Global import *


class Simulation :
    
    def __init__(self,PolyPool,Grid):
        self.SimPolymers=PolyPool
        self.SimSpace=Grid
    
    
    def TempControl(self,dt):#Only uniform temperture changer
        for i in range (0,len(self.SimSpace)-1):
            for j in range (0,len(self.SimSpace[i])):
                self.SimSpace[i][j].Tempr=self.SimSpace[i][j].Tempr+dt
                
    
    def MoveEstimation(self,PolyIndex,MonoIndex):
        
        Counter=0
        #dx,dy=NextPOT7Move(self.SimPolymers[PolyIndex].chain[MonoIndex].x,self.SimPolymers[PolyIndex].chain[MonoIndex].y,n)
        for dx in range (-1,2):
            for dy in range (-1,2):
                if dx==0 and dy==0 :
                    continue  
                else :
                    Freespace=True
                    Xn=(self.SimPolymers[PolyIndex].chain[MonoIndex].x+dx)
                    Yn=(self.SimPolymers[PolyIndex].chain[MonoIndex].y+dy)
                    if self.SimSpace[Xn][Yn].status == Occ: 
                        Freespace=False
                    
                    Bondlength=True 
                    if (MonoIndex > 0 and MonoIndex <len((self.SimPolymers[PolyIndex].chain))-1):
                           Xp=(self.SimPolymers[PolyIndex].chain[MonoIndex-1].x)
                           Yp=(self.SimPolymers[PolyIndex].chain[MonoIndex-1].y)
                           Xa=(self.SimPolymers[PolyIndex].chain[MonoIndex+1].x)
                           Ya=(self.SimPolymers[PolyIndex].chain[MonoIndex+1].y)
                           
                           a2=MT.pow(Xa-Xn,2)+MT.pow(Ya-Yn,2)
                           b2=MT.pow(Xp-Xn,2)+MT.pow(Yp-Yn,2)
                           if (a2>2 or b2>2):
                               
                                Bondlength= False
                           
                    Ei=self.Hamiltonian(PolyIndex,MonoIndex,0,0)
                    Ef=self.Hamiltonian(PolyIndex,MonoIndex,dx,dy)
                    dE=Ef-Ei
                   
                    D=MT.exp(dE/(kb*self.SimSpace[Xn][Yn].Tempr))
                    
                    if (Bondlength and Freespace ):
                        Counter=min(1,D)+Counter
        
        return (Counter/8)
    
    
    def Movmement (self):
        #this function calculate the mobility of each monomer
        Counter=0
        for i in range (0,len (self.SimPolymers)-1):
            for j in range  (0,len(self.SimPolymers[i].chain)-1):
                Counter=Counter+self.MoveEstimation(i,j)
                
        
        return (Counter/(len (self.SimPolymers)*len(self.SimPolymers[i].chain)))
                
                
                
            
            
              
        
        
        
    
    
    
    
    def DensityEstimationLocal(self,X1,Y1,X2,Y2):
        #X1<X2, Y1<Y2
        # Here we calculate each monomer possibility to move
        TotalinSpace=abs(X1-X2)*abs(Y1-Y2)
        
        counter=0
        for i in range (X1,X2):
            for j in range (Y1,Y2):
                if self.SimSpace[i][j].status==Occ:
                    counter=counter+1
        
        return counter/TotalinSpace
                    
   
    def Move_possiblilty(self,PolyIndex,MonoIndex):
         dx,dy=NextPOT7Move(self.SimPolymers[PolyIndex].chain[MonoIndex].x,self.SimPolymers[PolyIndex].chain[MonoIndex].y,n)
         #dx,dy= -1,0 # remove for the test 
#         print("dx", dx,"dy",dy )
         Xn=(self.SimPolymers[PolyIndex].chain[MonoIndex].x+dx)
         Yn=(self.SimPolymers[PolyIndex].chain[MonoIndex].y+dy)
         EmptySlot = False
#         print("Status",self.SimSpace[Xn][Yn].status)
         if  self.SimSpace[Xn][Yn].status == Occ:
              return False 
         
         
         if (MonoIndex > 0 and MonoIndex <len((self.SimPolymers[PolyIndex].chain))-1):
               Xp=(self.SimPolymers[PolyIndex].chain[MonoIndex-1].x)
               Yp=(self.SimPolymers[PolyIndex].chain[MonoIndex-1].y)
               Xa=(self.SimPolymers[PolyIndex].chain[MonoIndex+1].x)
               Ya=(self.SimPolymers[PolyIndex].chain[MonoIndex+1].y)
               
               a2=MT.pow(Xa-Xn,2)+MT.pow(Ya-Yn,2)
               b2=MT.pow(Xp-Xn,2)+MT.pow(Yp-Yn,2)
               if (a2>2 or b2>2):
                    
                    return False
               
         Ei=self.Hamiltonian(PolyIndex,MonoIndex,0,0)
         Ef=self.Hamiltonian(PolyIndex,MonoIndex,dx,dy)
         dE=Ef-Ei
        
         D=MT.exp(dE/(kb*self.SimSpace[Xn][Yn].Tempr))
         
         r=np.random.random(1)[0]
         
             # Slithering 
         
         if (r<=D):
             if (MonoIndex==0 ):
#                
                 Xo=self.SimPolymers[PolyIndex].chain[MonoIndex].x
                 Yo=self.SimPolymers[PolyIndex].chain[MonoIndex].y
                 self.SimSpace[self.SimPolymers[PolyIndex].chain[MonoIndex].x]\
                 [self.SimPolymers[PolyIndex].chain[MonoIndex].y].status=Free
                 self.SimSpace[Xn][Yn].status=Occ
                 self.SimPolymers[PolyIndex].chain[MonoIndex].x=Xn
                 self.SimPolymers[PolyIndex].chain[MonoIndex].y=Yn
                 for i in range (len((self.SimPolymers[PolyIndex].chain))-1,1,-1):
                     self.SimPolymers[PolyIndex].chain[i].x=self.SimPolymers[PolyIndex].chain[i-1].x
                     self.SimPolymers[PolyIndex].chain[i].y=self.SimPolymers[PolyIndex].chain[i-1].y
                     if (i==len(self.SimPolymers[PolyIndex].chain)-2):
                         self.SimSpace[self.SimPolymers[PolyIndex].chain[len(self.SimPolymers[PolyIndex].chain)-2].x]\
                         [self.SimPolymers[PolyIndex].chain[len(self.SimPolymers[PolyIndex].chain)-2].y].status=Free
                         
                 self.SimPolymers[PolyIndex].chain[1].x=Xo
                 self.SimPolymers[PolyIndex].chain[1].y=Yo
                 return True
                         

             elif ( MonoIndex == len((self.SimPolymers[PolyIndex].chain))-1) :
#
                 Xo=self.SimPolymers[PolyIndex].chain[MonoIndex].x
                 Yo=self.SimPolymers[PolyIndex].chain[MonoIndex].y
                 self.SimSpace[self.SimPolymers[PolyIndex].chain[MonoIndex].x][self.SimPolymers[PolyIndex].chain[MonoIndex].y].status=Free
                 self.SimSpace[Xn][Yn].status=Occ
                 self.SimPolymers[PolyIndex].chain[MonoIndex].x=Xn
                 self.SimPolymers[PolyIndex].chain[MonoIndex].y=Yn
#                 print("Index", MonoIndex,"New X",self.SimPolymers[PolyIndex].chain[MonoIndex].x,"New Y",self.SimPolymers[PolyIndex].chain[MonoIndex].y)
                 
                 for i in range (0,len((self.SimPolymers[PolyIndex].chain))-2):
                     self.SimPolymers[PolyIndex].chain[i].x=self.SimPolymers[PolyIndex].chain[i+1].x
                     self.SimPolymers[PolyIndex].chain[i].y=self.SimPolymers[PolyIndex].chain[i+1].y
                     
#                     print("Index", i,"New X",self.SimPolymers[PolyIndex].chain[i].x,"New Y",self.SimPolymers[PolyIndex].chain[i].y)
                     if (i==1):
                          self.SimSpace[self.SimPolymers[PolyIndex].chain[0].x][self.SimPolymers[PolyIndex].chain[0].y].status=Free
#                          print("Index", 0,"New X",self.SimPolymers[PolyIndex].chain[0].x,"New Y",self.SimPolymers[PolyIndex].chain[0].y)
                 self.SimPolymers[PolyIndex].chain[len((self.SimPolymers[PolyIndex].chain))-2].x=Xo
                 self.SimPolymers[PolyIndex].chain[len((self.SimPolymers[PolyIndex].chain))-2].y=Yo
                    
                 return True
#                     
                     

                 
             else:    
                 self.SimSpace[self.SimPolymers[PolyIndex].chain[MonoIndex].x][self.SimPolymers[PolyIndex].chain[MonoIndex].y].status=Free
                 self.SimSpace[Xn][Yn].status=Occ
                 self.SimPolymers[PolyIndex].chain[MonoIndex].x=Xn
                 self.SimPolymers[PolyIndex].chain[MonoIndex].y=Yn
#                 print("Accepted")
                 return True
             
         
              
              
         #case one check if the new position is free 
         
         
         return False
         
        
    def Hamiltonian (self,PolyIndex,MonoIndex,dx,dy,dz): # Energy calculation function
        # energy calcualtion
        Ebo=0 # the value of the bond energy
        Ebe=0 # the value of the bend energy
        Eth=0 # the value of the Trosion energy
        Elj=0 # LJ potential
        # bond potential calculation
        
        if MonoIndex < len(self.SimPolymers[PolyIndex].chain)-1 :

            lx= MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex].x+dx) -  (self.SimPolymers[PolyIndex].chain[MonoIndex+1].x),2)
            ly= MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex].y+dy) -  (self.SimPolymers[PolyIndex].chain[MonoIndex+1].y),2)
            lz= MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex].z+dz) -  (self.SimPolymers[PolyIndex].chain[MonoIndex+1].z),2)
        else:

            lx= MT.pow (((self.SimPolymers[PolyIndex].chain[MonoIndex].x)+dx) -(self.SimPolymers[PolyIndex].chain[MonoIndex-1].x),2)
            ly= MT.pow (((self.SimPolymers[PolyIndex].chain[MonoIndex].y)+dy) -(self.SimPolymers[PolyIndex].chain[MonoIndex-1].y),2)
            lz= MT.pow (((self.SimPolymers[PolyIndex].chain[MonoIndex].z)+dz) -(self.SimPolymers[PolyIndex].chain[MonoIndex-1].z),2)
        
        
     
            

        l=MT.sqrt(lx+ly+lz)
        Ebo=kbo*MT.pow((l0-l),2) # Ebond
        
        

        if  MonoIndex <len(self.SimPolymers[PolyIndex].chain)-3 :
            ax= MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex].x+dx) - (self.SimPolymers[PolyIndex].chain[MonoIndex+2].x),2)
          
            ay= MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex].y+dy) - (self.SimPolymers[PolyIndex].chain[MonoIndex+2].y),2)
            
            az= MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex].z+dz) - (self.SimPolymers[PolyIndex].chain[MonoIndex+2].z),2)

            
            
            bx =MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex].x+dx) - (self.SimPolymers[PolyIndex].chain[MonoIndex+1].x),2)

            by =MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex].y+dy) - (self.SimPolymers[PolyIndex].chain[MonoIndex+1].y),2)
            
            bz =MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex].z+dz) - (self.SimPolymers[PolyIndex].chain[MonoIndex+1].z),2)

            
            cx= MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex+1].x) - (self.SimPolymers[PolyIndex].chain[MonoIndex+2].x),2)

            cy= MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex+1].y) - (self.SimPolymers[PolyIndex].chain[MonoIndex+2].y),2)
            
            cz= MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex+1].z) - (self.SimPolymers[PolyIndex].chain[MonoIndex+2].z),2)

            
        else:
            
            ax= MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex].x+dx) - (self.SimPolymers[PolyIndex].chain[MonoIndex-2].x),2)
            ay= MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex].y+dy) - (self.SimPolymers[PolyIndex].chain[MonoIndex-2].y),2)
            az= MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex].z+dz) - (self.SimPolymers[PolyIndex].chain[MonoIndex-2].z),2)

            
            
            bx =MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex].x+dx) - (self.SimPolymers[PolyIndex].chain[MonoIndex-1].x),2)
            by =MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex].y+dy) - (self.SimPolymers[PolyIndex].chain[MonoIndex-1].y),2)
            bz =MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex].z+dz) - (self.SimPolymers[PolyIndex].chain[MonoIndex-1].z),2)

            
            cx= MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex-1].x) - (self.SimPolymers[PolyIndex].chain[MonoIndex-2].x),2)
            cy= MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex-1].y) - (self.SimPolymers[PolyIndex].chain[MonoIndex-2].y),2)
            cz= MT.pow ((self.SimPolymers[PolyIndex].chain[MonoIndex+1].z) - (self.SimPolymers[PolyIndex].chain[MonoIndex-2].z),2)
            
            
            
            
        a  =  MT.sqrt(ax+ay+az)
        b  =  MT.sqrt(bx+by+bz)
        c  =  MT.sqrt(cx+cy+cz)
       
        a2 =  ax+ay+az
        b2 =  bx+by+bz
        c2 =  cx+cy+cz
        
        
        
        
                
        if a == 0.0 or b==0.0 or c ==0.0:
             Ebe=0
                  
        else:

            

             A= MT.acos((b2+c2-a2)/(2*b*c))
             B= MT.acos((a2+c2-b2)/(2*a*c))
             C= MT.acos((b2+a2-c2)/(2*b*a))
             
             Ebe=kbe*(MT.cos(A)-MT.cos(th0))*(MT.cos(A)-MT.cos(th0)) # E bend
          
        
        
        
        
        Elj=0 # Lenord John energy 
        
        
        for k in [-1,0,1]:
            
            for j in [-1,0,1]:
                
                for i in [-1,0,1]:
                    
                    if i==0 and j==0 and k==0 :
                        continue   
                    else:
                    #print(i,j)
                        if self.SimSpace[(self.SimPolymers[PolyIndex].chain[MonoIndex].x+dx)+i][(self.SimPolymers[PolyIndex].chain[MonoIndex].y+dy)+j][(self.SimPolymers[PolyIndex].chain[MonoIndex].z+dz)+k].status == Occ :
                            Elj=Elj+1 
            
        if dx==0 and dy==0: # i have simplifed the calulcation in this regard
            Elj=Elj-1
       
       
        E=0.1*Elj+0.1*Ebe+0.1*Ebo # minmization of LJ potentioal
        return E
    




                
               
#    def Adanced_Density (self):
#        Counter=0
#        for i in range (0,len (self.SimPolymers)-1):
#            for j in range  (0,len(self.SimPolymers[i].chain)-1):
#                Counter=Counter+2*self.MoveEstimation(i,j)
#        
#        return (Counter)

     
     
     
     
