# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 17:02:29 2020

@author: mhdho
"""

### Simulation engin #####


import matplotlib.pyplot as plt
import os
import sys
from mpl_toolkits.mplot3d import Axes3D
import time
import random
import math
import numpy
from scipy import signal

from Energycalc import *

from Simulation import Polypool
import constants as consts

        
"""
Metropolis Monte-Carlo class derived from Simulation base class.
  
  Metropolis Monte-Carlo generates a set of configurations for a system by
  propogating an initial configuration forward in time through random
  displacements of coordinates which are accepted or rejected based on the
  relative potential energy Boltzmann factor.
  
  Args:
    infile_name (str): Path of simulation input file. May be relative or
    absolute, though absolute is safer.
      
  Attributes:
    totconf (int): Total number of steps.
    conf (int): Current simulation step number.
    dispmag (float): Magnitude of average random displacement [Angstrom].
    dispinc (float): Rate constant for 'dispmag' adjustment.
    n_accept (int): Number of accepted MMC trials.
    n_reject (int): Number of rejected MMC trials.
    NeiUpdat (int): Number of steps to update neighboors' list 
    dispconf (int): Number of configurations between adjusting 'dispmag' value.
    energyconf (int): Number of configurations between energy printing.
    geomconf (int): Number of configurations between geometry printing.
    rand_disp (float**): Random displacement vector of coordinates.
    cprintchar (int): Total characters for conf number printing.
"""



class MonteCarlo :
    
    def __init__(self):
        
        self.totconf = 1000
        self.conf = 0
        self.dispmag = 0.1
        self.dispinc = math.log(2.0)
        self.n_accept = 0
        self.n_reject = 0
        self.dispconf = 100
        self.energyconf = 100
        self.geomconf = 100
        self.melt=Polypool("def.txt")
        

    def setMC(self, infile_name, T,L):
        self.temperature = T
        self.Length = L
        self.totconf = consts.TS
        self.conf = 0
        self.dispmag = 0.1
        self.dispinc = math.log(2.0)
        self.n_accept = 0
        self.n_reject = 0
        self.dispconf = 100
        self.energyconf = 100
        self.geomconf = 100
        self.melt=Polypool(infile_name)
        #self.rand_disp = numpy.zeros((self.mol.n_atoms, const.NUMDIM))
        self.cprintchar = 7
        if (not os.path.exists(infile_name)):
             print("No DataFile was found or the file is unreadable Simulation will be stopped")
             sys.exit()
    def _GetRandDisp(self):
        """Generate random displacment vector for coordinates.
        
        Random trial displacements for MMC are selected from a Gaussian distribution
        of mu = 0.0 and sigma = 'dispmag' attribute for all 3N atomic coordinates."""
    
        self.rand_disp.fill(0.0)
        for i in range(self.mol.n_atoms):
          for j in range(const.NUMDIM):
            randval = numpy.random.normal(0.0, self.dispmag)
            self.rand_disp[i][j] = numpy.random.normal(0.0, self.dispmag)
    

    def _TempControl(self, Beta):
         """Update the system temperature based on the cooling rate.
    
         Args:
         Beta (float): Cooling/Heating rate."""
         self.temperature = consts.Temp - (self.conf*Beta)
        
        
    def _ChangeDisp(self):
        """Change root-mean-square magnitude of displacement vector.
        
        The MMC random displacement vector has mu = 0.0, and sigma = 'dispmag'
        chosen to best approach 50% acceptance ratio. Increase 'dispmag' when
        'p_accept' > 0.5 and vice versa. """
       
        p_accept = float(self.n_accept) / float(self.n_reject + self.n_accept)
        self.n_accept, self.n_reject = 0, 0
        self.dispmag *= math.exp(2.0 * self.dispinc * (p_accept - 0.5))
        
    def PrintSimulation(self,Filename):
        """Print the Polymer coorindatess.
        
        Args:
           Filename_out (str): Path of simulation Filename_out. May be relative or
           absolute, though absolute is safer."""
           
        RESULTSTEXT = str(Filename)+".txt"
        
        f=open(RESULTSTEXT, "w")
        
        f.write("ITEM: TIMESTEP")
        f.write("\n")
        f.write(str(self.conf))
        f.write("\n")
        f.write("ITEM: NUMBER OF ATOMS")
        f.write(str(len(self.melt.PMMA)*len(self.melt.PMMA[0].Chain)))
        f.write("\n")
        f.write("ITEM: BOX BOUNDS pp pp pp")
        f.write("\n")
        for i in range (0,3):
             f.write("0 "+str(self.Length))
             f.write("\n")
      
        f.write("ITEM: ATOMS id mol type x y z E" )
        f.write("\n")
        Counter = 0
        for i in range (len(self.melt.PMMA)):
            
            for j in range (0,len(self.melt.PMMA[i].Chain)):
                Counter += 1
                
                f.write(str(Counter))
                
                f.write("\t")
                
                f.write(str(i))
                
                f.write("\t")
                
                f.write(str(self.melt.PMMA[i].Chain[j].X))
                
                f.write("\t")
                
                f.write(str(self.melt.PMMA[i].Chain[j].Y))
                
                f.write("\t")
                
                f.write(str(self.melt.PMMA[i].Chain[j].Y))
                
                f.write(str(self.melt.PMMA[i].Chain[j].E))
                
                f.write("\t")
            
                f.write("\n")
            
        f.close()
    
    
    def Run(self):
        
             self.melt.UpdateNeighbors(3)#add the RC
             ####calculate the energy
             self.PrintSimulation("Data")
             for i in range (0,len(self.melt.PMMA)):
                 for j in range (0,len(self.melt.PMMA[i].Chain)):
                     self.melt.PMMA[i].Chain[j].E = self.melt.CalEnrMonv(i,j,[0,0,0])
           
             NPU = 0 # the Interval at which the neighbors to be updated
             
             Pri = 0 # print interval
             
             while self.conf < self.totconf:
                 for i in range (0,len(self.melt.PMMA)):
                     for j in range (0,len(self.melt.PMMA[i].Chain)):
                         
                          PI = random.randint(0,len(self.melt.PMMA)-1)
                         
                          MI = random.randint(0,len(self.melt.PMMA[PI].Chain)-1)
                         
                          DP = []
                          for K in range (0,3):
                              
                             DP.append(numpy.random.normal(0.0, self.dispmag))
                             
                          NewE= self.melt.CalEnrMonv(PI,MI,DP)
                          
                          if (NewE >= math.pow(10,16)):
                              self.n_reject += 1
                              continue
                  
                          OldE = self.melt.PMMA[PI].Chain[MI].E

                          delta_E = NewE - OldE
                          
                          bf = math.exp(min(1.0, -delta_E / (consts.KB * self.temperature )))
                          
                          if bf >= numpy.random.random():
                              
                              self.n_accept += 1
                              self.melt.PMMA[PI].Chain[MI].E = NewE
                              self.melt.PMMA[PI].Chain[MI].X += DP[0]
                              self.melt.PMMA[PI].Chain[MI].Y += DP[1]
                              self.melt.PMMA[PI].Chain[MI].Z += DP[2]
                              self.melt.PMMA[PI].Chain[MI].Coordinate = [self.melt.PMMA[PI].Chain[MI].X,self.melt.PMMA[PI].Chain[MI].Y,self.melt.PMMA[PI].Chain[MI].Z]
                              
                          else:
                              
                              self.n_reject += 1
                             
                 print("Time ",self.totconf-self.conf, "Temperature ", self.temperature)         
                 self._ChangeDisp()
                 self._TempControl(consts.CR)                 
                 self.conf += 1
                 
                 if (NPU < consts.NeighborUp):
                     NPU +=1
                 else:
                     self.melt.UpdateNeighbors(3)
                     NPU=0
                     
                 if (Pri < consts.PrintingT):
                     Pri +=1
                 else:
                     print("Printing")
                     self.melt.PrintSimulation(str(self.conf))
                     Pri=0                    
    
    
        