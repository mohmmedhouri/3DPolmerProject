# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 17:02:29 2020

@author: mhdho
"""

### Simulation engin #####


import matplotlib.pyplot as plt
import os.path
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
        

    def __init__(self, infile_name,T,L):
        self.temperature = T
        self.Length = L
        self.totconf = 100000
        self.conf = 0
        self.dispmag = 0.1
        self.dispinc = math.log(2.0)
        self.n_accept = 0
        self.n_reject = 0
        self.dispconf = 100
        self.energyconf = 100
        self.geomconf = 100
        self.Polypool.__init__(infile_name)
        #self.rand_disp = numpy.zeros((self.mol.n_atoms, const.NUMDIM))
        self.cprintchar = 7
        
    def _GetRandDisp(self):
        """Generate random displacment vector for coordinates.
        
        Random trial displacements for MMC are selected from a Gaussian distribution
        of mu = 0.0 and sigma = 'dispmag' attribute for all 3N atomic coordinates."""
    
        self.rand_disp.fill(0.0)
        for i in range(self.mol.n_atoms):
          for j in range(const.NUMDIM):
            randval = numpy.random.normal(0.0, self.dispmag)
            self.rand_disp[i][j] = numpy.random.normal(0.0, self.dispmag)
        
    def _ChangeDisp(self):
        """Change root-mean-square magnitude of displacement vector.
        
        The MMC random displacement vector has mu = 0.0, and sigma = 'dispmag'
        chosen to best approach 50% acceptance ratio. Increase 'dispmag' when
        'p_accept' > 0.5 and vice versa. """
       
        p_accept = float(self.n_accept) / float(self.n_reject + self.n_accept)
        self.n_accept, self.n_reject = 0, 0
        self.dispmag *= math.exp(2.0 * self.dispinc * (p_accept - 0.5))
        
    def PrintPolymer(self,Filename):
        """Print the Polymer coorindatess.
        
        Args:
           Filename_out (str): Path of simulation Filename_out. May be relative or
           absolute, though absolute is safer."""
           
        RESULTSTEXT = str(Filename)+".txt"
        
        
        f=open(RESULTSTEXT, "w")
        
#        f.write("The coordinates and the simulation data:")
#        f.write("\n")
#        f.write("Temperature in K:")
#        f.write("\n")
#        f.write(str(self.temperature))
#        f.write("\n")
#        f.write("Simulation volume in A:")
#        f.write("\n")
#        for i in range (0,3):
#            f.write(str(self.Length))
#            f.write("\t")
#        f.write("\n")
#        f.write("Monomers ID POID X Y Z E")
#        f.write("\n")
        
        
        f.write("ITEM: TIMESTEP")
        f.write("\n")
        f.write(str(self.conf))
        f.write("\n")
        f.write("ITEM: NUMBER OF ATOMS")
        f.write(str(len(self.Polypool.PMMA)*len(self.Polypool.PMMA[0].Chain)))
        f.write("\n")
        f.write("ITEM: BOX BOUNDS pp pp pp")
        f.write("\n")
        for i in range (0,3):
             f.write("0 "+str(self.Length))
             f.write("\n")
      
        f.write("ITEM: ATOMS id mol type x y z E" ) 
        Counter = 0
        for i in range (len(self.Polypool.PMMA)):
            
            for j in range (0,len(self.Polypool.PMMA[i].Chain)):
                Counter += 1
                
                f.write(str(Counter))
                
                f.write("\t")
                
                f.write(str(i))
                
                f.write("\t")
                
                f.write(str(self.Polypool.PMMA[i].Chain.X))
                
                f.write("\t")
                
                f.write(str(self.Polypool.PMMA[i].Chain.Y))
                
                f.write("\t")
                
                f.write(str(self.Polypool.PMMA[i].Chain.Y))
                
                f.write(str(self.Polypool.PMMA[i].Chain.E))
                
                f.write("\t")
            
                f.write("\n")
            
        f.close()
    
    
    def Run(self):
        
        
             self._OpenOutputFiles()
             #self._ZeroVels()
             numpy.random.seed(self.random_seed)
             self.Polypool.UpdateNeighbors(3)#add the RC
             ####calculate the energy
             for i in range (0,self.Polypool.PMMA):
                 for j in range (0,len(self.Polypool.PMMA[i].Chain)):
                     self.Polypool.PMMA[i].Chain[j].E = self.Polypool.CalEnrMonv(i,j,[0,0,0])
             
             
             
             NPU = 0 # the Interval at which the neighbors to be updated
             
             Pri = 0
             
             while self.conf < self.totconf:
                 for i in range (0,self.Polypool.PMMA):
                     for j in range (0,len(self.Polypool.PMMA[i].Chain)):
                         
                          PI = random.randint(0,len(self.Polypool.PMMA)-1)
                         
                          MI = random.randint(0,len(self.Polypool.PMMA[PI].Chain)-1)
                         
                          DP = []
                          for K in range (0,3):
                             DP.append(numpy.random.normal(0.0, self.dispmag))
                             
                          NewE= self.Polypool.CalEnrMonv(PI,MI,DP)
                  
                          OldE = self.Polypool.PMMA[PI].Chain[MI].E

                          delta_E = NewE - OldE
                          
                          bf = math.exp(min(1.0, -delta_E / (const.KB * self.temperature)))
                          
                          if bf >= numpy.random.random():
                              
                              self.n_accept += 1
                              self.Polypool.PMMA[PI].Chain[MI].E = NewE
                              self.Polypool.PMMA[PI].Chain[MI].X += DP[0]
                              self.Polypool.PMMA[PI].Chain[MI].Y += DP[1]
                              self.Polypool.PMMA[PI].Chain[MI].Z += DP[2]
                              self.Polypool.PMMA[PI].Chain[MI].Coordinate = [self.Polypool.PMMA[PI].Chain[MI].X,self.Polypool.PMMA[PI].Chain[MI].Y,self.Polypool.PMMA[PI].Chain[MI].Z]
                              
                          else:
                              
                              self.n_reject += 1
                              #self._CheckDisp()
                              #self._CheckPrint(0)
                              #self._CloseOutputFiles()
                          
                 self._CheckDisp()
                 self.conf += 1
                 self._CheckPrint(1)
                 if (NPU < 1000):
                     NPU +=1
                 else:
                     self.Polypool.UpdateNeighbors(3)
                     NPU=0
                     
                 if (Pri < 1000):
                     Pri +=1
                 else:
                     self.Polypool.PrintPolymer(str(self.conf))
                     Pri=0                    
    
    
        