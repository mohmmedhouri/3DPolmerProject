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
import math
from scipy import signal

from Energycalc import *




        
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
        

    def __init__(self, infile_name):
    
        self.totconf = 100000
        self.conf = 0
        self.dispmag = 0.1
        self.dispinc = math.log(2.0)
        self.n_accept = 0
        self.n_reject = 0
        self.dispconf = 100
        self.energyconf = 100
        self.geomconf = 100
        Simulation.__init__(self, infile_name)
        self.rand_disp = numpy.zeros((self.mol.n_atoms, const.NUMDIM))
        self.cprintchar = 7
        
    def Run():
        
        
             self._OpenOutputFiles()
             self._ZeroVels()
             numpy.random.seed(self.random_seed)
             self.mol.GetEnergy('standard')
             self._CheckPrint(0, print_all=True)
             previous_energy = self.mol.e_total
             
             
             while self.conf < self.totconf:
                  self._GetRandDisp()
                  self._DispCoords(self.rand_disp)
                  self.mol.GetEnergy('standard')
                  delta_e = self.mol.e_total - previous_energy
                  bf = math.exp(min(1.0, -delta_e / (const.KB * self.temperature)))
                  if bf >= numpy.random.random():
                      self._CheckPrint(1)
                      self.conf += 1
                      self.n_accept += 1
                      previous_energy = self.mol.e_total
                  else:
                      self._DispCoords(-self.rand_disp)
                      self.n_reject += 1
                      self._CheckDisp()
                      self._CheckPrint(0)
                      self._CloseOutputFiles()    
              
              
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
        
        
        
        
        def _PrintEnergyHeader(self):
            """Print header of energy output columns to file."""
            e = self.efile
            e.write('#\n# INPUTFILE %s' % (self.infile))
            e.write('\n#\n# -- INPUT DATA --\n#')
            e.write('\n# MOLFILE %s' % (self.mol.infile))
            e.write('\n# ENERGYOUT %s' % (self.energyout))
            e.write('\n# GEOMOUT %s' % (self.geomout))
            e.write('\n# RANDOMSEED %i' % (self.random_seed))
            e.write('\n# TEMPERATURE %.6f K' % (self.temperature))
            e.write('\n# BOUNDARY %.6f A' % (self.mol.boundary))
            e.write('\n# BOUNDARYSPRING %.6f kcal/(mol*A^2)' % (self.mol.k_box))
            e.write('\n# BOUNDARYTYPE %s' % (self.mol.boundary_type))
            e.write('\n# STATUSTIME %.6f s' % (self.statustime))
            e.write('\n# ENERGYCONF %i' % (self.energyconf))
            e.write('\n# GEOMCONF %i' % (self.geomconf))
            e.write('\n# TOTALCONF %i' % (self.totconf))
            e.write('\n#\n# -- ENERGY DATA --\n#')
            e.write('\n# energy terms [kcal/mol] vs. configuration\n')
            e.write('#  conf        e_pot  e_nonbond   ')
            e.write('e_bonded e_boundary      e_vdw     e_elst     e_bond    ')
            e.write('e_angle     e_tors      e_oop\n')
