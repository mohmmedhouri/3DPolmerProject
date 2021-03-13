# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 18:29:32 2020

@author: mhdho
"""

import constants as consts

#import Simulation

#import Energycalc

import SimEng

CR = 100000

Configuration = SimEng.MonteCarlo()

Configuration.setMC("PMMA_nc200_cl200_D200.dat", consts.Temp, consts.SimDimention)


Configuration.Run()
