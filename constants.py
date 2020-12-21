# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 18:02:00 2020

@author: mhdho
"""


"""Constant values for use in other modules within mmlib.

Contains physical constants, arbitrary values, and dictionary mappings for
plotting (listed in alphabetical order).
"""

import math

### 

N = 200 # Chain length
n = 20 # Number of polymer chains 

##All the constants in the simulation program

### LJ potential
BOND = True
BEND = True
Torsion = True

eps_ij = 0.03

ro_ij = 0.04

### Torsion
v_n = 1

gamma = 1

nfold = 1

paths = 1

##Bend
a_eq = 1

k_a = 1

## Bond
r_eq = 1

k_b = 1

## Simulation L

SimDimention = [30,30,30]

NUMDIM = 3

Coordinates = 50

Pres = 1 

Temp = 500 # start temperature

TempE = 200 #end temperature

k_box = 500 

T_Mode = "V"    #Char V is vector other values set the linear cooling rate ! note the equation mode is still under test

TempVector =[455,300,60]

TempInterval =[6,4]

TS = 20  ## steps in ns

CR = (TempE)/TS

NeighborUp = 10 # Neighbors list update

PrintingT = 2 #

KB = 0.001987204

# Threshold beyond covalent radii sum to determine bond cutoff.
BONDTHRESHOLD = 1.2

# Number of Cartesian dimensions
NUMDIM = 3

# Conversion from degrees to radians
DEG2RAD = math.pi / 180.0

# Boltzmann constant [kcal/(mol*K)].
KB = 0.001987204

# Conversion from [kcal*A^3/mol] to [Pa] for pressure.
KCALAMOL2PA = 69476.95

# Conversion of kinetic energy from [amu*A^2/ps^2] to [kcal/mol].
KIN2KCAL = 0.00239005736

# Conversion from radians to degrees.
RAD2DEG = 180.0 / math.pi

# Default optimization criteria keyword dictionary.
# [delta_e, grad_rms, grad_max, disp_rms, disp_max]

OPTCRITERIAREFS = {
    'loose':     [1.0E-4,  1.0E-3, 2.0E-3, 1.0E-2, 2.0E-2],
    'default':   [1.0E-6,  1.0E-4, 2.0E-4, 1.0E-3, 2.0E-3],
    'tight':     [1.0E-8,  1.0E-5, 2.0E-5, 1.0E-4, 2.0E-4],
    'verytight': [1.0E-10, 1.0E-6, 2.0E-6, 1.0E-5, 2.0E-5]}



# Factor by which to adjust the initial line search step size between steps.
OPTSTEPADJUSTOR = math.sqrt(2)



# Fraction of image width which is covered by the plot field.
PERCENTIMAGEPLOT = 0.75



# Unit conversion between points and inches
POINTSPERINCH = 72



# Legend labels, line colors, and plotting priority for properties.
# [energy_term, print_priority, line_color, index]
PROPERTYDICTIONARY = {
    'e_total':    ['Total',      12, '#000000', 1],
    'e_kin':      ['Kinetic',    11, '#007D34', 2],
    'e_pot':      ['Potential',   1, '#C10020', 3],
    'e_nonbond':  ['Non-bonded',  7, '#0000FF', 4],
    'e_bonded':   ['Bonded',      2, '#FF6800', 5],
    'e_boundary': ['Boundary',   10, '#551A8B', 6],
    'e_vdw':      ['Vdw',         9, '#00BFFF', 7],
    'e_elst':     ['Elst',        8, '#EEC900', 8],
    'e_bond':     ['Bonds',       3, '#F08080', 9],
    'e_angle':    ['Angles',      4, '#90EE90', 10],
    'e_tors':     ['Torsions',    6, '#FF83FA', 11],
    'e_oop':      ['Outofplanes', 5, '#A9A9A9', 12]}

# Physical property keys for output file data labels.
PROPERTYKEYS = [
    'e_total', 'e_kin', 'e_pot', 'e_nonbond', 'e_bonded', 'e_boundary',
    'e_vdw', 'e_elst', 'e_bond', 'e_angle', 'e_tors', 'e_oop', 'temperature',
    'pressure']




# Gas constant in units of [amu*A^2/(ps^2*K)].
RGAS = 0.83144598

# Dictionary of order-of-magnitude axis tick labels.
TICCHARS = {0: '', 1: 'k', 2: 'M', 3: 'B', 4: 'T', 5: 'P'}


