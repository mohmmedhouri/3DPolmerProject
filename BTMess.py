# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 09:25:41 2020

@author: mhdho
"""

import matplotlib.pyplot as plt
import os.path
from mpl_toolkits.mplot3d import Axes3D
import time
import random
import numpy as np

fig = plt.figure()
#fig.suptitle('The effect of  Bond potential parameters on the density of the cooled polymer and cooling rate is 1000 K/sec')

# First subplot
ax = fig.add_subplot(111)


t = np.linspace(10,1,10)

BTH=np.linspace(1.47094504,1.459728385,10)
BTH=BTH/max(BTH)
BTL=np.linspace(1.48095144,1.47844666,10)
BTL=BTL/max(BTL)

#BTH=[]
#BTL=[]
#for i in range (0,13):
#    BTH.append(random.randint(135,167))
#    BTL.append(random.randint(95,118))


plt.plot(t,BTH,marker="o")
plt.plot(t,BTL,marker="o")
plt.xlabel("distance x (0.1 mm)")
plt.ylabel("Normalized Refreactive index")
plt.legend(["10000 K/sec","1000 K/sec"])
plt.grid(True)


