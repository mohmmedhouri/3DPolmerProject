# -*- coding: utf-8 -*-
"""
Created on Tue May  7 15:52:21 2019

@author: alhouri
"""

# Here is the main Program 


from Global import*
from Allocation import*
from helper import*
from MCSim import*
from os import system, name 




Space =[[Lattice(Free,273) for i in range(n)] for j in range(n)]

PolymerPool= []

for i in range (Startup, Endlocation, MT.ceil(MT.sqrt(PLYlen))+1):
    
    for j in range (Startup, Endlocation, MT.ceil(MT.sqrt(PLYlen))+2):
        
        Space,polymerChain=allocatePolymer(i,j,Space,PLYlen)
        Data=polymer()
        Data.FillPolymer(polymerChain)
        PolymerPool.append(Data)


# Prouning extra  monomer
del PolymerPool[-1]   

for i in range (0,len(PolymerPool)-1):
     del PolymerPool[i].chain[-1]  
       
print("Wait here")    
Siml=Simulation(PolymerPool,Space)

Density=Siml.DensityEstimationLocal(600,600,900,900)
print("Density= ",Density)
Siml.TempControl(20 )
print(Siml.SimSpace[300][300].Tempr)
for i in range (0,10):
    Siml.TempControl(10)
    print(Siml.SimSpace[300][300].Tempr,' Density',Siml.Adanced_Density()," Movement_possibility",Siml.Movmement() )
    
    for j in range (0,300):
        for ii in range (0,len(Siml.SimPolymers)-1):
            for jj in range (0,len(Siml.SimPolymers[ii].chain)-1):
                Siml.Move_possiblilty(int(ii),int(jj)) 
                
#        PolyIndex=random.randint(0,len(Siml.SimPolymers)-1)
#        MonoInd=random.randint(0,len(Siml.SimPolymers[PolyIndex].chain)-1)
#        print("PolyIndex",PolyIndex,"MonoInd",MonoInd)
#        Siml.Move_possiblilty(int(PolyIndex),int(MonoInd))


Density=Siml.DensityEstimationLocal(Startup,Startup,Endlocation,Endlocation)
print("Density",Density)
Siml.TempControl(0)
for i in range (0,8):
    Siml.TempControl(-10)
    print(Siml.SimSpace[300][300].Tempr,' Density',Siml.Adanced_Density()," Movement_possibility",Siml.Movmement())
    
    for j in range (0,300):
        for ii in range (0,len(Siml.SimPolymers)-1):
            for jj in range (0,len(Siml.SimPolymers[ii].chain)-1):
                Siml.Move_possiblilty(int(ii),int(jj)) 

print(Siml.SimSpace[300][300].Tempr)
x=[]
y=[]

print(Siml.SimSpace[300][300].Tempr)
for i in range (0,len(Siml.SimPolymers)-1):
     for j in range (0,len (Siml.SimPolymers[i].chain)):
         x.append(Siml.SimPolymers[i].chain[j].x)
         y.append(Siml.SimPolymers[i].chain[j].y)

pylab.plot(x, y,"*") 
pylab.grid() 

pylab.show()    
         
        
        # bend potential calculation
        
        
            
        # Helper file should be added 
        
        
        
        
#    def MCMove(self):
#        for i in range (0,len(self.SimPolymers)-1):
#            
#        