# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 08:29:18 2019
@author: alhouri
"""

import numpy 
import pylab 
import random 


Free="f" # to determine status of the lattice cell
Occ="o"

# movment direction
North=[0,1]
South=[0,-1]
East=[1,0]
West=[-1,0]

#test=[North, North, North]



PLYlen=800  #polymerlength 

# defining the number of steps 
n = 1500

# Class lattice 2d

class Lattice :
    
    def __init__(self,status,Tempr):
        self.status=status
        self.Tempr=Tempr



class SAWI :   #this will create the 
    
    def __init__(self):
        self.POS=[] # the possibility 
        self.Choice=[] # the previous taken choice from the choice array
        self.Xc=0 # Current Xc location
        self.Yc=0 # Current Yc location


# a function to define the neighbor
def NeighborCheck(X,Y, Grid):
    Choice=[]
    
    if (Grid[int(X-1)][Y]).status == Free:
        Choice.append(West)
    if (Grid[int(X+1)][Y]).status == Free:
        Choice.append(East)
    if (Grid[int(X)][int(Y-1)]).status == Free:
        Choice.append(South)
    if (Grid[X][int(Y+1)]).status == Free:
        Choice.append(North)
    return Choice


class monomer:
    
    def __init__(self,INDEX,X,Y):
        self.Index=INDEX  # index of the monomer in the chain 
        self.x=X # the location of the monomer in the lattice
        self.y=Y # the location of the monomer in the lattice
    def setlocation(self,NewIndex,X,Y):
        self.Index=NewIndex
        self.x=X
        self.y=Y

class polymer:
    
    def __init__(self):
        self.chain=[monomer(INDEX=i,X=0,Y=0) for i in range(n)]


def random_walk(X,Y,SAW,Grid):
#check the neighbor possibility
    TempELe=SAWI()
    i=0
    j=0
    ChoicePool=[]
    Trails=[]
    x = numpy.zeros(PLYlen) 
    y = numpy.zeros(PLYlen) 
    #print(NeighborCheck(X,Y,Grid))
    while i< PLYlen and j<PLYlen*4 :
        TempELe=SAWI()
        ChoicePool = NeighborCheck(X,Y,Grid)
        #print(ChoicePool)
        print ('X= ', X, 'Y= ', Y  )
        if ChoicePool == [] and i==0 :
            print ('empty')
            return Grid
        
        elif ChoicePool == [] or Grid[X][Y]==Occ:
            #print("Better situation but still empty choice") was used for debugging
            
            #print (i)
            X=Trails[i-1].Xc
            Y=Trails[i-1].Yc
            Grid[X][Y].status= Free
            Trails.pop()
            i=i-1
            #print(i)
            #print('Length is  ',len(Trails[i-1].POS))
            j=j+1
            
            # remove all the elements which contains only one choice till we reached the old choice which led to the dead end

            while (len(Trails[i-1].POS)<=1):
                # print('You are in the while loop And the length is:  ',len(Trails[i-1].POS)) 
                X=Trails[i-1].Xc
                Y=Trails[i-1].Yc
                Grid[X][Y].status= Free
                Trails.pop()
                i=i-1
                j=j+1
                #print('New length is:  ' ,len(Trails[i-1].POS))
            
            #print('The Choice was  ', Trails[i-1].Choice)
            
            #print(Trails[i-1].POS)
            
            NEWINDI=Trails[i-1].POS.index(Trails[i-1].Choice) # remove the old used index
            
            #print('NEWINDI  ',NEWINDI)
            
            del Trails[i-1].POS[NEWINDI]
            
            #print('New array of choices  ',Trails[i-1].POS)
            if(len(Trails[i-1].POS)>1):
                Indicator=random.randint(0,len(Trails[i-1].POS)-1)
               
            
            else:
                Indicator=0
                
            
            Trails[i-1].Choice=Trails[i-1].POS[Indicator]
            # print('New Choice ', Trails[i-1].Choice)
                
            X=X+Trails[i-1].Choice[0]
            Y=Y+Trails[i-1].Choice[1]
                
            

            
        else :            
            # print("ChoicePool   ",ChoicePool)
            
            Grid[X][Y].status = Occ
            TempELe.POS=ChoicePool
            Indicator=random.randint(0,len(ChoicePool)-1)
            TempELe.Choice=ChoicePool[Indicator]
            
            TempELe.Xc=X
            TempELe.Yc=Y
            Trails.append(TempELe)
           # print(Trails[i-1].POS)
            X=X+TempELe.Choice[0]
            Y=Y+TempELe.Choice[1]
            i=i+1
            j=j+1
            print(j)
            # print(i)
           
            
           #pylab.title("Random Walk ($n = " + str(n) + "$ steps)")
    x=[]
    y=[]
    for i in range (0,len(Trails)-1):
        x.append(Trails[i].Xc)
        y.append(Trails[i].Yc)
        
    pylab.plot(x, y,'*-') 
            #pylab.savefig("rand_walk"+str(n)+".png",bbox_inches="tight",dpi=600) 
    pylab.show() 
        
    return Grid,Trails
             
            
#check neighbor    
    

Space =[[Lattice(Free,273) for i in range(n)] for j in range(n)]
#print(Space[0][0].status)


#print(NeighborCheck(3,5, Space))



Space,polymer = random_walk(500,750,2,Space)