# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 14:00:55 2019

@author: alhouri
"""

import random 




def NextPOT7Move(X,Y,Size):
    if X==0 and Y==0:
        Indicator=random.randint(0,2)
        if Indicator==0:
            dx=0+1
            dy=0
        elif Indicator==1:
            dx=0+1
            dy=0+1
        else:
            dx=0
            dy=0+1
            
    elif X==Size-1 and Y==Size-1:
        Indicator=random.randint(0,2)
        if Indicator==0:
            dx=0-1
            dy=0
        elif Indicator==1:
            dx=0-1
            dy=0-1
        else:
            dx=0
            dy=0-1
        
    elif X==0 and Y==Size-1:
        Indicator=random.randint(0,2)
        if Indicator==0:
            dx=0+1
            dy=0
        elif Indicator==1:
            dx=0+1
            dy=0-1
        else:
            dx=0
            dy=0-1
            
    elif X==Size-1 and Y==0: 
        Indicator=random.randint(0,2)
        if Indicator==0:
            dx=0-1
            dy=0
        elif Indicator==1:
            dx=0-1
            dy=0+1
        else:
            dx=0
            dy=0+1
            
    elif X==0 and Y>0 and Y<Size-1: # case for the edges but not at the cor
        Indicator=random.randint(0,4)
        if Indicator==0:
            dx=0
            dy=0+1
        elif Indicator==1:
            dx=0+1
            dy=0+1
        
        elif Indicator==2:
            dx=0+1
            dy=0
      
        elif Indicator==3:
            dx=0+1
            dy=0-1
    
        else:
            dx=0
            dy=0-1
        
        
    elif X==Size-1 and Y>0 and Y<Size-1: # case for the edges 
        Indicator=random.randint(0,4)
        if Indicator==0:
            dx=0
            dy=0+1
        elif Indicator==1:
            dx=0-1
            dy=0+1
        
        elif Indicator==2:
            dx=0-1
            dy=0
      
        elif Indicator==3:
            dx=0-1
            dy=0-1
    
        else:
            dx=0
            dy=0-1
            
    
    elif X==Size-1 and Y>0 and Y<Size-1: # case for the edges 
        Indicator=random.randint(0,4)
        if Indicator==0:
            dx=0
            dy=0+1
        elif Indicator==1:
            dx=0-1
            dy=0+1
        
        elif Indicator==2:
            dx=0-1
            dy=0
      
        elif Indicator==3:
            dx=0-1
            dy=0-1
    
        else:
            dx=0
            dy=0-1
    
    
    elif Y==0 and X>0 and X<Size-1: # case for the edges 
        Indicator=random.randint(0,4)
        if Indicator==0:
            dx=0-1
            dy=0
        elif Indicator==1:
            dx=0-1
            dy=0+1
        
        elif Indicator==2:
            dx=0
            dy=0+1
      
        elif Indicator==3:
            dx=0+1
            dy=0+1
    
        else:
            dx=0+1
            dy=0
      
    elif Y==Size-1 and X>0 and X<Size-1: # case for the edges 
        Indicator=random.randint(0,4)
        if Indicator==0:
            dx=0-1
            dy=0
        elif Indicator==1:
            dx=0-1
            dy=0-1
        
        elif Indicator==2:
            dx=0
            dy=0-1
      
        elif Indicator==3:
            dx=0-1
            dy=0+1
    
        else:
            dx=0+1
            dy=0
    
    else: #We are in the mid area
        
        Indicator=random.randint(0,7) 
        
        if Indicator == 0:
            dx=0+1
            dy=0
           
        elif Indicator == 1:
            dx=0+1
            dy=0+1
            
        elif Indicator == 2:
            dx=0
            dy=0+1
           
        elif Indicator == 3:
            dx=0-1
            dy=0+1
            
        elif Indicator == 4:
            dx=0-1
            dy=0
            
        elif Indicator == 5:
            dx=0-1
            dy=0-1
            
        elif Indicator == 6:
            dx=0
            dy=0-1
          
        else:
            dx=0+1
            dy=0-1
    return dx,dy