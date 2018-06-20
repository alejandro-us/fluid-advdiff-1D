#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 17:07:27 2018

@author: luiggi
modified by alex
"""

from Coefficients import Coefficients

class Diffusion1D(Coefficients):
    """
    Class that defines coefficients originated from difussion equation using Finite Volume method
    
    Parameters:
        nvx: Number of volumes\n
        dx: Size space step in a constant grid\n
        Gamma: Permeability coefficient
    """
    
    def __init__(self, nvx = None, Gamma = None, dx = None):
        super().__init__(nvx, dx)
        self._nvx = nvx
        self._Gamma = Gamma
        self._dx = dx

    def __del__(self):
        del(self._Gamma)
        del(self._dx)
    
    def calcCoef(self):
        """
            Calculates the proper coefficients for diffusive movement.
            
            Return:
                None
        """
        aE = self.aE()
        aW = self.aW()
        aP = self.aP()
        
        aE += self._Gamma / self._dx
        aW += self._Gamma / self._dx
        aP += aE + aW - self.Sp()
 
#        for i in range(self.__nvx):
#            aE[i] += self.__Gamma / self.__dx
#            aW[i] += self.__Gamma / self.__dx
#            aP[i] += aE[i] + aW[i]

if __name__ == '__main__':
    
    df1 = Diffusion1D(5, 5, 1)
    df1.alloc(5)
    df1.calcCoef()
    df1.source(100)

    print('-' * 20)  
    print(df1.aP(), df1.aE(), df1.aW(), df1.Su(), sep = '\n')
    print('-' * 20)  

    df1.bcDirichlet('LEFT_WALL', 2)
    df1.bcDirichlet('RIGHT_WALL', 1)
    print(df1.aP(), df1.aE(), df1.aW(), df1.Su(), sep = '\n')
    print('-' * 20)  
