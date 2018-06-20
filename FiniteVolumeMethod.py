#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  4 13:21:50 2018

@author: luiggi
"""

from Mesh import Mesh
from Coefficients import Coefficients
from Diffusion import Diffusion1D
from Advection import Advection1DUpwind
from Advection import Advection1DCDS
from Advection import Advection1DQuick
from Advection import Advection1DUpwind2
from Matrix import Matrix
from Temporal import Temporal1DBackEuler
import time

# def crono(f):
# 	"""
# 	Regresa el tiempo que toma en ejecutarse la funcion.
# 	"""
# 	def eTime():
# 		t1 = time.time()
# 		f()
# 		t2
# 		return 'Elapsed time: ' + str((t2 - t1)) + "\n"
# 	return eTime

def decorate(f):
    def nicePrint(**kargs):
        line = '-' * 70
        print('.'+ line + '.')
        print('|{:^70}|'.format('NoNacos : Numerical Objects for Natural Convection Systems'))
        print('.'+ line + '.')
        print('|{:^70}|'.format(' Ver. 0.1 Copyright LMCS 2018'))
        print('.'+ line + '.')
        f(**kargs)
        print('.'+ line + '.')
    return nicePrint
 
@decorate
def printData(**kargs):
	for (key,value) in kargs.items():
		print('|{:^70}|'.format('{0:>15s} = {1:10.5e}'.format(key, value)))


if __name__ == '__main__':
 
    Coefficients.alloc(5)
    m = Mesh(nodes = 5)
    d = Diffusion1D(m.volumes())
    ma = Matrix(m.volumes())
    a = Advection1DQuick(m.volumes())

    print(m.delta(), d.aP(), a.aP(), ma.mat(), sep='\n')

    printData(nvx =5, nx = 6, longitud = 1.3)

