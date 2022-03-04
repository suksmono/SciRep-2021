# -*- coding: utf-8 -*-
"""
Created on Wed May  6 11:48:25 2020

@author: Asus
"""
#
import matplotlib.pyplot as plt
from prb_turyn import *
#
def generate_XYZW_4():
    #W=np.zeros([4,1])
    #
    x=[1,-1,1,-1]
    y=[1,-1,-1,-1]
    z=[1,-1,-1,1]
    w=[1,1,1]
    return(x,y,z,w)
#
def generate_XYZW_6():
    #W=np.zeros([4,1])
    #
    x=[1,-1,1,1,-1,1]
    y=[1,1,1,-1,-1,1]
    z=[1,-1,1,1,1,-1]
    w=[1,-1,-1,-1,-1]
    return(x,y,z,w)
		
def generate_XYZW_8():
    #W=np.zeros([4,1])
    #
    x=[1,1,-1,1,-1,1,-1,1]
    y=[1,-1,-1,-1,-1,-1,-1,1]
    z=[1,-1,-1,1,1,1,1,-1]
    w=[1,1,1,-1,1,1,-1]
    return(x,y,z,w)
#
if __name__ == '__main__': 
#    X,Y,Z,W=generate_XYZW_8()
    X,Y,Z,W=generate_XYZW_4()
    H=construct_hadamard(X,Y,Z,W)
    DG=np.matmul(H, H.transpose())
    # display indicator
    plt.imshow(DG.astype(float), cmap="gray")
    plt.show()
    # display H-matrix
    plt.imshow(H.astype(float), cmap="gray")
    plt.show()
    if is_hadamard(H):
        print('Hadamard matrix is found ...!!!')
    else:
        print('Hadamard matrix is NOT found ...!!!')