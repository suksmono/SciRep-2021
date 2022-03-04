"""
-------------------------------------------------------------------------------
Created on April, 2020
author: Suksmono@{STEI-ITB, MDR Inc.}
Find a Hadamard matrix using QA Computer/DWave
#2020/April: cek williamson matrix
H = [  A   B   C   D ]
    [ -B   A  -D   C ]
    [ -C   D   A  -B ]
    [ -D  -C   B   A ]
where A, B,C, D is a symmetric circulant matrix
eg. C = [c0 c1    ... c2 c1]
        [c1 c0  c1 ... c2]
eg. 3x3 => [a b b]
           [b a b]
           [b b a]

For a 4kx4k H-matrix, we need 4*floor( (k-1)/2 +1) variables  

For a 12x12 H-matrix, we need 4*floor(( ((12/4)-1)/2) +1)=4x2= 8 variables  
A=> {q0, q1}---> [ [q0 q1 q1]; [q1 q0 q1]; [q1 q1 q0] ]
B=> {q2, q3}---> [ [q2 q3 q3]; [q3 q2 q3]; [q3 q3 q2] ]
C=> {q4, q5}---> [ [q4 q4 q5]; [q5 q4 q5]; [q5 q5 q4] ]
D=> {q6, q7}---> [ [q6 q7 q7]; [q7 q6 q7]; [q7 q7 q6] ]
-------------------------------------------------------------------------------
"""
from sympy import *
#from prb_symbolic import *
from prb_williamson import *
import neal
import numpy as np
###
import pandas as pd
import matplotlib.pyplot as plt
    
if __name__ == '__main__': 
    # order definition
    K =5# order = 4*K
    M=4*K
    K05=int(K/2)+1
    
    # ------------------------------------------------------------
    #    number of required qubits
    # ------------------------------------------------------------    
    # each of (first row of) A,B,C,D submatrix needs K05+1
    # extra 1 because the middle term repeated => 0,1,1; 0,1,2,2,1
    # 0,1,2, ..., K05-1, K05-1, K05-2, ....2,1 
    # ------------------------------------------------------------
    NQ0=4*K05 # no.of main qubits
#    NQ=NQ0 + int(NQ0*(NQ0-1)/2)
    NQ=int(1.5*NQ0) # + int(NQ0*(NQ0-1)/2)
    # 
#    Jij_name = strcat(strcat('..\\Jij_DATA\\Jij_MATLAB', int2str(M)),'.txt');
#    hi_name = strcat(strcat('..\\Jij_DATA\\hi_MATLAB', int2str(M)),'.txt');
#    b_name = strcat(strcat('..\\Jij_DATA\\b_MATLAB', int2str(M)),'.txt');

    print('Finding Williamson matrix of order:', M)
    print('Reading Ising coefficients from file ...')
    Jij_name='../Jij_DATA/'+'Jij_MATLAB'+str(M) + '.txt'
    hi_name='../Jij_DATA/'+'hi_MATLAB'+str(M) + '.txt'
    b_name='../Jij_DATA/'+'b_MATLAB'+str(M) + '.txt'
    Jij=readHMatrix(Jij_name)
    hi=readHMatrix(hi_name)
    b=readHMatrix(b_name)

    #normalize
    maxCoeff=np.max([np.max(abs(hi)), np.max(abs(Jij)),b])
    hi=hi/maxCoeff
    Jij=Jij/maxCoeff
    # measure of energy minimum
    b=b/maxCoeff

    '''
    -----------------------------------------------------------------------------
    Convert Ising coefficients into data of Neal's format
    -----------------------------------------------------------------------------
    '''
    #in dictionary format
    h={0:0}
    J={(0,1):1}
    
    for m in range(0,len(hi)):
        h[m]=hi[m]
        for n in range (m+1,len(hi)):
            J[m,n]=Jij[m,n]
        
    '''
    -----------------------------------------------------------------------------
    4. SOLVE THE PROBLEM
    -----------------------------------------------------------------------------
    '''
    #
    print('Solving the problem using Neal  ...')
    solver=neal.SimulatedAnnealingSampler()
    NSWEEPS=1*10*10*1000 #1000*1000
    NR=M*M#*1000 #*10 #25 #100
    response = solver.sample_ising(h, J, sweeps=NSWEEPS, num_reads=NR)
    #
    vE=response.data_vectors['energy']
    aSol=response.record['sample']
    #
    print('Configurations:\n', aSol)
    print('Energy:\n',vE)
    #
    idxMinE=np.argmin(vE)
    print('Minimum Energy:',vE[idxMinE])#, 'supposed to be', -b)
    print('Minimum Configurations:',aSol[idxMinE])
    tSol=aSol[idxMinE]
    vSol=tSol#[0]
    #
    # --- construct H-matrix --
    # you can replace symbols with numeric
    tt=vSol[0:NQ0]
    tA,tB,tC,tD=construct_ABCD(tt)
    W=construct_symm_williamson(tA,tB,tC,tD)
    DG=np.matmul(W, W.transpose())
#    print('W-Matrix \n',W)
#    print('Indicator matrix \n',DG)
    #
    if is_hadamard(W):
        hname='../H_MATRIX/'+'H_'+str(M) + '.txt'
        writeHMatrix(hname, W)
    # save if W a Hadamard matrix

    # display
    plt.imshow(DG.astype(float))
#    plt.imshow()

        