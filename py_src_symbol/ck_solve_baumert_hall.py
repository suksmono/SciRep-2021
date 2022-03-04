"""
-------------------------------------------------------------------------------
Created on April, 2020
author: Suksmono@{STEI-ITB, MDR Inc.}
Find a Hadamard matrix using QA Computer/DWave
    --- Baumert-Hall ---
    H=  [ A  A  A  B  -B  C -C -D   B  C -D -D]
        [ A -A  B -A  -B -D  D -C  -B -D -C -C]
        [ A -B -A  A  -D  D -B  B  -C -D  C -C]
        [ B  A -A -A   D  D  D  C   C -B -B -C]  
        
        [ B -D  D  D   A  A  A  C  -C  B -C  B] 
        [ B  C -D  D   A -A  C -A  -D  C  B -B]
        [ D -C  B -B   A -C -A  A   B  C  D -D] 
        [-C -D -C -D   C  A -A -A  -D  B -B -B]  
        
        [ D -C -B -B  -B  C  C -D   A  A  A  D]  
        [-D -B  C  C   C  B  B -D   A -A  D -A] 
        [ C -B -C  C   D -B -D -B   A -D -A  A] 
        [-C -D -D  C  -C -B  B  B   D  A -A -A]   
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
from prb_williamson import *
from prb_baumert_hall import *
import neal
import numpy as np
###
import pandas as pd
import matplotlib.pyplot as plt
   
def read_matrix(fname):
    data=pd.read_csv(fname,sep=',', header=None )
    return data.values #.tolist()
   
if __name__ == '__main__': 
    # order definition
    K =3# order = 4*K
    M=12*K
    K05=int(K/2)+1
    
    # ------------------------------------------------------------
    #    number of required qubits
    # ------------------------------------------------------------    
    # each of (first row of) A,B,C,D submatrix needs K05+1
    # extra 1 because the middle term repeated => 0,1,1; 0,1,2,2,1
    # 0,1,2, ..., K05-1, K05-1, K05-2, ....2,1 
    # ------------------------------------------------------------
    NQ0=4*K05 # no.of main qubits
    NQ=NQ0 + int(NQ0*(NQ0-1)/2)
    NQ=int(1.5*NQ0)
    # 
    print('Finding Baumert-Hall matrix of order:', M)
    print('Reading Ising coefficients from file ...')
#    fname='../Jij_DATA/'+'Jij_H'+str(M) + '.txt'
    Jij_name='../Jij_DATA/'+'Jij_BH_MATLAB'+str(M) + '.txt'
    hi_name='../Jij_DATA/'+'hi_BH_MATLAB'+str(M) + '.txt'
    b_name='../Jij_DATA/'+'b_BH_MATLAB'+str(M) + '.txt'
    Jij=read_matrix(Jij_name)
    hi=read_matrix(hi_name)
    b=read_matrix(b_name)

    #normalize
    maxCoeff=np.max([np.max(abs(hi)), np.max(abs(Jij))])
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
    NSWEEPS = 1*10*1000 # M*NQ  # 1*1*10*1000 #1000*1000
    NR= M*M #10*1000      #M #M#*1000 #*10 #25 #100
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
#    W=construct_symm_williamson(tA,tB,tC,tD)
    W=construct_baumert_hall(tA,tB,tC,tD)
    DG=np.matmul(W, W.transpose())
    #
    # save if W a Hadamard matrix
    if is_hadamard(W):
        print('Hadamard matrix is found ...!')
        hname='../H_MATRIX/'+'H_BH_'+str(M) + '.txt'
        writeHMatrix(hname, W)
    else:
        print('Hadamard matrix is NOT found ...!')
    #        
                  
    plt.imshow(DG.astype(float), cmap='gray')
    plt.show()
    plt.imshow(W.astype(float),  cmap='gray')
    plt.show()

        