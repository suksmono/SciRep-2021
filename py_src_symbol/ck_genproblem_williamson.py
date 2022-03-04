"""
-------------------------------------------------------------------------------
Created on April, 2020
author: Suksmono@{STEI-ITB, MDR Inc.}
Find a Hadamard matrix using QA Computer/DWave
#2020/April: Symmetric Williamson Matrix
H = [  A   B   C   D ]
    [ -B   A  -D   C ]
    [ -C   D   A  -B ]
    [ -D  -C   B   A ]
where A, B,C, D is a symmetric circulant matrix
eg. C = [c0 c1    ... c2 c1]
        [c1 c0  c1 .. c3 c2]
eg. 3x3 => [a b b]
           [b a b]
           [b b a]
A=> {q0, q1}---> [ [q0 q1 q1]; [q1 q0 q1]; [q1 q1 q0] ]
B=> {q2, q3}---> [ [q2 q3 q3]; [q3 q2 q3]; [q3 q3 q2] ]
C=> {q4, q5}---> [ [q4 q5 q5]; [q5 q4 q5]; [q5 q5 q4] ]
D=> {q6, q7}---> [ [q6 q7 q7]; [q7 q6 q7]; [q7 q7 q6] ]
-------------------------------------------------------------------------------
"""
#from sympy import *
#from symengine import *
import sympy
#from prb_symbolic import *
from prb_williamson import *
import neal
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == '__main__': 
    F_SIMULATE=0
    # order definition
    K =3
    M=4*K # order of H-matrix
    print('Hadamard matrix of William type, order:', M)
    '''
    #
    # --------------------------
    # number of required qubits
    # --------------------------
    
    # each of (first row of) A,B,C,D submatrix needs K05+1
    # extra 1 because the middle term repeated => 0,1,1; 0,1,2,2,1
    # 0,1,2, ..., K05-1, K05-1, K05-2, ....2,1 
    '''
    K05=int(K/2)+1 # 3: 0,1,1,0-> {0,1}, 5: 0,1,2,2,1-> {0,1,2}
    NQ0=4*K05 # no.of main qubits
    NQ=NQ0 + int(NQ0*(NQ0-1)/2) # main qubits + ancillaries 
    
    '''
    -------------------------------------------------------------------------------
    1. Formulate Hks
    -------------------------------------------------------------------------------
    '''
    """
    for 4x4 case, the arrangement of the qubits is
    -------------------------------------------------------------------------------
     |<---problem--->|<------ ancillas ------>|
    -------------------------------------------------------------------------------
    """
    ss=symbols('s0:%d'%NQ)
    #
    ss0=np.asarray(ss)
    #
    v_qmain=ss[0:NQ0]
    A,B,C,D=construct_ABCD(v_qmain)    
    H=construct_symm_williamson(A,B,C,D)
    #
    DS=np.matmul(H.transpose(),H)
    #
    Hks1 = calculate_Hks(DS)
    print('Hks-original\n',Hks1)
    # remove ID
    Hks = remove_squares(Hks1,ss)
    print('\nHks-simplified\n',Hks)
    
    # -- s->q
    qq=symbols('q0:%d'%NQ)
    qq0=np.asarray(qq)
    #
    '''
    -------------------------------------------------------------------------------
    2. Transform Hkq -> H2q -> H2s
    -------------------------------------------------------------------------------
    ''' 
    spair=gen_subs_pair(qq,NQ0)
    #
    [NPAIR, XX]=np.shape(spair)
    dijMax=M*M #**2
    #HMax=NPAIR*dijMax
    HMax=dijMax
    delta=4*HMax #6*16 #2*NPAIR*M #2*(M**2)
    print('Transform: Hkq->H2q->H2s ...')
    Hkq=Hks_to_Hkq(Hks,ss,qq)
    #
    H2q=Hkq_to_H2q(Hkq, spair,delta)
    H2s=H2q_to_H2s(H2q, qq0, ss0)
    print('H2s=', H2s)
    
    ''' check all variables --
    # v_list=H2s.free_symbols
    '''
    #
    ''' 
    ------------------------------------------------------------
    3. EXTRACT ISING PARAMETERS FROM SYMBOLIC SOLUTION
    ------------------------------------------------------------
    '''
    print('Obtaining Ising coefficients ...')
    b, hi, Jij = isingCoeffs(H2s,NQ)
    

        
    # normalize coefficients
    maxCoeff=np.max([np.max(abs(hi)), np.max(abs(Jij))])
    hi=hi/maxCoeff
    Jij=Jij/maxCoeff
    # measure of energy minimum
    b=b/maxCoeff
    #    
    '''
    -----------------------------------------------------------------------------
    3.b. save Ising Coefficients
    embed hi in the diagonal matrix of Jij
    -----------------------------------------------------------------------------
    '''
    for m in range(0,NQ):
        Jij[m,m]=hi[m]
    # save into file
    fname='../Jij_DATA/'+'Jij_H'+str(M) + '.txt'
    writeHMatrix(fname, Jij)
    
    if F_SIMULATE==1:
        # continue with simulation
        '''
        -----------------------------------------------------------------------------
        convert the problem into Ising coefficients
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
        select a solver
        > dimod: ExaxtSolver
        > neal:  SimulatedAnnealingSampler
        '''
        #
        print('Solving the problem using neal  ...')
        solver=neal.SimulatedAnnealingSampler()
        NSWEEPS=1*1*100*10*1000 #1000*1000
        NR=50 #100
        response = solver.sample_ising(h, J, sweeps=NSWEEPS, num_reads=NR)
        #
        vE=response.data_vectors['energy']
        #        aSol=response.samples_matrix
        aSol=response.record['sample']
        #
        print('Configurations:\n', aSol)
        print('Energy:\n',vE)
        #
        idxMinE=np.argmin(vE)
        print('Minimum Energy:',vE[idxMinE], 'supposed to be', -b)
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
        print('W-Matrix \n',W)
        print('Indicator matrix \n',DG)
        # display
        plt.imshow(DG.astype(float), cmap="gray")
        plt.imshow()

        