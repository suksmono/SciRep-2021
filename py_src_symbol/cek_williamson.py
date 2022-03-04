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
import pandas as pd
## define functions
#def writeHMatrix(fname, matx):
#   # df=pd.DataFrame(data=matx.astype(int))
#    df=pd.DataFrame(data=matx)
#    df.to_csv(fname,sep=',', header=False, index=False)
#    
#def readHMatrix(fname):
#    data=pd.read_csv(fname,sep=',', header=None )
#    return data.as_matrix() #.tolist()
##
#def construct_symm_williamson(A,B,C,D):
#    """
#    input  : submatrix A,B,C,D
#    output : williamson (hadamard) matrix
#    --- williamson ---
#    H = [  A   B   C   D ]
#        [ -B   A  -D   C ]
#        [ -C   D   A  -B ]
#        [ -D  -C   B   A ]
#    """
#    M,N=A.shape
#    W=np.empty([4*M,4*N], dtype=object)
#    '''
#    # row-1: [  A   B   C   D ]
#    '''
#    W[0:M,0:M]=     A[0:M,0:M]
#    W[0:M,M:2*M]=   B[0:M,0:M]
#    W[0:M,2*M:3*M]= C[0:M,0:M]
#    W[0:M,3*M:4*M]= D[0:M,0:M]
#
#    # row-2: [ -B   A  -D   C ]
#    W[M:2*M,0:M] =      -B[0:M,0:M]
#    W[M:2*M,M:2*M] =    A[0:M,0:M]
#    W[M:2*M,2*M:3*M] = -D[0:M,0:M]
#    W[M:2*M,3*M:4*M] =  C[0:M,0:M]
#
#    '''
#    # row-3: [ -D  -C   B   A ]
#    '''
#    W[3*M:4*M,0:M] =    -D[0:M,0:M]
#    W[3*M:4*M,M:2*M] =  -C[0:M,0:M]
#    W[3*M:4*M,2*M:3*M]= B[0:M,0:M]
#    W[3*M:4*M,3*M:4*M]= A[0:M,0:M]
#   
#    '''
#    # row-3: [ -C   D   A  -B ]
#    '''
#    W[2*M:3*M,0:M] =    -C[0:M,0:M]
#    W[2*M:3*M,M:2*M] =  D[0:M,0:M]
#    W[2*M:3*M,2*M:3*M]= A[0:M,0:M]
#    W[2*M:3*M,3*M:4*M]= -B[0:M,0:M]
#    #
#    return(W)
##
#def calculate_Hks(D):
#    '''
#    given D, calculate Hks
#    '''
#    M,N=D.shape
#    Hks=0
#    for m in range(0,M):
#        for n in range(m+1,M):
#            #print('m=',m,'n=',n)
#            t=expand(D[m,n]**2)
#            Hks=Hks+t
#    #
#    return(Hks)
##
#def remove_squares(Hks,ss):
#    '''
#    # substitute: s**2->1, s**3->s, s**4->1
#    '''
#    for m in range(len(ss)):
#        Hks= Hks.subs( {ss[m]**2:1})
#        Hks= Hks.subs( {ss[m]**3:ss[m]})
#        Hks= Hks.subs( {ss[m]**4:1}) 
#    return(simplify(Hks))
##
##
#def Hks_to_Hkq(Hks,ss,qq):
#    '''
#    # substitute: s**2->1, s**3->s, s**4->1
#    '''
#    Hkq=Hks
#    for m in range(len(ss)):
#        Hkq= Hkq.subs( {ss[m]:qq[m]})
#    return(simplify(Hkq))
#
##    #
#def gen_subs_pair(qq,M):
#    NR=int(M*(M-1)/2)
#    tbl_q=np.empty([NR,3], dtype=object)
#    idx=0
#    for m in range(M):
#        for n in range(m+1,M):
#            # print('idx:', idx+M, ':', qq[m],qq[n],'->',qq[idx+M])
#            tbl_q[idx,0]=qq[m]
#            tbl_q[idx,1]=qq[n]
#            tbl_q[idx,2]=qq[idx+M]
#            #
#            idx=idx+1
#        #--
#    #--
#    return(tbl_q)
#
##
#def C_boole(qi,qj,qk,dij):
#    '''
#    Boolean reduction qi*qj <- qk + C(qi,qj,qk;dij)
#    >> this one do: C=dij*(3qk + qi*qj - 2qi*qk-2qj*qk)
#    '''
#    return(dij*(3*qk + qi*qj -2*qi*qk - 2*qj*qk))
###
#def Hkq_to_H2q(Hkq,spair,delta):
#    H2q=Hkq
#    for tsp in spair:
#        #print(tsp[0],'*',tsp[1], '->' ,tsp[2])
#        H2q=H2q.subs({ \
#                   tsp[0]*tsp[1]:tsp[2]\
#                          }) \
#                + C_boole(tsp[0], tsp[1], tsp[2],delta)
#    # 
#    print('Simplifying H2q ...')
#    H2q=simplify(H2q)
#    #
#    return(H2q)
##
#def H2q_to_H2s(H2q,qq0,ss0):
#    H2s=H2q
#    N=len(qq0)
#    for m in range(N):
#        H2s = H2s.subs( {qq0[m]:q2s(ss0[m])} )
#    #
#    print('Simplifying H2s ...')
#    H2s=simplify(H2s)
#    #
#    return(H2s)
###
#def construct_ABCD(v):
#    #
#    """
##    A=> {q0, q1}---> [ [q0 q1 q1]; [q1 q0 q1]; [q1 q1 q0] ]
##    B=> {q2, q3}---> [ [q2 q3 q3]; [q3 q2 q3]; [q3 q3 q2] ]
##    C=> {q4, q5}---> [ [q4 q5 q5]; [q5 q4 q5]; [q5 q5 q4] ]
##    D=> {q6, q7}---> [ [q6 q7 q7]; [q7 q6 q7]; [q7 q7 q6] ]
#    #
#    ## USE: np.roll(vector,2)
#    """
#    vA=[ v[0], v[1], v[1]]
#    vB=[ v[2], v[3], v[3]]
#    vC=[ v[4], v[5], v[5]]
#    vD=[ v[6], v[7], v[7]]
#    
#    N=len(vA)
#    #init
#    A=np.empty([N,N], dtype=object)
#    B=np.empty([N,N], dtype=object)
#    C=np.empty([N,N], dtype=object)
#    D=np.empty([N,N], dtype=object)
#    #
#    A[:,0]=vA
#    B[:,0]=vB
#    C[:,0]=vC
#    D[:,0]=vD
#    for m in range(1,N):
#        A[:,m] = np.roll(vA,m)
#        B[:,m] = np.roll(vB,m)
#        C[:,m] = np.roll(vC,m)
#        D[:,m] = np.roll(vD,m)
#    #
#    return(A,B,C,D)

######
    

if __name__ == '__main__': 
    #
    K =3
    K05=int(K/2)+1
    M=4*(K05+1) 
    #
    F_SIMULATE=1
    # define matrix order
    # number of required qubits
    NQ0 = 4*(int(M/8)+1) # #M + int(M*(M-1)/2)
    NQ=NQ0 + int(NQ0*(NQ0-1)/2)
    
    
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
                 s-DOMAIN
      s0 s4 s8  s12  | s16 s20 s24 s28 s32 s36
      s1 s5 s9  s13  | s17 s21 s25 s29 s33 s37
      s2 s6 s10 s14  | s18 s22 s26 s30 s34 s38
      s3 s7 s11 s15  | s19 s23 s27 s31 s35 s39
    -------------------------------------------------------------------------------
    --- williamson ---
    H = [  A   B   C   D ]
        [ -B   A  -D   C ]
        [ -C   D   A  -B ]
        [ -D  -C   B   A ]

    """
#    NQ=8
    ss=symbols('s0:%d'%NQ)
    #
    ss0=np.asarray(ss)
#    s=ss0.reshape(int(NQ/M),M).tolist()
    """
    A=> {q0, q1}---> [ [q0 q1 q1]; [q1 q0 q1]; [q1 q1 q0] ]
    B=> {q2, q3}---> [ [q2 q3 q3]; [q3 q2 q3]; [q3 q3 q2] ]
    C=> {q4, q5}---> [ [q4 q5 q5]; [q5 q4 q5]; [q5 q5 q4] ]
    D=> {q6, q7}---> [ [q6 q7 q7]; [q7 q6 q7]; [q7 q7 q6] ]

    """
    #
    A,B,C,D=construct_ABCD(ss[0:NQ0])    
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
    

    '''    
    # normalize coefficients
    maxCoeff=np.max([np.max(abs(hi)), np.max(abs(Jij))])
    hi=hi/maxCoeff
    Jij=Jij/maxCoeff
    #
    b=b/maxCoeff
    '''
    #
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
    #fname='isingCoeff_H'+str(M) + '.txt'
    fname='../Jij_DATA/'+'Jij_H'+str(M) + '.txt'
    writeHMatrix(fname, Jij)
    
    F_SIMULATE=1    
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
        NSWEEPS=1*5*10*10*1000 #1000*1000
        NR=25 #100
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

        