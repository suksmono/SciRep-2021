"""
-------------------------------------------------------------------------------
Created on April,27, 2020
author: Suksmono@{STEI-ITB, MDR Inc.}
Find a Hadamard matrix using QA Computer/DWave
\-------------------------------------------------------------------------------
"""
# from sympy import *
# #from prb_symbolic import *
# from prb_williamson import *
import neal
import numpy as np
###
import pandas as pd
import matplotlib.pyplot as plt
##
# from prb_small_williamson import *
from prb_turyn import *
#
def writeHMatrix(fname, matx):
   # df=pd.DataFrame(data=matx.astype(int))
    df=pd.DataFrame(data=matx)
    df.to_csv(fname,sep=',', header=False, index=False)
    
def readHMatrix(fname):
    data=pd.read_csv(fname,sep=',', header=None )
    # return data.as_matrix() #.tolist()
    return data.values #.tolist()
#
def construct_XYZW_normal(v):
    # corr are taken from ABCD, not XYZW
    Nr=len(v)
    N=int((Nr+11)/4)
#    print('N=',N)
    #-- init XYZW
    X=np.zeros([N,1])
    Y=np.zeros([N,1])
    Z=np.zeros([N,1])
    W=np.zeros([N-1,1])
    #
    stX=0      
    edX=N-4
    #
    stY=edX #+1
    edY=2*N-7
    stZ=edY#+1
    edZ=3*N-9
    stW=edZ#+1
    edW=4*N-1-10
#    print('stX:',stX,'stY:',stY,'stZ:',stZ,'stW:',stW)
#    print('edX:',edX,'edY:',edY,'edZ:',edZ,'edW:',edW)
    #    
    if (edX>0):
        tX=v[stX: edX] # ;%get ss0(1) of  %X=[1; 1; 1; -1];
    else:
        tX=[]
    #--
    # --
    tY=v[stY:edY] #; %;  % get ss0(1) of Y=[1; ss0(0); -ss0(0); -1];
    tZ=v[stZ:edZ] #; %; % Z=[1; ss0(1); ss0(2); 1];
    tW=v[stW:edW] #; %;      %W =[1; ss0(3); ss0(4)];
    #
    if edX>0:
        NX=len(tX)
        for m in range(NX):
            X[m+2]=tX[m] #; !!!!!!!!
    else:
        tX=[]
    #--
    NY=len(tY)
    for m in range(NY):
        Y[m+1]=tY[m]
    #
    Y[N-2]= -Y[1] #enforcing constraint to Y
#    %
    
    NZ=len(tZ);
    for m in range(NZ):
        Z[m+1]=tZ[m]
#    #
    NW=len(tW)
    for m in range(NW):
        W[m+1]=tW[m]
    #
#    %---
    X[0]= 1
    Y[0]= 1
    Z[0]= 1
    W[0]= 1
    #
    X[N-1]= -1
    Y[N-1]= -1
    Z[N-1]= 1
    #
    X[1] = 1
    X[N-2]=1 
    #
    return(X,Y,Z,W)
###
def contruct_XXYY(vSol, X,Y,N_XTD):
    N_ORG=len(X)
    sSol=vSol[0:NQ0]
    N_ORG05=int(N_ORG/2)  
    # X->XX
    XX=np.zeros(N_ORG+N_XTD)
    XX[0:N_ORG05]=X[0:N_ORG05] #[x0,x1,*...]
    XX[N_ORG05:N_ORG05+N_XTD]=sSol[0:N_XTD]
    XX[N_ORG05+N_XTD:N_ORG+N_XTD]=X[N_ORG05:N_ORG]
    # Y->YY
    YY=np.zeros(N_ORG+N_XTD)
    YY[0:N_ORG05]=Y[0:N_ORG05] #[x0,x1,*...]
    YY[N_ORG05:N_ORG05+N_XTD]=sSol[N_XTD:2*N_XTD]
    YY[N_ORG05+N_XTD:N_ORG+N_XTD]=X[N_ORG05:N_ORG]
    return(XX,YY)
    
# ---
if __name__ == '__main__': 
    # order definition
    N_ORG=4
    N_XTD=8
    #--
    N=N_ORG+N_XTD # final order 
    #    P=N-1
    M=4*(3*N-1) # order of H-matrix
    #
    #    Nr=4*N-11 # length of random vector
    #
    NQ0=2*N_XTD
    idx=0
    ## 4->8: 4,7,18 corret
    idx_f=3 # idx file to read (+1 from the displaye idx in step-6)
    
    # --
    # read an entry of filled zw and incomplete xy
    path_name0='..\\py_src\\DATA\\'
    zzww_name= 'XYZW_star_xy_zzww_fill_'+str(N_ORG)+'_'+str(N_ORG+N_XTD)+ '.txt'
    fname_zzww=path_name0+zzww_name
    #
    in_file = open(fname_zzww, "r")
    #num_lines = sum(1 for line in open(fname_zzww))
    # read 
    for m in range(idx_f):
        t_hex=in_file.readline()
#        print(t_hex)
        m+=1
    #--
    in_file.close()
    # convert into spin
    NHO=int(np.ceil(N_ORG/4)) #hex length of  originalspin
    NHX=int(np.ceil((N_XTD+N_ORG)/4)) #hex length of  extended spin
    #
    hX=t_hex[0:NHO]
    dX=int(hX,16)
    X=int_to_spin(dX,N_ORG)
    #
    hY=t_hex[NHO:NHO+1]
    dY=int(hY,16)
    Y=int_to_spin(dY,N_ORG)
    #
    hZZ=t_hex[2*NHO:2*NHO+NHX]
    dZZ=int(hZZ,16)
    ZZ=int_to_spin(dZZ,N_ORG+N_XTD)
    #
    hWW=t_hex[2*NHO+NHX:2*NHO+2*NHX]
    dWW=int(hWW,16)
    WW0=int_to_spin(dWW,N_ORG+N_XTD)
    WW=WW0[1:N_ORG+N_XTD]
        
    # ------------------------------------------------------------
    #    number of required qubits
    # ------------------------------------------------------------    
    
    print('Finding Turyn Goethals-Siedel matrix of order:', M)
    print('Reading Ising coefficients from file ...')
    #    hi_complete_XY_4_2_1
    path_name='..\\Jij_DATA\\'
    Jij_name=path_name+'Jij_complete_XY_'+str(N_ORG)+\
    '_'+str(N_XTD)+ '_'+str(idx_f)+ '.txt'
    hi_name= path_name+'hi_complete_XY_'+str(N_ORG)+'_'+str(N_XTD)+ '_'+str(idx_f)+ '.txt'
    b_name= path_name +'b_complete_XY_'+str(N_ORG)+'_'+str(N_XTD)+ '_'+str(idx_f)+ '.txt'
    #
    Jij=readHMatrix(Jij_name)
    hiA=readHMatrix(hi_name)
    hi=hiA[:,0]
    bA=readHMatrix(b_name)
    b=bA[0][0]

    #normalize
    maxCoeff=np.max([np.max(abs(hi)), np.max(abs(Jij)),b])
    hi=hi/maxCoeff
    Jij=Jij/maxCoeff
    # measure of energy minimum
    b=b/maxCoeff

    '''
    # -----------------------------------------------------------------------------
    #     Convert Ising coefficients into data of Neal's format
    # -----------------------------------------------------------------------------
    '''
    #in dictionary format
    h={0:0}
    J={(0,1):1}
    
    for m in range(0,len(hi)):
        h[m]=hi[m]
        for n in range (m+1,len(hi)):
            J[m,n]=Jij[m,n]
        # --
    # --- end of filling h and J 
        
    '''
    -----------------------------------------------------------------------------
    4. SOLVE THE PROBLEM
    -----------------------------------------------------------------------------
    '''
    #
    print('Solving the problem using Neal  ...')
    solver=neal.SimulatedAnnealingSampler()
    NSWEEPS=1*M*M #1*1*10*1000 #1000*1000
    NR=10*1000 #*1000 #*10 #25 #100
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
    XX,YY=contruct_XXYY(vSol, X,Y,N_XTD)
    H=construct_hadamard(XX,YY,ZZ,WW)
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
        print('Hadamard matrix is NOT found in minimum Energy...!!!')
    ## check all solutions
    NR,NC=aSol.shape
    idx=0
    vIdx=[]
    for m in range(NR):
        tt=aSol[m,:]
        vtt=tt#[0]
        #--
        XX,YY=contruct_XXYY(tt[0:NQ0], X,Y,N_XTD)
        H=construct_hadamard(XX,YY,ZZ,WW)
        if is_hadamard(H):
            vIdx.append(m)
            idx=idx+1
    #--
    if(len(vIdx)>0):
        print('Solutions found: ',len(vIdx))
        tSol=aSol[vIdx[0],:]
        vSol=tSol#[0]
        #--
        XX,YY=contruct_XXYY(vSol[0:NQ0], X,Y,N_XTD)
        H=construct_hadamard(XX,YY,ZZ,WW)
        DG=np.matmul(H, H.transpose())
        plt.imshow(DG.astype(float), cmap="gray")
        plt.show()
