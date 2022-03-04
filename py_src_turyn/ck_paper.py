"""
-------------------------------------------------------------------------------
Created on April,27, 2020
author: Suksmono@{STEI-ITB, MDR Inc.}
Find a Hadamard matrix using QA Computer/DWave
\-------------------------------------------------------------------------------
"""
from sympy import *
# #from prb_symbolic import *
# from prb_williamson import *
#import neal
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
def auto_corr_sym(x):
    #non-periodic autocorrelation
    N=len(x)
    ML=N#-1
    idx=0
    r=np.empty([ML,1], dtype=object)
    # r=np.zeros([ML,1])
    for s in range(ML):
        t=0
        for m in range (N-s): 
            t=t+x[m]*x[m+s]
        #
        r[idx,0]=t
        idx=idx+1
    return(r)
    #--  
def remove_id(NX,s):
    #
    for m in range(len(NX)):
        tsym=NX[m][0]
        # print('m=',m,'->',tsym)
        if len(tsym.free_symbols)>0:
            for n in range(len(s)):
                #print(n)
                tsym=tsym.subs({s[n]**2:1})
            #
            NX[m][0]=simplify(tsym)
        # endif
        
    # --
    return(NX)
#
def remove_equ(EX,s):
    #
    # for m in range(len(NX)):
    #     tsym=NX[m][0]
    #     # print('m=',m,'->',tsym)
    # if len(tsym.free_symbols)>0:
    EXS=EX
    for n in range(len(s)):
        #print(n)
        EXS=EXS.subs({s[n]**2:1})
    #        
    # --
    return(simplify(EXS))

def sum_autocorr(NX):
    R=0    
    LN=len(NX)
    for m in range(1,LN):
        item=NX[m][0]
        R=R+item
    return(R)
#
if __name__ == '__main__': 
    # order definition
    N=4 #4,6, ...even number
    P=N-1
    M=4*(3*N-1) # order of H-matrix
    #
    Nr=4*N-5 # length of random vector
    NQ=6
    #
    ss=symbols('s0:%d'%NQ)
    #
    ss0=np.asarray(ss)
    A=[1,ss0[2],ss0[3],1,1,ss0[4],ss0[5] ]
    B=[1,ss0[2],ss0[3],1,-1,-ss0[4],-ss0[5] ]
    C=[1, ss0[0], 1, -1]
    D=[1,ss0[1], -ss0[1],-1] 
    #
    NA=auto_corr_sym(A)
    NAS=remove_id(NA,ss0)
    RA=sum_autocorr(NAS)
    #
    NB=auto_corr_sym(B)
    NBS=remove_id(NB,ss0)
    RB=sum_autocorr(NBS)

    #
    NC=auto_corr_sym(C)
    # NCS=remove_id(NC,ss0)
    RC=sum_autocorr(NC)

    #
    ND=auto_corr_sym(D)
    #NDS=remove_id(ND,ss0)
    #
    RD=sum_autocorr(ND)
    #energy
    E0=RA+RB+RC+RD
    E1=expand(E0**2)
    E2=simplify(E1)
    E=remove_equ(E2,ss0)
    
    
    
    # ------------------------------------------------------------
    #    number of required qubits
    # ------------------------------------------------------------    
    # .... 
    # ------------------------------------------------------------
 
    # print('Finding Turyn Goethals-Siedel matrix of order:', M)
    # print('Reading Ising coefficients from file ...')
    # Jij_name='Jij_MATLAB'+str(M) + '.txt'
    # hi_name= 'hi_MATLAB'+str(M) + '.txt'
    # b_name= 'b_MATLAB'+str(M) + '.txt'
    #
    # Jij=readHMatrix(Jij_name)
    # hiA=readHMatrix(hi_name)
    # hi=hiA[:,0]
    # bA=readHMatrix(b_name)
    # b=bA[0][0]

    # #normalize
    # maxCoeff=np.max([np.max(abs(hi)), np.max(abs(Jij)),b])
    # hi=hi/maxCoeff
    # Jij=Jij/maxCoeff
    # # measure of energy minimum
    # b=b/maxCoeff

    # #     '''
    # #     -----------------------------------------------------------------------------
    # #     Convert Ising coefficients into data of Neal's format
    # #     -----------------------------------------------------------------------------
    # #     '''
    # #in dictionary format
    # h={0:0}
    # J={(0,1):1}
    
    # for m in range(0,len(hi)):
    #     h[m]=hi[m]
    #     for n in range (m+1,len(hi)):
    #         J[m,n]=Jij[m,n]
        
    # '''
    # -----------------------------------------------------------------------------
    # 4. SOLVE THE PROBLEM
    # -----------------------------------------------------------------------------
    # '''
    # #
    # print('Solving the problem using Neal  ...')
    # solver=neal.SimulatedAnnealingSampler()
    # NSWEEPS=1*1*10*1000 #1000*1000
    # NR=M#*1000 #*10 #25 #100
    # response = solver.sample_ising(h, J, sweeps=NSWEEPS, num_reads=NR)
    # #
    # vE=response.data_vectors['energy']
    # aSol=response.record['sample']
    # #
    # print('Configurations:\n', aSol)
    # print('Energy:\n',vE)
    # #
    # idxMinE=np.argmin(vE)
    # print('Minimum Energy:',vE[idxMinE])#, 'supposed to be', -b)
    # print('Minimum Configurations:',aSol[idxMinE])
    # tSol=aSol[idxMinE]
    # vSol=tSol#[0]
    # [X,Y,Z,W]=construct_XYZW(vSol[0:Nr])
    # H=construct_hadamard(X,Y,Z,W)
    # DG=np.matmul(H, H.transpose())
    # # display indicator
    # plt.imshow(DG.astype(float), cmap="gray")
    # plt.show()
    # # display H-matrix
    # plt.imshow(H.astype(float), cmap="gray")
    # plt.show()
    # #
    # #     if is_hadamard(W):
    # #         hname='../H_MATRIX/'+'H_'+str(M) + '.txt'
    # #         writeHMatrix(hname, W)


        