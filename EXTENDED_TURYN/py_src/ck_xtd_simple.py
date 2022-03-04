# -*- coding: utf-8 -*-
"""
Created on Fri May 22 01:35:49 2020

@author: Alvin F
"""
import pandas as pd
from prb_turyn import *
from prb_generate_XYZW import *
from prb_london import *
##
#--
def construct_XXYY(sSol, X,Y,N_XTD):
    N_ORG=len(X)
    # sSol=vSol[0:NQ0]
    N_ORG05=int(N_ORG/2)  
    # N_XTD05=int(N_XTD/2)
    # X->XX
    XX=np.zeros(N_ORG+N_XTD)
    XX[0:N_ORG05]=X[0:N_ORG05] #[x0,x1,*...]
    XX[N_ORG05:N_ORG05+N_XTD]=sSol[0:N_XTD]
    XX[N_ORG05+N_XTD:N_ORG+N_XTD]=X[N_ORG05:N_ORG]
    # Y->YY
    YY=np.zeros(N_ORG+N_XTD)
    YY[0:N_ORG05]=Y[0:N_ORG05] #[x0,x1,*...]
    YY[N_ORG05:N_ORG05+N_XTD]=sSol[N_XTD:2*N_XTD]
    YY[N_ORG05+N_XTD:N_ORG+N_XTD]=Y[N_ORG05:N_ORG]
    return(XX,YY)
##

# x=[ 1,+1, 1,-1,-1,-1]
# y=[ 1,-1, 1,+1,-1,-1] 
# z=[ 1,+1, 1,-1, 1,+1]
# w=[ 1,+1,-1,+1,-1]

if __name__=='__main__':
    X=[ 1,+1,-1,-1]
    Y=[ 1,-1,-1,-1] 
    ZZ=[ 1,+1, 1,-1,+1,+1]
    WW=[ 1,+1,-1,+1,-1]
    #
    N=len(ZZ)
    N_XTD=2
    cntf=0
    for m in range(2**(2*N_XTD)):
        print('try',m,'of',2**(2*N_XTD))
        #
        ts=int_to_spin(m,2*N_XTD)
        #print(ts)
        XX,YY=construct_XXYY(ts, X,Y,N_XTD)
        H=construct_hadamard(XX,YY,ZZ,WW)
        if is_hadamard(H):
            print('H-matrix found')
            DG=np.matmul(H, H.transpose())
            # display indicator
            plt.imshow(DG.astype(float), cmap="gray")
            plt.show()
            cntf+=1
        #--
    #---
    print('found H:', cntf)
    #
    tbl_cos=gen_tbl_cos_t(N)
    zeta_ohmega={0,-1,1,-2,2,-3,3,-4,4}
    #
    #--- convert to array    
    Z=np.zeros([N,1])
    Z[:,0]=ZZ[0:N]
    #
    W=np.zeros([N-1,1])
    W[:,0]=WW[0:N-1]
    #
    TF_Z,TF_ZO, fz=is_zw_qualified(Z, zeta_ohmega, tbl_cos)
    TF_W,TF_WO, fw=is_zw_qualified(W, zeta_ohmega, tbl_cos)
    #
    print('3N-1 -> ', 3*N-1)
    print('vQZ=',max(max(fz)),'TF_Z:', TF_Z,'TF_ZO:', TF_ZO)
    print('vQW=',max(max(fw)),'TF_W:', TF_W,'TF_WO:', TF_WO)
    print('mean(vQW+vQZ)=',np.mean(np.mean(fz+fw)))
    print('max(vQW+vQZ)=',np.max(np.max(fz+fw)))


    