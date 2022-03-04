# -*- coding: utf-8 -*-
"""
final step:fillx,y
atlernative method
"""
import pandas as pd
from prb_turyn import *
from prb_generate_XYZW import *
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
if __name__=='__main__':
    #read data
    ##################### ZZ, WW: BUG-->FIXED ################
    # original length
    N_ORG=4
    N_BIT_ORG=4*N_ORG
    # extension length
    N_XTD=2 #extend number-of-bits each
    #--- file definition ---
    path_name='DATA//'
    fname_xtd=path_name+'XYZW_star_xy_zzww_fill_'+\
        str(N_ORG)+'_'+str(N_ORG+N_XTD)+'.txt'
    ## prepare open file
    in_file = open(fname_xtd, "r")
    ##
    num_lines = sum(1 for line in open(fname_xtd))
    #-- from X*,Y*,Z*,W*; insert spin-vector and check
    NQ=2*N_XTD
    # 
    NFIN=4*(N_ORG+N_XTD) #NR+4*NXTD
    start_lag=int(NFIN/2)
    MAX_ITER=2**NQ
    cnt=0
    idx=0
    ## estimate number of row
    NROW=0
    for tv in in_file:
        NROW+=1
    #--
    print('Total row in file=', NROW)
    in_file.close()
    ##-- reopen
    in_file = open(fname_xtd, "r")
    
    NHO=int(np.ceil(N_ORG/4)) #hex length of  originalspin
    NHX=int(np.ceil((N_XTD+N_ORG)/4)) #hex length of  extended spin
    cntH=0
    cnt_row=0
    for t_hex in in_file:
        # convert into spin
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
        #
        for m in range(2**(2*N_XTD)):
            ts=int_to_spin(m,2*N_XTD)
            #print(ts)
            XX,YY=construct_XXYY(ts, X,Y,N_XTD)
            H=construct_hadamard(XX,YY,ZZ,WW)
            if is_hadamard(H):
                cntH +=1
                # print('H-matrix found')
                DG=np.matmul(H, H.transpose())
                # display indicator
                plt.imshow(DG.astype(float), cmap="gray")
                plt.show()
            #--
        #-
        cnt_row+=1
        print('cek',cnt_row,'of',NROW)
    # --
    if cntH>0:
        print('H-matrix found=', cntH)
   