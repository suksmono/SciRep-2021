# -*- coding: utf-8 -*-
"""
Created on Tue May 19 19:13:27 2020
@author: Alvin F
# assume n=6
then x^2+y^2+ 2(z^2+w^2) =6*6-2 = 34
list of sqr numbers<=32: {0,1,4,9,16,25} ---rightterms^2<=18 {0,1,4,9,16}
--> {1,9,36}-> {1,9,2*9,2*9} => zw={3,3}
--> {6,3;1,0} -->zw={0,1}
xx--> {1^2,2^2;2*4^2,2*1^2} -->{12|34}
--------------------------------------------------------   
    x2      y2   2*z2          2*w         ttl   {z,w}
--------------------------------------------------------   
    0       0     2*4^2=32      2*1^2=2     34  {4,1}
    0       0     2*1^2=2       2*4^2=32    34  {1,4}
    3^2     5^2     0           0           34  {0,0}
    0       4^2     2*2^2=8      2*2^2=8    34  {2,2}
    4^2     0     2*2^2=8      2*2^2=8      34  {2,2}
    0       4     2*3^2        0            34  {3,0}  
    >> zeta={0,-1,1,-2,2,-4,4}
--------------------------------------------------------  
"""

import pandas as pd
from prb_turyn import *
from prb_generate_XYZW import *
##
#--
def gen_tbl_cos_t(N):
    #-generate table cos(j*theta), for theta= k*pi/100, k=1,..,100
    # so,it is a 2D array of [j,k]; j= 1, ..,(n-1)
    # N:number of element in X, or Y, ...
    # ---------------------------------------------------------------------
    #  lag theta->  0      1            2          ...    100  
    #   \/       unused  t1=1*pi/100  t2=2*pi/100 ... t100=100*pi/100   
    # ---------------------------------------------------------------------
    #   0           x       1           1       ...     1
    #   1           x    cos(1*t1)  cos(1*t2)            cos(1*t100)  
    #   2           x    cos(2*t1)  cos(2*t2)            cos(2*t100)  
    # ....
    #  N-1          x    cos((N-1)*t1) cos((N-2)*t2)        cos(1**t100)
    # ---------------------------------------------------------------------
    M=100
    tbl_cos_t=np.zeros([N+1,M+1])
    for j in range(1,N+1):     #lag: 1..N-1
        for k in range(1,M+1): #theta: 1 ... 100 *pi/100
            theta=k*np.pi/M
            tbl_cos_t[j,k]=np.cos(j*theta)
        #--end-k--
    #--end-j--
    return(tbl_cos_t)
#--

def fA_t(A,tbl_cos):
    N=len(A)
    TR,TC=tbl_cos.shape
    # fz(t)=Nz(0)+ 2*SUM (Nz(j)*cos(j*t))
    Na=auto_corr(A)
    #
    k=1 #index-theta
    vfA_t=np.zeros([TC,1])
    while (k<TC):
        fa_t=np.copy(Na[0])
        for m in range (1,len(Na)):
            fa_t += 2*Na[m]*tbl_cos[m,k] ## tblstart from1, not 0
        vfA_t[k,0]=fa_t
        #-end-for--
        k+=1
    #--
    return(vfA_t)
# #--
#
def is_zw_qualified(z, zeta_ohmega, tbl_cos):
    # evaluate z and w in step 2, and 3
    N=len(z)
    s_z=int(sum(sum(z)))
    f_zt=fA_t(z,tbl_cos)
    TR,TC=tbl_cos.shape
    TF_fzt=True
    for m in range(1,TC):
        TF_fzt=TF_fzt and f_zt[m]<(3*N-1)
    #
    # print('f_zt=', f_zt)
    TF_ZO=s_z in zeta_ohmega
    return(TF_fzt, TF_ZO, f_zt)
#--
if __name__=='__main__':
    #--
    X,Y,Z,W=generate_XYZW_6()
    N=len(X)
    #--- convert to array    
    ZZ=np.zeros([N,1])
    ZZ[:,0]=Z[0:N]
    #
    WW=np.zeros([N-1,1])
    WW[:,0]=W[0:N-1]
    #
    XX=np.zeros([N,1])
    XX[:,0]=X[0:N]
    YY=np.zeros([N,1])
    YY[:,0]=Y[0:N]
    #------
    H=construct_hadamard(X,Y,Z,W)
    if is_hadamard(H):
        print('H-matrix found')
        DG=np.matmul(H, H.transpose())
        # display indicator
        plt.imshow(DG.astype(float), cmap="gray")
        plt.show()
    #---
    zeta_ohmega={0,-1,1,-2,2,-3,3, -4,4}
    # zeta_ohmega={0,-1,1,-2,2,-3,3,-4,4}
    # createcosinus table
    tbl_cos=gen_tbl_cos_t(N)

    #    return(TF_fzt, TF_ZO, f_zt)
    TF_fzz,TF_ZOz,vQZ=is_zw_qualified(ZZ, zeta_ohmega, tbl_cos)
    TF_fzw,TF_ZOw,vQW=is_zw_qualified(WW, zeta_ohmega, tbl_cos)
    print('3N-1 -> ', 3*N-1)
    print('vQZ=',max(max(vQZ)),'TF_fzz:', TF_fzz,'TF_ZOz:', TF_ZOz)
    print('vQW=',max(max(vQW)),'TF_fzw:', TF_fzw,'TF_ZOw:', TF_ZOw)
    print('mean(vQW+vQZ)=',np.mean(np.mean(vQW+vQZ)))
    print('max(vQW+vQZ)=',np.max(np.max(vQW+vQZ)))
       
    plt.hist(vQZ+vQW)

    
