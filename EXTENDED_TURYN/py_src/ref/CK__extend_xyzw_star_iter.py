"""
----------------------------------------------------------------------------
Iteratively construct XYZW*; eg. X=#[x0,x1,*,*,*,*,x6,x7]
1. Run CK__gen_xyzw_iter.py for initial length |-> file: XYZW_star_iter_N.txt
2. Iteratively run CK__extend_xyzw_star_iter.py
3. Cek, there should be at least TT-seq, by calling 
   CK__extract_tt_sequence.py
--
 Extending XYZW* 
     X = [x0, x1, *   *   *  *  *    *, x8, x9] >> R8, R9
    XX = [x0, x1, x2, x3, *  *, x6, x7, x8, x9] >> R6, R7, R8, R9
----------------------------------------------------------------------------
"""
##
import pandas as pd
from prb_turyn import *
from prb_generate_XYZW import *
###
if __name__=='__main__':
    # --read data of XYZW* into array--
    # original length
    N_ORG=8
    # extension length
    NXTD=2 #extend number-of-bits each
    path_name='DATA//'
    fname_org=path_name+'XYZW_star_iter_'+str(N_ORG)+'.txt'
    fname_xtd=path_name+'XYZW_star_iter_'+str(N_ORG+NXTD)+'.txt'
    mydat=np.loadtxt(fname_org)
    #
    NR,NC=mydat.shape
    #
    N=int((NR+1)/4)
    #-- from X*,Y*,Z*,W*; insert spin-vector and check
    NQ=4*NXTD
    # 
    NFIN=NR+4*NXTD
    start_lag=int(NFIN/2)
    MAX_ITER=2**NQ
    cnt=0
    idx=0
    xtd_data=np.zeros([NFIN, NC*MAX_ITER])
    for m in range(NC):#NC
        v=mydat[:,m]
        x,y,z,w=extract_XYZW(v)
        #
        for itr in range(MAX_ITER):
            q=int_to_spin(itr,NQ)
            xx,yy,zz,ww=extend_XYZW(x,y,z,w,q)
            start_lag=int(len(xx)/2)
            if (tt_check_start_lag(xx,yy,zz,ww,start_lag) \
                 and is_XYZW_gt_rev(xx,yy,zz,ww)):
                tq=join_XYZW(xx,yy,zz,ww)
                xtd_data[:,cnt]=tq[0:NFIN,0]
                cnt=cnt+1
                # -------------
                # X=np.copy(xx)
                # Y=np.copy(yy)
                # Z=np.copy(zz)
                # W=np.copy(ww)
                # -------------
            #--end-if--
            idx=idx+1
            print('Processing ', idx,'of', NC*MAX_ITER)
        #---end-for-itr---
    #---end-for-m-----
    ##
    print('original',NC*MAX_ITER,'->',cnt)
    ## save data
    path_name='DATA//'
    np.savetxt(fname_xtd,xtd_data[:,0:cnt])

