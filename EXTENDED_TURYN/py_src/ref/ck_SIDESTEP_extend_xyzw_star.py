"""
STEP 2 – Find the possible strings for Z
We find all strings of length which satisfy the following.
1) Sum of the elements of is +/- z
2) fz(t) <= 3n-1; t in (k*pi/100);1<=k<=100
3) z0=1; z_{n-1}=1
4) Z>=Zr
"""
import pandas as pd
from prb_turyn import *
from prb_generate_XYZW import *
# ##
# ##
if __name__=='__main__':
    # --read data of XYZW* into array--
    path_name='DATA//'
    fname='XYZW_star_iter.txt'
    mydat=np.loadtxt(path_name+fname)
    #
    NR,NC=mydat.shape
    #
    N=int((NR+1)/4)
    #-- from X*,Y*,Z*,W*; insert spin-vector and check
    NXTD=4 #extend 4 bits each-->total 4*NXTD
    NQ=4*NXTD
    # 
    NFIN=NR+4*NXTD
    start_lag=int(NFIN/2)
    MAX_ITER=2**NXTD
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
            #            
            if (tt_check_start_lag(xx,yy,zz,ww,start_lag) \
            and is_XYZW_gt_rev(xx,yy,zz,ww)):                
            #if (is_XYZW_gt_rev(xx,yy,zz,ww)):   # xyzw_star oriented
            # if (tt_check(xx,yy,zz,ww)):       # hadamard-orientd
                tq=join_XYZW(xx,yy,zz,ww)
                xtd_data[:,cnt]=tq[0:NFIN,0]
                cnt=cnt+1
                # -------------
                # X=np.copy(xx)
                # Y=np.copy(yy)
                # Z=np.copy(zz)
                # W=np.copy(ww)
                # -------------
            #---
            idx=idx+1
        #-----
        print('Processing ', idx,'of', NC*MAX_ITER)
    #--------
    # #
    print('original',NC*MAX_ITER,'->',cnt)
    ## save data
    path_name='DATA//'
    fname='XYZW_extend_iter.txt'
    np.savetxt(path_name+fname,xtd_data[:,0:cnt])

