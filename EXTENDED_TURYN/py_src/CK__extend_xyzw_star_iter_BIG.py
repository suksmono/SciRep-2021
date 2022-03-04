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
##
if __name__=='__main__':
    # --read data of XYZW* into array--
    # original length
    N_ORG=6
    N_BIT_ORG=4*N_ORG
    # extension length
    NXTD=2 #extend number-of-bits each
    N_BIT_XTD=4*(N_ORG+NXTD)
    #--- file definition ---
    path_name='DATA//'
    fname_org=path_name+'XYZW_star_iter_'+str(N_ORG)+'.txt'
    fname_xtd=path_name+'XYZW_star_iter_'+str(N_ORG+NXTD)+'.txt'
    ## prepare open file
    in_file = open(fname_org, "r")
    out_file = open(fname_xtd, "w")
    ##
    num_lines = sum(1 for line in open(fname_org))
    #-- from X*,Y*,Z*,W*; insert spin-vector and check
    NQ=4*NXTD
    # 
    NFIN=4*(N_ORG+NXTD) #NR+4*NXTD
    start_lag=int(NFIN/2)
    MAX_ITER=2**NQ
    cnt=0
    idx=0
    for tv in in_file:
        # --
        v=hex_to_vec_xyzw0(tv, N_BIT_ORG)
        x,y,z,w=extract_XYZW(v)
        #
        for itr in range(MAX_ITER):
            q=int_to_spin(itr,NQ)
            xx,yy,zz,ww=extend_XYZW(x,y,z,w,q)
            start_lag=int(len(xx)/2)
            if (tt_check_start_lag(xx,yy,zz,ww,start_lag) \
                 and is_XYZW_gt_rev(xx,yy,zz,ww) \
                     and is_xx_yy_reflex_zero(xx,yy)):
                #--
                tv=join_XYZW(xx,yy,zz,ww)    
                tv_hex=vec_xyzw0_to_hex(tv, N_BIT_XTD)
                #
                out_file.write(str(tv_hex)+'\n')
                cnt=cnt+1
                # -------------
            #--end-if--
            idx=idx+1
            print('Processing ', idx,'of', num_lines*MAX_ITER)
        #---end-for-itr---
    #---end-for-m-----
    ##
    print('original',idx,'->',cnt)
    ## save data
    in_file.close()
    out_file.close()
    # path_name='DATA//'
    # np.savetxt(fname_xtd,xtd_data[:,0:cnt])

