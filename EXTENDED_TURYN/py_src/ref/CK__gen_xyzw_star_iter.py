"""
----------------------------------------------------------------------------
Iteratively construct XYZW*; eg. X=#[x0,x1,*,*,*,*,x6,x7]
1. Run CK__gen_xyzw_iter.py for initial length |-> file: XYZW_star_iter_N.txt
2. Iteratively run CK__extend_xyzw_star_iter.py
3. Cek, there should be at least TT-seq, by calling 
   CK__extract_tt_sequence.py
----------------------------------------------------------------------------
"""
#
import matplotlib.pyplot as plt
from prb_small_williamson import *
from prb_turyn import *
import numpy as np
import time
from prb_generate_XYZW import *
##
if __name__=='__main__':
    # N2=2
    N=6  #[x0,x1,*,..*,x_{n-2},x_{n-1}],etc
    #
    NQ0= 4*N-1-10
    MAX_ITR=2**NQ0
    s_data=np.zeros([4*N,MAX_ITR])
    idx=0
    itr=0
    # #
    while(itr<MAX_ITR):
        q=int_to_spin(itr,NQ0)
        start_lag=int(N/2) #assume symmetric division
        #
        X,Y,Z,W=generate_XYZW_normal(q,N)
        x=np.asarray(X)
        y=np.asarray(Y)
        z=np.asarray(Z)
        w=np.asarray(W)
        if (tt_check_start_lag(X,Y,Z,W,start_lag) \
             and is_XYZW_gt_rev(x,y,z,w)):
            s=join_XYZW(X,Y,Z,W)
            s_data[:,idx]=s[:,0]
            idx=idx+1
        #
        print('iteration->',itr+1,'of', MAX_ITR) #,'q=',q)
        itr=itr+1
    #          
    print('Total', itr,'->',idx)
    ## save data
    path_name='DATA//'
    fname='XYZW_star_iter_'+str(N)+'.txt'
    np.savetxt(path_name+fname,s_data[:,0:idx])