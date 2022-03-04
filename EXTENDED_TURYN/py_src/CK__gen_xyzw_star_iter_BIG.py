"""
----------------------------------------------------------------------------
BIG:file-mode, hexa
Iteratively construct XYZW*; eg. X=#[x0,x1,*,*,*,*,x6,x7]
1. Run CK__gen_xyzw_iter.py for initial length |-> file: XYZW_star_iter_N.txt
2. Iteratively run CK__extend_xyzw_star_iter.py
3. Cek, there should be at least TT-seq, by calling 
   CK__extract_tt_sequence.py
NB:
1. all saved-W to file, the last0-padded removed already
----------------------------------------------------------------------------
"""
#
import matplotlib.pyplot as plt
from prb_small_williamson import *
from prb_turyn import *
import numpy as np
import time
from prb_generate_XYZW import *
# ##
#def is_y1yN_2_minus(Y):
#    N=len(Y)
#    return((Y[1,0]*Y[N-2,0])<0)

if __name__=='__main__':
    N=4 #[x0,x1,*,..*,x_{n-2},x_{n-1}],etc
    #
    NQ0= 4*N-1-10
    MAX_ITR=2**NQ0
    N_BIT=4*N
    ## save data
    path_name='DATA//'
    fname=path_name+'XYZW_star_iter_'+str(N)+'.txt'
    #--
    idx=0
    itr=0
    # #
    out_file = open(fname, "w")
    while(itr<MAX_ITR):
        q=int_to_spin(itr,NQ0)
        start_lag=int(N/2) #assume symmetric division
        #
        X,Y,Z,W=generate_XYZW_normal(q,N)
        if (tt_check_start_lag(X,Y,Z,W,start_lag) \
             and is_XYZW_gt_rev(X,Y,Z,W) and \
                 is_xx_yy_reflex_zero(X,Y)):
            #-- save as hex-text file --
            tv=join_XYZW(X,Y,Z,W)    
            #print(np.transpose(X),np.transpose(Y),\
            #      np.transpose(Z),np.transpose(W))
            tv_hex=vec_xyzw0_to_hex(tv, N_BIT)
            # --
#            print(tv_hex)
            out_file.write(str(tv_hex)+'\n')
            idx=idx+1
        #
#        print('iteration->',itr+1,'of', MAX_ITR) #,'q=',q)
        itr=itr+1
    #
    out_file.close()          
    print('Total', itr,'->qualified =',idx)
# ------- END OF CODE ------