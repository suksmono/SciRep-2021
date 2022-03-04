# ck big-turyn
# 6,587,373 count/s ~2^22 count/s
# we will find a sequence of length 8
#[x0,x1,*,*,*,*,x6,x7]
#
# first generate all sequence[x0,x1,x6,x7]...
# such that (Nx+Ny+2Nz+2Nw)(s)=0;for s>=6 
# -------------------------------------------------
# -- to be considered later--
# normalized case
# x1=y1=z1=w1= 1 --->  #4
# x_n=y_n=-1,z_n= 1 --># 3
# x2=x_{n-1}=1, y2*y_{n-1}=-1 --> #2+1=#4
# total qubits= (4N-1) -11
# -------------------------------------------------

import matplotlib.pyplot as plt
from prb_small_williamson import *
# from prb_baumert_hall import *
from prb_turyn import *
import numpy as np
import time
from prb_generate_XYZW import *
# #
##
if __name__=='__main__':
    # N2=2
    N=8  #[x0,x1,*,..*,x_{n-2},x_{n-1}],etc
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
        #if (tt_check_start_lag(X,Y,Z,W,start_lag)):
        if (tt_check_start_lag(X,Y,Z,W,start_lag) \
            and is_XYZW_gt_rev(x,y,z,w)):
            s=join_XYZW(X,Y,Z,W)
            s_data[:,idx]=s[:,0]
            idx=idx+1
        #
        print('iteration->',itr,'of', MAX_ITR) #,'q=',q)
        itr=itr+1
    #          
    print('Total', itr,'->',idx)
    ## save data
    path_name='DATA//'
    fname='XYZW_star.txt'
    #write_matrix(path_name+fname, s_data[:,0:idx-1])
    np.savetxt(path_name+fname,s_data[:,0:idx-1])