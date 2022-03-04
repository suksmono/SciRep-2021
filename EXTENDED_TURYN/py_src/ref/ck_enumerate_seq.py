#ck Baumert Hall
import matplotlib.pyplot as plt
from prb_small_williamson import *
# from prb_baumert_hall import *
from prb_turyn import *
import time
#
def int_to_spin(q,N):
    # fill bitswith 0's or spin -1
    x=[-1 for x in range(N)]
    t=[int(x) for x in bin(q)[2:]]
    #
    NB=len(t)
    for n in range(NB):
        t[n]=int(2*(t[n]-0.5))
    #    
    x[N-len(t):N]=t[0:len(t)]
    return(x)
#
#
def generate_XYZW_4():
    #W=np.zeros([4,1])
    #
    x=[1,-1,1,-1]
    y=[1,-1,-1,-1]
    z=[1,-1,-1,1]
    w=[1,1,1]
    return(x,y,z,w)
#
def generate_XYZW_6():
    #W=np.zeros([4,1])
    #
    x=[1,-1,1,1,-1,1]
    y=[1,1,1,-1,-1,1]
    z=[1,-1,1,1,1,-1]
    w=[1,-1,-1,-1,-1]
    return(x,y,z,w)
		
def generate_XYZW_8():
    #W=np.zeros([4,1])
    #
    x=[1,1,-1,1,-1,1,-1,1]
    y=[1,-1,-1,-1,-1,-1,-1,1]
    z=[1,-1,-1,1,1,1,1,-1]
    w=[1,1,1,-1,1,1,-1]
    return(x,y,z,w)
#
def auto_corr_start_lag(x,K):
    # non-periodic autocorrelation
    # starting at lag K
    N=len(x)
    ML=N#-1
    idx=0
    r=np.zeros([ML,1])
    for s in range(K,ML):
        t=0
        for m in range (N-s): 
            t=t+x[m]*x[m+s]
        #
        r[idx,0]=t
        idx=idx+1
    return(r)
    #--    
def tt_check_start_lag(X,Y,Z,W,K):
# check Nx+Ny+2Nz+2Nw=0 on TT-seq
    N=len(X)
    W1=np.zeros([N,1])
    for m in range(len(W)):
        W1[m]=W[m]
    Nt=auto_corr_start_lag(X,K) + auto_corr_start_lag(Y,K) + \
        2*auto_corr_start_lag(Z,K)+2*auto_corr_start_lag(W1,K)
    # ttF= (sum(vF)==(N-1));
    SS=np.sum(abs(Nt[1:len(Nt)]));
    return(SS<1e-6)
    # return(np.abs(SS)<1e-6)
# 6,587,373 count/s ~2^22 count/s
if __name__=='__main__':
    # N2=2
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
    N=4
    NQ0=4*N-1
    MAX_ITR=2**NQ0
    q_data=np.zeros([NQ0,MAX_ITR])
    idx=0
    itr=0
    # #
    while(itr<MAX_ITR):
        q=int_to_spin(itr,NQ0)
        start_lag=2
        X=q[0:N]
        Y=q[N:2*N]
        Z=q[2*N:3*N]
        W=q[3*N:4*N-1]
        if tt_check_start_lag(X,Y,Z,W,start_lag):
            q_data[:,idx]=q
            idx=idx+1
        #
        print('iteration->',itr,'of', MAX_ITR,'q=',q)
        itr=itr+1
    #              
    ## save data
    path_name='DATA//'
    fname='XYZW_star.txt'
    write_matrix(path_name+fname, q_data[:,0:idx-1])
    
    
