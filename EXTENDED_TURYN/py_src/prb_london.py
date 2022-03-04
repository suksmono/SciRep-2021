from prb_turyn import *
from itertools import permutations 

def is_xyzw_root(tp,N):
    RHS=6*N-2
    x=tp[0]
    y=tp[1]
    z=tp[2]
    w=tp[3]
    tval=x**2+y**2+2*z**2+2*w**2
    return(tval==RHS)
#
def generate_zeta_ohmega(N):
    # generat a set of root satisfying
    RHS=6*N-2
    # list all perfect square less than RHS
    R=set()
    for k in range(int(np.sqrt(RHS))+1):
        R.add(k)
    #--
    # make all permutations of (R,4): x,y,z,w 
    LR=list(R)
    P=list(permutations(LR,4))
    #
    QR=set()
    for tp in P:
#        print(tc)
        if is_xyzw_root(tp,N):
            print('root',tp)
            for k in range(2,4):  
                # extract zeta an ohmega only
                QR.add(tp[k])
                QR.add(-1*tp[k])
            #--
        #---
    #----
    return(QR)
#
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