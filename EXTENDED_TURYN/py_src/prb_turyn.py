#
from prb_small_williamson import *
#
def construct_XYZW(v):
    Nr=len(v)
    N=int((Nr+5)/4)
    N_1=N-1
    X= np.zeros([N,1])
    Y= np.zeros([N,1])
    Z= np.zeros([N,1])
    W= np.zeros([N-1,1])
    ## insert 1 in first XYZW
    X[0]=1
    Y[0]=1
    Z[0]=1
    W[0]=1
    for m in range(N_1):
        X[m+1,0] = v[m]
        Y[m+1,0] = v[N_1+m]
        Z[m+1,0] = v[2*N_1+m]
    for m in range(N_1-1):    
        W[m+1,0] = v[3*N_1+m]
    #
    return(X,Y,Z,W)
#
def auto_corr(x):
    #non-periodic autocorrelation
    N=len(x)
    ML=N#-1
    idx=0
    r=np.zeros([ML,1])
    for s in range(ML):
        t=0
        for m in range (N-s): 
            t=t+x[m]*x[m+s]
        #
        r[idx,0]=t
        idx=idx+1
    return(r)
    #--    

def tt_check(X,Y,Z,W):
# check Nx+Ny+2Nz+2Nw=0 on TT-seq
    N=len(X)
    W1=np.zeros([N,1])
    for m in range(len(W)):
        W1[m]=W[m]
    Nt=auto_corr(X)+auto_corr(Y)+2*auto_corr(Z)+2*auto_corr(W1)
    # ttF= (sum(vF)==(N-1));
    SS=np.sum(abs(Nt[1:len(Nt)]));
    return(SS<1e-6)
    # return(np.abs(SS)<1e-6)
#
def tt_energy(X,Y,Z,W):
# check Nx+Ny+2Nz+2Nw=0 on TT-seq
    N=len(X)
    W1=np.zeros([N,1])
    for m in range(len(W)):
        W1[m]=W[m]
    Nt=auto_corr(X)+auto_corr(Y)+2*auto_corr(Z)+2*auto_corr(W1)
    # ttF= (sum(vF)==(N-1));
    SS=np.sum(abs(Nt[1:len(Nt)]));
    return(SS)
#
def construct_hadamard(X,Y,Z,W):
    #  construct H-matrix from XYZW
    # 1. create base sequence 
    NX=len(X)
    NW=len(W)
    #
    A=np.zeros([NX+NW,1]);
    B=np.zeros([NX+NW,1]);
    C=np.zeros([NX,1]);
    D=np.zeros([NX,1]);
    
    for m in range(NX):
        A[m]=Z[m]
        B[m]=Z[m]
        C[m]=X[m]
        D[m]=Y[m]
    #
    for m in range(NW):
        A[NX+m]= W[m]
        B[NX+m]=-W[m]
    #-- end
    # 2. Create T-sequences
    T1=np.zeros([2*NX+NW,1])
    T2=np.zeros([2*NX+NW,1])
    T3=np.zeros([2*NX+NW,1])
    T4=np.zeros([2*NX+NW,1])
    #
    for m in range(NX+NW):
        T1[m]=0.5*(A[m]+B[m])
        T2[m]=0.5*(A[m]-B[m])
    #
    for m in range (NX):
        T3[NX+NW+m]=0.5*(C[m]+D[m])
        T4[NX+NW+m]=0.5*(C[m]-D[m])
    #
    # 3. Construct A1,A2,A3,A4
    A1=  T1 + T2 + T3 + T4
    A2= -T1 + T2 + T3 - T4
    A3= -T1 - T2 + T3 + T4
    A4= -T1 + T2 - T3 + T4
    #
    # # 4. Construct circulant matrices
    XA1=construct_circular_matrix(A1)
    XA2=construct_circular_matrix(A2)
    XA3=construct_circular_matrix(A3)
    XA4=construct_circular_matrix(A4)
    # #
    H=construct_hadamard_turyn(XA1,XA2,XA3,XA4)
    # #---
    return(H)
# ---
    
def construct_circular_matrix(y):        
    # Construct a circular matrix for a given v-vector
    N=len(y)
    X=np.zeros([N,N])
    X[:,0]=y[0:N,0]
    # MM,NN=y.shape
    # print('y',MM,NN)
    for m in range (N-1):
        t=np.roll(y,m+1)
        X[:, m+1]=t[0:N,0]
    #
    return(X)
#
def construct_hadamard_turyn(A1,A2,A3,A4):
    #given 4-circulant matrices, construct a hadamard matrix
    M,N=A1.shape
    # construct a back identity matrix
    R=np.zeros([M,M])
    for m in range(M): 
        R[m,M-m-1]=1
    #
    # %%
    # H=[  A1     A2*R    A3*R    A4*R;
    #     -A2*R   A1      A4'*R  -A3'*R;
    #     -A3*R  -A4'*R   A1      A2'*R;
    #     -A4*R   A3'*R  -A2'*R   A1];
    # %% -- end --
    H=np.zeros([4*M,4*M])
    #first row: A1, A2*R, A3*R, A4*R
    H[0:M,0:M]= A1
    H[0:M,M:2*M]= np.matmul(A2,R)
    H[0:M,2*M:3*M]= np.matmul(A3,R)
    H[0:M,3*M:4*M]= np.matmul(A4,R)
   
    # second: -A2*R, A1, A4'*R, -A3'*R;
    H[M:2*M,0:M]= -np.matmul(A2,R)
    H[M:2*M,M:2*M]= A1
    H[M:2*M,2*M:3*M]= np.matmul(np.transpose(A4),R)
    H[M:2*M,3*M:4*M]= -np.matmul(np.transpose(A3),R)
    
    # third:-A3*R  -A4'*R   A1  A2'*R; 
    H[2*M:3*M,0:M]= -np.matmul(A3,R)
    H[2*M:3*M,M:2*M]= -np.matmul(np.transpose(A4),R)
    H[2*M:3*M,2*M:3*M]= A1
    H[2*M:3*M,3*M:4*M]= np.matmul(np.transpose(A2),R)
    
    # fourth:-A4*R   A3'*R  -A2'*R   A1]; 
    H[3*M:4*M,0:M]= -np.matmul(A4,R)
    H[3*M:4*M,M:2*M]= np.matmul(np.transpose(A3),R)
    H[3*M:4*M,2*M:3*M]= -np.matmul(np.transpose(A2),R)
    H[3*M:4*M,3*M:4*M]= A1
    #
    return(H)
#   
def write_XYZW(X,Y,Z,W,path_name):
    N=X.size
    M_XYZW=[X[:,0], Y[:,0], X[:,0], W[:,0]]
    XYZW_name=path_name+'XYZW_'+str(4*(3*N-1))+'.txt'
    write_matrix(XYZW_name, M_XYZW)
    #
#
def random_n_flip(q,N):
    # single spin-flip -1-> +1,v.v.
    rd=np.random.permutation(len(q))
    #
    q[rd[0:N]]=-q[rd[0:N]];
    return(q)
##
#
def pmc_schedule_turyn(N,pStart, c_speed):
    '''
    # Create probabilistic schedule
    # N: number of iteration
    # pStart: starting probability
    # c_speed: asymptotic rate
    '''
    tSched=np.ones([N,1])
   #
    c=c_speed
    #c=10.0
    tt1= np.arange(0,N,1)
    dGap=1.0-pStart
    xt1=tt1/N
    st1=np.exp(-c*xt1)
    st1=pStart+dGap*(1-st1)/np.max(1-st1)
    tSched[0:N,0] = st1[0:N]
    return(tSched)
#
def flip_schedule(N, N_start, c_speed):
    '''
    # Create probabilistic schedule
    # N: number of iteration
    # pStart: starting probability
    # c_speed: asymptotic rate
    '''
    v_flip=np.ones([N,1])
    v_sch=1-pmc_schedule_turyn(N, 0, c_speed)
    for m in range(N):
        v_flip[m,0]=1+int(N_start*v_sch[m,0])
    #
    return(v_flip)
## --
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
# #
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
    # check Nx+Ny+2Nz+2Nw=0 on TT-seq, but start at lag-K
    N=len(X)
    W1=np.zeros([N,1])
    # W1[1:N-2,0] = W[1:N,0] ##??
    for m in range(len(W)):
        W1[m]=W[m]
    Nt=auto_corr_start_lag(X,K) + auto_corr_start_lag(Y,K) + \
        2*auto_corr_start_lag(Z,K)+2*auto_corr_start_lag(W1,K)
    # ttF= (sum(vF)==(N-1));
    SS=np.sum(abs(Nt[1:len(Nt)]));
    return(SS<1e-6)
# 
def generate_XYZW_normal(q,N):
    x=np.zeros([N,1])
    y=np.zeros([N,1])
    z=np.zeros([N,1])
    w=np.zeros([N,1])
    #
    # x1=y1=z1=w1= 1 --->  #4
    x[0]=1
    y[0]=1
    z[0]=1
    w[0]=1
    #
    # x_n=y_n=-1,z_n= 1 --># 3
    x[N-1]=-1
    y[N-1]=-1
    z[N-1]= 1
    # x2=x_{n-1}=1, y2*y_{n-1}=-1 --> #2+1=#4
    x[1]= 1
    x[N-2]= 1
    # -- X--
    NX=N-4
    stX=0
    edX=stX+NX
    x[2:2+NX,0]=q[stX:edX]
    # -- Y--
    NY=N-3
    stY=edX
    edY=stY+NY
    y[1:1+NY,0]=q[stY:edY]
    y[N-2]= -y[1]
    # print(NY, stY, edY)
    # -- Z--
    NZ=N-2
    stZ=edY
    edZ=stZ+NZ
    z[1:1+NZ,0]=q[stZ:edZ]
    # print(NZ, stZ, edZ)
    # -- W--
    NW=N-1-1
    stW=edZ
    edW=stW+NW
    w[1:1+NW,0]=q[stW:edW]
    # print(NW, stW, edW)
    #
    return(x,y,z,w)
##
def join_XYZW(x,y,z,w):
    N=len(x)
    v=np.zeros([4*N,1])
    for m in range(N):
        v[m,0]=x[m]
        v[N+m,0]=y[m]
        v[2*N+m,0]=z[m]
        v[3*N+m,0]=w[m]
    #
    return(v)
##
def spin_to_bin(v,K):
    #convert vector of{-1,+1} into vector of{0,1}
    N=len(v)
    b=np.zeros([K,1])
    if K<N:
        return(0)
    else:
        #
        for m in range(N):
            # b[K-1-m]=round((v[N-1-m]+0.5)/2)
            b[K-1-m]=round((v[N-1-m].item()+0.5)/2)
        #
    #
    return(b)
##
def spin_to_int(v,K):
    # convert spin vector to integer
    s=spin_to_bin(v,K)
    val=0
    for m in range(K):
        val=val+s[K-1-m]*2**(m)
    #
    return(int(val[0]))
#
def bin_to_int(s):
    # convert binary vector to integer
    K=len(s)
    val=0
    for m in range(K):
        val=val+s[K-1-m]*2**(m)
    #
    return(int(val))
#
def is_spin_gt_rev(v,K):
    # is a spin vector A is greater-or-equal than its reverse Ar?
    # A>=Ar ??
    s_v=spin_to_bin(v,K)
    s_vr=np.zeros([K,1])
    for m in range(K):
        s_vr[m]=s_v[K-m-1]
    #
    v_val=bin_to_int(s_v)
    vr_val=bin_to_int(s_vr)
    #
    return(v_val>=vr_val)
 
#
def is_spin_gt_pm_rev(v,K):
    # is a spin vector A is greater-or-equal than its reverse +/- Ar?
    # A>=Ar ??
    s_v=spin_to_bin(v,K)
    s_vr=np.zeros([K,1])
    v_neg=np.zeros([K,1])
    
    for m in range(len(v)):
        v_neg[m]=-1*v[m]
    #
    s_v_neg=spin_to_bin(v_neg,K)
    s_vr_neg=np.zeros([K,1])
    for m in range(K):
        s_vr[m]=s_v[K-m-1]
        s_vr_neg[m]=s_v_neg[K-m-1]
    #
    v_val=bin_to_int(s_v)
    vr_val=bin_to_int(s_vr)
    vr_val_neg=bin_to_int(s_vr_neg)
    #
    return( (v_val>=vr_val) and (v_val>=vr_val_neg))
 
##
def is_xx_yy_reflex_zero(x,y):
    #check, whether xi*x_{n-1-i}+yi*y_{n-1-i}=0
    N=len(x)
    N05=int(N/2)
    for m in range(1,N05):
        tval=x[m]*x[N-1-m]+y[m]*y[N-1-m]
    #
    return(abs(tval)<0.001)
    
##
def is_XYZW_gt_rev(x,y,z,w):
    #check wether XYZW >= +/ XYZW_reversed?is_spin_gt_pm_rev
    N=len(x)
#    TFx=is_spin_gt_rev(x,N)
#    TFy=is_spin_gt_rev(y,N)
#    TFz=is_spin_gt_rev(z,N)
#    TFw=is_spin_gt_rev(w[0:N-1],N-1)

    TFx=is_spin_gt_pm_rev(x,N)
    TFy=is_spin_gt_pm_rev(y,N)
    TFz=is_spin_gt_pm_rev(z,N)
    TFw=is_spin_gt_pm_rev(w[0:N-1],N-1)
    #
    return(TFx and TFy and TFz and TFw)
##
def ft(A,t):
    Na=auto_corr(A)
    t_val=Na[0]
    N=length(Na)
    for jj in range(1,N):
        t_val=t_val+Na(jj)*np.cos(jj*t)
    #
    return(t_val)
##
def join_XYZW(x,y,z,w):
    N=len(x)
    v=np.zeros([4*N,1])
    for m in range(N):
        v[m,0]=x[m]
        v[N+m,0]=y[m]
        v[2*N+m,0]=z[m]
        v[3*N+m,0]=w[m]
    #
    return(v)
###
#
def extract_XYZW(v):
    # EXTRACT x,y,z,w from 4*N length vector, w padded with 0
    N=round(len(v)/4)
    x=v[0:N]
    y=v[N:2*N]
    z=v[2*N:3*N]
    w=v[3*N:4*N]
    #
    return(x,y,z,w)
#
##
def extend_XYZW(x,y,z,w,v):
    # -----------------------------------------------------------
    # extend xyzw in the middle by len(v)/4 each
    # x=[1,2,3,4]; y=[5,6,7,8], z=[9,10,11,12]; w=[13,14,15,16]
    # q=[11,22,33,44,55,66,77,88]
    # exten_XYZW(x,y,z,w,q)-> x=[1,2,11,22,3,4]; ... etc
    # -----------------------------------------------------------
    Nv=len(v)
    Nxtd=int(Nv/4)
    Nx=len(x)
    Nx05=int(Nx/2)
    #
    Nfin=Nx+Nxtd
    xx=np.zeros([Nx+Nxtd,1])
    yy=np.zeros([Nx+Nxtd,1])
    zz=np.zeros([Nx+Nxtd,1])
    ww=np.zeros([Nx+Nxtd,1])
    # split-n-fill
    #x
    xx[0:Nx05,0]=x[0:Nx05]
    xx[Nx05+Nxtd:Nfin,0]=x[Nx05:Nx]
    xx[Nx05:Nx05+Nxtd,0]=v[0:Nxtd]
    #y
    yy[0:Nx05,0]=y[0:Nx05]
    yy[Nx05+Nxtd:Nfin,0]=y[Nx05:Nx]
    yy[Nx05:Nx05+Nxtd,0]=v[Nxtd:2*Nxtd]
    #z
    zz[0:Nx05,0]=z[0:Nx05]
    zz[Nx05+Nxtd:Nfin,0]=z[Nx05:Nx]
    zz[Nx05:Nx05+Nxtd,0]=v[2*Nxtd:3*Nxtd]
    #w
    ww[0:Nx05,0]=w[0:Nx05]
    ww[Nx05+Nxtd:Nfin,0]=w[Nx05:Nx]
    ww[Nx05:Nx05+Nxtd,0]=v[3*Nxtd:4*Nxtd]
    #
    return(xx,yy,zz,ww)
##
##
def extend_spin_vector(x,v):
    # -----------------------------------------------------------
    # extend x or y..zw in the middle by len(v)/4 
    # x=[1,2,3,4];
    # q=[11,22]
    # exten_spin_vector(x,q)-> x=[1,2,11,22,3,4]; ... etc
    # -----------------------------------------------------------
    Nv=len(v)
    Nxtd=int(Nv)
    Nx=len(x)
    Nx05=int(Nx/2)
    #
    Nfin=Nx+Nxtd
    xx=np.zeros([Nx+Nxtd,1])
    # split-n-fill
    #x
    xx[0:Nx05,0]=x[0:Nx05]
    xx[Nx05+Nxtd:Nfin,0]=x[Nx05:Nx]
    xx[Nx05:Nx05+Nxtd,0]=v[0:Nxtd]
    #
    return(xx)
##
##
def vec_to_hex(tv,N_BIT):
    tv_bin=spin_to_bin(tv,N_BIT)
    tv_int=bin_to_int(tv_bin)
    tv_hex=np.base_repr(tv_int,base=16)
    #
    return(tv_hex)
##
##
def vec_xyzw0_to_hex(tv, N_BIT):
    # prinsip: padding should be removed when it is saved
    # convert XYZW-vector; where W last position was-padded with 0,
    # to hex with removed 0 in last position of W
    tv_bin=spin_to_bin(tv[0:N_BIT-1],N_BIT)
    tv_int=bin_to_int(tv_bin)
    tv_hex=np.base_repr(tv_int,base=16)
    return(tv_hex)
##
def hex_to_vec_xyzw0(tv_hex, N_BIT):
    # prinsip: padding should be removed when it is saved
    # prinsip: padding can be inserted in reading
    # convert hex-to-XYZW-vector 
    # the W's last position is padded with 0
    v0=int(tv_hex, 16) # directly to integer
    v1=int_to_spin(v0, N_BIT)
    v2=np.copy(v1)
    v2[0:N_BIT-1]=v1[1:N_BIT]
    v2[N_BIT-1]=0
    return(v2)
##--
def is_y1yN_2_minus(Y):
    N=len(Y)
    return((Y[1,0]*Y[N-2,0])<0)