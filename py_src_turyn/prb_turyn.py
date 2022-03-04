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