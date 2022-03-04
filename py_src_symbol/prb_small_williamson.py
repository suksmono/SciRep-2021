#--------------------------------------------------------------
# Probabilistic Construction of Williamson Hadamard Matrix
# ** Simulated Annealing **
# When the energy formula is correct, SA should perform better
#--------------------------------------------------------------
# Williamson H-matrix
# H = [ A  B  C  D]
#     [-B  A -D  C]
#     [-C  D  A -B]
#     [-D -C  B  A]
# A=> {q0, q1}---> [ [q0 q1 q1]; [q1 q0 q1]; [q1 q1 q0] ]
# B=> {q2, q3}---> [ [q2 q3 q3]; [q3 q2 q3]; [q3 q3 q2] ]
# C=> {q4, q5}---> [ [q4 q5 q5]; [q5 q4 q5]; [q5 q5 q4] ]
# D=> {q6, q7}---> [ [q6 q7 q7]; [q7 q6 q7]; [q7 q7 q6] ]
# ---- E.g. ---
#  q=[1 -1 1 -1 -1 -1 1 -1];
#--------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
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
def write_matrix(fname, matx):
   # df=pd.DataFrame(data=matx.astype(int))
    df=pd.DataFrame(data=matx)
    df.to_csv(fname,sep=',', header=False, index=False)
#
def is_hadamard(H):
    #check wether H is Hadamard/orthogonal matrix
    M,N=H.shape
    DD=np.matmul(np.transpose(H),H)
    TF=((np.sum(np.sum(abs(DD)))-M**2)==0)
    return(TF)
#
def gen_seed(M):
    x=np.random.randint(2, size=M)
    x=x.astype(float)
    x=2*(x-0.5)
    # print('Test')
    return(x.astype(int))
#
def write_ABCD(A,B,C,D,path_name):
    M,N=A.shape
    M_ABCD=[A[0,:], B[0,:], C[0,:], D[0,:]]
    ABCD_name=path_name+'ABCD_'+str(4*N)+'.txt'
    write_matrix(ABCD_name, M_ABCD)
    #
#
def construct_symm_williamson(A,B,C,D):
    """
    input  : submatrix A,B,C,D
    output : williamson (hadamard) matrix
    --- williamson ---
    H = [  A   B   C   D ]
        [ -B   A  -D   C ]
        [ -C   D   A  -B ]
        [ -D  -C   B   A ]
    """
    M,N=A.shape
    W=np.empty([4*M,4*N], dtype=object)
    '''
    # row-1: [  A   B   C   D ]
    '''
    W[0:M,0:M]=     A[0:M,0:M]
    W[0:M,M:2*M]=   B[0:M,0:M]
    W[0:M,2*M:3*M]= C[0:M,0:M]
    W[0:M,3*M:4*M]= D[0:M,0:M]

    # row-2: [ -B   A  -D   C ]
    W[M:2*M,0:M] =      -B[0:M,0:M]
    W[M:2*M,M:2*M] =    A[0:M,0:M]
    W[M:2*M,2*M:3*M] = -D[0:M,0:M]
    W[M:2*M,3*M:4*M] =  C[0:M,0:M]

    '''
    # row-3: [ -D  -C   B   A ]
    '''
    W[3*M:4*M,0:M] =    -D[0:M,0:M]
    W[3*M:4*M,M:2*M] =  -C[0:M,0:M]
    W[3*M:4*M,2*M:3*M]= B[0:M,0:M]
    W[3*M:4*M,3*M:4*M]= A[0:M,0:M]
   
    '''
    # row-3: [ -C   D   A  -B ]
    '''
    W[2*M:3*M,0:M] =    -C[0:M,0:M]
    W[2*M:3*M,M:2*M] =  D[0:M,0:M]
    W[2*M:3*M,2*M:3*M]= A[0:M,0:M]
    W[2*M:3*M,3*M:4*M]= -B[0:M,0:M]
    #
    return(W)
#
def construct_baumert_hall(A,B,C,D):
    """
    input  : submatrix A,B,C,D
    output : williamson (hadamard) matrix
    --- williamson ---
    H = [  A   B   C   D ]
        [ -B   A  -D   C ]
        [ -C   D   A  -B ]
        [ -D  -C   B   A ]
    --- Baumert-Hall ---
    H=  [ A  A  A  B  -B  C -C -D   B  C -D -D]
        [ A -A  B -A  -B -D  D -C  -B -D -C -C]
        [ A -B -A  A  -D  D -B  B  -C -D  C -C]
        [ B  A -A -A   D  D  D  C   C -B -B -C]  
        
        [ B -D  D  D   A  A  A  C  -C  B -C  B] 
        [ B -C -D  D   A -A  C -A  -D  C  B -B]
        [ D -D  B -B   A -C -A  A   B  C  D -D] 
        [-C -C -C -D   C  A -A -A  -D  B -B -B]  
        
        [ D -D -B -B  -B  C  C -D   A  A  A  D]  
        [-D -C  C  C   C  B  B -D   A -A  D -A] 
        [ C -B -C  C   D -B -D -B   A -D -A  A] 
        [-C -D -D -C  -C -B  B  B   D  A -A -A]   
    """
    M,N=A.shape
    W=np.empty([4*M,4*N], dtype=object)

    '''
    # row-1: [ A  A  A  B  -B  C -C -D   B  C -D -D]
    '''
    # 1.block-1: A  A  A  B
    W[0:M,0:M]=     A[0:M,0:M]
    W[0:M,M:2*M]=   A[0:M,0:M]
    W[0:M,2*M:3*M]= A[0:M,0:M]
    W[0:M,3*M:4*M]= B[0:M,0:M]

    # 1.block-2: -B  C -C -D
    W[0:M,4*M:5*M]= -B[0:M,0:M]
    W[0:M,5*M:6*M]=  C[0:M,0:M]
    W[0:M,6*M:7*M]= -C[0:M,0:M]
    W[0:M,7*M:8*M]= -D[0:M,0:M]

    # 1.block-3: B  C -D -D
    W[0:M,8*M:9*M]=    B[0:M,0:M]
    W[0:M,9*M:10*M]=   C[0:M,0:M]
    W[0:M,10*M:11*M]= -D[0:M,0:M]
    W[0:M,11*M:12*M]= -D[0:M,0:M]

    '''
    # row-2: [ A -A  B -A  -B -D  D -C  -B -D -C -C]
    '''
  


    # #
    return(W)
#
#
def construct_ABCD(q):
#
    """
#    A=> {q0, q1}---> [ [q0 q1 q1]; [q1 q0 q1]; [q1 q1 q0] ]
#    B=> {q2, q3}---> [ [q2 q3 q3]; [q3 q2 q3]; [q3 q3 q2] ]
#    C=> {q4, q5}---> [ [q4 q4 q5]; [q5 q4 q5]; [q5 q5 q4] ]
#    D=> {q6, q7}---> [ [q6 q7 q7]; [q7 q6 q7]; [q7 q7 q6] ]
    #
    ## USE: np.roll(vector,2)
    """
    N=len(q)
    K05=int(N/4)
#    print('K05=', K05)

    #
    tA = [ q[0:K05], q[(K05-1):0:-1] ]   
    tB = [ q[K05:2*K05], q[(2*K05-1):K05:-1] ]   
    tC = [ q[2*K05:3*K05], q[(3*K05-1):2*K05:-1] ]   
    tD = [ q[3*K05:4*K05], q[(4*K05-1):3*K05:-1] ]   
    qA=[]
    for item in tA:
        qA.extend(item)
    qB=[]
    for item in tB:
        qB.extend(item)
    qC=[]
    for item in tC:
        qC.extend(item)
    qD=[]
    for item in tD:
        qD.extend(item)  
    #init
    NN=2*K05-1
    A=np.empty([NN,NN], dtype=object)
    B=np.empty([NN,NN], dtype=object)
    C=np.empty([NN,NN], dtype=object)
    D=np.empty([NN,NN], dtype=object)
    #
    A[:,0]=qA
    B[:,0]=qB
    C[:,0]=qC
    D[:,0]=qD
    for m in range(1,NN):
        A[:,m] = np.roll(qA,m)
        B[:,m] = np.roll(qB,m)
        C[:,m] = np.roll(qC,m)
        D[:,m] = np.roll(qD,m)
    #
    return(A,B,C,D)
#
#
def H_energy(H):
    M,N=H.shape
    #check wether H is Hadamard/orthogonal matrix
    DD=np.matmul(np.transpose(H),H)
    E= np.sum(np.sum(abs(DD)))-M*N
    return(E)
#
def random_flip(q):
    rd=np.random.permutation(len(q))
    #
    q[rd[1]]=-q[rd[1]];
    return(q)
######
def pmc_schedule__org (N,pStart):
    '''
    # Create probabilistic schedule
    # 1/4 part: pStart-0.99
    # 1/4 part: 0.999
    # 1/4 part: 0.9999
    # 1/4 part: 1.0
    '''
    N025=int(N/4)
    tSched=np.ones([N,1])
    tSched[N025:2*N025,0]=0.999
    tSched[2*N025:3*N025,0]=0.9999
    tSched[3*N025:N,0]=1.0
    #
    tt1= np.arange(0,N025+1,1)
    dGap=1.0-pStart
    xt1=tt1/N025
    st1=np.exp(-xt1)
    st1=pStart+dGap*(1-st1)/np.max(1-st1)
    tSched[0:N025,0] = st1[0:N025]
    return(tSched)
#
def pmc_schedule (N,pStart):
    '''
    # Create probabilistic schedule
    '''
    N05=int(N/2)
    tSched=np.ones([N,1])
    tSched[N05-2:int(1.5*N05),0]=0.999
    tSched[int(1.5*N05):int(1.9*N05),0]=0.9999
    #tSched[3*N025:N,0]=1.0
    #
    tt1= np.arange(0,N05,1)
    dGap=1.0-pStart
    xt1=tt1/N05
    st1=np.exp(-xt1)
    st1=pStart+dGap*(1-st1)/np.max(1-st1)
    tSched[0:N05,0] = st1[0:N05]
    return(tSched)
#