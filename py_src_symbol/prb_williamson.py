# prb_williamson
from sympy import *
from prb_symbolic import *
#import neal
#import numpy as np
###
import pandas as pd
#
def writeHMatrix(fname, matx):
   # df=pd.DataFrame(data=matx.astype(int))
    df=pd.DataFrame(data=matx)
    df.to_csv(fname,sep=',', header=False, index=False)
    
def readHMatrix(fname):
    data=pd.read_csv(fname,sep=',', header=None )
    return data.as_matrix() #.tolist()
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
def calculate_Hks(D):
    '''
    given D, calculate Hks
    '''
    M,N=D.shape
    Hks=0
    for m in range(0,M):
        for n in range(m+1,M):
            #print('m=',m,'n=',n)
            t=expand(D[m,n]**2)
            Hks=Hks+t
    #
    return(Hks)
#
def remove_squares(Hks,ss):
    '''
    # substitute: s**2->1, s**3->s, s**4->1
    '''
    for m in range(len(ss)):
        Hks= Hks.subs( {ss[m]**2:1})
        Hks= Hks.subs( {ss[m]**3:ss[m]})
        Hks= Hks.subs( {ss[m]**4:1}) 
    return(simplify(Hks))
#
#
def Hks_to_Hkq(Hks,ss,qq):
    '''
    # substitute: s**2->1, s**3->s, s**4->1
    '''
    Hkq=Hks
    for m in range(len(ss)):
        Hkq= Hkq.subs( {ss[m]:qq[m]})
    return(simplify(Hkq))

#    #
def gen_subs_pair(qq,M):
    NR=int(M*(M-1)/2)
    tbl_q=np.empty([NR,3], dtype=object)
    idx=0
    for m in range(M):
        for n in range(m+1,M):
            # print('idx:', idx+M, ':', qq[m],qq[n],'->',qq[idx+M])
            tbl_q[idx,0]=qq[m]
            tbl_q[idx,1]=qq[n]
            tbl_q[idx,2]=qq[idx+M]
            #
            idx=idx+1
        #--
    #--
    return(tbl_q)

#
def C_boole(qi,qj,qk,dij):
    '''
    Boolean reduction qi*qj <- qk + C(qi,qj,qk;dij)
    >> this one do: C=dij*(3qk + qi*qj - 2qi*qk-2qj*qk)
    '''
    return(dij*(3*qk + qi*qj -2*qi*qk - 2*qj*qk))
##
def Hkq_to_H2q(Hkq,spair,delta):
    H2q=Hkq
    for tsp in spair:
        #print(tsp[0],'*',tsp[1], '->' ,tsp[2])
        H2q=H2q.subs({ \
                   tsp[0]*tsp[1]:tsp[2]\
                          }) \
                + C_boole(tsp[0], tsp[1], tsp[2],delta)
    # 
    print('Simplifying H2q ...')
    H2q=simplify(H2q)
    #
    return(H2q)
#
def H2q_to_H2s(H2q,qq0,ss0):
    H2s=H2q
    N=len(qq0)
    for m in range(N):
        print(' q<-s substitution:',m, 'of',N)
        H2s = H2s.subs( {qq0[m]:q2s(ss0[m])} )
    #
    print('Simplifying H2s ...')
    H2s=simplify(H2s)
    #
    return(H2s)
##
def test_ABCD(q):
    N=len(q)
    K05=int(N/4)
    print('K05=', K05)
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
    #
    return(qA,qB,qC,qD)

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
def is_hadamard(W):
    M,N=W.shape
    DD=np.matmul(W,W.transpose())
    SS=sum(sum(abs(W)))
    if abs(SS-M*N)<1:
        return(True)
    else:
        return(False)
######