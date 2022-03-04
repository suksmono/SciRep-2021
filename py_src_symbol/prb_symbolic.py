"""
prb_symbolic
symbolic manipulation module
"""

import numpy as np
from sympy import *
import itertools as itr

"""
-------------------------------------------------------------------------------
Formulate/Symbolically-calculate Hks for a given H-order
order can be recognized from sList
-------------------------------------------------------------------------------
"""
def formHks(N_VECT,sList):
    s=sList
    M=N_VECT
    Hks=0;
    for m in range(0,M):
        ts=s[m]
        for nn in range(m+1,M):
            Hks = Hks + (np.dot(s[m],s[nn]))**2 
    Hks=expand(Hks)
    return(Hks)

"""
for the 4x4 case, the arrangement of of the qubits is
             s-DOMAIN
-------------------------------------------------------------------------------
 |<---problem-->|<------ ancillas ------>|
-------------------------------------------------------------------------------
  s0 s4   1  1  |s8 
  s1 s5  -1  1  |s9 
  s2 s6   1  1  |s10 
  s3 s7  -1  1  |s11 
-------------------------------------------------------------------------------
"""
# generate s_matrix with constraint
def gen_csqmatrix(M,NQ,sk):
    NC=len(sk)
    NO=M-NC
    ss=symbols('s0:%d'%NQ)
    ss0=np.asarray(ss)
    ss1=ss0.reshape(int(NQ/M),M).tolist()
    '''
    -------------------------------------------------------------------------------
    insert known vectors here
    -------------------------------------------------------------------------------
    '''
    s=[ss1[0]]
    # put variables of searched vect0r
    for m in range(1,NO):
        s.append(ss1[m])
    # insert known vevtor to the rest of H
    for m in range(0,NC):
        s.append(sk[m])
    
    # insert ancillary variables
    for m in range(NO,int(NQ/M)):
        #print(m)
        s.append(ss1[m])
    '''
    -------------------------------------------------------------------------------
    gen q
    -------------------------------------------------------------------------------
    '''
    qq=symbols('q0:%d'%NQ)

    qq0=np.asarray(qq)
    qq1=qq0.reshape(int(NQ/M),M).tolist()
    #s0=np.asarray(ss)
    #ss1=ss0.reshape(int(NQ/M),M).tolist()
    '''
    -------------------------------------------------------------------------------
    insert known vectors here
    -------------------------------------------------------------------------------
    '''
    qk=ones(1,M).tolist() #[ [1,1,1,1] ]
    oneM=qk[0]
    for m in range(1,NC):
        qk.append(oneM)
    
    q=[qq1[0]]
    # put variables of searched vect0r
    for m in range(1,NO):
        q.append(qq1[m])
    # insert known vevtor to the rest of H
    for m in range(0,NC):
        q.append(qk[m])
    
    # insert ancillary variables
    for m in range(NO,int(NQ/M)):
        #print(m)
        q.append(qq1[m])
    return(s,q)

'''
------------------------------------------------------------
 generate list (vector) of substitution pair
 i.e., for qi*qj->qk, the pair is [i, j, k]
------------------------------------------------------------
'''
def genSubsPair(M):    
    w, h = 3, int(M*(M-1)/2)
    spair = [[0 for x in range(w)] for y in range(h)]
    listCol=list(range(0, M))
    lc=list(itr.combinations(listCol,2))
    for m in range(0,len(lc)):
        spair[m]=[lc[m][0], lc[m][1], M+m]
    return(spair)
#
def genSubsPairNEW(M,NO): 
    w, h = 3, int(NO*(NO-1)/2)
    spair = [[0 for x in range(w)] for y in range(h)]
    listCol=list(range(0, NO))
    lc=list(itr.combinations(listCol,2))
    for m in range(0,len(lc)):
        spair[m]=[lc[m][0], lc[m][1], M+m]
    return(spair)

'''
------------------------------------------------------------
 define Hks->Hkq transform as function
------------------------------------------------------------
'''
def Hks2Hkq(Hks,s,q):
    [NCOL, NROW]=np.shape(s)
    Hkq=Hks
    k=0
    for m in range(0,NCOL):
        for n in range(0,NROW):
            #Hkq= expand(Hkq.subs( {s[m][n]:s2q(q[m][n])})) 
            Hkq= Hkq.subs( {s[m][n]:s2q(q[m][n])}) 
            k=k+1
            print('Processing si->qi: ', k, 'of', NROW*NCOL)
    Hkq=expand(Hkq)
    #Hkq=simplify(Hkq)
    return(Hkq)

#
def Hks2HkqNEW(Hks,ss,qq):
    Hkq=Hks
    for m in range(0,len(ss)):
        Hkq= Hkq.subs( {ss[m]:s2q(qq[m])}) 
        print('Processing si->qi: ', m, 'of', len(ss)-1)
    Hkq=expand(Hkq)
    return(Hkq)
'''
------------------------------------------------------------
REMOVE ALL SELF_SQUARED TERMS: si**2 -> 1
------------------------------------------------------------
'''
def rmvIdSquare(Hks,s):
    [NCOL, NROW]=np.shape(s)
    k=0
    for m in range(0,NCOL):
        for n in range(0,NROW):
            #Hks= simplify(Hks.subs( {s[m][n]**2:1})) 
            Hks= Hks.subs( {s[m][n]**2:1}) 
            k=k+1
            print('Processing qi^2->1: ', k, 'of', NROW*NCOL)
    return(simplify(Hks))
#
def rmvIdSquareNEW(Hks,ss):
    for m in range(0,len(ss)):
        Hks= Hks.subs( {ss[m]**2:1}) 
        print('Processing qi^2->1: ', m, 'of', len(ss)-1)
    return(simplify(Hks))

'''
------------------------------------------------------------
maximum value of Hq, used to estimate delta > max(Hq)
WHAT IS THE MIN for threshold? it should be in s-domain
------------------------------------------------------------
'''
def minMaxHq(Hq, NQ):
    q=symbols('q0:%d'%NQ)
    vE=zeros(2**NQ,1)
    for m in range(0,2**NQ):
        vq=np.binary_repr(m, width=NQ)
        tE=Hq
        for n in range(0,NQ):
            tE=tE.subs({q[n]:float(vq[n]) })
        # put energy spectra to vE
        vE[m]= float(tE)
    return(min(vE), max(vE))
'''
------------------------------------------------------------
k-body to 2-body transform
------------------------------------------------------------
'''
def H2sub(x1,x2,y,d12):
    return(d12*(3*y+x1*x2-2*x1*y-2*x2*y))

'''
There are two binary representation
    s={-1,+1} and q={0,1}
The default in Ising simulation is s-domain
------------------------------------------------------------
symbolic: q-to-s and s-t-q transform
------------------------------------------------------------
'''
def q2s(x):
    return(1/2 -x/2)

def s2q(x):
    return(1 - 2*x)

'''
------------------------------------------------------------
numerical: q-to-s and s-to-q transforms
------------------------------------------------------------
'''
 # define function vq2s
def vq2s(x):
    return(1-2*x)
# define function vs2q
def vs2q(x):
    return(1/2-x/2)

'''
------------------------------------------------------------
EXTRACT ISING PARAMETERS THE RESULTS FROM SYMBOLIC SOLUTION
copy coefficients into {b, hi, Jij} of Ising parameters
input: Hamiltonian in s-domain H(s0, s1, ...)
output: Ising parameters {b, hi, Jij}
------------------------------------------------------------
'''
def isingCoeffs(Hs, NQ):
    '''
    since we work in s-domain, assume the symbols are
    s0, s1, ...., s_(NQ-1)
    '''
    hi=np.zeros(NQ)
    Jij=np.zeros((NQ,NQ))
    dc=Hs.as_coefficients_dict()
    # list all symbols: s0, s1, ...
    ss=symbols('s0:%d'%NQ)

    # extract b
    b=dc[1]
    # extract hi
    for m in range(NQ):
        hi[m]=dc[ss[m]];
    # extract Jij
    for m in range(NQ):
        for n in range(m+1,NQ):
            #print('m=',m,'n=',n)
            Jij[m,n]=dc[ss[m]*ss[n]]
    # return the results
    return(b, hi, Jij)
#
'''
------------------------------------------------------------
    SYMBOLIC COMPUTATION
------------------------------------------------------------
  >input:  Hkq in q-domain
  >output: H2s in s-domain
 ALGORITHM
  1> Expand H
  2> Simplify-1: substitute  qi**2 <- qi
  3> Simplify-2: transform k-body terms-> 2-body+compenste
  4> Transform Hq-> Hs
Process:
    Hq0 -> Hq -> Hq2b -> H2s
 Hq0: original input
 Hq1: expanded Hq0
 Hq2b: Hq with only 2-body interaction
 H2s: H in s-domain with 2-body interaction
------------------------------------------------------------
'''

def q2s_symbolig(Hkq,q,s,spair,delta):
    #dc=Hkq.as_coefficients_dict()    
    [NCOL, NROW]=np.shape(s)
    if len(spair)>0:
        [NPAIR, XX]=np.shape(spair)
    else:
        NPAIR=0
    d=delta; #2*vMax 
    # do substitution iteratively 
    H2q=Hkq
    # ------------------------------------------------------------------------
    # row -0-1
    # ------------------------------------------------------------------------
    for k in range(0,NPAIR):
        for m in range(0, NROW):
            ii=spair[k][0]
            jj=spair[k][1]
            kk=spair[k][2]
            print(ii,'*' , jj, '->',kk )
            H2q=H2q.subs({ \
                  q[ii][m]*q[jj][m]:q[kk][m]\
                          }) \
                + H2sub(q[ii][m], q[jj][m],q[kk][m],d)
    
    print('Simplifying H2q ...')
    H2q=simplify(H2q)
    '''
    ------------------------------------------------------------
    TRANSFORM H in q-domain TO s-DOMAIN FOR SIMULATION
    ------------------------------------------------------------
    '''
    print('Transform H2q->H2s ...')
    H2s=H2q;
    for m in range(0,NCOL):
        for n in range(0,NROW):
            print('substitute col-', m,'/',NCOL, 'row-',n,'/',NROW)
            H2s = H2s.subs( {q[m][n]:q2s(s[m][n])} )
    #
    #print('Simplify H2s ...')
    #H2s=simplify(H2s)
    H2s=expand(H2s)
    return(H2q, H2s)   