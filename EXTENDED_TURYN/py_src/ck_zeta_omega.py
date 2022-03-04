# -*- coding: utf-8 -*-
"""
calculate zw, the root of z and w that satisfy
x^2 + y^2 + 2z^2 + 2w^2=6*N-2 
"""
#
import numpy as np
#from itertools import combinations 
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
if __name__=='__main__':
    #
    N=24
    QR=generate_zeta_ohmega(N)
    print('Root found of order N',N, 'are',QR)




