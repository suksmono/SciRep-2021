#  ------------------------------------------------------------------------
#  Work-flow
#  ------------------------------------------------------------------------
#  1. Find Turyn-Type sequences: TT_algo-> {X,Y,Z,W}
#  2. Create base-sequences: {X,Y,Z,W} -> {A,B,C,D}
#  3. Create T-sequences: {A,B,C,D}->{T1,T2,T3,T4}
#  4. Create seed sequences: {T1,T2,T3,T4} -> {A1,A2,A3,A4}
#  5. Create circular matrices: {A1,A2,A3,A4}->{XA1,XA2, XA3, XA4}
#  6. Construct Hadamard matrix:
#  H=[  A1     A2*R    A3*R    A4*R;
#      -A2*R   A1      A4'*R  -A3'*R;
#      -A3*R  -A4'*R   A1      A2'*R;
#      -A4*R   A3'*R  -A2'*R   A1];
#  where R is back-diagonal identity matrix, eg.: 
#   R = [ 0 0 1; 
#         0 1 0;
#         1 0 0]
#  ------------------------------------------------------------------------
from prb_small_williamson import *
from prb_turyn import *
from prb_generate_XYZW import *
# #
def construct_XY_antisymm(q,xo,yo,N):
    NO=len(xo)
    NO05=int(NO/2)
    NX=N-2*NO05
    NX05=int(NX/2)
    #--
    x=np.zeros([N,1])
    y=np.zeros([N,1])
    # prefill
    x[0:NO05,0]=xo[0:NO05,0]
    x[NX+NO05:N,0]=xo[NO05:NO,0]
    #
    y[0:NO05,0]=yo[0:NO05,0]
    y[NX+NO05:N,0]=yo[NO05:NO,0]
    
    # fill with SA-vector-q    
    x[NO05:NX+NO05,0]=q[0:NX]
    ty=q[NX:NX+NX05]
    # --[first half of y]
    y[NO05:NX05+NO05,0]=ty[0:NX05]
    for m in range(NX05):
        tq=-1*x[NO05+m]*x[N-1-NO05-m]*y[NO05+m]
        y[N-1-NO05-m,0]=tq
    #--
    return(x,y)
##
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
##
if __name__=='__main__':
    N=36 #4,6, ...even number
    # generate XYZW sequence, delete some sy entries, try to recover by SA
    X,Y,Z,W=generate_XYZW_36()
    #
    N05=int(N/2)
    P=N-1
    M=4*(3*N-1) # order of H-matrix
    # assume original-> extended 
    NX=10
    NO=N-NX
    #
    NX05=int(NX/2)
    NO05=int(NO/2)
    xo=np.zeros([NO,1])
    yo=np.zeros([NO,1])
    # simulate original x and y, before extension
    # [NO05;NX05;NX05;NO05]
    xo[0:NO05,0]=X[0:NO05]
    xo[NO05:NO,0]=X[NO05+NX:N]
    #
    yo[0:NO05,0]=Y[0:NO05]
    yo[NO05:NO,0]=Y[NO05+NX:N]
    #
    # #-- lost string
    # x1=np.zeros([NX,1])
    # y1=np.zeros([NX,1])
    # # y1=np.zeros([NX05,1])
    # x1[0:NX,0]=X[NO05:NO05+NX]
    # y1[0:NX,0]=Y[NO05:NO05+NX]
    # use property xi*x_{n-1-i}+yi*y_{n-1-i}=0
    # or, since xi,yi={-,+}=> y_{n-1-i}= -(xi*x_{n-1-i})*yi 
    NQ=int(3*NX/2) # number of unknown elements
    
    idx=0
    q=int_to_spin(idx,NQ)
    X,Y=construct_XY_antisymm(q,xo,yo,N)
    E=tt_energy_rev(X,Y,Z,W) 
    epsilon=1e-6
    idx=0
    MAX_ITER=2**(NQ-1) # nb:1 million, for N=12->E=4.0
    while(not(E<epsilon) and (idx<MAX_ITER)):
        #try-transition to a lower energy:: random-flip q
        q=int_to_spin(idx,NQ)
        X,Y=construct_XY_antisymm(q,xo,yo,N)
        E=tt_energy_rev(X,Y,Z,W)
        #--
        # print('H-',M,'q=',q,'; E=', E,'; iteration:',idx,'of',MAX_ITER)
        print('H-',M,'; iteration:',idx,'of',MAX_ITER)
        idx+=1
    # -- end
    # print('Lost q=',np.transpose(x1),np.transpose(y1))
    # Save XYZW vectors
    if(E<epsilon):
        print('Hadamard matrix is found !')
        #path_name='..//XYZW_VECTOR//'
        #write_XYZW(X,Y,Z,W,path_name)
    else:
        print('Hadamard matrix is NOT found !')        
    #--
    H=construct_hadamard(X,Y,Z,W)
    DG=np.matmul(H, H.transpose())
    # display indicator
    plt.imshow(DG.astype(float), cmap="gray")
    plt.show()
    # display H-matrix
    plt.imshow(H.astype(float), cmap="gray")
    plt.show()
     #plot energy curve
