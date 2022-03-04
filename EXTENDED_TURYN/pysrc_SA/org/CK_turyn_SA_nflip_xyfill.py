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
if __name__=='__main__':
    N=12 #4,6, ...even number
    P=N-1
    M=4*(3*N-1) # order of H-matrix
    # generate XYZW sequence, delete some sy entries, try to recover by SA
    #x,y,z,w=generate_XYZW_10()
    Nr=4*N-11 # length of random vector <<-- (-11)
    #
    idx=0
    MAX_ITER=2*1000*1000 #1 million, for N=12->E=4.0
    # generate schedule
    # T_schedule=pmc_schedule(MAX_ITER,0.5)
    # param: max-iteration, starting value,decay speed
    T_schedule=pmc_schedule_turyn(MAX_ITER,0.5,2*10)
    v_flip=flip_schedule(MAX_ITER, Nr,1.0)
    v_E=np.zeros([MAX_ITER,1])
    t=np.zeros([MAX_ITER,1])
    #--
    q_curr=gen_seed(Nr)
#    [X,Y,Z,W]=construct_XYZW(q_curr)
    [X,Y,Z,W]=generate_XYZW_normal_SA(q_curr,N)
    # E_curr=tt_energy(X,Y,Z,W)
    E_curr=tt_energy_rev(X,Y,Z,W) 
#    E_curr=tt_energy_xyzw_sqr(X,Y,Z,W) 
    #
    epsilon=1e-6
    idx=0
    # while( (tt_check(X,Y,Z,W)==False) and (idx<MAX_ITER)) :
    while(not(E_curr<epsilon) and (idx<MAX_ITER)):
        #try-transition to a lower energy:: random-flip q
        q_next=np.copy(q_curr)
        # N_FLIP= 1+int(0.25*Nr)-int(0.25*Nr*(idx/(MAX_ITER)));
        N_FLIP=int(v_flip[idx])
        q_next=random_n_flip(q_next,N_FLIP)
#        [X,Y,Z,W]=construct_XYZW(q_next)
        [X,Y,Z,W]=generate_XYZW_normal_SA(q_next,N)

        # E_next=tt_energy(X,Y,Z,W)
        E_next=tt_energy_rev(X,Y,Z,W)
#        E_next=tt_energy_xyzw_sqr(X,Y,Z,W) 
        #
        if (E_next<E_curr) or (np.random.rand()>T_schedule[idx]):
           #accept transition
            q_curr=np.copy(q_next)
            E_curr=E_next
        #
        print('H-',M,'; E=', E_curr,';n_flip=',N_FLIP,'; iteration:',idx,'of',MAX_ITER)
        t[idx]=idx
        v_E[idx]=E_curr
        idx=idx+1  
    # -- end
    # Save XYZW vectors
    if(E_curr<epsilon):
        print('Hadamard matrix is found !')
        path_name='..//XYZW_VECTOR//'
        write_XYZW(X,Y,Z,W,path_name)
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
    plt.plot(v_E[0:idx])
    plt.show()
    # #
    # plt.plot(v_flip[0:idx])
    # plt.show()
    # #
    # #
    # plt.plot(T_schedule)
    # plt.show()
    

