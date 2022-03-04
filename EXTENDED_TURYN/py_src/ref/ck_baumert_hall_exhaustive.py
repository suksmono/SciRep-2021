#ck Baumert Hall
import matplotlib.pyplot as plt
from prb_small_williamson import *
from prb_baumert_hall import *

#
if __name__=='__main__':
    K=3
    K05=int(K/2)+1
    NQ0=4*K05 # no.of main qubits
    M=4*K# order of H-matrix
    TF_FLAG=False
    idx=0
    MAX_ITR=2**(NQ0-1) #10*1000
    #
    while(not(TF_FLAG) and idx<MAX_ITR):
        q=int_to_spin(idx,NQ0)
        A,B,C,D=construct_ABCD(q)
        W=construct_baumert_hall(A,B,C,D)
        print('H-',12*K,'; iteration->',idx,'of',MAX_ITR)
        TF_FLAG=is_hadamard(W)
        idx=idx+1
    ##whend
    # Save vector of first row when H is found
    if(TF_FLAG):
        #path_name='..//ABCD_VECTOR//'
        path_name='DATA//'
        write_ABCD(A,B,C,D,path_name)
    #--
    DG=np.matmul(W, W.transpose())
    # print('W-Matrix \n',W)
    print('Indicator matrix \n',DG)
    # display
    plt.imshow(DG.astype(float), cmap="gray")
    
