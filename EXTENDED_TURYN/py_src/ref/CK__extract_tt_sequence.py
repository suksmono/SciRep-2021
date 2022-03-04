"""
----------------------------------------------------------------------------
Iteratively construct XYZW*; eg. X=#[x0,x1,*,*,*,*,x6,x7]
1. Run CK__gen_xyzw_iter.py for initial length |-> file: XYZW_star_iter_N.txt
2. Iteratively run CK__extend_xyzw_star_iter.py
3. Cek, there should be at least TT-seq, by calling 
   CK__extract_tt_sequence.py
----------------------------------------------------------------------------
"""
### ---
import pandas as pd
from prb_turyn import *
from prb_generate_XYZW import *
### ---
if __name__=='__main__':
    # --read data of XYZW* into array--
    path_name='DATA//'
    N_length=10
    fname=path_name+'XYZW_star_iter_'+str(N_length)+'.txt'
    mydat=np.loadtxt(fname)
    #
    NR,NC=mydat.shape
    #
    cnt=0
    idx=0
    tts_data=np.zeros([NR,NC]) # worst case: they all qualified
    for m in range(NC):
        v=mydat[:,m]
        #
        x,y,z,w=extract_XYZW(v)
        if tt_check(x,y,z,w):
            tts_data[:,cnt]=v
            cnt=cnt+1
            # -------------
            X=np.copy(x)
            Y=np.copy(y)
            Z=np.copy(z)
            W=np.copy(w)
            # -------------
        #
        print('Processing ', idx,'of', NC)
        idx=idx+1
    #--------
    ##
    print('original',NC,'->',cnt)
    ## save data
    path_name='DATA//'
    fname='XYZW_TTS_'+str(N_length)+'.txt'
    np.savetxt(path_name+fname,tts_data[:,0:cnt])
    ## -- check one with constructing H-matrix
    if cnt>0:
        print('Hadamard matrix is found')
        NX=len(X)
        H=construct_hadamard(X,Y,Z,W[0:NX-1])
        DD=np.matmul(np.transpose(H),H)
        plt.imshow(DD)
    else:
        print('H-matrix NOT found!!!')
