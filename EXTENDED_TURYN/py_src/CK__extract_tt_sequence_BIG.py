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
    N_length=6
    in_fname=path_name+'XYZW_star_iter_'+str(N_length)+'.txt'
    out_fname=path_name+'XYZW_TTS_'+str(N_length)+'.txt'
    ## prepare open file
    in_file = open(in_fname, "r")
    out_file = open(out_fname, "w")
    #
    cnt=0
    idx=0
    N_BIT=4*N_length
    ##
    num_lines = sum(1 for line in open(in_fname))
    # for m in range(NC):
    for tv_hex in in_file:
        v=hex_to_vec_xyzw0(tv_hex, N_BIT)
        #
        x,y,z,w=extract_XYZW(v)
        if tt_check(x,y,z,w):
            out_file.write(str(tv_hex)) # no '\n' required, already there
            cnt=cnt+1
            # -------------
            X=np.copy(x)
            Y=np.copy(y)
            Z=np.copy(z)
            W=np.copy(w)
            # -------------
        #---end if---
        print('Processing ', idx,'of', num_lines)
        idx=idx+1
    #--------
    in_file.close()
    out_file.close()
    ##
    print('original',idx,'->',cnt)
    ## -- check one with constructing H-matrix
    if cnt>0:
        print('Hadamard matrix is found')
        NX=len(X)
        H=construct_hadamard(X,Y,Z,W[0:NX-1])
        DD=np.matmul(np.transpose(H),H)
        plt.imshow(DD)
    else:
        print('H-matrix NOT found!!!')
    #--
