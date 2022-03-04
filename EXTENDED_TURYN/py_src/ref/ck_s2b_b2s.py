from prb_turyn import *
from prb_generate_XYZW import *
## cek bin to spin to bin to hex

def spin_to_bin1(v,K):
    #convert vector of{-1,+1} into vector of{0,1}
    N=len(v)
    b=np.zeros([K,1])
    if K<N:
        return(0)
    else:
        #
        for m in range(N):
             b[K-1-m]=abs(round((v[N-1-m]+0.5)/2))
            #b[K-1-m]=round((v[N-1-m].item()+0.5)/2)
        #
    #
    return(b)
##
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
##
#
def bin_to_int(s):
    # convert binary vector to integer
    K=len(s)
    val=0
    for m in range(K):
        val=val+s[K-1-m]*2**(m)
    #
    return(int(val))
#
##
if __name__=='__main__':    
    # --read data of XYZW* into array--
    path_name='DATA//'
    N_length=6
    fname=path_name+'XYZW_star_iter_'+str(N_length)+'.txt'
    mydat=np.loadtxt(fname)
    ##
    [NR,NC]=mydat.shape
    #spin->bin->int-hex
    idx=0
    for m in range(NC):
        tv0=mydat[:,m]
        tv=tv0[0:NR-1]
        tv_bin=spin_to_bin1(tv,NR-1)
        tv_int=bin_to_int(tv_bin)
        tv_hex=np.base_repr(tv_int,base=16)
#        print(idx,':', np.transpose(tv_bin),'->',tv_int,'->', tv_hex)
        print(idx,':', tv_int,'->', tv_hex)
        idx+=1
        
#
#    a=10
#    a_spin=int_to_spin(a,8)
#    print(a, a_spin)
#    
#    a_bin=spin_to_bin1(a_spin,8)
#    print(a, a_spin, a_bin)
#    #
#    v=[ 1,0,0,0,1,1,1,1, 1,0,0,0,1,1,1,1, 1,0,0,0,1,1,1,1, 1,0,0,0,1,1,1,1, \
#        1,0,0,0,1,1,1,1, 1,0,0,0,1,1,1,1, 1,0,0,0,1,1,1,1, 1,0,0,0,1,1,1,1, \
#        1,0,0,0,1,1,1,1, 1,0,0,0,1,1,1,1]
#    vint=bin_to_int(v)
#    print('v=', vint)
#    #
#    v2=int_to_spin(vint,80)
#    print('v2=',v2)
#    #
#    v_bin=np.binary_repr(vint, width=80)
#    print('vbin=', v_bin)
#    ## 
#    vhex=np.base_repr(vint, base=16)
#    print('vhex=',vhex)