"""
STEP 2: Find the possible strings for Z
We find all strings of length which satisfy the following.
1) Sum of the elements of is +/- z
2) fz(t) <= 3n-1; t in (k*pi/100);1<=k<=100
3) z0=1; z_{n-1}=1
4) Z>=Zr
# assume n=6
then x^2+y^2+ 2(z^2+w^2) =6n-2 =34
list of sqr numbers: {1,4,9,16,25}
--------------------------------------------------------   
    x2      y2   2*z2          2*w         ttl   {z,w}
--------------------------------------------------------   
    0       0     2*4^2=32      2*1^2=2     34  {4,1}
    0       0     2*1^2=2       2*4^2=32    34  {1,4}
    3^2     5^2     0           0           34  {0,0}
    0       4^2     2*2^2=8      2*2^2=8    34  {2,2}
    4^2     0     2*2^2=8      2*2^2=8      34  {2,2}
    0       4     2*3^2        0            34  {3,0}  
    >> zeta={0,-1,1,-2,2,-3,3-4,4}
--------------------------------------------------------   

possible extended z->zz; and w->ww
put the original xy at the front
"""
import pandas as pd
from prb_turyn import *
from prb_generate_XYZW import *
from prb_london import *
##
# --
if __name__=='__main__':
    # --read data of XYZW* into array--
    # original length
    N_ORG=4
    N_BIT_ORG=4*N_ORG
    #
    # extension length
    N_XTD=8 #extend number-of-bits each
    N=N_ORG+N_XTD
    N_BIT_XTD=4*(N_ORG+N_XTD)
    #--- file definition ---
    path_name='DATA//'
    fname_org=path_name+'XYZW_star_iter_'+str(N_ORG)+'.txt'
    fname_xtd=path_name+'XYZW_star_xy_zzww_fill_'+\
        str(N_ORG)+'_'+str(N_ORG+N_XTD)+'.txt'
    ## prepare open file
    in_file = open(fname_org, "r")
    out_file_xy_zzww = open(fname_xtd, "w")
    ##
    num_lines = sum(1 for line in open(fname_org))
    #-- from X*,Y*,Z*,W*; insert spin-vector and check
    NQ=2*N_XTD
    # 
    NFIN=4*(N_ORG+N_XTD) #NR+4*N_XTD
    start_lag=int(NFIN/2)
    MAX_ITER=2**NQ
    cnt=0
    idx=0
    ## prepare all z and w, such that x**2+2y**2+z**2 + 2w**2 <6N-2
#    zeta_ohmega={0,-1,1,-2,2,-3,3,-4,4}  ## N=6
    zeta_ohmega=generate_zeta_ohmega(N)
    #-----
    # assume n=8 
    #then x^2+y^2+ 2(z^2+w^2) =6*8-2 = 46
    #list of sqr numbers < 46: {1,2,3,4,5,6}->{25,36; 1,4,9,16}-->{25,36; 2,8}
    #46=1+9+36(=2*9+2*0)  --->0,1,3,3
    #46=36+10=6^2+0+2+8 = 6^2+0+2*1^2+2*2^2 --{6,0,1,2}
    #==> zw={0,1,2,3,6}
    #-----
#    zeta_ohmega={0,-1,1,-2,2,-3,3,-6,6}  ## N=8
    # createcosinus table
    tbl_cos=gen_tbl_cos_t(N)
    ## estimate number of row
    NROW=0
    for tv in in_file:
        NROW+=1
    #--
    print('Total row in file=', NROW)
    in_file.close()
    ##-- reopen
    in_file = open(fname_org, "r")
    
    cnt_row=0
    for tv in in_file:
        # --
        v=hex_to_vec_xyzw0(tv, N_BIT_ORG)
        x,y,z,w0=extract_XYZW(v)
        ## extend z by applying some restrictions
        for itr in range(MAX_ITER):
            print('Extending zz,ww',itr,'of', MAX_ITER)
            q=int_to_spin(itr,NQ)
            zz=extend_spin_vector(z,q[0:int(NQ/2)])
            ww0=extend_spin_vector(w0,q[int(NQ/2):NQ])
            ww=ww0[0:len(ww0)-1]
            #
            TF_Z,TF_ZO, fz=is_zw_qualified(zz, zeta_ohmega, tbl_cos)
            TF_W,TF_WO, fw=is_zw_qualified(ww, zeta_ohmega, tbl_cos)
            if (TF_Z and TF_W) and (np.max(fz+fw)<3*N-1):
                # print('qualified')
                x_hex = vec_to_hex(x,len(x))
                y_hex = vec_to_hex(y,len(y))
                zz_hex= vec_to_hex(zz,len(zz))
                ww_hex= vec_to_hex(ww,len(ww))
                out_file_xy_zzww.write(str(x_hex)+str(y_hex)+\
                                 str(zz_hex)+str(ww_hex)+'\n')
                cnt=cnt+1
            #--end-if --
            idx=idx+1
        #-- end for itr --
        print('Processing', cnt_row, 'of',NROW)
        cnt_row +=1
    #--end for tv --
    print('qualified=', cnt, 'of', idx)

    #--
    in_file.close()
    out_file_xy_zzww.close()
