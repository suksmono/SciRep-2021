import numpy as np
#

from prb_small_williamson import *
#
def write_ABCD_baumert_hall(A,B,C,D,path_name):
    M,N=A.shape
    M_ABCD=[A[0,:], B[0,:], C[0,:], D[0,:]]
    ABCD_name=path_name+'ABCD_'+str(12*N)+'.txt'
    write_matrix(ABCD_name, M_ABCD)
#
"""
    --- Baumert-Hall ---
    H=  [ A  A  A  B  -B  C -C -D   B  C -D -D]
        [ A -A  B -A  -B -D  D -C  -B -D -C -C]
        [ A -B -A  A  -D  D -B  B  -C -D  C -C]
        [ B  A -A -A   D  D  D  C   C -B -B -C]  
        
        [ B -D  D  D   A  A  A  C  -C  B -C  B] 
        [ B  C -D  D   A -A  C -A  -D  C  B -B]
        [ D -C  B -B   A -C -A  A   B  C  D -D] 
        [-C -D -C -D   C  A -A -A  -D  B -B -B]  
        
        [ D -C -B -B  -B  C  C -D   A  A  A  D]  
        [-D -B  C  C   C  B  B -D   A -A  D -A] 
        [ C -B -C  C   D -B -D -B   A -D -A  A] 
        [-C -D -D  C  -C -B  B  B   D  A -A -A]   
"""
def construct_baumert_hall(A,B,C,D):
    """
    input  : submatrix A,B,C,D
    output : Baumert-Hall matrix
    """
    M,N=A.shape
    W=np.empty([12*M,12*N], dtype=object)

    '''
     # 1.[ A  A  A  B    -B  C -C -D     B  C -D -D]
   '''
    # 1.block-1: A  A  A  B
    W[0:M,0:M]=     A[0:M,0:M]
    W[0:M,M:2*M]=   A[0:M,0:M]
    W[0:M,2*M:3*M]= A[0:M,0:M]
    W[0:M,3*M:4*M]= B[0:M,0:M]

    # 1.block-2: -B  C -C -D
    W[0:M,4*M:5*M]= -B[0:M,0:M]
    W[0:M,5*M:6*M]=  C[0:M,0:M]
    W[0:M,6*M:7*M]= -C[0:M,0:M]
    W[0:M,7*M:8*M]= -D[0:M,0:M]

    # 1.block-3: B  C -D -D
    W[0:M,8*M:9*M]=    B[0:M,0:M]
    W[0:M,9*M:10*M]=   C[0:M,0:M]
    W[0:M,10*M:11*M]= -D[0:M,0:M]
    W[0:M,11*M:12*M]= -D[0:M,0:M]

    '''
    # 2.[ A -A  B -A    -B -D  D -C    -B -D -C -C]
    '''
    #: A -A  B -A
    W[M:2*M,0:M]=      A[0:M,0:M]
    W[M:2*M,M:2*M]=   -A[0:M,0:M]
    W[M:2*M,2*M:3*M]=  B[0:M,0:M]
    W[M:2*M,3*M:4*M]= -A[0:M,0:M]

    #: -B -D  D -C 
    W[M:2*M,4*M:5*M]= -B[0:M,0:M]
    W[M:2*M,5*M:6*M]= -D[0:M,0:M]
    W[M:2*M,6*M:7*M]=  D[0:M,0:M]
    W[M:2*M,7*M:8*M]= -C[0:M,0:M]

    #: -B -D -C -C
    W[M:2*M,8*M:9*M]=   -B[0:M,0:M]
    W[M:2*M,9*M:10*M]=  -D[0:M,0:M]
    W[M:2*M,10*M:11*M]= -C[0:M,0:M]
    W[M:2*M,11*M:12*M]= -C[0:M,0:M]


    # 3.[ A -B -A  A    -D  D -B  B    -C -D  C -C]
    #: A -B -A  A 
    W[2*M:3*M,0:M]=       A[0:M,0:M]
    W[2*M:3*M,M:2*M]=    -B[0:M,0:M]
    W[2*M:3*M,2*M:3*M]=  -A[0:M,0:M]
    W[2*M:3*M,3*M:4*M]=   A[0:M,0:M]

    #: -D  D -B  B  
    W[2*M:3*M,4*M:5*M]= -D[0:M,0:M]
    W[2*M:3*M,5*M:6*M]=  D[0:M,0:M]
    W[2*M:3*M,6*M:7*M]= -B[0:M,0:M]
    W[2*M:3*M,7*M:8*M]=  B[0:M,0:M]

    #: -C -D  C -C
    W[2*M:3*M,8*M:9*M]=   -C[0:M,0:M]
    W[2*M:3*M,9*M:10*M]=  -D[0:M,0:M]
    W[2*M:3*M,10*M:11*M]=  C[0:M,0:M]
    W[2*M:3*M,11*M:12*M]= -C[0:M,0:M]
    
    
    # 4.[ B  A -A -A     D  D  D  C     C -B -B -C]  
    #:B  A -A -A 
    W[3*M:4*M,0:M]=      B[0:M,0:M]
    W[3*M:4*M,M:2*M]=    A[0:M,0:M]
    W[3*M:4*M,2*M:3*M]= -A[0:M,0:M]
    W[3*M:4*M,3*M:4*M]= -A[0:M,0:M]

    #: D  D  D  C  
    W[3*M:4*M,4*M:5*M]=  D[0:M,0:M]
    W[3*M:4*M,5*M:6*M]=  D[0:M,0:M]
    W[3*M:4*M,6*M:7*M]=  D[0:M,0:M]
    W[3*M:4*M,7*M:8*M]=  C[0:M,0:M]

    #:  C -B -B -C
    W[3*M:4*M,8*M:9*M]=    C[0:M,0:M]
    W[3*M:4*M,9*M:10*M]=  -B[0:M,0:M]
    W[3*M:4*M,10*M:11*M]= -B[0:M,0:M]
    W[3*M:4*M,11*M:12*M]= -C[0:M,0:M]  
    #
    
    # 5.[ B -D  D  D     A  A  A  C    -C  B -C  B] 
    #: B -D  D  D
    W[4*M:5*M,0:M]=      B[0:M,0:M]
    W[4*M:5*M,M:2*M]=   -D[0:M,0:M]
    W[4*M:5*M,2*M:3*M]=  D[0:M,0:M]
    W[4*M:5*M,3*M:4*M]=  D[0:M,0:M]

    #: A  A  A  C 
    W[4*M:5*M,4*M:5*M]=  A[0:M,0:M]
    W[4*M:5*M,5*M:6*M]=  A[0:M,0:M]
    W[4*M:5*M,6*M:7*M]=  A[0:M,0:M]
    W[4*M:5*M,7*M:8*M]=  C[0:M,0:M]

    #:  -C  B -C  B
    W[4*M:5*M,8*M:9*M]=   -C[0:M,0:M]
    W[4*M:5*M,9*M:10*M]=   B[0:M,0:M]
    W[4*M:5*M,10*M:11*M]= -C[0:M,0:M]
    W[4*M:5*M,11*M:12*M]=  B[0:M,0:M]  
    #
    
    # 6.[ B  C -D  D     A -A  C -A    -D  C  B -B]
    #: B  C -D  D
    W[5*M:6*M,0:M]=      B[0:M,0:M]
    W[5*M:6*M,M:2*M]=    C[0:M,0:M]
    W[5*M:6*M,2*M:3*M]= -D[0:M,0:M]
    W[5*M:6*M,3*M:4*M]=  D[0:M,0:M]

    #:  A -A  C -A  
    W[5*M:6*M,4*M:5*M]=  A[0:M,0:M]
    W[5*M:6*M,5*M:6*M]= -A[0:M,0:M]
    W[5*M:6*M,6*M:7*M]=  C[0:M,0:M]
    W[5*M:6*M,7*M:8*M]= -A[0:M,0:M]

    #: -D  C  B -B
    W[5*M:6*M,8*M:9*M]=   -D[0:M,0:M]
    W[5*M:6*M,9*M:10*M]=   C[0:M,0:M]
    W[5*M:6*M,10*M:11*M]=  B[0:M,0:M]
    W[5*M:6*M,11*M:12*M]= -B[0:M,0:M]  
    #    

    # 7.[ D -C  B -B     A -C -A  A     B  C  D -D] 
    #: D -C  B -B
    W[6*M:7*M,0:M]=      D[0:M,0:M]
    W[6*M:7*M,M:2*M]=   -C[0:M,0:M]
    W[6*M:7*M,2*M:3*M]=  B[0:M,0:M]
    W[6*M:7*M,3*M:4*M]= -B[0:M,0:M]

    #:  A -C -A  A   
    W[6*M:7*M,4*M:5*M]=  A[0:M,0:M]
    W[6*M:7*M,5*M:6*M]= -C[0:M,0:M]
    W[6*M:7*M,6*M:7*M]= -A[0:M,0:M]
    W[6*M:7*M,7*M:8*M]=  A[0:M,0:M]

    #: B  C  D -D
    W[6*M:7*M,8*M:9*M]=    B[0:M,0:M]
    W[6*M:7*M,9*M:10*M]=   C[0:M,0:M]
    W[6*M:7*M,10*M:11*M]=  D[0:M,0:M]
    W[6*M:7*M,11*M:12*M]= -D[0:M,0:M]  
    #    

    # 8.[-C -D -C -D     C  A -A -A    -D  B -B -B]  
    #: -C -D -C -D 
    W[7*M:8*M,0:M]=     -C[0:M,0:M]
    W[7*M:8*M,M:2*M]=   -D[0:M,0:M]
    W[7*M:8*M,2*M:3*M]= -C[0:M,0:M]
    W[7*M:8*M,3*M:4*M]= -D[0:M,0:M]

    #: C  A -A -A   
    W[7*M:8*M,4*M:5*M]=  C[0:M,0:M]
    W[7*M:8*M,5*M:6*M]=  A[0:M,0:M]
    W[7*M:8*M,6*M:7*M]= -A[0:M,0:M]
    W[7*M:8*M,7*M:8*M]= -A[0:M,0:M]

    #: -D  B -B -B
    W[7*M:8*M,8*M:9*M]=   -D[0:M,0:M]
    W[7*M:8*M,9*M:10*M]=   B[0:M,0:M]
    W[7*M:8*M,10*M:11*M]= -B[0:M,0:M]
    W[7*M:8*M,11*M:12*M]= -B[0:M,0:M]  
    #    
            
    # 9.[ D -C -B -B    -B  C  C -D     A  A  A  D]  
    #: D -C -B -B
    W[8*M:9*M,0:M]=      D[0:M,0:M]
    W[8*M:9*M,M:2*M]=   -C[0:M,0:M]
    W[8*M:9*M,2*M:3*M]= -B[0:M,0:M]
    W[8*M:9*M,3*M:4*M]= -B[0:M,0:M]

    #: -B  C  C -D   
    W[8*M:9*M,4*M:5*M]= -B[0:M,0:M]
    W[8*M:9*M,5*M:6*M]=  C[0:M,0:M]
    W[8*M:9*M,6*M:7*M]=  C[0:M,0:M]
    W[8*M:9*M,7*M:8*M]= -D[0:M,0:M]

    #: A  A  A  D
    W[8*M:9*M,8*M:9*M]=    A[0:M,0:M]
    W[8*M:9*M,9*M:10*M]=   A[0:M,0:M]
    W[8*M:9*M,10*M:11*M]=  A[0:M,0:M]
    W[8*M:9*M,11*M:12*M]=  D[0:M,0:M]  
    #    
    
    #10.[-D -B  C  C     C  B  B -D     A -A  D -A] 
    #: -D -B  C  C 
    W[9*M:10*M,0:M]=     -D[0:M,0:M]
    W[9*M:10*M,M:2*M]=   -B[0:M,0:M]
    W[9*M:10*M,2*M:3*M]=  C[0:M,0:M]
    W[9*M:10*M,3*M:4*M]=  C[0:M,0:M]

    #: C  B  B -D   
    W[9*M:10*M,4*M:5*M]=  C[0:M,0:M]
    W[9*M:10*M,5*M:6*M]=  B[0:M,0:M]
    W[9*M:10*M,6*M:7*M]=  B[0:M,0:M]
    W[9*M:10*M,7*M:8*M]= -D[0:M,0:M]

    #: A -A  D -A
    W[9*M:10*M,8*M:9*M]=    A[0:M,0:M]
    W[9*M:10*M,9*M:10*M]=  -A[0:M,0:M]
    W[9*M:10*M,10*M:11*M]=  D[0:M,0:M]
    W[9*M:10*M,11*M:12*M]= -A[0:M,0:M]  
    #  
    
    #11.[ C -B -C  C     D -B -D -B     A -D -A  A] 
    #: C -B -C  C 
    W[10*M:11*M,0:M]=      C[0:M,0:M]
    W[10*M:11*M,M:2*M]=   -B[0:M,0:M]
    W[10*M:11*M,2*M:3*M]= -C[0:M,0:M]
    W[10*M:11*M,3*M:4*M]=  C[0:M,0:M]

    #: D -B -D -B   
    W[10*M:11*M,4*M:5*M]=  D[0:M,0:M]
    W[10*M:11*M,5*M:6*M]= -B[0:M,0:M]
    W[10*M:11*M,6*M:7*M]= -D[0:M,0:M]
    W[10*M:11*M,7*M:8*M]= -B[0:M,0:M]

    #: A -D -A  A
    W[10*M:11*M,8*M:9*M]=    A[0:M,0:M]
    W[10*M:11*M,9*M:10*M]=  -D[0:M,0:M]
    W[10*M:11*M,10*M:11*M]= -A[0:M,0:M]
    W[10*M:11*M,11*M:12*M]=  A[0:M,0:M]  
    #    
    
    #12.[-C -D -D C    -C -B  B  B     D  A -A -A]   
    #: -C -D -D C
    W[11*M:12*M,0:M]=     -C[0:M,0:M]
    W[11*M:12*M,M:2*M]=   -D[0:M,0:M]
    W[11*M:12*M,2*M:3*M]= -D[0:M,0:M]
    W[11*M:12*M,3*M:4*M]=  C[0:M,0:M]

    #: -C -B  B  B  
    W[11*M:12*M,4*M:5*M]= -C[0:M,0:M]
    W[11*M:12*M,5*M:6*M]= -B[0:M,0:M]
    W[11*M:12*M,6*M:7*M]=  B[0:M,0:M]
    W[11*M:12*M,7*M:8*M]=  B[0:M,0:M]

    #:  D  A -A -A
    W[11*M:12*M,8*M:9*M]=    D[0:M,0:M]
    W[11*M:12*M,9*M:10*M]=   A[0:M,0:M]
    W[11*M:12*M,10*M:11*M]= -A[0:M,0:M]
    W[11*M:12*M,11*M:12*M]= -A[0:M,0:M]  
    #    

    # #
    return(W)
#