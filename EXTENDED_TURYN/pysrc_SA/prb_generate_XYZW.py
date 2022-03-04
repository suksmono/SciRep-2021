#
def generate_XYZW_4():
    #W=np.zeros([4,1])
    #
    x=[1,-1,1,-1]
    y=[1,-1,-1,-1]
    z=[1,-1,-1,1]
    w=[1,1,1]
    return(x,y,z,w)
#
def generate_XYZW_6():
    #W=np.zeros([4,1])
    #
    # x=[1,-1,1,1,-1,1]
    # y=[1,1,1,-1,-1,1]
    # z=[1,-1,1,1,1,-1]
    # w=[1,-1,-1,-1,-1]

    #-- normalize --
    x=[ 1,+1, 1,-1,-1,-1]
    y=[ 1,-1, 1,+1,-1,-1] 
    z=[ 1,+1, 1,-1, 1,+1]
    w=[ 1,+1,-1,+1,-1]

    return(x,y,z,w)
		
def generate_XYZW_8():
    #W=np.zeros([4,1])
    #
    x=[1, 1,-1, 1, -1, 1,-1,1]
    y=[1,-1,-1,-1, -1,-1,-1,1]
    z=[1,-1,-1, 1,  1, 1, 1,-1]
    w=[1, 1, 1,-1,  1, 1,-1]
    return(x,y,z,w)
#
def generate_XYZW_10():
    x=[+1,+1,-1,-1,+1, +1,-1,+1,+1,-1]
    y=[+1,-1,+1,+1,-1, +1,-1,+1,+1,-1]
    z=[+1,+1,+1,+1,+1, +1,-1,-1,-1,+1]
    w=[+1,+1,-1,+1,-1, +1,+1,+1,-1]
    return(x,y,z,w)
#
def generate_ZW_36():
    #from London's thesis
    #x=[+++++-***-++++-]
    #y=[++---+***-+++--]
    z=[+1,+1,+1,-1,+1,-1,+1,-1,-1,+1,+1,+1,-1,+1,-1,+1,+1,+1,\
       +1,-1,-1,+1,-1,+1,+1,+1,-1,-1,-1,-1,-1,+1,-1,-1,+1,+1]
    w=[+1,+1,+1,+1,+1,-1,-1,-1,+1,+1,-1,-1,-1,+1,-1,-1,-1,+1,\
       -1,-1,+1,-1,-1,-1,+1,-1,+1,+1,-1,+1,+1,-1,+1,-1,-1]
    return(z,w)
#
def generate_XYZW_36():
    # from kharagani-Rezaie paper428
    X=[+1,+1,+1,-1,-1,-1,-1,+1,+1,-1,+1,-1,+1,-1,-1,-1,-1,-1,\
       +1,+1,+1,+1,-1,+1,+1,-1,+1,+1,+1,+1,-1,-1,-1,-1,+1,-1]
    #
    Y=[+1,-1,+1,+1,+1,+1,+1,-1,-1,+1,-1,+1,-1,-1,+1,-1,-1,+1,\
       +1,-1,-1,+1,+1,+1,+1,-1,+1,+1,+1,+1,-1,-1,-1,+1,+1,-1]
    #
    Z = [+1, -1, +1, +1, +1, +1, +1, -1, +1, -1, -1,\
         +1, +1, +1, +1, -1, +1, +1, +1, -1, +1, +1,\
             -1, -1, +1, +1, +1, -1, +1, -1, -1, +1,\
                 -1, -1, -1, +1]
    W = [+1, +1, +1, -1, +1, -1, -1, -1, -1, -1, +1,\
         +1, -1, -1, +1, -1, +1, +1, +1, -1, -1, +1,\
             -1, +1, -1, +1, +1, +1, -1, +1, +1, +1,\
                 +1, -1, +1]
    return(X,Y,Z,W)