###########################################################################
# File: CreateChebyshevD.py
# 
# Forms the Chebyshev differentiation matrix, cosine-spaced vector of independent variables
# 
#
# Created by: Kim Shish
# Created on: 02/01/2011
#
# Major Modifications by: Francisco Capristan
# Major Modification on : 07/21/2012
#
#
#           ---------- Inputs ----------
# N          = [0 N] Number of points desired in trajectory vector/ order
#              of Chebyshev polynomial
#
#           ---------- Outputs ----------
# x      = cosine-spaced vector
# D[N,N] = Derivative matrix
# err        = Error code
#
###########################################################################

import numpy as np

def ChebyshevMatrices(N):

    err = 0
    x = np.cos(np.arange(N)*np.pi/(N-1.0))

    if N<1:
        print '\n'
        print 'Error:N out of bounds'
        err = 1
        return
        

    D=np.zeros((N,N))

    for i in np.arange(N):
        for j in np.arange(N):
            
            ci = 1.
            cj = 1.
    
            if (i==0) or (i==N-1):
                ci = 2.
            elif (j==0)or(j==N-1):  
                cj = 2.
            
        
            if (i==j):
                if (i!=0) and (i!=N-1):
                    
                    D[i,j] = -0.5*x[j]/(1.-x[j]**2)
            else : 
                D[i,j] = (ci/cj)*((-1.)**(i+j))/(x[i]-x[j])
        
    D[0,0] = (2.*(N-1.)**2 + 1.)/6.
    D[N-1,N-1] = -D[0,0]
    D[0,N-1] = .5*D[0,N-1]
    D[N-1,0] = .5*D[N-1,0]

# remap x
    x = .5 * (1 - x)
    D = -2.*D
            


    return (x,D,err)


