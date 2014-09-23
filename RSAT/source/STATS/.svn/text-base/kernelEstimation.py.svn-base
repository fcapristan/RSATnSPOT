# python script to calculate the KERNEL DENSITY ESTIMATION
# Base on RANGE SAFETY APPLICATIONN OF KERNEL DENSITY ESTIMATION
# from the AUSTRALIAN GOVERNMENT-DEPARTMENT OF DEFENCE
import numpy as np



def calculateS(xy,nSamples):
    
    xybar = np.mean(xy,axis=0)
    xbar = xybar[0]
    ybar = xybar[1]
  #  print xbar,ybar
    S11 = 0.0 # initializing S values
    S12 = 0.0
    S22 = 0.0
   # print xy
    for index in range(0,nSamples):
    #    print index
        S11 = S11 + (xy[index,0]-xbar)**2
        S12 = S12 + (xy[index,0]-xbar)*(xy[index,1]-ybar)
        S22 = S22 + (xy[index,1]-ybar)**2
    S = np.array([[S11,S12],[S12,S22]])
    S = 1.0/(nSamples-1.0)*S
    return S

def calculateU(xy,nSamples):
    S = calculateS(xy,nSamples)
    D,U = np.linalg.eig(S)
        # print D
        # print U
    return U

def PQframe(xy,Umatrix):
    xyMatrix = np.matrix(xy)
    PQ = np.array(xyMatrix*Umatrix)
    return PQ


def calculateH2(xy):
    nSamples = np.size(xy,axis=0)
    U = calculateU(xy,nSamples)
    Umatrix = np.matrix(U)
    PQ = PQframe(xy,Umatrix)
    x = PQ[:,0]
    y = PQ[:,1]
    sigma1 = np.std(x)
    sigma2 = np.std(y)
    IQR1 = IQR(x)
    IQR2 = IQR(y)
    # print U
    h1 = 1.06*min(sigma1,IQR1/1.34)*nSamples**(-1.0/5.0)
    h2 = 1.06*min(sigma2,IQR2/1.34)*nSamples**(-1.0/5.0)
#    print h1,h2
#    print Umatrix
    H2 = Umatrix*np.matrix([[h1**2,0],[0,h2**2]])*np.linalg.inv(Umatrix)
#    print H2
    return H2,nSamples
    

def IQR(x):
    import scipy.stats.mstats as st
    quarts = st.idealfourths(x)
    val = quarts[1]-quarts[0]
    return val
'''
# This section commented out because other routine is more efficient
def fhatKernel3(xVector,yVector,xy):
    H2,nSamples = calculateH2(xy) #H2 is already of matrix type
    H2inv = np.linalg.inv(H2)
    Xi = xy[:,0]
    Yi = xy[:,1]
    xLen = len(xVector)
    yLen = len(yVector)
    fMatrix = np.zeros((yLen,xLen))
    for indexX in range(0,xLen):
        for indexY in range(0,yLen):
            fMatrix[indexY,indexX] = fsumhat(xVector[indexX],yVector[indexY],Xi,Yi,H2inv,nSamples)
    fMatrix = 1.0/(2.0*np.pi*nSamples*(np.linalg.det(H2))**0.5)*fMatrix
    return fMatrix
'''
   





def fsumhat(xl,yl,Xi,Yi,H2inv,nSamples):
    fsum = 0.0
    h11 = H2inv[0,0]
    h22 = H2inv[1,1]
    h12 = H2inv[0,1] # h12 is the same as h21

    for index in range(0,nSamples):
    #    print Xi[index]
       # xMatrix = np.matrix([xl-Xi[index],yl-Yi[index]]) # transpose becuase of the way python reads the values inside Xi and Yi
        x1 = xl-Xi[index]
        x2 = yl-Yi[index]
      #  matrixMult = xMatrix*H2inv*xMatrix.T
        matrixMult = x1*x1*h11+2*x1*x2*h12+x2*x2*h22
      #  print matrixMult
        fsum = np.exp(-0.5*matrixMult)+fsum
#    print fsum
    return fsum

def fsumhatFAST(xl,yl,nhits,Xi,Yi,H2inv):
    fsum = 0.0
    h11 = H2inv[0,0]
    h22 = H2inv[1,1]
    h12 = H2inv[0,1] # h12 is the same as h21

    for index in range(0,len(nhits)):
        x1 = xl-Xi[index]
        x2 = yl-Yi[index]
  
        matrixMult = x1*x1*h11+2*x1*x2*h12+x2*x2*h22
        fsum = nhits[index]*np.exp(-0.5*matrixMult)+fsum
    return fsum
    

def fhatKernelFAST(xVector,yVector,xy,dx,dy):
    xmin = np.min(xy[:,0])-.1*dx
    xmax = np.max(xy[:,0])+.1*dx
    ymin = np.min(xy[:,1])-.1*dy
    ymax = np.max(xy[:,1])+.1*dy
    xybins = hits2bins(xy,xmin,xmax,ymin,ymax,dx,dy)
    H2,nSamples = calculateH2(xy) #H2 is already of matrix type
    H2inv = np.linalg.inv(H2)
    Xi = xybins[:,0]
    Yi = xybins[:,1]
    nhits = xybins[:,2]
    XX,YY=np.meshgrid(xVector,yVector)
    fMatrix = fsumhatFAST(XX,YY,nhits,Xi,Yi,H2inv)
    fMatrix = 1.0/(2.0*np.pi*nSamples*(np.linalg.det(H2))**0.5)*fMatrix
    return fMatrix



def fhatKernel(xVector,yVector,xy):
    H2,nSamples = calculateH2(xy) #H2 is already of matrix type
    H2inv = np.linalg.inv(H2)
    Xi = xy[:,0]
    Yi = xy[:,1]
    XX,YY=np.meshgrid(xVector,yVector)
    fMatrix = fsumhat(XX,YY,Xi,Yi,H2inv,nSamples)
    fMatrix = 1.0/(2.0*np.pi*nSamples*(np.linalg.det(H2))**0.5)*fMatrix
    return fMatrix


def fhatKernelV2(XX,YY,xy):
    H2,nSamples = calculateH2(xy) #H2 is already of matrix type
    H2inv = np.linalg.inv(H2)
    Xi = xy[:,0]
    Yi = xy[:,1]
    fMatrix = fsumhat(XX,YY,Xi,Yi,H2inv,nSamples)
    fMatrix = 1.0/(2.0*np.pi*nSamples*(np.linalg.det(H2))**0.5)*fMatrix
    return fMatrix


def hits2bins(xy,xmin,xmax,ymin,ymax,dx,dy):
    xvector = np.arange(xmin,xmax,dx)
    yvector = np.arange(ymin,ymax,dy)
    bins = np.zeros((len(xy),3))
    xindex = np.floor((xy[:,0]-xmin)/dx)
    yindex = np.floor((xy[:,1]-ymin)/dy)
    index = len(xy)-1
    xylen = len(xy)
    binsIndex =0
    print xindex
    while index>=0:
        xval = xindex[index]
        yval = yindex[index]

        xindex = np.delete(xindex,index)
        yindex = np.delete(yindex,index)
        print xval
        
        bins[binsIndex,:]=[xvector[xval],yvector[yval],1]
        index = index - 1
        for index2 in range(index,-1,-1):

            if xindex[index2]==xval and yindex[index2]==yval:
 
                bins[binsIndex,2]=bins[binsIndex,2]+1
                xindex = np.delete(xindex,index2)
                yindex = np.delete(yindex,index2)
                index = index - 1
        binsIndex = binsIndex+1
    bins = bins[0:binsIndex,:]
    return bins

        
        
            


    
