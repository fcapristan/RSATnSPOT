import numpy as np
import sys
sys.path.append("../STATS/kernelQuadCpp")
import meancov
import pdfCoordTrans as pdf
import numpy as np
import lonlatChecks
import serialKernel as SK
import meancov
import statsPython
import matplotlib.pyplot as plt


def calculate(lonlat1,lonlat2,n1,n2,delta,nsigma,pdfoption):
    # n2 is the number of samples in the second probability...n2>n1

    lonMesh,latMesh,ZZpdf2,XX,YY,xlen,ylen,transformDetails = pdf.getPDF(lonlat2,n2,delta,nsigma,pdfoption)
        #since all the pieces are independent from each other, the first n1 samples can be simply used
        
    R1 = transformDetails[0]
    R2 = transformDetails[1]
    R3 = transformDetails[2]
    ind = np.random.permutation(n1)
    ZZpdf1 = getPDF4KLDiv(lonlat1[ind,:],n1,pdfoption,XX,YY,R1,R2,R3,xlen,ylen)
    
    A = delta**2
    Log1 = np.log(ZZpdf2)-np.log(ZZpdf1)
    #print ZZpdf2
    #print ZZpdf1
    #print ZZpdf2/ZZpdf1
    intMult = ZZpdf2*Log1
    KLDist = np.sum(intMult)
    return(KLDist)



def calculateKLMain(refLonLat,lonlatList,delta,nsigma,pdfoption):
    
    # NEED RECTANGULAR MESH IN ORIGINAL FRAME
    
    xMesh,yMesh,xlen,ylen= statsPython.areaOfInterestKDE(refLonLat,delta) # 

    A = delta**2

    lonMesh,latMesh,ZZpdfNMain = calculatePDF(refLonLat,xMesh,yMesh,pdfoption)
    nlist = []
    KLlist = []
    ZZpdfNMain = ZZpdfNMain*A
    ZZpdfNMain[ZZpdfNMain==0] = 10.**-50.
    KLsurf = []
    #print 'part1 done'
    for index in range(len(lonlatList)):
        lonlat = lonlatList[index]
        nSamples,colvals = np.shape(lonlat)
        lonMesh,latMesh,ZZpdfC = calculatePDF(lonlat,xMesh,yMesh,pdfoption)
        ZZpdfCarea = A*ZZpdfC
        #KL dist
        ZZpdfCarea[ZZpdfCarea==0] = 10.**-50.
        Log1 = np.log(ZZpdfNMain) - np.log(ZZpdfCarea)
        intMult = ZZpdfNMain*Log1
        KLDist = np.sum(intMult)
        KLlist.append(KLDist)
        nlist.append(nSamples)
        KLsurf.append([lonMesh,latMesh,ZZpdfC])
    return (nlist,KLlist,KLsurf)


def calculatePDF(refLonLat,xMesh,yMesh,pdfoption):
    
    '''
    plt.figure()
    plt.plot(refLonLat[:,0],refLonLat[:,1],'x')
    plt.show()
    exit()
    '''
    xyz,x,y,z = pdf.getxyz(refLonLat)
    xmean = np.mean(x)
    ymean = np.mean(y)
    zmean = np.mean(z)
    mag = (xmean**2+ymean**2+zmean**2)#**.5
    lonlat0 = pdf.getlonlat(xmean/mag,ymean/mag,zmean/mag)
    lon0 = lonlat0[0,0]
    lat0 = lonlat0[0,1]
    beta = 0 # initial guess
    nSamples,colvals = np.shape(refLonLat)
    row,col = np.shape(xMesh)
    beta = pdf.optimizeHeadingAngle(refLonLat,nSamples,lon0,lat0,beta)
    R1 = pdf.Rotz(lon0*np.pi/180.) # generating rotation matrices
    R2 = pdf.Roty(-lat0*np.pi/180.)
    R3 = pdf.Rotx(beta*np.pi/180.)
    transformDetails = [R1,R2,R3,lon0,lat0,beta]
    
    Uint = np.dot(R2,R1)
    U  = np.dot(R3,Uint)
    lonlatPframe = pdf.originalFrame2Pframe(refLonLat,U)
    
    lonlatMesh = np.zeros((np.size(xMesh),2))
    lonlatMesh[:,0] = np.reshape(xMesh,np.size(xMesh))
    lonlatMesh[:,1] = np.reshape(yMesh,np.size(yMesh))

    lonlatPMesh = pdf.originalFrame2Pframe(lonlatMesh,U)
 
    ylen,xlen = np.shape(xMesh)
    lonMeshP = np.reshape(lonlatPMesh[:,0],[ylen,xlen])
    latMeshP = np.reshape(lonlatPMesh[:,1],[ylen,xlen])

    if pdfoption.lower()=='kde' or pdfoption.lower()=='kernel':        
        ZZpdfPframe = SK.serialK(int(nSamples),lonlatPframe,row,col,lonMeshP,latMeshP,.0,0)
    elif pdfoption=='normal' or pdfoption=='gaussian' or pdfoption=='gauss':
        ZZpdfPframe = meancov.normalbivariate(mean,covar,lonMeshP,latMeshP,col,row) #fortran version    
    '''
    plt.figure()
    plt.contour(lonMeshP,latMeshP,ZZpdfPframe,30)
    plt.plot(lonlatPframe[:,0],lonlatPframe[:,1],'x')
    
    plt.figure()
    plt.plot(lonlatPMesh[:,0],lonlatPMesh[:,1],'x')
    exit()
    '''
    
    
    pqrVec ,pVec,qVec,rVec= pdf.getxyz(lonlatPMesh) # mesh locations in P frame
    
    pdfN = pdf.transformPDF(U,pqrVec,np.pi/180*lonlatMesh[:,0],np.pi/180*lonlatMesh[:,1],np.reshape(ZZpdfPframe,xlen*ylen))
    ZZpdfN = np.reshape(pdfN,[ylen,xlen])
    lonMesh = np.reshape(lonlatMesh[:,0],[ylen,xlen])
    latMesh = np.reshape(lonlatMesh[:,1],[ylen,xlen])
    
    '''
    plt.figure()
    plt.contour(lonMeshP,latMeshP,ZZpdfPframe,30)
    plt.plot(refLonLat[:,0],refLonLat[:,1],'x')
    plt.show()
    
    plt.figure()
    plt.contour(lonMesh,latMesh,ZZpdfN,30)
    plt.plot(refLonLat[:,0],refLonLat[:,1],'x')
    plt.show()
    '''
    return (lonMesh,latMesh,ZZpdfN)




def calculateMLS(lonlat1,lonlat2,n1,n2,delta,nsigma,pdfoption):
    # n2 is the number of samples in the second probability...n2>n1
    
    lonMesh,latMesh,ZZpdf2,XX,YY,xlen,ylen,transformDetails = pdf.getPDF(lonlat2,n2,delta,nsigma,pdfoption)
    #since all the pieces are independent from each other, the first n1 samples can be simply used
    
    R1 = transformDetails[0]
    R2 = transformDetails[1]
    R3 = transformDetails[2]
    n1Total,cols = np.shape(lonlat1)
    ind = np.random.permutation(n1Total)

    ZZpdf1 = getPDF4KLDiv(lonlat1[ind[0:n1],:],n1,pdfoption,lonMesh,latMesh,R1,R2,R3,xlen,ylen)
    A = delta**2.0
    
    #intSum = (ZZpdf1-ZZpdf2)**2/(ZZpdf2*(1./A - ZZpdf2))
    #MLSval = 1./(xlen*ylen)*np.sum(np.log10(intSum))
    Ptest = ZZpdf1*A
    Pover = ZZpdf2*A
    maxval = np.max(Pover)
    distmin = 1e-10 # range of probabilities from distmin*maxval to maxval. This ensures divide by zero errors
    minval = distmin*maxval
    Pover[Pover<minval] = minval
    Ptest[Ptest<minval] = minval
    Log1 = np.log(Pover)-np.log(Ptest)
    intMult = Pover*Log1
    KLDist = np.sum(intMult)

    #print 'int',np.sum(Pover)
    print KLDist
    #intSum = ((Ptest-Pover)**2)/(Pover*(1.-Pover))
    #print np.min(intSum)
    #MLSval = 1./(xlen*ylen)*np.sum(np.log10(intSum))
    #MLSval = 1./(xlen*ylen)*np.log10(np.sum(intSum))
    #return(MLSval)
    return (KLDist)
def getPDF4KLDiv(lonlat,nSamples,pdfoption,lonMesh,latMesh,R1,R2,R3,xlen,ylen):
    # lon0, lat0 : main vehicle lat and lon used as initial conditions for debris propagation
    # beta : main vehicle heading angle (counterclockwise starting East) [deg]
    # this routine should only be used for KL Divergence check. 
    # SIMPLE MODIFICATION FROM pdfCoordTrans.py
    
    xyz,x,y,z = pdf.getxyz(lonlat)

    
    
    Uint = np.dot(R2,R1)
    U  = np.dot(R3,Uint)
    #U = Uint
    UT = U.T
    pqr = np.dot(U,xyz) # rotating point to new frame...(P frame)
    
    p = pqr[0,:] 
    q = pqr[1,:]
    r = pqr[2,:]
    lonlatPframe = pdf.getlonlat(p,q,r)
    lonlatPframe = lonlatChecks.fixlon4pdf(lonlatPframe,nSamples)
    '''
    mean,covar = meancov.meancovmatrix(lonlatPframe,nSamples)
    #print 'Mean in P frame =',mean
    D,Ulocal = np.linalg.eig(covar)
    
    
    
    lonlatAframe = np.dot(lonlatPframe,Ulocal) 
    #lonlatAframe = lonlatChecks.fixlon4pdf(lonlatAframe,nSamples)
    
    mean,covar = meancov.meancovmatrix(lonlatAframe,nSamples)
  
    
    
    if pdfoption=='KDE':
        ZZpdf = meancov.kde(lonlatAframe,XX,YY,[nSamples, xlen,ylen])
    elif pdfoption=='normal':
        ZZpdf = meancov.normalbivariate(mean,covar,XX,YY,xlen,ylen) #fortran version
    '''
    
    
    
    lonVec = np.reshape(lonMesh,np.size(lonMesh))
    latVec = np.reshape(latMesh,np.size(latMesh))
    lonlatVec = np.zeros((len(lonVec),2))
    lonlatVec[:,0] = lonVec
    lonlatVec[:,1] = latVec
    lonlatVecP = pdf.originalFrame2Pframe(lonlatVec,U)
    ylen,xlen = np.shape(lonMesh)
    lonMeshP = np.reshape(lonlatVecP[:,0],[ylen,xlen])
    latMeshP = np.reshape(lonlatVecP[:,1],[ylen,xlen])
    
    
    pqrVec ,pVec,qVec,rVec= pdf.getxyz(lonlatVecP) # mesh locations in P frame
    

    if pdfoption=='KDE':
        ZZpdf = meancov.kde(lonlatPframe,lonMeshP,latMeshP,[nSamples, xlen,ylen])
    elif pdfoption=='normal' or pdfoption=='gaussian' or pdfoption=='gauss':
        mean,covar = meancov.meancovmatrix(lonlatPframe,nSamples)
        ZZpdf = meancov.normalbivariate(mean,covar,lonMeshP,latMeshP,xlen,ylen) #fortran version

    pdfN = pdf.transformPDF(U,pqrVec,np.pi/180*lonlatVec[:,0],np.pi/180*lonlatVec[:,1],np.reshape(ZZpdf,xlen*ylen))
    
    ZZpdfN = np.reshape(pdfN,[ylen,xlen])

    
 
    return (ZZpdfN)



