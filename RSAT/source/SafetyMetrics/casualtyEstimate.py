# routines to calculate a gaussian distribution (from debris) in a circular planet. 


import numpy as np
import meancov
import statsPython
import lonlatChecks
import serialKernel as SK
import safetyMetrics as SMF


def getCorrelation2InRotatedFrame(beta,vals):
    # lon0, lat0 : main vehicle lat and lon used as initial conditions for debris propagation
    # beta : main vehicle heading angle (counterclockwise starting East) [deg]
    
    lonlat = vals[0]
    nSamples = vals[1]
    lon0 = vals[2]
    lat0 = vals[3]
    R1 = Rotz(lon0*np.pi/180.) # generating rotation matrices
    R2 = Roty(-lat0*np.pi/180.)
    R3 = Rotx(beta*np.pi/180.)
    
    Uint = np.dot(R2,R1)
    U  = np.dot(R3,Uint)
    UT = U.T
    xyz,x,y,z = getxyz(lonlat)
    pqr = np.dot(U,xyz) # rotating point to new frame...(P frame)
    p = pqr[0,:] 
    q = pqr[1,:]
    r = pqr[2,:]
    lonlatPframe = getlonlat(p,q,r)
    lonlatPframe = lonlatChecks.fixlon4pdf(lonlatPframe,nSamples)
    mean,covar = meancov.meancovmatrix(lonlatPframe,nSamples)

    sigma1 = (covar[0,0])**.5
    sigma2 = (covar[1,1])**.5
    #print sigma1,sigma2
    correlation = covar[0,1]/(sigma1*sigma2)
    
    return abs(correlation)*sigma2 # objective function for optimization part (minimizing rho) 

def optimizeHeadingAngle(lonlat,nSamples,lon0,lat0,beta0):
    #optimizing heading angle for PDF computation using KDE
    from scipy.optimize import fmin
    vals = [lonlat,nSamples,lon0,lat0]
    funcLambda = lambda beta : getCorrelation2InRotatedFrame(beta,vals)
    res = fmin(funcLambda,beta0,disp=0)
    return res


def pdfSetup(lonlat,pdfoption='kde',delta=None,nMesh=50):
    nSamples,dims = np.shape(lonlat)
    if dims!=2:
        print 'Error in lonlat, dimensions are not correct'
        print 'Check inputs to pdfSetup in casualtyEstimate.py'
        exit()
    '''
    #checks if KDE is necessary...no rotation. Debris always lands in a given location
    lonCheck = lonlat[:,0]
    latCheck = lonlat[:,1]
    lonSTD = np.std(lonCheck)
    latSTD = np.std(latCheck)
    '''
    
    xyz,x,y,z = getxyz(lonlat)
    xmean = np.mean(x)
    ymean = np.mean(y)
    zmean = np.mean(z)
    mag = (xmean**2+ymean**2+zmean**2)
    lonlat0 = getlonlat(xmean/mag,ymean/mag,zmean/mag)
    lon0 = lonlat0[0,0]
    lat0 = lonlat0[0,1]
    beta = 0 # initial guess
    beta = optimizeHeadingAngle(lonlat,nSamples,lon0,lat0,beta)

    R1 = Rotz(lon0*np.pi/180.) # generating rotation matrices
    R2 = Roty(-lat0*np.pi/180.)
    R3 = Rotx(beta*np.pi/180.)
    transformDetails = [R1,R2,R3,lon0,lat0,beta]
    
    Uint = np.dot(R2,R1)
    U  = np.dot(R3,Uint)
    lonlatPframe = originalFrame2Pframe(lonlat,U)
    
    
    if pdfoption.lower()=='kde' or pdfoption.lower()=='kernel':
        xMeshP,yMeshP,xlen,ylen,areaInt= statsPython.areaOfInterestKDE(lonlatPframe,nMesh = nMesh)
    #xMeshP,yMeshP,xlen,ylen,areaInt = statsPython.areaOfInterestKDE(lonlatPframe,delta) # A frame could also be used
    elif pdfoption=='normal' or pdfoption=='gaussian' or pdfoption=='gauss':
        #mean,covar = meancov.meancovmatrix(lonlatPframe,nSamples)
        #xMeshP,yMeshP,xlen,ylen = statsPython.areaOfInterest(mean,covar,delta,nsigma)
        xMeshP,yMeshP,xlen,ylen,areaInt= statsPython.areaOfInterestKDE(lonlatPframe,nMesh = 50) # A frame could also be used 
    lonlatPframeMesh = np.zeros((np.size(xMeshP),2))
    lonlatPframeMesh[:,0] = np.reshape(xMeshP,np.size(xMeshP))
    lonlatPframeMesh[:,1] = np.reshape(yMeshP,np.size(yMeshP))
    
    
    lonlatOrFrameMesh = Pframe2originalFrame(lonlatPframeMesh,U)
    lonOrMesh = np.reshape(lonlatOrFrameMesh[:,0],np.shape(xMeshP))
    latOrMesh = np.reshape(lonlatOrFrameMesh[:,1],np.shape(xMeshP))
    
    PframeVals = [xMeshP,yMeshP,lonlatPframe]
    OrFrameVals = [lonOrMesh,latOrMesh]
    rotVals = [U,transformDetails,areaInt]

    return (PframeVals,OrFrameVals,rotVals)



def getPDF(lonlat,pdfoption='kde',delta=None,nMesh=50):
    #lonlat is a nx2 matrix with lon lat locations
    # nsigma is only used for getting PDF area when using normal distribution
    # pdf option 'kde' or 'normal' dist 
    # delta is spacing used for mesh in PDF
    # if delta not specified,then it will use nMesh which is the number of points used in one dim...total mesh will be nMesh x nMesh
    
    #this routine rotates the SPHERICAL planet to get a better 'flat planet' approximation by using local lat and lons and then applies the KDE or normal PDF routines
    
    xyz,x,y,z = getxyz(lonlat)
    xmean = np.mean(x)
    ymean = np.mean(y)
    zmean = np.mean(z)
    mag = (xmean**2+ymean**2+zmean**2)**.5
    lonlat0 = getlonlat(xmean/mag,ymean/mag,zmean/mag)
    lon0 = lonlat0[0,0]
    lat0 = lonlat0[0,1]
    beta = 0 # initial guess
    nSamples,dims = np.shape(lonlat)
    if dims!=2:
        print 'Error in lonlat, dimensions are not correct'
        print 'Check inputs to getPDF in casualtyEstimate.py'
        exit()
    
    beta = optimizeHeadingAngle(lonlat,nSamples,lon0,lat0,beta)
    
    R1 = Rotz(lon0*np.pi/180.) # generating rotation matrices
    R2 = Roty(-lat0*np.pi/180.)
    R3 = Rotx(beta*np.pi/180.)
    transformDetails = [R1,R2,R3,lon0,lat0,beta]
    
    Uint = np.dot(R2,R1)
    U  = np.dot(R3,Uint)
    lonlatPframe = originalFrame2Pframe(lonlat,U)

    if pdfoption.lower()=='kde' or pdfoption.lower()=='kernel':
        xMeshP,yMeshP,xlen,ylen,areaInt = statsPython.areaOfInterestKDE(lonlatPframe,delta=delta,nMesh=nMesh) # A frame could also be used
    elif pdfoption=='normal' or pdfoption=='gaussian' or pdfoption=='gauss':
        #mean,covar = meancov.meancovmatrix(lonlatPframe,nSamples)
        #xMeshP,yMeshP,xlen,ylen = statsPython.areaOfInterest(mean,covar,delta,nsigma)
        xMeshP,yMeshP,xlen,ylen,areaInt = statsPython.areaOfInterestKDE(lonlatPframe,delta=delta,nMesh=nMesh) # A frame could also be used

    
    
    lonlatPframeMesh = np.zeros((np.size(xMeshP),2))
    lonlatPframeMesh[:,0] = np.reshape(xMeshP,np.size(xMeshP))
    lonlatPframeMesh[:,1] = np.reshape(yMeshP,np.size(yMeshP))
    lonlatOrFrameMesh = Pframe2originalFrame(lonlatPframeMesh,U)
    lonMesh,latMesh = getLonLatMesh(lonlatOrFrameMesh[:,0],lonlatOrFrameMesh[:,1],delta)
    
    lonVec = np.reshape(lonMesh,np.size(lonMesh))
    latVec = np.reshape(latMesh,np.size(latMesh))
    lonlatVec = np.zeros((len(lonVec),2))
    lonlatVec[:,0] = lonVec
    lonlatVec[:,1] = latVec
    lonlatVecP = originalFrame2Pframe(lonlatVec,U)
    ylen,xlen = np.shape(lonMesh)
    lonMeshP = np.reshape(lonlatVecP[:,0],[ylen,xlen])
    latMeshP = np.reshape(lonlatVecP[:,1],[ylen,xlen])
    
    
    pqrVec ,pVec,qVec,rVec= getxyz(lonlatVecP) # mesh locations in P frame
    
    
    if pdfoption.lower()=='kde' or pdfoption.lower()=='kernel':
        #ZZpdf = meancov.kde(lonlatPframe,lonMeshP,latMeshP,[nSamples, xlen,ylen])
        print 'Starting quad KDE'
        ZZpdf = SK.serialK(nSamples,lonlatPframe,ylen,xlen,lonMeshP,latMeshP,.1,1)
        print 'Done with quad KDE'
    elif pdfoption=='normal' or pdfoption=='gaussian' or pdfoption=='gauss':
        ZZpdf = meancov.normalbivariate(mean,covar,lonMeshP,latMeshP,xlen,ylen) #fortran version
    
    
    pdfN = transformPDF(U,pqrVec,np.pi/180*lonlatVec[:,0],np.pi/180*lonlatVec[:,1],np.reshape(ZZpdf,xlen*ylen))
    ZZpdfN = np.reshape(pdfN,[ylen,xlen])

    PframeVals = [lonPMesh,latPMesh,ZZpdf]
    OrFrameVals = [lonMesh,latMesh,ZZpdfN]
    return (PframeVals,OrFrameVals,transformDetails,areaInt)
 



def getPDFWeighted(lonlat,delta=None,nMesh=50,weightArray=1.0):
    #this routine returns a weighted KDE formulation...useful for Ec calculation
    #each sample get weighted by its corresponded weightArray
    nSamples,dims = np.shape(lonlat)
    if dims!=2:
        print 'Error in lonlat, dimensions are not correct'
        print 'Check inputs to getPDF in casualtyEstimate.py'
        exit()
    
    xyz,x,y,z = getxyz(lonlat)
    xmean = np.mean(x)
    ymean = np.mean(y)
    zmean = np.mean(z)
    mag = (xmean**2+ymean**2+zmean**2)
    lonlat0 = getlonlat(xmean/mag,ymean/mag,zmean/mag)
    lon0 = lonlat0[0,0]
    lat0 = lonlat0[0,1]
    beta = 0 # initial guess
    beta = optimizeHeadingAngle(lonlat,nSamples,lon0,lat0,beta)
    
    R1 = Rotz(lon0*np.pi/180.) # generating rotation matrices
    R2 = Roty(-lat0*np.pi/180.)
    R3 = Rotx(beta*np.pi/180.)
    transformDetails = [R1,R2,R3,lon0,lat0,beta]
    
    Uint = np.dot(R2,R1)
    U  = np.dot(R3,Uint)
    lonlatPframe = originalFrame2Pframe(lonlat,U)
 
    xMeshP,yMeshP,xlen,ylen,areaInt = statsPython.areaOfInterestKDE(lonlatPframe,delta=delta,nMesh=nMesh)     
    
    lonlatPframeMesh = np.zeros((np.size(xMeshP),2))
    lonlatPframeMesh[:,0] = np.reshape(xMeshP,np.size(xMeshP))
    lonlatPframeMesh[:,1] = np.reshape(yMeshP,np.size(yMeshP))
    row,col = np.shape(xMeshP)
    
    ZZpdfPframe = SK.serialWeightKDE(int(nSamples),lonlatPframe,row,col,xMeshP,yMeshP,weightArray)

    lonlatOrFrameMesh = Pframe2originalFrame(lonlatPframeMesh,U)
    lonOrMesh = np.reshape(lonlatOrFrameMesh[:,0],[ylen,xlen])
    latOrMesh = np.reshape(lonlatOrFrameMesh[:,1],[ylen,xlen])
    OrFrameVals = [lonOrMesh,latOrMesh]
    PframeVals = [xMeshP,yMeshP,ZZpdfPframe]
    return (PframeVals,OrFrameVals,transformDetails,areaInt)






def getPDFfromSetup(pdfoption,lonlatPframe,xMeshP,yMeshP):
    
    row,col = np.shape(xMeshP)
    if pdfoption.lower()=='kde' or pdfoption.lower()=='kernel':
        ZZpdfPframe = SK.serialK(int(nSamples),lonlatPframe,row,col,xMeshP,yMeshP,.0,0)
    elif pdfoption=='normal' or pdfoption=='gaussian' or pdfoption=='gauss':
        nSamples,dims = np.shape(lonlatPframe)
        if dims!=2:
            print 'Error in lonlatPframe, dimensions are not correct'
            print 'Check inputs to getPDFfromSetup in casualtyEstimate.py'
            exit() 
        mean,covar = meancov.meancovmatrix(lonlatPframe,nSamples)
        ZZpdfPframe = meancov.normalbivariate(mean,covar,xMeshP,yMeshP,col,row) #fortran version    
    return ZZpdfPframe


def getWeightedPDFfromSetup(lonlatPframe,xMeshP,yMeshP,weightArray):
    #simply calculated the weighted KDE. note that this is not a proper PDF...does not integrate to 1
    row,col = np.shape(xMeshP)
    nSamples,dims = np.shape(lonlatPframe)
    if dims!=2:
        print 'Error in lonlatPframe, dimensions are not correct'
        print 'Check inputs to getWeightedPDFfromSetup in casualtyEstimate.py'
        exit() 
  
    ZZpdfPframe = SK.serialWeightKDE(nSamples,lonlatPframe,row,col,xMeshP,yMeshP,weightArray)
  
    return ZZpdfPframe



def calculateEcBlast(lon,lat,populationClass,boundsOption,rad_Exp_Blast,angleDegBlast,polygon=None,desval=0.0 ):
    # TO DO...update to use classes instead of list for population
    import safetyMetrics as SMF
    import orbitTools
    #km2_to_m_2 = (1000.0)**2.0
    m2_to_km2 = 1/(1000.0)**2.0 
    #print lon,lat,angleDegBlast,rad_Exp_Blast
    # starting a mesh for blast overpressure
    lonVec = np.linspace(lon - 1.2*angleDegBlast,lon + 1.2*angleDegBlast,50)
    latVec = np.linspace(lat - 1.2*angleDegBlast,lat + 1.2*angleDegBlast,52)
    RplanetCirc = 6371000.8 # meters spherical Earth
    dlon = ( lonVec[1] - lonVec[0]) * np.pi/180.
    dlat = (latVec[1] - latVec[0]) * np.pi/180.

    lonMesh,latMesh = np.meshgrid(lonVec,latVec)

    popdensity,xMatLoc,yMatLoc = SMF.agsgetdensity(populationClass.keyPop,populationClass.keyArea,
                                                               populationClass.xllcorner,populationClass.yllcorner,
                                                               populationClass.cellsize,lonMesh,latMesh,
                                                               populationClass.xMax,populationClass.yMax,boundsOption)  
    if polygon!=None:
        xpol = polygon[0]
        ypol = polygon[1]
        matout = SMF.checkpolygon(lonMesh,latMesh,xpol,ypol)
        popdensity = SMF.updatematpolygon(lonMesh,latMesh,xpol,ypol,popdensity,desval)

    

    Ac = (dlon*dlat)*( RplanetCirc**2) * np.cos(latMesh*np.pi/180.)  
    r0= orbitTools.latlonalt2ECEF(lat,lon,0,0)
    rX,rY,rZ = orbitTools.latlonalt2ECEF(latMesh,lonMesh,0.,0)

    distMat = ((rX-r0[0])**2 + (rY-r0[1])**2 + (rZ-r0[2])**2)**.5
  
    indexSphere =  distMat > rad_Exp_Blast
    Ac[indexSphere] = 0.0

    casualties = popdensity  * (Ac*m2_to_km2)
    #print np.sum(Ac)*m2_to_km2
    #currPopulation = np.asarray(deepcopy(populationClass.keyPop),order='F') # making sure it is Fortran ordered
    #print 'P1',populationClass.keyPop[xMatLoc,yMatLoc]
    row,col = np.shape(xMatLoc)

    # population was being updated every single time
    #for i in range(row):
    #    for j in range(col):
    #        populationClass.keyPop[xMatLoc[i,j],yMatLoc[i,j]] = max(populationClass.keyPop[xMatLoc[i,j],yMatLoc[i,j]]- casualties[i,j],0.0) # updating population counts

    #print 'P2',populationClass.keyPop[xMatLoc,yMatLoc]

    #print 'Popden',popdensity
    #print np.sum(casualties),np.sum(Ac),rad_Exp_Blast**2*np.pi
    return np.sum(casualties) , populationClass




def calculateEc(lonlat,nSamples,delta,nsigma,key,xllcorner,yllcorner,cellsize,xMax,yMax,boundsOption,ncols,nrows,nPieces,arefMean,pdfoption):
    # this version uses the rotated pdf (actual lon lat coordinates). Gets pdf in regular frame, and integrates it. Inneficient for big ranges of lat lon (memory intensive)
    tempval = .3048 #1ft to meters
    PframeVals,OrFrameVals,transformDetails,areaInt = getPDF(lonlat,pdfoption=pdfoption,delta=delta)
    [lonPMesh,latPMesh,ZZpdf] = PframeVals
    [lonMesh,latMesh,ZZpdfN] = OrFrameVals
    popMatrix = SMF.agsgetvals(key,xllcorner,yllcorner,cellsize,lonMesh,latMesh,xMax,yMax,boundsOption)
    Ec = nPieces*(SMF.calculateec(ZZpdf,areaInt,arefMean,tempval,popMatrix))
    return (Ec,lonMesh,latMesh,ZZpdf)

'''
def calculateEcMatrixSheltering(lonlat,nSamples,delta,nsigma,keyPop,keyArea,xllcorner,yllcorner,cellsize,xMax,yMax,boundsOption,ncols,nrows,nPieces,casualtyArea,pdfoption,roofFraction):
    # calculates Ec by calculating pdf in P frame, then get corresponding population density for P frame.
    #print 'Npieces IN',nPieces
    if nPieces>0:
        import safetyMetrics as SMF 
        # this version uses the rotated pdf (actual lon lat coordinates)
        popOrigArea = cellsize**2
        areapop = delta*delta
        tempval = .3048 #1ft to meters
        
               
        PframeVals,OrFrameVals,rotVals = pdfSetup(lonlat,pdfoption=pdfoption,delta=delta,nMesh=50)
        
        
        [xMeshP,yMeshP,lonlatPframe] = PframeVals 
        [lonOrMesh,latOrMesh] = OrFrameVals 
        [U,transformDetails,areaInt] = rotVals
        
        ylen,xlen = np.shape(xMeshP)
        popMatrix,xmatlocations,ymatlocations = SMF.agsgetdensity(keyPop,keyArea,xllcorner,yllcorner,cellsize,lonOrMesh,latOrMesh,xMax,yMax,boundsOption,[ncols,nrows,xlen,ylen])
        
        if np.sum(popMatrix)>0.0:
            ZZpdfPframe= getPDFfromSetup(pdfoption,lonlatPframe,xMeshP,yMeshP)
            
            Ec,Ecmatrix = SMF.calculateecmatrixshelteringpopdensity(ZZpdfPframe,areapop,casualtyArea,roofFraction,tempval,popMatrix,ylen,xlen,len(casualtyArea))
            Ec = nPieces*Ec
            Ecmatrix = nPieces*Ecmatrix
        else:
            #print 'No pdf calc',np.sum(popMatrix)
            Ec = 0.0
            #Ec2 = 0.0
            xmatlocations = [1]
            ymatlocations = [1]
            Ecmatrix = np.zeros((1,1))
    
    else:
        Ec = 0.0
        #Ec2 = 0.0
        xmatlocations = [1]
        ymatlocations = [1]
        Ecmatrix = np.zeros((1,1))
    return Ec,Ecmatrix,xmatlocations,ymatlocations#,Ec2
'''
def getEc(popMatrix,weightedPDF,area):
    if weightedPDF==None:
        return 0,[]
    else:
        km2_to_m_2 = (1000.0)**2.0
        #popMatrix = popMatrix/km2_to_m_2
        Ecmatrix = (area/km2_to_m_2) *(popMatrix*weightedPDF)
        Ec = np.sum(Ecmatrix)
        return Ec,Ecmatrix

def getEcVc(popMatrix,weightedPDF,area,weightedPDF_sq):
    if weightedPDF==None:
        return [[0.,[]],0.]
    else:
        km2_to_m_2 = (1000.0)**2.0
        #popMatrix = popMatrix/km2_to_m_2
        Ecmatrix = (area/km2_to_m_2) *(popMatrix*weightedPDF)
        Ec = np.sum(Ecmatrix)
        # computing variance terms
        popMatrix2 = (popMatrix/km2_to_m_2)**2.0

        VcmatrixTerm = popMatrix2*weightedPDF_sq
        Vc = Ec - Ec**2.0 + area*np.sum(VcmatrixTerm)
        return [[Ec,Ecmatrix],Vc]





#def calculateEcMatrixShelteringWeighted(lonlat,populationClass,weightArray,
#                                        massSum,massTotalConserved,boundsOption=1,delta = None,nMesh = None,
#                                        polygon=None,ArefList = None,massList=None,CDref = 1.0,
#                                        debrisClass = None , variance = False):
def calculateEcMatrixShelteringWeighted(lonlat,populationClass,weightArray,
                                        E_N_simulated,V_N_simulated,fractionOnGround=1.0,boundsOption=1,delta = None,nMesh = None,
                                        polygon=None,ArefList = None,massList=None,CDref = 1.0,
                                        debrisClass = None , variance = False):
    # calculates Ec by calculating pdf in P frame, then get corresponding population density for P frame.
    # E_N_simulated is the expected number of debris pieces from the debris group
    # V_N_simulated is the variance of the number of debris pieces from the debris group
    # fractionOnground is the fraction of the debris pieces that impact the ground (pieces could reach orbit)
    
    # this version uses the rotated pdf (actual lon lat coordinates)
    nSamples,otherdim = np.shape(lonlat)
    E_N_ground = fractionOnGround*E_N_simulated
    V_N_ground = (fractionOnGround**2.0) * V_N_simulated
    if nSamples==0 or np.sum(weightArray)==0:
        return 0,[],[],[],0
    ################################################################################################
    #making sure we have distributions and not just point landing in the same location
    #KDE cannot handle points landing in the same location
    lonCheck = lonlat[:,0]
    latCheck = lonlat[:,1]
    maxlon = lonCheck.max()
    minlon = lonCheck.min()
    maxlat = latCheck.max()
    minlat = latCheck.min()
    dlon = maxlon - minlon
    dlat = maxlat - minlat

    #special case....deterministic point.. KDE not necessary
    if (dlon*dlat<=10**-25):
        lonOrMesh = [lonCheck.mean()]
        latOrMesh = [latCheck.mean()]
        if populationClass.keyDen ==None:
            popMatrix,xmatlocations,ymatlocations = SMF.agsgetdensity(populationClass.keyPop,populationClass.keyArea,
                                                                      populationClass.xllcorner,populationClass.yllcorner,
                                                                      populationClass.cellsize,lonOrMesh,latOrMesh,
                                                                      populationClass.xMax,populationClass.yMax,
                                                                      boundsOption)#,[populationClass.ncols,populationClass.nrows,xlen,ylen])
        else:
            popMatrix = SMF.agsgetvals(populationClass.keyDen,populationClass.xllcorner,populationClass.yllcorner,
                                       populationClass.cellsize,lonOrMesh,latOrMesh,populationClass.xMax,
                                       populationClass.yMax,boundsOption)
            
            xmatlocations = []
            ymatlocations = []
        #popMatrix = SMF.agsgetvals(populationClass.keyDen,populationClass.xllcorner,populationClass.yllcorner,
        #                           populationClass.cellsize,lonOrMesh,latOrMesh,populationClass.xMax,
        #                           populationClass.yMax,boundsOption)
        if polygon!=None:
            xpol = polygon[0]
            ypol = polygon[1]
            matout = SMF.checkpolygon(lonOrMesh,latOrMesh,xpol,ypol)
            desval = 0.0
            popMatrix = SMF.updatematpolygon(lonOrMesh,latOrMesh,xpol,ypol,popMatrix,desval)
        km2_to_m_2 = (1000.0)**2.0   
        Ec = np.sum(weightArray)/nSamples *(popMatrix[0][0]/km2_to_m_2)*E_N_ground
        return Ec,[],[],[],0 
    ################################################################################################ 

   
    popOrigArea = populationClass.cellsize**2
    
    
    (PframeVals,OrFrameVals,rotVals) = pdfSetup(lonlat,pdfoption='kde',delta=delta,nMesh=nMesh)
    
    [xMeshP,yMeshP,lonlatPframe] = PframeVals
    [lonOrMesh,latOrMesh] = OrFrameVals
    [U,transformDetails,areaInt] = rotVals

    
    
    
    ylen,xlen = np.shape(xMeshP)
    
    if populationClass.keyDen ==None:
        popMatrix,xmatlocations,ymatlocations = SMF.agsgetdensity(populationClass.keyPop,populationClass.keyArea,
                                                                  populationClass.xllcorner,populationClass.yllcorner,
                                                                  populationClass.cellsize,lonOrMesh,latOrMesh,
                                                                  populationClass.xMax,populationClass.yMax,
                                                                  boundsOption)#,[populationClass.ncols,populationClass.nrows,xlen,ylen])
    else:
        popMatrix = SMF.agsgetvals(populationClass.keyDen,populationClass.xllcorner,populationClass.yllcorner,
                                   populationClass.cellsize,lonOrMesh,latOrMesh,populationClass.xMax,
                                   populationClass.yMax,boundsOption)
        
        xmatlocations = []
        ymatlocations = []
   
    if polygon!=None:
        xpol = polygon[0]
        ypol = polygon[1]
        matout = SMF.checkpolygon(lonOrMesh,latOrMesh,xpol,ypol)
        desval = 0.0
        popMatrix = SMF.updatematpolygon(lonOrMesh,latOrMesh,xpol,ypol,popMatrix,desval)
        
   
    if np.sum(popMatrix)>0.0:
        

        
        # now divide into ballistic coefficient groups
        nAref,nmass,nlonlat,nweightArray =groupBallisticCoeff(ArefList=ArefList,massList=massList,lonlat=lonlat,weight=weightArray,CDref=CDref)
        nSamplesTotal  = 0.0
        ZZpdfPframe = 0.0
        ZZpdfPframe_sq = 0.0
        Vc = 0.0
        for index in range(len(nAref)):
            lonlatPframe = originalFrame2Pframe(nlonlat[index],U)
            weightArrayLocal = nweightArray[index]
            nSamplesloc = len(weightArrayLocal)
            nSamplesTotal = nSamplesloc + nSamplesTotal
            if nSamplesloc<=10:
                print 'Error: Adjust debris catalog. Current samples for this subgroup (rearranged by ballistic coeff) <10 samples '
                print 'Try subdiviging this debris group or add more samples'
                print 'Check debris catalog',debrisClass.name
                exit(2)

            ZZpdfPframetemp = getWeightedPDFfromSetup(lonlatPframe,xMeshP,yMeshP,weightArrayLocal)
            percentPDF = float(nSamplesloc)/float(nSamples)
            assert (percentPDF<=1.0)
            ZZpdfPframe = percentPDF*ZZpdfPframetemp + ZZpdfPframe
          
            ZZpdfPframetemp_sq = getWeightedPDFfromSetup(lonlatPframe,xMeshP,yMeshP,weightArrayLocal**2.0)
            ZZpdfPframe_sq = percentPDF*ZZpdfPframetemp_sq + ZZpdfPframe_sq

        
       
        [[Ecsingle,Ecmatrixsingle],Vcsingle] = getEcVc(popMatrix,ZZpdfPframe,areaInt,ZZpdfPframe_sq)
        Ec = Ecsingle*E_N_ground
        Vc = (Vcsingle*E_N_ground) + (V_N_ground * Ecsingle**2.0) # variance computation
        Ecmatrix = Ecmatrixsingle*E_N_ground
    
      
    else:
        print 'No pdf calc',np.sum(popMatrix)
        Ec = 0.0
        Vc = 0.0
        #Ec2 = 0.0
        xmatlocations = [1]
        ymatlocations = [1]
        Ecmatrix = np.zeros((1,1))


    return Ec,Ecmatrix,xmatlocations,ymatlocations,Vc

'''

def calculateWeightedKDE(lonlat,nSamplesModeled,delta,weightArray,massSum,massTotal,ArefList = None,massList=None,CDref = 1.0,debrisClass = None):
    # calculates Ec by calculating pdf in P frame, then get corresponding population density for P frame.
    #massSum is the added mass for all sampled debris pieces
    #massTotal is the actual physical mass for the debris group...constrained by the vehicle size
    # nSamplesModeled is not necessary equal to the number of samples in lonlat..because
    #some pieces might be in orbit, or would not cause a casualty. It is needed to calculate the weights in the KDE
    
    import safetyMetrics as SMF 
    # this version uses the rotated pdf (actual lon lat coordinates)
    nSamples,otherdim = np.shape(lonlat)
    
    if nSamples==0:
        return 0,[],[],0
    
    areapop = delta*delta
    
    mc = massSum/float(nSamplesModeled)
    kt = massTotal/mc
    
    
    
    
    print 'setting kde'
    
    
    lonlatPframe,lonOrMesh,latOrMesh,U,xMeshP,yMeshP,transformDetails = pdfSetup(lonlat=lonlat,nSamples=nSamples,
                                                                                 delta=delta,pdfoption='kde')
    ylen,xlen = np.shape(xMeshP)
    
    # now divide into ballistic coefficient groups
    nAref,nmass,nlonlat,nweightArray =groupBallisticCoeff(ArefList=ArefList,massList=massList,lonlat=lonlat,weight=weightArray,CDref=CDref)
    nSamplesTotal  = 0.0
    ZZpdfPframe = 0.0
    for index in range(len(nAref)):
        lonlatPframe = originalFrame2Pframe(nlonlat[index],U)
        weightArray = nweightArray[index]
        nSamplesloc = len(weightArray)
        nSamplesTotal = nSamplesloc + nSamplesTotal
        #print 'nSamples in minigroup',nSamplesloc
        if nSamplesloc<=10:
            print 'Error: Adjust debris catalog. Current samples for this subgroup (rearranged by ballistic coeff) <10 samples '
            print 'Try subdiviging this debris group or add more samples'
            print 'Check debris catalog',debrisClass.name
            exit(2)
        print 'running kde'
        #print lonlatPframe
        lonPMesh,latPMesh,ZZpdfPframetemp,lonOrMesh,latOrMesh = getWeightedPDFfromSetup(nSamplesloc,lonlatPframe,lonOrMesh,latOrMesh,U,xMeshP,yMeshP,weightArray)
        ZZpdfPframe = float(nSamplesloc)*ZZpdfPframetemp + ZZpdfPframe
    
    ZZpdfPframe = (1./float(nSamplesTotal))*ZZpdfPframe
    

    return lonOrMesh,latOrMesh,ZZpdfPframe,kt

def getEcfromWeightedKDE(lonOrMesh=None,latOrMesh=None,ZZpdfPframe=None,delta=None,
                         kt=None,populationClass=None,polygonClass=None,desval=0.0,boundsOption=1):
    areapop = delta**2.0
    
    if populationClass.keyDen ==None:
        popMatrix,xmatlocations,ymatlocations = SMF.agsgetdensity(populationClass.keyPop,populationClass.keyArea,
                                                                  populationClass.xllcorner,populationClass.yllcorner,
                                                                  populationClass.cellsize,lonOrMesh,latOrMesh,
                                                                  populationClass.xMax,populationClass.yMax,
                                                                  boundsOption,[populationClass.ncols,populationClass.nrows,xlen,ylen])
    else:
        popMatrix = SMF.agsgetvals(populationClass.keyDen,populationClass.xllcorner,populationClass.yllcorner,
                                   populationClass.cellsize,lonOrMesh,latOrMesh,populationClass.xMax,
                                   populationClass.yMax,boundsOption)
    
    
    if polygonClass!=None:
        
        xpol = polygonClass.xy[0]
        ypol = polygonClass.xy[1]
        matout = SMF.checkpolygon(lonOrMesh,latOrMesh,xpol,ypol)
        #print matout
        popMatrix = SMF.updatematpolygon(lonOrMesh,latOrMesh,xpol,ypol,popMatrix,desval)

    Ec,Ecmatrix = getEc(popMatrix,ZZpdfPframe,areapop)
    Ec = Ec*kt
    Ecmatrix = Ecmatrix*kt
    return Ec,Ecmatrix
'''


def getxyz(lonlat):
    deg2rad = np.pi/180.
    lon = deg2rad*lonlat[:,0]
    lat = deg2rad*lonlat[:,1]
    #x = np.zeros(nSamples)
    #y = np.zeros(nSamples)
    #z = np.zeros(nSamples)
    
    
    x = np.cos(lat)*np.cos(lon)
    y = np.cos(lat)*np.sin(lon)
    z = np.sin(lat)
    xyz = np.zeros((3,len(x)))
    xyz[0,:] = x
    xyz[1,:] = y
    xyz[2,:] = z
    
    return(xyz,x,y,z)
def singlexyz(lonlat):
    
    deg2rad = np.pi/180.
    lon = deg2rad*lonlat[0]
    lat = deg2rad*lonlat[1]
    
    
    x = np.cos(lat)*np.cos(lon)
    y = np.cos(lat)*np.sin(lon)
    z = np.sin(lat)
    
    xyz = np.array([x,y,z])
    return xyz




def getlonlat(x,y,z):
    if isinstance(x,float):
        lonlat = np.zeros((1,2))
    else:
        lonlat = np.zeros((len(x),2))
    lonlat[:,0] = 180./np.pi*np.arctan2(y,x)
    lonlat[:,1] = 180./np.pi*np.arctan2(z,np.sqrt(x**2+y**2))
    
    return lonlat

def Rotz(theta):
    R = np.array([[np.cos(theta),np.sin(theta),0],
                  [-np.sin(theta),np.cos(theta),0],
                  [0,0,1]])
    
    return R



def Roty(theta):
    R = np.array([[np.cos(theta),0,-np.sin(theta)],
                  [0,1,0],
                  [np.sin(theta),0,np.cos(theta)]])
    
    return R


def Rotx(theta):
    R = np.array([[1.,0,0],
                  [0.,np.cos(theta),np.sin(theta)],
                  [0.,-np.sin(theta),np.cos(theta)]])
    
    return R


def transformPDF(ROT,pqrVec,uVec,vVec,pdfVec):
    # form needed is fuv(u,v) = fxy(h1(x,y),h2(x,y))*abs(det(Jacobian(u,v))), det(Jacobian(u,v)) = (det(Jacobian(x,y)))^-1
    # routine to transform pdf in actual lat lon
    
    pdfVecN = np.zeros((len(pdfVec),1))
    for index in range(0,len(pdfVec)):
        pqr = pqrVec[:,index]
        detJXY = getdetJacobianXY(ROT,pqr,uVec[index],vVec[index])
        pdfVecN[index,0] = pdfVec[index]*abs(detJXY)
    return (pdfVecN)




# this section just for derivative calculation (Jacobian) to translate pdf to new frame


def getdetJacobianXY(ROT,pqr,u,v):
    # form needed is fuv(u,v) = fxy(h1(x,y),h2(x,y))*abs(det(Jacobian(u,v))), det(Jacobian(u,v)) = (det(Jacobian(x,y)))^-1
    # Derivative deriavation in Francisco Capristan's research notes
    dxdpqr,dydqpr = getdxdydqpr(pqr)
    dpqrdu,dpqrdv = getdpqrduv(ROT,u,v)
    dxdu = dxdpqr[1]*dpqrdu[1,0] + dxdpqr[0]*dpqrdu[0,0]
    dxdv = dxdpqr[1]*dpqrdv[1,0] + dxdpqr[0]*dpqrdv[0,0]
    dydu = dydqpr[0]*dpqrdu[0,0] + dydqpr[1]*dpqrdu[1,0] + dydqpr[2]*dpqrdu[2,0]
    dydv = dydqpr[0]*dpqrdv[0,0] + dydqpr[1]*dpqrdv[1,0] + dydqpr[2]*dpqrdv[2,0]
    detJXY = dxdu*dydv-dxdv*dydu
    #print 'tm',dxdpqr[1]*dpqrdv[1,0],dxdpqr[0]*dpqrdv[0,0]
    #print dxdu,dxdv,dydu,dydv
    return (detJXY)

def getdxdydqpr(pqr):
    pqr1 = pqr[0]
    pqr2 = pqr[1]
    pqr3 = pqr[2]
    pqrSQ = pqr1**2+pqr2**2
    dxdpqr1 = -pqr2/(pqrSQ)
    dxdpqr2 = pqr1/pqrSQ
    dydqpr1 = -pqr1*pqr3/(np.sqrt(pqrSQ)*(pqr1**2+pqr2**2+pqr3**2))
    dydqpr2 = -pqr2*pqr3/(np.sqrt(pqrSQ)*(pqr1**2+pqr2**2+pqr3**2))
    dydqpr3 = np.sqrt(pqrSQ)/(pqr3**2 + pqrSQ)
    
    dxdpqr = [dxdpqr1,dxdpqr2]
    dydqpr = [dydqpr1,dydqpr2,dydqpr3]
    return (dxdpqr,dydqpr)

def getdpqrduv(ROT,u,v):
    duvec = np.array([[-np.cos(v)*np.sin(u)],[np.cos(v)*np.cos(u)],[0.0]])
    dvvec = np.array([[-np.sin(v)*np.cos(u)],[-np.sin(v)*np.sin(u)],[np.cos(v)]])
    
    dpqrdu = np.dot(ROT,duvec)
    dpqrdv = np.dot(ROT,dvvec)
    #print dpqrdv
    return (dpqrdu,dpqrdv)

def getlonlatSimple(lon0,lat0,ROT):
    lonlat = np.array([lon0,lat0])
    xyz0=singlexyz(lonlat)
    xyz = np.dot(ROT,xyz0)
    
    lonlat2 = getlonlat(xyz[0],xyz[1],xyz[2],1)
    return lonlat2

def getNumDeriv(lon0,lat0,ROT):
    h = .0001
    lonlat0 = np.array([lon0,lat0])
    reflonlat1 = getlonlatSimple(lon0,lat0,ROT)
    lonlat1 = getlonlatSimple(lon0+h,lat0,ROT)
    dxdu = (lonlat1[0,0]-reflonlat1[0,0])/h
    dxdv = (lonlat1[0,1]-reflonlat1[0,1])/h
    lonlat2 = getlonlatSimple(lon0,lat0+h,ROT)
    dydu = (lonlat2[0,0]-reflonlat1[0,0])/h
    dydv = (lonlat2[0,1]-reflonlat1[0,1])/h
    #print dxdu,dxdv,dydu,dydv
    
    return (dxdu)

def Pframe2originalFrame(lonlatvecP,U):
    
    pqrVec ,pVec,qVec,rVec= getxyz(lonlatvecP) # mesh locations in P frame
    
    xyzVec = np.dot(U.T,pqrVec) # transforming mesh back to original frame
    lonlatVec = getlonlat(xyzVec[0,:],xyzVec[1,:],xyzVec[2:])
    return lonlatVec

def originalFrame2Pframe(lonlat,U):
    
    xyz,x,y,z = getxyz(lonlat)
    
    pqr = np.dot(U,xyz) # rotating point to new frame...(P frame)
    p = pqr[0,:] 
    q = pqr[1,:]
    r = pqr[2,:]
    lonlatPframe = getlonlat(p,q,r)

    lonlatPframe = lonlatChecks.fixlon4pdf(lonlatPframe,len(x)) # P frame should be technically in the principal axis
    
    return lonlatPframe



def getLonLatMesh(LonMesh,LatMesh,delta):
    minLon = np.min(LonMesh)
    maxLon = np.max(LonMesh)
    minLat = np.min(LatMesh)
    maxLat = np.max(LatMesh)
    
    lonvec = np.arange(minLon,maxLon,delta)
    latvec = np.arange(minLat,maxLat,delta)
    XX,YY = np.meshgrid(lonvec,latvec)
    return (XX,YY)


def groupBallisticCoeff(ArefList=None,massList=None,lonlat=None,weight=None,CDref = 1.0):
    #this function returns a group of lists with the appropriate debris group for PDF calculation
    #groups based on ballistic coefficient.... grouping suggested in from Chp 4 pg 126 Safety Design for Space Operation
    # note that it is expected that all pieces have the same CD. Also...we are looking at a ratio, so CD cancels out.
    #It is only useful for the lower bound....which is really low
    #beta is defined as = W/(CD*Aref) ...notice that weight is in SI (Newtons) and Aref in (m^2)
    deltaMax = 500000.
    g = 9.81 #gravity accel
    validBeta = 3.0 #lower bound for grouping
    
    arefOut = []
    massOut = []
    lonlatOut = []
    weightOut = []
    ballistic = g*(np.array(massList))/(CDref*np.array(ArefList))
    
    indexOutBounds = ballistic<validBeta
    if np.sum(indexOutBounds)>0:#handling cases that are below the lower bound threshold
        arefOut.append(ArefList[indexOutBounds])
        massOut.append(massList[indexOutBounds])
        lonlatOut.append(lonlat[indexOutBounds])
        weightOut.append(weight[indexOutBounds])
        indexInBounds = ballistic>=validBeta
        ballistic = ballistic[indexInBounds]
        ArefList = ArefList[indexInBounds]
        lonlat = lonlat[indexInBounds,:]
        weight = weight[indexInBounds]
    
    ballMin = ballistic.min()-.001 # added small value to make sure all debris belongs to a group when they are located in the extreme boundaries
    ballMax = ballistic.max()+.001
    log10BallMin = np.log10(ballMin)
    log10BallMax = np.log10(ballMax)
    
    
    
    deltabeta = log10BallMax-log10BallMin
    kval = deltabeta/deltaMax
    kdes = np.ceil(kval)
    deltaDes = deltabeta/kdes
    ballBins = np.linspace(log10BallMin,log10BallMax,kdes+1)
    ballBins = 10.0**(ballBins)
    
    
    
    
    for index in range(len(ballBins)-1):
        ballMinRange = ballBins[index]
        ballMaxRange = ballBins[index+1]
        location = (ballistic>=ballMinRange)==(ballistic<ballMaxRange) #returns index locations of ballistic coeff in range
        arefOut.append(ArefList[location])
        massOut.append(massList[location])
        lonlatOut.append(lonlat[location,:])
        weightOut.append(weight[location])
    return arefOut,massOut,lonlatOut,weightOut


