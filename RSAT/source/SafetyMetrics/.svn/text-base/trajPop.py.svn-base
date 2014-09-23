import numpy as np
import sys

# Importing f2py debrisPropagation routine folder
sys.path.append("../DebrisPropagation/SourceCode")
sys.path.append("../SafetyMetrics")
import debrisPropagation as DP # Debris Propagator. Fortran routine wrapped in python

import population
import safetyMetrics as SMF
def checkStateVector(Vmag,gamma,beta,lat,long,alt,lonRef,latRef,mainVelRelOption,Vdebris,ballMin,ballMax,altitudeList,densitylist,key,xllcorner,yllcorner,cellsize,nrows, ncols,xMax,yMax):

    #################################################### MAIN CALL TO DEBRIS PROPAGATOR #############################################################################
     
    #print ballMax,ballMin
     
     
    initialposition = np.array([lat,long,alt,lonRef,latRef])
    mainVel = np.array([Vmag+Vdebris,gamma,beta])
    #mainVelRelOption = 1 # 0 for main Vehicle inertial velocity, 1 for main vehicle relative to rotating planet velocity
    debrisRelVel = np.array([0,0,0])
    loverd = 0
                     
     
    ## Debris Modeling Options
    cloption = 0# 1 =>using CL value instead of LoverD for computing Lift force, 0 => use LoverD value for lift calculation
    planetModel = 0#0 for spherical, 1 for oblate
    geoptions = 0# 
    atmosoption = -1# -1 Vacuum, 0 exp density, 1 gram density, 2 gram wind and density
    cl = 0.
    aref = 1.
    mass = 1.
    
    cd = mass/(ballMin*aref)
    
    minfcd = 1.
    minfcl = 1.
    #altitudeList = 1.
    #densitylist = 1.
    #ulist = 1.
    #vlist = 1.
    #wlist = 1.
    filename = 'piece1'
    #print ballMin,ballMax
    ncd = 1
    ncl = 1
    nlist = len(densitylist)
    zerolist = np.zeros((nlist))
    ulist = zerolist
    vlist = zerolist
    wlist = zerolist
    #print nlist
    finalConditions = DP.debrispropagation(initialposition, mainVel,mainVelRelOption,
                                            debrisRelVel, mass,aref, # debris piece state vector wrt to main vehicle
                                            minfcd,cd,cloption,minfcl,cl,loverd, # aerodynamic inputs
                                            atmosoption,altitudeList, # dt for RK4, angrate for tumbling debris. Atmospheric options
                                            densitylist,ulist,vlist,wlist, # atmospheric profile
                                            geoptions,filename,planetModel,[ncd,ncl,nlist]) # geoptions & filename used for Google Earth visualiation. ncd => length of CD array
     # ncl => length of CL array. nlist=> length of atmospheric data array (altitude)
     
    filename = 'piece2'

     
    altitudefinal = finalConditions[0]
    latitudefinal = finalConditions[1]
    longitudefinal = finalConditions[2]
    flag = finalConditions[8]
    vrelmag = finalConditions[3]
    mainVel = np.array([Vmag-Vdebris,gamma,beta])
    
    
    cd = mass/(ballMax*aref)
    atmosoption = 1# -1 Vacuum, 0 exp density, 1 gram density, 2 gram wind and density

    
    lowerConditions = DP.debrispropagation(initialposition, mainVel,mainVelRelOption,
                                               debrisRelVel, mass,aref, # debris piece state vector wrt to main vehicle
                                               minfcd,cd,cloption,minfcl,cl,loverd, # aerodynamic inputs
                                               atmosoption,altitudeList, # dt for RK4, angrate for tumbling debris. Atmospheric options
                                               densitylist,ulist,vlist,wlist, # atmospheric profile
                                               geoptions,filename,planetModel,[ncd,ncl,nlist]) # geoptions & filename used for Google Earth visualiation. ncd => length of CD array
    # ncl => length of CL array. nlist=> length of atmospheric data array (altitude)
    
    
    
    altitudeLower = lowerConditions[0]
    latitudeLower = lowerConditions[1]
    longitudeLower = lowerConditions[2]
    flagLower = lowerConditions[8]
    vrelmagLower = lowerConditions[3]
    
    
    
    
    if (altitudefinal<0)and(altitudeLower<0):
#   calculating radius
        
        
        thetaLength = 20
        rLength = 20
        deltaLon = longitudefinal - longitudeLower
        deltaLat = latitudefinal - latitudeLower
        rmax = max(1.,(deltaLat**2+deltaLat**2)**.5) #+ .5 # .5 added to be conservative
        rmin = 0.0
        centerLon = .5*(longitudefinal + longitudeLower) 
        centerLat = .5*(latitudefinal + latitudeLower)
        theta = np.linspace(0,2*np.pi,thetaLength)
        r = np.linspace(rmin,rmax,rLength)
        rMesh,thetaMesh = np.meshgrid(r,theta)
        xlon = rMesh*np.cos(thetaMesh) + centerLon
        ylat = rMesh*np.sin(theta) + centerLat
        boundsOption = 1
        popMatrix = SMF.agsgetvals(key,xllcorner,yllcorner,cellsize,xlon,ylat,xMax,yMax,boundsOption,[ncols,nrows,rLength,thetaLength])
        val = np.sum(popMatrix)
        if val < 1e-16:
            dangerZone = 0
        else :
            dangerZone = 1
    elif (altitudeLower>0)and(altitudefinal>0):
        dangerZone = 0 # debris will be in orbit
    elif altitudeLower<0 :
        thetaLength = 20
        rLength = 20
        rmax = 2.5 # conservative bound for upper level case
        rmin = 0.0
        centerLon =  longitudeLower 
        centerLat = latitudeLower
        theta = np.linspace(0,2*np.pi,thetaLength)
        r = np.linspace(rmin,rmax,rLength)
        rMesh,thetaMesh = np.meshgrid(r,theta)
        xlon = rMesh*np.cos(thetaMesh) + centerLon
        ylat = rMesh*np.sin(theta) + centerLat
        boundsOption = 1
        popMatrix = SMF.agsgetvals(key,xllcorner,yllcorner,cellsize,xlon,ylat,xMax,yMax,boundsOption,[ncols,nrows,rLength,thetaLength])
        val = np.sum(popMatrix)
        if val < 1e-16:
            dangerZone = 0
        else :
            dangerZone = 1
    print rmax,val,centerLon,centerLat
# print xlon
#   print ylat
    return dangerZone
        
    



        