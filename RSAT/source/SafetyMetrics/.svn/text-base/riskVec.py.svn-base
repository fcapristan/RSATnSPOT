import numpy as np
import sys

# Importing f2py debrisPropagation routine folder
sys.path.append("../DebrisPropagation/SourceCode")
sys.path.append("../SafetyMetrics")
import debrisPropagation as dp # Debris Propagator. Fortran routine wrapped in python

import population
import safetyMetrics as SMF

def checkStateVector(stateVec,Vdeb,planetModel,dt,thetag,populationData):
    key,xllcorner,yllcorner,cellsize,xMax,yMax,ncols,nrows = populationData
    #    def checkStateVector(Vmag,gamma,beta,lat,long,alt,lonRef,latRef,mainVelRelOption,Vdebris,ballMin,ballMax,altitudeList,densitylist,key,xllcorner,yllcorner,cellsize,nrows, ncols,xMax,yMax):
    massval = 1.0
    arefval = 1.0
    minfcdval = 1.0
    CDval = 1.0
    minfclval = 1.0
    CLval = 0.0
    altitudeList = [1.0]
    densityList = 1.0
    uList = 1.0
    vList = 1.0
    wList = 1.0
    ncdval = 1
    nclval = 1
    nlist = 1
    

    finalConditions = dp.debrispropagation(initialstate = stateVec,
                                        debrisvel = Vdeb,
                                        mass=massval,sref=arefval,minfcd=minfcdval,cd=CDval,cloption=1,
                                        minfcl=minfclval,cl=CLval,loverd = 0,atmosoption=-1,altitudelist=altitudeList,
                                        densitylist=densityList,ulist=uList,vlist=vList,
                                        wlist=wList,geoptions=0,filename='none',planetmodel=planetModel,dtinterval = dt,thetag0=thetag,
                                        ncd=ncdval,
                                        ncl=nclval,
                                        nlist=len(altitudeList))
    
    
    
    

     

     
    altitudefinal = finalConditions[0]
    latitudefinal = finalConditions[1]
    longitudefinal = finalConditions[2]
    
    
    
    
    ballMin = 0.1
    cd = massval/(ballMin*arefval)
    atmosoption = 1# -1 Vacuum, 0 exp density, 1 gram density, 2 gram wind and density

    
    lowerConditions = dp.debrispropagation(initialstate = stateVec,
                                           debrisvel = Vdeb,
                                           mass=massval,sref=arefval,minfcd=minfcdval,cd=CDval,cloption=1,
                                           minfcl=minfclval,cl=CLval,loverd = 0,atmosoption=0,altitudelist=altitudeList,
                                           densitylist=densityList,ulist=uList,vlist=vList,
                                           wlist=wList,geoptions=0,filename='none',planetmodel=planetModel,dtinterval = dt,thetag0=thetag,
                                           ncd=ncdval,
                                           ncl=nclval,
                                           nlist=len(altitudeList))
    
    
    altitudeLower = lowerConditions[0]
    latitudeLower = lowerConditions[1]
    longitudeLower = lowerConditions[2]

    
    
    
    if (altitudefinal<0)and(altitudeLower<0):
#   calculating radius
        
        
        thetaLength = 20
        rLength = 20
        deltaLon = longitudefinal - longitudeLower
        deltaLat = latitudefinal - latitudeLower
        rmax = max(1.,(deltaLon**2+deltaLat**2)**.5) #+ .5 # .5 added to be conservative
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
#print rmax,val,centerLon,centerLat
# print xlon
#   print ylat
    return dangerZone
        
    



        