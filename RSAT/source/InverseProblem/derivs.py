#routine to calculate partial derivatives 
import copy
import trajectoryPropagationTools as TPT



# basic imports 

import numpy as np
# safety assessment tool path
#mainPATH = '../../../source/'
### importing safety assessment tool
# path and import for debris propagation path

### importing safety metrics and population reader
import simpleArcReader as SAR
import safetyMetrics as SF

import debrisPropagation as dp # debris propagation only tracks values at impact with planet


def EcderivsFiniteDiff(vals,params,options,step):
    import copy
    #calculating referece value
    
    valsRef = copy.deepcopy(vals)
    valsRef[11] = 0.0
    valsRef[12] = 0.0
    valsRef[13] = 0.0
    
    ec0 = runEc(valsRef,params,options)
    
    dec_dx = np.zeros((len(vals),1))
    for index in range(len(vals)):
        diffVals = copy.deepcopy(vals)
        #if (isinstance(diffVals[index], (int, long, float, complex))==False):
        #diffVals[index] = np.array(diffVals[index])
        diffVals[index] = diffVals[index] + step
        
        ecvals = runEc(diffVals,params,options)
        dec_dx[index,0] = (ecvals-ec0)/step
    return dec_dx,ec0

def EcderivsFiniteDiffNormVars(normvars,params,options,step,lowerbound,upperbound):
    import copy
    #calculating referece value
    
    vals = normVars2Vars(normvars,lowerbound,upperbound)
    
    valsRef = copy.deepcopy(vals)
    valsRef[11] = 0.0
    valsRef[12] = 0.0
    valsRef[13] = 0.0
    
    lonlat,altitudeFinal =runTraj(valsRef,params,options)
    lonlat = np.reshape(lonlat,(1,2))
    ec0 = runEc(lonlat,[altitudeFinal],params)
    
    dec_dxnorm = np.zeros((len(vals),1))
    lonlatMat = np.zeros((len(vals),2))
    altVec = np.zeros((len(vals)))
    for index in range(len(vals)):
        diffValsNorm = copy.deepcopy(normvars)
        #if (isinstance(diffVals[index], (int, long, float, complex))==False):
        #diffVals[index] = np.array(diffVals[index])
        diffValsNorm[index] = diffValsNorm[index] + step
        
        diffVals = normVars2Vars(diffValsNorm,lowerbound,upperbound)
        lonlat,altitudeFinal =runTraj(diffVals,params,options)
        lonlatMat[index,:] = lonlat
        altVec[index] = altitudeFinal
    
    
    ecvals = runEc(lonlatMat,altVec,params)
    
    dec_dxnorm = (ecvals-ec0)/step
    dec_dxnorm = np.reshape(dec_dxnorm,(len(dec_dxnorm),1))
    return dec_dxnorm,ec0

def EcderivsFiniteDiffVars(normvars,params,options,step,lowerbound,upperbound):
    import copy
    #calculating referece value
    
    vals = normVars2Vars(normvars,lowerbound,upperbound)
    
    valsRef = copy.deepcopy(vals)
    valsRef[11] = 0.0
    valsRef[12] = 0.0
    valsRef[13] = 0.0
    
    lonlat,altitudeFinal =runTraj(valsRef,params,options)
    lonlat = np.reshape(lonlat,(1,2))
    ec0 = runEc(lonlat,[altitudeFinal],params)
    
    dec_dxnorm = np.zeros((len(vals),1))
    lonlatMat = np.zeros((len(vals),2))
    altVec = np.zeros((len(vals)))
    for index in range(len(vals)):
        diffValsNorm = copy.deepcopy(normvars)
        #if (isinstance(diffVals[index], (int, long, float, complex))==False):
        #diffVals[index] = np.array(diffVals[index])
        diffValsNorm[index] = diffValsNorm[index] + step
        
        diffVals = normVars2Vars(diffValsNorm,lowerbound,upperbound)
        lonlat,altitudeFinal =runTraj(diffVals,params,options)
        lonlatMat[index,:] = lonlat
        altVec[index] = altitudeFinal
    
    
    ecvals = runEc(lonlatMat,altVec,params)
    
    dec_dxnorm = (ecvals-ec0)/step
    dec_dxnorm = np.reshape(dec_dxnorm,(len(dec_dxnorm),1))
    delta = upperbound - lowerbound
    step_reg = .5*delta*step#regular step for regular partial derivatives
    dec_dx = (ecvals-ec0)/step_reg
    dec_dx = np.reshape(dec_dx,(len(dec_dx),1))
    return dec_dx,dec_dxnorm,ec0





def EcCalc(normvars,params,options,step,lowerbound,upperbound):
    import copy
    #calculating referece value
    
    vals = normVars2Vars(normvars,lowerbound,upperbound)
    
    valsRef = copy.deepcopy(vals)
 
    
    lonlat,altitudeFinal =runTraj(valsRef,params,options)
    lonlat = np.reshape(lonlat,(1,2))
    ec0 = runEc(lonlat,[altitudeFinal],params)
    
    
    
    
    
    
 
    return ec0




def EcCalcMalFunc(normvars,params,options,optionsMalFunc,lowerbound,upperbound):
    vals = normVars2Vars(normvars,lowerbound,upperbound)
    
    tLaunchDesired = vals[0] #time of launch within launch window
    tfail = vals[1] 
    deltatfail = vals[2]
    Toffset = vals[3]
    ThrustOffsetAngDeg=[vals[4],vals[5]]
    
    missionList = optionsMalFunc[0]
    dt = options[5]
    
    paramDebris = params[0]
    atmProfile = params[1]
    
    paramsPopulation = params[2]
    filenameGE = None
    
    
    
    #obtaining the required trajectory parameters
    stateVec,thetag,mass,indexEvent,Vmag,trefFail = TPT.generateRandomStateVectorMalFunc(missionList,tLaunchDesired,dt,tfail,deltatfail,Toffset,ThrustOffsetAngDeg,atmProfile,filenameGE)
    
    diffVals = np.concatenate((stateVec,vals[6:]))
    #########    
    
    valsRef = copy.deepcopy(diffVals)
    
    
    
    
    lonlat,altitudeFinal =runTraj(valsRef,params,options,thetag)
    lonlat = np.reshape(lonlat,(1,2))
    casArea = paramDebris[5]
    ec0 = runEc(lonlat,[altitudeFinal],paramsPopulation,casArea)


    return ec0

def runEc(lonlatMat,altVec,paramsPopulation,casualtyArea=1.0):
    rp = .3048
    
    keyPop = paramsPopulation[0]
    keyArea = paramsPopulation[1]
    xllcorner = paramsPopulation[2]
    yllcorner = paramsPopulation[3]
    cellsize = paramsPopulation[4]
    xMax = paramsPopulation[5]
    yMax = paramsPopulation[6]
    ignorebounds = 1
    ecVec = np.zeros((len(altVec)))
    
    
    ecVec= SF.calculateecsinglepiece(keyPop,keyArea,xllcorner,yllcorner,cellsize,xMax,yMax,lonlatMat,ignorebounds,rp,casualtyArea)
    
    return ecVec


def runTraj(diffVals,params,option,thetag0=0.0):
    
    initialstate = diffVals[0:6]
    debrisvel = diffVals[6:9]
    massval = diffVals[9]
    sref = diffVals[10]
    cdoff = diffVals[11]
    cloff = diffVals[12]
    loverdoff = diffVals[13]
    
    cloption = option[0]
    atmosoption = option[1]
    geoptions = option[2]
    filename = option[3]
    planetmodel = option[4]
    dtinterval = option[5]
    
    
    paramsDebris = params[0]
    paramsAtm = params[1]
    
    
    minfcd = paramsDebris[0]
    cd = np.array(paramsDebris[1]) + cdoff
    minfcl = paramsDebris[2]
    cl = np.array(paramsDebris[3]) + cloff
    loverd = np.array(paramsDebris[4]) + loverdoff
    altitudelist = paramsAtm[0]
    densitylist = paramsAtm[1]
    ulist = paramsAtm[2]
    vlist = paramsAtm[3]
    wlist = paramsAtm[4]
    
    
    
    
    debrisResults= dp.debrispropagation(initialstate,debrisvel,massval,sref,minfcd,cd,cloption,
                                        minfcl,cl,loverd,atmosoption,altitudelist,
                                        densitylist,ulist,vlist,wlist,geoptions,filename,planetmodel,dtinterval,thetag0)
    
    altitudeFinal = debrisResults[0]
    latitudeFinal = debrisResults[1]
    longitudeFinal = debrisResults[2]
    Vfinal = debrisResults[3] # final velocity magnitude relative to Earth
    
    
    xylocations = [longitudeFinal,latitudeFinal]
    
    return xylocations,altitudeFinal


def normalizeVars(vars,lowerbound=None,upperbound=None):
    # function that normalizes input variables to be within [-1,1]
    print vars
    dim,varjN = np.shape(vars)
    
    if (lowerbound==None)and(upperbound==None):
        upperbound = np.max(vars,axis=1)
        lowerbound = np.min(vars,axis=1)
    elif (lowerbound==None)or(upperbound==None):
        print 'Only upper or lower bound were specified. Both or none needed'
        print 'Error in normVars in derivs.py'
        error()
    
    
    upperbound = np.reshape(upperbound,(dim,1))
    lowerbound = np.reshape(lowerbound,(dim,1))
    
    deltaBound = upperbound - lowerbound
    deltaBound[deltaBound==0.] = 1e-15 # to avoid divide by zero error...has no effect on calculation if this is the case
    newVars = 2.0*((vars-lowerbound)/deltaBound) -1.0
    return newVars,lowerbound,upperbound



def normVars2Vars(normvars,lowerbound,upperbound):
    lowerbound = np.reshape(lowerbound,(len(normvars)))
    upperbound = np.reshape(upperbound,(len(normvars)))
    deltaBound = upperbound - lowerbound
 
    vals = lowerbound + 0.5*deltaBound*(normvars+ 1.0)
    return vals

def EcderivsMalFunc(normvars,params,options,optionsMalFunc,step,lowerbound,upperbound):
    vals = normVars2Vars(normvars,lowerbound,upperbound)
    
    tLaunchDesired = vals[0] #time of launch within launch window
    tfail = vals[1] 
    deltatfail = vals[2]
    Toffset = vals[3]
    ThrustOffsetAngDeg=[vals[4],vals[5]]
    
    missionList = optionsMalFunc[0]
    dt = options[5]
    
    paramDebris = params[0]
    atmProfile = params[1]

    paramsPopulation = params[2]
    filenameGE = None

    
    
    #obtaining the required trajectory parameters
    stateVec,thetag,mass,indexEvent,Vmag,trefFail = TPT.generateRandomStateVectorMalFunc(missionList,tLaunchDesired,dt,tfail,deltatfail,Toffset,ThrustOffsetAngDeg,atmProfile,filenameGE)

    diffVals = np.concatenate((stateVec,vals[6:]))
    #########    
    
    valsRef = copy.deepcopy(diffVals)



    
    lonlat,altitudeFinal =runTraj(valsRef,params,options,thetag)
    lonlat = np.reshape(lonlat,(1,2))
    casArea = paramDebris[5]
    ec0 = runEc(lonlat,[altitudeFinal],paramsPopulation,casArea)
    if np.isnan(ec0):
        print 'ValsRef',valsRef
        print 'vals',vals
        print lonlat,altitudeFinal,casArea,ec0
        exit()
    dec_dxnorm = np.zeros((len(vals),1))
    lonlatMat = np.zeros((len(vals),2))
    altVec = np.zeros((len(vals)))
    for index in range(len(vals)):
        diffValsNorm = copy.deepcopy(normvars)
        #if (isinstance(diffVals[index], (int, long, float, complex))==False):
        #diffVals[index] = np.array(diffVals[index])
        diffValsNorm[index] = diffValsNorm[index] + step
        
        diffValsTEMP = normVars2Vars(diffValsNorm,lowerbound,upperbound)
        
        tLaunchDesired = diffValsTEMP[0] #time of launch within launch window
        tfail = diffValsTEMP[1] 
        deltatfail = diffValsTEMP[2]
        Toffset = diffValsTEMP[3]
        ThrustOffsetAngDeg=[diffValsTEMP[4],diffValsTEMP[5]]
        
        stateVec,thetag,mass,indexEvent,Vmag,trefFail = TPT.generateRandomStateVectorMalFunc(missionList,tLaunchDesired,dt,tfail,deltatfail,Toffset,ThrustOffsetAngDeg,atmProfile,filenameGE)    
        
        diffVals = np.concatenate((stateVec,diffValsTEMP[6:]))

        
        
        lonlat,altitudeFinal =runTraj(diffVals,params,options,thetag)
        lonlatMat[index,:] = lonlat
        altVec[index] = altitudeFinal

    ecvals = runEc(lonlatMat,altVec,paramsPopulation,casArea)
    print 'ec',ecvals
    dec_dxnorm = (ecvals-ec0)/step
    dec_dxnorm = np.reshape(dec_dxnorm,(len(dec_dxnorm),1))
    return dec_dxnorm,ec0

