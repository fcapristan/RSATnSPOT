# this routine is used to generate trajectories to provide the state vectors for safety analysis
import malFunc

def generateTraj(missionList,atmProfile,nProp):
    #import sys
    #sys.path.append('../../../safetyAssessmentTool/DebrisPropagation/SourceCode')
    
    import numpy as np
    import orbitProp as op
    from scipy.io import loadmat
    from scipy.interpolate import UnivariateSpline
    from scipy import interpolate
    import AtmosProfile as AP
    import data2GoogleEarth as GE
    import copy
    #import 
    
    altitudelist,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,nlist = atmProfile
    
    planetmodel = 0 # 0 for spherical planet
    
    #missionList = loadmat(trajectoryFile)['missionList']
    
    cloption = 1
    cl = 0.
    minfcl = [1]
    loverd = 0.0
    atmosoption = 2

    geoptions = 0
    dtinterval = 1 # [sec]
    indval = 0
    
    stages = getNumberOfStages(missionList[0])
    # iterating through each optimal trajectory
    retTrajs = [[0]*stages for i in range(len(missionList)*nProp)]
    #exit()
    for index in range(len(missionList)):
        
        filename = 'debrisTrash'+str(index)
        mission = missionList[index]# getting mission parameters for current mission
        omegaE = getAngVelEarth(mission)

        for randindex in range(nProp):
            stateVectors = []
            densitylist, ulist,vlist,wlist = AP.generateAtmosProfile(densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList)
            fileNameGE = 'TrajProp'+str(indval) + '.kml'
            indval = indval + 1
            for staging in range(stages):
                # setting initial condition

                if staging ==0: # working with first stage
                    xM,yM,zM,VxM,VyM,VzM,ux,uy,uz,timeM,mc,eta,thetag = nominalTrajFromMission(mission,staging)
                    initialstate = np.array([xM[0],yM[0],zM[0],VxM[0],VyM[0],VzM[0]])

                    offset =  1.0*np.random.uniform(.95,1.05)
                
                else : 
                    xM,yM,zM,VxM,VyM,VzM,ux,uy,uz,timeM,mc,eta,thetag = nominalTrajFromMission(mission,staging)

                    thetag = thetag + omegaE*t[-1]
                    offset = 1.0*np.random.uniform(.95,1.05)
                    initialstate = np.array([x[-1],y[-1],z[-1],Vx[-1],Vy[-1],Vz[-1]])

                mf = mc[-1]
                m0,Fmax,isp,minfcd,cd,sref = vehicleParamFromMission(mission,eventIndex)
                Fmax = Fmax*offset
            
                time = time - time[0] # shifting time
                timeVec = np.linspace(time0[0],time0[-1],1000) 
            
            
                fcd = UnivariateSpline(minfcd,cd)
                fux = interpolate.interp1d(time,ux)
                fuy = interpolate.interp1d(time,uy)
                fuz = interpolate.interp1d(time,uz)
                feta = interpolate.interp1d(time,eta)
                
                uxN = fux(timeVec)
                uyN = fuy(timeVec)
                uzN = fuz(timeVec)
                etaN = feta(timeVec)
                minfcdN = np.linspace(minfcd[0],minfcd[-1],100)
                cdN = fcd(minfcdN)
                Fx = np.array([Fmax*etaN*uxN]).T
                Fy = np.array([Fmax*etaN*uyN]).T
                Fz = np.array([Fmax*etaN*uzN]).T
                Farray = np.concatenate([Fx,Fy,Fz],axis=1)
                ndtinterval = int(np.ceil(2.*timelist[-1]/dtinterval))
                ntime = len(timelist)
                finalconditions,derivVals = op.propagate(initialstate,m0,mf,sref,minfcdN,cdN,cloption,minfcl,cl,loverd,atmosoption,altitudelist,densitylist,ulist,vlist,wlist,Farray,timelist,isp,geoptions,filename,planetmodel,dtinterval,ndtinterval,thetag,1,ncd=len(minfcd),ncl=len(minfcl),ntime=ntime,nlist=len(altitudelist))
                
                x,y,z,Vx,Vy,Vz,m,t,newfinalConditions = fixPropagateResults(finalconditions)
                stateVectors.append(newfinalConditions)
                retTrajs[indval-1][staging] = newfinalConditions

    return retTrajs







def nominalTrajFromMission(mission,eventIndex):
    x = mission['solutionStage']['x'][eventIndex]
    y = mission['solutionStage']['y'][eventIndex]
    z = mission['solutionStage']['z'][eventIndex]
    Vx = mission['solutionStage']['Vx'][eventIndex]
    Vy = mission['solutionStage']['Vy'][eventIndex]
    Vz = mission['solutionStage']['Vz'][eventIndex]
    ux = mission['solutionStage']['ux'][eventIndex]
    uy = mission['solutionStage']['uy'][eventIndex]
    uz = mission['solutionStage']['uz'][eventIndex]
    time = mission['solutionStage']['t'][eventIndex]
    m = mission['solutionStage']['m'][eventIndex]
    eta = mission['solutionStage']['eta'][eventIndex]
    thetag = mission['solutionStage']['thetag'] 
    return (x,y,z,Vx,Vy,Vz,ux,uy,uz,time,m,eta,thetag)

def thetagFromMission(mission):
    return mission['solutionStage']['thetag'] 

def vehicleParamFromMission(mission,eventIndex):
    m0 = mission['vehicle']['mass']['m0'][eventIndex]
    Fmax = mission['vehicle']['propulsion']['F'][eventIndex]
    isp = mission['vehicle']['propulsion']['ISP'][eventIndex]
    minfcd = mission['vehicle']['aero']['MinfCD'][eventIndex]
    cd = mission['vehicle']['aero']['CD'][eventIndex]
    sref = mission['vehicle']['aero']['Aref'][eventIndex]

    return (m0,Fmax,isp,minfcd,cd,sref)

def latlonaltFromMission(mission):
    latitude = mission['launchSite']['lat']
    longitude = mission['launchSite']['long']
    height = mission['launchSite']['h']

    return (latitude,longitude,height)

def getAngVelEarth(mission):
    return mission['planet']['omega']

def getNumberOfStages(mission):
    return mission['vehicle']['stages']

def getEtaMinFromMission(mission,indexEvent):
    return mission['vehicle']['propulsion']['minThrottle'][indexEvent]
    
def getTimeFinalFromMission(mission,indexEvent):

    return mission['solution']['tf'][indexEvent]
def getFinalMassFromMission(mission,indexEvent):
    massf = mission['vehicle']['mass']['mf'][indexEvent]
    return massf


def getTimeVectorFromMission(mission,indexEvent):
    #print  len(mission['solutionStage']['t'])
    #print mission['solutionStage']['t'][indexEvent]
    #print indexEvent

    return mission['solutionStage']['t'][indexEvent]
    
def getRfromMission(mission):    
    return mission['planet']['R']



def fixPropagateResults(finalconditions):

    newIndex = finalconditions[:,10]>=0.0
    newfinalConditions = finalconditions[newIndex,:]
    # updating parameters for propagation in next stage
    x = newfinalConditions[:,0]
    y = newfinalConditions[:,1]
    z = newfinalConditions[:,2]
    Vx = newfinalConditions[:,3]
    Vy = newfinalConditions[:,4]
    Vz = newfinalConditions[:,5]      
    m = newfinalConditions[:,9]
    t = newfinalConditions[:,10] 

    return (x,y,z,Vx,Vy,Vz,m,t,newfinalConditions)

def fixPropagateResultsGetTime(finalconditions):
    
    newIndex = finalconditions[:,10]>=0.0
    newfinalConditions = finalconditions[newIndex,:]
    # updating parameters for propagation in next stage

    t = newfinalConditions[:,10] 
    
    return t

def pickTraj(trajectory,stage,indTraj,indTime):
    # trajectory => trajectory list
    # stage
    # indTime => index corresponding to the time since beginning of current stage
    # indTraj => index of trajectory to look at
    
    stateVector = trajectory[indTraj][stage][indTime,0:6]
    time = trajectory[indTraj][stage][indTime,10]
    thetag = trajectory[indTraj][stage][indTime,11]
    return stateVector,time,thetag
#generateTraj('missionList.mat')

def pickOptTraj(mission,eventIndex,indTraj,dt):
    
    import numpy as np
    from scipy.interpolate import UnivariateSpline
    x,y,z,Vx,Vy,Vz,ux,uy,uz,time,m,eta,thetag = nominalTrajFromMission(mission,eventIndex) # obtaining nominal trajectory
    timeVec = np.arange(time[0],time[-1],dt)
    fx = UnivariateSpline(time,x,s=0)
    fy = UnivariateSpline(time,y,s=0)
    fz = UnivariateSpline(time,z,s=0)
    fVx = UnivariateSpline(time,Vx,s=0)
    fVy = UnivariateSpline(time,Vy,s=0)
    fVz = UnivariateSpline(time,Vz,s=0)
    
    
    #print 'time',time
    if (len(eta)==1 and eta[0]== 0):
        ux = np.zeros((len(time)))
        uy = np.zeros((len(time)))
        uz = np.zeros((len(time)))
        eta = np.zeros((len(time)))
        m = np.zeros((len(time)))
    
    #print 'eta',m
    feta = UnivariateSpline(time,eta,s=0)
    
    fm = UnivariateSpline(time,m,s=0)
    
    fux = UnivariateSpline(time,ux,s=0)
    fuy = UnivariateSpline(time,uy,s=0)
    fuz = UnivariateSpline(time,uz,s=0)
    
    xn = np.array([fx(timeVec)]).T
    yn = np.array([fy(timeVec)]).T
    zn = np.array([fz(timeVec)]).T
    Vxn =np.array([fVx(timeVec)]).T
    Vyn = np.array([fVy(timeVec)]).T
    Vzn = np.array([fVz(timeVec)]).T
    uxn = np.array([fux(timeVec)]).T
    uyn = np.array([fuy(timeVec)]).T
    uzn = np.array([fuz(timeVec)]).T
    etan = np.array([feta(timeVec)]).T
    mn = np.array([fm(timeVec)]).T
    stateList = np.concatenate((xn,yn,zn,Vxn,Vyn,Vzn),1)
    controlList = np.concatenate((uxn,uyn,uzn,etan),1)
    #print controlList
    #print np.shape(controlList)
    return (stateList,controlList,m,timeVec,thetag)



def getALLTraj(missionFileName,dt):
    from scipy.io import loadmat
    import copy
    #trajectoryList is expected to have a set of trajectories from SPOT in increasing time of launch
    # trajectoryList[0] should have a trajectory for the beginning of launch window...trajectoryList[-1] should have the end of launch window
    missionList = loadmat(missionFileName)['missionList']
    mission0 = missionList[0]# getting mission parameters for current mission
    ntraj = len(missionList)   # getting number of trajectories
    omegaE = getAngVelEarth(mission0) # planet's angular velocity
    retList = []
    events = getNumberOfStages(mission0)
    for indTraj in range(ntraj):
        mission = missionList[indTraj]# getting mission parameters for current mission
        stateList = []
        timeList = []
        mList = []
        controlList = []
        
        for eventIndex in range(events):
            stateVector,controlVector,m,timeVec,thetag=pickOptTraj(mission,eventIndex,indTraj,dt)
            if eventIndex==0:
                thetaLaunch = thetag# getting Launch thetag
                if indTraj==0: # beginning of launch window...reference time   
                    thetaLaunchRef = copy.deepcopy(thetag)
                tLaunch = (thetaLaunch - thetaLaunchRef)/omegaE
            stateList.append(stateVector)
            timeList.append(timeVec)
            controlList.append(controlVector)
        
        retList.append([stateList,controlList,mList,timeList,thetag,tLaunch])
    return retList



def makeFunctions(retList,tdes):
    
    # tdes includes the time from launch... this value is simply used to make a 2D interpolation a 1D interpolation with a plane at tdes
    from scipy.interpolate import UnivariateSpline
    import numpy as np
    epsilon = 1e-9
    uxList = []
    uyList = []
    uzList = []
    etaList = []
    timeList = []  
    for indexTraj in range(len(retList)):
        
        [stateVec,controlVec,mVec,timeVec,thetag,timeLaunch] = retList[indexTraj]
        
        timeList = np.concatenate((timeList,[timeLaunch]),1)
        for indexEvent in range(len(stateVec)):
            timeLocal = timeVec[indexEvent]
            for indexTime in range(len(timeLocal)):
                if abs(timeLocal[indexTime]-tdes)<=epsilon:
                    uxn,uyn,uzn,etan= controlVec[indexEvent][indexTime,:]
                    uxList.append(uxn)
                    uyList.append(uyn)
                    uzList.append(uzn)
                    etaList.append(etan)  

    fux = UnivariateSpline(timeList,uxList,s=0)
    fuy = UnivariateSpline(timeList,uxList,s=0)
    fuz = UnivariateSpline(timeList,uxList,s=0)
    feta = UnivariateSpline(timeList,uxList,s=0)
    return (fux,fuy,fuz,feta,indexEvent)




def getTFevent(missionList,eventIndex):
    import numpy as np
    nevents = getNumberOfStages(missionList[0]) #nummber of stages or events

    tfList = np.zeros((nevents))

    for indexEvent in range(nevents): 
        mission = missionList[indexEvent]# getting mission parameters for current mission
        tfList[indexEvent] =getTimeFinalFromMission(mission,indexEvent)
    return tfList

def getTLaunch(missionList):
    import numpy as np
    import copy
    mission0 = missionList[0]# getting mission parameters for current mission
    ntraj = len(missionList)   # getting number of trajectories
    omegaE = getAngVelEarth(mission0)# planet's angular velocity
    events = getNumberOfStages(mission0)
    tLaunchArray = np.zeros((ntraj))-999999.9
    thetagArray = np.zeros((ntraj))-999999.9
    for indTraj in range(ntraj):
        mission = missionList[indTraj]# getting mission parameters for current mission
        thetaLaunch = thetagFromMission(mission)
        if indTraj==0: # beginning of launch window...reference time   
            thetaLaunchRef = copy.deepcopy(thetaLaunch)
        tLaunch = (thetaLaunch - thetaLaunchRef)/omegaE
        tLaunchArray[indTraj] = tLaunch
        thetagArray[indTraj] = thetaLaunch
    return tLaunchArray,thetagArray

def interpolateTrajectories(missionList,tlaunchDesired,dt):
    # getting necessary parameters for trajectory generation for a given launch window  
    import numpy as np
    from scipy import interpolate
    #from mpl_toolkits.mplot3d import Axes3D
    #from matplotlib import cm
    #from matplotlib.ticker import LinearLocator, FormatStrFormatter
    #import matplotlib.pyplot as plt
    from scipy.interpolate import griddata
    retList = []
    mission0 = missionList[0]# getting mission parameters. Mainly used to obtain constant values throughout all trajectories
    nevents = getNumberOfStages(mission0) #nummber of stages or events
    
    tLaunchVector,thetagVector = getTLaunch(missionList) # getting launch time from beginning of launch window
    
    tdeltaArray = np.zeros((len(missionList)))
    tMaxEvent = getTFevent(missionList,2)
    
    fthetag = interpolate.interp1d(tLaunchVector,thetagVector,kind='linear')
    thetag = fthetag(tlaunchDesired) # thetag required ECI frame parameters


    retList = [[]]*nevents
    
    for indexEvent in range(nevents):
        etaMin = getEtaMinFromMission(mission0,indexEvent)
    	time0 = getTimeVectorFromMission(mission0,indexEvent)
        timeMax0 = time0[-1]
        timeMin0 = time0[0]
    	timeNonDimVector = (time0-timeMin0)/(timeMax0-timeMin0) # nondimensional for all cases (between 0 and 1)....spacing should be identical between stages...needs to be calculated only once
    	timeNonDimMatrix,tLaunchMatrix = np.meshgrid(timeNonDimVector,tLaunchVector) # should be the same for all events in multiple trajectories
    	row,col = np.shape(timeNonDimMatrix)
        # initializing Matrices for fitting
    	uxMatrix = np.zeros((row,col)) - 9999999.9
    	uyMatrix = np.zeros((row,col)) - 9999999.9
    	uzMatrix = np.zeros((row,col)) - 9999999.9
    	etaMatrix = np.zeros((row,col)) - 9999999.9
        tfArray = np.zeros((len(missionList)))
        
    	for indexTraj in range(len(missionList)):
            mission = missionList[indexTraj]# getting mission parameters for current mission
            
            xM,yM,zM,VxM,VyM,VzM,ux,uy,uz,time,mc,eta,thetagc = nominalTrajFromMission(mission,indexEvent)
    	    uxMatrix[indexTraj,:] = ux 
    	    uyMatrix[indexTraj,:] = uy 
    	    uzMatrix[indexTraj,:] = uz 
    	    etaMatrix[indexTraj,:] = eta 
            tMax = time[-1]
            tMin = time[0]
            tfArray[indexTraj] = tMax - tMin # each time for each stage/event starts at zero
        
        # fit for tlaunch and t ...1D interpolation
        ftf = interpolate.interp1d(tLaunchVector,tfArray)
        tffit = ftf(tlaunchDesired)
        
        if tffit>dt and etaMin>=0:
            dtNonDim = dt/tffit
            tfitNonDim = np.arange(0.,1.,dtNonDim)
    	    tLaunchTemp = tlaunchDesired*np.ones((len(tfitNonDim)))
            
            
            timeNonDimMatrixR = np.reshape(timeNonDimMatrix,[row*col,1]) #reshaping for griddata
    	    tLaunchMatrixR = np.reshape(tLaunchMatrix,[row*col,1])
            points = np.concatenate((timeNonDimMatrixR,tLaunchMatrixR),1)
    	    uxMatrixR = np.reshape(uxMatrix,row*col)
    	    uyMatrixR = np.reshape(uyMatrix,row*col)
    	    uzMatrixR = np.reshape(uzMatrix,row*col)
    	    etaMatrixR = np.reshape(etaMatrix,row*col)
            
            uxFit = griddata(points,uxMatrixR,(tfitNonDim,tLaunchTemp),method='linear')
            uyFit = griddata(points,uyMatrixR,(tfitNonDim,tLaunchTemp),method='linear')
            uzFit = griddata(points,uzMatrixR,(tfitNonDim,tLaunchTemp),method='linear')
            etaFit = griddata(points,etaMatrixR,(tfitNonDim,tLaunchTemp),method='linear')
            
            tfit = tffit*tfitNonDim # actual time for t [sec]
            retList[indexEvent] = [indexEvent,tfit,uxFit,uyFit,uzFit,etaFit]
        

        elif etaMin<0. :
            
            retList[indexEvent] = [indexEvent,tffit]

    return (retList,thetag)
    
    
def getMatricesTraj(missionList,tlaunchDesired,dt):
    # getting necessary parameters for trajectory generation for a given launch window...this funcition populates matrices 
    import numpy as np
    
    retList = []
    mission0 = missionList[0]# getting mission parameters. Mainly used to obtain constant values throughout all trajectories
    nevents = mission0['vehicle']['stages'] #nummber of stages or events
    
    tLaunchVector,thetagVector = getTLaunch(missionList) # getting launch time from beginning of launch window
    
    
    retList = [[]]*nevents
    
    for indexEvent in range(nevents):
        etaMin = getEtaMinFromMission(mission,indexEvent)
    	time0 = getTimeVectorFromMission(mission0,indexEvent)
        timeMax0 = time0[-1]
        timeMin0 = time0[0]
    	timeNonDimVector = (time0-timeMin0)/(timeMax0-timeMin0) # nondimensional for all cases (between 0 and 1)....spacing should be identical between stages...needs to be calculated only once
    	timeNonDimMatrix,tLaunchMatrix = np.meshgrid(timeNonDimVector,tLaunchVector) # should be the same for all events in multiple trajectories
    	row,col = np.shape(timeNonDimMatrix)
        # initializing Matrices for fitting
    	uxMatrix = np.zeros((row,col)) - 9999999.9
    	uyMatrix = np.zeros((row,col)) - 9999999.9
    	uzMatrix = np.zeros((row,col)) - 9999999.9
    	etaMatrix = np.zeros((row,col)) - 9999999.9
        tfArray = np.zeros((len(missionList)))
        
    	for indexTraj in range(len(missionList)):
            mission = missionList[indexTraj]# getting mission parameters for current mission
                    
            
            xM,yM,zM,VxM,VyM,VzM,ux,uy,uz,time,mc,eta,thetag = nominalTrajFromMission(mission,staging)
    	    uxMatrix[indexTraj,:] = ux 
    	    uyMatrix[indexTraj,:] = uy 
    	    uzMatrix[indexTraj,:] = uz 
    	    etaMatrix[indexTraj,:] = eta 

            tMax = time[-1]
            tMin = time[0]
            tfArray[indexTraj] = tMax - tMin # each time for each stage/event starts at zero
        
        
        if  etaMin>=0:

            timeNonDimMatrixR = np.reshape(timeNonDimMatrix,[row*col,1]) #reshaping for griddata
    	    tLaunchMatrixR = np.reshape(tLaunchMatrix,[row*col,1])
            points = np.concatenate((timeNonDimMatrixR,tLaunchMatrixR),1)
    	    uxMatrixR = np.reshape(uxMatrix,row*col)
    	    uyMatrixR = np.reshape(uyMatrix,row*col)
    	    uzMatrixR = np.reshape(uzMatrix,row*col)
    	    etaMatrixR = np.reshape(etaMatrix,row*col)
            
            
            retList[indexEvent] = [indexEvent,points,uxMatrixR,uyMatrixR,uzMatrixR,etaMatrixR,tfArray]
        
        
        elif etaMin<0. :
            
            retList[indexEvent] = [indexEvent,tfArray]
    
    return [retList,tLaunchVector,thetagVector]

def interpolateTrajectoriesFromMatrices(trajMatrixList,tlaunchDesired,dt):
    # getting necessary parameters for trajectory generation for a given launch window  
    import numpy as np
    from scipy import interpolate
    from scipy.interpolate import griddata
    
    trajList,tLaunchVector,thetagVector = trajMatrixList
    
    retList = []
    nevents = len(trajList)
    
    fthetag = interpolate.interp1d(tLaunchVector,thetagVector,kind='linear')
    thetag = fthetag(tlaunchDesired) # thetag required ECI frame parameters
    
    retList = [[]]*nevents
    
    for indexEvent in range(nevents):
        
        
        if len(trajList[indexEvent])>=3:
            
            [indexEventVal,points,uxMatrixR,uyMatrixR,uzMatrixR,etaMatrixR,tfArray] = trajList[indexEvent]
        else:
            [indexEvent,tfArray] = trajList[indexEvent]
        # fit for tlaunch and t ...1D interpolation
        
        ftf = interpolate.interp1d(tLaunchVector,tfArray)
        tffit = ftf(tlaunchDesired)
        if tffit>dt and len(trajList[indexEvent])>=3:
            
            
            dtNonDim = dt/tffit
            tfitNonDim = np.arange(0.,1.,dtNonDim)
    	    tLaunchTemp = tlaunchDesired*np.ones((len(tfitNonDim)))
            
            
            uxFit = griddata(points,uxMatrixR,(tfitNonDim,tLaunchTemp),method='linear')
            uyFit = griddata(points,uyMatrixR,(tfitNonDim,tLaunchTemp),method='linear')
            uzFit = griddata(points,uzMatrixR,(tfitNonDim,tLaunchTemp),method='linear')
            etaFit = griddata(points,etaMatrixR,(tfitNonDim,tLaunchTemp),method='linear')
            
            tfit = tffit*tfitNonDim # actual time for t [sec]
            retList[indexEvent] = [indexEvent,tfit,uxFit,uyFit,uzFit,etaFit]
        
        
        else :
            
            
            retList[indexEvent] = [indexEvent,tffit]
    
    
    return (retList,thetag)





def generateRandomMultiTraj(missionList,tLaunchDesired,dt,atmProfile=[],TOffset = 1.0,ThrustOffsetAngDeg=[0,0]):
    # this function generates entire trajectories
    import numpy as np
    #calculating ECI coordinates for launch 
    import orbitTools # library for coordinate transformation
    import orbitProp as op # trajectory propagation routine
    import data2GE
    import copy
    
    if len(atmProfile)>0:
        atmosoption =2
        [altitudelist,densitylist,ulist,vlist,wlist] = atmProfile
    else:
        atmosoption = 0
        densitylist = [0]
        ulist = [0]
        vlist = [0]
        wlist = [0]
        altitudelist = [0]
    
    mission0 = missionList[0]# getting mission parameters for current mission
    omegaE = getAngVelEarth(mission0)# planet's angular velocity
    retList = []
    
    latitude,longitude,height = latlonaltFromMission(mission0)

    planetmodel = 0# 0 for Circ ...1 for ellipt
    trajList,thetag0 = interpolateTrajectories(missionList,tLaunchDesired,dt)
    
    r0 = orbitTools.latlonalt2ECI(latitude,longitude,height,thetag0,planetmodel)
    Rearth = getRfromMission(mission0)
    omega = getAngVelEarth(mission0)    

    nEvents = getNumberOfStages(mission0)
    cloption = 1
    cl = 0.
    minfcl = [1]
    loverd = 0.0
    
    geoptions = 0
    
    dtinterval = dt
    
    ndtinterval = 20000
    fileList = []
    currentTime = 0.0
    retList= []
    for indexEvent in range(nEvents):
        if indexEvent ==0 : #first stage...setting initial conditions for rest of cases
            
            [vS0,vE0,vZ0] = [0,0,0] # initial velocity in SEZ frame
            inertial = 0
            thetag = copy.deepcopy(thetag0)
            # inertial = 1. input velocities are inertial velocities in SEZ frame
            # inertial =0. input velocities are local velocities (earth rotation not accounted)
            #Veci0 = orbitTools.SEZ2ECI(latitude,longitude,r0,vS0,vE0,vZ0,inertial,thetag)
            Veci0 = orbitTools.SEZ2ECI(latitude,longitude,height,vS0,vE0,vZ0,inertial,thetag,planetmodel)
            initialstate = np.array([r0[0],r0[1],r0[2],Veci0[0],Veci0[1],Veci0[2]])
        m0,Fmax,isp,minfcd,cd,sref = vehicleParamFromMission(mission0,indexEvent)
        
        massf = getFinalMassFromMission(mission0,indexEvent)# m0[-1]        
        Fmax = TOffset*Fmax
        etaMin = getEtaMinFromMission(mission0,indexEvent)
        filename = 'pTraj'+str(indexEvent)+'.txt'
        fileList.append(filename)
        if etaMin>=0:
            
            propOption = 1
            currTraj = trajList[indexEvent]
            [indexEvent,tfit,uxFit,uyFit,uzFit,etaFit] = currTraj
            # ensuring that unit vector is actually a unit vector (due to interpolation)
          
            uMag = (uxFit**2 + uyFit**2 + uzFit**2)**.5
            uxetaFit = np.array([etaFit*uxFit/uMag]).T
            uyetaFit = np.array([etaFit*uyFit/uMag]).T
            uzetaFit = np.array([etaFit*uzFit/uMag]).T
            uetaMat = np.concatenate((uxetaFit,uyetaFit,uzetaFit),1)
            Fmat = Fmax*uetaMat 
            timelist = tfit
            propCond = massf
            ntime = len(timelist)
            

            
            
            finalconditions,finalderivs = op.propagate(initialstate, m0, propCond,sref,minfcd,cd,cloption, minfcl,cl,loverd,atmosoption, altitudelist,densitylist,ulist,vlist ,wlist,Fmat,timelist,isp,geoptions,filename,planetmodel,dtinterval,ndtinterval,thetag,propOption,ThrustOffsetAngDeg,ncd=len(minfcd),ncl=len(minfcl),ntime=ntime,nlist=len(altitudelist))
            newIndex = finalconditions[:,10]>=0.0
            newfinalConditions = finalconditions[newIndex,:]
            # updating parameters for propagation in next stage
            xc = newfinalConditions[:,0]
            yc = newfinalConditions[:,1]
            zc = newfinalConditions[:,2]
            Vxc = newfinalConditions[:,3]
            Vyc = newfinalConditions[:,4]
            Vzc = newfinalConditions[:,5]      
            mc = newfinalConditions[:,9]
            tc = newfinalConditions[:,10] + currentTime
            initialstate = np.array([xc[-1],yc[-1],zc[-1],Vxc[-1],Vyc[-1],Vzc[-1]])
            mass0 = mc[-1]   
            
            currentTime = tc[-1] # updating time
            thetagVector = omegaE*tc + thetag0
            thetag = omegaE*currentTime + thetag0 # updating thetag for next event
        
        elif etaMin<0:
            
            propOption = 2
            currTraj = trajList[indexEvent]
            [indexEvent,tffit] = currTraj
            dtlocal = tffit/5.
            Tmax = 0.0
            uetaMat = np.array([[1,1,1]])
            Fmat = Fmax*uetaMat 
            timelist = [tffit]
            propCond = tffit
            ntime = len(timelist)
            
            
            
            finalconditions,finalderivs = op.propagate(initialstate, mass0, propCond,sref,minfcd,cd,cloption, minfcl,cl,loverd,atmosoption, altitudelist,densitylist,ulist,vlist ,wlist,Fmat,timelist,isp,geoptions,filename,planetmodel,dtlocal,ndtinterval,thetag,propOption,ThrustOffsetAngDeg,ncd=len(minfcd),ncl=len(minfcl),ntime=ntime,nlist=len(altitudelist))
            
            xc,yc,zc,Vxc,Vyc,Vzc,mc,tc,newfinalConditions = fixPropagateResults(finalconditions)
            
            tc = tc + currentTime


            mass0 = mc[-1]   
            currentTime = tc[-1] # updating time
            thetagVector = omegaE*tc + thetag0
            thetag = omegaE*currentTime + thetag0# updating thetag for next event
        
        retList.append([xc,yc,zc,Vxc,Vyc,Vzc,mc,tc,thetagVector,indexEvent])
    return retList
#data2GE.convertMultiple(fileList)




def generateRandomMultiTrajMalFunc(missionList,tLaunchDesired,dt,timeMalfunc,deltatfail,TOffset,ThrustOffsetAngDegMalFunc,atmProfile):
    # this function generates entire trajectories
    import numpy as np
    #calculating ECI coordinates for launch 
    import orbitTools # library for coordinate transformation
    import orbitProp as op # trajectory propagation routine
    import data2GE
    import copy
    
    
    if len(atmProfile)>0:
        atmosoption =1
        [altitudelist,densitylist,ulist,vlist,wlist] = atmProfile
    else:
        atmosoption = 0
        densitylist = [0]
        ulist = [0]
        vlist = [0]
        wlist = [0]
        altitudelist = [0]
    
    mission0 = missionList[0]# getting mission parameters for current mission
    omegaE = getAngVelEarth(mission0)# planet's angular velocity
    retList = []
    
    latitude,longitude,height = latlonaltFromMission(mission0)
    
    planetmodel = 0# 0 for Circ ...1 for ellipt
    trajList,thetag0 = interpolateTrajectories(missionList,tLaunchDesired,dt)
    
    r0 = orbitTools.latlonalt2ECI(latitude,longitude,height,thetag0,planetmodel)
    Rearth = getRfromMission(mission0)
    omega = getAngVelEarth(mission0)    
    
    nEvents = getNumberOfStages(mission0)
    cloption = 1
    cl = 0.
    minfcl = [1]
    loverd = 0.0
    
    geoptions = 0
    
    dtinterval = dt
    
    ndtinterval = 20000
    fileList = []
    currentTime = 0.0
    retList= []
    malFuncBool = False
    for indexEvent in range(nEvents):
        ThrustOffsetAngDeg= [0,0]#ThrustOffsetAngDegMalFunc

        if indexEvent ==0 : #first stage...setting initial conditions for rest of cases
            
            [vS0,vE0,vZ0] = [0,0,0] # initial velocity in SEZ frame
            inertial = 0
            thetag = copy.deepcopy(thetag0)
            # inertial = 1. input velocities are inertial velocities in SEZ frame
            # inertial =0. input velocities are local velocities (earth rotation not accounted)
            #Veci0 = orbitTools.SEZ2ECI(latitude,longitude,r0,vS0,vE0,vZ0,inertial,thetag)
            Veci0 = orbitTools.SEZ2ECI(latitude,longitude,height,vS0,vE0,vZ0,inertial,thetag,planetmodel)
            initialstate = np.array([r0[0],r0[1],r0[2],Veci0[0],Veci0[1],Veci0[2]])
        m0,Fmax,isp,minfcd,cd,sref = vehicleParamFromMission(mission0,indexEvent)
        massf = getFinalMassFromMission(mission0,indexEvent)# m0[-1]        
        Fmax = TOffset*Fmax
        etaMin = getEtaMinFromMission(mission0,indexEvent)
        filename = 'pTraj'+str(indexEvent)+'.txt'
        fileList.append(filename)
        if etaMin>=0:
            
            propOption = 1
            currTraj = trajList[indexEvent]
            [indexEvent,tfit,uxFit,uyFit,uzFit,etaFit] = currTraj
            # ensuring that unit vector is actually a unit vector (due to interpolation)
            
            uMag = (uxFit**2 + uyFit**2 + uzFit**2)**.5
            uxetaFit = np.array([etaFit*uxFit/uMag]).T
            uyetaFit = np.array([etaFit*uyFit/uMag]).T
            uzetaFit = np.array([etaFit*uzFit/uMag]).T
            uetaMat = np.concatenate((uxetaFit,uyetaFit,uzetaFit),1)
            Fmat = Fmax*uetaMat 
            timelist = tfit
            propCond = massf
            ntime = len(timelist)
            
            
            # letting vehicle propagate until propOption is reached (final mass from input parameter from SPOT!! not from propagated trajectory or desired time). 
            finalconditions,finalderivs = op.propagate(initialstate, m0, propCond,sref,minfcd,cd,cloption, minfcl,cl,loverd,atmosoption, altitudelist,densitylist,ulist,vlist ,wlist,Fmat,timelist,isp,geoptions,filename,planetmodel,dtinterval,ndtinterval,thetag,propOption,ThrustOffsetAngDeg,ncd=len(minfcd),ncl=len(minfcl),ntime=ntime,nlist=len(altitudelist))
            
            tfTemp = fixPropagateResultsGetTime(finalconditions)
        
            if (timeMalfunc - currentTime )<=tfTemp[-1]:
   
                geoptions=0
                #print 'TOff',ThrustOffsetAngDegMalFunc

                timeMalfuncLocal = timeMalfunc -currentTime
                finalconditions,finalderivs = malFunc.propagate(initialstate,m0,propCond,sref,minfcd,cd,cloption,minfcl,cl,loverd,atmosoption,altitudelist,densitylist,ulist,vlist,wlist,Fmat,timeMalfuncLocal,deltatfail,timelist,isp,geoptions,filename,planetmodel,dtinterval,ndtinterval,thetag,propOption,ThrustOffsetAngDegMalFunc)
                malFuncBool = True
                
            xc,yc,zc,Vxc,Vyc,Vzc,mc,tc,newfinalConditions = fixPropagateResults(finalconditions)
            #print 'FC',newfinalConditions[-1]

            #print 'new',newfinalConditions[-1]
            # print 'tc',tc,currentTime

            #print 'tc',tc[0],tc[-1]
            tc = tc + currentTime
            initialstate = np.array([xc[-1],yc[-1],zc[-1],Vxc[-1],Vyc[-1],Vzc[-1]])
            mass0 = mc[-1]   
            
            currentTime = tc[-1] # updating time
            #print 'Current',currentTime,tc[-1]
            thetagVector = omegaE*tc + thetag0
            thetag = omegaE*currentTime + thetag0 # updating thetag for next event



        
        elif etaMin<0:
            
            propOption = 2
            currTraj = trajList[indexEvent]
            [indexEvent,tffit] = currTraj
            dtlocal = tffit/5.
            Tmax = 0.0
            uetaMat = np.array([[1,1,1]])
            Fmat = Fmax*uetaMat 
            timelist = [tffit]
            propCond = tffit
            ntime = len(timelist)
            
            
            
            finalconditions,finalderivs = op.propagate(initialstate, mass0, propCond,sref,minfcd,cd,cloption, minfcl,cl,loverd,atmosoption, altitudelist,densitylist,ulist,vlist ,wlist,Fmat,timelist,isp,geoptions,filename,planetmodel,dtlocal,ndtinterval,thetag,propOption,ThrustOffsetAngDeg,ncd=len(minfcd),ncl=len(minfcl),ntime=ntime,nlist=len(altitudelist))
            
            xc,yc,zc,Vxc,Vyc,Vzc,mc,tc,newfinalConditions = fixPropagateResults(finalconditions)
            tc = tc + currentTime
            
            
            mass0 = mc[-1]   
            currentTime = tc[-1] # updating time
            thetagVector = omegaE*tc + thetag0
            thetag = omegaE*currentTime + thetag0# updating thetag for next event

    
        

        retList.append([xc,yc,zc,Vxc,Vyc,Vzc,mc,tc,thetagVector,indexEvent])
        if malFuncBool==True:
            break

    #print 'tt', timelist[-1],currentTime,timeMalfunc
    if (timelist[-1] + currentTime - tc[-1])< timeMalfunc:
        print 'Warning when propagating malfunction turns. Failure time seems to be greater than thrust profile time bounds' 
    return retList


def generateRandomStateVector(missionList,tLaunchDesired,tfail,dt,atmProfile=[],TOffset = 1.0,ThrustOffsetAngDeg=[0,0],fileNameGE=None):
    from scipy import interpolate
    import numpy as np
    import orbitTools
    import data2GoogleEarth as DG
    retList = generateRandomMultiTraj(missionList,tLaunchDesired,dt,atmProfile,TOffset,ThrustOffsetAngDeg)
    mission0 = missionList[0]# getting mission parameters for current mission
    omegaE = getAngVelEarth(mission0)  # planet's angular velocity
    for index in range(len(retList)):
        [xc,yc,zc,Vxc,Vyc,Vzc,mc,tc,thetagVector,indexEvent] = retList[index]
        if (index==(len(retList)-1)) and tc[-1]<tfail:
            # ensuring that values are bounded...this should not be an issue because late in the trajectory the IIP should be in orbit....no chance of hitting ground
            tfail = tc[-1]
            print 'Warning in getPropTraj.py. Failure time is greater than modeled trajectory time.'
            print 'Adjusting failure time to maximum trajectory time available. New failure time is ',tfail,' sec.'
        
        if tc[-1]>=tfail : # look in this event (stage)
  
            fx = interpolate.interp1d(tc,xc,kind='linear')
            fy = interpolate.interp1d(tc,yc,kind='linear')
            fz = interpolate.interp1d(tc,zc,kind='linear')
            fVx = interpolate.interp1d(tc,Vxc,kind='linear')
            fVy = interpolate.interp1d(tc,Vyc,kind='linear')
            fVz = interpolate.interp1d(tc,Vzc,kind='linear')
            fm = interpolate.interp1d(tc,mc,kind='linear')
            fthetag = interpolate.interp1d(tc,thetagVector,kind='linear')
            [x,y,z,Vx,Vy,Vz] = [fx(tfail),fy(tfail),fz(tfail),fVx(tfail),fVy(tfail),fVz(tfail)]
            stateVec = [x,y,z,Vx,Vy,Vz]
            mass = fm(tfail)
            thetag = fthetag(tfail)

            Vrot = np.cross([0,0,omegaE],[x,y,z])
            Vinf = np.array([Vx,Vy,Vz]) - Vrot
            Vmag = (Vinf[0]**2 + Vinf[1]**2 + Vinf[2]**2)**.5  
            trefFail = tfail - tc[0] # fail time from beginning of event
            
            if fileNameGE!=None:
                # temporary addition
                lon = []
                lat = []
                h = []
                for index2 in range(len(xc)):
                    templatlon = orbitTools.ECI2latlonalt([xc[index2],yc[index2],zc[index2]],thetagVector[index2],0)
                    lon.append(templatlon[1])
                    lat.append(templatlon[0])
                    h.append(templatlon[2])
                 
                
                DG.convertPropagate(lat,lon,h,fileNameGE)
            
            
            return (stateVec,thetag,mass,indexEvent,Vmag,trefFail)






def generateRandomStateVectorMalFunc(missionList,tLaunchDesired,dt,timeMalfunc,deltatfail,TOffset = 1.0,ThrustOffsetAngDeg=[0,0],atmProfile=[],fileNameGE=None):
    from scipy import interpolate
    import numpy as np
    import orbitTools
    import data2GoogleEarth as DG
    #print 'Toff1',ThrustOffsetAngDeg
    retList = generateRandomMultiTrajMalFunc(missionList,tLaunchDesired,dt,timeMalfunc,deltatfail,TOffset,ThrustOffsetAngDeg,atmProfile)
    mission0 = missionList[0]# getting mission parameters for current mission
    omegaE = getAngVelEarth(mission0)  # planet's angular velocity
    tfail = timeMalfunc + deltatfail
    for index in range(len(retList)):
        [xc,yc,zc,Vxc,Vyc,Vzc,mc,tc,thetagVector,indexEvent] = retList[index]
        if (index==(len(retList)-1)) and tc[-1]<tfail:
            # ensuring that values are bounded...this should not be an issue because late in the trajectory the IIP should be in orbit....no chance of hitting ground
            tfail = tc[-1]
            print 'Warning in getPropTraj.py. Failure time is greater than modeled trajectory time.'
            print 'Adjusting failure time to maximum trajectory time available. New failure time is ',tfail,' sec.'
        
        if tc[-1]>=tfail : # look in this event (stage)
            
            fx = interpolate.interp1d(tc,xc,kind='linear')
            fy = interpolate.interp1d(tc,yc,kind='linear')
            fz = interpolate.interp1d(tc,zc,kind='linear')
            fVx = interpolate.interp1d(tc,Vxc,kind='linear')
            fVy = interpolate.interp1d(tc,Vyc,kind='linear')
            fVz = interpolate.interp1d(tc,Vzc,kind='linear')
            fm = interpolate.interp1d(tc,mc,kind='linear')
            fthetag = interpolate.interp1d(tc,thetagVector,kind='linear')
            [x,y,z,Vx,Vy,Vz] = [fx(tfail),fy(tfail),fz(tfail),fVx(tfail),fVy(tfail),fVz(tfail)]
            stateVec = [x,y,z,Vx,Vy,Vz]
            mass = fm(tfail)
            thetag = fthetag(tfail)
            
            Vrot = np.cross([0,0,omegaE],[x,y,z])
            Vinf = np.array([Vx,Vy,Vz]) - Vrot
            Vmag = (Vinf[0]**2 + Vinf[1]**2 + Vinf[2]**2)**.5  
            trefFail = tfail - tc[0] # fail time from beginning of event
            
            if fileNameGE!=None:
                # temporary addition
                lon = []
                lat = []
                h = []
                for index2 in range(len(tc)):
                    templatlon = orbitTools.ECI2latlonalt([xc[index2],yc[index2],zc[index2]],thetagVector[index2],0)
                    lon.append(templatlon[1])
                    lat.append(templatlon[0])
                    h.append(templatlon[2])

                
                DG.convertPropagate(lat,lon,h,fileNameGE)
            return (stateVec,thetag,mass,indexEvent,Vmag,trefFail)


def generateRandomStateVectorMalFunc_SIMPLE(missionList,tLaunchDesired,dt,timeMalfunc,deltatfail,TOffset = 1.0,ThrustOffsetAngDeg=[0,0],atmProfile=[],fileNameGE=None):
    from scipy import interpolate
    import numpy as np
    import orbitTools
    import data2GoogleEarth as DG
    #print 'Toff1',ThrustOffsetAngDeg
    retList = generateRandomMultiTrajMalFunc(missionList,tLaunchDesired,dt,timeMalfunc,deltatfail,TOffset,ThrustOffsetAngDeg,atmProfile)
    mission0 = missionList[0]# getting mission parameters for current mission
    omegaE = getAngVelEarth(mission0)  # planet's angular velocity
    tfail = timeMalfunc + deltatfail
    for index in range(len(retList)):
        [xc,yc,zc,Vxc,Vyc,Vzc,mc,tc,thetagVector,indexEvent] = retList[-1]
        if (index==(len(retList)-1)) and tc[-1]<tfail:
            # ensuring that values are bounded...this should not be an issue because late in the trajectory the IIP should be in orbit....no chance of hitting ground
            tfail = tc[-1]
            print 'Warning in getPropTraj.py. Failure time is greater than modeled trajectory time.'
            print 'Adjusting failure time to maximum trajectory time available. New failure time is ',tfail,' sec.'
        
        if tc[-1]>=tfail : # look in this event (stage)
            
  
            [x,y,z,Vx,Vy,Vz] =   [xc[-1],yc[-1],zc[-1],Vxc[-1],Vyc[-1],Vzc[-1]] 
            stateVec = [x,y,z,Vx,Vy,Vz]
            mass = mc[-1]
            thetag = thetagVector[-1]
            
            Vrot = np.cross([0,0,omegaE],[x,y,z])
            Vinf = np.array([Vx,Vy,Vz]) - Vrot
            Vmag = (Vinf[0]**2 + Vinf[1]**2 + Vinf[2]**2)**.5  
            trefFail = tfail - tc[0] # fail time from beginning of event
            
            if fileNameGE!=None:
                # temporary addition
                lon = []
                lat = []
                h = []
                for index2 in range(len(tc)):
                    templatlon = orbitTools.ECI2latlonalt([xc[index2],yc[index2],zc[index2]],thetagVector[index2],0)
                    lon.append(templatlon[1])
                    lat.append(templatlon[0])
                    h.append(templatlon[2])
                
                
                DG.convertPropagate(lat,lon,h,fileNameGE)
            return (stateVec,thetag,mass,indexEvent,Vmag,trefFail)




def generateRandomTrajfromMatrix(missionList,trajMatrixList,tLaunchDesired,dt,atmProfile=[],TOffset = 1.0,ThrustOffsetAngDeg=[0,0]):
    # this function generates entire trajectories
    import numpy as np
    #calculating ECI coordinates for launch 
    import orbitTools # library for coordinate transformation
    import orbitProp as op # trajectory propagation routine
    import data2GE
    import copy
    
    if len(atmProfile)>0:
        atmosoption =2
        [altitudelist,densitylist,ulist,vlist,wlist] = atmProfile
    else:
        atmosoption = 0
        densitylist = [0]
        ulist = [0]
        vlist = [0]
        wlist = [0]
        altitudelist = [0]
    
    mission0 = missionList[0]# getting mission parameters for current mission
    omegaE = getAngVelEarth(mission0)# planet's angular velocity
    retList = []
    
    latitude,longitude,height = latlonaltFromMission(mission0)

            
    
    planetmodel = 0# 0 for Circ ...1 for ellipt
    #trajList,thetag0 = interpolateTrajectories(missionList,tLaunchDesired,dt)
    trajList,thetag0 = interpolateTrajectoriesFromMatrices(trajMatrixList,tLaunchDesired,dt)
    r0 = orbitTools.latlonalt2ECI(latitude,longitude,height,thetag0,planetmodel)
    
    nEvents = getNumberOfStages(mission0)
    cloption = 1
    cl = 0.
    minfcl = [1]
    loverd = 0.0
    
    geoptions = 0
    
    dtinterval = dt
    
    ndtinterval = 20000
    fileList = []
    currentTime = 0.0
    retList= []
    for indexEvent in range(nEvents):
        if indexEvent ==0 : #first stage...setting initial conditions for rest of cases
            
            [vS0,vE0,vZ0] = [0,0,0] # initial velocity in SEZ frame
            inertial = 0
            thetag = copy.deepcopy(thetag0)
            # inertial = 1. input velocities are inertial velocities in SEZ frame
            # inertial =0. input velocities are local velocities (earth rotation not accounted)
            #Veci0 = orbitTools.SEZ2ECI(latitude,longitude,r0,vS0,vE0,vZ0,inertial,thetag)
            Veci0 = orbitTools.SEZ2ECI(latitude,longitude,height,vS0,vE0,vZ0,inertial,thetag,planetmodel)
            initialstate = np.array([r0[0],r0[1],r0[2],Veci0[0],Veci0[1],Veci0[2]])
        
        
        m0,Fmax,isp,minfcd,cd,sref = vehicleParamFromMission(mission,eventIndex)
        Fmax = Fmax*offset
        
        etaMin = getEtaMinFromMission(mission,indexEvent)
        filename = 'pTraj'+str(indexEvent)+'.txt'
        fileList.append(filename)
        if etaMin>=0:
            
            propOption = 1
            currTraj = trajList[indexEvent]
            [indexEvent,tfit,uxFit,uyFit,uzFit,etaFit] = currTraj
            # ensuring that unit vector is actually a unit vector (due to interpolation)
            
            uMag = (uxFit**2 + uyFit**2 + uzFit**2)**.5
            uxetaFit = np.array([etaFit*uxFit/uMag]).T
            uyetaFit = np.array([etaFit*uyFit/uMag]).T
            uzetaFit = np.array([etaFit*uzFit/uMag]).T
            uetaMat = np.concatenate((uxetaFit,uyetaFit,uzetaFit),1)
            Tmat = Tmax*uetaMat 
            timelist = tfit
            propCond = massf
            ntime = len(timelist)
            
            Rearth = getRfromMission(mission0)
            omega = getAngVelEarth(mission0)
            
            
            finalconditions,finalderivs = op.propagate(initialstate, m0, propCond,sref,minfcd,cd,cloption, minfcl,cl,loverd,atmosoption, altitudelist,densitylist,ulist,vlist ,wlist,Tmat,timelist,isp,geoptions,filename,planetmodel,dtinterval,ndtinterval,thetag,propOption,ThrustOffsetAngDeg,ncd=len(minfcd),ncl=len(minfcl),ntime=ntime,nlist=len(altitudelist))
            xc,yc,zc,Vxc,Vyc,Vzc,mc,tc,newfinalConditions = fixPropagateResults(finalconditions)

            newIndex = finalconditions[:,10]>=0.0
            newfinalConditions = finalconditions[newIndex,:]
            # updating parameters for propagation in next stage

            tc = tc + currentTime
            initialstate = np.array([xc[-1],yc[-1],zc[-1],Vxc[-1],Vyc[-1],Vzc[-1]])
            mass0 = mc[-1]   
            
            currentTime = tc[-1] # updating time
            thetagVector = omegaE*tc + thetag0
            thetag = omegaE*currentTime + thetag0 # updating thetag for next event
        
        elif etaMin<0:
            
            propOption = 2
            currTraj = trajList[indexEvent]
            [indexEvent,tffit] = currTraj
            dtlocal = tffit/5.
            Tmax = 0.0
            uetaMat = np.array([[1,1,1]])
            Tmat = Tmax*uetaMat 
            timelist = [tffit]
            propCond = tffit
            ntime = len(timelist)
            
            
            
            finalconditions,finalderivs = op.propagate(initialstate, mass0, propCond,sref,minfcd,cd,cloption, minfcl,cl,loverd,atmosoption, altitudelist,densitylist,ulist,vlist ,wlist,Tmat,timelist,isp,geoptions,filename,planetmodel,dtlocal,ndtinterval,thetag,propOption,ThrustOffsetAngDeg,ncd=len(minfcd),ncl=len(minfcl),ntime=ntime,nlist=len(altitudelist))
            newIndex = finalconditions[:,10]>=0.0
            newfinalConditions = finalconditions[newIndex,:]
            # updating parameters for propagation in next stage
            xc = newfinalConditions[:,0]
            yc = newfinalConditions[:,1]
            zc = newfinalConditions[:,2]
            Vxc = newfinalConditions[:,3]
            Vyc = newfinalConditions[:,4]
            Vzc = newfinalConditions[:,5]      
            mc = newfinalConditions[:,9]
            tc = newfinalConditions[:,10] + currentTime
            initialstate = np.array([xc[-1],yc[-1],zc[-1],Vxc[-1],Vyc[-1],Vzc[-1]])
            mass0 = mc[-1]   
            currentTime = tc[-1] # updating time
            thetagVector = omegaE*tc + thetag0
            thetag = omegaE*currentTime + thetag0# updating thetag for next event
        
        retList.append([xc,yc,zc,Vxc,Vyc,Vzc,mc,tc,thetagVector,indexEvent])
    return retList
#data2GE.convertMultiple(fileList)

def generateRandomStateVectorfromMatrix(missionList,trajMatrixList,tLaunchDesired,tfail,dt,atmProfile=[],TOffset = 1.0,ThrustOffsetAngDeg=[0,0]):
    from scipy import interpolate
    #retList = generateRandomTraj(missionList,tLaunchDesired,dt,atmProfile,TOffset,ThrustOffsetAngDeg)  
    retList = generateRandomTrajfromMatrix(missionList,trajMatrixList,tLaunchDesired,dt,atmProfile,TOffset ,ThrustOffsetAngDeg)
    for index in range(len(retList)):
        [xc,yc,zc,Vxc,Vyc,Vzc,mc,tc,thetagVector,indexEvent] = retList[index]
        
        if (index==(len(retList)-1)) and tc[-1]<tfail:
            # ensuring that values are bounded...this should not be an issue because late in the trajectory the IIP should be in orbit....no chance of hitting ground
            tfail = tc[-1]
            print 'Warning in getPropTraj.py. Failure time is greater than modeled trajectory time.'
            print 'Adjusting failure time to maximum trajectory time available... tfail is now ',tfail,' secs'
        
        if tc[-1]>=tfail : # look in this event (stage)
            fx = interpolate.interp1d(tc,xc,kind='linear')
            fy = interpolate.interp1d(tc,yc,kind='linear')
            fz = interpolate.interp1d(tc,zc,kind='linear')
            fVx = interpolate.interp1d(tc,Vxc,kind='linear')
            fVy = interpolate.interp1d(tc,Vyc,kind='linear')
            fVz = interpolate.interp1d(tc,Vzc,kind='linear')
            fm = interpolate.interp1d(tc,mc,kind='linear')
            fthetag = interpolate.interp1d(tc,thetagVector,kind='linear')
            [x,y,z,Vx,Vy,Vz] = [fx(tfail),fy(tfail),fz(tfail),fVx(tfail),fVy(tfail),fVz(tfail)]
            stateVec = [x,y,z,Vx,Vy,Vz]
            mass = fm(tfail)
            thetag = fthetag(tfail)
            return (stateVec,thetag,mass,indexEvent)

