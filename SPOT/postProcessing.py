import variables
from constraints import computeNumericMass
import numpy as np
import matplotlib.pyplot as plt 

def RecordSolution(x_opt,mission):
    x,y,z,Vx,Vy,Vz,theta1,theta2,eta,tf = variables.decomposeVariables(x_opt,mission['numerics']['range'])
    mission['solution']['x'] = mission['scale']['R']*x
    mission['solution']['y'] = mission['scale']['R']*y
    mission['solution']['z'] = mission['scale']['R']*z
    mission['solution']['Vx'] = mission['scale']['V']*Vx
    mission['solution']['Vy'] = mission['scale']['V']*Vy
    mission['solution']['Vz'] = mission['scale']['V']*Vz
    
    mission['solution']['theta1'] = theta1
    mission['solution']['theta2'] = theta2

    
    ux = np.sin(theta1)*np.cos(theta2)
    uy = np.sin(theta1)*np.sin(theta2)
    uz = np.cos(theta1)
    mission['solution']['ux'] = ux
    mission['solution']['uy'] = uy
    mission['solution']['uz'] = uz
    mission['solution']['eta'] = eta
    tState,tControl = getTimeVector(mission,tf)
    mission['solution']['tState'] = mission['scale']['t']*tState
    mission['solution']['tControl'] = mission['scale']['t']*tControl
    mission['solution']['m'] = computeNumericMass(eta,tf,mission,'real')*mission['scale']['m']
    mission['solution']['mf'] = mission['solution']['m'][-1]
    mission['solution']['tf'] = tf*mission['scale']['t']
    mission['solution']['margin'] = mission['solution']['mf'] - mission['vehicle']['mass']['mf']
    mission['solution']['N'] = mission['numerics']['N']


                
    return mission

def RecordSolutionPerStage(x_opt,mission):
    import orbitPropSPOT as op

    x,y,z,Vx,Vy,Vz,theta1,theta2,eta,tf = variables.decomposeVariables(x_opt,mission['numerics']['range'])
    mission['solution']['x'] = mission['scale']['R']*x
    mission['solution']['y'] = mission['scale']['R']*y
    mission['solution']['z'] = mission['scale']['R']*z
    mission['solution']['Vx'] = mission['scale']['V']*Vx
    mission['solution']['Vy'] = mission['scale']['V']*Vy
    mission['solution']['Vz'] = mission['scale']['V']*Vz
    mission['solution']['theta1'] = theta1
    mission['solution']['theta2'] = theta2
    
    
    ux = np.sin(theta1)*np.cos(theta2)
    uy = np.sin(theta1)*np.sin(theta2)
    uz = np.cos(theta1)
    
    
    mission['solution']['ux'] = ux
    mission['solution']['uy'] = uy
    mission['solution']['uz'] = uz
    mission['solution']['eta'] = eta
    tState,tControl = getTimeVector(mission,tf)
    mission['solution']['tState'] = mission['scale']['t']*tState
    mission['solution']['tControl'] = mission['scale']['t']*tControl
    mission['solution']['m'] = computeNumericMass(eta,tf,mission,'real')*mission['scale']['m']
    mission['solution']['mf'] = mission['solution']['m'][-1]
    mission['solution']['tf'] = tf*mission['scale']['t']
    mission['solution']['margin'] = mission['solution']['mf'] - mission['vehicle']['mass']['mf']
    mission['solution']['N'] = mission['numerics']['N']
    
    print 'Time for each stage is ',mission['solution']['tf'],' [sec]'
    mission['solutionStage']['thetag'] =  mission['launchSite']['thetaLaunch']
    newTime = mission['solutionStage']['t'][0][-1]+ mission['solution']['tState']
    
    # storing everything in a nice list so that each stage can be access quickly
    for stage in range(int(mission['vehicle']['stages'])):
        rangeState = mission['numerics']['range']['stageState'][stage]
        rangeControl = mission['numerics']['range']['stageControl'][stage]
        
        if stage ==0: # adding points from rk45 straight up trajectory
            
            mission['solutionStage']['x'][stage] = np.concatenate((mission['solutionStage']['x'][stage][0:-1],mission['solution']['x'][rangeState]),2)
            mission['solutionStage']['y'][stage] = np.concatenate((mission['solutionStage']['y'][stage][0:-1],mission['solution']['y'][rangeState]),2)
            mission['solutionStage']['z'][stage] = np.concatenate((mission['solutionStage']['z'][stage][0:-1],mission['solution']['z'][rangeState]),2)
            mission['solutionStage']['Vx'][stage] = np.concatenate((mission['solutionStage']['Vx'][stage][0:-1],mission['solution']['Vx'][rangeState]),2)
            mission['solutionStage']['Vy'][stage] = np.concatenate((mission['solutionStage']['Vy'][stage][0:-1],mission['solution']['Vy'][rangeState]),2)
            mission['solutionStage']['Vz'][stage] = np.concatenate((mission['solutionStage']['Vz'][stage][0:-1],mission['solution']['Vz'][rangeState]),2)        
            mission['solutionStage']['ux'][stage] = np.concatenate((mission['solutionStage']['ux'][stage][0:-1],mission['solution']['ux'][rangeControl]),2)
            mission['solutionStage']['uy'][stage] = np.concatenate((mission['solutionStage']['uy'][stage][0:-1],mission['solution']['uy'][rangeControl]),2)
            mission['solutionStage']['uz'][stage] = np.concatenate((mission['solutionStage']['uz'][stage][0:-1],mission['solution']['uz'][rangeControl]),2)
            mission['solutionStage']['eta'][stage] = np.concatenate((mission['solutionStage']['eta'][stage][0:-1],mission['solution']['eta'][rangeControl]),2)
            mission['solutionStage']['t'][stage] = np.concatenate((mission['solutionStage']['t'][stage][0:-1],newTime[rangeControl]),2)
            mission['solutionStage']['m'][stage] = np.concatenate((mission['solutionStage']['m'][stage][0:-1],mission['solution']['m'][rangeControl]),2)
            print 'Calculated Mass at the end of stage ',stage+1,' = ',mission['solutionStage']['m'][stage][-1] ,' [kg]'

        
        else:
            mission['solutionStage']['x'][stage] = mission['solution']['x'][rangeState]
            mission['solutionStage']['y'][stage] = mission['solution']['y'][rangeState]
            mission['solutionStage']['z'][stage] = mission['solution']['z'][rangeState]
            mission['solutionStage']['Vx'][stage] = mission['solution']['Vx'][rangeState]
            mission['solutionStage']['Vy'][stage] = mission['solution']['Vy'][rangeState]
            mission['solutionStage']['Vz'][stage] = mission['solution']['Vz'][rangeState]
            mission['solutionStage']['t'][stage] = newTime[rangeState]

            if len(rangeControl)>0:
                mission['solutionStage']['eta'][stage] = mission['solution']['eta'][rangeControl]    
                mission['solutionStage']['theta1'][stage] = mission['solution']['theta1'][rangeControl] 
                mission['solutionStage']['theta2'][stage] = mission['solution']['theta2'][rangeControl]    
                theta1 = mission['solutionStage']['theta1'][stage]
                theta2 = mission['solutionStage']['theta2'][stage]
            
                ux = np.sin(theta1)*np.cos(theta2)
                uy = np.sin(theta1)*np.sin(theta2)
                uz = np.cos(theta1)
                mission['solutionStage']['ux'][stage] = ux
                mission['solutionStage']['uy'][stage] = uy
                mission['solutionStage']['uz'][stage] = uz
                
                
                mission['solutionStage']['m'][stage] = mission['solution']['m'][rangeControl]
                print 'Calculated Mass at the end of stage ',stage+1,' = ',mission['solutionStage']['m'][stage][-1] ,' [kg]'
    print 'Physical Mass at the end of each stage =' ,mission['vehicle']['mass']['mf'],' [kg]'

    print 'Final Mass is ',mission['solution']['mf'],' [kg]'
    print 'Final Physical Mass (structural + payload) is ', mission['vehicle']['mass']['mf'][-1],' [kg]'
    
                    
    # propagating any desired trajectory forward 
    if mission['options']['propagateStage'][0]>0:
        mission['solutionStage']['propagateLast']['x'] = [-999]*len(mission['options']['propagateStage'])
        mission['solutionStage']['propagateLast']['y'] = [-999]*len(mission['options']['propagateStage'])
        mission['solutionStage']['propagateLast']['z'] = [-999]*len(mission['options']['propagateStage'])
        mission['solutionStage']['propagateLast']['Vx'] = [-999]*len(mission['options']['propagateStage'])
        mission['solutionStage']['propagateLast']['Vy'] = [-999]*len(mission['options']['propagateStage'])
        mission['solutionStage']['propagateLast']['Vz'] = [-999]*len(mission['options']['propagateStage'])
        mission['solutionStage']['propagateLast']['t'] = [-999]*len(mission['options']['propagateStage'])
        mission['solutionStage']['propagateLast']['latitude'] = [-999]*len(mission['options']['propagateStage'])
        mission['solutionStage']['propagateLast']['longitude'] = [-999]*len(mission['options']['propagateStage'])
        mission['solutionStage']['propagateLast']['height'] = [-999]*len(mission['options']['propagateStage'])
        mission['solutionStage']['propagateLast']['mass'] = [-999]*len(mission['options']['propagateStage'])
        

        indval = -1
        for stageval in mission['options']['propagateStage']:
            stage = int(stageval)-1
            indval = indval+1 
            thetag = mission['launchSite']['thetaLaunch'] + mission['planet']['omega']*mission['solutionStage']['t'][stage][-1]
            #getting last known state vectors from optimal trajectory
            propOption = 2
            propCond = 8000.0 # [seconds] to stop propagation




            #mass = mission.vehicle.mass.structural[stage]
            mass = mission['solutionStage']['m'][stage][-1]
            xD = mission['solutionStage']['x'][stage][-1]
            yD = mission['solutionStage']['y'][stage][-1]
            zD = mission['solutionStage']['z'][stage][-1]
            VxD = mission['solutionStage']['Vx'][stage][-1]
            VyD = mission['solutionStage']['Vy'][stage][-1]
            VzD = mission['solutionStage']['Vz'][stage][-1]
            initialstate = np.array([xD,yD,zD,VxD,VyD,VzD])
            
        
            planetmodel = 0
            
            
            #debrisvel = np.array([0,0,0])
            
            sref = mission['vehicle']['aero']['Aref'][stage]
            minfcd = mission['vehicle']['aero']['MinfCD'][stage]
            cd = mission['vehicle']['aero']['CD'][stage]#1.0
            cloption = 1
            cl = 0.
            minfcl = [1]
            loverd = 0.0
            atmosoption = 0
            densitylist = [0]
            ulist = [0]
            vlist = [0]
            wlist = [0]
            altitudelist = [0]
            geoptions = 0
            filename = 'PropagationStage'+str(stage)+'.txt'
            dtinterval = 1
            #print 'Thrust' ,mission.vehicle.propulsion.Ftotal[0]
            thrust = mission['vehicle']['propulsion']['Ftotal'][0]
            #print 'ratio',thrust/(mass*9.81) 
            
            thrustList = np.array([[0,0,0]])
            timelist = np.array([0])
            ntime = len(timelist)
            isp = 1.0 # does not matter since no thrust!
            ndtinterval = 10000
            '''
            #print thrustList
            print cd
            print minfcd
            print 'R val ',np.sqrt(xD**2+yD**2+zD**2)
            print 'V val ',np.sqrt(VxD**2+VyD**2+VzD**2)
            print mass
            '''
            finalconditions,finalderivs = op.propagate(initialstate,mass,propCond,sref,minfcd,cd,cloption,minfcl,cl,loverd,atmosoption,altitudelist,densitylist,ulist,vlist,wlist,thrustList,timelist,isp,geoptions,filename,planetmodel,dtinterval,ndtinterval,thetag,propOption,ncd=len(minfcd),ncl=len(minfcl),ntime=ntime,nlist=len(altitudelist))
            
            newIndex = finalconditions[:,10]>=-.01
            newfinalConditions = finalconditions[newIndex,:]
            mission['solutionStage']['propagateLast']['x'][indval] = newfinalConditions[:,0]
            mission['solutionStage']['propagateLast']['y'][indval] = newfinalConditions[:,1]
            mission['solutionStage']['propagateLast']['z'][indval] = newfinalConditions[:,2]
            mission['solutionStage']['propagateLast']['Vx'][indval] = newfinalConditions[:,3]
            mission['solutionStage']['propagateLast']['Vy'][indval] = newfinalConditions[:,4]
            mission['solutionStage']['propagateLast']['Vz'][indval] = newfinalConditions[:,5]
            mission['solutionStage']['propagateLast']['longitude'][indval] = newfinalConditions[:,6]
            mission['solutionStage']['propagateLast']['latitude'][indval] = newfinalConditions[:,7]
            mission['solutionStage']['propagateLast']['height'][indval] = newfinalConditions[:,8]
            mission['solutionStage']['propagateLast']['mass'][indval] = newfinalConditions[:,9]
            mission['solutionStage']['propagateLast']['t'][indval] = newfinalConditions[:,10]

    
    return mission

def getTimeVector(mission,tf):
    tState = np.zeros(mission['numerics']['Nstate'])
    tControl = np.zeros(mission['numerics']['Ncontrol'])
    rangeState = mission['numerics']['range']['stageState']
    rangeControl = mission['numerics']['range']['stageControl']
    tflast = 0.0
    for stage in range(mission['vehicle']['stages']):
        N = mission['numerics']['N'][stage]
        t = mission['numerics']['t'][stage]
        tcurrent = tf[stage]*t+tflast
        tState[rangeState[stage]] = tcurrent
        tControl[rangeControl[stage]] = tcurrent

        tflast = tcurrent[-1]
    #print 'Tcont ',tControl
    #print 'Tstate ',tState
    return (tState,tControl)





def plotting(x_opt,mission):
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

# added some simple plots
    
    x,y,z=variables.decomposeLocationVariables(x_opt,mission['numerics']['range'])


    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.plot(x,y,z)
    u, v = np.mgrid[0:2*np.pi:200j, 0:np.pi:200j]
    x1=np.cos(u)*np.sin(v)
    y1=np.sin(u)*np.sin(v)
    z1=np.cos(v)
    ax.plot_surface(x1, y1, z1, color="r")

    plt.figure()
    plt.plot(mission['solution']['m'])
    
    plt.figure()
    plt.plot(mission['solution']['eta'])
    plt.show()



def convert2LatLon(mission):
    xT = mission['solution']['x']
    yT = mission['solution']['y']
    zT = mission['solution']['z']

    VxT = mission['solution']['Vx']
    VyT = mission['solution']['Vy']
    VzT = mission['solution']['Vz']

    tT = mission['solution']['tState']
    thetag = mission['launchSite']['thetag']
    revPerDay = mission['VernalEquinox']['revPerDay']
    Nstages = mission['vehicle']['stages']
    Ntotal = mission['numerics']['Ncontrol']
    Nstages = mission['vehicle']['stages']
    Ntotal = mission['numerics']['Ncontrol']
    xGreen = np.zeros(Ntotal)
    yGreen = np.zeros(Ntotal)
    zGreen = np.zeros(Ntotal)
    VxGreen = np.zeros(Ntotal)
    VyGreen = np.zeros(Ntotal)
    VzGreen = np.zeros(Ntotal)
    latitude = np.zeros(Ntotal)
    longitude = np.zeros(Ntotal)
    Rdistance = np.zeros(Ntotal)
    Vdown = np.zeros(Ntotal)
    Veast = np.zeros(Ntotal)
    Vnorth = np.zeros(Ntotal)

    maincount = -1
    omega = mission['planet']['omega']
    for stage in range(Nstages):
        rangeState = mission['numerics']['range']['stageState'][stage]
        x = xT[rangeState]
        y = yT[rangeState]
        z = zT[rangeState]
        Vx = VxT[rangeState]
        Vy = VyT[rangeState]
        Vz = VzT[rangeState]
        tloc = tT[rangeState]

        for counter in range(len(tloc)):
            t = tloc[counter]
            xyzIJK = np.array([[x[counter]],[y[counter]],[z[counter]]])
            VxVyVzIJK = np.array([[Vx[counter]],[Vy[counter]],[Vz[counter]]])
            theta = thetag + omega*t
            Rz = Rotatez(theta)
            xyzGreen = np.dot(Rz,xyzIJK)
            VxVyVzGreen = np.dot(Rz,VxVyVzIJK)
            Vcirc = np.array([[-xyzGreen[1]*omega],[xyzGreen[0]*omega],[0]])

            maincount = maincount + 1
            xGreen[maincount] = xyzGreen[0]
            yGreen[maincount] = xyzGreen[1]
            zGreen[maincount] = xyzGreen[2]
            lon,lat,Rloc = cart2sph(xyzGreen[0],xyzGreen[1],xyzGreen[2])
            longitude[maincount] = lon*180./np.pi
            latitude[maincount] = lat*180./np.pi
            Rdistance[maincount] = Rloc
    mission['solution']['latitude'] = latitude
    mission['solution']['longitude'] = longitude
    mission['solution']['height'] = Rdistance - mission['scale']['R']

    return mission

def convert2LatLonPerStage(mission):
    Nstages = mission['vehicle']['stages']
    thetag = mission['solutionStage']['thetag']
    omega = mission['planet']['omega']

    for stage in range(Nstages):
        x = mission['solutionStage']['x'][stage]
        y = mission['solutionStage']['y'][stage]
        z = mission['solutionStage']['z'][stage]
        tloc = mission['solutionStage']['t'][stage]
        latitude = np.zeros(len(tloc))
        longitude = np.zeros(len(tloc))
        height = np.zeros(len(tloc))
        for counter in range(len(tloc)):
            t = tloc[counter]

            theta = thetag + omega*t
            xyzIJK = np.array([[x[counter]],[y[counter]],[z[counter]]])
            Rz = Rotatez(theta)
            xyzGreen = np.dot(Rz,xyzIJK)
            lon,lat,Rloc = cart2sph(xyzGreen[0],xyzGreen[1],xyzGreen[2])
            latitude[counter] = lat*180.0/np.pi
            longitude[counter] = lon*180.0/np.pi
            height[counter] = Rloc - mission['planet']['R']
        mission['solutionStage']['height'][stage] = height
        mission['solutionStage']['latitude'][stage] = latitude
        mission['solutionStage']['longitude'][stage] = longitude
        
            #print 'Theta is ',thetag
        #plt.figure()
        #plt.plot(x)
        #plt.plot(y)
        #plt.plot(z)

    #plt.show()
    return mission






def cart2sph(x,y,z):
    theta = np.arctan2(y,x)
    phi = np.arctan2(z,np.sqrt(x**2+y**2))
    r = np.sqrt(x**2+y**2+z**2)

    return (theta,phi,r)

def Rotatez(theta):
    R = np.array([[np.cos(theta),np.sin(theta),0.0],[-np.sin(theta),np.cos(theta),0],[0,0,1.]])
    return R

def Rotaty(theta):
    R = np.array([[np.cos(theta),0,-np.sin(theta)],[0,1,0],[np.sin(theta),0,np.cos(theta)]])
    return R
    




def UpdateMissionAndInitialGuess(mission):
    import numpy as np
    import variables
    from LaunchDataProcessing import NumericData
    from scipy.interpolate import interp1d

    # First, need to update the relevant parameters in mission according the the iteration schedule
    curIter = mission['numerics']['iterationSchedule']['curIter']
    mission['numerics']['N'] = mission['numerics']['iterationSchedule']['N'][curIter]
    mission = NumericData(mission)
    # ----------- Done Updating Misison --------------

    # Generate an updated initial guess based on previous solution
    Nstate = mission['numerics']['Nstate']
    Ncontrol = mission['numerics']['Ncontrol']
        
    tf = mission['solution']['tf']/mission['scale']['t']
    tState,tControl = getTimeVector(mission,tf)
    scaleR = mission['scale']['R']
    scaleV = mission['scale']['V']
    tStateSol = mission['solution']['tState']/mission['scale']['t']
    tControlSol = mission['solution']['tControl']/mission['scale']['t']
    
    # making sure end points are the same, sometimes there is a 1e-15 difference, enough to show an extrapolation error
    tState[-1] = tStateSol[-1] 
    tControl[-1] = tControlSol[-1]

    x = interp1d(tStateSol, mission['solution']['x']/scaleR)(tState)
    y = interp1d(tStateSol, mission['solution']['y']/scaleR)(tState)
    z = interp1d(tStateSol, mission['solution']['z']/scaleR)(tState)
    Vx = interp1d(tStateSol, mission['solution']['Vx']/scaleV)(tState)
    Vy = interp1d(tStateSol, mission['solution']['Vy']/scaleV)(tState)
    Vz = interp1d(tStateSol, mission['solution']['Vz']/scaleV)(tState)

    theta1 = interp1d(tControlSol, mission['solution']['theta1'])(tControl)
    theta2 = interp1d(tControlSol, mission['solution']['theta2'])(tControl)
    
    eta = interp1d(tControlSol, mission['solution']['eta'])(tControl)
    tf = tf; #Keep the old tf, no change needed

    vars0 = variables.assembleVariables(x,y,z,Vx,Vy,Vz,theta1,theta2,eta,tf)
    
    
    return mission, vars0

def writeTextFiles(missionList):
    for index in range(len(missionList)):
        fileOut = open('Traj'+str(index+1)+'.txt','w')
        mission = missionList[index]
        nstages = mission['vehicle']['stages']
        fileOut.write('Date for this Trajectory\n')
        fileOut.write(str(mission['launchSite']['month']) + '/'+str(mission['launchSite']['day'])+'/'+str(mission['launchSite']['year'])+'\n')
        fileOut.write('Local Time (hour-minute): '+str(mission['launchSite']['localTime'][0])+'-'+str(mission['launchSite']['localTime'][1])+'\n')
        fileOut.write('Thetag0 : '+ str(mission['launchSite']['thetaLaunch'])+'\n')
        stateVectors = ['t','m','longitude','latitude','height','x','y','z','Vx','Vy','Vz','ux','uy','uz','eta']
        for stage in range(nstages):
            fileOut.write('Stage '+str(stage+1)+ ' Results\n')
            fileOut.write('\n   Results\n')
            for stateIndex in range(len(stateVectors)):
                fileOut.write('\n      *'+stateVectors[stateIndex]+'*\n')
                states = mission['solutionStage']
                for indexXYZ in range(len(states[stateVectors[stateIndex]][stage])):
                    fileOut.write('        '+str(states[stateVectors[stateIndex]][stage][indexXYZ])+'\n')
        fileOut.close()
            

