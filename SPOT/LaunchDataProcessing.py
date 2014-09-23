#
#   Francisco Capristan, ADL, Stanford University 
#   Last Modified : Francisco Capristan 
#   Started:        7/15/11
#
#   Inputs:         mission                     class structure
#
#   Outputs:        mission                     updated structure class 
#

#################################################################
import numpy as np
from scipy import interpolate

def processAllMissionData(mission,flagthetag = False):
    mission = massPropsData(mission,'derived') # initial test passed
    mission = planetaryData(mission) # initial test passed

    mission = scaleFactors(mission) # initial test passed
    mission = massPropsData(mission,'scale') # initial test passed
    mission = NumericData(mission) # initial test passed
    mission = PropulsionData(mission) # initial test passed

    mission = LaunchSiteData(mission,flagthetag) # initial test passed

    mission = OrbitData(mission) # initial test passed

    mission = aeroData(mission)
    return mission
    
def massPropsData(mission,flag):
    # compute derived mass property and non-dimensionalize
    #   based on MassPropsData.m by Michael Colonno


    mass = mission['vehicle']['mass']

    if flag=='derived':
        
        # stage masses
        nend = mission['vehicle']['stages']
        mass['m0'] = np.zeros((mission['vehicle']['stages']))
        mass['mf'] = np.zeros((mission['vehicle']['stages']))
        for index in range(0,mission['vehicle']['stages']):
            mass['m0'][index] = np.sum(mass['propellant'][index:nend]) + np.sum(mass['structural'][index:nend]) + mass['payload']         # kg
            mass['mf'][index] = np.sum(mass['propellant'][index+1:nend]) + np.sum(mass['structural'][index:nend]) + mass['payload']       # kg
        mass['nonDim']['MR'] = mass['m0']/mass['mf']                                                                   # non-dim

    elif flag=='scale':
        
        mass['nonDim']['m0'] = mass['m0']/mission['scale']['m']
        mass['nonDim']['mf'] = mass['mf']/mission['scale']['m']
        mass['nonDim']['mS'] = mass['structural']/mission['scale']['m']
        mass['nonDim']['mP'] = mass['propellant']/mission['scale']['m']
        mass['nonDim']['mPL'] = mass['payload']/mission['scale']['m']

    else:
        
        print 'ERROR in LaunchDataProcessing.py: Unknown flag'
        exit()




    mission['vehicle']['mass'] = mass 

    return mission
    #################################################################



def planetaryData(mission):
    #computes requried planetary data

    
    G = 6.67384e-11                                                # universal gravity (m^3/kg-s^2)
    mission['planet']['mu'] = G*mission['planet']['m']                        # planet grav. parameter (m^3/s^2)                     
    mission['planet']['g0'] = mission['planet']['mu']/(mission['planet']['R']**2)      # reference gravity (m/s^2)
    mission['planet']['omega'] = 2*np.pi/mission['planet']['T']                   # planetary rotation rate (rad/s)
    #print mission['planet']['omega

    return mission




def scaleFactors(mission):
    
    # Calculate quantities used for non-dimensionalization
    
    mission['scale']['R'] = mission['planet']['R']                                 # m
    mission['scale']['m'] = mission['vehicle']['mass']['m0'][0]                # kg picking total initial first stage mass to non-dimensionalize
    mission['scale']['t'] = (mission['planet']['R']**3/mission['planet']['mu'])**.5      # s
    mission['scale']['V'] = mission['scale']['R']/mission['scale']['t']                  # m/s

    return mission

#################################################################




#   NumericData.m: Compute Chebyshev differentiation & integration matricies, non-dimensional time vector, 
#                   variable and constraint counts & ranges
#   Modified : Francisco Capristan
#
#   Michael Colonno, ADL, Stanford University 
#
#   Modified : Francisco Capristan
#   Started:        7/2/11
#   Last Updated:   5/01/12
#
#   Inputs:         mission               
#
#   Outputs:        t           cosine-spaced time vector
#                   D           differentiation matrix 
#                   I           integration matrix 
#                   nvars       total number of variables 
#                   ['range']      indices for variable specific variables 
#   Notes : it creates t,D,I vectors/matrices in 3D to handle different
#   stages
#################################################################

def NumericData(mission):
    import CreateChebyshevD as CC
    
    etavec = mission['vehicle']['propulsion']['minThrottle'] 

    Ntotal0 = int(np.sum(mission['numerics']['N'])) # total number of nodes, including intersections
    print 'Cheb points used = ',mission['numerics']['N']
    numberofStages = int(mission['vehicle']['stages'])
    # accounting
    Nstate = int(Ntotal0-(numberofStages-1))
    Ntotal = int(np.sum(mission['numerics']['N'][etavec>0])) # total number of nodes, including intersections

    mission['numerics']['Ncontrol'] = int(Ntotal)
    mission['numerics']['Nstate'] = int(Nstate)
    
    mission['numerics']['nvars'] = int( 6*Nstate+3*Ntotal + mission['vehicle']['stages'])
    

    Nstate = int(Nstate)
    # variable ranges
    
    # setting state vector indexing
    rangevar = np.arange(0,Nstate)
    mission['numerics']['range']['x'] = rangevar
    rangevar = rangevar + Nstate
    mission['numerics']['range']['y'] = rangevar
    rangevar = rangevar + Nstate
    mission['numerics']['range']['z'] = rangevar
    rangevar = rangevar + Nstate
    mission['numerics']['range']['Vx']= rangevar
    rangevar = rangevar+Nstate
    mission['numerics']['range']['Vy'] = rangevar
    rangevar = rangevar+Nstate
    mission['numerics']['range']['Vz']  = rangevar
    # done with state vector indexing
    # starting indexing for control variable ux uy uz eta
    rangevar = np.arange(rangevar[Nstate-1]+1,rangevar[Nstate-1]+1+Ntotal) # make that this is equivalent to matlab indexing
    mission['numerics']['range']['theta1'] = rangevar
    rangevar = rangevar+Ntotal
    mission['numerics']['range']['theta2'] = rangevar
    rangevar = rangevar+Ntotal
    #mission['numerics']['range']['uz'] = rangevar
    #rangevar = rangevar+Ntotal
    mission['numerics']['range']['eta'] = rangevar
    mission['numerics']['range']['tf'] = np.arange((rangevar[len(rangevar)-1]+1),(rangevar[len(rangevar)-1]+mission['vehicle']['stages']+1))# need to fix this line and below
    
    
    upperRangeEQMotion = -1
    mission['numerics']['t'] = range(numberofStages)
    mission['numerics']['D'] = range(numberofStages)
    mission['numerics']['dthalf'] = range(numberofStages)
    mission['numerics']['range']['stageState'] =  range(numberofStages)
    mission['numerics']['range']['stageControl'] =  range(numberofStages)
    mission['numerics']['range']['EQ'] = range(numberofStages)
    Nold = 0
    for stageNumber in  range(0,numberofStages):
        N = int(mission['numerics']['N'][stageNumber])
        
        mission['numerics']['t'][stageNumber], mission['numerics']['D'][stageNumber],err = CC.ChebyshevMatrices(N)
        t = mission['numerics']['t'][stageNumber]
        dthalf = np.array(.5*(t[1:N]-t[0:N-1]))
        mission['numerics']['dthalf'][stageNumber] = np.matrix(np.concatenate(([0.0], dthalf),1)).T
        
        
        if stageNumber==0:
            lowerRangeState = 0
            upperRangeState = N-1       
        
        else :
            lowerRangeState = upperRangeState
            upperRangeState = lowerRangeState+N-1
                    
        rangeSt = np.arange(lowerRangeState,upperRangeState+1)
        mission['numerics']['range']['stageState'][stageNumber] = rangeSt#range for state vectors
        if etavec[stageNumber]>=0:
            mission['numerics']['range']['stageControl'][stageNumber] = np.arange(0,N) + Nold
            Nold = N + Nold

            #mission['numerics']['range'].stageControl[stageNumber] = rangeSt + stageNumber#range for other control parameters
        else:
            mission['numerics']['range']['stageControl'][stageNumber] = []    
        #Nold = N + Nold

        lowerRangeEQMotion = upperRangeEQMotion+1
        upperRangeEQMotion = lowerRangeEQMotion+6*N-1
        rangeEQ = np.arange(lowerRangeEQMotion,upperRangeEQMotion+1)

        mission['numerics']['range']['EQ'][stageNumber]=rangeEQ
    

    ## This section is just calculated to speed up gradient in the objective function
    ## used tp calculate change in final mass wrt to eta
    N = int(mission['numerics']['N'][-1])
    dmfdeta = np.zeros((N))
    tn = mission['numerics']['t'][-1]
    index1 = range(1,N-1) # equivalent to 1:(N-2) matlab type
    index2 = range(2,N)
    index3 = range(0,N-2)
    dmfdeta[index1] = -.5*(tn[index2]-tn[index3])
    dmfdeta[0] = -.5*(tn[1]-tn[0])
    dmfdeta[-1] = -.5*(tn[N-1]-tn[N-2])
  
    mission['numerics']['dmfdetaNondim'] = dmfdeta
    #print mission['numerics']['range'].stageState



    #exit()
    return mission

#################################################################





def LaunchSiteData(mission,flagthetag=False):
    from datetime import datetime
    import orbitPropSPOT as op
    import data2GE

#   FROM LaunchSiteData.m: Compute derives launch site data and non-dimensionalize 
#
#   Michael Colonno, ADL, Stanford University 
#   Last Modified : Francisco Capristan (Added Vernal Equinox data)
#   Started:        7/1/11
#   Last Updated:   7/1/11
#
#   Inputs:         mission.launch_site         user launch site data   
#
#   Outputs:        mission.launch_site         updated / non-dimensional data 
#


#################################################################


    mission['launchSite']['nonDim']['r'] = 1 + mission['launchSite']['h']/mission['scale']['R']  # non-dim
    lat = np.pi*mission['launchSite']['lat']/180                                       # rad
    long = np.pi*mission['launchSite']['long']/180                                     # rad

    omega = mission['planet']['omega']



    if flagthetag==False:
        local_time = mission['launchSite']['localTime'] # local launch time [hours, minutes]
        local_time_to_UT = mission['launchSite']['shiftUT'] #hour shift to convert local time to UT [hour]
        year = mission['launchSite']['year']
        month = mission['launchSite']['month']
        day = mission['launchSite']['day']
        mission = RefVernalEquinox(mission)
        refYear = mission['VernalEquinox']['year']
        refDay = mission['VernalEquinox']['day']
        refMonth = mission['VernalEquinox']['month']
        
        a = datetime(year,month,day)
        b = datetime(refYear,refMonth,refDay)
        days = (a-b).days
        timeSEC = getUTsec(local_time,local_time_to_UT)
    #    days = getDaysFromRefVernalEquinox(month,day,year,mission)

    #    timeSEC = getUTsec(local_time,local_time_to_UT) #covert the local time [hours minutes] to UT seconds
        D = days + timeSEC/(3600.0*24.0) # converting to fractions of DAYS ...ref Fund of ASTRO. Bates
        thetag0 = mission['VernalEquinox']['thetag0']
        revPerDay = mission['VernalEquinox']['revPerDay']
        thetag = thetag0 + revPerDay*2.*np.pi*D

    else : # ignore previous thetag and use previously stored thetag
        thetag = mission['launchSite']['thetaLaunch']
        #print 'here'
    #thetag = -2.202536422172285e+04 # temporary...MUST CHANGE!!!!!

  

    ## x0 y0 z0 in IJK frame (see Vernal Equinox)
    
    [x0, y0, z0] = sph2cart(long+thetag,lat,mission['launchSite']['nonDim']['r'])            # non-dim
    #print flagthetag,thetag,x0,y0
    mission['launchSite']['nonDim']['r0'] = np.array([x0, y0, z0])

    ## Using Earth's rotation at Planet's surface as default

    # in V0 in IJK frame (see Vernal Equinox)    
    mission['launchSite']['nonDim']['V0'] = (mission['planet']['omega']*mission['planet']['R']*np.cos(lat)/mission['scale']['V'])*np.array([-np.sin(long+thetag), np.cos(long+thetag),0.])
                                    # non-dim

    ## For Vertical or Horizontal Launches
    finalHeight = mission['launchSite']['straightHeight']
    finalTime = mission['launchSite']['straightTime']
    propagateSwitch = 0
            


    #print 'mass'
    #print newfinalConditions[:,6]/mission['scale']['m']
    #exit()
    # setting up flight section where vehicle goes straight up. And rk type routine is used for propagation straight up
    

    mission['solutionStage']['x'] = int(mission['vehicle']['stages'])*[0]
    mission['solutionStage']['y'] = int(mission['vehicle']['stages'])*[0]
    mission['solutionStage']['z']= int(mission['vehicle']['stages'])*[0]
    mission['solutionStage']['Vx']= int(mission['vehicle']['stages'])*[0]
    mission['solutionStage']['Vy'] = int(mission['vehicle']['stages'])*[0]
    mission['solutionStage']['Vz']  = int(mission['vehicle']['stages'])*[0]
    mission['solutionStage']['theta1']  = int(mission['vehicle']['stages'])*[0]
    mission['solutionStage']['theta2']  = int(mission['vehicle']['stages'])*[0]
    mission['solutionStage']['ux'] = int(mission['vehicle']['stages'])*[0]
    mission['solutionStage']['uy'] = int(mission['vehicle']['stages'])*[0]
    mission['solutionStage']['uz'] = int(mission['vehicle']['stages'])*[0]
    mission['solutionStage']['uz'] = int(mission['vehicle']['stages'])*[0]
    mission['solutionStage']['eta'] = int(mission['vehicle']['stages'])*[0]
    mission['solutionStage']['t'] = int(mission['vehicle']['stages'])*[0]
    mission['solutionStage']['m'] = int(mission['vehicle']['stages'])*[0]
    mission['solutionStage']['latitude'] = int(mission['vehicle']['stages'])*[0]
    mission['solutionStage']['longitude'] = int(mission['vehicle']['stages'])*[0]
    mission['solutionStage']['height'] = int(mission['vehicle']['stages'])*[0]



    if finalTime>0 and finalHeight>0:
        print 'Both, straight time and straight height are specified. Please pick one'
        exit()
    elif finalHeight>0:
        propOption = 3
        propCond = finalHeight
        propagateSwitch = 1
    elif finalTime > 0:
        propOption = 2
        propCond = finalTime
        propagateSwitch = 1

    thetag2 = 1.0*thetag
    if (mission['launchSite']['orientation']=='vertical') and propagateSwitch==1:
        mass = mission['vehicle']['mass']['m0'][0]
        xD = mission['scale']['R']*x0
        yD = mission['scale']['R']*y0
        zD = mission['scale']['R']*z0
        VD = mission['scale']['V']*mission['launchSite']['nonDim']['V0']
        initialstate = np.array([xD,yD,zD,VD[0],VD[1],VD[2]])
        
        #finalHeight = mission['launchSite'].straightHeight
        #finaltime = 45.0
        #print 'FinalHeight',finalHeight
        rmag = np.sqrt(x0**2+y0**2+z0**2)
        u0 = np.array([x0,y0,z0])/rmag
        mission['launchSite']['launchDirection'] = u0
        uMat,etavec,timelist = get_uarray(u0,thetag,2.*omega)
        fux = interpolate.interp1d(timelist,uMat[:,0])
        fuy = interpolate.interp1d(timelist,uMat[:,1])
        fuz = interpolate.interp1d(timelist,uMat[:,2])

        planetmodel = 0


        #debrisvel = np.array([0,0,0])

        sref = mission['vehicle']['aero']['Aref'][0]
        minfcd = mission['vehicle']['aero']['MinfCD'][0]
        cd = mission['vehicle']['aero']['CD'][0]#1.0
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
        filename = 'debrisTrash10'
        dtinterval = .1
        #print 'Thrust' ,mission['vehicle']['propulsion']['Ftotal'][0]
        thrust = mission['vehicle']['propulsion']['Ftotal'][0]
        #print 'ratio',thrust/(mass*9.81) 

        thrustList = thrust*uMat
        ntime = len(timelist)
        isp = mission['vehicle']['propulsion']['ISP'][0]
        ndtinterval = 2000
        
        #print thrustList

        finalconditions,finalderivs = op.propagate(initialstate,mass,propCond,sref,minfcd,cd,cloption,minfcl,cl,loverd,atmosoption,altitudelist,densitylist,ulist,vlist,wlist,thrustList,timelist,isp,geoptions,filename,planetmodel,dtinterval,ndtinterval,thetag,propOption,ncd=len(minfcd),ncl=len(minfcl),ntime=ntime,nlist=len(altitudelist))
        
        newIndex = finalconditions[:,10]>=0.0
        newfinalConditions = finalconditions[newIndex,:]
   
        # updating new initial conditions
        
        mission['launchSite']['nonDim']['r0'] = newfinalConditions[-1,0:3]/mission['scale']['R']
        mission['launchSite']['nonDim']['V0'] = newfinalConditions[-1,3:6]/mission['scale']['V']
        #print mission['launchSite']['nonDim']['r0']
        #print mission['launchSite']['nonDim']['V0']
        thetag2 = thetag + omega*newfinalConditions[-1,7]
        mission['vehicle']['mass']['nonDim']['m0'][0] = newfinalConditions[-1,9]/mission['scale']['m']
        mfmin = mission['vehicle']['mass']['nonDim']['mf'][0]
        
        if (mfmin>=mission['vehicle']['mass']['nonDim']['m0'][0]):
            print 'Initial mass after propagation is less than possible minimum mass at the given stage. Try decreasing straight up altitude'
            exit()
        


        mission['solutionStage']['x'][0] = newfinalConditions[:,0]
        mission['solutionStage']['y'][0] = newfinalConditions[:,1]
        mission['solutionStage']['z'][0] = newfinalConditions[:,2]
        mission['solutionStage']['Vx'][0] = newfinalConditions[:,3]
        mission['solutionStage']['Vy'][0] = newfinalConditions[:,4]
        mission['solutionStage']['Vz'][0] = newfinalConditions[:,5]      
        mission['solutionStage']['m'][0] = newfinalConditions[:,9]
        mission['solutionStage']['t'][0] = newfinalConditions[:,10]
        mission['solutionStage']['ux'][0] = fux(mission['solutionStage']['t'][0])
        mission['solutionStage']['uy'][0] = fuy(mission['solutionStage']['t'][0])
        mission['solutionStage']['uz'][0] = fuz(mission['solutionStage']['t'][0])


        mission['solutionStage']['eta'][0] = np.ones(len(newfinalConditions[:,10]))

    elif (mission['launchSite']['orientation']=='vertical') and propagateSwitch==0:
        
        mission['solutionStage']['x'][0] = []
        mission['solutionStage']['y'][0] = []
        mission['solutionStage']['z'][0] = []
        mission['solutionStage']['Vx'][0] = []
        mission['solutionStage']['Vy'][0] = []
        mission['solutionStage']['Vz'] [0] = []      
        mission['solutionStage']['m'][0] = []
        mission['solutionStage']['t'][0] = []
        mission['solutionStage']['ux'][0] = []
        mission['solutionStage']['uy'][0] = []
        mission['solutionStage']['uz'][0] = []
        
        
    else:
        print 'ERROR: Specify launch type. Currently only vertical launch supported'
    mission['launchSite']['thetaLaunch'] = thetag
    mission['launchSite']['thetag'] = thetag2

    return mission


def get_uarray(u0,thetag,omega):

    maxtime = 500
    maxint = 10000
    u0 = np.array([u0]).T # making u a column vector
    timearray = np.linspace(0,maxtime,maxint) # harcoded, this is just temporary!!!
    etavec = np.ones((maxint))
    uMat = np.zeros((maxint,3))
    for index in range(len(timearray)):
        
        R = Rotatez(omega*timearray[index]).T #since we need the inertial frame coordinates
        uc = np.dot(R,u0)
        uMat[index,0] = uc[0,0]
        uMat[index,1] = uc[1,0]
        uMat[index,2] = uc[2,0]
    return (uMat,etavec,timearray)
    


def RefVernalEquinox(mission):
    if mission['planet']['name']=='EARTH':
        refYear = 2012
        refDay = 1
        refMonth = 1
        g0Hours = 6.6706801 # sidereal Hours (from Tom Colvin) data obtained from US Naval Observatory Website  
        thetag0 = g0Hours*2.*np.pi/24.#thetag0 referenced in BATES,MUELLER,White Fund of Astrodynamics pg 104
        revPerDay = 1.0027379093
    else:
        print 'ERROR in dateProcessing.py: Planet Ref Equinox not specified'
        exit(1)
    mission['VernalEquinox']['year'] = refYear
    mission['VernalEquinox']['day'] = refDay
    mission['VernalEquinox']['month'] = refMonth
    mission['VernalEquinox']['revPerDay'] = revPerDay
    mission['VernalEquinox']['thetag0'] = thetag0
    return mission

def getUTsec(localtime,timeShift):
    UTsec = localtime[0]*3600. + localtime[1]*60. + timeShift*3600
    return UTsec


#################################################################



def PropulsionData(mission):

    prop = mission['vehicle']['propulsion'] 

    # thrust, burn times, mass flow rates (powered phases)


    # dimensional quantities
    prop['Ftotal'] = prop['Nengines']*prop['F']                                     # N
    prop['mdotTotal'] = prop['Ftotal']/(prop['ISP']*9.8)                             # kg/s

    # non-dimensional quantities
    prop['nonDim']['mdot'] = prop['mdotTotal']*(mission['scale']['t']/mission['scale']['m'])      # non-dim
    prop['nonDim']['Fmg0'] = prop['Ftotal']/(mission['scale']['m']*mission['planet']['g0'])      # non-dim
    prop['nonDim']['tb'] = mission['vehicle']['mass']['nonDim']['mP']/prop['nonDim']['mdot']        # non-dim

    #print prop['nonDim'].mdot
    #print prop['nonDim'].Fmg0
    #print prop['nonDim'].tb
    mission['vehicle']['propulsion'] = prop 

    return mission


def OrbitData(mission):
# computes the eccentricity and angular momentum vectors to be used in forming the non-linear constraints
    import OrbitConversions as OC
    
    rp = mission['orbit']['hp'] + mission['planet']['R']
    e = mission['orbit']['e']
    a = rp/(1.-e)
    i = mission['orbit']['i']
    Om = mission['orbit']['Om']
    w = mission['orbit']['w']
    mu = mission['planet']['mu']
    h = mission['orbit']['h']
    v = mission['orbit']['vel']
    fpa = mission['orbit']['flightPathAngle']
    #print a,e,i,Om,w,mu
    
    if mission['options']['orbitConstraint']==1:
    
        hvec,evec = OC.OrbitalVectors(a,e,i,Om,w,0.0,mu)
        mission['orbit']['nonDim']['hvec'] = hvec/(mission['scale']['R']*mission['scale']['V'])
        mission['orbit']['nonDim']['evec'] = evec
    elif mission['options']['stateConstraint'] ==1:
        mission['orbit']['nonDim']['h'] = h/mission['scale']['R']
        mission['orbit']['nonDim']['vel'] = v/mission['scale']['V']
        mission['orbit']['nonDim']['flightPathAngle'] = fpa # in radians...['nonDim'] will be the same
    
    return mission


def aeroData(mission):
    from scipy import interpolate
    MinfCDList = mission['vehicle']['aero']['MinfCD']
    yList = mission['vehicle']['aero']['CD']
    flist = []
    for index in range(len(yList)):
        x = MinfCDList[index]
        y = yList[index]
        f = interpolate.interp1d(x,y,kind='linear',bounds_error=False,fill_value=y[-1]) # python doesn't seem to like complex values and interp1d
        flist.append(f)
    mission['vehicle']['aero']['CDinterp'] = flist
    return mission

def sph2cart(azimuth,elevation,r):
    # spherical 2 cartesian routine...same as Matlab's routine
    x = r*np.cos(elevation)*np.cos(azimuth)
    y = r*np.cos(elevation)*np.sin(azimuth)
    z = r*np.sin(elevation)
    return x,y,z


def Rotatez(theta):
    R = np.array([[np.cos(theta),np.sin(theta),0.0],[-np.sin(theta),np.cos(theta),0],[0,0,1.]])
    return R
















