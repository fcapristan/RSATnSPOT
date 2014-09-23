
#****************************************************************************************
# File: propagate.py
# 
# Reads the input file
# Based on a sample input reader routine provided by Thomas Economon
# 
# Created by:          Francisco C.

#
#*****************************************************************************************		



import numpy as np
#import dictDefinitions
import readAeroData  as rA
import readThrust as rT
planetModel =0





def RefVernalEquinox():

    refYear = 2012
    refDay = 1
    refMonth = 1
    g0Hours = 6.6706801 # sidereal Hours (from Tom Colvin) data obtained from US Naval Observatory Website  
    thetag0 = g0Hours*2.*np.pi/24.#thetag0 referenced in BATES,MUELLER,White Fund of Astrodynamics pg 104
    revPerDay = 1.0027379093
    dateListRef = [refMonth,refDay,refYear]

    return thetag0,revPerDay,dateListRef

def getUTsec(localtime,timeShift):
    UTsec = localtime[0]*3600. + localtime[1]*60. + timeShift*3600
    return UTsec


def calcThetag(dateList):
    from datetime import datetime
    thetag0,revPerDay,dateListRef = RefVernalEquinox()
    [refMonth,refDay,refYear] = dateListRef
    [month,day,year,UTshift,localTime] = dateList
    a = datetime(int(year),int(month),int(day))
    b = datetime(refYear,refMonth,refDay)
    days = (a-b).days
    timeSEC = getUTsec(localTime,UTshift)
    D = days + timeSEC/(3600.0*24.0) # converting to fractions of DAYS ...ref Fund of ASTRO. Bates
    thetag = thetag0 + revPerDay*2.*np.pi*D
    return thetag


# function to check if number is an integer
# returns 1 if it is an integer, 0 otherwise
def isint(number):
    x = number % 1
    b = 0
    if x == 0:
        b = 1
    return b
# end of integer function

def makeVector(array,stages,checkPos,checkint):
    err = 0
    if len(array)!=int(stages):
        print ' \n  Error: check number of stages agreement'
        err = 1
    ret = []
    for index in range(0,int(stages)):
        val = float(array[index])
        if (checkPos==1 and val<0):
            print ' \n Warning: value must be positive'
            err = 2
        if ((checkint==1) and (isint(val)==0)):
            err = 3
        ret.append(float(array[index]))
    return (np.array(ret),err)
def makeVectorSimple(array,checkPos,checkint):
    err = 0
    
    ret = []
    for index in range(0,len(array)):
        val = float(array[index])
        if (checkPos==1 and val<0):
            print ' \n Warning: value must be positive'
            err = 2
        if ((checkint==1) and (isint(val)==0)):
            err = 3
        ret.append(float(array[index]))
    return (np.array(ret),err)     


def readInput(fileName):
    #mission = dictDefinitions.createMain() # defining dictionaries
    
    try:
        inputFile = open(fileName,'r')
    except:
        if len(sys.argv) == 1:
            print '\n!!! Error: No input file specified !!!\n' \
            + ' Proper command line usage: $ readInputTrajectory.py [inputFileName] \n \n'
        else:
            print '\n!!! Error: Could not open input file: ' + sys.argv[1] + ' !!!\n'
        exit(1)
    thetag0 = 'NaN'
    for line in inputFile:
        if line[0] == '#' or line.find("#") > 0:
            iComment = line.index('#')
            line = line[:iComment]
        if line.find(":") > 0:
            iColon = line.index(":")
            key = (line[:iColon].lower()).split()
            value = (line[iColon+1:]).split()
            if len(value) < 1:
                print '\n!!! Error: Invalid input arguments !!!\n' \
                + ' At line: ' + line.strip() + '\n'
                exit(1)
            

            # Vehicle information section  

            elif key[0] == 'vehicle' and key[1] == 'stages':
                VehicleStages = float(value[0])
                if isint(VehicleStages) == 0 or VehicleStages <= 0:
                    print '\n!!! Error: Invalid number of stages, it must be a positive integer \n'\
                        + ' At line: ' + line.strip() + '\n'
                    exit(1)
                VehicleStages = int(VehicleStages)
            elif key[0] == 'vehicle' and key[1] == 'payload' and key[2]=='mass':
                VehiclePayloadMass = float(value[0])
                if VehiclePayloadMass <= 0:
                    print '\n!!! Error: Invalid Payload Mass \n'\
                        + ' At line: ' + line.strip() + '\n'
                    exit(1)
            
            
            
            elif key[0] == 'vehicle' and key[1] == 'structural' and key[2] == 'mass':
                checkpos = 1 # making sure values are postive.
                checkint = 0# if integer check is needed
                VehicleStructuralMass,err= makeVector(value,VehicleStages,checkpos,checkint)
                if err!=0:
                    print '\n!!! Error: Invalid Structural Mass \n'\
                    + ' At line: ' + line.strip() + '\n'
                    exit(1)
            elif key[0] == 'vehicle' and key[1] == 'propellant' and key[2] == 'mass':
                checkpos = 1 # making sure values are postive
                checkint = 0# if integer check is needed
                
                VehiclePropellantMass,err= makeVector(value,VehicleStages,checkpos,checkint)
                if err!=0:
                    print '\n!!! Error: Invalid Propellant Mass \n'\
                    + ' At line: ' + line.strip() + '\n'
                    exit(1)
            
            
            
            
            elif key[0] == 'vehicle' and key[1] == 'isp':
                checkpos = 1 # making sure values are postive
                checkint = 0# if integer check is needed
                
                VehicleIsp,err= makeVector(value,VehicleStages,checkpos,checkint)
                if err!=0:
                    print '\n!!! Error: Invalid ISP \n'\
                    + ' At line: ' + line.strip() + '\n'
                    exit(1)



            elif key[0] == 'vehicle' and key[1] == 'cd':
                #reading CD text files
                MinfCDlist = []
                CDList = []
                ArefList = []
                if len(value)!= VehicleStages:
                    print '\n!!! ERROR in readInputTrajectory.py. Incorrect number of CD files'
                    exit(1)
                for index in range(len(value)):
                    #print value[index]
                    Minf,CD,Aref = rA.readInputCD(value[index])
                    MinfCDlist.append(Minf)
                    CDList.append(CD)
                    ArefList.append(Aref)
            
            
            elif key[0] == 'latitude' and key[1] == 'initial':
                LatitudeInitial = float(value[0])
                if LatitudeInitial < -90 or LatitudeInitial > 90:
                    print '\n!!! Error: Latitude must be within -90 to 90 degrees \n'\
                    + ' At line: ' + line.strip() + '\n'
                    exit(1)
            elif key[0] == 'longitude' and key[1] == 'initial':
                LongitudeInitial = float(value[0])
                if LongitudeInitial < -180 or LongitudeInitial > 180:
                    print '\n!!! Error: Longitude must be within -180 to 180 degrees \n'\
                    + ' At line: ' + line.strip() + '\n'
                    exit(1)
            elif key[0] == 'altitude' and key[1] == 'initial':
                ElevationInitial = float(value[0])
                if ElevationInitial < 0:
                    print '\n!!! Error: Invalid Initial Elevation, it must be positive \n'\
                    + ' At line: ' + line.strip() + '\n'
                    exit(1)

            elif key[0] == 'thetag' and key[1] == 'initial':
                thetag0 = float(value[0])


            elif key[0] == 'thrust' and key[1] == 'file':
                #reading CD text files
                timelist = []
                thrustList = []
                stageList = []
                if len(value)!= VehicleStages:
                    print '\n!!! ERROR in readInputTrajectory.py. Incorrect number of CD files'
                    exit(1)
                for index in range(len(value)):
                    #print value[index]
                    stageVal,timeVal,Tx,Ty,Tz = rT.readInputThrust(value[index])
                    timelist.append(timeVal)
                    thrustList.append([Tx,Ty,Tz])
                    stageList.append(stageVal)
                        
            
            elif key[0]=='atmospheric' and key[1]=='option':
                atmoOption = float(value[0])
            elif key[0]=='atmospheric' and key[1]=='file':
                atmoFile = value[0]
            elif key[0] == 'launch' and key[1] == 'time':
                checkpos = 1 # making sure values are postive
                checkint = 1# if integer check is needed
                timeLength = 2 # hours and mins
                localTime,err= makeVector(value,timeLength,checkpos,checkint)
                if err!=0:
                    print '\n!!! Error: Invalid Local Time \n'\
                    + ' At line: ' + line.strip() + '\n'
                    exit(1)
            elif key[0] == 'launch' and key[1] == 'ut':
                UTshift = float(value[0])
            elif key[0] =='dt':
                dtval = float(value[0])

            elif key[0] =='launch' and key[1]=='date':
                dateval = value[0].split('/')
                month = float(dateval[0])
                day = float(dateval[1])
                year = float(dateval[2])
                
                if (isint(month)==0 or isint(day)==0 or isint(year)==0):
                    print '\n!!! Error: Invalid date, it must be positive integers \n'\
                    + ' At line: ' + line.strip() + '\n'
                    exit(1)                                 

        elif len(line.strip()) != 0:
            print '\n!!! Error: Unrecognized input parameter !!!!\n' \
            + ' At line: ' + line.strip() + '\n'
            
            exit(1)
    
    inputFile.close()
                
    
    dateList = [month,day,year,UTshift,localTime]
    initialLocation = [LongitudeInitial,LatitudeInitial,ElevationInitial]
    # stage masses
    nend = VehicleStages
    m0 = np.zeros((VehicleStages))
    mf = np.zeros((VehicleStages))
    for index in range(0,VehicleStages):
        m0[index] = np.sum(VehiclePropellantMass[index:nend]) + np.sum(VehicleStructuralMass[index:nend]) + VehiclePayloadMass         # kg
        mf[index] = np.sum(VehiclePropellantMass[index+1:nend]) + np.sum(VehicleStructuralMass[index:nend]) + VehiclePayloadMass      # kg
        print VehiclePropellantMass[index+1:nend]
    massList = [m0,mf]
 
    if thetag0=='NaN':
        thetag0 = calcThetag(dateList)
    else:
        print 'Using specified thetag in input file. Date is ignored.'
    outputs = [initialLocation,massList,ArefList,MinfCDlist,CDList,atmoOption,atmoFile,timelist,thrustList,VehicleIsp,thetag0,dateList,dtval]
        
    '''finalconditions,finalderivs = propagate(initialstate,mass,mass_time_alt_final,sref,minfcd,cd,cloption,minfcl,cl,loverd,atmosoption,altitudelist,densitylist,ulist,vlist,wlist,tlist,timelist,isp,geoptions,filename,planetmodel,dtinterval,ndtinterval,thetag,mass_time_alt_opt,thrustoffangledeg,ncd=len(minfcd),ncl=len(minfcl),ntime=shape(tlist,0),nlist=len(altitudelist))
    '''
                
    return(outputs)



def getStateVector(fileName):
    import orbitTools
    from scipy.interpolate import UnivariateSpline
    from scipy import interpolate
    import orbitProp as op
    import sys
    sys.path.append('../PythonScripts')
    import AtmosProfile as AP
    import matplotlib.pyplot as plt
    planetModel = 0

    omegaE = 2.*np.pi/(86164.0906)
    inputs = readInput(fileName)
    [lon0,lat0,h0] = inputs[0]
    thetag = inputs[10]
    ispList = inputs[9]
    
    [x0,y0,z0] = orbitTools.latlonalt2ECI(lat0,lon0,h0,thetag,planetModel)
    #SEZ2ECI(lat,lon,hs,Vsouth,Veast,Vzenith,inertial,thetag,planetModel)
    V0 = orbitTools.SEZ2ECI(lat0,lon0,h0,0,0,0,0,thetag,planetModel)
    Vx0 = V0[0]
    Vy0 = V0[1]
    Vz0 = V0[2]
    massList = inputs[1]
    m0Vec = massList[0]
    mfVec = massList[1]
    SrefList = inputs[2]
    MinfCDlist = inputs[3]
    CDList = inputs[4]
    atmosoption = inputs[5]
    atmoFile = inputs[6]
    timelist = inputs[7]
    thrustList = inputs[8]
    retTrajs = []
    #atmosoption = 2
    if atmosoption>0:
        print 'using given atmo data'
        altitudeList,densityMeanList,uMeanList,vMeanList,wMeanList,densitySDList,uSDList,vSDList,wSDList,nPoints= AP.readGramAtmos(atmoFile)
    else:
        
        altitudeList = [0]
        densityMeanList = [0]
        uMeanList = [0]
        vMeanList = [0]
        wMeanList = [0]
    
    densitylist = densityMeanList
    ulist = uMeanList
    vlist = vMeanList
    wlist = wMeanList
    dtinterval = inputs[12]
    
    cloption = 0
    cl = [0]
    loverd = 0
    minfcl = [1]
    geoptions = 1
    fileGE = 'trajProp'
    thrustoffangledeg =[0,0]
    mass_time_alt_opt = 1
    dtinterval = 0.01
    stateVectors = []
    for index in range(0,len(m0Vec)):
        mass = m0Vec[index]
        mf = mfVec[index]
        initialstate = [x0,y0,z0,Vx0,Vy0,Vz0]
        thrust = thrustList[index]
        Tx = thrust[0]
        Ty = thrust[1]
        Tz = thrust[2]

        timeVec = timelist[index]
        '''
        fTx = UnivariateSpline(timeVec,Tx)
        fTy = UnivariateSpline(timeVec,Ty)
        fTz = UnivariateSpline(timeVec,Tz)

        timeVec = np.linspace(timeVec[0],timeVec[-1],1000)
        Tx = fTx(timeVec)
        Ty = fTy(timeVec)
        Tz = fTz(timeVec)
        '''    
        cd0  = CDList[index]
        minfcd0 = MinfCDlist[index]
        sref = SrefList[index]
        fcd0 = UnivariateSpline(minfcd0,cd0)
        
        isp = ispList[index]
        
     
        minfcd = np.linspace(minfcd0[0],minfcd0[-1],100)
        cd = fcd0(minfcd)
        Tlist = np.array([Tx,Ty,Tz]).T#np.concatenate([[Tx],[Ty],[Tz]],axis=1)
      
        
        ndtinterval = int(np.ceil(2*timeVec[-1]/dtinterval))
        ntime = len(timeVec)
        #finalconditions,derivVals = op.propagate(initialstate,mass,mf,sref,minfcd,cd,cloption,minfcl,cl,loverd,atmosoption,altitudeList,densitylist,ulist,vlist,wlist,Tlist,timeVec,isp,geoptions,fileGE,planetModel,dtinterval,ndtinterval,thetag,1,ncd=len(minfcd),ncl=len(minfcl),ntime=ntime,nlist=len(altitudeList))
        
        #finalconditions,finalderivs = propagate(initialstate,mass,mass_time_alt_final,sref,minfcd,cd,cloption,minfcl,cl,loverd,atmosoption,altitudelist,densitylist,ulist,vlist,wlist,tlist,timelist,isp,geoptions,filename,planetmodel,dtinterval,ndtinterval,thetag,mass_time_alt_opt,thrustoffangledeg,ncd=len(minfcd),ncl=len(minfcl),ntime=shape(tlist,0),nlist=len(altitudelist))


        if index == len(m0Vec)-1:
            mass_time_alt_opt = 2
            fcond = timeVec[-1]
        else:
            fcond = mf
        finalconditions,finalderivs = op.propagate(initialstate,mass,fcond,sref,minfcd,cd,cloption,minfcl,cl,loverd,atmosoption,altitudeList,densitylist,ulist,vlist,wlist,Tlist,timeVec,isp,geoptions,fileGE+str(1+index),planetModel,dtinterval,ndtinterval,thetag,mass_time_alt_opt)
        
        newIndex = finalconditions[:,10]>=0.0
        newfinalConditions = finalconditions[newIndex,:]
        stateVectors.append(newfinalConditions)


        x0 = newfinalConditions[-1,0]
        y0 = newfinalConditions[-1,1]
        z0 = newfinalConditions[-1,2]
        Vx0 = newfinalConditions[-1,3]
        Vy0 = newfinalConditions[-1,4]
        Vz0 = newfinalConditions[-1,5]
        massCheck = newfinalConditions[-1,9]
        timeCheck = newfinalConditions[-1,7]
 
        thetag = thetag + omegaE*newfinalConditions[-1,10]
        retTrajs.append(newfinalConditions)

import data2GE as GE

getStateVector('nominalParam.txt')
GE.convert('trajProp',2)
