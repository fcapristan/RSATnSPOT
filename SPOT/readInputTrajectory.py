
#****************************************************************************************
# File: ReadInput.py
# 
# Reads the input file
# Based on a sample input reader routine provided by Thomas Economon
# 
# Created by:          Francisco C.
# Created on:          04/09/2011
# $Rev: 60 $  
# $LastChangedDate: 2012-07-06 15:34:41 -0700 (Fri, 06 Jul 2012) $ 
# $LastChangedBy: fcapristan $     
#
#*****************************************************************************************		



import numpy as np
import dictDefinitions
import readAeroData  as rA

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
    mission = dictDefinitions.createMain() # defining dictionaries

    try:
           inputFile = open(fileName,'r')
    except:
           if len(sys.argv) == 1:
                   print '\n!!! Error: No input file specified !!!\n' \
                   + ' Proper command line usage: $ readInputTrajectory.py [inputFileName] \n \n'
           else:
                   print '\n!!! Error: Could not open input file: ' + sys.argv[1] + ' !!!\n'
           exit(1)

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

                   if key[0] == 'planet' and key[1] == 'name':
                           PlanetName = value[0]
                   elif key[0] == 'planet' and key[1] == 'mu':
                           PlanetMu = float(value[0])
                           if PlanetMu <= 0:
                              print '\n!!! Error: Invalid Gravitational Parameter mu, it must be positive \n'\
                              + ' At line: ' + line.strip() + '\n'
                              exit(1)
                   elif key[0] == 'planet' and key[1] == 'omega':
                           PlanetOmega = float(value[0])
                   elif key[0] == 'planet' and key[1] == 'hs':
                           PlanetHs = float(value[0])
                           if PlanetHs <= 0:
                              print '\n!!! Error: Invalid Exponential Model Height, it must be positive \n'\
                              + ' At line: ' + line.strip() + '\n'
                              exit(1)
                   elif key[0] == 'planet' and key[1] == 'r':
                      PlanetR = float(value[0])
                      if PlanetR <= 0:
                          print '\n!!! Error: Invalid Radius, it must be positive \n'\
                          + ' At line: ' + line.strip() + '\n'
                          exit(1)
        
                                  
                                  

                  # Vehicle information section  
                   elif key[0] == 'vehicle' and key[1] == 'name':
                           VehicleName = value[0]
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
                   elif key[0] == 'vehicle' and key[1] == 'max' and key[2] == 'thrust':
                           checkpos = 1 # making sure values are postive
                           checkint = 0 # if integer check is needed

                           VehicleMaxThrust,err= makeVector(value,VehicleStages,checkpos,checkint)
                           if err!=0:
                              print '\n!!! Error: Invalid  Max Thrust \n'\
                              + ' At line: ' + line.strip() + '\n'
                              exit(1)
                   elif key[0] == 'vehicle' and key[1] == 'num' and key[2] == 'engines':
                           checkpos = 1 # making sure values are postive
                           checkint = 1# if integer check is needed
                           VehicleNumEngines,err= makeVector(value,VehicleStages,checkpos,checkint)
                           if err!=0:
                              print '\n!!! Error: Invalid Number of Engines \n'\
                              + ' At line: ' + line.strip() + '\n'
                              exit(1)
                   elif key[0] == 'vehicle' and key[1] == 'min' and key[2] == 'throttle':
                           checkpos = 1 # making sure values are postive
                           checkint = 0# if integer check is needed
                           VehicleMinThrottle,err= makeVector(value,VehicleStages,checkpos,checkint)
                           if err!=0:
                              print '\n!!! Using coasting period. Min Throttle <0 \n'\
                              + ' At line: ' + line.strip() + '\n'
                              ''' 
                              print '\n!!! Error: Invalid Minimum Throttle \n'\
                              + ' At line: ' + line.strip() + '\n'
                              exit(1)
                              '''
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
                   elif key[0] == 'elevation' and key[1] == 'initial':
                           ElevationInitial = float(value[0])
                           if ElevationInitial < 0:
                              print '\n!!! Error: Invalid Initial Elevation, it must be positive \n'\
                              + ' At line: ' + line.strip() + '\n'
                              exit(1)
                   elif key[0] == 'launch' and key[1] == 'orientation':
                           LaunchOrientation = value[0]
        
                   #constraints to determine what orbit/state is desired at the end of optimization
                               
                               
                   elif key[0] == 'orbit' and key[1] == 'constraints':
                       mission['options']['orbitConstraint'] = int(value[0])
                   elif key[0] == 'state' and key[1] == 'constraints':
                       mission['options']['orbitConstraint'] = 0
                       mission['options']['stateConstraint'] = int(value[0])
                   elif key[0] == 'hp' and key[1] == 'final':
                           hpFinal = float(value[0])
                           if hpFinal <= 0:
                              print '\n!!! Error: Invalid altitude of periapsis, it must be greater than zero \n'\
                              + ' At line: ' + line.strip() + '\n'
                              exit(1)
                   elif key[0] == 'e' and key[1] == 'final':
                           eFinal = float(value[0])
                   elif key[0] == 'inc' and key[1] == 'final':
                           incFinal = float(value[0])
                   elif key[0] == 'komega' and key[1] == 'final':
                           komegaFinal = float(value[0])
                   elif key[0] == 'omega' and key[1] == 'final':
                           omegaFinal = float(value[0])
                   elif key[0] == 'nu' and key[1] == 'final':
                           nuFinal = float(value[0])

                   # conditions if a desired final state is required...currently considering sounding rockts
                   elif key[0] == 'altitude' and key[1] == 'final':
                        mission['orbit']['h'] = float(value[0])
                   elif key[0] == 'velocity' and key[1] == 'final':
                        mission['orbit']['v'] = float(value[0])
                   elif key[0] == 'flight' and key[1] == 'path' and key[2]=='final':
                        mission['orbit']['flightPathAngle'] = float(value[0])*np.pi/180. # converting to radians
                                       

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
                                  
                    

                   elif key[0] =='launch' and key[1]=='date':
                        dateval = value[0].split('/')
                        month = float(dateval[0])
                        day = float(dateval[1])
                        year = float(dateval[2])

                        if (isint(month)==0 or isint(day)==0 or isint(year)==0):
                            print '\n!!! Error: Invalid date, it must be positive integers \n'\
                            + ' At line: ' + line.strip() + '\n'
                            exit(1)                         
                   elif key[0].lower() == 'straight' and key[1].lower() =='height':
                           standardHeight = float(value[0])
                           mission['launchSite']['straightHeight'] = standardHeight
                   elif key[0].lower() == 'straight' and key[1].lower() =='time':
                           standardTime = float(value[0])
                           mission['launchSite']['straightTime'] = standardTime
                        
                   elif key[0] == 'n' and key[1] == 'points':
                           checkpos = 1 # making sure values are postive
                           checkint = 1# if integer check is needed
                           Npoints,err= makeVector(value,VehicleStages,checkpos,checkint)
                           if err!=0:
                              print '\n!!! Error: Invalid Number of Points \n'\
                              + ' At line: ' + line.strip() + '\n'
                              exit(1)
                   elif key[0] == 'optimizer':
                        optimizer = value[0]
                   elif key[0] =='deltatlaunch':
                        dtLaunch = float(value[0])
                   elif key[0].lower() =='ntrajectories':
                        ntraj = float(value[0])

                   elif key[0] == 'objective' and key[1] == 'function':
                        mission['options']['objective'] = value[0]
        
                   elif key[0] == 'propagate' and key[1] == 'stage':
                       checkpos = 1 # making sure values are postive
                       checkint = 1# if integer check is needed
                       propStage,err= makeVectorSimple(value,checkpos,checkint)
                       if err!=0:
                           print '\n!!! Error: Invalid Number of Stage to Propagate \n'\
                           + ' At line: ' + line.strip() + '\n'
                           exit(1)
                       mission['options']['propagateStage']=propStage
                   else:
                       print '\n!!! Error: Unrecognized input parameter !!!\n' \
                       + ' At line: ' + line.strip() + '\n'
                       exit(1)    
            
           elif len(line.strip()) != 0:
                   print '\n!!! Error: Unrecognized input parameter !!!!\n' \
                   + ' At line: ' + line.strip() + '\n'
                
                   exit(1)

    inputFile.close()
 
    mission['optimizer'] = optimizer
# setting planet parameters
    mission['planet']['name'] = PlanetName 

# setting vehicle information
    mission['vehicle']['name'] = VehicleName
    mission['vehicle']['stages'] = np.array(VehicleStages)
    mission['vehicle']['mass']['structural'] = np.array(VehicleStructuralMass)
    mission['vehicle']['mass']['propellant'] = np.array(VehiclePropellantMass)
    mission['vehicle']['mass']['payload'] = np.array(VehiclePayloadMass)

    mission['vehicle']['propulsion']['ISP'] = np.array(VehicleIsp) 
    mission['vehicle']['propulsion']['F'] = np.array(VehicleMaxThrust)
    mission['vehicle']['propulsion']['Nengines'] = np.array(VehicleNumEngines)
    mission['vehicle']['propulsion']['minThrottle'] = np.array(VehicleMinThrottle)
    mission['vehicle']['aero']['MinfCD'] = MinfCDlist
    mission['vehicle']['aero']['CD'] = CDList
    mission['vehicle']['aero']['Aref'] = ArefList
    
# launch site  information
    mission['launchSite']['lat'] = np.array(LatitudeInitial)
    mission['launchSite']['long'] = np.array(LongitudeInitial)
    mission['launchSite']['h'] = np.array(ElevationInitial)
    mission['launchSite']['orientation'] = np.array(LaunchOrientation)
    mission['launchSite']['localTime'] = np.array(localTime)
    mission['launchSite']['shiftUT'] = np.array(UTshift)
    mission['launchSite']['month'] = np.array(month)
    mission['launchSite']['day'] = np.array(day)
    mission['launchSite']['year'] = np.array(year)
    #print standardHeight,ElevationInitial
    #mission['launchSite']['straightHeight = ElevationInitial + standardHeight
                       
    mission['orbit']['hp'] = np.array(hpFinal)
    mission['orbit']['e'] = np.array(eFinal)    
    mission['orbit']['i'] = np.array(incFinal)*np.pi/180.
    mission['orbit']['Om'] = np.array(komegaFinal)*np.pi/180.
    mission['orbit']['w'] = np.array(omegaFinal)*np.pi/180.
#   nuFinal currently ignored
    #print 'checking', mission['numerics']
    mission['numerics']['N'] = Npoints

    mission['multiTraj']['ntraj'] = ntraj
    mission['multiTraj']['dt'] = dtLaunch #  given in minutes
    return(mission)

    




