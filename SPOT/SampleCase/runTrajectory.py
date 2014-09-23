SPOT_folder = '../' # location for SPOT 
import sys
sys.path.append(SPOT_folder)


import copy_reg
import types

def reduce_method(m):
    return (getattr, (m.__self__, m.__func__.__name__))

copy_reg.pickle(types.MethodType, reduce_method)



def initialRun(mission):
    
    numIters = len(mission0['numerics']['iterationSchedule']['N'])
    for curIter in range(numIters):
        mission['numerics']['iterationSchedule']['curIter'] = curIter
        mission['numerics']['N'] = mission['numerics']['iterationSchedule']['N'][curIter]
        #mission = LDP.processAllMissionData(mission)
        info = -999
        infoOpt = -999
        if curIter ==0:
            includeDrag = True
            mission = LDP.processAllMissionData(mission)
            vars0 = IC.CreateInitialGuess(mission)  
        else:
            #mission.optimizer='snopt'
            if curIter==numIters-1:
                includeDrag = True
            else:
                includeDrag = True
            mission, vars0 = postProcessing.UpdateMissionAndInitialGuess(mission)
        
        if mission['optimizer'].lower() =='ipopt' and curIter==0:
            #print vars0
            x,info = optimSetup.IPOPTrun(vars0,mission,includeDrag)
            #x,info = optimSetup.IPOPTrunINIT(vars0,mission,includeDrag)
            infoOpt = info['status_msg'] 
            print 'status is ',info['status']
        elif mission['optimizer'].lower() =='ipopt':
            x,info = optimSetup.IPOPTrun(vars0,mission,includeDrag)
            infoOpt = info['status_msg'] 
            #mission.optimizer ='snopt'
            print 'status is ',info['status']
        elif mission['optimizer'].lower() =='snopt' and curIter==0:
            #exitval = optimSetup.SNOPTrunINIT(vars0,mission,includeDrag)
            exitval = optimSetup.SNOPTrun(vars0,mission,includeDrag)
            f_opt = exitval[0]
            x = exitval[1]
            #print 'cal to read ',exitval[2]['text']
            infoOpt = exitval[2]['text']
            
            print infoOpt
        #print 'objective function = ',f_opt
        
        elif mission['optimizer'].lower() =='snopt':
            exitval = optimSetup.SNOPTrun(vars0,mission,includeDrag)
            f_opt = exitval[0]
            x = exitval[1]
            #print 'cal to read ',exitval[2]['text']
            infoOpt = exitval[2]['text']
            print infoOpt
        #print 'objective function = ',f_opt
        else:
            print 'Error: Optimizer '+mission['optimizer']+' not found'
        if (curIter==0):
            if (mission['optimizer'].lower() =='ipopt' and info['status'] !=0) or (mission['optimizer'].lower() =='snopt' and infoOpt!='finished successfully'):
	        x = IC.CreateInitialGuess(mission)# previous solution is useless..better to start fresh
        mission = postProcessing.RecordSolution(x,mission)
    return mission,x,info,infoOpt







import LaunchDataProcessing as LDP
import readInputTrajectory as RI
import optimSetup
import initialConditions as IC
import planetLibrary as PL
import numpy as np
import postProcessing
import os
import data2GoogleEarthSPOT as GE
import shelve
from scipy.io import savemat
import cPickle
import copy
import sys

fileName = sys.argv[-1]#'input_SPOT.dat'

fileWarningName = 'optWarnings.txt'
fileWarning = open(fileWarningName,'w')
fileWarning.write('')
# Load the mission data

mission0 = RI.readInput(fileName)# reading input file
mission0 = PL.planetParameters(mission0) # getting planet parameters from library


Noriginal = 1*mission0['numerics']['N']
Nnew = np.zeros((len(Noriginal)))+7 # 10 nodes per stage for any initial guess
N14 = (1.0/2.0*(Noriginal+Nnew))
N14= np.floor(N14)

mission0['numerics']['N'] = Nnew
mission0 = LDP.processAllMissionData(mission0)
np.set_printoptions(precision=15)

theta0 = mission0['launchSite']['thetaLaunch']
omega = mission0['planet']['omega']
#print 'theta0',theta0
#print theta0
# First run -> 5 nodes per stage + NO DRAG
#numIters = 4

mission0['numerics']['iterationSchedule']['N'] = np.array([Nnew,N14,Noriginal])

mission0['numerics']['iterationSchedule']['N'][mission0['numerics']['iterationSchedule']['N']<5] = 5
numIters = len(mission0['numerics']['iterationSchedule']['N'])
refmin = 1.0*mission0['launchSite']['localTime'][1]

if mission0['multiTraj']['ntraj']>1:

    minvec = np.linspace(refmin - mission0['multiTraj']['dt'],refmin + mission0['multiTraj']['dt'],mission0['multiTraj']['ntraj'])
else: 
    minvec = np.array([refmin])

#missionList = int(mission0.multiTraj.ntraj)*[mission0]
lenMiss =  int(mission0['multiTraj']['ntraj'])

missionList = [mission0 for i in range(lenMiss)]


nrun = 0
for minvecindex in range(len(minvec)):
    fileNameGE = 'GEperStage'+str(nrun)+'.kml'
    mission = missionList[minvecindex]
    mission['launchSite']['localTime'] = [mission0['launchSite']['localTime'][0],minvec[minvecindex]]
    #mission.launchSite.thetaLaunch = 1.0*thetavec[thetaindex]
    #mission = LDP.processAllMissionData(mission,flagthetag=True)
    mission = LDP.LaunchSiteData(mission) # initial test passed
    #print mission.launchSite.localTime

    
    #exit()
    
    if nrun==0:
        mission,x,info,infoOpt = initialRun(mission)
        nrun = nrun + 1
    elif nrun>0:
        #mission = LDP.processAllMissionData(mission,flagthetag=True)

        curIter = 2#numIters - 1 
        for curIter in range(2,len(mission0['numerics']['iterationSchedule']['N'])):
            mission['numerics']['iterationSchedule']['curIter'] = curIter
            includeDrag = True
            mission, vars0 = postProcessing.UpdateMissionAndInitialGuess(mission)
            if mission['optimizer'].lower() =='ipopt':
                x,info = optimSetup.IPOPTrun(vars0,mission,includeDrag)
                infoOpt = info['status_msg'] 
                print infoOpt
            elif mission['optimizer'].lower() =='snopt':
                flagVal = 0
                if curIter==len(mission0['numerics']['iterationSchedule']['N'])-1:
                    flagVal = 1
                
                exitval = optimSetup.SNOPTrun(vars0,mission,includeDrag,flagFORCE=flagVal)
                f_opt = exitval[0]
                x = exitval[1]
                #print exitval[2]['text']
                #print 'objective function = ',f_opt
                infoOpt = exitval[2]['text']
                print infoOpt
            else:
                print 'Error: Optimizer '+mission['optimizer']+' not found'
            mission = postProcessing.RecordSolution(x,mission)
            nrun = nrun + 1
    '''
    # comment the 3 lines below if the above commented is desired
    mission,x,info,infoOpt = initialRun(mission)
    mission = postProcessing.RecordSolution(x,mission)
    nrun = 0    
    '''

    if nrun >0   : # recording potential error in the optimization process
        if (mission['optimizer'].lower() =='ipopt' and info['status'] !=0) or (mission['optimizer'].lower() =='snopt' and infoOpt!='finished successfully'):
            
            fileWarning.write('Check run ' + str(nrun) +'\n')
            fileWarning.write(infoOpt + '\n')
            fileWarning.write('----------------------------------\n')
            print 'Rerunning'
            
            mission,x,info,infoOpt = initialRun(mission)
            if (mission['optimizer'].lower() =='ipopt' and info['status'] ==0) or (mission['optimizer'].lower() =='snopt' and infoOpt=='finished successfully'):
                
                fileWarning.write('run ' + str(nrun) +' seems to be fine after rerunning \n')
                fileWarning.write(infoOpt + '\n')
                fileWarning.write('----------------------------------\n')



    mission = postProcessing.RecordSolutionPerStage(x,mission)
    #mission = postProcessing.convert2LatLon(mission)
    mission = postProcessing.convert2LatLonPerStage(mission)
    missionList[minvecindex] = copy.deepcopy(mission)
    #print 'thetag is ',mission.launchSite.thetaLaunch

    #postProcessing.plotting(x,mission)
    #GE.convert(mission)
    GE.convertPerStage(mission,fileNameGE) # generating KML file for Google Earth
# using cpickle to save the results to read later in python
fileWarning.close()
#savemat('./missionSolSnopt.mat', dict(mission=mission,x_opt=x))
#savemat('./missionSolMulti.mat',dict(mission=mission))
#savemat('./missionList.mat',dict(missionList=missionList))
postProcessing.writeTextFiles(missionList)

#print missionList
#print missionList[0]


pickleFileNameList = 'solList.txt'
missionListPickle = missionList
pickleFile = open(pickleFileNameList,'w')
cPickle.dump(missionListPickle,pickleFile)
pickleFile.close()


#pickleFileName = 'solution.txt'
#missionPickle = [mission['solution']]
#pickleFile = open(pickleFileName,'w')
#cPickle.dump(missionPickle,pickleFile)
#pickleFile.close()/Users/fcapristan/Desktop/Capristan2.jpg

