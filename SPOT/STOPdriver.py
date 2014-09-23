import readInputTrajectory as RI
import LaunchDataProcessing as LDP
import planetLibrary as PL
import initialConditions as IC
import constraints
import costFunction
import scipy.optimize as SO
import numpy as np
import pyOpt
from scipy.io import savemat
import time
import variables
import pdb
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
fileName = 'input_STOP.dat'
import postProcessing

# Define the functions the optimizer will call
def objectiveFunction(inVars,mission):
    x=costFunction.ObjectiveMass(inVars,mission)
    eq=constraints.equality(inVars,mission)
    ineq=constraints.inequality(inVars,mission)

def objectiveFunction(inVars,mission):
    x=costFunction.ObjectiveMass(inVars,mission)
    eq=constraints.equality(inVars,mission)
    ineq=constraints.inequality(inVars,mission)
    
    g = np.concatenate((eq,ineq),1)
    fail = 0
    return x,g,fail

def sensitivityFunction(inVars,f,g,mission):
    x = costFunction.fprimeObjectiveMass(inVars,mission)
    eq = constraints.fprimeequality(inVars,mission)
    ineq = constraints.fprimeinequality(inVars,mission)
    
    g = np.concatenate((eq,ineq),0)
    fail = 0
    return x,g,fail


# Load the mission data
mission = RI.readInput(fileName)# reading input file
mission = PL.planetParameters(mission) # getting planet parameters from library

# Simulate loading in an iteration schedule from input.dat
mission['numerics']['iterationSchedule']['numIterations'] = 2
mission['numerics']['iterationSchedule']['N'] = np.array([[5, 5],[ 10, 10]])

numIters = mission['numerics']['iterationSchedule']['numIterations']

for curIter in range(numIters):
    # Create the initial guess for the current optimization problem
    mission['numerics']['iterationSchedule']['curIter'] = curIter

    if curIter == 0:
        mission = LDP.processAllMissionData(mission)
        vars0 = IC.CreateInitialGuess(mission)
    else:
        mission, vars0 = postProcessing.UpdateMissionAndInitialGuess(mission)
        #vars0 = IC.UpdateInitialGuess(mission)

    numEquality = len(constraints.equality(vars0,mission))
    numInequality = len(constraints.inequality(vars0,mission))

    # Find the upper and lower bounds
    #boundsCase = constraints.bounds(mission)
    lb,ub = constraints.getXLXU(vars0,mission)


    #TJC
    opt_prob = pyOpt.Optimization('Trajectory Optimization',lambda x: objectiveFunction(x,mission))
    opt_prob.addObj('Objective Mass')

    # Specify all of the variables
    print 'Setting up variables in a hackish way.  MUST CHANGE!!!'
    for curVar in range(len(vars0)):
        opt_prob.addVar('var'+str(curVar), 'c', value=vars0[curVar], lower=lb[curVar], upper=ub[curVar])

    # Now add in equality constraints
    for curCon in range(numEquality):
        opt_prob.addCon('g' + str(curCon), 'e')

    # Now add in inequality constraints
    for curCon in range(numEquality,numEquality + numInequality):
        opt_prob.addCon('g' + str(curCon), 'i')

    # Confirm that everything is correct
    print opt_prob

    # Set up the optimizer
    snopt = pyOpt.pySNOPT.SNOPT()
    snopt.setOption('Major feasibility tolerance',value=1e-4)
    snopt.setOption('Major optimality tolerance',value=1e-4)
    snopt.setOption('Minor feasibility tolerance',value=1e-4)

    # Optimize and save results
    startTime = time.time()
    sens2 = lambda x,f,g:sensitivityFunction(x,f,g,mission)
    #exitVal = snopt(opt_prob,sens_type= sens2)
    exitVal = snopt(opt_prob,sens_type= 'FD')

    elapsedTime = time.time() - startTime
    f_opt = exitVal[0]
    x_opt = exitVal[1]
    print exitVal[2]['text']
    print elapsedTime
    savemat('./missionSolSnopt.mat', dict(mission=mission,f_opt=f_opt,x_opt=x_opt))
    mission = postProcessing.RecordSolution(x_opt,mission)






# added some simple plots
makePlot = False

if makePlot:
    x,y,z=variables.decomposeLocationVariables(x_opt,mission['numerics']['range'])

    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.plot(x,y,z)
    u, v = np.mgrid[0:2*np.pi:200j, 0:np.pi:200j]
    x1=np.cos(u)*np.sin(v)
    y1=np.sin(u)*np.sin(v)
    z1=np.cos(v)
    ax.plot_surface(x1, y1, z1, color="r")

    plt.show()




