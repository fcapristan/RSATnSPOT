# this set of routines have the initial guess, and the interpolated guess (once a previous solution has been calculated)

def CreateInitialGuess(mission):
    import numpy as np
    import variables

#    curIter = mission.numerics.iterationSchedule.curIter 
#    if curIter == 0

    Nstate = mission['numerics']['Nstate']
    Ncontrol = mission['numerics']['Ncontrol']
    
    # generate lousy guess
    r0 = mission['launchSite']['nonDim']['r0']
    V0 = mission['launchSite']['nonDim']['V0']
    onesVector = np.ones(Nstate)
    
    x = np.linspace(r0[0],1.001*r0[0],Nstate)
    y = np.linspace(r0[1],1.001*r0[1],Nstate)
    z = np.linspace(r0[2],1.001*r0[2],Nstate)
    Vx = np.linspace(V0[0],1.001*V0[0],Nstate)
    Vy = np.linspace(V0[1],1.001*V0[1],Nstate)
    Vz = np.linspace(V0[2],1.001*V0[2],Nstate)
    
    
    #Vx = V0[0]*onesVector
    #Vy = V0[1]*onesVector
    #Vz = V0[2]*onesVector
    
    onesVector = np.ones(Ncontrol)
    zerosVector = np.zeros((Ncontrol))
    norm_r0 = np.linalg.norm(r0)
    u0 = r0/norm_r0
    
    theta1_0 = np.arccos(u0[2])
    theta2_0 = np.arctan(u0[1]/u0[0])
    if theta2_0<0:
        theta2_0 = theta2_0 + 2.0*np.pi
    '''
    ux = u0[0]*onesVector
    uy = u0[1]*onesVector
    uz = u0[2]*onesVector*0
    '''
    theta1  = theta1_0*onesVector
    theta2 = theta2_0*onesVector
    
    print 'thetas',theta1,theta2
    eta = onesVector
    #stages = mission.vehicle.stages
    #tfTotal = np.sum(mission.vehicle.propulsion.nonDim.tb)
    
    # mdot = mission.vehicle.propulsion.mdotTotal
    #mpropellant = mission.vehicle.mass.propellant
    tf = mission['vehicle']['propulsion']['nonDim']['tb']    
    #tf = np.ones(stages)*tfTotal/stages
    vars0 = variables.assembleVariables(x,y,z,Vx,Vy,Vz,theta1,theta2,eta,tf)
    
   
    return vars0


