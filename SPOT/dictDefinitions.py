#****************************************************************************************
# File: dictDefinitions.py
# 
# Defines the dictionaries to be used by python
# 
# Created by:          Francisco C.
# Created on:          06/07/2012
# $Rev: 43 $  
# $LastChangedDate: 2011-06-02 16:15:15 -0700 (Thu, 02 Jun 2011) $ 
# $LastChangedBy: fcapristan $     
#
#*****************************************************************************************		

def createMain():
    mission = {'vehicle':createVehicle(),
                'launchSite':createLaunchSite(),
                'orbit':createOrbit(),
                'nonDim':createNonDimMass(),
                'scale':createScale(),
                'numerics':createNumerics(),
                'planet':createPlanet(),
                'VernalEquinox':createVernalEquinox(),
                'solution':createSolution(),
                'solutionStage':createSolutionStage(),
                'optimizer':0,
                'multiTraj':createMultiTraj(),
                'tol':5e-7,
                'options':createOptions()}
    return mission

def createMultiTraj():
    multiTraj = {'ntraj':1,'dt':0}
    return multiTraj

def createOptions():
    options = {'propagateStage' :[-9999], # negative by default
                'orbitConstraint' : 1,
                'stateConstraint' : 0,
                'objective' : 'mf'} # final mass is the default
    return options
def createSolution():
    solution = {'x' : 0,
                'y' : 0,
                'z' : 0,
                'Vx' : 0,
                'Vy' : 0,
                'Vz' : 0,
                'ux' : 0,
                'uy' : 0,
                'uz' : 0,
                'eta' : 0,    
                'tState' : 0,
                'tControl' : 0,
                'm' : 0,
                'mf' : 0,
                'tf' : 0,
                'margin' : 0,
                'N' : 0,
                'latitude' : 0,
                'longitude' : 0,
                'height' : 0}
    return solution


def createSolutionStage():
    solutionStage={'thetag' : 'NaN',
                    'x' : 0,
                    'y' : 0,
                    'z' : 0,
                    'Vx' : 0,
                    'Vy' : 0,
                    'Vz' : 0,
                    'ux' : 0,
                    'uy' : 0,
                    'uz' : 0,
                    'eta': 0,   
                    't' : 0,
                    'm' : 0,
                    'latitude' : 0,
                    'longitude' : 0,
                    'height' : 0,
                    'propagateLast' : createPropagateLast()}
    return solutionStage
def createPropagateLast():
    propagateLast = {'thetag' : 'NaN',
                    'x' : 0,
                    'y' : 0,
                    'z' : 0,
                    'Vx' : 0,
                    'Vy' : 0,
                    'Vz' : 0,
                    'ux' : 0,
                    'uy' : 0,
                    'uz' : 0,
                    'eta' : 0,    
                    't': 0,
                    'm' : 0,
                    'latitude' : 0,
                    'longitude' : 0,
                    'height' : 0}
    return propagateLast
def createVernalEquinox():
    VernalEquinox = {'thetag0' : 0,
                    'thetag' : 0,
                    'year' : 0,
                    'day' : 0,
                    'month' : 0,
                    'revPerDay' : 0}

    return VernalEquinox
def createVehicle():
    vehicle = {'stages' : 0,
                'propulsion' : createPropulsion(),
                'name' : 'none',
                'mass' : createVehicleMass(),
                'aero' : createAero()}
    return vehicle
def createVehicleMass():
    mass = {'propellant' : 0,
            'structural' : 0,
            'payload' : 0,
            'm0' : 0,
            'mf' : 0,
            'nonDim' : createNonDimMass()}
    return mass
def createAero():
    aero = {'MinfCD' : 0,
            'MinfCL' : 0,
            'CD' : 0,
            'CL' : 0,
            'Aref' : 0,
            'CDinterp' : 0}

    return aero
        

def createLaunchSite():
    launchSite = {'h' : 0, #height [m]
                    'lat' : 0, # [deg]
                    'long' : 0, # [deg]
                    'localTime' : 0, 
                    'shiftUT' : 0,
                    'month' : 0,
                    'day' : 0,
                    'year' : 0,
                    'thetag' : 0,
                    'orientation' : 0,
                    'nonDim' : createNonDimLaunchSite(),
                    'launchDirecion' : 0,
                    'straightHeight' : 0,
                    'initialTraj' : 0,
                    'thetaLaunch' : 'NaN'}
    return launchSite

def createOrbit():
    orbit = {'hp' : 0, # altitude of periapsis
                'e' : 0, # eccentricity
                'i' : 0, # inclination [deg]
                'Om' : 0, # ascending node [deg]
                'w' : 0, # argument of periapsis [deg]
                # parameters below only used if stateConstraints is set to 1
                'h' : 0,# regular altitude desired [m]
                'vel' : 0, # velocity magnitude desired [m/sec]
                'flightPathAngle' : 0, # flight path angle desired [rad]
                'nonDim' : createNonDimOrbit()}
    return orbit

def createPlanet():
    planet = {'mu' : 0, # planet's grav parameter
                'g0' : 0, # ref gravity
                'omega' : 0, # rotation rate
                'name' : 'none',
                'HS' : 0, #planet scale height for exponential density 
                'R' : 0, # planet's radius
                'T' : 0, # sidereal day
                'm' : 0} # mass [kg]

    return planet
def createPropulsion():
    propulsion = {'ISP' : 0,
                    'F' : 0,
                    'Ftotal' : 0,
                    'Nengines' : 0,
                    'minThrottle' : 0,
                    'mdotTotal': 0,
                    'nonDim' : createNonDimPropulsion()}   
    return propulsion

def createNonDimMass():
    nonDimMass = {'m0' : 0, # initial mass
                    'mf' : 0, # final mass
                    'mS' : 0, # structural mass
                    'mP' : 0, # propellant mass
                    'mPL' : 0, # payload mass
                    'MR' : 0} # mass ratio
    return nonDimMass
def createNonDimLaunchSite():
    nonDimLaunch = {'r' : 0,
                    'r0' : 0,
                    'V0' : 0,
                    'r0_ref' : 0,
                    'V0_ref' : 0}
    return nonDimLaunch

def createNonDimPropulsion():
    nonDimPropulsion = {'mdot' : 0,
                        'Fmg0': 0,
                        'tb' : 0}
    return nonDimPropulsion
def createNonDimOrbit():
    nonDimOrbit = {'hvec' : 0,
                    'evec' : 0,
                    'h' : 0, # regular altitude desired
                    'vel' : 0, # velocity magnitude desired
                    'flightPathAngle' : 0} # flight path angle desired [rad]
    return nonDimOrbit

def createScale():
    scale = {'m': 0,
                'R' : 0,
                't' : 0,
                'V' : 0}
    return scale

def createNumerics():
    numerics = {'range' : createRange(),
                'iterationSchedule' : createIterationSchedule(),
                'Nstate' : 0,
                'Ncontrol' : 0,
                'nvars' : 0,
                'dthalf' : 0,
                'N' : 0,
                'dmfdetaNONDIM' : 0} 
    return numerics
def createRange():
    range = {'x' : 0,
                'y' : 0,
                'z' : 0,
                'Vx' : 0,
                'Vy' : 0,
                'Vz' : 0,
                'ux' : 0,
                'uy' : 0,
                'uz' : 0,
                'eta' : 0,
                'tf' : 0}
    return range

def createIterationSchedule():
    iterSch = {'numIterations' : 0,
                'N' : 0,
                'dragOn' : 0}
    return iterSch
