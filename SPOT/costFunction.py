
import numpy as np

def ObjectiveMass(vars,mission):
    eta = vars[mission['numerics']['range']['eta']]
    tf = vars[mission['numerics']['range']['tf']]
    m = computeFinalNumericMass(eta,tf,mission)
    ret = -m 
    #ret = 5.0
    #print m
    #print 'Cost function is ',m,tf
    return (ret)# minimizing this value

def ObjectiveMassINIT(vars,mission):
    import variables 
    x,y,z,Vx,Vy,Vz,ux,uy,uz,eta,tf = variables.decomposeVariables(vars,mission['numerics']['range'])

    #eta = vars[mission.numerics.range.eta]
    tf = vars[mission['numerics']['range']['tf']][-1]
    #m = computeFinalNumericMass(eta,tf,mission)
    #ret = -(x[-1]**2 + y[-1]**2 + z[-1]**2)**.5
    ret = tf
    #ret = tf
    #ret = 5.0
    #print m
    #print 'Cost function is ',m,tf
    return (ret)# minimizing this value

def computeFinalNumericMass(etaTotal,tfTotal,mission):
    
    Ncontrol = mission['numerics']['Ncontrol']
    numberofStages = mission['vehicle']['stages']
    
    stage = numberofStages-1
    
    rangeControl = mission['numerics']['range']['stageControl'][stage]
    eta = etaTotal[rangeControl]
    t = tfTotal[stage]*mission['numerics']['t'][stage]
    #dthalf = tfTotal[stage]*mission.numerics.dthalf[stage]
    #etasum = eta + np.concatenate(([0.0],eta[0:len(eta)-1]),1)
    #inteta = np.sum(np.multiply(etasum,dthalf))
    
    inteta = np.trapz(eta,t)
    m0 = mission['vehicle']['mass']['nonDim']['m0'][stage]
    mdot = mission['vehicle']['propulsion']['nonDim']['mdot'][stage]
    m = m0-mdot*inteta
    
    return (m)


def fprimeObjectiveMass(vars,mission): # careful when using this derivative. ONLY VALID FOR TRAPEZOIDAL INTEGRATION
    h = 1.e-100
    tf = vars[mission['numerics']['range']['tf'][-1]]
    eta1= vars[mission['numerics']['range']['eta']]
    rangeControl = mission['numerics']['range']['stageControl'][-1]
    eta = eta1[rangeControl]
    mdot = mission['vehicle']['propulsion']['nonDim']['mdot'][-1]
    dmfdeta = tf*mdot*mission['numerics']['dmfdetaNondim']
    t = np.complex(tf,h)*mission['numerics']['t'][-1]
    m0 = mission['vehicle']['mass']['nonDim']['m0'][-1]

    x = integrateFinalMass(t,eta,mdot,m0)
    dmfdtf = x.imag/h
    dmdx = np.zeros((1,mission['numerics']['nvars']))
    dmdx[0,mission['numerics']['range']['eta'][rangeControl]] = dmfdeta
    dmdx[0,mission['numerics']['range']['tf'][-1]] = dmfdtf
    dmdx = -dmdx

    return dmdx

def fprimeObjectiveMassINIT(vars,mission): # for no cost function

    ncons = 1
    nvars = len(vars)
    h = 1e-100
    x0 = np.array(vars,dtype=complex)
    h1i = complex(0,h)
    dceqdx = np.zeros((ncons,nvars))
    for n in range(0,nvars):
        x = 1.*x0
        
        x[n] = np.complex(vars[n],h)
        
        ceq = ObjectiveMassINIT(x,mission)
        #if n ==0:
        #    print ceq.imag/h
        dceqdx[:,n] = ceq.imag/h
    # dceqdx = np.reshape(dceqdx,ncons*nvars)
    
    return dceqdx

def integrateFinalMass(t,eta,mdot,m0):
    inteta = np.trapz(eta,t)
    mfinal = m0-mdot*inteta
    return mfinal
'''

def fprimeObjectiveMass(vars,mission): # VALID FOR ALL...TIME CONSUMING
    ncons = 1
    nvars = len(vars)
    h = 1e-100
    x0 = np.array(vars,dtype=complex)
    h1i = complex(0,h)
    dceqdx = np.zeros((ncons,nvars))
    for n in range(0,nvars):
        x = 1.*x0
        
        x[n] = np.complex(vars[n],h)
        
        ceq = ObjectiveMass(x,mission)
        #if n ==0:
        #    print ceq.imag/h
        dceqdx[:,n] = ceq.imag/h
   # dceqdx = np.reshape(dceqdx,ncons*nvars)

    return dceqdx
'''
