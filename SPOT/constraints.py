import numpy as np
import variables 
import scipy.integrate as SI
from AeroCalc import AeroCalcFast
import cmath

def getXLXU(vars,mission):
    r0 = mission['launchSite']['nonDim']['r0']

    XL = np.zeros(len(vars)) - 100.
    XU = np.zeros(len(vars)) + 100.

    
    XL[mission['numerics']['range']['theta1']] = -np.inf
    XL[mission['numerics']['range']['theta2']] = -np.inf


    
    XU[mission['numerics']['range']['theta1']] = np.inf
    XU[mission['numerics']['range']['theta2']] = np.inf
    numberofStages = mission['vehicle']['stages']
    etaListDown = []
    for stage in range(0,numberofStages):
        rangeControl = mission['numerics']['range']['stageControl'][stage]
        etaboundDown = len(rangeControl)*[mission['vehicle']['propulsion']['minThrottle'][stage]]
        etaListDown = etaListDown + etaboundDown
    XL[mission['numerics']['range']['eta']] = etaListDown
    XU[mission['numerics']['range']['eta']] = 1.
    N = mission['numerics']['N'][0]
    if N<=6:
        XL[mission['numerics']['range']['tf']] = .5/mission['scale']['t']
    else:
        XL[mission['numerics']['range']['tf']] = 5/mission['scale']['t']
    XU[mission['numerics']['range']['tf']] = np.inf
    #XU[mission['numerics']['range']['tf'][1]] = .5/mission['scale']['t']

    return (XL,XU)
    
        
def eval_g(vars,mission,includeDrag=True):
    ceq = equality(vars,mission,'real',includeDrag) # must be = 0
    cieq = inequality(vars,mission,'real') # must be >=0, 
    
    return np.concatenate((ceq,cieq),1)



def eval_jac_g(vars,mission,flag='2d',includeDrag=True):

    # this function should only be used with CYIPOPT
    dceqdx = fprimeequality(vars,mission,flag,includeDrag)
    dcineqdx = fprimeinequality(vars,mission,flag)
    ret = np.concatenate((dceqdx,dcineqdx),1)

    return ret
    




def equality(vars,mission,flag='real',includeDrag=True):
    
##### First SIMPLE equality constraints
    etavec = mission['vehicle']['propulsion']['minThrottle']
    nvars = mission['numerics']['nvars']
    numberofStages = mission['vehicle']['stages']
    x,y,z,Vx,Vy,Vz,theta1,theta2,eta,tf = variables.decomposeVariables(vars,mission['numerics']['range'])
    
    u0 = mission['launchSite']['launchDirection']
    r0 = mission['launchSite']['nonDim']['r0']
    V0 = mission['launchSite']['nonDim']['V0']
    
    theta1_0 = np.arccos(u0[2])
    theta2_0 = np.arctan(u0[1]/u0[0])
    
    if mission['options']['objective']=='mf':
        nLen = numberofStages-1
    elif mission['options']['objective']=='none':
        nLen = numberofStages
    
    if flag=='complex': # recasting to complex
        linearVals = np.zeros(8,dtype=complex)
        x = np.array(x,dtype=complex)
        y = np.array(y,dtype=complex)
        z = np.array(z,dtype=complex)
        Vx = np.array(Vx,dtype=complex)
        Vy = np.array(Vy,dtype=complex)
        Vz = np.array(Vz,dtype=complex)
        theta1 = np.array(theta1,dtype=complex)
        theta2 = np.array(theta2,dtype=complex)
        eta = np.array(eta,dtype=complex)
        tf = np.array(tf,dtype=complex)
        d = np.zeros(nLen,dtype=complex)

    elif flag=='real':
        #print theta1,theta2

        d = np.zeros(nLen)
        linearVals = np.zeros(8)

    # launch direction
    linearVals[0] = 0#theta1_0 - theta1[0]
    linearVals[1] = 0#theta2_0 - theta2[0] 
    # launch location
    linearVals[2] = vars[mission['numerics']['range']['x'][0]] - r0[0]
    linearVals[3] = vars[mission['numerics']['range']['y'][0]] - r0[1]
    linearVals[4] = vars[mission['numerics']['range']['z'][0]] - r0[2]

    # Initial Velocity
    linearVals[5] = vars[mission['numerics']['range']['Vx'][0]] - V0[0]
    linearVals[6] = vars[mission['numerics']['range']['Vy'][0]] - V0[1]
    linearVals[7] = vars[mission['numerics']['range']['Vz'][0]] - V0[2]




    ceq_ODEs,mT = equationsOfMotion(x,y,z,Vx,Vy,Vz,theta1,theta2,eta,tf,mission,flag,includeDrag)
    #ceq_umag = 1. - (ux**2.+ uy**2.+uz**2.)**.5
    AngMom = mission['orbit']['nonDim']['hvec']
    evec = mission['orbit']['nonDim']['evec']
    Rvec = np.array([x[-1],y[len(y)-1],z[-1]])
    Rmag = (np.sum(Rvec**2))**.5
    Vvec = np.array([Vx[-1],Vy[-1],Vz[-1]])
    
    if mission['options']['orbitConstraint']==1 and mission['options']['stateConstraint']==1:
        print 'Orbital Constraint and State Constraint cannot be both activated'
        exit()
    
    if mission['options']['orbitConstraint']==1:
       
        h = np.cross(Rvec,Vvec)
        e = np.cross(Vvec,h) - Rvec/Rmag
        angdiff = (AngMom - h)
        evecdiff = (evec - e)
        #print np.array([angdiff,evecdiff]).T
        ceq_orbit = np.concatenate((angdiff,evecdiff),1)
    elif mission['options']['stateConstraint']==1:
        
        Vx_atm = -mission['planet']['omega']*y[-1]; # (m/s) Velocity of atmosphere
        Vy_atm = mission['planet']['omega']*x[-1];
        #Vz_atm = 0.*z;
        Vx_inf = Vx[-1] - Vx_atm; # (m/s) Velocity relative to atmosphere
        Vy_inf = Vy[-1] - Vy_atm;
        Vz_inf = Vz[-1]; #- Vz_atm;
        V_inf = np.array([Vx_inf,Vy_inf,Vz_inf])
            
        
        Rf = mission['orbit']['nonDim']['h'] + 1.
        fpaF = mission['orbit']['nonDim']['flightPathAngle']
        altitude2 = x[-1]**2 + y[-1]**2 + z[-1]**2 
        hfdiff = altitude2-Rf**2
       
        zdiff = []# z[-1] - r0[2]
        #fpa = np.pi - np.arccos(np.dot(Rvec,Vvec)/(np.linalg.norm(Rvec)*np.linalg.norm(Vvec)))
        #fpa = np.pi - cmath.acos(np.dot(Rvec,Vvec)/(np.linalg.norm(Rvec)*np.linalg.norm(Vvec)))
        if flag=='complex':
            fpa = .5*np.pi - cmath.acos(mydot(Rvec,V_inf)/(mynorm(Rvec)*mynorm(V_inf)))
        elif flag=='real':
            fpa = .5*np.pi - np.arccos(mydot(Rvec,V_inf)/(mynorm(Rvec)*mynorm(V_inf)))

        diffFPA = fpa - fpaF
        ceq_orbit = np.concatenate(([],[diffFPA]),1)

        #ceq_orbit = np.concatenate(([hfdiff],[zdiff]),1)
        #ceq_orbit = np.concatenate(([hfdiff],[]),1)

    if numberofStages>1:
        for stage in range(0,nLen):
            rangeControl = mission['numerics']['range']['stageControl'][stage]
            if len(rangeControl)>0:
                mu = mT[rangeControl[len(rangeControl)-1]]
                mf = mission['vehicle']['mass']['nonDim']['mf'][stage]
                d[stage] = (mf-mu)*10

    #print mu,mf,rangeControl
    #print np.shape(linearVals),np.shape(ceq_orbit)
    #ceq = np.concatenate((ceq_orbit,linearVals,ceq_ODEs,ceq_umag,d),1)
    ceq = np.concatenate((ceq_orbit,linearVals,ceq_ODEs,d),1)

    #ceq = np.concatenate((ceq_ODEs,ceq_umag,ceq_orbit,d),1)
    #ceq = np.concatenate((ceq_ODEs,ceq_orbit,d),1) # NOT INCLUDING UNIT VECTOR CONSTRAINT

    return ceq

def mynorm(x):
    return (x[0]**2 + x[1]**2 + x[2]**2)**.5
def mydot(x,y):
    return (x[0]*y[0] + x[1]*y[1] + x[2]*y[2])


def equationsOfMotion(xT,yT,zT,VxT,VyT,VzT,theta1T,theta2T,etaT,tfT,mission,flag='real',includeDrag=True):
    numberofStages = mission['vehicle']['stages']
    etavec = mission['vehicle']['propulsion']['minThrottle']
    if flag=='complex':
        RT = np.zeros((6*np.sum(mission['numerics']['N'])),dtype=complex)
    elif flag=='real':
        RT = np.zeros((6*np.sum(mission['numerics']['N'])))
    mT = computeNumericMass(etaT,tfT,mission,flag)
    #mT[mT<0] = 0.00001
    mT_m1 =  1./mT
    rangeEQ = mission['numerics']['range']['EQ']
    #if flag=='real':
        #print mT
    for stage in range(0,numberofStages):
        if etavec[stage]>=0:
            N = mission['numerics']['N'][stage]
            rangeState = mission['numerics']['range']['stageState'][stage]
            rangeControl = mission['numerics']['range']['stageControl'][stage]
            x = xT[rangeState]
            y = yT[rangeState]
            z = zT[rangeState]
            Vx = VxT[rangeState]
            Vy = VyT[rangeState]
            Vz = VzT[rangeState]
            theta1 = theta1T[rangeControl]
            theta2 = theta2T[rangeControl]
            eta = etaT[rangeControl]
            tf = tfT[stage]
            
            ux = np.sin(theta1)*np.cos(theta2)
            uy = np.sin(theta1)*np.sin(theta2)
            uz = np.cos(theta1)
            if tf==0.: # to avoid divide by zero errors
                tf = .000001/mission['scale']['t']
            
            ddt = mission['numerics']['D'][stage]/tf      
            F = mission['vehicle']['propulsion']['nonDim']['Fmg0'][stage]

            m_m1 = mT_m1[rangeControl]
            c = F*eta*m_m1
            rm3 = (x**2.+y**2.+z**2.)**(-1.5)



            #print np.dot(ddt,Vz),'ww'
            #x=RT[rangeEQ[stage]]   
            if includeDrag:
     
                DX,DY,DZ = AeroCalcFast(x,y,z,Vx,Vy,Vz,mission['vehicle']['aero'],stage,mission['scale'],mission['planet'],flag)
            
            #Non-dimensionalize the drag and find the acceleration
                dx = DX*m_m1/(mission['scale']['m']*mission['planet']['g0'])
                dy = DY*m_m1/(mission['scale']['m']*mission['planet']['g0'])
                dz = DZ*m_m1/(mission['scale']['m']*mission['planet']['g0'])
                
                
                
                RT[rangeEQ[stage]]= np.concatenate((Vx - np.dot(ddt,x),\
                                   Vy-np.dot(ddt,y),\
                                   Vz-np.dot(ddt,z),\
                                   c*ux + dx - x*rm3 - np.dot(ddt,Vx),\
                                   c*uy + dy - y*rm3 - np.dot(ddt,Vy),\
                                   c*uz + dz - z*rm3 - np.dot(ddt,Vz)))
            else:
                RT[rangeEQ[stage]]= np.concatenate((Vx - np.dot(ddt,x),\
                                                    Vy-np.dot(ddt,y),\
                                                    Vz-np.dot(ddt,z),\
                                                    c*ux - x*rm3 - np.dot(ddt,Vx),\
                                                    c*uy - y*rm3 - np.dot(ddt,Vy),\
                                                    c*uz - z*rm3 - np.dot(ddt,Vz)))
            #ddt = mission['numerics'].D[stage]/tf      
           #F = mission['vehicle']['propulsion']['nonDim'].Fmg0[stage]
                   
            #m_m1 = mT_m1[rangeControl]
            #c = F*eta*m_m1
            #rm3 = (x**2.+y**2.+z**2.)**(-1.5)

            #print np.dot(ddt,Vz),'ww'
            #x=RT[rangeEQ[stage]]   
            ''' 
            RT[rangeEQ[stage]]= np.concatenate((Vx - np.dot(ddt,x),\
                                                Vy-np.dot(ddt,y),\
                                                Vz-np.dot(ddt,z),\
                                                c*ux - x*rm3 - np.dot(ddt,Vx),\
                                                c*uy - y*rm3 - np.dot(ddt,Vy),\
                                                c*uz - z*rm3 - np.dot(ddt,Vz)))
                
            '''
        else:
            N = mission['numerics']['N'][stage]
            rangeState = mission['numerics']['range']['stageState'][stage]
            rangeControl = mission['numerics']['range']['stageControl'][stage]
            x = xT[rangeState]
            y = yT[rangeState]
            z = zT[rangeState]
            Vx = VxT[rangeState]
            Vy = VyT[rangeState]
            Vz = VzT[rangeState]
            tf = tfT[stage]
            #if tf<=0.: # to avoid divide by zero errors
                #tf = .0001/mission['scale']['t']
            
            ddt = mission['numerics']['D'][stage]/tf      
            
            m_m10= 1./mission['vehicle']['mass']['nonDim']['m0'][stage]
            #m_m1 = mT_m1[rangeControl]
            rm3 = (x**2.+y**2.+z**2.)**(-1.5)
            
            
            
            #print np.dot(ddt,Vz),'ww'
            #x=RT[rangeEQ[stage]]   
            if includeDrag:
                
                DX,DY,DZ = AeroCalcFast(x,y,z,Vx,Vy,Vz,mission['vehicle']['aero'],stage,mission['scale'],mission['planet'],flag)
                
                #Non-dimensionalize the drag and find the acceleration
                dx = DX*m_m10/(mission['scale']['m']*mission['planet']['g0'])
                dy = DY*m_m10/(mission['scale']['m']*mission['planet']['g0'])
                dz = DZ*m_m10/(mission['scale']['m']*mission['planet']['g0'])
                
                
                
                RT[rangeEQ[stage]]= np.concatenate((Vx - np.dot(ddt,x),\
                                                    Vy-np.dot(ddt,y),\
                                                    Vz-np.dot(ddt,z),\
                                                    dx - x*rm3 - np.dot(ddt,Vx),\
                                                    dy - y*rm3 - np.dot(ddt,Vy),\
                                                    dz - z*rm3 - np.dot(ddt,Vz)))
            else:
                RT[rangeEQ[stage]]= np.concatenate((Vx - np.dot(ddt,x),\
                                                    Vy-np.dot(ddt,y),\
                                                    Vz-np.dot(ddt,z),\
                                                    - x*rm3 - np.dot(ddt,Vx),\
                                                    - y*rm3 - np.dot(ddt,Vy),\
                                                    - z*rm3 - np.dot(ddt,Vz)))
    return (RT,mT)



def computeNumericMass(etaTotal,tf,mission,flag):
    
    Ncontrol = mission['numerics']['Ncontrol']
    numberofStages = mission['vehicle']['stages']
    etavec = mission['vehicle']['propulsion']['minThrottle']
    if flag=='complex':
        mret = np.zeros((Ncontrol),dtype=complex)
    elif flag=='real':
        mret = np.zeros((Ncontrol))
    for stage in range(numberofStages):
        if etavec[stage]>=0:
            rangeControl = mission['numerics']['range']['stageControl'][stage]
            eta = etaTotal[rangeControl]
            #dthalf = tf[stage]*mission['numerics'].dthalf[stage]
            #etasum = eta + np.concatenate(([0.0],eta[0:len(eta)-1]),1)
            t = tf[stage]*mission['numerics']['t'][stage]
            #inteta = np.cumsum(np.multiply(etasum,dthalf))
            
            m0 = mission['vehicle']['mass']['nonDim']['m0'][stage]
            inteta = np.concatenate(([0],SI.cumtrapz(eta,t)),1)
            mdot = mission['vehicle']['propulsion']['nonDim']['mdot'][stage]
            m = m0-mdot*inteta
            mret[rangeControl] = m
            #mMin = .001*mission.vehicle.mass['nonDim'].mPL
            #mret[mret<mMin] = mMin
    
    return (mret)




def inequality(vars,mission,flag='real'):
    x,y,z = variables.decomposeLocationVariables(vars,mission['numerics']['range'])
    #x,y,z,Vx,Vy,Vz,ux,uy,uz,eta,tf = variables.decomposeVariables(vars,mission['numerics'].range)
    c = -(x**2+y**2+z**2)+1. #c(x)<=0
    
    #m=computeNumericMass(eta,tf,mission)
    
    
    #ret = np.concatenate((-c,[-m[-1]]),1)
    return c


def fprimeequality(vars,mission,flag='2d',includeDrag=True):
    ret=equality(vars,mission,'real',includeDrag=False)
    ncons = len(ret)
    nvars = len(vars) # this can be avoided by doing a calculation before hand....to be fixed
    h = 1e-40
    x0 = np.array(vars,dtype=complex)
    h1i = complex(0,h)
    dceqdx = np.zeros((ncons,nvars))
    for n in range(0,nvars):
        x = 1*x0
        
        x[n] = np.complex(vars[n],h)
        #x[n] = vars[n]+h
        ceq = equality(x,mission,'complex',includeDrag)
        #ceq = equality(x,mission)
        #if n ==0:
        #    print ceq.imag/h
        dceqdx[:,n] = ceq.imag/h
    if flag == '1d':
        dceqdx = np.reshape(dceqdx,ncons*nvars)
    #   dceqdx = np.concatenate((dceqdx,ceq.imag/h))
    return dceqdx



def fprimeinequality(vars,mission,flag='2d'):
    ret=inequality(vars,mission,'real')
    ncons = len(ret)
    nvars = len(vars)
    h = 1e-40
    x0 = np.array(vars,dtype=complex)
    h1i = complex(0,h)
    dceqdx = np.zeros((ncons,nvars))
    for n in range(0,nvars):
        x = 1.*x0
        
        x[n] = np.complex(vars[n],h)
        
        ceq = inequality(x,mission,'complex')
        #if n ==0:
        #    print ceq.imag/h
        dceqdx[:,n] = ceq.imag/h
    if flag=='1d':
        dceqdx = np.reshape(dceqdx,ncons*nvars)

    return dceqdx


def bounds(vars,mission):
    # only called once so additions here won't affect the running time dramatically
    
    ceq = equality(vars,mission,'real') # must be = 0
    cieq = inequality(vars,mission,'real') # must be >=0,
    gu = np.zeros((len(ceq)+len(cieq)))
    gl1 = np.zeros((len(ceq)))
    gl2 = np.zeros((len(cieq))) - np.inf
    gl = np.concatenate((gl1,gl2),1)
    '''
        gvals=eval_g(vars,mission)
        nvars = mission['numerics'].nvars
        numberofStages = mission['vehicle']['stages']
        x,y,z,Vx,Vy,Vz,ux,uy,uz,eta,tf = variables.decomposeVariables(vars,mission['numerics'].range)
        neta = len(eta)
        gleta = np.ones((neta))
        gueta = np.ones((neta))
        gl = np.zeros((len(gvals)))
        gu = np.zeros((len(gvals)))
        ineq = inequality(vars,mission)
        ineqIndex = np.arange(len(ineq)) + len(gvals)  - len(ineq)
        gu[ineqIndex] = 1e15
        '''
    return (gl,gu)

