import numpy as np
from scipy import interpolate 
from planetLibrary import CalculateAtmosphere,getSpeedOfSound
def AeroCalcFast(x,y,z,Vx,Vy,Vz,aero,stage,scale,planet,flag):
   #pdb.set_trace()
    A_ref = aero['Aref'][stage]
    #CDarray = aero.CD[stage]
    #MinfCD = aero.MinfCD[stage]
    fCD = aero['CDinterp'][stage]
    #print MinfCD
    #print CDarray
    #print A_ref

    Vx = Vx*scale['V']
    Vy = Vy*scale['V']
    Vz = Vz*scale['V']

    Vx_atm = -planet['omega']*y*scale['R'] # (m/s) Velocity of atmosphere
    Vy_atm = planet['omega']*x*scale['R']
    #Vz_atm = 0.*z;
    Vx_inf = Vx - Vx_atm # (m/s) Velocity relative to atmosphere
    Vy_inf = Vy - Vy_atm
    Vz_inf = Vz; #- Vz_atm;
    V_inf = (Vx_inf**2 + Vy_inf**2 + Vz_inf**2)**(0.5);


    alt = (x**2 + y**2 + z**2)**(0.5) - 1.; # non-dim altitude
    alt2 = alt*scale['R']; # altitude in meters
    c_sound = getSpeedOfSound(alt2,planet,flag) # this will return speed of sound array
    Minf = V_inf/c_sound
    
    Cd = fCD(Minf)

    atm_rho = CalculateAtmosphere(alt,planet);

    #q_inf = .5.*atm_rho.*V_inf.**2; % (Pa) Dynamic pressure

    #D = Cd.*q_inf.*A_ref;
    #DVinf = -D./(V_inf+epsilon);
    DVinf = (-.5*A_ref*Cd)*atm_rho*V_inf;

    Dx = DVinf*Vx_inf;
    Dy = DVinf*Vy_inf;
    Dz = DVinf*Vz_inf;

    return Dx, Dy, Dz

     
    # Dx = -D.*Vx_inf./V_inf;
    # Dy = -D.*Vy_inf./V_inf;
    # Dz = -D.*Vz_inf./V_inf;


            

