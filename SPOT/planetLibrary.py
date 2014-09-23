#   planetLibrary.py: Script that returns required planet parameters given
#                    desired planet
#
#   Francisco Capristan, ADL
#
#   Started:        5/02/12
#
#
#   Inputs:     mission             mission data structure
#
#
#   Outputs:    mission             output data structure containing collected input data,
#                                   derived data, and solution data 

#################################################################
import numpy as np

def planetParameters(mission):

    planetName = mission['planet']['name']
    planet = mission['planet']
    if planetName == 'MARS':
        planet['m'] = .64185e24 # (kg) Mass of Mars
        planet['R'] = 3389.5e3 # (m) Volumetric mean radius of Mars
        planet['T'] = 88642 # (s) Mars sidereal day
    elif planetName =='EARTH':
        planet['m'] = 5.9742e24
        planet['R'] = 6371.0008e3#(m) Volumetric mean radius
        planet['T'] = 86164.0906#(s) Earth sidereal day. NOTE sidereal day != Solar day
    else:   
        print 'ERROR in PlanetLibrary.py : Planet not in current library'
        exit()


    mission['planet'] = planet

#mission = RefVernalEquinox(mission)# adding Vernal Equinox data for IJK frame

    return mission

def CalculateAtmosphere(alt,planet):
    
    planetName = planet['name'];
    
    if planetName == 'EARTH':
        rho = 1.225*np.exp(-754.97611931455*alt) #alt is non dim. From Cantwell's reader
    elif planetName == 'MARS': # this model is from Ross Allen's AA 290 work
        #[~,~,atm_rho]= MarsAtmV3(alt*3389.5e3);
        #[~, ~, ~, max_rho] = MarsAtmV3(0);
        #for i = 1:length(atm_rho) 
        #    if atm_rho(i) > max_rho || atm_rho(i) < 0
        #        atm_rho(i)=0;
        #    end
        #    if isnan(atm_rho(i)) || isinf(atm_rho(i))
        #        atm_rho(i) = 0;
        #    end
        
        #end
        #rho = atm_rho;
        error('ERROR in planetLibrary.py. MARS is not yet a valid option. ')
        exit(1)
    else:
        error('ERROR in planetLibrary.py. Atmosphere model cannot be identified');
        exit(1)
    return rho   




def getSpeedOfSound(hvec,planet,flag):
    #F2PY can be easily done in this part. FOR LOOP CAN BE A KILLER...F2PY soon to come
    # speed of sound derived by looking at the linear behavior as a function of altitude
    # hvec is the altitude array in meters
    # c is im m/sec  
    
    if flag=='real':
        c = np.zeros(len(hvec))
    elif flag=='complex':
        c = np.zeros(len(hvec),complex)
    
    
    if planet['name'] =='EARTH':
        
        
        for index in range(len(hvec)):
            h = hvec[index]
            hkm0 = h*.001
            hkm = hkm0.real # making sure we are comparing real numbers
            

            if (hkm<0.):
                c[index] = 340.3
            elif(hkm>86.) :
                c[index] = 274.1
            else:

                if (hkm<=11):
                    c2 = 295.2
                    c1 = 340.3
                    h2 = 11.
                    h1 = 0.
                elif (hkm<=20):
                    c2 = 295.11
                    c1 = 295.1
                    h2 = 20.
                    h1 = 11.
                elif (hkm<=32):
                    c2 = 303.0
                    c1 = 295.1
                    h2 = 32.
                    h1 = 20.
                elif (hkm<=47):
                    c2 = 329.2
                    c1 = 303.
                    h2 = 47.
                    h1 = 32.
                elif (hkm<=48):
                    c2 = 329.8
                    c1 = 329.2
                    h2 = 48.
                    h1 = 47.
                elif (hkm<=51):
                    c2 = 329.81
                    c1 = 329.8
                    h2 = 51.
                    h1 = 48.
                elif (hkm<=52):
                    c2 = 328.8
                    c1 = 329.82
                    h2 = 52.
                    h1 = 51.
                elif (hkm<=72):
                    c2 = 293.4
                    c1 = 328.8
                    h2 = 72.
                    h1 = 52.
                elif (hkm<=86):
                    c2 = 274.1
                    c1 = 293.4
                    h2 = 86.
                    h1 = 72.

                m = (c2-c1)/(h2-h1) #slope calculation
                c[index] = m*(hkm0-h1) + c1 # speed of sound calculation
        return c
    else:
        error('ERROR in planetLibrary.py. Speed of sound model cannot be identified');
        exit(1)
        


