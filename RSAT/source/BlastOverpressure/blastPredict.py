import numpy as np
import scipy
# set of routines for blast overpressure calculation





def getDangerArea(alt=None,pcritical=None,Et=None,Mf=None):
    # returns the casualty area centered at lon0 and lat0
    REarth = 6378145. #[m]
    p0 = stdatmf(alt)
    Rblast = 0.0
    pbarcrit = 0.0
    if p0>0.:
        Rblast,pbarcrit = getRcritical(pcritical,p0,Et,Mf)
    if alt>Rblast:
        area = 0
    else:
        
        #area = (Rblast**2-alt**2)*np.pi # this assumes a flat EARTH...GOOD APPROXIMATION!!
        
        # assumes spherical Earth
        d = REarth + alt
        a = (REarth**2-Rblast**2 + d**2)/(2.*d)
        alpha = np.arccos(a/REarth)
        area = (-REarth**2)*(np.cos(alpha) - 1.)*2.*np.pi
    
    return area,pbarcrit


def solidPropellantYield(mass = None,S='soft',v = None):
    kg_to_lb = 2.204622622
    ft_to_m = 0.3048
    m_to_ft = 1./ft_to_m
    
    mass_lb = mass*kg_to_lb
    v_ft_sec = v*m_to_ft
    if S.lower()=='soft':#soft for soft soil
        Sval = 1.81
    elif S.lower()=='water':
        Sval = 2.92
    elif S.lower()=='concrete':
        Sval = 1.41
    elif S.lower()=='steel':
        Sval = 1.0
    
    fracTNT = 1.28*(1.0 + 192000./(mass_lb**0.156) *(Sval/v_ft_sec)**1.55)**-1.0
    #print 'fracTNT',fracTNT,mass_lb,v_ft_sec
    massTNT = mass*fracTNT
    return massTNT


def getDangerAreaGround(pcritical,Et,Mf):
    p0 = stdatmf(0.0)
    Rblast,pbarcrit = getRcritical(pcritical,p0,Et,Mf)
    area = (Rblast**2.)*np.pi
    return area,Rblast

def getSimpleRadiusTNT(massTNT=None,TNT_K_factor=45.):
    #equation from Flight Safety Analysis Handbook 
    #K factor of 45 for 1 psi overpressure. K factor of 20 for 3psi
    lb_to_kg = 0.453592
    kg_to_lb = 2.204622622
    ft_to_m = 0.3048
    R = TNT_K_factor*(massTNT*kg_to_lb)**(1./3.)
    R = ft_to_m*R
    return R

def getDangerZone(lon0=None,lat0=None,alt=None,pcritical=None,
                  p0=None,Et=None,Mf=None,lonMesh=None,latMesh=None,TNT_mass = None,TNT_Kfactor=45.,tangCurve=True):
    # import matplotlib.pyplot as plt    
    # from pylab import *
    REarth = 6378145. #[m]
    Rlocal = REarth+alt
    kg_to_lb = 2.204622622
    ft_to_m = 0.3048
    # origin for blast overpressure
    x0 = Rlocal*np.cos(lat0*np.pi/180.)*np.cos(lon0*np.pi/180.)
    y0 = Rlocal*np.cos(lat0*np.pi/180.)*np.sin(lon0*np.pi/180.)
    z0 = Rlocal*np.sin(lat0*np.pi/180.)
    
    if tangCurve==True:
        Rblast = getRcritical(pcritical,p0,Et,Mf)
    else :
        #eqaution from Flight Safety Analysis Handbook 
        Rblast = TNT_Kfactor*(TNT_mass*kg_to_lb)**(1./3.)
        Rblast = Rblast*ft_to_m
    #print Rblast
    if Rblast <alt:
        print Rblast,alt     
        # no casualties...blast critical never reaches the ground
    	return 0
    else : # in case it reaches the ground
        
    	lonMeshrad = np.pi/180.*lonMesh
    	latMeshrad = np.pi/180.*latMesh
        xMesh = REarth*np.cos(latMeshrad)*np.cos(lonMeshrad)
    	yMesh = REarth*np.cos(latMeshrad)*np.sin(lonMeshrad)
    	zMesh = REarth*np.sin(latMeshrad)
        
       	# Equation of a sphere
       	Rvals2 =  (xMesh-x0)**2 + (yMesh-y0)**2 + (zMesh-z0)**2
        Rcheck = Rblast**2 - Rvals2
        
        Rcheck[Rcheck>=0] = 1.
        Rcheck[Rcheck<0] = 0.
        #print Rcheck
        #plt.contour(lonMesh,latMesh,Rcheck)
        #imshow(Rcheck)
        #colorbar()
        #plt.show()
        #print np.sum(Rcheck)
        
        return Rcheck
    
    
    ''' # Originally started calculation for 2 spheres intersecting,
    # but decided to avoid it (avoid loss of generality in case we want to use an oblate Earth)
    # essentially Earth curvature does not have a big impact as long as the explosive sphere is a lot smallet than the planet
    x = Rlocal*np.cos(lat0*np.pi/180.)*np.cos(lon0*np.pi/180.)
    y = Rlocal*np.cos(lat0*np.pi/180.)*np.sin(lon0*np.pi/180.)
    z = Rlocal*np.sin(lat0*np.pi/180.)
    vec0 = [x,y,z]
    d = np.sqrt(x**2+y**2+z**2)
    Rcrit = getRcritical(pcritical,p0,Et,Mf)
    print Rcrit,REarth**2-((REarth**2 + d**2 - Rcrit**2)/(2.0*d))**2
    a = np.sqrt(REarth**2-((REarth**2 + d**2 - Rcrit**2)/(2.0*d))**2) # from sphere intersection, radius of intersected area
    theta = np.arctan2(a,REarth)
    # must add next part here
    lonMesh = np.pi/180.*lonMesh
    latMesh = np.pi/180.*latMesh
    xMesh = np.cos(latMesh)*np.cos(lonMesh)
    yMesh = np.cos(latMesh)*np.sin(lonMesh)
    zMesh = np.sin(latMesh)
    
    iN,jN = np.shape(xMesh)
    retMesh = np.zeros((iN,jN))
    for i in range(0,iN):
    for j in range(0,jN):
    vec = [xMesh[i,j],yMesh[i,j],zMesh[i,j]]
    angle = np.dot(vec,vec0)
    if angle<= theta:
    retMesh[i,j] = 1
    return retMesh
    '''
def getlatlon(x,y,z):
    
    #radEarth = 6378145. #[m]
    rmag = np.sqrt(x**2+y**2+z**2)
    longitude = np.arctan2(y,x)*180./np.pi
    latitude = 90. - np.acos(z/rmag)*180./np.pi
    return (latitude,longitude)




def getRcritical(pcritical,p0,Et,Mf):
    #inputs critrial pressure, ambient p, Energy, flame Mach number
    pbarcrit = (pcritical-p0)/p0
    Rbar,pbar = getRbarPbar(Mf)
    Rbarcrit = np.interp(pbarcrit,pbar,Rbar)
    Rcrit = (Et/p0)**(1./3.)*Rbarcrit
    
    return Rcrit,pbarcrit


def getPfromVCE(R,Et,p0,Mf):
    # set of curves obtained from "A NEW SET OF BLAST CURVES FROM VAPOR CLOUD EXPLOSION"
    #  M.J Tang and Q.A Baker
    
    RbarCalc = R/(Et/p0)**(1.0/3.0)
    (Rbar,PbarMf5p2) = getRbarPbar(5.2)
    
    Pbar = np.interp(RbarCalc,Rbar,PbarMf5p2)
    p = Pbar*p0+p0
    
    return p

def getRbarPbar(Mf):
    
    if Mf == 5.2:
        
        Rbar = np.array([0.104398,0.118786,0.137111,0.156008,0.170032,0.181371,0.192083,0.207853,0.224918,0.243384,\
                         0.276928,0.312843,0.348382,0.407936,0.498677,0.622863,0.724124,0.910965,1.13782,1.41102,\
                         1.74981,2.12375,2.63367,3.06183,3.82432,4.98677,6.22863,7.61412,8.97989,9.71715])
        
        Pbar = np.array([20.1279,20.1278,20.1277,20.1276,20.6227,14.3258,9.47967,7.43546,5.97546,4.57444,3.33585,\
                         2.43263,1.95497,1.35804,1.01468,0.722187,0.566454,0.41308,0.324003,0.248036,0.194549,\
                         0.152596,0.125648,0.100976,0.0773011,0.0563709,0.0411078,0.0355331,0.0285559,0.0246834])
        indexing = np.arange(len(Pbar)-1,-1,-1)
        Pbar = Pbar[indexing]#making sure this is in ascensing order for interpolation
        Rbar = Rbar[indexing]
    else:
        print "Error in blastPredict.py, Mf not in current library \n"
        exit() 
    return (Rbar,Pbar)


class yieldLiquid:
    def __init__(self):
        self.name = None

    def propellantType(self,namePropellant,surface):

        ft_to_m = 0.3048
        if (namePropellant.lower()=='lo2/lh2' and surface.lower()=='soft'):
            x = [0,82.7471,528.85,800.]
            y = [.17,.17,1.52,1.52]
        elif (namePropellant.lower()=='lo2/rp-1' and surface.lower()=='soft'):
            x = [0,767.369,800.]
            y = [0.0,1.02424,1.02424]
        elif (namePropellant.lower()=='hypergols' and surface.lower()=='soft'):
            x = [0,195.863,296.176,495.702,800.]
            y = [0.0277695,0.439647,0.556279,0.610864,0.623251]
        elif (namePropellant.lower()=='lo2/lh2' and surface.lower()=='hard'):
            x = [0,82.7471,738.68,800.]
            y = [0.17,0.17,1.52368,1.52368]
        elif (namePropellant.lower()=='lo2/rp-1' and surface.lower()=='hard'):
            x= [0,165.693,800.]
            y = [0.00397017,0.224803,0.224803]        
        elif (namePropellant.lower()=='hypergols' and surface.lower()=='hard'):
            x = [1.43326,
                 40.7127,
                 85.8206,
                 130.945,
                 183.346,
                 248.862,
                 318.742,
                 391.544,
                 493.481,
                 569.237,
                 652.27,
                 796.504]
            y = [0.0238427,
                 0.0561898,
                 0.0886237,
                 0.113117,
                 0.141689,
                 0.170457,
                 0.20326,
                 0.232136,
                 0.265416,
                 0.274484,
                 0.287631,
                 0.30169]
        else:
            print 'Liquid Propellant name does not match'
            print 'Error in blastPredic.py'
            exit()
        self.fun = scipy.interpolate.interp1d(ft_to_m*np.array(x),np.array(y),kind='linear',fill_value=y[-1])

    def yieldFactor(self,vel):
        if vel <0:
            print 'Error: Impact Velocity cannot be less than zero. Error in blastPredict.py'
            exit()
        
        return self.fun(vel)
# this module translated from matlab version


#   *********** 1976 STANDARD ATMOSPHERE SUBROUTINE **********
#
#     Mason's BASIC program, converted to FORTRAN - Sept. 1, 1989
#     converted to MATLAB, by Paul Buller, 1998
#     converted to a MATLAB function by w.h. mason, 2001
#
#
#     W.H. Mason
#     Department of Aerospace and Ocean Engineering
#     Virginia Tech, Blacksburg, VA 24061
#     email: mason@aoe.vt.edu
#
#     k  -  = 0 - metric units
#          <> 0 - English units
#
#     KK -   0 - good return
#            1 - error: altitude out of table,
#                 do not use output (max altitude for this
#                 routine is 84.85 km or 282,152 ft.)
#
#     Z  - input altitude, in feet or meters (depending on k)
#
#     output:
#                      units: metric        English
#     T  - temp.               deg K         deg R
#     P  - pressure            N/m**2         lb/ft**2
#     R  - density (rho)       Kg/m**3        slug/ft**3
#     A  - speed of sound      m/sec         ft/sec
#     MU - viscosity           Kg/(m sec)    slug/<ft sec)
#     
#     TS - t/t at sea level
#     RR - rho/rho at sea level
#     PP - p/p at sea level
#
#     RM - Reynolds number per Mach per unit of length
#     QM - dynamic pressure/Mach**2
#

def stdatmf(Z,k=0):
    #function [T,R,P,A,MU,TS,RR,PP,RM,QM,KK] = stdatmf(Z,k)
    KK = 0
    K = 34.163195
    C1 = 3.048e-4
    T = 1
    PP = 0
    
    if (k==0):
        TL = 288.15
        PL = 101325
        RL = 1.22152281 #1.225 I switched.
        C1 = 0.001
        AL = 340.294
        ML = 1.7894e-5
        BT = 1.458e-6
    else:
        TL = 518.67
        PL = 2116.22
        RL = 0.0023769
        AL = 1116.45
        ML = 3.7373e-7
        BT = 3.0450963e-8
    
    H = C1*Z/(1 + C1*Z/6356.766)
    
    if (H<11.):
        T = 288.15 - 6.5*H
        PP = (288.15/T)**(-K/6.5)
    elif (H<20):
        T = 216.65
        PP = 0.22336*np.exp(-K*(H-11)/216.65)
    elif (H<32):
        T = 216.65 + (H-20)
        PP = 0.054032*(216.65/T)**K
    elif (H<47):
        T = 228.65 + 2.8*(H-32)
        PP = 0.0085666*(228.65/T)**(K/2.8)
    elif (H<51):
        T = 270.65
        PP = 0.0010945*np.exp(-K*(H-47)/270.65)
    elif (H<71):
        T = 270.65 - 2.8*(H-51)
        PP = 0.00066063*(270.65/T)**(-K/2.8)
    elif (H<84.852):
        T = 214.65 - 2*(H-71)
        PP = 3.9046e-5*(214.65/T)**(-K/2)
    else:
        # chosen altitude too high
        KK = 1
    
    
    M1 = np.sqrt(1.4*287*T)
    RR = PP/(T/288.15)
    MU = BT*T**1.5/(T+110.4)
    TS = T/288.15
    A = AL*np.sqrt(TS)
    T = TL*TS
    R = RL*RR
    P = PL*PP
    RM = R*A/MU
    QM = 0.7*P
    #    return [T,R,P,A,MU,TS,RR,PP,RM,QM,KK]
    return P



'''
    
# test driver
lon0 = 30.1
lat0 = 45.3
alt = 100
p0 = 1e5
pcritical = 1.025*p0
Et = 10000000000.
Mf = 5.2
lonVec = np.linspace(lon0-.005,lon0+.005)
latVec = np.linspace(lat0-.005,lat0+.005)

lonMesh,latMesh = np.meshgrid(lonVec,latVec)

#print getDangerZone(lon0,lat0,alt,pcritical,p0,Et,Mf,lonMesh,latMesh)

print getDangerArea(pcritical,Et,Mf)






'''

