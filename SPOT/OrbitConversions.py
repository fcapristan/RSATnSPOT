import numpy as np


def OrbitalVectors(a,e,inc,komega,omega,nu,mu):

#****************************************************************************************#* 										        *
#* Converts Standard Orbital Elements to Cartesian			    		*
#* 						                                        *
#* Created by:          Francisco Capristan                               			*
#* Created on:          07/21/2012 							* 
#											*
#											*
# Inputs - 										*
#	 a      = semimajor axis							*
#	 e      = eccentricity								*
#	 inc    = inclination (radians)							*
#	 komega = Longitude of the ascending node (radians)				*
#	 omega  = argument of periapsis (radians)					*
#	 nu	= true anomaly (radians) 						*
#	 mu 	= G*M (usually 1-->non-dimensionalized)					*
# Outputs 										*
#	 r  	= location in cartesian coordinates (x y z) 				*
# 	 V  	= velocity in cartesian cordinates (Vx Vy Vz)				*
#	 err	= error code								*
#											*
#											*
# 											*
# Equations from Fundamentals of Astrodynamics by Bate, Mueller, and White....pg ~82	*
#****************************************************************************************

    p = a*(1-e**2)
    radius = p/(1.+e*np.cos(nu))
    rv = np.zeros((3))
    vv = np.zeros((3))
    rotation = np.zeros((3,3))

                                
    rv[0]=radius*np.cos(nu)
    rv[1]=radius*np.sin(nu)
    rv[2]=0.
    vv[0]=-np.sqrt(mu/p)*np.sin(nu)
    vv[1]=np.sqrt(mu/p)*(e+np.cos(nu))
    vv[2]=0.
    
    rotation[0,0] = np.cos(komega)*np.cos(omega)-np.sin(komega)*np.sin(omega)*np.cos(inc)
    rotation[0,1] = -np.cos(komega)*np.sin(omega)-np.sin(komega)*np.cos(omega)*np.cos(inc)
    rotation[0,2] = np.sin(komega)*np.sin(inc)
    
    rotation[1,0] = np.sin(komega)*np.cos(omega)+np.cos(komega)*np.sin(omega)*np.cos(inc)
    rotation[1,1] = -np.sin(komega)*np.sin(omega)+np.cos(komega)*np.cos(omega)*np.cos(inc)
    rotation[1,2] = -np.cos(komega)*np.sin(inc)
    
    rotation[2,0] = np.sin(omega)*np.sin(inc)
    rotation[2,1] = np.cos(omega)*np.sin(inc)
    rotation[2,2] = np.cos(inc)
    
    
    r = np.dot(rotation,rv).T
    V = np.dot(rotation,vv).T
    angularMomentum = np.cross(r,V)
    eVector = np.cross(V,angularMomentum)/mu-r/np.linalg.norm(r)

    return angularMomentum,eVector





