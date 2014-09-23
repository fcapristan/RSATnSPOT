# routines that will assemble to be passed to the optimizer

import numpy as np

def assembleVariables(x,y,z,Vx,Vy,Vz,theta1,theta2,eta,tf):
    vars = np.concatenate((x,y,z,Vx,Vy,Vz,theta1,theta2,eta,tf),1)
    return vars

def decomposeVariables(vars,rangeVar):  
    x = vars[rangeVar['x']]
    y = vars[rangeVar['y']]
    z = vars[rangeVar['z']]
    Vx = vars[rangeVar['Vx']]
    Vy = vars[rangeVar['Vy']]
    Vz = vars[rangeVar['Vz']]
    theta1 = vars[rangeVar['theta1']]
    theta2 = vars[rangeVar['theta2']]
    eta = vars[rangeVar['eta']]
    tf = vars[rangeVar['tf']]
    return (x,y,z,Vx,Vy,Vz,theta1,theta2,eta,tf)


def decomposeDirectionVariables(vars,rangeVar):
    x = vars[rangeVar['x']]
    y = vars[rangeVar['y']]
    z = vars[rangeVar['z']]
    Vx = vars[rangeVar['Vx']]
    Vy = vars[rangeVar['Vy']]
    Vz = vars[rangeVar['Vz']]
    theta1 = vars[rangeVar['theta1']]
    theta2 = vars[rangeVar['theta2']]
    return (x,y,z,Vx,Vy,Vz,theta1,theta2)

def decomposeLocationVariables(vars,rangeVar):
    x = vars[rangeVar['x']]
    y = vars[rangeVar['y']]
    z = vars[rangeVar['z']]
    return (x,y,z)