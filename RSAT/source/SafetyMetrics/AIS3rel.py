import numpy


def AIS(V,m):
    # takes the impact velocity relative to rotating planet in m/s and mass in kg
    # Last check Nov 12,2014. Equations seems to work as intended
    casualty = 0
    W = m*9.81

    # This values obtained from --Estimation of Space Shuttle Orbiter Reentry Debris Casulaty Area-- by Jon D. Collins et all
    # Plot was digitized and a power fitt applied to get the constants (everything then converted to SI Units)
    a = 29.63
    b = -0.8934
    c = 3.6
    #*****
    Vlim = 29.63*W**b+c
    if V>=Vlim:
        casualty = 1
    return casualty
