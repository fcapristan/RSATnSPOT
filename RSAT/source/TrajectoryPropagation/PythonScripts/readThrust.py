import numpy as np
import sys
# these set of functions to be used to read CD and CL as function of Mach numner...--> to include AOA in later versions

def readInputThrust(fileName):
    try:
        inputFile = open(fileName,'r')
    except:
        print '\n!!! Error: Could not open input file: ' + str(fileName) + ' !!!\n'
        exit(1)
    time = np.array([])
    Tx = np.array([])
    Ty = np.array([])
    Tz = np.array([])

    for line in inputFile:
        if line[0] == '#' or line.find("#") > 0:
            iComment = line.index('#')
            line = line[:iComment]
        if line.find(":") > 0:
            iColon = line.index(":")
            key = (line[:iColon].lower()).split()
            value = (line[iColon+1:]).split()
            if key[0].lower() == 'stage':
                if len(value) < 1:
                    print '\n!!! Error: Invalid input arguments !!!\n' \
                    + ' At line: ' + line.strip() + '\n'
                    exit(1)
                stage = float(value[0])
        elif len(line.split())==4:
            key = line.split()
            time = np.concatenate((time,[float(key[0])]),1)
            Tx = np.concatenate((Tx,[float(key[1])]),1)
            Ty = np.concatenate((Ty,[float(key[2])]),1)
            Tz = np.concatenate((Tz,[float(key[3])]),1)
    time = time - time[0]
    return (stage,time,Tx,Ty,Tz)