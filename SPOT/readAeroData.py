import numpy as np
import sys
# these set of functions to be used to read CD and CL as function of Mach numner...--> to include AOA in later versions

def readInputCD(fileName):
    try:
        inputFile = open(fileName,'r')
    except:
        print '\n!!! Error: Could not open input file: ' + str(fileName) + ' !!!\n'
        exit(1)
    arefstring = 'Aref'
    Aref = np.array([])
    CD = np.array([])
    Minf = np.array([])
    for line in inputFile:
        if line[0] == '#' or line.find("#") > 0:
            iComment = line.index('#')
            line = line[:iComment]
        if line.find(":") > 0:
            iColon = line.index(":")
            key = (line[:iColon].lower()).split()
            value = (line[iColon+1:]).split()
            if key[0].lower() == arefstring.lower():
                if len(value) < 1:
                    print '\n!!! Error: Invalid input arguments !!!\n' \
                    + ' At line: ' + line.strip() + '\n'
                    exit(1)
                Aref = float(value[0])
        elif len(line.split())==2:
            key = line.split()
            Minf = np.concatenate((Minf,[float(key[0])]),1)
            CD = np.concatenate((CD,[float(key[1])]),1)
    return (Minf,CD,Aref)