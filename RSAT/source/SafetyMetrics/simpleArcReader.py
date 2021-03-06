#import sys
import numpy as np
import math
#xvec = [22,33,60]
#yvec = [19,40,60]


    #filename = 'trial.asc'
def arcReader(xvec,yvec,totalLength,filename):
    key,xllcorner,yllcorner,cellsize,nrows,xMax = getRawVals(filename)
    ValueVector = getVal(key,xllcorner,yllcorner,cellsize,nrows,xMax,xvec,yvec,totalLength)
    return ValueVector
              
def getRawVals(filename):
    try:
       inputFile = open(filename,'r')
    except:
       print 'Cannot open file==> Error in simpleArcReader.arcReader'


    #Initializing values for all parameters
    ncols = -1
    nrows = -1
    xllcorner = -1
    yllcorner = -1
    cellsize = -1
    noDataValue = 0

    ycount = 0
    xcount = 0

    for line in inputFile:
       key = line.split()
       
       if len(key)==2:
          if key[0] == 'ncols':
             ncols = float(key[1])
          elif key[0] == 'nrows':
             nrows = float(key[1])
          elif key[0] == 'xllcorner':
             xllcorner = float(key[1])
          elif key[0] == 'yllcorner':
             yllcorner = float(key[1])
             ylower = yllcorner
          elif key[0] == 'cellsize':
             cellsize = float(key[1])
          elif key[0] == 'NODATA_value':
             noDataValue = float(key[1])
             break
       elif len(key)>2 and cellsize>0:
          print 'ERROR: Please specify a NODATA_Value after the CELLSIZE==> simpleArcReader.arcReader'
          exit(1)

    ySpan = nrows*cellsize
    xSpan = ncols*cellsize
      
    xMax = xllcorner + xSpan
    yMax = yllcorner + ySpan

    key =np.zeros([nrows,ncols])
#    print line
    counter = 0
    for line2 in inputFile:
       key[counter][:]=np.fromstring(line2,dtype=float,sep=' ')
       counter = counter + 1

    return (key,xllcorner,yllcorner,cellsize,nrows,ncols,xMax,yMax)

def getRawValsPopDensity(filename,flag='complete'):
    try:
        inputFile = open(filename,'r')
    except:
        print 'Cannot open file==> Error in simpleArcReader.arcReader'
        exit(1)
    
    #Initializing values for all parameters
    ncols = -1
    nrows = -1
    xllcorner = -1
    yllcorner = -1
    cellsize = -1
    noDataValue = 0
    
    ycount = 0
    xcount = 0
    
    for line in inputFile:
        key = line.split()
        
        if len(key)==2:
            if key[0] == 'ncols':
                ncols = int(key[1])
            elif key[0] == 'nrows':
                nrows = int(key[1])
            elif key[0] == 'xllcorner':
                xllcorner = float(key[1])
            elif key[0] == 'yllcorner':
                yllcorner = float(key[1])
                ylower = yllcorner
            elif key[0] == 'cellsize':
                cellsize = float(key[1])
            elif key[0] == 'NODATA_value':
                noDataValue = float(key[1])
                break
        elif len(key)>2 and cellsize>0:
            print 'ERROR: Please specify a NODATA_Value after the CELLSIZE==> simpleArcReader.arcReader'
            exit(1)
    
    ySpan = nrows*cellsize
    xSpan = ncols*cellsize
    
    xMax = xllcorner + xSpan
    yMax = yllcorner + ySpan
    

    if flag == 'simple':
                
        return (xllcorner,yllcorner,cellsize,nrows,ncols,xMax,yMax)
    else:
        key =np.zeros([nrows,ncols])
        #    print line
        counter = 0
        for line2 in inputFile:
            key[counter][:]=np.fromstring(line2,dtype=float,sep=' ')
            counter = counter + 1
        if counter!=int(nrows):
            print 'Error reading Population data'
            exit(1)
        
        return (key,xllcorner,yllcorner,cellsize,nrows,ncols,xMax,yMax)
def getIndex(xllcorner,yllcorner,cellsize,nrows,xMax,xvec,yvec,totalLength):
    valueVector = np.zeros(totalLength,2)

    for index in range(len(xvec)):
        x = xvec[index]
        y = yvec[index]
        nRow = int(nrows) - int(math.floor((y - yllcorner)/cellsize))
        nColumn = 1 + int(math.floor((x-xllcorner)/cellsize))
        if x ==xMax:
            nColumn = nColumn -1
        elif x>xMax:
            print 'ERROR: Check X bounds in ASC file'
            exit(1)
        
        valueVector[index,:]=np.array([nRow-1,nColumn-1])
    return valueVector

#def getValfromIndex(line,valueVector):
# line is simply the file line obtained from line = open(filename,'r')
# valueVector obtained from getIndex


def getVal(key,xllcorner,yllcorner,cellsize,nrows,xMax,xvec,yvec,totalLength):
    valueVector = np.zeros(totalLength)

    counter = 0
    for index in range(len(xvec)):
       x = xvec[index]
       y = yvec[index]
       nRow = int(nrows) - int(math.floor((y - yllcorner)/cellsize))
       nColumn = 1 + int(math.floor((x-xllcorner)/cellsize))
       if x ==xMax:
          nColumn = nColumn -1
       elif x>xMax:
          print 'ERROR: Check X bounds in ASC file'
          exit(1)

       valueVector[counter]=key[nRow-1][nColumn-1]
       counter = counter+1
    return valueVector


