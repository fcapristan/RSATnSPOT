#this routine will wrap different functions to return a weighted sobol calculation in OT
# developed by Francisco Capristan
from openturns import *
import sys
import numpy as np
def wrapSobol(wInList=None,OTmodelListIn=None,OTinputList=None,seedList1=None,seedList2=None,nCases=None,blocksize=1):
  
    ##myDistList = []
    #myVectXList = []
    #setting up input parameters the way OT accepts them
    firstOrderList = []
    totalOrderList = []
    varList = []
    varTotal = 0.0
    for index in range(len(OTinputList)):
        OTinput = OTinputList[index]
        myDist,VectX = setInputRandomVector(OTinput)
        #myDistList.append(myDist)
        #myVectXlist.append(VectX)
        if seedList1[index]==seedList2[index]:
            print 'seeds cannot be the same value in OT'
            sys.exit(1)
        RandomGenerator.SetSeed(seedList1[index])
        myRandomExp1 = LHSExperiment(myDist,nCases)
        inp1 = myRandomExp1.generate()
        RandomGenerator.SetSeed(seedList2[index])
        myRandomExp2 = LHSExperiment(myDist,nCases)
        inp2 = myRandomExp2.generate()
        model = OTmodelListIn[index]
        sensAnalysis = mySensitivityAnalysis(inp1,inp2,model)
        sensAnalysis.setBlockSize(blocksize)
        firstOrderIndices = sensAnalysis.getFirstOrderIndices()

        
        totalOrderIndices = sensAnalysis.getTotalOrderIndices()
        firstOrderList.append(firstOrderIndices)
        totalOrderList.append(totalOrderIndices)
        out1 = sensAnalysis.getOutVals1()
        std1 = out1.computeStandardDeviation()[0,0]
        var = std1**2.0
        varList.append(var)
        varTotal = (wInList[index]**2.0)*var + varTotal
    # done doing all major computations. Now to rescale the sobol indices to get the actual values
    #First order indices first
    Xfirst = []
    XTotal = []
    print 'VarTotal',varTotal
    print 'STD',varTotal**.5
    for index1 in range(len(OTinputList)):
        OTinput = OTinputList[index1]
        for index2 in range(len(OTinput)):
            sobolFirst = firstOrderList[index1][index2]
            sobolTotal = totalOrderList[index1][index2]
            sobolNewFirst = (wInList[index1]**2.0 * varList[index1]*sobolFirst)/varTotal
            sobolNewTotal = (wInList[index1]**2.0 * varList[index1]*sobolTotal)/varTotal
            Xfirst.append(sobolNewFirst)
            XTotal.append(sobolNewTotal)
            
    return (Xfirst,XTotal)
    
def setInputRandomVector(inList):
    #setting the input random vector X
    myCollection = DistributionCollection(len(inList))
    for index in range(len(inList)):
        myCollection[index] = Distribution(inList[index])
    myDistribution = ComposedDistribution(myCollection)
    VectX = RandomVector(Distribution(myDistribution))
    return myDistribution,VectX
