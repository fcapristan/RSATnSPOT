#this routine will wrap different functions to return a weighted sobol calculation in OT
# developed by Francisco Capristan
from openturns import *
import sys
import numpy as np
import openTurnsTools as OTtools
import time
def wrapSobol(wInList=None,OTmodelListIn=None,OTinputList=None,seedList1=None,seedList2=None,nCases=None,blocksize=1,varTotal = None):
  
    ##myDistList = []
    #myVectXList = []
    #setting up input parameters the way OT accepts them
    firstOrderList = []
    totalOrderList = []
    varList = []
    varTotalC = 0.0
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
        varTotalC = (wInList[index]**2.0)*var + varTotalC
    if varTotal==None:
        varTotal = varTotalC
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


def fTotal(X,indexVariables,functionList,weightList):
    X = np.array(X)

    valOut = 0.0
    for index in range(len(weightList)):
        Xin = X[indexVariables[index]]

        valC = functionList[index](Xin)
        valOut = weightList[index]*valC + valOut

    return valOut


'''
def fTotal(X,indexArray,funList,weightArray):
    Xvars = getXsingle(X,indexArray)
    valOut = 0.0
    for index in range(len(Xvars)):
        cval = funList[index](Xvars[index])

        valOut =  valOut + weightArray[index]*cval

    return [valOut] 

def getXsingle(X,indexList):
    Xvar = []
    X = np.array(X)

    for index in range(len(indexList)):

        Xvar.append(np.array(X[indexList[index]]))

    return Xvar

'''


def sobolMod(weightList = None, functionList = None,nCases=None,blocksize=None,sharedIndex = None,ImportantBoolIndex = None,indexVariable = None,seedList = None,varTotal=-9999,multiPro=1):

    MainIndexOut = np.zeros(len(ImportantBoolIndex)) -9999
    TotalIndexOut = np.zeros(len(ImportantBoolIndex)) -9999
    if blocksize==None:
        blocksize = nCases
    # first let's do the shared Input Variables
    [seed1,seed2,seed3] = seedList
    sharedIndex = np.array(sharedIndex)

    sharedBoolIndex = ImportantBoolIndex[sharedIndex]
    
    sharedIndexImp = sharedIndex[sharedBoolIndex]

    checkBool = np.sum(sharedBoolIndex)
    varTotal = -999
    time0 = time.time()
    if checkBool >=1: # do sobol indices over shared function
        sharedBoolIndex = ImportantBoolIndex
        sharedBoolSobol = np.array([False]*len(ImportantBoolIndex))
        sharedBoolSobol[sharedIndexImp] = True
        f = lambda X:fTotal(X,indexVariable,functionList,weightList) # full function
        inList = [Uniform(0.0,1.0)]*len(ImportantBoolIndex)
     
        inp1,inp2,inpOther = OTtools.LHSExperimentWithIndexForSobol(inList,sharedBoolSobol,seed1,seed2,seed3,nCases)
        funcXimp = lambda Ximp:OTtools.getValOutfromSubset(Ximp,inpOther,sharedBoolSobol,f)

        newSamp,newDim = np.shape(inp1)
        OTfunXimp = OTtools.myfunction(multiPro=multiPro,dim=newDim,fun=funcXimp)
        model = NumericalMathFunction(OTfunXimp)
        sensAnalysis = mySensitivityAnalysis(inp1,inp2,model)

        sensAnalysis.setBlockSize(nCases)
        print 'starting2'

        firstOrderIndicesShared = sensAnalysis.getFirstOrderIndices()
        totalOrderIndicesShared = sensAnalysis.getTotalOrderIndices()
        print 'done'
        out1 = sensAnalysis.getOutVals1()
        varTotal  = (out1.computeStandardDeviation()[0,0])**2

        counter = 0
        for index in range(len(sharedBoolSobol)):
            if sharedBoolSobol[index] == True:
                MainIndexOut[index] = firstOrderIndicesShared[counter + 1]
                TotalIndexOut[index] = totalOrderIndicesShared[counter + 1]
                counter = counter + 1
    # Now doing the rest of the sobol indices
    varList = []

    for index in range(len(functionList)):
        currentIndexBool = ImportantBoolIndex[indexVariable[index]]
        currentIndexBool[sharedIndex] = False
        checkBoolSum = np.sum(currentIndexBool)
    
        if checkBoolSum>=1: # important index in given 
            inList = [Uniform(0.0,1.0)]*len(currentIndexBool)

            inp1,inp2,inpOther = OTtools.LHSExperimentWithIndexForSobol(inList,currentIndexBool,seed1+index,seed2+index,seed3+index,nCases)
            flocal = functionList[index]
            funcXimp = lambda Ximp:OTtools.getValOutfromSubset(Ximp,inpOther,currentIndexBool,flocal)

            newSamp,newDim = np.shape(inp1)
            OTfunXimp = OTtools.myfunction(multiPro=multiPro,dim=newDim,fun=funcXimp)
            model = NumericalMathFunction(OTfunXimp)
            sensAnalysis = mySensitivityAnalysis(inp1,inp2,model)

            sensAnalysis.setBlockSize(nCases)
            firstOrderIndices = sensAnalysis.getFirstOrderIndices()
            totalOrderIndices = sensAnalysis.getTotalOrderIndices()
            out1 = sensAnalysis.getOutVals1()
            std1 = out1.computeStandardDeviation()[0,0]
            var = std1**2.0
            sobolNewFirstList = []
            sobolNewTotalList = []
            for index1 in range(len(firstOrderIndices)-1):
                sobolFirst = firstOrderIndices[index1+1]
                sobolTotal = totalOrderIndices[index1+1]
                sobolNewFirst = (weightList[index]**2.0 * var*sobolFirst)/varTotal
                sobolNewTotal = (weightList[index]**2.0 * var*sobolTotal)/varTotal
                sobolNewFirstList.append(sobolNewFirst)
                sobolNewTotalList.append(sobolNewTotal)
            counter = 0
            for index2 in range(len(indexVariable[index])):
                if currentIndexBool[index2]==True:
                    MainIndexOut[indexVariable[index][index2]] = sobolNewFirstList[counter]
                    TotalIndexOut[indexVariable[index][index2]] = sobolNewTotalList[counter]
                    counter = counter + 1
    time2 = time.time()
    print 'time2',time2- time0
    print 'M',MainIndexOut
    print 'T',TotalIndexOut


def funFix(myfun,X,indexArray):
    X = np.array(X)
    return myfun(X[indexArray])











    
def setInputRandomVector(inList):
    #setting the input random vector X
    myCollection = DistributionCollection(len(inList))
    for index in range(len(inList)):
        myCollection[index] = Distribution(inList[index])
    myDistribution = ComposedDistribution(myCollection)
    VectX = RandomVector(Distribution(myDistribution))
    return myDistribution,VectX
