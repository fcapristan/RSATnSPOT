import numpy as np
from openturns import *
import openTurnsTools as OTtools
import OTsobolWrapper as OTwrap
import time

def selectImportant(firstOrderIndices,threshold =.05):
    firstOrderIndices = np.array(firstOrderIndices)
    firstOrderIndices[firstOrderIndices<0.] = 0.0
    maxVal = max(firstOrderIndices)
    boolList = [False]*len(firstOrderIndices)
    for index in range(len(firstOrderIndices)):
        ratio = firstOrderIndices[index]/maxVal

        if ratio>=threshold:
            boolList[index] = True

    return np.array(boolList)
def sobolFunc(myfunction=None,nCases=None,dim=None):
    # time of failure is already defined in Physicsal model

    inList = [Uniform(0.0,1.0)]*dim # making a list of uniform distributions of dimension "dim"

    myDistribution,VectX = OTtools.setInputRandomVector(inList)



    #f = lambda X : EcCalc(X,tfail)
    
    #model = PythonFunction(dim,1,f)
  
    func = OTtools.myfunction(dim=dim,fun=myfunction,multiPro=12)
    model = NumericalMathFunction(func)


    RandomGenerator.SetSeed(10)
    myRandomExp1 = LHSExperiment(myDistribution,nCases)
    inp1 = myRandomExp1.generate()

    RandomGenerator.SetSeed(75)
    myRandomExp2 = LHSExperiment(myDistribution,nCases)
    inp2 = myRandomExp2.generate()

    #myRandomExp3 = LHSExperiment(myDistribution,nCases)
    #inp3 = myRandomExp3.generate()
    
    try:
        sensAnalysis = mySensitivityAnalysis(inp1,inp2,model)

        sensAnalysis.setBlockSize(nCases)
        print '########################'
        print 'Starting Computations '
        #secondOrderIndices = sensAnalysis.getSecondOrderIndices()

        firstOrderIndices = sensAnalysis.getFirstOrderIndices()
        totalOrderIndices = sensAnalysis.getTotalOrderIndices()
        print 'Done withSobol'
        out1 = sensAnalysis.getOutVals1()
        print 'Mean', out1.computeMean()
        print 'STD', out1.computeStandardDeviation()
        #out2 = sensAnalysis.getSecondValsOut()
        
        print 'First',firstOrderIndices
        print 'Total', totalOrderIndices
        #print 'Second',tfail,secondOrderIndices

        return [firstOrderIndices,totalOrderIndices]
    except Exception as inst:

        print inst.args
        print 'Check error message. It is possible that all monte carlo runs return 0 casualties'
        return [[],[]]




def fTotal(X,indexArray,funList,weightArray):
    Xvars = getXsingle(X,indexArray)
    valOut = 0.0
    for index in range(len(Xvars)):
        cval = funList[index](Xvars[index])

        valOut =  valOut + weightArray[index]*cval

    return valOut 

def getXsingle(X,indexList):
    Xvar = []
    X = np.array(X)

    for index in range(len(indexList)):

        Xvar.append(np.array(X[indexList[index]]))

    return Xvar


# the first 3 variables are shared among all funs
def fun1(X1):
    #time.sleep(.2)
    ret = np.array([0.0*X1[0]**2 + X1[1] + 0.8*X1[2] + X1[3]**3 + 0.1*X1[4]] + .4*X1[5])
    return ret



def fun2(X2):
    #time.sleep(.2)
    ret = np.array([(0.4*X2[0] + X2[1])**2 + 0.8*X2[2] + X2[3]**3 + 0.1*X2[4] + 0.1*X2[5]])
    return ret


def fun3(X3):
    #time.sleep(.2)
    ret = np.array([ 0.0*X3[0]**2 - X3[1] + (0.8*X3[2] + X3[3]**3 + 0.1*X3[4])**2])
    return ret

def fun4(X4):
    #time.sleep(.2)
    ret = np.array([ 0.1*X4[0]**1.3 - X4[1] + (0.008*X4[2] + X4[3]**3 + 0.0*X4[4])**2])
    return ret

def getIndexArray():
    sharedVals = 3
    funDim = [3,3,2,2]
    arrayOut = []
    offset = sharedVals

    sharedArray = np.arange(sharedVals)
    for index in range(len(funDim)):
        arrayOut.append(np.concatenate((sharedArray, offset + np.arange(funDim[index])),1))
        offset = arrayOut[index][-1] +1
    return arrayOut,sharedArray

indexArray,sharedArray = getIndexArray()

#print indexArray
#print sharedArray
#exit()
weightArray = [1.,1.,0.,0.]
f = lambda X:fTotal(X,indexArray,[fun1,fun2,fun3,fun4],weightArray)
dims = indexArray[-1][-1] + 1




nCases = 1000
print 'Stating Timer'
t0 = time.time()
[firstOrderIndices,totalOrderIndices] = sobolFunc(myfunction=f,nCases=nCases,dim=dims)
tf =  time.time()
print 'time',tf - t0

importantIndex = selectImportant(firstOrderIndices)
print importantIndex
print indexArray
nCases = 2000

#exit()


OTwrap.sobolMod(weightList = weightArray,functionList=[fun1,fun2,fun3,fun4],nCases=nCases,sharedIndex = sharedArray,
                ImportantBoolIndex = importantIndex , indexVariable=indexArray, seedList = [1,2,3],multiPro=12)                             






















