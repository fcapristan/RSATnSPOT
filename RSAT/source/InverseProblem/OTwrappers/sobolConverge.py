import numpy as np

from openturns import *
import openTurnsTools as OTtools
import OTsobolWrapper as OTwrap
import time
import scipy.io
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
def sobolFunc(myfunction=None,nCases=None,dim=None,seedList = None,multiPro=12):
    # time of failure is already defined in Physicsal model

    inList = [Uniform(0.0,1.0)]*dim # making a list of uniform distributions of dimension "dim"

    myDistribution,VectX = OTtools.setInputRandomVector(inList)



    #f = lambda X : EcCalc(X,tfail)
    
    #model = PythonFunction(dim,1,f)
  
    func = OTtools.myfunction(dim=dim,fun=myfunction,multiPro=multiPro)
    model = NumericalMathFunction(func)


    RandomGenerator.SetSeed(seedList[0])
    myRandomExp1 = LHSExperiment(myDistribution,nCases)
    inp1 = myRandomExp1.generate()

    RandomGenerator.SetSeed(seedList[1])
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

def gfun(Xi,ai):
    gi = (abs(4.*Xi - 2) + ai)/(1. + ai)
    return gi

def fSobol(Xarr,a_arr):
    assert(len(Xarr)==len(a_arr))
    ftemp = 1.0
    for i in range(len(Xarr)):
        ftemp = ftemp * gfun(Xarr[i],a_arr[i])
    return np.array([ftemp])


def sobolTheoretical(a_arr):
    Si = np.zeros(len(a_arr))
    for index in range(len(a_arr)):
        Si[index] = 1./(3.*(1.+a_arr[index])**2)
    temp = 1.0    
    for index in range(len(a_arr)):
        temp = temp*(1.+Si[index])
    Vy = -1 + temp

    Si = Si/Vy
    return Si,Vy

def sobolTheoreticalAll(Si_List,Var_List,weightList):
    VarTotal = np.sum(Var_List)
    SiFinal = []
    for index in range(len(weightList)):
        SiFinal.append((Si_List[index] * Var_List[index] * (weightList[index]**2))/VarTotal)
    return SiFinal




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
    sharedVals = 0
    funDim = [6,5,4]
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

runCase= ['simple','complex','total']
runCaseOption = runCase[1]



weights = [1.0,2.0,0.5]
a1_arr = np.array([78.,12.,0.5,2.,97.,33.])
a2_arr = np.array([20,12.,10,2.,9.])
a3_arr = np.array([45,0.5,97.,33.])
indexArray,sharedArray = getIndexArray()

aList = [a1_arr,a2_arr,a3_arr]
f1 = lambda X:fSobol(X,a1_arr)
f2 = lambda X:fSobol(X,a2_arr)
f3 = lambda X:fSobol(X,a3_arr)
dims = len(a1_arr) + len(a2_arr) + len(a3_arr)
nBoot = 1000
MainSiMatrix = np.zeros((nBoot,dims))
TotalSiMatrix = np.zeros((nBoot,dims))

for indBoot in range(nBoot):


    if runCaseOption=='simple':



        a_arr = np.array([78.,12.,0.5,2.,97.,33.])


        f = lambda X:fSobol(X,a_arr)
        dims = len(a_arr)



        nCases = 500000
        print 'Stating Timer'
        t0 = time.time()
        [firstOrderIndices,totalOrderIndices] = sobolFunc(myfunction=f,nCases=nCases,dim=dims)
        tf =  time.time()



        print 'time',tf - t0
        Si,V = sobolTheoretical(a_arr)
        print Si

    elif runCaseOption=='complex':


        nCases = 21000# 17000

        importantIndex = np.array([True]*(len(a1_arr) + len(a2_arr) + len(a3_arr)))
        
        seedList = [indBoot,2*indBoot+3,2*indBoot+5]
        firstOrderIndices,totalIndices = OTwrap.sobolMod(weightList=weights,functionList=[f1,f2,f3],nCases=nCases,sharedIndex=sharedArray,ImportantBoolIndex = importantIndex,indexVariable=indexArray,seedList=seedList,multiPro=12) 
        '''
        Si_List = []
        Var_List = []
    
        for index in range(len(weights)):
            Si,V = sobolTheoretical(aList[index])
            Si_List.append(Si)
            Var_List.append(V)
        SiTotal =  sobolTheoreticalAll(Si_List,Var_List,weights)
        '''
        MainSiMatrix[indBoot,:] = firstOrderIndices
        TotalSiMatrix[indBoot,:] = totalIndices
        #scipy.io.savemat('sobolMod.mat',dict(mainMod=MainSiMatrix,totalMod=TotalSiMatrix))
    elif runCaseOption=='total':
        nCases = 21000
        funList = [f1,f2,f3]
        fsimp = lambda X:fTotal(X,indexArray,funList,weights)

        seedList = [indBoot+1,2*indBoot+2]
        [firstOrderIndices,totalOrderIndices] = sobolFunc(myfunction=fsimp,nCases=nCases,dim=dims,seedList=seedList,multiPro=12)
        MainSiMatrix[indBoot,:] = firstOrderIndices
        TotalSiMatrix[indBoot,:] = totalOrderIndices
        #scipy.io.savemat('sobolReg.mat',dict(mainReg=MainSiMatrix,totalReg=TotalSiMatrix))
        '''
        Si_List=[];
        Var_List=[];
        for index in range(len(weights)):
            Si,V = sobolTheoretical(aList[index])
            Si_List.append(Si)
            Var_List.append(V)
        SiTotal =  sobolTheoreticalAll(Si_List,Var_List,weights)
        print 'SS',SiTotal
        exit()
        '''
        print 'Main100000',MainSiMatrix
        print 'Total',TotalSiMatrix
        exit()
if runCaseOption=='complex':
    scipy.io.savemat('sobolMod.mat',dict(mainMod=MainSiMatrix,totalMod=TotalSiMatrix))
elif runCaseOption=='total':
    scipy.io.savemat('sobolReg.mat',dict(mainReg=MainSiMatrix,totalReg=TotalSiMatrix))


print 'MM',MainSiMatrix
print 'TT',TotalSiMatrix
'''
importantIndex = selectImportant(firstOrderIndices)
print importantIndex

nCases = 2000

OTwrap.sobolMod(weightList = weightArray,functionList=[fun1,fun2,fun3,fun4],nCases=nCases,sharedIndex = sharedArray,
                ImportantBoolIndex = importantIndex , indexVariable=indexArray, seedList = [1,2,3],multiPro=12)                             

'''




















