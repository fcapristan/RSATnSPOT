import numpy as np
from openturns import *
import pathos.multiprocessing as mp

#===================================================================
#                         Marginal PDFs
#===================================================================

def setInputRandomVector(inList):
    #setting the input random vector X
    myCollection = DistributionCollection(len(inList))
    for index in range(len(inList)):
        myCollection[index] = Distribution(inList[index])
    myDistribution = ComposedDistribution(myCollection)
    VectX = RandomVector(Distribution(myDistribution))
    return myDistribution,VectX

def divideInputs(inList,importantIndex,subsetIndex=False):
    # importantIndex is a list of True or False variables
    # True if variable is important
    # False if is not
    # True indexes will be sampled separately
    # False variables will be a bunched together

    inList = np.array(inList)
    assert(len(inList)==len(importantIndex))

    impVars = inList[importantIndex]
    otherVars = inList[~importantIndex]
    if subsetIndex==True:

        #impVars.append(Uniform(0.0,1.0))
        # last marginal will be updated to just be indexes. That distribution does not really matter
        impVars = np.concatenate(([Uniform(-99999.0001,-99999.)],impVars),axis=0)
    impDist,impVectX = setInputRandomVector(impVars)
    otherDist,otherVectX = setInputRandomVector(otherVars)

    return impDist,impVectX,otherDist,otherVectX

def wrapVariables(Ximp,Xother,importantIndex):
    X = np.zeros((len(importantIndex)))
    X[importantIndex] = Ximp
    X[~importantIndex] = Xother
    return X

def LHSExperimentWithIndexForSobol(inList,importantIndex,seed1,seed2,seed3,nCases):

    impDist,impVectX,otherDist,otherVectX = divideInputs(inList,importantIndex,subsetIndex=True) # dividing variables and setting them in openturns format
    RandomGenerator.SetSeed(seed1)    
    myRandomExp1 = LHSExperiment(impDist,nCases)
    inp1 = myRandomExp1.generate()
    RandomGenerator.SetSeed(seed2)    
    myRandomExp2 = LHSExperiment(impDist,nCases)
    inp2 = myRandomExp2.generate()


    for index in range(nCases) : #couldn't figure out how to do it with regular numpy assignments
        inp1[index,0] = index# assigning the indexes. THis is not really a variable, but a set of variables.
        inp2[index,0] = index + nCases
    RandomGenerator.SetSeed(seed3)
    otherRandomExp = LHSExperiment(otherDist,2*nCases)
    inpOther = otherRandomExp.generate()
    return inp1,inp2,inpOther

def getXfromSubset(Ximp,VarOther,importantIndex):
    # index 0 of Ximp is just the index for the subset vars
    Xindex = int(Ximp[0])
    assert (Xindex==Ximp[0]) #making sure Ximp is an integer, it is simply expressed as a double
    Ximp = np.array(Ximp)
    Ximp = Ximp[1:len(Ximp)] # keeping the actual variables
    Xother = VarOther[Xindex,:]
    X = wrapVariables(Ximp,Xother,importantIndex)
    return X

def getValOutfromSubset(Ximp,VarOther,importantIndex,fun):
    X = getXfromSubset(Ximp,VarOther,importantIndex)
    valOut = fun(X)
    return valOut
# wrapper function for running samples in parallel
class myfunction(OpenTURNSPythonFunction):
    def __init__(self,multiPro=12,dim=None,fun=None):
        OpenTURNSPythonFunction.__init__(self, dim, 1)
        self.Nproc = multiPro
        self.myfun = fun
    def _exec(self,in_point):
        return (self._exec_sample([in_point]))

    def _exec_sample(self, in_sample):
        Xlist=[[s for s in p] for p in in_sample]
        p = mp.Pool(self.Nproc)
        myfun = self.myfun 
        samp = p.map(myfun,Xlist)
        return samp



    



 
        
