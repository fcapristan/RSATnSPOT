#routines to calculate active subspaces

#import derivs
import numpy as np
import time
import os

def normalizeX(distribution,X):
    ndim,nSamples = np.shape(X)
    Xnormalized = np.zeros((ndim,nSamples))
    deltaVec = np.zeros((ndim))
    
    for index in range(ndim):
        delta = distribution.upperLimit[index] - distribution.lowerLimit[index]
        deltaVec[index] = delta
        Xnormalized[index,:] = 2.0*(X[index,:] - distribution.lowerLimit[index])/delta -1.0
    return Xnormalized




def getJacobian(fun,X,step):
    import copy
    X = np.array(X)
    #print X
    #fun-< function, X->vector of inputs,step->step vector for finite diff
    f0 = fun(X)
    n = len(X)
    J = np.zeros((n,1))
    #print 'Xor',X
    for index in range(n):
        Xnew = copy.copy(X)
        Xnew[index] = X[index] + step[index]
        fp1 = fun(Xnew)
        J[index,0] = (fp1-f0)/step[index]
    
    
    return J,f0


def steps(end,n):
    n=n+1
    start = 0
    if n<2:
        raise Exception("behaviour not defined for n<2")
    step = (end-start)/float(n-1)
    return [int(round(start+x*step)) for x in range(n)]




def parallelHelper(fun,Xmat,procNum):
    os.chdir('./')
    if os.path.isfile('activeSubThread'+str(procNum)+'.npy'):
        os.remove('activeSubThread'+str(procNum)+'.npy')
    ndim,nSamples = np.shape(Xmat)
    sumJ = 0.0
    for index in range(nSamples):
        J,f0 = fun(Xmat[:,index])
        sumJ = np.dot(J,J.T) + sumJ
    f0 = np.array(f0)
    print 'In Thread ',procNum,f0
    np.save('activeSubThread'+str(procNum),(sumJ,f0))

def checkParallelFiles(nproc):
    boolCheck = []  
    print os.getcwd()
    
    for index in range(nproc):
        boolCheck.append(os.path.isfile('activeSubThread'+str(index)+'.npy'))
    print boolCheck
    if np.sum(boolCheck)==nproc:
        return True
    else:
        return False




def simpleSampleRandom(dist,mu,sigma,baseVal=None):
    if dist.lower()=='uniform':
        fval = sampleUniform(mu,sigma,baseVal)
    elif dist.lower()=='normal':
        fval = sampleGaussian(mu,sigma,baseVal)
    return fval

def sampleGaussian(mu,sigma,baseVal = None):
    if baseVal==None:
        baseVal = np.random.normal(0.0,1.0)
    fval = baseVal*sigma + mu
    return fval

def sampleUniform(mu,sigma,baseVal=None):
    if baseVal==None:
        baseVal = np.random.uniform()
    b = mu + np.sqrt(3.0)*sigma
    a = 2.0*mu - b
    
    fval = (b-a)*baseVal + a
    return fval

def uniformVals(a,b):
    assert(b>a)
    mu = .5*(b+a)
    var = (1.0/12.0)*(b-a)**2.0
    sigma = (b-a)/(2.0*np.sqrt(3.))
    
    return mu,var,sigma





class varsClass:
    def __init__(self,distName=None):
        self.X = Xclass(uniform=0,distName=distName)
        self.Xuniform = Xclass(uniform=1,distName=distName)






class distribution:
    
    def __init__(self,name=None,mu=None,std=None,stdmin=None,stdmax=None):
        self.name = name
        self.mu = np.array(mu)
        self.std = np.array(std)
        self.dim = len(self.name)
        self.stdmax = np.array(stdmax)
        assert(len(self.name)==len(self.mu))
        assert(len(self.name)==len(self.std))
        
        self.minVal = np.zeros((len(self.name)))
        self.maxVal = np.zeros((len(self.name)))
        
        
        for index in range(len(self.name)):
            name =self.name[index]
            assert((name.lower()=='normal')or(name.lower()=='uniform'))
            if name.lower()=='uniform':
                b = self.mu[index] + np.sqrt(3.0)*self.stdmax[index]
                a = 2.0*self.mu[index] - b
                delta = b-a #positive value 
                self.maxVal[index] = b
                self.minVal[index] = a
            
            
            elif name.lower()=='normal':
                self.maxVal[index] = self.mu[index] + 7.0*self.stdmax[index]
                self.minVal[index] = self.mu[index] - 7.0*self.stdmax[index]





class samples:
    def __init__(self):
        self.uniform = []
        self.base = [] #generateNew(self,distributionClass,nSamples)
        self.values = []#sampleRandom(self)
    
    
    def generateNew(self,distributionClass,nSamples=None,uniform=0):
        self.nSamples = nSamples
        self.generateNewBase(distributionClass,nSamples,uniform)
        self.updateValues(distributionClass)
    
    def generateNewBase(self,distributionClass=None,nSamples=None,uniform=0):
        self.uniform = uniform
        n = len(distributionClass.name)
        baseVal = np.zeros((n,nSamples))
        
        if self.uniform==0:
            for index in range(n):
                dist = distributionClass.name[index]
                if dist.lower() =='uniform':
                    baseVal[index,:] = np.random.uniform(low=0.0,high=1.0,size=nSamples)
                elif dist.lower()=='normal':
                    baseVal[index,:] = np.random.normal(loc=0.0,scale=1.0,size=nSamples)
        else:
            for index in range(n):
                baseVal[index,:] = np.random.uniform(low=0.0,high=1.0,size=nSamples)
        self.base= baseVal
    
    def updateValues(self,distributionClass=None):
        row,col = np.shape(self.base)
        vals = np.zeros((row,col))
        for index in range(len(distributionClass.name)):
            dist = distributionClass.name[index]
            mu = distributionClass.mu[index]
            sigma = distributionClass.std[index]
            vals[index,:]=simpleSampleRandom(dist,mu,sigma,self.base[index,:])
        
        self.values= vals
    def updateValuesActiveSubspace(self,distributionClass,activeSubspaceClass):
        lowLim = distributionClass.minVal
        upLim = distributionClass.maxVal
        delta = upLim-lowLim
        activeSubspaceVar = np.zeros((distributionClass.dim,self.nSamples))
        values = self.values
        for index in range(len(distributionClass.name)):
            activeSubspaceVar[index,:] = 2.0*(values[index,:]-lowLim[index])/delta[index]-1.
        
        if activeSubspaceVar.max()>1.0 or activeSubspaceVar.min()<-1.0:
            print 'Error: New Sample Value exceeds limits [-1,1]' 
            exit()
        XS = np.dot((activeSubspaceClass.U).T,activeSubspaceVar)
        self.activeSubspaceVals = XS
    def updateOutputs(self,fun=None,multiPro=1):
        
        if multiPro==1:
            outputs = np.zeros((self.nSamples))
            for index in range(self.nSamples):
                X = self.values[:,index]
                outputs[index] = fun(X)
        else:
            import pathos.multiprocessing as mp
            
            valsMulti = ((self.values).T).tolist()
            p = mp.Pool(multiPro)
            outputs = p.map(fun,valsMulti)
            outputs = np.array(outputs)
        self.outputs = outputs



def localPool(X,funPrime,multiPro):
    valsMulti = (X.T).tolist()
    
    pool = Pool(processes=multiPro)
    J_vec,funcOut_list,sumJList = pool.map(funPrime,valsMulti)
    return J_vec,funcOut_list,sumJList



class activeSubspace:
    def __init__(self):
        self.eigenvals = []
        self.eigenvecs = []
        self.Nsubspace = []
        self.outputs = []
    def getSubspace(self,fun=None,funPrime=None,distributionClass=None,nSamples=None,nsubspace=None,step=None,controlIndex=None,multiPro=1):
        ndim = len(distributionClass.name)
        self.base = np.random.uniform(size=[ndim,nSamples])
        self.values = np.zeros((ndim,nSamples))
        self.controlIndex = controlIndex
        lowLim = distributionClass.minVal
        upLim = distributionClass.maxVal
        deltaVec = upLim - lowLim
        
        for index in range(ndim):
            self.values[index,:] = self.base[index,:]*(deltaVec[index])+lowLim[index]
        self.base = 2.0*self.base-1. # adjusting bounds so that it is bounded by [-1,1] instead of [0,1]
        
        #vars is a list of input vectors
        if step==None:
            step = 0.1*np.ones((ndim))
        self.step = step
        #nondim  and rescaled desired values new vars within [-1,1]
        
        A = np.diag(.5*deltaVec)
        #calculating Jacobian
        sumJ = 0.0
        self.outputs = np.zeros((nSamples))-999999.
        if funPrime==None:
            funPrime = lambda X:getJacobian(fun,X,step)
        
        
        if multiPro==1:
            for index in range(nSamples):
                
                J,funcOut = funPrime(self.values[:,index])
                self.outputs[index] = funcOut
                sumJ = np.dot(J,J.T) + sumJ
        else:
            import pathos.multiprocessing as mp
            valsMulti = ((self.values).T).tolist()
            p = mp.Pool(multiPro)
            resMulti= p.map(funPrime,valsMulti)
            for index in range(nSamples):
                Jlocal = resMulti[index][0]
                funcOut = resMulti[index][1]
                self.outputs[index] = funcOut
                sumJ = np.dot(Jlocal,Jlocal.T) + sumJ
        
        
        
        sumJ = np.dot(A,sumJ)#could change this to make matrix calc faster (elem*row) 
        sumJ = np.dot(sumJ,A)
        
        sumJ = (1./np.float(nSamples))*sumJ
        eigenvals,eigenvecs = np.linalg.eigh(sumJ)
        
        idx = eigenvals.argsort()[::-1]   # sorting from highest to lowest
        eigenvals = eigenvals[idx]
        eigenvecs = eigenvecs[:,idx]  
        self.eigenvals = eigenvals
        
        if nsubspace==None: #let the code pick the number of eigenvalues
            self.selectSubspaces()
        else:
            self.Nsubspace = nsubspace
        self.eigenvecs = eigenvecs
        self.U = eigenvecs[:,0:self.Nsubspace]
        self.subspaceBase = np.dot((self.U).T,self.base)
        self.importantVariables(distributionClass)
    
    def selectSubspaces(self):
        eigenvals = self.eigenvals
        #eigrnvalues must be in decreasing order
        totalSum = np.sum(eigenvals)
        currSum = 0.0
        threshold = 0.95
        for index in range(len(eigenvals)):
            currSum = currSum + eigenvals[index]
            #print 'temp',currSum/totalSum,index
            
            if (currSum/totalSum)>=threshold:
                #print 'Final',currSum/totalSum,index
                break
        self.Nsubspace = index + 1
    
    
    def importantVariables(self,distributionClass):#Any other method such as 
        n = self.Nsubspace
        eigenvals = self.eigenvals[0:n]
        ndim = len(distributionClass.name)
        eigenvecs = self.U
        rel_imp = 0.0
        for index in range(n):
            rel_imp = np.abs(eigenvecs[:,index])*eigenvals[index] + rel_imp
        self.variable_importance = rel_imp/rel_imp.max()
        imp_index = self.variable_importance>.05
        self.important_index = imp_index
        if self.controlIndex!=None:
            check_index = (np.array(imp_index)&np.array(self.controlIndex))
            if np.sum(check_index)==0:
                print 'Nothing to optimize'
            self.important_control_index = check_index
        
        self.important_var_index = imp_index 






