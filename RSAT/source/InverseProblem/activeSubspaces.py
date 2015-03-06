#routines to calculate active subspaces

#import derivs
import numpy as np
import time
import os
import copy

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
    X = np.array(X)
    #print X
    #fun-< function, X->vector of inputs,step->step vector for finite diff
    f0 = fun(X)
    n = len(X)
    J = np.zeros((n,1))

    for index in range(n):
        Xnew = copy.copy(X)
        Xnew[index] = X[index] + step[index]
        fp1 = fun(Xnew)
        J[index,0] = (fp1-f0)/step[index]

    return J,f0

def getJacobianMultiOut(fun,X,step,ndimOut=1):
    X = np.array(X)
    #print X
    #fun-< function, X->vector of inputs,step->step vector for finite diff
    print 'Starting Jacobian',len(X)
    f0_array = fun(X)
    f0_array = np.array(f0_array) # making sure it is a numpy array
    assert(len(f0_array)==ndimOut)
    n = len(X)
    J_array = np.zeros((n,ndimOut))# [np.zeros((n,1))]*ndimOut
    for index in range(n):
        Xnew = copy.copy(X)
        Xnew[index] = X[index] + step[index]
        fp1_array = np.array(fun(Xnew))
        
        J_array[index,:] = (fp1_array -f0_array)/step[index]
    print J_array
    return J_array,f0_array




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
            print 'Warning: New Sample Value exceeds limits [-1,1]' 
            print 'This might be OK if finite difference being used on the bounds'
            print 'Current values should be [-1,1] but are  ',activeSubspaceVar.max(),activeSubspaceVar.min()
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
            threshold
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
    def getSubspace(self,fun=None,funPrime=None,distributionClass=None,nSamples=None,nsubspace=None,step=None,controlIndex=None,multiPro=1,threshold=.95,LHS=False):
        ndim = len(distributionClass.name)
        self.ndim = ndim
        
        if LHS ==True:
            self.base = sampling.getLHSuniform(ndim,nSamples)
        elif LHS==False:
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
            step = 0.01*np.ones((ndim))
        self.step = step
        #nondim  and rescaled desired values new vars within [-1,1]
        
        A = np.diag(.5*deltaVec)
        self.AJA = A #storing A matrix to be used in computations outside this routine (giving the user more control)
        #calculating Jacobian
        sumJ = 0.0
        self.outputs = np.zeros((nSamples))-999999.
        self.J = np.zeros((ndim,nSamples))
        if funPrime==None:
            funPrime = lambda X:getJacobian(fun,X,step)
        
        
        if multiPro==1:
            for index in range(nSamples):
                
                J,funcOut = funPrime(self.values[:,index])
                self.J[:,index] = J[:,0]
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
                self.J[:,index] = Jlocal[:,0]
        
        
        sumJ = np.dot(A,sumJ)#could change this to make matrix calc faster (elem*row) 
        sumJ = np.dot(sumJ,A)
        
        sumJ = (1./np.float(nSamples))*sumJ
        eigenvals,eigenvecs = np.linalg.eigh(sumJ)
        
        idx = eigenvals.argsort()[::-1]   # sorting from highest to lowest
        eigenvals = eigenvals[idx]
        eigenvecs = eigenvecs[:,idx]  
        self.eigenvals = eigenvals
        
        if nsubspace==None: #let the code pick the number of eigenvalues
            self.selectSubspaces(threshold=threshold)
        else:
            self.Nsubspace = nsubspace
        self.eigenvecs = eigenvecs
        self.U = eigenvecs[:,0:self.Nsubspace]
        self.subspaceBase = np.dot((self.U).T,self.base)
        self.importantVariables()

    def updateSubspace(self,nsubspace=None,threshold=None):
        if nsubspace==None:
            if threshold==None:
                print 'Error in updateSubspace. No specified nsubspace or threshold'
                exit()
            self.selectSubspaces(threshold=threshold)
        else :
            if threshold!=None:
                print 'Error in updateSubspace. Cannot specify both nsubspace and threshold'
                exit()
            self.Nsubspace = nsubspace

        self.U = self.eigenvecs[:,0:self.Nsubspace]
        self.subspaceBase = np.dot((self.U).T,self.base)
        self.importantVariables()
    
    def selectSubspaces(self,threshold=None):
        eigenvals = self.eigenvals
        #eigrnvalues must be in decreasing order
        totalSum = np.sum(eigenvals)
        currSum = 0.0
        
        for index in range(len(eigenvals)):
            currSum = currSum + eigenvals[index]
            #print 'temp',currSum/totalSum,index
            
            if (currSum/totalSum)>=threshold:
                #print 'Final',currSum/totalSum,index
                break
        self.Nsubspace = index + 1
    
    
    def importantVariables(self):#Any other method such as 
        n = self.Nsubspace
        eigenvals = self.eigenvals[0:n]
        #self.values = np.zeros((ndim,nSamples))
        ndim , nSamples = np.shape(self.values)
        #ndim = len(distributionClass.name)
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



def addActiveSubspace_multiOuput(active1,active2,threshold=.99,nsubspace=None):
    
    ch1 = active1.ndim == active2.ndim
    ch2 = (active1.AJA == active2.AJA).all()

    assert ((ch1==ch2)==True)
    '''
    if len(resMulti[index])>2:
            otherVals = resMulti[index][2]
            otherList.append(otherVals)
    '''

    # adding classes   
    activeOut = copy.copy(active1)
    funDim = activeOut.funDim
    activeOut.values = np.concatenate((active1.values,active2.values),axis=1)
    activeOut.base = np.concatenate((active1.base,active2.base),axis=1)
    A = activeOut.AJA

    activeOut.outputs   
    nSamples =  len(active1.outputs[0]) + len(active2.outputs[0])
    activeOut.outputs = np.zeros((funDim,nSamples))
    ndim = active1.ndim
    activeOut.J = np.zeros((funDim,ndim,nSamples))
    for index in range(activeOut.funDim):

        activeOut.outputs[index,:] = np.concatenate((active1.outputs[index],active2.outputs[index]),axis=1)
        activeOut.J[index,:,:] = np.concatenate((active1.J[index,:,:],active2.J[index,:,:]),axis=1)


    
    sumJ = np.zeros((funDim,ndim,ndim))
    for indexDim in range(funDim):
        for index in range(nSamples):
            J = activeOut.J[indexDim,:,index,None]# resMulti[index][0]
            sumJ[indexDim,:,:] = np.dot(J,J.T) + sumJ[indexDim,:,:]


    activeOut.eigenvals = []#np.zeros((ndim,funDim))
    activeOut.eigenvecs = []#np.zeros((funDim,ndim,ndim))
    activeOut.Nsubspace = [0]*funDim
    activeOut.U = []
    activeOut.subspaceBase = []



    for index in range(funDim):

        sumJ[index,:,:] = np.dot(A,sumJ[index,:,:])#could change this to make matrix calc faster (elem*row) 
        sumJ[index,:,:] = np.dot(sumJ[index,:,:],A)
    
        sumJ[index,:,:] = (1./np.float(nSamples))*sumJ[index,:,:]
        eigenvals,eigenvecs = np.linalg.eigh(sumJ[index,:,:])
        idx = eigenvals.argsort()[::-1]   # sorting from highest to lowest
        eigenvals = eigenvals[idx]
        eigenvecs = eigenvecs[:,idx]  
        activeOut.eigenvals.append(eigenvals)# [:,index] = eigenvals
        activeOut.eigenvecs.append(eigenvecs)# [index,:,:] = eigenvecs

        if nsubspace==None: #let the code pick the number of eigenvalues
            nsubspace = activeOut.selectSubspaces(eigenvals,threshold=threshold)

        activeOut.Nsubspace[index] = nsubspace
        activeOut.U.append(eigenvecs[:,0:nsubspace])
        activeOut.subspaceBase.append(np.dot((activeOut.U[index]).T,activeOut.base))
    return activeOut

class activeSubspace_multiOutput:
    def __init__(self):
        self.eigenvals = []
        self.eigenvecs = []
        self.Nsubspace = []
        self.outputs = []
    def getSubspace(self,fun=None,funPrime=None,funDim=1,distributionClass=None,nSamples=None,nsubspace=None,step=None,controlIndex=None,multiPro=1,threshold=.95,LHS=False):
        ndim = len(distributionClass.name)
        self.ndim = ndim
        
        if LHS ==True:
            self.base = sampling.getLHSuniform(ndim,nSamples)
        elif LHS==False:
            self.base = np.random.uniform(size=[ndim,nSamples])
        self.values = np.zeros((ndim,nSamples))
        self.controlIndex = controlIndex
        self.funDim = funDim
        lowLim = distributionClass.minVal
        upLim = distributionClass.maxVal
        deltaVec = upLim - lowLim
        
        for index in range(ndim):
            self.values[index,:] = self.base[index,:]*(deltaVec[index])+lowLim[index]
        self.base = 2.0*self.base-1. # adjusting bounds so that it is bounded by [-1,1] instead of [0,1]
        
        #vars is a list of input vectors
        if step==None:
            step = 0.01*np.ones((ndim))
        self.step = step
        #nondim  and rescaled desired values new vars within [-1,1]
        
        A = np.diag(.5*deltaVec)
        self.AJA = A #storing A matrix to be used in computations outside this routine (giving the user more control)
        #calculating Jacobian
        sumJ = np.zeros((funDim,ndim,ndim))
        self.outputs = np.zeros((funDim,nSamples))-999999.
        self.J = np.zeros((funDim,ndim,nSamples))
        otherList = []
        if funPrime==None:
            funPrime = lambda X:getJacobianMultiOut(fun,X,step,ndimOut=funDim)# getJacobian(fun,X,step)
        
        
        if multiPro==1:
            for index in range(nSamples):
                #J,funcOut = funPrime(self.values[:,index])
                resMulti = funPrime(self.values[:,index])
                J = resMulti[0]
                funcOut = resMulti[1]
                if len(resMulti)>2:
                    otherVals = resMulti[2]
                    otherList.append(otherVals)
                for indexDim in range(funDim):
                    Jtemp = J[:,indexDim,None]                
                    self.J[indexDim,:,index] = Jtemp[:,0]
                    self.outputs[indexDim,index] = funcOut[indexDim]
                    sumJ[indexDim,:,:] = np.dot(Jtemp,Jtemp.T) + sumJ[indexDim,:,:]
        else:
            import pathos.multiprocessing as mp
            valsMulti = ((self.values).T).tolist()
            p = mp.Pool(multiPro)
            resMulti= p.map(funPrime,valsMulti)

            for index in range(nSamples):
                J = resMulti[index][0]
                funcOut = resMulti[index][1]
                if len(resMulti[index])>2:
                    otherVals = resMulti[index][2]
                    otherList.append(otherVals)

                for indexDim in range(funDim):
                    Jtemp = J[:,indexDim,None]                
                    self.J[indexDim,:,index] = Jtemp[:,0]

                    self.outputs[indexDim,index] = funcOut[indexDim]
                    sumJ[indexDim,:,:] = np.dot(Jtemp,Jtemp.T) + sumJ[indexDim,:,:]
        self.otherVals = otherList
        self.eigenvals = []#np.zeros((ndim,funDim))
        self.eigenvecs = []#np.zeros((funDim,ndim,ndim))
        self.Nsubspace = [0]*funDim
        self.U = []
        self.subspaceBase = []
        for index in range(funDim):

            sumJ[index,:,:] = np.dot(A,sumJ[index,:,:])#could change this to make matrix calc faster (elem*row) 
            sumJ[index,:,:] = np.dot(sumJ[index,:,:],A)
        
            sumJ[index,:,:] = (1./np.float(nSamples))*sumJ[index,:,:]


            eigenvals,eigenvecs = np.linalg.eigh(sumJ[index,:,:])
        
            idx = eigenvals.argsort()[::-1]   # sorting from highest to lowest
            eigenvals = eigenvals[idx]
            eigenvecs = eigenvecs[:,idx]  
            self.eigenvals.append(eigenvals)# [:,index] = eigenvals
            self.eigenvecs.append(eigenvecs)# [index,:,:] = eigenvecs

            if nsubspace==None: #let the code pick the number of eigenvalues
                nsubspace = self.selectSubspaces(eigenvals,threshold=threshold)

            self.Nsubspace[index] = nsubspace

            self.U.append(eigenvecs[:,0:nsubspace])
            self.subspaceBase.append(np.dot((self.U[index]).T,self.base))
            

    def updateSubspace(self,nsubspace=None,threshold=None):
        self.U = []
        self.subspaceBase = []
        self.Nsubspace = [0]*self.funDim
        if nsubspace==None:
            if threshold==None:
                print 'Error in updateSubspace. No specified nsubspace or threshold'
                exit()
            for index in range(self.funDim):
                eigenvals = self.eigenvals[index]
                eigenvecs = self.eigenvecs[index]
                nsubspace = self.selectSubspaces(eigenvals,threshold=threshold)
                self.Nsubspace[index] = nsubspace

                self.U.append(eigenvecs[:,0:nsubspace])
                self.subspaceBase.append(np.dot((self.U[index]).T,self.base))

        else :
            if threshold!=None:
                print 'Error in updateSubspace. Cannot specify both nsubspace and threshold'
                exit()
            for index in range(self.funDim):
                eigenvecs = self.eigenvecs[index]
                self.Nsubspace[index] = nsubspace
                self.U.append(eigenvecs[:,0:nsubspace])
                self.subspaceBase.append(np.dot((self.U[index]).T,self.base))

    
    def selectSubspaces(self,eigenvals,threshold=None):
        #eigenvals = self.eigenvals
        #eigrnvalues must be in decreasing order
        totalSum = np.sum(eigenvals)
        currSum = 0.0
        
        for index in range(len(eigenvals)):
            currSum = currSum + eigenvals[index]
            #print 'temp',currSum/totalSum,index
            
            if (currSum/totalSum)>=threshold:
                #print 'Final',currSum/totalSum,index
                break
        Nsubspace = index + 1
        return Nsubspace
    
    def selectSubspacesWeighted(self,fraction=None,indexDim=None,nsubspace=None,threshold=None):
        assert(0.<=fraction<=1.)
        #first determine max and min values of Y outputs
        self.Uweighted = [0]*self.funDim
        self.subspaceBaseWeighted = [0]*self.funDim
        ymax = np.max(self.outputs[indexDim,:])
        ymin = np.min(self.outputs[indexDim,:])
        deltay = ymax - ymin
        ynormalized = (self.outputs[indexDim,:] - ymin)/deltay # outputs now range from 0 to 1
        
        #get index of values by using fraction to divide
    
        indexDesired = ynormalized>=fraction
        NweightedSamples = np.sum(indexDesired)
        nSamples = len(self.outputs[indexDim,:])
        print 'nn',NweightedSamples
        weight =  float(NweightedSamples)/(float(nSamples)-float(NweightedSamples) )  # simple ratio of samples
        print 'w',weight
        #recalculate Covariance matrix by adding weight
        sumJweighted = 0.0
        for index in range(nSamples):
            J =  np.reshape(self.J[indexDim,:,index],(self.ndim,1))# not sure of NONE!!!!!!! NEED TO CHECK THIS!!!!!!!!!!!!!!!!!!!!1
            #print 'J has a shape of', np.shape(J)
            if indexDesired[index]==True:
                sumJweighted = np.dot(weight*J,weight*J.T) + sumJweighted
                #print J
            else:
                sumJweighted = 0* np.dot(J,J.T) + sumJweighted
            #print np.dot(J,J.T)
        A = self.AJA
        sumJweighted = np.dot(A,sumJweighted)
        sumJweighted = np.dot(sumJweighted,A)
        sumJweighted = (1./np.float(nSamples))*sumJweighted
        
        eigenvals,eigenvecs = np.linalg.eigh(sumJweighted)

        idx = eigenvals.argsort()[::-1]   # sorting from highest to lowest
        eigenvals = eigenvals[idx]
        eigenvecs = eigenvecs[:,idx]  
        #print eigenvals
        #print self.eigenvals[0]
        #self.eigenvalsWeig.append(eigenvals)# [:,index] = eigenvals
        #self.eigenvecs.append(eigenvecs)# [index,:,:] = eigenvecs
        
        if nsubspace==None: #let the code pick the number of eigenvalues
            nsubspace = self.selectSubspaces(eigenvals,threshold=threshold)
        
        #self.Nsubspace[index] = nsubspace
        
        #self.Uweighted[indexDim] = eigenvecs[:,0:nsubspace]
        #self.subspaceBaseWeighted[indexDim] = np.dot((self.U[index]).T,self.base)
        self.UW = [0]*self.funDim
        self.subspaceBaseW = [0]*self.funDim
        self.UW[indexDim] = eigenvecs[:,0:nsubspace]
        #print np.shape(self.U[indexDim])
        #self.subspaceBase = [0]*self.funDim
        
        self.subspaceBaseW[indexDim] = np.dot((self.UW[indexDim]).T,self.base)
    
    '''
    
    def importantVariables(self):#Any other method such as 
        n = self.Nsubspace
        eigenvals = self.eigenvals[0:n]
        #self.values = np.zeros((ndim,nSamples))
        ndim , nSamples = np.shape(self.values)
        #ndim = len(distributionClass.name)
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

    '''


class ActiveSubspacesLinear:
    def __init__(self,Xmatrix,Jmatrix,rangeVar,Amatrix,threshold = 0.99 ):
        #Xmatrix is the subspace.base [-1,1]
        Jsum = [0]*len(rangeVar)
        AJA = [0]*len(rangeVar)
        UList = [0]*len(rangeVar)
        ndim , nSamples = np.shape(Jmatrix)   
        Adiag = np.diag(Amatrix)  
        subspaceBaseList=[]
        for in0 in range(len(rangeVar)):   
            for in1 in range(nSamples):
                Jtemp = decomposeVariable(Jmatrix[:,in1],rangeVar,in0,flag_2d=True)
                Jsum[in0] = np.dot(Jtemp,Jtemp.T) + Jsum[in0]
            A = decomposeVariable(Adiag,rangeVar,in0)
            A = np.diag(A)
            AJA[in0] = np.dot(A,Jsum[in0])
            AJA[in0] = np.dot(AJA[in0],A)
            eigenvals,eigenvecs = np.linalg.eigh(AJA[in0])
            idx = eigenvals.argsort()[::-1]   # sorting from highest to lowest
            eigenvals = eigenvals[idx]
            eigenvecs = eigenvecs[:,idx]  

            Nsubspace = selectSubspaces(eigenvals,threshold = threshold)
            print eigenvals
            UList[in0] = eigenvecs[:,0:Nsubspace]
            subspaceBaseList.append(np.dot((UList[in0]).T,decomposeVariableMatrix(Xmatrix,rangeVar,in0)))
        self.UList = UList
        self.subspaceBaseList = subspaceBaseList



    def updateSubspace(self,nsubspace=None,threshold=None):
        if nsubspace==None:
            if threshold==None:
                print 'Error in updateSubspace. No specified nsubspace or threshold'
                exit()
            self.selectSubspaces(threshold=threshold)
        else :
            if threshold!=None:
                print 'Error in updateSubspace. Cannot specify both nsubspace and threshold'
                exit()
            self.Nsubspace = nsubspace

        self.U = self.eigenvecs[:,0:self.Nsubspace]
        self.subspaceBase = np.dot((self.U).T,self.base)
        self.importantVariables()
    
def selectSubspaces(eigenvals,threshold=None):
    #eigrnvalues must be in decreasing order
    totalSum = np.sum(eigenvals)
    currSum = 0.0
    
    for index in range(len(eigenvals)):
        currSum = currSum + eigenvals[index]
        #print 'temp',currSum/totalSum,index
        
        if (currSum/totalSum)>=threshold:
            #print 'Final',currSum/totalSum,index
            break
    return index + 1




def decomposeVariable(J_X,rangeVar,index,flag_2d = False):
    #J_X input vector
    # rangeVar is a list of indexes that contain the required grouping for inputs
    if flag_2d==False:
        return J_X[rangeVar[index]]
    else:
        return J_X[rangeVar[index],None]

def decomposeVariableMatrix(J_X,rangeVar,index):
    return J_X[rangeVar[index],:]





