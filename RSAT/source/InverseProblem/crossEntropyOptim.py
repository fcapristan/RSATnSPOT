import numpy as np

# a modified cross entropy method suited for Ec calculations

def sample(bounds=None,dist=None,nsamples=1):
    #bounds contain a list with lower and upper bounds for all the variables
    # len(bounds) = number of variables. bounds[index] = [lower,upper] for variable in index (uses uniform dist)
    # dist is similar to bounds, but contains [mulist,stdlist] instead for use in a normal distribution
    X = [-999]*nsamples
    for index0 in range(nsamples):
        if dist==None:
            nvars = len(bounds)
            X0 = np.zeros((nvars))
            for index in range(nvars):
                X0[index] = np.random.uniform(low=bounds[index][0],high=bounds[index][1])
        else:
            assert(len(dist[0])==len(dist[1]))
            nvars = len(dist[0])
            muvals = dist[0]
            stdvals = dist[1]
            X0 = np.zeros((nvars))
            for index in range(nvars):
                if bounds!=None:
                    if bounds[index][0]==bounds[index][1]:
	               curval = bounds[index][0]
                    else:
                       curval = accrej(muvals[index],stdvals[index],bounds[index][0],bounds[index][1])
                else:
                    curval= np.random.normal(loc=muvals[index],scale=stdvals[index])
                X0[index] = curval
        X[index0] = X0
    return X


def accrej(mu,std,low,high):
    #simple acceptance rejection to sample within the given bounds
    while True:
        sampVal = np.random.normal(loc=mu,scale=std)
        if low<=sampVal<=high:
            break
    return sampVal


def optimMin(bounds=None,dist=None,fun=None,epsilonVar=0.0000001,nMax = 5000,nsamples=100,option='min'):
    import copy
    #this works better if everything is scaled appropriately
    varOut = np.inf
    Yold = np.array([])
    Xoldn = None
    counter = 1
    while (varOut > epsilonVar)and(counter<= nMax):
        X = sample(bounds=bounds,dist=dist,nsamples=nsamples) # samples X
        Xmat = np.array(X)
        Y = np.zeros((nsamples))
 
        for index in range(nsamples):
            Y[index] = fun(X[index])
            #print X[index],Y[index]
        #print Xmat,Xoldn
        Y = np.concatenate((Y,Yold))
        if Xoldn!=None:
            Xmat = np.concatenate((Xmat,Xoldn),axis=0)

        if option=='min':
            indexSorted = np.argsort(Y)
        elif option=='max':
            indexSorted =  np.argsort(-Y)
        else:
            print 'Option in optimMin does not exist'
            print 'It can only be min or max'
            exit()
        indexSorted = indexSorted[0:20]
        #print 'XN'
        #print Xmat
        Xn = Xmat[indexSorted]
        #print 'Xn',Xn
        newMean = np.mean(Xn,axis=0)
        newvar = np.var(Xn,axis=0)
        dist = [newMean,newvar**.5]
        varOut = np.max(newvar)
        YvarOut = np.var(Y)
        varOut = np.min([varOut,10.*YvarOut])
        counter = counter + 1
        Xoldn = copy.copy(Xn)
	Yold = Y[indexSorted[0:10]]

    if counter>=nMax:
        message='Maximum number of iterations reached'
        success = False
    else:
        message='Optimal value is possible'
        success=True
    return (newMean,newvar,success,message)
    


if __name__ == '__main__':
    myfun = lambda X:(X[0] + .5)**6+(X[1]-.5)**2
    
    bounds = [[-100,100],[-100,100]]
    (Xopt,varOut,success,message) =optimMin(bounds=bounds,fun=myfun)
    #print Xopt,varOut
    #print message
