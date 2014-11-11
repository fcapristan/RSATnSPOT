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


def optimMin(bounds=None,dist=None,fun=None,epsilonVar=0.000001,nMax = 5000):
    #this works better if everything is scaled appropriately
    varOut = np.inf
    nsamples = 60

    counter = 1
    while (varOut > epsilonVar)and(counter<= nMax):
        X = sample(bounds=bounds,dist=dist,nsamples=nsamples) # samples X
        Xmat = np.array(X)
        Y = np.zeros((nsamples))
        for index in range(nsamples):
            Y[index] = fun(X[index])
        indexSorted = np.argsort(Y)
        indexSorted = indexSorted[0:10]
        Xn = Xmat[indexSorted]
        newMean = np.mean(Xn,axis=0)
        newvar = np.var(Xn,axis=0)
        dist = [newMean,newvar**.5]
        varOut = np.max(newvar)
        YvarOut = np.var(Y)
        varOut = np.min([varOut,YvarOut])
        counter = counter + 1
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
    print Xopt,varOut
    print message