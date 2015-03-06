# routine to bin the inputs to find input combinations that lead to a certain region
import numpy as np

def binFunction(inputs,outputs,nbins = 5):
    assert ( (np.max(inputs)<=1.0) and (np.min(inputs)>=0.))
    assert (np.max(outputs)<= 1.0 and np.min(outputs)>=0.)
    
    bin_limits = np.linspace(0,1,nbins+1) # setting size for bins
    dim, Nsamples = np.shape(outputs)
    bins = [bin_limits]*dim

    bin_out = [bin_limits]*1
    
    # looping over each limit
    
    for index in range(nbins + 1):
        i0 = index
        ip1 = index + 1
        if index == nbins:
            ides = (outputs <= bin_limits[i0]) & (outputs >= bin_limits[ip1])
        else:
            ides = (outputs <= bin_limits[i0]) & (outputs > bin_limits[ip1])
        X_ides = inputs[:,ides,None] # inputs in the desired range

'''

def binInputs(inputs,bin_limits):
    dim,Nsamples = np.shape(inputs)
    # fixing bin_limits to include bounds at 1
    #assert ( bin_limits[-1] ==1.)
    bin_limits[-1] = 1.00001
    ret = []#[None]*(len(bin_limits)-1)
    
    
    for index in range(len(bin_limits)-1):
        ides = (inputs[0,:] >= bin_limits[index] ) & (inputs[0,:] <  bin_limits[index + 1 ] )
        
        if np.sum(ides) != 0.0  :
            #print 'ii',index,ides
            ndim,nvar = np.shape(inputs[1:,ides])
            if ndim>0:
                print inputs[1:,ides]
                print 'b',binInputs(inputs[1:,ides],bin_limits)

                ret.append([index,binInputs(inputs[1:,ides],bin_limits)])
            else:
                print 'ii', index
                ret.append(index)
    return ret
'''
def binInputs(X,Y,bin_limits):

    ndim,nSamples = np.shape(X)
    assert(len(Y)==nSamples)
    xMax = np.amax(X,axis=1)
    xMin = np.amin(X,axis=1)
    deltax = xMax - xMin
    yMax = np.amax(Y)
    yMin = np.amin(Y)
    #print yMin,yMax
    Xnorm = np.zeros((ndim,nSamples))
    Ynorm = (Y - yMin)/(yMax-yMin)
 
    for indim in range(ndim):
        Xnorm[indim,:] = (X[indim,:] - xMin[indim])/deltax[indim]


    
    comboArr = np.concatenate((Xnorm,[Ynorm]),0)
    refInput = np.zeros((ndim+1,nSamples))
    bin_limits_Local = 1.0*bin_limits
    bin_limits_Local[-1] = 1.001
    for index in range(len(bin_limits)-1):

        ides = (comboArr >= bin_limits_Local[index] ) & (comboArr <  bin_limits_Local[index + 1 ] )
        refInput[ides] = index
    refList = (refInput.T).tolist()
    bins = np.array(purge_dublicates(refList))
            
    inputsRange = bins[:,:-1].T
    outputsRange = bins[:,-1]
    #print outputsRange
    # Now eliminate duplicates
    return inputsRange,outputsRange,[xMin,xMax],[yMin,yMax]#np.array(bins)#(refInput.T).tolist()# set(list(refInput))


def purge_dublicates(X):
    unique_X = []
    for i, row in enumerate(X):
        if row not in X[i + 1:]:
            unique_X.append(row)
    return unique_X


def index2sample(X,Y,yLimits,NnewSamples = 50,bin_limits=None):
    assert (yLimits[0]<yLimits[1])   

    # binning inputs and outputs
    if bin_limits==None:
        bin_limits = np.linspace(0,1.0,10)

    indexVals = np.arange(0,len(bin_limits)-1)

    inputsIndexRange,outputsIndexRange,[xMin,xMax],[yMin,yMax] = binInputs(X,Y,bin_limits)

    yLimitsScaled = (yLimits - yMin) / (yMax - yMin) # scaling to [0,1]
    if yLimitsScaled[1]==1.:
        yLimitsScaled[1] = 1.0001
    #desired_index_val = indexVals[ (yLimitsScaled[0] <= bin_limits) & (yLimitsScaled[1] > bin_limits)]

    #for index in range(len(bin_limits)-1):
    #    if  bin_limits[index] >= yLimitScaled[0] and bin_limits[index +1] <= yLimitScaled[1]


    desired_index_val = indexVals[ (yLimitsScaled[0] <= bin_limits[1:] ) & (yLimitsScaled[1] >= bin_limits[:-1])]
    indexTotal = np.arange(len(outputsIndexRange))
    indexFinal = []

    for index in range(len(desired_index_val)):

        indexMain = outputsIndexRange == desired_index_val[index]

        indexFinal  = np.concatenate((indexFinal,indexTotal[indexMain]),1)

    return indexFinal.astype(int),inputsIndexRange,outputsIndexRange,[xMin,xMax],[yMin,yMax]

def reesample(X,Y,yLimits,NnewSamples = 50,bin_limits=None):
    # binning inputs and outputs
    if bin_limits==None:
        bin_limits = np.linspace(0,1.0,40)

    index_desired,inputsIndexRange,outputsIndexRange,[xMin,xMax],[yMin,yMax] = index2sample(X,Y,yLimits,NnewSamples = NnewSamples,bin_limits=bin_limits)

    dim = len(xMin)
    Xnew = np.random.uniform(0,1,size=(dim,NnewSamples))

    # creating one sample at a time, from the desired range.
    binIndex = 0
    Xout = np.zeros((dim,NnewSamples))
    
    for indexSample in range(NnewSamples):

        if binIndex >= len(index_desired):
            binIndex = 0

        Xout[:,indexSample] = sampleOneAtaTime(inputsIndexRange[:,index_desired[binIndex]],xMin,xMax,bin_limits)
        binIndex = binIndex + 1

    return Xout

def sampleOneAtaTime(index_desired,xMin,xMax,bin_limits):
    assert (len(index_desired)==len(xMin))
    bin_scaled =  bin_limits
    xout = np.random.uniform(0,1,size=len(xMin))
    for index in range(len(xMin)):

        bin_bounds = [bin_limits[index_desired[index]],bin_limits[index_desired[index] + 1]]
        upperVal = bin_bounds[1] * (xMax[index]-xMin[index]) + xMin[index]
        lowerBound = bin_bounds[0] * (xMax[index]-xMin[index]) + xMin[index]
        xout[index] =  xout[index] * (upperVal - lowerBound) + lowerBound
 
    return xout


def resampleAcceptanceRejection(U,X,Y,yLimits,NnewSamples=50,bin_limits=None):    
    # binning inputs and outputs
    if bin_limits==None:
        bin_limits = np.linspace(0,1.0,4)


    index_desired,inputsIndexRange,outputsIndexRange,[xMin,xMax],[yMin,yMax] = index2sample(X,Y,yLimits,NnewSamples = NnewSamples,bin_limits=bin_limits)
    orig_dim,sam = np.shape(U)
    dim = len(xMin)
    Xnew = np.random.uniform(0,1,size=(dim,NnewSamples))

    # creating one sample at a time, from the desired range.
    binIndex = 0
    Xout = np.zeros((orig_dim,NnewSamples))

    for index in range(len(xMin)):

        if binIndex >= len(index_desired):
            binIndex = 0
        Xout[:,index] = sampleOneAcceptanceRejection(U,inputsIndexRange[:,index_desired[binIndex]],xMin,xMax,bin_limits)
        binIndex = binIndex +  1

    return Xout

def sampleOneAcceptanceRejection(U,index_desired,xMin,xMax,bin_limits):
  
    assert (len(index_desired)==len(xMin))
    bin_scaled =  bin_limits

    ori_dim,red_dim = np.shape(U)


    index_desired = np.array(index_desired).astype(int)
    bin_bounds = [bin_limits[index_desired],bin_limits[index_desired + 1]]
    lowerBounds = bin_bounds[0] * (xMax - xMin) + xMin
    upperBounds = bin_bounds[1] * (xMax - xMin) + xMin



    while True:
        
        Xc = np.random.uniform(-1.,1.,size=ori_dim)
        X_red = np.dot(U.T,Xc)
        indexCheck = (lowerBounds <= X_red) & (upperBounds >= X_red)
   
        
        if np.sum(indexCheck) == red_dim :
 
            break
    return Xc




def sampleFun(X):
    ndim,nSamples = np.shape(X)
    Y = np.zeros((nSamples))
    for index in range(nSamples):
        x0 = X[0,index]
        x1 = X[1,index]
        x2 = X[2,index]
        x3 = X[3,index]
        x4 = X[4,index]
        Y[index] = x0**2 + x1 + 0.1*x2 + x3**2 + x4**3
    return Y
        
'''
np.random.seed(1)
nSamples = 100
x0 = np.random.uniform(-50,1,size=(1,nSamples))
x1 = np.random.uniform(0,30,size=(1,nSamples))
x2 = np.random.uniform(0,25,size=(1,nSamples))
x3 = np.random.uniform(56,58,size=(1,nSamples))
x4 = np.random.uniform(0,1,size=(1,nSamples))


X = np.concatenate((x0,x1,x2,x3,x4),0)
#print X
#print np.shape(X)

yout = sampleFun(X)


import matplotlib.pyplot as plt
import cobWebPlotting as CW




yrange = np.array([.7,1])*(max(yout) - min(yout)) + min(yout)
fig = CW.plotCobWeb(X,yout,yrange)


#fig2 = CW.plotCobWeb(X ,yout,[.7,1])





#X = np.array([[0.1,0.1],[0.3,0.4]])

bin_limits = np.linspace(0.0,1.0,4)


delta = bin_limits[1]-bin_limits[0]
#inputsRange,outputsRange = binInputs(X,yout,bin_limits)

#indexFinal = index2sample(X,yout,[.3,.8])
Xout = reesample(X,yout,yLimits=yrange,NnewSamples = 200)
#print Xout

nX = np.concatenate((X,Xout),1)
nY = np.concatenate((yout,sampleFun(Xout)),1)

fig2 = CW.plotCobWeb(nX,nY,yrange)
plt.show()
'''
'''
import cobWebPlotting as CW

fig = CW.plotCobWeb(X,yout,[.7,1])


fig2 = CW.plotCobWeb(inputsRange*delta + .5*delta ,outputsRange*delta + .5*delta,[.7,1])

plt.show()
'''

