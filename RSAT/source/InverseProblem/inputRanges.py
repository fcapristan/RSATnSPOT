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

def binInputs2(X,Y,bin_limits):

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
    bin_limits[-1]= 1.001
    for index in range(len(bin_limits)-1):
        ides = (comboArr >= bin_limits[index] ) & (comboArr <  bin_limits[index + 1 ] )
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



def printBinInputs(treeList):
    
    
    
    #for index in range
    '''
    if len(treeList)<=1:
        print treeList ,'-' 
    else :
        print len(treeList[0]), treeList[0][0],printBinInputs(treeList[0][1]),'r-' 
        printBinInputs(treeList[1:len(treeList)])

    '''

    return None
'''

x0 = np.array([.67,.12,.67,.97,.567,1.,.3])
x1 = np.array([.1,.1,.2,.6,.8,.4,.1])
x2 = np.array([.2,.3,.22,.467,.6,.32,.4])
x3 = np.array([1.,1.,1.,1.,1.,1.,1.])
x4 = np.array([1.,.9,1.,1.,1.,1.,1.])

yout = np.array([.5,.4,.1,.7,.9])
X = np.concatenate(([x0],[x1],[x2],[x3],[x4])).T

#X = np.array([[0.1,0.1],[0.3,0.4]])

bin_limits = np.linspace(0.0,1.0,4)
print X

delta = bin_limits[1]-bin_limits[0]
inputsRange,outputsRange = binInputs2(X,yout,bin_limits)

import cobWebPlotting as CW

fig = CW.plotCobWeb(X,yout,[.7,1])


fig2 = CW.plotCobWeb(inputsRange*delta + .5*delta ,outputsRange*delta + .5*delta,[.7,1])

plt.show()
'''

