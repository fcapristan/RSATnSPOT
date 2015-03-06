import matplotlib.pyplot as plt
import numpy as np

def plotCobWeb(X,Y,Ybounds):
   
    # assuming that X is a matrix with nSamples vector inputs
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
    xlabel = ['X']*(ndim+1)
    fig = plt.figure()
    for indim in range(ndim):
        Xnorm[indim,:] = (X[indim,:] - xMin[indim])/deltax[indim]
        xlabel[indim] = xlabel[indim]+str(indim)
        plt.plot([indim,indim],[0,1],'k',figure=fig)
    varOut = np.concatenate((Xnorm,[Ynorm]),axis=0)
    xlabel[ndim]  = 'Y'
    index = (Y<=Ybounds[1])&(Y>=Ybounds[0])

    #print Ybounds
    try:
        plt.plot(varOut[:,~index],'c',figure=fig)
    except:
        print 'Warning in CobWeb Plot: All samples are inside desired region'
    try:
        plt.plot(varOut[:,index],'r',figure=fig,linewidth=2)
    except:
        print 'Warning in CobWeb Plot: All samples outside desired region'
    plt.xticks(range(ndim+1),xlabel,figure=fig)
    return fig


if __name__ == "__main__":
    X = np.random.uniform(low=0.,high=1.,size=(2,300))
    Y = X[0,:]+ X[1,:]
    #Y = X[0,:]**2 + X[3,:]**2 + X[1,:]
         
    #Ybounds = [2.5,3.]
    Ybounds = [1.5,5.]
    plotCobWeb(X,Y,Ybounds)

