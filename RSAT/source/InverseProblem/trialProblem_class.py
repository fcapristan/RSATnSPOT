import numpy as np
import activeSubspaces as active
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import pyOpt

from matplotlib import cm


import VyPy
from VyPy.regression import gpr
import random
from warnings import warn, simplefilter
simplefilter('error')

np.random.seed(2)
random.seed(2)






def setupSNOPT(distributionClass,subspaceClass,X0,Optfun):
    print 'SNOPT'


    opt_prob = pyOpt.Optimization('Trial UQ',Optfun)
    curvar = 1
    upperLimit = distributionClass.stdmax
    for index in range(distributionClass.dim):
        if subspaceClass.important_index[index]==True:
            opt_prob.addVar('var'+str(curvar),'c',lower=0.001,upper=upperLimit[index],value=X0[curvar-1])
            curvar = curvar + 1


    opt_prob.addObj('Objective Mass')

    opt_prob.addCon('g1','i')
    opt_prob.addCon('g2','i')





    print opt_prob
    snopt = pyOpt.pySNOPT.SNOPT()


    snopt.setOption('Major feasibility tolerance',value=5e-6)
    snopt.setOption('Major optimality tolerance',value=1e-5)
    snopt.setOption('Minor feasibility tolerance',value=5e-6)
    snopt.setOption('Major iterations limit',500)


    exitVal = snopt(opt_prob,sens_type= 'FD')
    print exitVal






def cost(X,distributionClass,subspaceClass,Scaling,Model,samplesClass):
    
    distributionClass.std[subspaceClass.important_index] = X
    
    samplesClass.updateValues(distributionClass=distributionClass)#updating samples with pdated distribution
    FI = surroFunc(distributionClass,subspaceClass,Scaling,Model,samplesClass)
    
    E = np.mean(FI)
    std = np.std(FI)
    g = np.array([E-83000.,std-20000.])
    fail = 0
    vals = -np.sum(X)
    '''
    distribution.std = xVec
    vars = active.sampleRandom(distribution,basevals)
    
    
    fvals = func(vars)
    vals = -np.std(fvals)
    print 'v',vals
    g=np.mean(fvals)
    print 'g',g
    g  = g-21
    fail = 0
    '''
    print vals,g,E,std
    return vals,g,fail


def func(x):
    myval = 10.0*np.power(x[0],2.8) + 10.*np.power(x[1],3.0) + x[0]*x[0]*x[2]**2 + x[0]*x[1] + x[0]+x[1]+x[2]+x[3]+x[4]*x[5]
    
    return  myval


def surroFunc(distributionClass,subspaceClass,Scaling,Model,samplesClass):
    activeSubspaceVar = np.zeros((distributionClass.dim,samplesClass.nSamples))

    lowLim = distributionClass.minVal
    upLim = distributionClass.maxVal
    delta = upLim-lowLim
    for index in range(distributionClass.dim):
        activeSubspaceVar[index,:] = 2.0*((samplesClass.values[index,:] - lowLim[index])/delta[index]) - 1.0
    
    if activeSubspaceVar.max()>1.0 or activeSubspaceVar.min()<-1.0:
        print 'Error: New Sample Value exceeds limits [-1,1]' 
        exit()
    
    XS = np.dot((subspaceClass.U).T,activeSubspaceVar)
      
    XI = XS.T
    
    XI_scl = XI / Scaling.XI

    # ---------------------------------------------------------    
    # evaluate the model with the scaled features
    The_Data_Scl = Model.predict(XI_scl)
    # since we constructed a model with scaled training data,
    # we must provide scaled feature locations to predict    
    # ---------------------------------------------------------    

    # un-scale the estimated targets
    The_Data = Scaling.unset_scaling(The_Data_Scl)

    # pull the estimated target function values
    FI = The_Data.YI

    #E = np.mean(FI)
    #std = np.std(FI)
    '''
    Fcheck = np.zeros((samplesClass.nSamples))
    for index in range(samplesClass.nSamples):
        Fcheck[index] = func(samplesClass.values[:,index])
    print 'Surr',E,std
    print 'Actual',np.mean(Fcheck),np.std(Fcheck)
    '''        

    return FI




name = ['uniform','uniform','uniform','uniform','normal','normal']

mu = [10,15.0,-10.,15.,5,.2]
std = [3.,5.,10.,10.,10.,10.]
nSamples = 1000
#stdmax = [3.,5.,10.,10.,10.,10.]
stdmax = [5.,7.,30.,20.,20.,20.]
inputsDistribution = active.distribution(name=name,mu=mu,std=std,stdmax=stdmax) #setting distribution class


subspace = active.activeSubspace()
subspace.getSubspace(fun = func,distributionClass=inputsDistribution,nSamples=nSamples,nsubspace=2) # getting active subspace

print inputsDistribution.minVal,inputsDistribution.maxVal
print subspace.values.max(),subspace.values.min()


samples = active.samples()
samples.generateNew(distributionClass=inputsDistribution,nSamples=nSamples,uniform=0)#setting samples and base for current distribution

print inputsDistribution.std, subspace.important_index
X0 = inputsDistribution.std[subspace.important_index] # getting initial guess for important variables only
#lowerLimit = inputsDistribution.lowerLimit[important_index]
#upperLimit = inputsDistribution.upperLimit[important_index]


#generating surrogate for the main function

FS = np.zeros((1,nSamples))
XS = (subspace.subspaceBase).T


for index in range(nSamples):
    FS[0,index] = func(subspace.values[:,index])

XB = np.zeros((subspace.Nsubspace,2)) 
for index in range(subspace.Nsubspace):
    XB[:,0] = XS.min()# FIX THIS AND FINISH
    XB[:,1] = XS.max()
FS = FS.T


Train = gpr.training.Training(XB,XS,FS)
# find scaling factors for normalizing data
Scaling = gpr.scaling.Linear(Train)

# scale the training data
Train_Scl = Scaling.set_scaling(Train)
# In this case the scaling is performed by normalizing on 
# the experimental range ( max()-min() ) for input feature
# samples XS (by dimension) and output target samples FS.

# choose a kernel 
Kernel = gpr.kernel.Gaussian(Train_Scl,probNze=-2.0)

# choose an inference model (only one for now)
Infer  = gpr.inference.Gaussian(Kernel)

# choose a learning model (only one for now)
Learn  = gpr.learning.Likelihood(Infer)

# start a gpr modling object
Model  = gpr.modeling.Regression(Learn)
# This object holds all the model's assumptions and 
# data, in order to expose a simple interface.

# learn on the training data
Model.learn()


FI =surroFunc(inputsDistribution,subspace,Scaling,Model,samples)

Optfun = lambda X:cost(X,inputsDistribution,subspace,Scaling,Model,samples)


setupSNOPT(inputsDistribution,subspace,X0,Optfun)

'''

# ---------------------------------------------------------
#  Evaluate the Model
# ---------------------------------------------------------    
# Run a grid sample of the model for plotting

# grid number of points per dimension
nx = 70 

# generate a grid, in the original design space bounds
xi = np.linspace(XB[0,0], XB[0,1], nx) 
yi = np.linspace(XB[1,0], XB[1,1], nx)
xi,yi = np.meshgrid(xi,yi)

# zip data into a design matrix
XI = np.zeros([nx*nx,2])
it = 0
for ix in range(nx):
    for iy in range(nx):
        XI[it,:] = [xi[ix,iy],yi[ix,iy]]
        it = it+1

# scale feature locations to match the model
XI_scl = XI / Scaling.XI

# ---------------------------------------------------------    
# evaluate the model with the scaled features
The_Data_Scl = Model.predict(XI_scl)
# since we constructed a model with scaled training data,
# we must provide scaled feature locations to predict    
# ---------------------------------------------------------    

# un-scale the estimated targets
The_Data = Scaling.unset_scaling(The_Data_Scl)

# pull the estimated target function values
FI = The_Data.YI


# ---------------------------------------------------------
#  Evaluate the Truth Function
# ---------------------------------------------------------      
# Run the same grid on the truth function




# ---------------------------------------------------------
#  Plotting
# ---------------------------------------------------------
# plot the estimated and truth surface, evaluate rms error

print 'Plot Response Surface ...'

# unzip data for plotting a gridded surface
fi = xi*0. # estimated targets
ft = xi*0. # truth targets
it = 0
for ix in range(nx):
    for iy in range(nx):
        fi[ix,iy] = FI[it]
        it = it+1

# start the plot
fig = plt.figure()
ax = fig.gca(projection='3d')

# plot the estimated response surface
surf = ax.plot_surface(xi,yi,fi, cmap=cm.jet, rstride=1,cstride=1, linewidth=0, antialiased=False)
pnts = ax.plot(XS[:,0],XS[:,1],FS[:,0],'bx',mew=3,ms=10,label='Test Samples')
print 'min and max',XS.min(),XS.max()
plt.show(block=True)


'''










