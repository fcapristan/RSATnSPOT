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

np.random.seed(3)
random.seed(3)


def cost(distribution,U,Scaling,Model,Xdesign,basevals):
    vals = -np.sum(Xdesign) #maximize
    
    distribution.std = [Xdesign[0],Xdesign[1],Xdesign[2],10.,10.,10.]
    E,std = surroFunc(distribution,U,Scaling,Model,Xdesign,basevals)
    g = np.array([E-73000.,std-41000.])
    fail = 0
    vals = vals
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


def surroFunc(distribution,U,Scaling,Model,Xdesign,basevals):
    #distribution.std = Xdesign
    distribution.std = [Xdesign[0],Xdesign[1],Xdesign[2],10.,10.,10.]
    X = active.sampleRandom(distribution,basevals)
    Xnormalized = active.normalizeX(distribution,X)
    XS = np.dot(U.T,Xnormalized) # new inputs for surrogate model
    
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

    E = np.mean(FI)
    std = np.std(FI)
        

    return E,std



distList = active.activeDistribution()

distList.name = ['uniform','uniform','uniform','uniform','normal','normal']
distList.mu = [10,15.0,-10.,15.,5,.2]
xVec = [5.,5.,10.]

distList.std = [5.,5.,10.,10.,10.,10.]
distList.UpdateMinMaxLimits()



nSamples = 1000
basevals = active.baseSample(distList,nSamples)
vars = active.sampleRandom(distList,basevals)


step = 0.0001*np.ones((len(distList.name)))
U,XS,Xnormalized,X = active.getSubspaceSimple(func,distList,step,nSamples,basevals,nsubspace=2)


print 'U is'
print U


fvals = np.zeros((nSamples))
for index in range(nSamples):
    fvals[index] = func(vars[:,index])

for index in range(len(vars)):
    plt.figure()
    plt.plot(vars[index],fvals,'x')




The_Func = func
ns = nSamples
FS = np.zeros((1,ns))

for index in range(ns):
    FS[0,index] = The_Func(vars[:,index])


activeDims,samples = np.shape(XS)

XB = np.zeros((activeDims,2))

for index in range(activeDims):
    XB[index,0] = -1.1
    XB[index,1] =  1.1

XS = XS.T
FS = FS.T
#DFS = DFS.T

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




nSamples = 1000
basevals = active.baseSample(distList,nSamples)


fun = lambda Xdesign:cost(distList,U,Scaling,Model,Xdesign,basevals)

print fun(xVec)

#bnds = ((0.001,10),(0.001,10))

#res = minimize(fun, (2, 0), method='SLSQP', bounds=bnds)
#print res
print 'SNOPT'


opt_prob = pyOpt.Optimization('Trial UQ',fun)
opt_prob.addVar('var1','c',lower=0.1,upper=7.,value=xVec[0])
opt_prob.addVar('var2','c',lower=0.1,upper=7.,value=xVec[1])
opt_prob.addVar('var3','c',lower=0.1,upper=100.,value=xVec[2])
#opt_prob.addVar('var4','c',lower=0.1,upper=1000.,value=xVec[3])
#opt_prob.addVar('var5','c',lower=0.1,upper=10.,value=xVec[1])
#opt_prob.addVar('var6','c',lower=0.1,upper=10.,value=xVec[1])



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

slsqp_prob = pyOpt.SLSQP()
exit2Val = slsqp_prob(opt_prob)
print exit2Val





'''
Xdesign = xVec


nSamples = 1000
basevals = active.baseSample(distList,nSamples)
X = active.sampleRandom(distList,basevals)

E,V,FI=surroFunc(distList,U,Scaling,Model,Xdesign,basevals)

fcheck = np.zeros((nSamples))
for index in range(nSamples):
    fcheck[index] = func(X[:,index])
print 'Mean',np.mean(fcheck)
print 'STD',np.std(fcheck)




print E,V


exit()
'''


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

# show the plot

# plot the training data
pnts = ax.plot(XS[:,0],XS[:,1],FS[:,0],'bx',mew=3,ms=10,label='Test Samples')




ax.set_xlabel('Reduced Coordinate 1')
ax.set_ylabel('Reduced Coordinate 2')
ax.set_zlabel('Casualty Measure')
# plot the truth response surface
#truh = ax.plot_surface(xi,yi,ft, cmap=cm.autumn, rstride=1,cstride=1, linewidth=0, antialiased=False)



plt.show(block=True)


exit()


'''



