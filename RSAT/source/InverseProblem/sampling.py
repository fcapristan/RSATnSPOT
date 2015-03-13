import numpy as np
#import matplotlib.pyplot as plt


#simple LHS for uniform distributions [0,1)


def getLHSuniform(dim,samples,seed=None):
    if seed!=None:
        np.random.seed(seed)
    pij  = np.zeros((dim,samples))
    for index in range(dim):
        pij[index,:] =  np.random.permutation(samples)
    Uij = np.random.uniform(low=0,high=1,size=(dim,samples))
    Xij = (pij + Uij)/float(samples)

    return Xij

'''
nSamples = 1000
dim = 2

X = getLHSuniform(dim,nSamples)

print X

fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(X[0,:],X[1,:],'o')
ax.xaxis.set_ticks(np.linspace(0,1,nSamples+1))
ax.yaxis.set_ticks(np.linspace(0,1,nSamples + 1))

ax.grid(True)


X2 = np.random.uniform(low=0,high=1,size=(dim,nSamples))

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
plt.plot(X2[0,:],X2[1,:],'o')
ax2.xaxis.set_ticks(np.linspace(0,1,nSamples+1))
ax2.yaxis.set_ticks(np.linspace(0,1,nSamples + 1))

ax2.grid(True)
plt.show()



'''
