import numpy as np
import serialKernel as SK

# Driver for kernelEstimation, developed for testing

import matplotlib.pyplot as plt

# first generate random bivariate gaussian points
np.random.seed(1)
mean  =[10,-5]

cov = [[1,-.6],[-.6,1]]

x,y = np.random.multivariate_normal(mean,cov,5000).T
xnew = np.array(np.matrix(x).T)
ynew = np.array(np.matrix(y).T)

xvec = np.linspace(0,20,600)
yvec = np.linspace(-15,5,600)
dx = xvec[2]-xvec[1]
dy = yvec[2]-yvec[1]

xy = np.concatenate((xnew,ynew),1)
#print xy
#print np.std(ynew)
X,Y = np.meshgrid(xvec,yvec)




fhatFAST = SK.serialK(5000,xy,600,600,X,Y,2,1)
print fhatFAST

print np.sum(fhatFAST*dx*dy)

plt.figure(1)
plt.plot(x,y,'x'); plt.axis('equal')



plt.figure(2)
plt.contour(X,Y,fhatFAST,50)
plt.show()

