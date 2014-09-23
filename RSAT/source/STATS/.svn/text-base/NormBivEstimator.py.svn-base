import numpy as np
from scipy.optimize import fmin
from scipy.integrate import dblquad


# Takes Bivariate Normal estimator
# calculates the mean by looking for the centroid

def normalEstimator(yin):
    y = yin
    y = np.matrix(y)
    nSamples = np.size(y,axis=0)
    muCalc = np.mean(y,axis=0).T
    covar = 0
    for index in range(0,nSamples):
	Xi = y[index,0:2]
	Xi = Xi.T
	prod = (Xi-muCalc)
	covar = covar + prod*(prod.T)
#	print covar
    covar = (1.0/(nSamples-1.0))*covar

    mle_mean = muCalc
    mle_cov = covar
    return (mle_mean,mle_cov)


def MLEestimator(yin):
    y = yin
    y = np.matrix(y)
    nSamples = np.size(y,axis=0)
    muCalc = np.mean(y,axis=0).T
    covar = 0
    for index in range(0,nSamples):
        Xi = y[index,0:2]
        Xi = Xi.T
        prod = (Xi-muCalc)
        covar = covar + prod*(prod.T)
    #	print covar
    covar = (1.0/(nSamples-1.0))*covar
    
    mle_mean = muCalc
    mle_cov = covar
    return (mle_mean,mle_cov)

def normal_bivariate(mu,covar,X,Y):
#	print mu
	mu1 = mu[0,0]
	mu2 = mu[1,0]
#	print mu1
#	print mu2
	sigma1 = np.sqrt(covar[0,0])
	sigma2 = np.sqrt(covar[1,1])
	rho = covar[0,1]/(sigma1*sigma2)
	X1 = np.array(X)
	X2 = np.array(Y)

	z = ((X1-mu1)/sigma1)**2-2*rho*(X1-mu1)*(X2-mu2)/(sigma1*sigma2)+((X2-mu2)/sigma2)**2
	pdf = (np.exp(-z/(2*(1-rho**2))))/(np.pi*2*sigma1*sigma2*np.sqrt(1-rho**2))
	return pdf 
def surfIntegral(mu,covar,xin,yin,delta):
   # print xin,yin,delta
   # print mu,covar
    area = dblquad(lambda Y,X:normal_bivariate(mu,covar,X,Y),xin-delta,xin+delta, lambda X: yin-delta,lambda X:yin+delta)
    return area
'''
def myDoubleIntegral(mu,covar,xin,yin,delta):
    error = 1.0

    area = (mu,covar,xin,yin,delta)
    while error>.01
'''       

def mySimpleIntegral(mu,covar,xin,yin,delta):
    xL = xin - delta/4
    xU = xin + delta/4
    yU = yin + delta/4
    yL = yin - delta/4

    f1 = normal_bivariate(mu,covar,xU,yU)
    f2 = normal_bivariate(mu,covar,xL,yU)
    f3 = normal_bivariate(mu,covar,xL,yL)
    f4 = normal_bivariate(mu,covar,xU,yL)

    area = .25*(f1+f2+f3+f4)*delta**2
    return area
    
    


'''
plt.figure
plt.plot(y[:,0],y[:,1],'x')
#plt.axis('equal')
#plt.show()

delta = 0.0001
xx = np.arange(mle_mean[0]-5*sigma1,mle_mean[0]+5*sigma1,delta)
yy = np.arange(mle_mean[1]-5*sigma2,mle_mean[1]+5*sigma2,delta)
XX, YY =np.meshgrid(xx,yy)
ZZ = normal_bivariate(mle_mean,mle_cov,XX,YY)
#print ZZ
#ZZ = mlab.bivariate_normal(XX,YY,mle_  

CS = plt.contour(XX,YY,ZZ)
plt.axis('equal')
plt.show()
'''


