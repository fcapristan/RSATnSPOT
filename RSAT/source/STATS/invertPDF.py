# Bisection-based inversion of a CDF
import numpy as np

'''
def invert(xaxis,P):
    #P is an array length N. P = Prob(X>=k) cum prob
    N = len(P)
    indexVec = range(N)
    u = np.random.uniform()
    if u==0.0:
        validIndex = P>0.0
        return xaxis[validIndex][0]
    L = 0
    R = N
    while L<(R-1):
        k = int(np.floor((L+R)/2.0))
        if u > P[k]:
            L=k
        else:
            R=k
    return xaxis[R],R
# add min
'''
import matplotlib.pyplot as plt
from scipy import stats
xk = [2.98787,8.28055,17.842,22.9655,28.172,33.1236,37.9975,43.2043,48.0719,53.0297,58.1514,63.1024,67.9776,73.1773,78.1287,82.9201,88.2168,92.9825,98.0902,103.123,107.988,112.939,118.232,124]
pk = (0.00056085,
      0.000988932,
      0.002432342,
      0.009169074,
      0.005251101,
      0.006727511,
      0.041845671,
      0.037016317,
      0.045145089,
      0.074129003,
      0.070848365,
      0.068966885,
      0.113473836,
      0.076933808,
      0.07520832,
      0.130627497,
      0.133642938,
      0.080685328,
      0.020303041,
      0.002633194,
      0.00133886,
      0.000850314,
      0.001221724,
      0)
custm = stats.rv_discrete(name='custm', values=(xk, pk))
print custm
print dir(custm)
h = plt.plot(xk, custm.pmf(xk))
plt.show()