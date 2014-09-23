import blastPredict as BP
import numpy as np
import matplotlib.pyplot as plt


x = np.linspace(.1,10,1000)

y = BP.VCEcurves(x,3)

plt.loglog(x,y)
plt.show()


