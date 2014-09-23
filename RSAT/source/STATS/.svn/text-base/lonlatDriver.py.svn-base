# driver for lonlatChecks
import lonlatChecks
import numpy as np


nsamples = 100
lon = np.linspace(-10,220,nsamples)
lat = np.linspace(33,70,nsamples)


lonlat = np.zeros((nsamples,2))

lonlat[:,0] = lon
lonlat[:,1] = lat

print lonlat

lonlat2 = lonlatChecks.boundlon(lonlat,nsamples)

lonlat3 = lonlatChecks.fixlon4pdf(lonlat2,nsamples)

print lonlat-lonlat3
