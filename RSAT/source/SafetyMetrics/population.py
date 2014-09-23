import simpleArcReader
import numpy as np
import simpleArcReader as SAR

# This functions are used to read an 'asc' type file. In this case the 'asc' file has the population density in [Persons/km^2]

class populationClass:
    def __init__(self,PopFile=None,AreaFile=None):
        
        keyPop,xllcorner0,yllcorner0,cellsize0,nrows0,ncols0,xMax0,yMax0= SAR.getRawValsPopDensity(PopFile)
        keyArea,xllcorner,yllcorner,cellsize,nrows,ncols,xMax,yMax= SAR.getRawValsPopDensity(AreaFile)

        if ((xllcorner0-xllcorner)**2+(yllcorner-yllcorner0)**2 +(cellsize-cellsize0)+(nrows-nrows0)**2>0.00000001):
            print 'Population Count and Population Area parameters do not match'
            error(1)

        self.keyPop = keyPop
        self.keyArea = keyArea
        self.xllcorner = xllcorner
        self.yllcorner = yllcorner
        self.cellsize = cellsize
        self.xMax = xMax
        self.yMax = yMax
        self.nrows = nrows
        self.ncols = ncols
        self.keyDen = None


class polygon:
    def __init__(self,xy=None):
        #xy is a list with point pairs that define the polygon e.g xy = [[4,5],[1,2],[3,3]]
        self.xy = xy
    
    def square(self,length=1.0,center=None):
        if center==None:
            center = self.center
        else:
            self.center = center
        
        xMax = center[0] + 0.5*length
        xMin = center[0] - 0.5*length
                
        yMax = center[1] + 0.5*length
        yMin = center[1] - 0.5*length
        x = [xMin,xMin,xMax,xMax]
        y = [yMin,yMax,yMax,yMin]
        self.xy = [x,y]

    def max(self,xy):
        xy = np.array(self.xy)
        self.xmax = np.max(xy[:,0])
        self.xmin = np.min(xy[:,0])
        self.ymax = np.max(xy[:,1])
        self.ymin = np.min(xy[:,1])


    
        


# determine if a point is inside a given polygon or not
# Polygon is a list of (x,y) pairs.
# code obtained from http://www.ariel.com.au/a/python-point-int-poly.html
def point_inside_polygon(x,y,poly):
    
    n = len(poly)
    inside =False
    
    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y
    
    return inside



class populationDensityClass:

    def __init__(self,PopDenFile=None):
        
        keyDen,xllcorner,yllcorner,cellsize,nrows,ncols,xMax,yMax= SAR.getRawValsPopDensity(PopDenFile)
        self.keyDen = keyDen
        self.xllcorner = xllcorner
        self.yllcorner = yllcorner
        self.cellsize = cellsize
        self.xMax = xMax
        self.yMax = yMax
        self.nrows = nrows
        self.ncols = ncols


def Density(X,Y,xlen,ylen,filename):
# Takes in the X, Y which are equivalent to latitude and longitude

  	#filename = 'usads05ag.asc'
	totalLength = int(xlen*ylen)
	
	xvec = np.reshape(X,totalLength) # making X a vector
	yvec = np.reshape(Y,totalLength) # making Y a vector

	popDensity = simpleArcReader.arcReader(xvec,yvec,totalLength,filename) # this result is a vector
	popDensity = np.asarray(popDensity)
	popDensity = popDensity.reshape(ylen,xlen)

	return popDensity


def DensityfromRaw(X,Y,xlen,ylen,key,xllcorner,yllcorner,cellsize,nrows,xMax):
	totalLength = int(xlen*ylen)
	
	xvec = np.reshape(X,totalLength) # making X a vector
	yvec = np.reshape(Y,totalLength) # making Y a vector

	popDensity = simpleArcReader.getVal(key,xllcorner,yllcorner,cellsize,nrows,xMax,xvec,yvec,totalLength)
	popDensity = np.asarray(popDensity)
	popDensity = popDensity.reshape(ylen,xlen)
	return popDensity

def DensityRaw(filename):
# simply making a call to simpleArcReader, check simpleArcReader for information
# This function is here just to keep all the population calls together, this function could be avoided all together
	key,xllcorner,yllcorner,cellsize,nrows,ncols,xMax,yMax = simpleArcReader.getRawVals(filename)
	return (key,xllcorner,yllcorner,cellsize,nrows,xMax)




def gridSize(filename):
# this function return the cellsize for an asc file
	try:
		inputFile = open(filename,'r')
	except:
		print 'Cannot open file==> ERROR in population.gridSize'

	for line in inputFile:
		key = line.split()
		if key[0] == 'cellsize':
			cellsize = float(key[1])
			break
	return cellsize
def cellArea(filename):
	import math
# Calculates the area of each grid in meter square by assuming a perfect spherical earth
# this section could improve later on to include Earth's oblateness
	cellsize = gridSize(filename)
	cellsize = cellsize*(math.pi)/180.0 # converting cellsize to radians
	earthRadius = 6378.145*1000.0 # [m]
	deltaL = earthRadius*cellsize
	areaCell = deltaL*deltaL
	return areaCell # in meters^2

def deltaArea(delta):
	import math
# calculates the area of a selected delta[deg] by assuming a spherical Earth
	earthRadius = 6378.145*1000.0 # [m]
	delta = delta*(math.pi)/180.0
	deltaL = delta*earthRadius
	area = deltaL*deltaL
	return area

def interestArea(meanList,covarList):
# Calculates the area of interest by taking in all the results from all groups
# the result is a big rectangle that contains the areas of interest
	
	nGroups = len(meanList)

	XU = -999.0
	XL =  999.0
	YU = -999.0
	YL = 999.0
	nSigma = 5.0
	for index in range(0,nGroups):
		mu = meanList[index]
		covar = covarList[index]
		if mu!='NAN':
			sigma1 = np.sqrt(covar[0,0])
			sigma2 = np.sqrt(covar[1,1])
			xmin = mu[0]-nSigma*sigma1
			xmax = mu[0]+nSigma*sigma1
			ymin = mu[1]-nSigma*sigma2
			ymax = mu[1]+nSigma*sigma2
			if XL >xmin:
				XL = xmin
			if XU < xmax:
				XU = xmax
			if YL > ymin:
				YL = ymin
			if YU < ymax :
				YU = ymax

	return (XL,XU,YL,YU)

def interestAreaGauss(mu,covar):
	XU = -999.0
	XL =  999.0
	YU = -999.0
	YL = 999.0
	nSigma = 5.0
	sigma1 = np.sqrt(covar[0,0])
	sigma2 = np.sqrt(covar[1,1])
	xmin = mu[0]-nSigma*sigma1
	xmax = mu[0]+nSigma*sigma1
	ymin = mu[1]-nSigma*sigma2
	ymax = mu[1]+nSigma*sigma2
	if XL >xmin:
		XL = xmin
	if XU < xmax:
		XU = xmax
	if YL > ymin:
		YL = ymin
	if YU < ymax :
		YU = ymax

	return (XL,XU,YL,YU)
	

def interestAreaKernel(y):
# y is just a nx2 numpy array with lon lat
	noffset = .25 #1.5 bu default...need to work on this parameter
	xb = y[:,0]
	yb = y[:,1]
	XL = min(xb)
	XU = max(xb)
	YL = min(yb)
	YU = max(yb)
	
	dx = XU-XL
	dy = YU-YL
	XL = XL-noffset*dx
	XU = XU+noffset*dx
	YL = YL-noffset*dy
	YU = YU+noffset*dy
	return (XL,XU,YL,YU)
		
	
	
