import numpy as np

def areaOfInterest(mean,covar,delta,nsigma):
# delta = spacing in between points
# nsgima = amount to move from mean

    #minvals = np.min(lonlat,axis=0) # getting min values for lat and lon
    #maxvals = np.max(lonlat,axis=0) # geting max values for lat and lon
    
    #minlon = minvals[0]
    #minlat = minvals[1]
    #maxlon = maxvals[0]
    #maxlat = maxvals[1]
    
    sigma1 = np.sqrt(covar[0,0])
    sigma2 = np.sqrt(covar[1,1])
    xmax = mean[0]+nsigma*sigma1
    xmin = mean[0]-nsigma*sigma1
    ymin = mean[1]-nsigma*sigma2
    ymax = mean[1]+nsigma*sigma2
    
    if (xmax-xmin >250):
        xmax = mean[0]+ 125
        xmin = mean[0] -125
    
    
            #print xmin,xmax
    '''
    
    if xmin<-250.: # 200 for extreme cases. Just trying to give some continuity to Earth. Population densities should adjust to handle angles out of this bounds
        xmin = -250
    if xmax > 250.: 
        xmax = 250 
    
    if ymin<-90. : # no need for continuity. Population density are neglected near the poles, so values in these areas will be neglected
        ymin = -90.
    if ymax>90.:
        ymax = 90.
    print sigma1,xmin
    '''
    
    xx = np.arange(xmin,xmax,delta)
    yy = np.arange(ymin,ymax,delta)
    xlen = len(xx)
    ylen = len(yy)
    XX, YY =np.meshgrid(xx,yy)

    
    return XX,YY,xlen,ylen
    #return xx,yy,xlen,ylen

def areaOfInterestKDE(lonlat,delta=None,nMesh = 50):
    x = lonlat[:,0]
    y = lonlat[:,1]
    
    xmin = min(x) 
    xmax = max(x)
    ymin = min(y)
    ymax = max(y)
    dx = np.min([xmax-xmin,105])
    dy = np.min([ymax-ymin,105])
    

    xmax = xmax+.3*dx
    xmin = xmin-.3*dx
    ymax = ymax+.3*dy
    ymin = ymin-.3*dy
    

    
    if xmin<-200.: # 200 for extreme cases. Just trying to give some continuity to Earth. Population densities should adjust to handle angles out of this bounds
        xmin = -200
    if xmax > 200.: 
        xmax = 200
    if ymin<-90. : # no need for continuity. Population density are neglected near the poles, so values in these areas will be neglected
        ymin = -90.
    if ymax>90.:
        ymax = 90.
    
 
    if delta!=None:
        xx = np.arange(xmin,xmax,delta)
        yy = np.arange(ymin,ymax,delta)
        areaIntegral = delta**2.0
    elif delta==None:
        xx = np.linspace(xmin,xmax,nMesh)
        yy = np.linspace(ymin,ymax,nMesh)
        dx = xx[1]-xx[0]
        dy = yy[1]-yy[0]
        areaIntegral = dx*dy
    xlen = len(xx)
    ylen = len(yy)
    XX, YY =np.meshgrid(xx,yy)

    #print xlen,ylen
    #print np.shape(XX)
    #exit()
    
    return XX,YY,xlen,ylen,areaIntegral 
    #return xx,yy,xlen,ylen

    
