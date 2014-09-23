import population
import NormBivEstimator as NBE
import numpy as np
import kernelEstimation as KE

def expectationKernel(y,filename,fragArea,Npieces,kSpeed):
    nGroups = len(Npieces)
    gridArea = population.cellArea(filename)
    cfraction = .1
    cKernel = .1
    deltaGrid = population.gridSize(filename)
    delta = cfraction*deltaGrid
    area = delta*delta
    rp = .3048
    key,xllcorner,yllcorner,cellsize,nrows,xMax = population.DensityRaw(filename)
    E = np.zeros(nGroups)
    XPlist = []
    YPlist = []
    fkernelList = []
    for index in range(0,len(y)):
        if Npieces[index]>.1:
            yin = y[index]
            XL,XU,YL,YU = population.interestAreaKernel(yin)
       #     print XL,XU,YL,YU
            xpop = np.arange(XL,XU,delta)
            ypop = np.arange(YL,YU,delta)
            xlen = len(xpop)
            ylen = len(ypop)
            XP,YP = np.meshgrid(xpop,ypop)
            popDensity = population.DensityfromRaw(XP,YP,xlen,ylen,key,xllcorner,yllcorner,cellsize,nrows,xMax)
      #      print popDensity
            if kSpeed=='fast':
                print 'Using Fast KDE'
                fkernel = KE.fhatKernelFAST(xpop,ypop,yin,deltaGrid*cKernel,deltaGrid*cKernel)
            else:
                print 'Using Regular KDE'
                fkernel = KE.fhatKernel(xpop,ypop,yin)
            Ec=calculateE(fkernel,area,fragArea[index],rp,popDensity)

            E[index]=Ec
            XPlist.append(XP)
            YPlist.append(YP)
            fkernelList.append(fkernel)
        else:
            E[index]=0.0
    Etotal = E*Npieces

    Efinal = np.sum(Etotal)
    return XPlist,YPlist,fkernelList,Efinal
    
def fkernel(xy):
    XL,XU,YL,YU = population.interestAreaKernel(xy)
    #delta = .001*(XU-XL)
    delta = .01
    xpop = np.arange(XL,XU,delta)
    ypop = np.arange(YL,YU,delta)
    xlen = len(xpop)
    ylen = len(ypop)
    XP,YP = np.meshgrid(xpop,ypop)
   # f = KE.fhatKernelFAST(xpop,ypop,xy,delta,delta)
    f = KE.fhatKernel(xpop,ypop,xy)
    return XP,YP,f

def fgauss(xy):
    # xy is nx2 containing all lat lon locations
    mu,covar = NBE.MLEestimator(xy)
    print covar
    print covar[1,1]
    XL,XU,YL,YU = population.interestAreaGauss(mu,covar)# adding index zero because mu in MLEestimator is returned as a list
    delta = .001*(XU-XL)
    xpop = np.arange(XL,XU,delta)
    ypop = np.arange(YL,YU,delta)
    xlen = len(xpop)
    ylen = len(ypop)
    XP,YP = np.meshgrid(xpop,ypop)
    f = NBE.normal_bivariate(mu,covar,XP,YP)
    return XP,YP,f
    


def expectation(yList,filename,fragArea,Npieces):
    nGroups = len(Npieces)
    gridArea = population.cellArea(filename)
    cfraction = .1
    delta = cfraction*population.gridSize(filename)
    muList=[]
    covarList=[]
    for index in range(0,len(yList)):
        if Npieces[index]>.1:
            mu,covar = NBE.normalEstimator(yList[index])
        else:
            mu = 'NAN'
            covar = 'NAN'
        muList.append(mu)
        covarList.append(covar)
    XL,XU,YL,YU = population.interestArea(muList,covarList)
    print XL,XU,YL,YU
    xpop = np.arange(XL,XU,delta)
    ypop = np.arange(YL,YU,delta)
    xlen = len(xpop)
    ylen = len(ypop)
    XP,YP = np.meshgrid(xpop,ypop)
    print 'Reading Population Density from asc file'
    key,xllcorner,yllcorner,cellsize,nrows,xMax = population.DensityRaw(filename)
    popDensity = population.DensityfromRaw(XP,YP,xlen,ylen,key,xllcorner,yllcorner,cellsize,nrows,xMax)
    #print np.shape(XP)
   # print np.shape(YP)
   # print xlen
    print 'Done reading Population Density from asc file'
    area = delta*delta
    print delta,area
#    Ap = 1.0
    rp = .3048

    E = np.zeros(nGroups)
    Egrid=np.zeros([ylen,xlen])

    for index in range(0,nGroups):
        mu = muList[index]
        covar = covarList[index]

        if mu =='NAN':
            f = 0.0*XP # MUST CHANGE THIS
        else:
            f = NBE.normal_bivariate(mu,covar,XP,YP)
        Ec=calculateE(f,area,fragArea[index],rp,popDensity)

        '''
        Ec = 0
        I=10
        J=10
        delta=10.0
        P = NBE.surfIntegral(mu,covar,xpop[J],ypop[I],0.5*delta)
        print P
        '''
        '''
        Ec = 0
        Pt =0.0
        for J in range(0,xlen):
            for I in range(0,ylen):
                P = (f[I][J])*area
                #P,error = NBE.surfIntegral(mu,covar,xpop[J],ypop[I],0.5*delta)
                #print P
                if P>=1:
                    print 'Error:Probability of debris hit greater than 1==>casualtyEstimate'
                    print P
                    exit(1)
                Pt = P+Pt
                Ac = math.pi*(((Ap/math.pi)**0.5)+rp)**2
                Ec = Ec+(P*Ac*popDensity[I][J])*(0.001)**2
        '''
        E[index]=Ec

        #print Pt,Ec
   # print E
   # print Npieces
    Etotal = E*Npieces

    Efinal = np.sum(Etotal)    

  #  print 'Etotal =',Etotal
    print 'Efinal =',Efinal
        
    return (Efinal,XP,YP)
    


def expectationRAW(yList,filename,fragArea,Npieces,key,xllcorner,yllcorner,cellsize,nrows,xMax,XP,YP,xlen,ylen):
    nGroups = len(yList)
    gridArea = population.cellArea(filename)
    cfraction = .1
    delta = cfraction*population.gridSize(filename)
    popDensity = population.DensityfromRaw(XP,YP,xlen,ylen,key,xllcorner,yllcorner,cellsize,nrows,xMax)
    area = delta*delta
#    Ap = 1.0
    rp = .3048

    E = np.zeros(nGroups)
    Egrid=np.zeros([ylen,xlen])
 

    for index in range(0,nGroups):
        mu,covar = NBE.MLEestimator(yList[index])
        if mu =='NAN':
            f = 0.0*XP
        else:
            f = NBE.normal_bivariate(mu,covar,XP,YP)
        Ec=calculateE(f,area,fragArea[index],rp,popDensity)

        Egrid = Egrid+Ec*Npieces[index]
      #  Ec = Ac*P*popDensity
        Ec = np.sum(Ec)

        E[index]=Ec

    Etotal = E*Npieces

    Efinal = np.sum(Etotal)    

    print 'Efinal =',Efinal
        
    return (Efinal,Egrid)
    
def expectationKernelRAW(xpoplist,ypoplist,xlenlist,ylenlist,filename,y,key,xllcorner,yllcorner,cellsize,nrows,xMax,fragArea,Npieces):
    nGroups = len(Npieces)
    gridArea = population.cellArea(filename)
    cfraction = .1
    delta = cfraction*population.gridSize(filename)
    area = delta*delta
    rp = .3048
    E = np.zeros(nGroups)
 
    fkernelList = []
    for index in range(0,len(y)):
        if Npieces[index]>.1:
            yin = y[index]
          #  XL,XU,YL,YU = population.interestAreaKernel(yin)
          #  xpop = np.arange(XL,XU,delta)
          #  ypop = np.arange(YL,YU,delta)
          #  xlen = len(xpop)
          #  ylen = len(ypop)
            xpop=xpoplist[index]
            ypop=ypoplist[index]
            XP,YP = np.meshgrid(xpop,ypop)
            xlen = xlenlist[index]
            ylen = ylenlist[index]
     
            popDensity = population.DensityfromRaw(XP,YP,xlen,ylen,key,xllcorner,yllcorner,cellsize,nrows,xMax)
            #print popDensity
 
            fkernel = KE.fhatKernel(xpop,ypop,yin)
            Ec=calculateE(fkernel,area,fragArea[index],rp,popDensity)

            E[index]=Ec
  
            fkernelList.append(fkernel)
        else:
            E[index]=0.0
    Etotal = E*Npieces

    Efinal = np.sum(Etotal)
    print Efinal
   # return XPlist,YPlist,fkernelList,Efinal
    return Efinal
        
def calculateE(f,area,fragArea,rp,popDensity):
    P = f*area
    print np.sum(P)
    Ap = fragArea
    Ac = np.pi*(((Ap/np.pi)**0.5)+rp)**2
    Ec = Ac*(P*popDensity)*(0.001)**2
    Ec = np.sum(Ec)
    return Ec
