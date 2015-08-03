import numpy as np
# roof models for calculating casualty areas from debris weight

# Wood 
def getModel():
    from scipy import interpolate
    # data digitized from "Large Region Population Sheltering Models for Space Debris Risk Analysis" - Erik W. F. Larson
    # returns a list with the different roof models..ALL UNITS IN SI
    # data has been expanded to return logical values in case some extrapolation is needed when doing the expectation calculation
    woodWeight = [-10000.0,
                  0.0,
                  0.43,
                  0.43025,
                0.438297,
                0.438444,
                0.446672,
                0.454998,
                0.507711,
                0.576914,
                0.885167,
                1.42143,
                2.41076,
                2.97266,
                3.40845,
                6.82344,
                8.89053,
                11.0659,
                14.9491,
                18.7743,
                23.1548,
                31.57,
                39.6554,
                54.0648,
                69.7838,
                117.3,
                1323.68,
                  1e6,]
    woodCasualtyArea = [0.0,
                        0.0,
                        0.0,
                        0.103118,
                    0.130885,
                    0.168707,
                    0.224251,
                    0.271798,
                    0.347621,
                    0.427813,
                    0.43427,
                    0.447641,
                    0.457863,
                    0.472077,
                    0.542114,
                    2.31852,
                    3.37909,
                    4.22271,
                    5.03855,
                    5.61026,
                    6.69476,
                    8.69349,
                    11.1174,
                    13.8919,
                    16.1984,
                    18.7384,
                    33.8049,
                    33.8049]
    tileWeight = [-10000.0,
                  0.0,
                  0.56,
                  0.560304,
                  0.560532,
                  0.57627,
                  0.614399,
                  0.730918,
                  1.01474,
                  1.51506,
                  2.2619,
                  2.97288,
                  4.05176,
                  5.37443,
                  8.56994,
                  9.73856,
                  11.7972,
                  15.5054,
                  19.8274,
                  25.8257,
                  44.2355,
                  87.64,
                  1348.01,
                  1e6]
    tileCasualtyArea = [0.,
                        0.0,
                        0.,
                        0.102302,
                        0.13916,
                        0.183557,
                        0.232974,
                        0.369558,
                        0.434215,
                        0.461601,
                        0.464989,
                        0.498193,
                        0.482957,
                        0.546065,
                        2.70354,
                        3.45767,
                        5.0786,
                        5.44125,
                        5.31589,
                        6.10381,
                        11.2021,
                        18.0363,
                        33.8043,
                        33.8043]
    lightMetalWeight = [-10000.0,
                        0.0,
                        1.23,
                        1.2375,
                        1.55624,
                        2.10309,
                        3.73342,
                        5.99554,
                        6.6883,
                        9.03318,
                        11.6595,
                        18.759,
                        23.1409,
                        29.0639,
                        55.5127,
                        95.9119,
                        141.946,
                        208.2,
                        353.318,
                        588.587,
                        1231.17,
                        1e6]
    lightMetalCasualtyArea = [0.,
                              0.,
                              0.,
                              0.102226,
                              0.314208,
                              0.486987,
                              0.537913,
                              0.576216,
                              0.598748,
                              0.593992,
                              0.692614,
                              3.03196,
                              4.25236,
                              4.95852,
                              7.05931,
                              9.2356,
                              12.3667,
                              18.4424,
                              29.4704,
                              37.6775,
                              46.3426,
                              46.3426]
    compositeWeight = [-10000.,
                       0.0,
                       0.41,
                       0.414854,
                       0.438551,
                       0.592841,
                       0.869207,
                       1.81794,
                       4.6032,
                       14.1168,
                       16.9368,
                       22.6683,
                       27.701,
                       38.1261,
                       65.9359,
                       79.1347,
                       102.148,
                       131.827,
                       164.056,
                       225.733,
                       588.294,
                       680.694,
                       1300.31,
                       1e6]
    compositeCasualtyArea = [0.0,
                             0.0,
                             0.0,
                             0.102331,
                             0.202914,
                             0.399187,
                             0.441011,
                             0.494605,
                             0.517515,
                             0.681917,
                             0.666246,
                             0.681613,
                             0.753157,
                             1.24138,
                             3.37271,
                             4.28024,
                             5.18659,
                             5.38867,
                             5.95418,
                             7.91224,
                             25.8455,
                             29.2265,
                             45.9851,
                             45.9851]
    steelWeight = [-1000.,
                   0.0,
                   3.3,
                   3.31557,
                   4.16573,
                   5.62678,
                   6.81194,
                   9.80581,
                   30.3479,
                   37.4181,
                   49.186,
                   83.5684,
                   106.926,
                   178.157,
                   372.992,
                   508.535,
                   712.452,
                   1324.24,
                   1e6]
    steelCasualtyArea = [0.0,
                         0.0,
                         0.0,
                         0.430407,
                         0.662015,
                         0.70927,
                         0.646613,
                         0.641438,
                         0.864932,
                         0.832132,
                         0.985319,
                         3.87271,
                         5.91094,
                         8.61281,
                         20.8463,
                         26.4525,
                         30.6057,
                         46.6972,
                         46.6972]
    concreteWeight = [-1000.0,
                      0.0,
                      14.4,
                      14.4949,
                      19.2389,
                      23.9432,
                      31.1825,
                      64.6299,
                      84.9577,
                      107.652,
                      140.311,
                      176.443,
                      225.792,
                      372.745,
                      831.903,
                      1324.42,
                      1e6]
    concreteCasualtyArea = [0.0,
                            0.0,
                            0.0,
                            0.346524,
                            0.636118,
                            0.719284,
                            0.747302,
                            0.87099,
                            1.04732,
                            1.00758,
                            1.89283,
                            5.64166,
                            9.66405,
                            12.6439,
                            34.8765,
                            51.6085,
                            51.6085]
    wood = [woodWeight,woodCasualtyArea]
    tile = [tileWeight,tileCasualtyArea]
    lightMetal = [lightMetalWeight,lightMetalCasualtyArea]
    composite = [compositeWeight,compositeCasualtyArea]
    steel = [steelWeight,steelCasualtyArea]
    concrete = [concreteWeight,concreteCasualtyArea]
    
    fwood = interpolate.interp1d(wood[0],wood[1])
    ftile = interpolate.interp1d(tile[0],tile[1])
    flightMetal = interpolate.interp1d(lightMetal[0],lightMetal[1])
    fcomposite = interpolate.interp1d(composite[0],composite[1])
    fsteel = interpolate.interp1d(steel[0],steel[1])
    fconcrete = interpolate.interp1d(concrete[0],concrete[1])
    #return (fwood,ftile,flightMetal,fcomposite,fsteel,fconcrete)
    return (fwood,fsteel,fconcrete,fcomposite,flightMetal,ftile)

def Expectation(model,dist,distparam1,distparam2,Aproj,distparamextra = 0):
    from scipy import integrate
    import numpy as np
# this routine takes in the name of the roof model desired, the distibution of mass to be integrated with respect to,
# and parameters needed to define the distribution
    #fwood,fsteel,fconcrete,fcomposite,flightMetal,ftile= getModel()
    functions = getModel()
    if model.lower()=='all':
    # HERE COMES AN ASSUMPTION TO GIVE ROOF CASUALTY AREA TO ROOF TYPES FOR WHICH DATA IS NOT AVAILABLE:
        # assuming that all roof of similar types have the same casualty areas. e.g all wood roofs are the same
    # the order of models is assumed to be as follows
    # Wood Roof, Wood 1st, Wood 2nd, Steel Roof, Steel 1st, Steel 2nd, Concrete, Concrete 1st, Concrete 2nd, Composite
    # Light Metal, Tile Roof, Tile 1st, Tile 2nd, Car, Open
        
        if dist.lower() =='uniform':
            lowbound = distparam1
            highbound = distparam2
            denc = highbound-lowbound
            if denc>0.:
                c = 1./denc
            elif denc<0.:
                print 'Check bounds in upper and lower bound. roofPenetration.Expectation'
                exit(1)

        elif dist.lower() =='gaussian' or dist.lower() == 'normal' or dist.lower()=='gauss':
            mu = distparam1
            sigma = distparam2
            lowbound = max(-1e3,mu-5.*sigma)
            highbound = min(1e6,mu+5.*sigma)
            fgauss = lambda x: np.exp(-.5*((x-mu)/sigma)**2)/(sigma*(2.*np.pi)**.5)
        
        EretMain = np.zeros((len(functions)))
        for index in range(len(functions)):
            froof = functions[index]
        
            if dist.lower() =='uniform':
                if denc ==0.0:
                    retval = froof(lowbound)
                else:
                    invals = integrate.fixed_quad(froof,lowbound,highbound,n=50)
                    retval = invals[0]*c
            elif dist.lower() =='gaussian' or dist.lower() == 'normal' or dist.lower()=='gauss':
                fvals = lambda x: fgauss(x)*froof(x)
                invals = integrate.fixed_quad(fvals,lowbound,highbound,n=50)
                retval = invals[0]
                
            EretMain[index] = retval
        Eret = np.zeros((16))
        Eret = np.concatenate((3*[EretMain[0]],3*[EretMain[1]],3*[EretMain[2]],[EretMain[3]],[EretMain[4]],3*[EretMain[5]],[EretMain[4]],[Aproj]),1)
        return Eret
        
#Expectation('all','uniform',10,100)

def getAcSheltering(model=None,mass=None,Aref=None,AF = 1.0):
    #this routine returns the casualty area for different roofs for a given debris with a certain mass and reference area. Model contains all the required functions to interpolate. Model can be obtained frm the function getModel() in this file.
    rp = .3048 #1 foot 
    Aproj = np.pi*(np.sqrt(AF*Aref/np.pi) + rp )**2.0
    AcArray = np.concatenate((3*[model[0](mass)],3*[model[1](mass)],3*[model[2](mass)],[model[3](mass)],[model[4](mass)],3*[model[5](mass)],[model[4](mass)],[Aproj]),1)
    return AcArray

        
class casualtyAreaClass:
    def __init__(self,shelter=True):
        
        self.shelter = False
        if shelter==True:
            self.shelter=True
            self.myShelterClass = shelteringClass()
        
    def getAreaInertDebris(self,mass=None,Af=None,Aref=None,CD=None,AFF=1.0):
        #assuming that fragment projected area is the same as reference area for ballistic coeff
        #AFF is a multiplier to the fragment projected area. Mainly used to model secondary effects
        if Af==None:
            Af = Aref
        rp = 0.3048
        Aopen = np.pi*(np.sqrt(AFF*Af/np.pi)+rp)**2.0
        Acas = []
        if self.shelter==True:
            Acas = self.myShelterClass.sampleAllGroups(mass,Aref,CD)
        return np.array(np.concatenate((Acas,[Aopen]),1))
#    def getAreaExplosiveDebris(self,massTNT


def sampling(modelClass,mass,ballisticCoeff):
    
    index = np.arange(len(modelClass.bound))
    
    index = index[modelClass.bound<ballisticCoeff]
    indexMax = len(modelClass.bound)-1
    nlen = len(index)
    if nlen==0:
        return 0.0
    elif nlen==len(modelClass.bound):
        index = index[-2]
    else:
        index = index[-1]
    Acas = modelClass.fun[index](mass)
    return Acas

class shelteringClass:
    def __init__(self):
        self.getAllModels()
            
    def getAllModels(self):
        gpA = ballClass()
        gpB = ballClass()
        gpC = ballClass()
        gpD = ballClass()
        
        gpA.getModelLightMetal()
        gpB.getModelComposite()
        gpC.getModelWood()
        gpD.getModelConcreteSteel()       
        self.groupA = gpA
        ####################
        #changing these lines becuase they are not in the right order in FAA handbook
        #self.groupB = gpB
        #self.groupC = gpC
        self.groupB = gpC
        self.groupC = gpB
        #####################

        self.groupD = gpD
    
    def sampleAllGroups(self,massIN,ArefIN,CD):
        #handling the sampling so that SI is used for inputs outputs
        lb_to_kg = 0.453592
        kg_to_lb = 1./lb_to_kg
        ft_to_m = 0.3048
        m_to_ft = 1.0/ft_to_m
        m2_to_ft2 = m_to_ft**2.0
        ft2_to_m2 = ft_to_m**2.0
        
        
        mass = massIN*kg_to_lb
        Aref = ArefIN*m2_to_ft2
        ballisticCoeff = mass/(CD*Aref)
        
        AreaA = sampling(self.groupA,mass,ballisticCoeff)
        #if mass<100.:
        AreaB = sampling(self.groupB,mass,ballisticCoeff)
        AreaC = sampling(self.groupC,mass,ballisticCoeff)
        AreaD = sampling(self.groupD,mass,ballisticCoeff)
    
        Acas = ft2_to_m2*np.array([AreaA,AreaB,AreaC,AreaD])
        return Acas
        #return np.array([AreaA])
        
def myInterp(x,y,leftVal=0.0,rightVal=0.0):
    # this function to avoid making mistakes with lambda functions in for loops
    
    #fun = lambda x0:np.interp(x0,x,y,leftVal,right=rightVal)
    fun = lambda x0:mylogInterp(x0,x,y)
    return fun

def mylogInterp(x0,x,y):
    # function to do a simple log log linear interpolation
    newx = np.log10(x)
    newy = np.log10(y)
    nx0 = np.log10(x0)
    leftVal = -9.
    rightVal = newy[-1]
    val = np.interp(nx0,newx,newy,left=leftVal,right=rightVal)

    if val<-8.:
        return 0.0
    val = 10.**val
    return val




class ballClass:
    def __init__(self):
        self.bound = []
        self.mass = []
        self.casArea = []
        self.name = []
        self.fun = []
        
    def getModelLightMetal(self):
        # returns the required list of values for sheltering. It assumes 4 sheltering categories. 
        #data obtained from FLigh Safety Analysis Handbook  Pg 116
        self.bound = []
        self.mass = []
        self.casArea = []
        self.name = []
        self.fun = []
        #group A (Light Metal Roof)
        self.name = 'Light Metal'
        self.bound = [3,10,17.5,30,55,69.9]
        self.mass.append([289.427,
                          991.397,
                          3048.26,
                          9372.54])
        
        self.casArea.append([11.7125,
                             150.038,
                             468.253,
                             542.692])
        
        self.mass.append([9.95689,
                          29.9604,
                          102.626,
                          302.204,
                          1013.04,
                          2919.38,
                          5119.1])
            
            
        self.casArea.append([0.465792,
                             3.75293,
                             21.1321,
                             181.377,
                             327.246,
                             578.114,
                             958.723])

        self.mass.append([2.9068,
                           9.95689,
                           29.9604,
                           100.433,
                           302.204,
                           991.397,
                           2919.38,
                           9786.31])

        self.casArea.append([0.994745,
                             5.14838,
                             6.35629,
                             27.793,
                             197.331,
                             348.605,
                             566.057,
                             919.151])
        self.mass.append([0.300901,
                          0.987123,
                          9.74413,
                          29.9604,
                          100.433,
                          271.266,
                          1013.04,
                          3048.26,
                          9577.19])
        self.casArea.append([0.446566,
                             4.17001,
                             4.83293,
                             7.8476,
                             35.0447,
                             98.4317,
                             307.195,
                             478.227,
                             760.338])
        self.mass.append([0.294471,
                          1.00868,
                          9.95689,
                          29.9604,
                          98.2868,
                          308.803,
                          991.397,
                          3048.26,
                          9577.19])
        self.casArea.append([3.83287,
                             4.08304,
                             4.83293,
                             7.8476,
                             37.3321,
                             100.528,
                             282.359,
                             458.487,
                             728.954])

    
        for index in range(len(self.mass)):
            self.fun.append(myInterp(self.mass[index],self.casArea[index],leftVal=0.0,rightVal=self.casArea[index][-1]))
        
        

        #self.fun =[ lambda index=index, x:myInterp(x,xin,self.casArea[index],leftVal=0.0,rightVal=self.casArea[index][-1])]

  
        self.mass = np.array(self.mass)
        self.casArea = np.array(self.casArea)
        self.bound = np.array(self.bound)
    def getModelComposite(self):
        self.bound = []
        self.mass = []
        self.casArea = []
        self.name = []
        self.fun = []
        #group B (Composite Roof)
        self.name = 'Composite Roof'
        self.bound = [3,10,17.5,30,55,100.]
        self.mass.append([10.0866,
                          30.2882,
                          99.1413])
        self.casArea.append([0.142856,
                             5.05676,
                             15.7])

        self.mass.append([1.02621,
                          3.08151,
                          10.0866,
                          32.312,
                          101.302,
                          310.82])
        self.casArea.append([3.12108,
                             5.49945,
                             6.92702,
                             8.72518,
                             11.4611,
                             16.7199])
        self.mass.append([0.300281,
                          1.07143,
                          3.01579,
                          10.3064,
                          99.1413,
                          310.82,
                          1039.57,
                          3121.63])
                        
        self.casArea.append([0.262503,
                             4.55317,
                             5.05676,
                             5.85671,
                             10.3197,
                             24.9089,
                             78.9759,
                             87.7107])
        self.mass.append([0.300281,
                          0.9829,
                          3.08151,
                          10.3064,
                          29.6422,
                          101.302,
                          304.191,
                          1039.57,
                          3189.67,
                          10440.6])
        self.casArea.append([4.95177,
                             4.09974,
                             4.95177,
                             5.05676,
                             6.23719,
                             9.89565,
                             27.0896,
                             91.4694,
                             190.626,
                             643.659])

        self.mass.append([0.300281,
                          1.04857,
                          3.21729,
                          10.0866,
                          32.312,
                          101.302,
                          304.191,
                          1062.23,
                          3121.63,
                          10000.])
        self.casArea.append([4.74829,
                              4.45864,
                              4.55317,
                              5.16397,
                              6.23719,
                              9.89565,
                              27.0896,
                              84.1064,
                              302.439,
                              685.474])
    
        for index in range(len(self.mass)):
            self.fun.append(myInterp(self.mass[index],self.casArea[index],leftVal=0.0,rightVal=self.casArea[index][-1]))
    
        self.mass = np.array(self.mass)
        self.casArea = np.array(self.casArea)
        self.bound = np.array(self.bound)
        
        
        
    def getModelWood(self):
        self.bound = []
        self.mass = []
        self.casArea = []
        self.name = []
        self.fun = []
        
        
        #group C (Wood Roof)
        self.name = 'Wood Roof'
        self.bound = [3,10,17.5,30,55,100.]
        self.mass.append([9.70911,
                          30.3751,
                          67.041,
                          96.3126,
                          289.427,
                          1021.7,
                          3070.29])
        self.casArea.append([4.01281,
                             3.2044,
                             9.61078,
                             16.1026,
                             66.3497,
                             233.247,
                             412.043])
        self.mass.append([1.01896,
                          3.02125,
                          9.9732,
                          30.7854,
                          100.269,
                          309.511,
                          1035.5,
                          2949.15,
                          4069.65])
        self.casArea.append([2.27155,
                             5.0921,
                             4.82957,
                             17.9008,
                             98.6854,
                             193.801,
                             412.043,
                             777.687,
                             986.854])
        self.mass.append([0.393012,
                          1.00538,
                          2.94125,
                          9.84027,
                          30.3751,
                          98.9323,
                          1008.08,
                          3029.37,
                          6960.77])
        self.casArea.append([0.0986854,
                             4.58057,
                             4.76608,
                             6.04796,
                             23.6354,
                             120.354,
                             356.225,
                             509.21,
                             1000])
        self.mass.append([0.308681,
                          1.01896,
                          2.90204,
                          9.84027,
                          29.9702,
                          101.623,
                          309.511,
                          994.647,
                          2949.15,
                          9735.2])
        self.casArea.append([4.17532,
                             3.85662,
                             4.23094,
                             7.57374,
                             26.2748,
                             88.772,
                             188.739,
                             295.981,
                             428.73,
                             747.417])
        self.mass.append([0.304567,
                          1.01896,
                          3.02125,
                          9.84027,
                          29.9702,
                          100.269,
                          305.386,
                          1021.7,
                          3070.29,
                          9605.45])
        self.casArea.append([3.80592,
                             3.85662,
                             4.17532,
                             7.98544,
                             28.0722,
                             77.7687,
                             193.801,
                             255.886,
                             380.592,
                             529.832])
        
        for index in range(len(self.mass)):
            self.fun.append(myInterp(self.mass[index],self.casArea[index],leftVal=0.0,rightVal=self.casArea[index][-1]))
        self.mass = np.array(self.mass)
        self.casArea = np.array(self.casArea)
        self.bound = np.array(self.bound)



    def getModelConcreteSteel(self):
        self.bound = []
        self.mass = []
        self.casArea = []
        self.name = []
        self.fun = []
        
        
        #group D (concrete reinforced Roof)
        self.name = 'concrete reinforced with steel'
        self.bound = [10,17.5,30,55,100.]
        self.mass.append([99.7186,
                          299.318,
                          991.581,
                          3018.59,
                          9722.1])
        self.casArea.append([7.61477,
                             8.72626,
                             15.2557,
                             25.9538,
                             152.557])


        self.mass.append([29.6797,
                          102.569,
                          299.318,
                          991.581,
                          3061.43,
                          8933.89])

        self.casArea.append([8.49169,
                             16.3312,
                             16.5552,
                             38.53,
                             229.587,
                             741.007])

        self.mass.append([2.94299,
                          10.1705,
                          29.6797,
                          98.3232,
                          295.129,
                          1005.65,
                          3018.59,
                          9860.07])
        self.casArea.append([1.9232,
                             5.27101,
                             7.01704,
                             9.46959,
                             17.7223,
                             100,
                             369.868,
                             884.597])
        self.mass.append([0.304422,
                          1.0228,
                          3.07007,
                          10.3149,
                          30.5281,
                          98.3232,
                          303.565,
                          991.581,
                          3104.87,
                          9860.07])
        self.casArea.append([0.114597,
                             3.64863,
                             4.35565,
                             4.85725,
                             6.04039,
                             8.37678,
                             20.0343,
                             97.3118,
                             423.857,
                             837.678])
        for index in range(len(self.mass)):
            self.fun.append(myInterp(self.mass[index],self.casArea[index],leftVal=0.0,rightVal=self.casArea[index][-1]))

        self.mass = np.array(self.mass)
        self.casArea = np.array(self.casArea)
        self.bound = np.array(self.bound)





def groundBlastFunctions():
    # interpolated models  from Flight Safety Analysis Handbook Version 1.0, FAA. This returns the casualty area for people in 4 different shelter categories
    # handles casualty areas for explosive pieces of debris
    from scipy import interpolate
    lb_to_kg = 0.453592
    sqft2sqm = 0.092903
    
    Ayield = lb_to_kg*np.array([0.0,9.9858,21.0768,43.4598,124.259,359.45,1052.05,3008.29,8805.42,25475.4,72848.4,306233,1243040,7417370])
    classA = sqft2sqm*np.array([0.0,388.113,1101.46,3047.39,11738,45213.9,157333,470049,1369150,3695260,9973080,33855800,103817000,400428000])
    fA = interpolate.interp1d(Ayield,classA)
    
    Byield = lb_to_kg*np.array([0.0,9.98432,12.3186,17.0801,27.5632,71.7819,172.274,540.811,1510.69,4269.58,10991.2,38778.9,183196,1019130,2446250,9476860])
    classB = sqft2sqm*np.array([0.0,568.166,771.005,1250.18,2361.67,8010.13,23923.4,85400.6,268432,801951,2056740,7159410,29041200,127171000,252940000,700615000])
    fB = interpolate.interp1d(Byield,classB)
    
    Cyield = lb_to_kg*np.array([0.0,9.99003,24.8286,80.7106,319.969,1210.67,2496.51,5208.52,10250.8,21638,43592.6,85796.3,157438,328492,685407,1413540,2949450,7246990])
    classC = sqft2sqm*np.array([0.0,130.159,476.389,2489.69,14779.5,75323.8,178930,425054,889171,2007640,3991880,7737760,13547500,26262900,48390400,84741500,148403000,295184000])
    fC = interpolate.interp1d(Cyield,classC)
    
    Dyield = lb_to_kg*np.array([0.0,10.,13.0665,23.1491,85.5603,327.527,1239.25,4970.95,19940.3,83815.6,328496,1413530,7247060])
    classD = sqft2sqm*np.array([0.0,120.609,185.854,441.369,2756.33,15160.5,79253.8,393825,1813360,7737420,25604100,86922200,287778000])
    fD = interpolate.interp1d(Dyield,classD)
    
    return (fA,fB,fC,fD)


class casAreaExplosiveDebris:

    # this class is to determine casualty areas due to different propellant types or inert debris hitting the ground
    def __init__(self):
        import blastPredict as bp

        self.inertCasAreaClass = casualtyAreaClass()
        
        propLiquidTypes = ['lo2/lh2','lo2/rp-1','hypergols']
        surfaceTypes = ['soft','hard']
        (fA,fB,fC,fD) = groundBlastFunctions()
        liqFun = [[],[]]

        self.solidProp  = bp.solidPropellantYield
        self.getSimpleRadiusTNT = bp.getSimpleRadiusTNT
        yieldLiquidClass1 = bp.yieldLiquid()
        yieldLiquidClass2 = bp.yieldLiquid()
        yieldLiquidClass3 = bp.yieldLiquid()
        yieldLiquidClass4 = bp.yieldLiquid()
        yieldLiquidClass5 = bp.yieldLiquid()
        yieldLiquidClass6 = bp.yieldLiquid()
        yieldLiquidClass = [yieldLiquidClass1,yieldLiquidClass2,yieldLiquidClass3,
                            yieldLiquidClass4,yieldLiquidClass5,yieldLiquidClass6]
        counter = 0
        for index1 in range(len(surfaceTypes)):
            for index2 in range(len(propLiquidTypes)):

                #print 'LQ', yieldLiquidClass.propellantType(propLiquidTypes[index2],surfaceTypes[index1])
                yieldLiquidClass[counter].propellantType(propLiquidTypes[index2],surfaceTypes[index1])
                liqFun[index1].append(yieldLiquidClass[counter])
                counter = counter + 1

        self.liqFun = liqFun
        self.propLiquidTypes = propLiquidTypes
        self.fA = fA
        self.fB = fB
        self.fC = fC
        self.fD = fD
    
    def getCasArea(self,propellant=None,surface=None,v=None,mass=None,aref=None,CD=None,TNT_K_factor=45.):
        if propellant.lower()=='solid':
            #getting yield mass
            massTNT = self.solidProp(mass = mass,S=surface,v = v)
        elif propellant.lower()=='inert':
            Aout = self.inertCasAreaClass.getAreaInertDebris(mass=mass,Aref=aref,CD=CD)
            return np.array(Aout) #returning cas area if inert debris
        else:
            index1 = 0
            if surface=='hard':
                index1 = 1
            index2 = np.arange(len(self.propLiquidTypes))
            index2 = index2[self.propLiquidTypes==propellant.lower()]

            yieldTNT = self.liqFun[index1][index2].yieldFactor(v)
            massTNT = mass*yieldTNT
        R_TNT = self.getSimpleRadiusTNT(massTNT=massTNT,TNT_K_factor=TNT_K_factor)
        AcasOpen = np.pi*R_TNT**2.0
        if massTNT>=4.53592: # this check needed because actual model does not show results for TNT mass <10 lb (~4.5kg)
            #using inert debris if less than 10 lbs for TNT mass
            Aout = np.array([self.fA(massTNT),self.fB(massTNT),self.fC(massTNT),self.fD(massTNT),AcasOpen])
        else:
            Aout = self.inertCasAreaClass.getAreaInertDebris(mass=mass,Aref=aref,CD=CD)
            Aout[-1] = AcasOpen #making sure explosive debris is used for people in the open instead of inert debris    
        return np.array(Aout)










'''


import matplotlib.pyplot as plt

myShelterClass = shelteringClass()

#massArray = np.array([.5,1,5,10,15,100,500,1000,10000])# np.linspace(1.,10000.,10000)
massArray = np.linspace(np.log10(.1),np.log10(10000.))
massArray =.453592*10**massArray
ball = [5.,15.,20.,50.,60.]

A1 = []
A2 = []
A3 = []
A4 = []
A5 = []
pltindex = 3

CD = 1.0

for index in range(len(massArray)):
    Ared = ball
    Array0 = myShelterClass.sampleAllGroups(massArray[index],ball[0])
    Array1 = myShelterClass.sampleAllGroups(massArray[index],ball[1])
    Array2 = myShelterClass.sampleAllGroups(massArray[index],ball[2])
    Array3 = myShelterClass.sampleAllGroups(massArray[index],ball[3])
    Array4 = myShelterClass.sampleAllGroups(massArray[index],ball[4])

    A1.append(Array0[pltindex])
    A2.append(Array1[pltindex])
    A3.append(Array2[pltindex])
    A4.append(Array3[pltindex])
    A5.append(Array4[pltindex])
plt.figure()
plt.plot(massArray,A1,label='BC Class 2: 3 to 10 psf')
plt.figure()
plt.loglog(massArray,A2,label='BC Class 3: 10 to 17.5 psf')
plt.loglog(massArray,A3,label='BC Class 4: 17.5 to 30 psf')

plt.loglog(massArray,A4,label='BC Class 5: 30 to 55 psf')
plt.loglog(massArray,A5,label='BC Class 6: 55 to 69.9 psf')
plt.legend(loc='upper left')
plt.grid()

plt.show()

'''


