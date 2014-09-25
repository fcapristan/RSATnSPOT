import openturns as ot


from numpy import empty, argmin, array, arange, floor, sqrt, linspace
#
#===================================================================
#                         Marginal PDFs
#===================================================================
weights = [-6.2,.5,10.]
def setInputRandomVector(inList):
    #setting the input random vector X
    myCollection = ot.DistributionCollection(len(inList))
    for index in range(len(inList)):
        myCollection[index] = ot.Distribution(inList[index])
    myDistribution = ot.ComposedDistribution(myCollection)
    VectX = ot.RandomVector(ot.Distribution(myDistribution))
    return myDistribution,VectX

def myfun1(X):
    [X0,X1] = X
    return [.01*(10*X0 + 5*X1)**2.]
  
def myfun2(X):
    [X2,X3,X4] = X
    return [(X2+0.*X3-.1*X2+.7*X4)**3.0]
def myfun3(X):
    [X5,X6]=X
    return [X5 + X6]
    
def myfunTotal(X):
    [X0,X1,X2,X3,X4,X5,X6] = X
    v1 = myfun1([X0,X1])
    v2 = myfun2([X2,X3,X4])
    v3 = myfun3([X5,X6])
    #print 'Eval'
    
    return [ weights[0]*v1[0]+ weights[1]*v2[0]+weights[2]*v3[0]]

def sobolTotal(nCases):

    X0 = ot.Uniform(0.0,1.0)
    X1 = ot.Uniform(0.0,1.0)
    X2 = ot.Uniform(0.0,1.0)
    X3 = ot.Uniform(0.0,1.0)
    X4 = ot.Uniform(0.0,1.0)
    X5 = ot.Uniform(0.0,1.0)
    X6 = ot.Uniform(0.0,1.0)

    inList = [X0,X1,X2,X3,X4,X5,X6]

    dim = len(inList) #number of random variables
    myDistribution,VectX = setInputRandomVector(inList)



    f = myfunTotal
    model = ot.PythonFunction(dim,1,f)





    ot.RandomGenerator.SetSeed(10)
    myRandomExp1 = ot.LHSExperiment(myDistribution,nCases)
    inp1 = myRandomExp1.generate()

    ot.RandomGenerator.SetSeed(75)
    myRandomExp2 = ot.LHSExperiment(myDistribution,nCases)
    inp2 = myRandomExp2.generate()

    #myRandomExp3 = LHSExperiment(myDistribution,nCases)
    #inp3 = myRandomExp3.generate()
    sensAnalysis = ot.mySensitivityAnalysis(inp1,inp2,model)

    sensAnalysis.setBlockSize(1)
    print '########################'
    print 'Starting Computations'
    #secondOrderIndices = sensAnalysis.getSecondOrderIndices()
    firstOrderIndices = sensAnalysis.getFirstOrderIndices()
    totalOrderIndices = sensAnalysis.getTotalOrderIndices()
    print 'Done withSobol'
    out1 = sensAnalysis.getOutVals1()
    print 'Mean', out1.computeMean()
    print 'STD', out1.computeStandardDeviation()
    out2 = sensAnalysis.getOutVals2()
    
    print 'First',firstOrderIndices
    print 'Total',totalOrderIndices
    #print 'Second',secondOrderIndices
    return [firstOrderIndices,totalOrderIndices]#,secondOrderIndices,out1]
    
def sobolFun2(nCases):
    X2 = ot.Uniform(0.0,1.0)
    X3 = ot.Uniform(0.0,1.0)
    X4 = ot.Uniform(0.0,1.0)
    inList=[X2,X3,X4]
    dim = len(inList) #number of random variables
    myDistribution,VectX = setInputRandomVector(inList)
    ot.RandomGenerator.SetSeed(2)
    myRandomExp1 = ot.LHSExperiment(myDistribution,nCases)
    inp1 = myRandomExp1.generate()

    ot.RandomGenerator.SetSeed(5)
    myRandomExp2 = ot.LHSExperiment(myDistribution,nCases)
    inp2 = myRandomExp2.generate()
    
    f = myfun2
    model = ot.PythonFunction(dim,1,f)
    sensAnalysis = ot.mySensitivityAnalysis(inp1,inp2,model)
    #secondOrderIndices = sensAnalysis.getSecondOrderIndices()
    firstOrderIndices = sensAnalysis.getFirstOrderIndices()
    totalOrderIndices = sensAnalysis.getTotalOrderIndices()
    out1 = sensAnalysis.getOutVals1()
    out2 = sensAnalysis.getOutVals2()
    return [firstOrderIndices,totalOrderIndices,out1]
def sobolFun1(nCases):
    X0 = ot.Uniform(0.0,1.0)
    X1 = ot.Uniform(0.0,1.0)
    inList=[X0,X1]
    dim = len(inList) #number of random variables
    myDistribution,VectX = setInputRandomVector(inList)
    ot.RandomGenerator.SetSeed(1)
    myRandomExp1 = ot.LHSExperiment(myDistribution,nCases)
    inp1 = myRandomExp1.generate()

    ot.RandomGenerator.SetSeed(4)
    myRandomExp2 = ot.LHSExperiment(myDistribution,nCases)
    inp2 = myRandomExp2.generate()
    
    f = myfun1
    model = ot.PythonFunction(dim,1,f)
    sensAnalysis = ot.mySensitivityAnalysis(inp1,inp2,model)
    #secondOrderIndices = sensAnalysis.getSecondOrderIndices()
    firstOrderIndices = sensAnalysis.getFirstOrderIndices()
    totalOrderIndices = sensAnalysis.getTotalOrderIndices()
    out1 = sensAnalysis.getOutVals1()
    out2 = sensAnalysis.getOutVals2()
    return [firstOrderIndices,totalOrderIndices,out1]#,secondOrderIndices,out1]    
def sobolFun3(nCases):
    X0 = ot.Uniform(0.0,1.0)
    X1 = ot.Uniform(0.0,1.0)
    inList=[X0,X1]
    dim = len(inList) #number of random variables
    myDistribution,VectX = setInputRandomVector(inList)
    ot.RandomGenerator.SetSeed(3)
    myRandomExp1 = ot.LHSExperiment(myDistribution,nCases)
    inp1 = myRandomExp1.generate()

    ot.RandomGenerator.SetSeed(6)
    myRandomExp2 = ot.LHSExperiment(myDistribution,nCases)
    inp2 = myRandomExp2.generate()
    
    f = myfun3
    model = ot.PythonFunction(dim,1,f)
    sensAnalysis = ot.mySensitivityAnalysis(inp1,inp2,model)
    #secondOrderIndices = sensAnalysis.getSecondOrderIndices()
    firstOrderIndices = sensAnalysis.getFirstOrderIndices()
    totalOrderIndices = sensAnalysis.getTotalOrderIndices()
    out1 = sensAnalysis.getOutVals1()
    out2 = sensAnalysis.getOutVals2()
    return [firstOrderIndices,totalOrderIndices,out1]   
    
    
#sobolTotal(1000000)
'''
[ind1,indtot1,out1] = sobolFun1(1000)
[ind2,indtot2,out2] = sobolFun2(1000)
[ind3,indtot3,out3] = sobolFun3(1000)

stdTotal = (out1.computeStandardDeviation()[0,0]**2 + out2.computeStandardDeviation()[0,0]**2 + out3.computeStandardDeviation()[0,0]**2)**.5
std1 = out1.computeStandardDeviation()[0,0]
print stdTotal

in0new = ind1[0]*std1**2/(stdTotal**2)
in0Tot = indtot1[0]*std1**2/(stdTotal**2) #(stdTotal**2.0 + (std1**2.0)*indtot1[0] - std1**2.0)/(stdTotal**2.0)
print 'X0',in0new,in0Tot


in1new = ind1[1]*std1**2/(stdTotal**2)
in1Tot = indtot1[1]*std1**2/(stdTotal**2)
print 'X1', in1new,in1Tot

std2 = out2.computeStandardDeviation()[0,0]
in2new = ind2[0]*std2**2/(stdTotal**2)
in2tot = indtot2[0]*std2**2/(stdTotal**2)
in3new = ind2[1]*std2**2/(stdTotal**2)
in4new = ind2[2]*std2**2/(stdTotal**2)
print 'X2', in2new,in2tot
print 'X3', in3new
print 'X4', in4new


std3 = out3.computeStandardDeviation()[0,0]
in5new = ind3[0]*std3**2/(stdTotal**2)
in6new = ind3[1]*std3**2/(stdTotal**2)
print 'X5',in5new
print 'X6',in6new

'''
import OTsobolWrapper as OTw
wInList = weights
X0 = ot.Uniform(0.0,1.0)
X1 = ot.Uniform(0.0,1.0)
X2 = ot.Uniform(0.0,1.0)
X3 = ot.Uniform(0.0,1.0)
X4 = ot.Uniform(0.0,1.0)
X5 = ot.Uniform(0.0,1.0)
X6 = ot.Uniform(0.0,1.0)

OTinputList = [[X0,X1],[X2,X3,X4],[X5,X6]]
seedList1 = [1,2,3]
seedList2 = [4,5,6]
nCases = 50000
blocksize = 1

model1 = ot.PythonFunction(2,1,myfun1)
model2 = ot.PythonFunction(3,1,myfun2)
model3 = ot.PythonFunction(2,1,myfun3)
modelList = [model1,model2,model3]
(X1,X2) =OTw.wrapSobol(wInList=wInList,OTmodelListIn=modelList,OTinputList=OTinputList,
              seedList1=seedList1,seedList2=seedList2,nCases=nCases,blocksize=blocksize)
           
print X1
print X2
print len(X1),len(X2)






