## these routines used for sheltering set up and computations
import numpy as np
def read(fileName):
    table1 = False
    table2 = False
    table3 = False
    vector1 = False
    vector2 = False
    school = False
    fraction = False
    fileIN = open(fileName,'r')
    table1Matrix = np.array([])
    table2Matrix = np.array([])
    table3Matrix = np.array([])

    
    for line in fileIN:
        key = line.split()
        if len(key)>0:
            if key[0].lower()=='$endtable1':
                table1 = False
            elif key[0].lower()=='$endtable2':
                table2 = False
            elif key[0].lower()=='$endtable3':
                table3 = False       
            elif key[0].lower()=='$endvector1':
                vector1 = False          
            elif key[0].lower()=='$endvector2':
                vector2 = False          
            elif key[0].lower() =='$endschoolvec':
                school = False
            elif key[0].lower() =='$endfractionvals':
                fraction = False
        
        if table1 or table2 or table3 or vector1 or vector2 or school or fraction:
            # table with inputs MUST FOLLOW GIVEN FORMAT
            if line.find(":") > 0:
                iColon = line.index(":")
                key = (line[:iColon].lower()).split()
                value = (line[iColon+1:]).split()
            if table1:
                lentable1 = len(value)
                table1Matrix = np.concatenate((table1Matrix,np.array(value,dtype=float)),1)
            elif table2:
                lentable2 = len(value)
                table2Matrix = np.concatenate((table2Matrix,np.array(value,dtype=float)),1)
            elif table3:
                lentable3 = len(value)
                table3Matrix = np.concatenate((table3Matrix,np.array(value,dtype=float)),1)
            elif vector1 :
                vector1Array = np.array(value,dtype=float)
            elif vector2 :
                vector2Array = np.array(value,dtype=float)
            elif school:
                schoolvec = np.array(value,dtype=float)
            elif fraction :
                fractionvals = np.array(value,dtype=float)
                
        if len(key)>0:
            if key[0].lower()=='$table1':
                table1 = True
            elif key[0].lower()=='$table2':
                table2 = True
            elif key[0].lower()=='$table3':
                table3 = True
            elif key[0].lower()=='$vector1':
                vector1 = True          
            elif key[0].lower()=='$vector2':
                vector2 = True     
            elif key[0].lower() =='$schoolvec':
                school = True
            elif key[0].lower() == '$fractionvals':
                fraction = True
                
    mat1 = np.reshape(table1Matrix,[len(table1Matrix)/lentable1,lentable1])
    mat2 = np.reshape(table2Matrix,[len(table2Matrix)/lentable2,lentable2])
    mat3 = np.reshape(table3Matrix,[len(table3Matrix)/lentable3,lentable3])

    return (mat1,mat2,mat3,vector1Array,vector2Array,schoolvec,fractionvals)


def getVals(fileName,seasonAndTime):
    # Season and time options:
    # "weekday-daytime-summer" => 0
    # "weekday-daytime-winter" => 1
    # "weekday-night"          => 2
    # "weekend-daytime-summer" => 3
    # "weekend-daytime-winter" => 4
    # "weekend-night"          => 5
    mat1,mat2,mat3,vector1,vector2,schoolvec,fractionvals = read(fileName)
    q = schoolvec
    # getting e2,s2,,d,v from seasonAndTime
    s2 = mat3[seasonAndTime,0]
    e2 = mat3[seasonAndTime,1]
    d = mat3[seasonAndTime,2]
    v = mat3[seasonAndTime,3]
    
    e1 = fractionvals[0]
    s1 = fractionvals[1]
    term1 = e1*e2*np.dot(mat2.T,vector2[np.newaxis].T)
    term2 = s1*s2*q
    term2 = term2[np.newaxis].T
    term3 = (1.-d-v)*np.dot(mat1.T,vector1[np.newaxis].T)
    term3 = np.concatenate((term3,[[0]]),0)
    #print np.shape(term3)
    term3[-1] = term3[-1] + 100.*d
    term3[-2] = term3[-2] + 100.*v
    term3 = (1.-e1*e2-s1*s2)*term3
    c = term1+term2+term3
    #print np.sum(np.dot(mat2.T,vector2[np.newaxis].T)) , np.sum(q), np.sum(np.dot(mat1.T,vector1[np.newaxis].T)),np.sum(c),e1,e2,s1,s2,np.sum(term3)
    #print c
    #print np.sum(c)
    #print mat1
    #print mat2
    #print mat3
    #print vector1
    #print vector2
    #print schoolvec
    #print term1,np.sum(term1),e1*e2,e1,e2
#print term3
    return c

def getVals5Groups(c):
    #the new groups are
    # the order of models is assumed to be as follows
    # Wood Roof, Wood 1st, Wood 2nd, Steel Roof, Steel 1st, 
    # Steel 2nd, Concrete, Concrete 1st, Concrete 2nd, Composite
    # Light Metal, Tile Roof, Tile 1st, Tile 2nd, Car,
    #Open
    
    # Wood Roof, Wood 1st, Wood 2nd,Tile 1st, Tile 2nd -> assumes group 3
    #  Steel Roof, Steel 1st, Steel 2nd, Concrete, Concrete 1st, Concrete 2nd -> assumes group 4
    # composite -> group 2
    # light metal , cars-> group 1
    # group 5...people in the open
    c1 = c[10] + c[14]#light metal and cars
    c2 = c[9] #composite
    c3 = c[0] + c[1] + c[2] + c[11] + c[12] + c[13] #wood
    c4 = c[3] + c[4] + c[5] + c[6] + c[7] +c[8] #steel and concrete
    c5 = c[15] #in the open considered to be outside
    #print 'c1',c1
    return np.array([c1,c2,c3,c4,c5])
#fileName = 'shelterTables.txt'
#getVals(fileName,0)


class shelter:
    def __init__(self,fileName=None,seasonAndTime=None):
        if seasonAndTime==None:
            print 'Season and time options:'
            print 'weekday-daytime-summer => 0'
            print 'weekday-daytime-winter => 1'
            print 'weekday-night          => 2'
            print 'weekend-daytime-summer => 3'
            print 'weekend-daytime-winter => 4'
            print 'weekend-night          => 5'
            exit()
        
        c = getVals(fileName,seasonAndTime)
        c_sum = np.sum(c)
        #print 'Value for sheltering fraction',c_sum,'This value will be used to normalize the sheltering fraction'
        c = c/c_sum
        c5array = getVals5Groups(c)

        if (np.sum(c5array) - 1.0)**2 >.000001:
            print 'Error: sheltering fraction does not add up to 1.'
            print 'Check sheltering.py'
            exit
        self.vals = c5array
        self.myInfo = 'self.vals = [light Metals,Composite,Wood,Concrete and steel,open]'
    def update(self,percentInOpen):
        if percentInOpen > 1.0:
            print 'Invalid input for shelter.update, value must be less than 1'
            print 'Percent in Open',percentInOpen
            print 'Adjusting to upper bound of 1' 
            percentInOpen = 1.0

        sumVals = np.sum(self.vals) - self.vals[4,0] # sum of percentages except for people in the open
        weight = (1.0 - percentInOpen)/sumVals
        self.vals = weight*self.vals
  

        self.vals[4,0] = percentInOpen
        if np.abs(np.sum(self.vals)-1.)>.000000001:
            print 'Error: New sheltering percentage does not equal 1'
            print 'Check inputs to sheltering.py shelter class'
            exit()
    def __doc__(self):
        print 'Season and time options:'
        print 'weekday-daytime-summer => 0'
        print 'weekday-daytime-winter => 1'
        print 'weekday-night          => 2'
        print 'weekend-daytime-summer => 3'
        print 'weekend-daytime-winter => 4'
        print 'weekend-night          => 5'



