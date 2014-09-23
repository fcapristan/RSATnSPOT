# Python routine to read CD/CL profile from a text file
# Created by Francisco Capristan


def reader(filename):
    inputFile = open(filename,'r')
    # The format expected is Minf -> first column. Cd or Cl -> second column
    
    M = [] # Minf
    Cx = []# coefficient CL or CD 
    for line in inputFile:
        key = line.split()
        if len(key)>0:
            if (line[0]!='#'):# this way "#" can be used to comment an entire line 
                M.append(float(key[0]))
                Cx.append(float(key[1]))
    n = len(M) 
    inputFile.close()

    return (M,Cx,n)