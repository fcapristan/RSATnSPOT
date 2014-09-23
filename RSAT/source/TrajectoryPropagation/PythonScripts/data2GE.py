
# File to get Latitude, Longitude and Altitude from debris.plt
# and then output it to Google Earth

def convert(baseFile,nSamples):

    N = nSamples

    altitude =0
    fileNameGE = 'GE.kml'
    outFileGE = open(fileNameGE,'w')

    line1GE  = '<?xml version="1.0" encoding="UTF-8"?>\n'
    line2GE  = '<kml xmlns="http://www.opengis.net/kml/2.2">\n'
    line3GE  = ' <Document>\n'
    line4GE  = '  <Style id="style1">\n'
    line5GE  = '   <LineStyle>\n'
    line6GE  = '    <colorMode>random</colorMode>\n'
    line7GE  = '    <width>4</width>\n'
    line8GE  = '   </LineStyle>\n'
    line9GE  = '  </Style>\n'
    line10GE = ' <name>'+baseFile+'</name>\n\n'  
    outFileGE.write(line1GE)
    outFileGE.write(line2GE)
    outFileGE.write(line3GE)
    outFileGE.write(line4GE)
    outFileGE.write(line5GE)
    outFileGE.write(line6GE)
    outFileGE.write(line7GE)
    outFileGE.write(line8GE)
    outFileGE.write(line9GE)
    outFileGE.write(line10GE)



    for index in range(1,N+1):

       fileName = baseFile+str(index)
       print fileName
       try:
          inputFile = open(fileName,'r')
    #   outputFile = open(fileNameOut, 'w')
       except:
          print 'Cant open file'

       line1GE = '  <Placemark>\n'
       line2GE = '   <name>Piece'+str(index)+'</name>\n'
       line3GE = '   <styleUrl>#style1</styleUrl>\n'
       line4GE = '   <MultiGeometry>\n'
       line5GE = '    <LineString>\n'
       line6GE = '      <altitudeMode>absolute</altitudeMode>\n'
       line7GE = '      <coordinates>\n'

       outFileGE.write(line1GE)
       outFileGE.write(line2GE)
       outFileGE.write(line3GE)
       outFileGE.write(line4GE)
       outFileGE.write(line5GE)
       outFileGE.write(line6GE)
       outFileGE.write(line7GE)





       lineNumber = 0
       sample = 1
       counter = 0
       altitude = -500

       for line in inputFile:
          counter = counter+1
          lineNumber = lineNumber + 1
          key = line.split() 
        
          if lineNumber > 2 :
             alt = float(key[2])
         #    print alt,index
             if abs(altitude - alt)>2: 
                 # if counter%sample==0:
                newLine = key[0]+','+key[1]+','+key[2]+'\n'
                outFileGE.write(newLine)
                altitude = float(key[2])

       line1GE = '    </coordinates>\n'
       line2GE = '   </LineString>\n'
       line3GE = '  </MultiGeometry>\n'
       line4GE = ' </Placemark>\n\n'
      
       outFileGE.write(line1GE)
       outFileGE.write(line2GE)
       outFileGE.write(line3GE)
       outFileGE.write(line4GE)


    line1GE = '  </Document>\n'
    line2GE = '</kml>'


    outFileGE.write(line1GE)
    outFileGE.write(line2GE)

    inputFile.close()
    outFileGE.close()
