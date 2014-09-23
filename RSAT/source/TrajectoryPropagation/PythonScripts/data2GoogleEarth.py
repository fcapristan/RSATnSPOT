

# File to get Latitude, Longitude and Altitude from debris.plt
# and then output it to Google Earth

def convert(mission):
    
    N=1
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
    line10GE = ' <name>Debris-Path</name>\n\n'  
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
        
        
        
        
        
        sample = 1
        counter = 0
        lonvals = mission.solution.longitude
        latvals = mission.solution.latitude
        height = mission.solution.height
        for indexW in range(len(lonvals)):
            newLine = str(lonvals[indexW])+','+str(latvals[indexW])+','+str(height[indexW])+'\n'
            outFileGE.write(newLine)
            altitude = height[indexW]
        
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
    
    outFileGE.close()




def convertPerStage(mission,fileNameGE):
    
    N = mission.vehicle.stages
    altitude =0
    #fileNameGE = 'GEperStage.kml'
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
    line10GE = ' <name>TrajPath</name>\n\n'  
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
        
        
        line1GE = '  <Placemark>\n'
        line2GE = '   <name>Stage'+str(index)+'</name>\n'
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
        
        
        
        
        
        sample = 1
        counter = 0
        lonvals = mission.solutionStage.longitude[index-1]
        latvals = mission.solutionStage.latitude[index-1]
        height = mission.solutionStage.height[index-1]
        for indexW in range(len(lonvals)):
            newLine = str(lonvals[indexW])+','+str(latvals[indexW])+','+str(height[indexW])+'\n'
            outFileGE.write(newLine)
            altitude = height[indexW]
        
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
    
    outFileGE.close()

def convertPerStagePropagate(stateList,fileNameGE):
    
    N = len(stateList)
    altitude =0
    #fileNameGE = 'GEperStage.kml'
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
    line10GE = ' <name>TrajPath</name>\n\n'  
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
        
        
        line1GE = '  <Placemark>\n'
        line2GE = '   <name>Stage'+str(index)+'</name>\n'
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
        
        
        
        
        
        sample = 1
        counter = 0
        lonvals = stateList[index-1][:,6]
        latvals = stateList[index-1][:,7]
        height = stateList[index-1][:,8]
        for indexW in range(len(lonvals)):
            newLine = str(lonvals[indexW])+','+str(latvals[indexW])+','+str(height[indexW])+'\n'
            outFileGE.write(newLine)
            altitude = height[indexW]
        
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
    
    outFileGE.close()

def convertPropagate(latvals,lonvals,height,fileNameGE):
    
    N = 1
    altitude =0
    #fileNameGE = 'GEperStage.kml'
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
    line10GE = ' <name>TrajPath</name>\n\n'  
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
        
        
        line1GE = '  <Placemark>\n'
        line2GE = '   <name>Stage'+str(index)+'</name>\n'
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
        
        
        
        
        
        sample = 1
        counter = 0

        for indexW in range(len(lonvals)):
            newLine = str(lonvals[indexW])+','+str(latvals[indexW])+','+str(height[indexW])+'\n'
            outFileGE.write(newLine)
            altitude = height[indexW]
        
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
    
    outFileGE.close()


def convertSimplePerStage(latList,lonList,heightList,fileNameGE):
    N = len(latList)
    altitude =0
    #fileNameGE = 'GEperStage.kml'
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
    line10GE = ' <name>TrajPath</name>\n\n'  
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
        
        
        line1GE = '  <Placemark>\n'
        line2GE = '   <name>Event'+str(index)+'</name>\n'
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
        
        
        
        
        
        sample = 1
        counter = 0
        lonvals = lonList[index-1]
        latvals = latList[index-1]
        height = heightList[index-1]
        for indexW in range(len(lonvals)):
            newLine = str(lonvals[indexW])+','+str(latvals[indexW])+','+str(height[indexW])+'\n'
            outFileGE.write(newLine)
            altitude = height[indexW]
        
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
    
    outFileGE.close()

