
def calculatesingleEc(lonlat,nSamples,delta,nsigma,key,xllcorner,yllcorner,cellsize,xMax,yMax,boundsOption,ncols,nrows,nPieces,arefMean,pdfoption):
    # calculates Ec by calculating pdf in P frame, then get corresponding population density for P frame.
    sys.path.append("../SafetyMetrics")
    import safetyMetrics as SMF
    
    # this version uses the rotated pdf (actual lon lat coordinates)
    areapop = delta*delta
    tempval = .3048 #1ft to meters
    lonPMesh,latPMesh,ZZpdfPframe,lonOrMesh,latOrMesh,xlen,ylen,transformDetails  = getPDF2(lonlat,nSamples,delta,nsigma,pdfoption)
    popMatrix = SMF.agsgetvals(key,xllcorner,yllcorner,cellsize,lonOrMesh,latOrMesh,xMax,yMax,boundsOption,[ncols,nrows,ylen,xlen])
    Ec = nPieces*(SMF.calculateec(ZZpdfPframe,areapop,arefMean,tempval,popMatrix,[ylen,xlen]))
    return (Ec,lonPMesh,latPMesh,ZZpdfPframe)
