subroutine fixLON4PDF(lonlat,nSamples,lonlatFIXED)
! Fix for lonlat array to ensure that there is continuity in the groundtrack
! e.g 180 and - 180 should be the same value

  implicit none
  integer , intent(in) :: nSamples
  double precision , dimension(nSamples,2) ,intent(in) :: lonlat
  double precision , dimension(nSamples,2) , intent(out) :: lonlatFIXED
  double precision , dimension(nSamples) :: longitude , latitude
  double precision :: lonmax,lonmin,deltalon,deltadesired
  integer :: i
  
  deltadesired = 300.0 ! value to activate fix. If difference between max lon and min lon is greater than this val, then fix value to be in same local plane
                       ! mainly to be used for PDF calculation
  longitude = lonlat(:,1)
  latitude = lonlat(:,2)
  
  lonmin = minval(longitude)
  lonmax = maxval(longitude)
  deltalon  = lonmax - lonmin

  if (deltalon > deltadesired) then
     do i = 1,nSamples
        !if (longitude(i)<-150.-(deltadesired - 360.)) then
        if (longitude(i)<0) then
 
           lonlatFIXED(i,1) = longitude(i) + 360.
           lonlatFIXED(i,2) = latitude(i)
        else
           lonlatFIXED(i,1) = longitude(i)
           lonlatFIXED(i,2) = latitude(i)
        end if
     end do
  else
     lonlatFIXED = lonlat
  end if
end subroutine fixLON4PDF

subroutine boundLON(lonlat,nSamples,lonlatFIXED)
  ! ensure that longitude is between -180 and 180 degrees

  implicit none
  integer , intent(in) :: nSamples
  double precision , dimension(nSamples,2) ,intent(in) :: lonlat
  double precision , dimension(nSamples,2) , intent(out) :: lonlatFIXED
  double precision , dimension(nSamples) :: longitude , latitude
  double precision :: lonmax,lonmin,fixedlon
  integer :: i           

  longitude = lonlat(:,1)
  latitude = lonlat(:,2)
  
  lonmin = minval(longitude)
  lonmax = maxval(longitude)

  if ((lonmin<-180).or.(lonmax>180)) then
     do i=1,nSamples
        fixedlon = bounder(longitude(i))
        lonlatFIXED(i,1) = fixedlon
        lonlatFIXED(i,2) = latitude(i)
     end do
  else
     lonlatFIXED = lonlat
  end if
contains

 double precision function bounder(longitude)
    ! bounder to ensure values are within -180 and 180 degrees
    implicit none
    double precision , intent(in) :: longitude
    double precision :: n
    bounder = 0.0
    n = longitude/360.
    if (longitude>180) then
       bounder = longitude - floor(n)*360
    elseif (longitude<180) then
       bounder = longitude - ceiling(n)*360
    end if

    if (bounder>180) then
       bounder = bounder - 360
    elseif (bounder<-180) then
       bounder = bounder + 360
    end if
    
  end function bounder
end subroutine boundLON

subroutine rectangularGrid(lonVec,latVec,ZVec,nPoints,delta,lonCorner,latCorner,nlon,nlat,ZMesh)
  implicit none
  integer , intent(in) :: nPoints,nlon,nlat
  double precision , dimension(nPoints) , intent(in) :: lonVec,latVec,Zvec
  double precision , intent(in) :: delta
  !integer , dimension(nPoints) :: lonLocation,latLocation
  double precision ,  dimension(nlat,nlon) ,intent(out):: Zmesh
  double precision ,  dimension(nlat,nlon) :: nZmesh
  double precision , intent(in) :: lonCorner,latCorner
  !double precision ,  dimension(nlon) , intent(in) :: lonVecMesh
  !double precision ,  dimension(nlat) , intent(in) :: latVecMesh
  integer :: i,j,lonLocation,latLocation
  !double precision :: minlon,maxlon,minlat,maxlat, lonCorner,latCorner
  !minlon = minval(lonVec)-2*delta
  !maxlon = maxval(lonVec)+2*delta
  !minlat = minval(latVec)-2*delta
  !maxlat = maxval(latVec)+2*delta

  !lonCorner = minlon - .5*delta
  !latCorner = minlat - .5*delta
  !nlon = ceiling((maxlon-minlon)/delta)
  !nlat = ceiling((maxlat-minlat)/delta)
  !allocate(lonVecMesh(nlon))
  !allocate(latVecMesh(nlat))
  !allocate(lonMesh(nlat,nlon))
  !allocate(latMesh(nlat,nlon))
  !allocate(Zmesh(nlat,nlon))
  !allocate(nZmesh(nlat,nlon))
  
  !lonVecMesh(1) = minlon
  !do i=2,nlon
  !   lonVecMesh(i) = lonVecMesh(i-1) + delta
  !end do
  !latVecMesh(1) = minlat
  !do i=2,nlat
  !   latVecMesh(i) = latVecMesh(i-1) + delta
  !end do
  
  do i=1,nlat
     do j = 1,nlon
        !lonMesh(i,j) = lonVecMesh(j)
        !latMesh(i,j) = latVecMesh(i)
        Zmesh(i,j) = 0.0
        nZmesh(i,j) = 0.0
     end do
  end do
  


       ! nRow = keyrows - floor((y - yllcorner)/cellsize)
       ! nColumn = 1 + floor((x-xllcorner)/cellsize)

  do i=1,nPoints
     lonLocation = ceiling((lonVec(i) - lonCorner)/delta)
     latLocation = nlat - floor((latVec(i) - latCorner)/delta)
     nZmesh(latLocation,lonLocation) = nZmesh(latLocation,lonLocation) + 1
     !Zmesh(latLocation,lonLocation) = Zvec(i) + Zmesh(latLocation,lonLocation)
  !   if (nZmesh(latLocation,lonLocation)>1) then
        Zmesh(latLocation,lonLocation) = ((nZmesh(latLocation,lonLocation)-1)*Zmesh(latLocation,lonLocation)+Zvec(i))&
             /nZmesh(latLocation,lonLocation) ! taking the average inside the cell
   !  else
    !    Zmesh(latLocation,lonLocation) = Zvec(i)
    ! end if
  end do
end subroutine rectangularGrid
  
  
  
        
