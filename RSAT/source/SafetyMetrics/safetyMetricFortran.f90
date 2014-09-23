
  subroutine calculateEc(f,area,frows,fcols,Aproj,rp,popDensity,Ec)
! simple rectangular integration to obtain Prob from a bivariate pdf
! popDensity in Npeople per km2
    implicit none
    double precision , parameter :: Pi = 2*ACOS(0.0)
    integer , intent(in) :: frows,fcols
    double precision , dimension(frows,fcols) , intent(in):: f,popDensity
    double precision , intent(in) :: Aproj , rp , area
    double precision :: Ac,cdf
    !double precision , dimension(frows,fcols) :: tempCalculation
    double precision ,intent(out):: Ec
    integer :: i,j
   ! tempCalculation = f*popDensity
    Ec = 0
    cdf = 0
    do i=1,frows
       do j = 1,fcols
          cdf = cdf + f(i,j)
          Ec = Ec + f(i,j)*popDensity(i,j)
       end do
    end do

    !Ec = sum(tempCalculation)
    cdf = cdf* area
    print *, 'Cdf ',cdf
    Ac = Pi*((Aproj/Pi)**.5D0 + rp )**2.0D0
    Ec = Ac*Ec*area*(.001D0)**2.0D0 ! .001 **2 added to convert km2 to m2
  end subroutine calculateEc

!subroutine getCasualtySinglePiece(keyPop,keyArea,keycols,keyrows,xllcorner,yllcorner,&
!     cellsize,xMax,yMax,xyLocations,EcVals,ignoreBounds,rp,Aproj,samples)
!     implicit none
!    double precision , parameter :: Pi = 2*ACOS(0.0)

!  integer , intent(in) :: keycols,keyrows,ignoreBounds,samples
!  double precision , dimension(keyrows,keycols) , intent(in) :: keyPop,keyArea
!  double precision , dimension(samples,2) , intent(in) :: xyLocations ! essentially lattitude and longitude locations

!  double precision , intent(in) :: xllcorner, yllcorner, cellsize , xMax , yMax
!  double precision  , dimension(samples),intent(out) :: casualties
!  double precision ,dimension(samples) :: popDensity
!  double precision , dimension(samples) ,intent(in) :: Aproj
!  double precision  :: Ac,km2_to_m2
!  double precision ,intent(in) :: rp
  
!  km2_to_m2 = (0.001D0)**(2.0D0)
!  !call agsGetDensitySingle(keyPop,keyArea,keycols,keyrows,xllcorner,yllcorner,&
!  !   cellsize,xMax,yMax,xyLocations,popDensity,ignoreBounds,samples)

!  call agsGetDensitySingleCont(keyPop,keyArea,keycols,keyrows,xllcorner,yllcorner,&
!     cellsize,xMax,yMax,xyLocations,popDensity,ignoreBounds,samples)
!  Ac = Pi*((Aproj/Pi)**.5D0 + rp)**2.0D0
!  casualties = popDensity*(Ac*km2_to_m2)
  
!end subroutine getCasualtySinglePiece


subroutine calculateEcSinglePiece(keyPop,keyArea,keycols,keyrows,xllcorner,yllcorner,&
     cellsize,xMax,yMax,xyLocations,EcVals,ignoreBounds,rp,Aproj,samples)
     implicit none
    double precision , parameter :: Pi = 2*ACOS(0.0)

  integer , intent(in) :: keycols,keyrows,ignoreBounds,samples
  double precision , dimension(keyrows,keycols) , intent(in) :: keyPop,keyArea
  double precision , dimension(samples,2) , intent(in) :: xyLocations ! essentially lattitude and longitude locations

  double precision , intent(in) :: xllcorner, yllcorner, cellsize , xMax , yMax
  double precision  , dimension(samples),intent(out) :: EcVals
  double precision ,dimension(samples) :: popDensity
  double precision  :: Aproj,AC
  double precision ,intent(in) :: rp
  

  !call agsGetDensitySingle(keyPop,keyArea,keycols,keyrows,xllcorner,yllcorner,&
  !   cellsize,xMax,yMax,xyLocations,popDensity,ignoreBounds,samples)

  call agsGetDensitySingleCont(keyPop,keyArea,keycols,keyrows,xllcorner,yllcorner,&
     cellsize,xMax,yMax,xyLocations,popDensity,ignoreBounds,samples)
  Ac = Pi*((Aproj/Pi)**.5D0 + rp)**2.0D0
    !print *,'PopDen', popDensity
    !print *,'Ac',Ac,Aproj,rp
  EcVals = popDensity*(Ac*(.001D0)**2.0D0)!popdensity in km2 Ac in m2
end subroutine calculateEcSinglePiece



  subroutine calculateEcMatrix(f,area,frows,fcols,Aproj,rp,popDensity,Ec,EcMatrix)
! simple rectangular integration to obtain Prob from a bivariate pdf
! popDensity in Npeople per km2
    implicit none
    double precision , parameter :: Pi = 2*ACOS(0.0)
    integer , intent(in) :: frows,fcols
    double precision , dimension(frows,fcols) , intent(in):: f,popDensity
    double precision , intent(in) :: Aproj , rp , area
    double precision :: Ac,cdf,tempval
    !double precision , dimension(frows,fcols) :: tempCalculation
    double precision , dimension(frows,fcols), intent(out) :: EcMatrix
    double precision ,intent(out):: Ec
    integer :: i,j
   ! tempCalculation = f*popDensity
    Ec = 0.0
    cdf = 0.0
    Ac = Pi*((Aproj/Pi)**.5 + rp )**2

    do i=1,frows
       do j = 1,fcols
          cdf = cdf + f(i,j)
          !Ec = Ec + f(i,j)*popDensity(i,j)
          tempval= f(i,j)*popDensity(i,j)*Ac*area*(.001)**2! .001 **2 added to convert km2 to m2
          Ec = Ec + tempval
          EcMatrix(i,j) = tempval
       end do
    end do

    !Ec = sum(tempCalculation)
    !cdf = cdf* area
    !print *, 'Cdf ',cdf
    !Ec = Ac*Ec*area*(.001)**2 ! .001 **2 added to convert km2 to m2
  end subroutine calculateEcMatrix






    subroutine calculateEcBlast(popDensity,area,Ec)
      implicit none
      double precision , intent(in) :: popDensity, area
      double precision , intent(out) :: Ec
      
      Ec = popDensity*area*(.001)**2
    end subroutine calculateEcBlast

    subroutine calculateEcShelteringpopDensity(f,area,frows,fcols,nRoofs,CasualtyArea,roofFraction,rp,popDensity,Ec)
! simple rectangular integration to obtain Prob from a bivariate pdf
! this version takes into account sheltering information
    implicit none
    double precision , parameter :: Pi = 2*ACOS(0.0)
    integer , intent(in) :: frows,fcols,nRoofs
    double precision , dimension(frows,fcols) , intent(in):: f,popDensity
    double precision , dimension(nRoofs) :: CasualtyArea,roofFraction
    double precision , intent(in) :: rp , area
    double precision :: Ac,intSum,popKdensity
    !double precision , dimension(frows,fcols) :: tempCalculation
    double precision ,intent(out):: Ec
    integer :: i,j,k
    Ec = 0

    do i=1,frows
       do j = 1,fcols
         ! cdf = cdf + f(i,j) 
          intSum = 0.0
          if (popDensity(i,j)>0.0) then
             do k = 1,(nRoofs-1) ! only to n-1, because last shelter category is considered to be in the open, thus other calculations is required
                popKdensity = roofFraction(k)*popDensity(i,j)
                intSum = intSum + popKdensity*CasualtyArea(k)
             end do
             Ac = Pi*((CasualtyArea(nRoofs)/Pi)**.5 + rp )**2.0 ! Casualty Area for calculation in the open
             intSum = intSum + Ac*popDensity(i,j)*roofFraction(nRoofs)
             Ec = Ec + f(i,j)*area*intSum*(.001D0)**2.0   ! .001 **2 added to convert km2 to m2
       
          end if
       end do
    end do




  end subroutine calculateEcShelteringpopDensity

    subroutine calculateEcMatrixShelteringpopDensity(f,area,frows,fcols,nRoofs,CasualtyArea,&
      roofFraction,rp,popDensity,Ec,EcMatrix)
! simple rectangular integration to obtain Prob from a bivariate pdf
! this version takes into account sheltering information
! population contains the number of people inside popArea. popArea is in units of degrees ^2 (lat lon map)
    implicit none
    double precision , parameter :: Pi = 2*ACOS(0.0)
    integer , intent(in) :: frows,fcols,nRoofs
    double precision , dimension(frows,fcols) , intent(in):: f,popDensity
    double precision , dimension(nRoofs) :: CasualtyArea,roofFraction
    double precision , intent(in) :: rp , area
    double precision :: Ac,intSum,popKdensity
    !double precision , dimension(frows,fcols) :: tempCalculation
    double precision , dimension(frows,fcols),intent(out) :: EcMatrix
    double precision ,intent(out):: Ec
    double precision :: EcVal,cdf
    integer :: i,j,k
    Ec = 0.0
    cdf = 0.0
    do i=1,frows
       do j = 1,fcols
          ! cdf = cdf + f(i,j) 
          intSum = 0.0
          cdf = cdf + f(i,j)

          if (popDensity(i,j)>0.0) then
             do k = 1,(nRoofs-1) ! only to n-1, because last shelter category is considered to be in the open, thus other calculations is required
                popKdensity = roofFraction(k)*popDensity(i,j)
                intSum = intSum + popKdensity*CasualtyArea(k)
             end do
             Ac = Pi*((CasualtyArea(nRoofs)/Pi)**.5 + rp )**2 ! Casualty Area for calculation in the open
             intSum = intSum + Ac*popDensity(i,j)*roofFraction(nRoofs)
             EcVal = f(i,j)*area*intSum*(.001D0)**2   ! .001 **2 added to convert km2 to m2
             Ec = Ec + EcVal
             EcMatrix(i,j) = Ecval 
          else
             EcMatrix(i,j) = 0.0
          end if
       end do
    end do
    print *,'cdf',cdf*area


  end subroutine calculateEcMatrixShelteringpopDensity




    subroutine calculateEcSheltering(f,area,frows,fcols,nRoofs,CasualtyArea,roofFraction,rp,population,popArea,Ec)
! simple rectangular integration to obtain Prob from a bivariate pdf
! this version takes into account sheltering information
! population contains the number of people inside popArea. popArea is in units of degrees ^2 (lat lon map)
    implicit none
    double precision , parameter :: Pi = 2*ACOS(0.0)
    integer , intent(in) :: frows,fcols,nRoofs
    double precision , dimension(frows,fcols) , intent(in):: f,population
    double precision , dimension(nRoofs) :: CasualtyArea,roofFraction
    double precision , intent(in) :: rp , area,popArea
    double precision :: Ac,Npeople,intSum,populationFixed
    !double precision , dimension(frows,fcols) :: tempCalculation
    double precision ,intent(out):: Ec
    integer :: i,j,k
    Ec = 0
    !cdf = 0
    print *,'AreaRatio', area/popArea
    print *, 'CasArea',CasualtyArea
    print *, 'RoofFraction',roofFraction
    do i=1,frows
       do j = 1,fcols
         ! cdf = cdf + f(i,j) 
          intSum = 0.0
          populationFixed = population(i,j)*area/popArea ! area/popArea returns the correct number of people inside the area used for integration
          do k = 1,(nRoofs-1) ! only to n-1, because last shelter category is considered to be in the open, thus other calculations is required
             Npeople = populationFixed*roofFraction(k) ! number of people in the shelter group inside the area used for integration (area).
             intSum =  intSum + Npeople*CasualtyArea(k)
          end do
          intSum = 0.0
          Ac = Pi*((CasualtyArea(nRoofs)/Pi)**.5 + rp )**2 ! Casualty Area for calculation in the open
          intSum = intSum + Ac*populationFixed*roofFraction(nRoofs)

          if (intSum>populationFixed) then
             intSum = populationFixed
           end if
           Ec = Ec + f(i,j)*intSum
       end do
    end do


  end subroutine calculateEcSheltering

subroutine checkPolygon(Xmat,Ymat,vertx,verty,MatOut,rows,cols,lenvert)
!checks if a given Mesh has values that are inside a polygon
implicit none
integer , intent(in) :: rows,cols,lenvert
double precision, dimension(rows,cols), intent(in) :: Xmat,Ymat
double precision , dimension(lenvert), intent(in):: vertx,verty
integer, dimension(rows,cols) ,intent(out)::MatOut
integer i,j,l,m
double precision :: testx,testy


do i=1,rows
   do j=1,cols
      testx = Xmat(i,j)
      testy = Ymat(i,j)
      call locpt (testx, testy, vertx, verty, lenvert, l, m)
      if (l.GE.0) then
         MatOut(i,j) = 1
      else
         MatOut(i,j)= 0
      end if
   end do
end do


end subroutine checkPolygon



subroutine updateMatPolygon(Xmat,Ymat,vertx,verty,MatIn,MatOut,rows,cols,lenvert,desVal)
!UPDATES VALUES IN A MATRIX WITH DESVAL IF VALUES IN THE CORRESPONDING MESH ARE INSIDE THE DESIRED POLYGON
implicit none
integer , intent(in) :: rows,cols,lenvert
double precision, dimension(rows,cols), intent(in) :: Xmat,Ymat
double precision , dimension(lenvert), intent(in):: vertx,verty
double precision , intent(in) :: desVal
double precision, dimension(rows,cols) ,intent(in)::MatIn
double precision, dimension(rows,cols) ,intent(out)::MatOut
integer i,j,l,m
double precision :: testx,testy


do i=1,rows
   do j=1,cols
      testx = Xmat(i,j)
      testy = Ymat(i,j)
      call locpt (testx, testy, vertx, verty, lenvert, l, m)
      if (l.GE.0) then
         MatOut(i,j) = desVal
      else
         MatOut(i,j)= MatIn(i,j)
      end if
   end do
end do


end subroutine updateMatPolygon


subroutine agsGetVals(key,keycols,keyrows,xllcorner,yllcorner,cellsize,Xmat,Ymat,XYcols,XYrows,xMax,yMax,popMatrix,ignoreBounds)
  implicit none

  integer , intent(in) :: keycols,keyrows,XYcols,XYrows,ignoreBounds
  double precision , dimension(keyrows,keycols) , intent(in) :: key
  double precision , dimension(XYrows,XYcols) , intent(in) :: Xmat,Ymat
  double precision , intent(in) :: xllcorner, yllcorner, cellsize , xMax , yMax
  double precision , dimension(XYrows,XYcols) , intent(out) :: popMatrix
  integer :: i,j , nRow, nColumn
  double precision :: x,y
  do i = 1,XYrows
     do j = 1,XYcols
        x = Xmat(i,j)
        y = Ymat(i,j)


        if (x>180.) then !making sure there are some upper bounds on longitude

           x = x - 360.
          
        elseif (x<-180.) then
           x = x + 360.
        end if

        if (y > 90.) then
           y = 180. - y
        elseif (y < -90.) then
           y = -180 - y
        end if


        
        nRow = keyrows - floor((y - yllcorner)/cellsize)
        nColumn = 1 + floor((x-xllcorner)/cellsize)




        if (x==xMax) then
           nColumn = nColumn - 1
        end if
       

        if ((x>xMax).or.(y>yMax).or.(x<xllcorner).or.(y<yllcorner)) then

           if (ignoreBounds.NE.1) then


              print *, 'Error in agsGetVals.f90, check X Y bounds in ASC file \n'
              print *,'Xmax is ',Xmax,'. X given is ',x
              print *,'Ymax is ',Ymax,'. X given is ',y

              print *,'Xmin is ',xllcorner,'. X given is ',x
              print *,'Ymin is ',yllcorner,'. X given is ',y

              STOP
           else
              popMatrix(i,j) = 0.0 !ignoring out of bounds values and setting them to zero
             ! print *, 'Ignoring out of bounds val'

           end if
        else 
          ! print *,x,xMax,y,yMax

           if (key(nRow,nColumn)<0) then
              popMatrix(i,j) = 0.0 !ignoring negative values for population
              !print *, 'setting neg vals to zero'
           else
              
              popMatrix(i,j) = key(nRow,nColumn)


           end if
        end if

     end do
  end do

end subroutine agsGetVals





subroutine updatePopulation(keyPop,EcMatrix,xMatLocations,yMatLocations,keyRows,keyCols,EcRows,EcCols)
   implicit none

  integer,intent(in) :: keyRows, keyCols, EcRows, EcCols

  double precision, dimension(keyRows,keyCols),intent(inout) :: keyPop
  double precision, dimension(EcRows,EcCols), intent(in) :: EcMatrix
  integer, dimension(EcRows,EcCols), intent(in) :: xMatLocations,yMatLocations
  integer :: i,j

  do i=1,EcRows
     do j=1,EcCols
        if (EcMatrix(i,j)>0.0) then
           ! in theory xMatLocations will be valid as long as Ec >0. When out of bounds Ec = 0, so no update is necessary 
           keyPop(xMatLocations(i,j),yMatLocations(i,j)) = max(keyPop(xMatLocations(i,j),yMatLocations(i,j)) - EcMatrix(i,j),0.0)
        end if
     end do
  end do
end subroutine updatePopulation



subroutine pnpoly(nvert,vertx,verty,testx,testy,cval)
  implicit none
  integer, intent(in):: nvert
  integer, intent(out):: cval
  double precision , intent(in) :: testx,testy
  double precision , dimension(nvert),intent(in)::vertx,verty
  integer ::i,j
  cval = 0
  do i=1,nvert
     if (i.eq.1) then
        j = nvert
     else
        j = i-1
     end if
     if (((verty(i) > testy).NEQV.(verty(j)>testy)).AND.(testx <((vertx(j) &
        -vertx(i))*(testy-verty(i))/ &
         (verty(j)-verty(i)) + vertx(i)))) then
     
        cval = 1
        exit
     end if
  end do
end subroutine pnpoly


subroutine agsGetDensity(keyPop,keyArea,keycols,keyrows, &
xllcorner,yllcorner,cellsize,Xmat,Ymat,XYcols,XYrows,xMax,&
yMax,popDensity,XMatLocations,YMatLocations,ignoreBounds)
  implicit none

  integer , intent(in) :: keycols,keyrows,XYcols,XYrows,ignoreBounds
  double precision , dimension(keyrows,keycols) , intent(in) :: keyPop,keyArea
  double precision , dimension(XYrows,XYcols) , intent(in) :: Xmat,Ymat
  double precision , intent(in) :: xllcorner, yllcorner, cellsize , xMax , yMax
  double precision , dimension(XYrows,XYcols) , intent(out) :: popDensity
  integer , dimension(XYrows,XYcols) , intent(out) :: XMatLocations,YMatLocations
  integer :: i,j , nRow, nColumn
  double precision :: x,y
  do i = 1,XYrows
     do j = 1,XYcols
        x = Xmat(i,j)
        y = Ymat(i,j)


        if (x>180.) then !making sure there are some upper bounds on longitude

           x = x - 360.
          
        elseif (x<-180.) then
           x = x + 360.
        end if

        if (y > 90.) then
           y = 180. - y
        elseif (y < -90.) then
           y = -180 - y
        end if


        
        nRow = keyrows - floor((y - yllcorner)/cellsize)
        nColumn = 1 + floor((x-xllcorner)/cellsize)




        if (x==xMax) then
           nColumn = nColumn - 1
        end if
       

        if ((x>xMax).or.(y>yMax).or.(x<xllcorner).or.(y<yllcorner)) then

           if (ignoreBounds.NE.1) then

              
              print *, 'Error in agsGetVals.f90, check X Y bounds in ASC file \n'
              print *,'Xmax is ',Xmax,'. X given is ',x
              print *,'Ymax is ',Ymax,'. X given is ',y

              print *,'Xmin is ',xllcorner,'. X given is ',x
              print *,'Ymin is ',yllcorner,'. X given is ',y
              
              STOP
           else
              popDensity(i,j) = 0.0 !ignoring out of bounds values and setting them to zero
             ! print *, 'Ignoring out of bounds val'
              xMatLocations(i,j) = -999 ! locations out of bounds
              yMatLocations(i,j) = -999
           end if
        else 
          ! print *,x,xMax,y,yMax
              xMatLocations(i,j) = nRow-1 ! -1 becuase it needs to agree with python indexing
              yMatLocations(i,j) = nColumn-1
           if ((keyPop(nRow,nColumn)<=0).or.(keyArea(nRow,nColumn)<=0)) then
              popDensity(i,j) = 0.0 !ignoring negative values for population
              !print *, 'setting neg vals to zero'
           else
              
              popDensity(i,j) = keyPop(nRow,nColumn)/keyArea(nRow,nColumn)


           end if
        end if




     end do
  end do

end subroutine agsGetDensity





subroutine agsGetDensitySingle(keyPop,keyArea,keycols,keyrows,xllcorner,yllcorner,&
     cellsize,xMax,yMax,xyLocations,popDensity,ignoreBounds,samples)
  implicit none

  integer , intent(in) :: keycols,keyrows,ignoreBounds,samples
  double precision , dimension(keyrows,keycols) , intent(in) :: keyPop,keyArea
  double precision , dimension(samples,2) , intent(in) :: xyLocations ! essentially lattitude and longitude locations

  double precision , intent(in) :: xllcorner, yllcorner, cellsize , xMax , yMax
  double precision , dimension(samples),  intent(out) :: popDensity
  integer :: nRow, nColumn,iloc,jloc,index
  double precision :: x,y,popdensitytemp
  double precision :: x0,y0,x1,y1,s,popX,popY
  

  do index = 1,samples


     x = xyLocations(index,1)
     y = xyLocations(index,2)


        if (x>180.) then !making sure there are some upper bounds on longitude

           x = x - 360.
          
        elseif (x<-180.) then
           x = x + 360.
        end if

        if (y > 90.) then
           y = 180. - y
        elseif (y < -90.) then
           y = -180 - y
        end if


        
        nRow = keyrows - floor((y - yllcorner)/cellsize)
        nColumn = 1 + floor((x-xllcorner)/cellsize)




        if (x==xMax) then
           nColumn = nColumn - 1
        end if
       

        if ((x>xMax).or.(y>yMax).or.(x<xllcorner).or.(y<yllcorner)) then

           if (ignoreBounds.NE.1) then

              
              print *, 'Error in agsGetVals.f90, check X Y bounds in ASC file \n'
              print *,'Xmax is ',Xmax,'. X given is ',x
              print *,'Ymax is ',Ymax,'. X given is ',y

              print *,'Xmin is ',xllcorner,'. X given is ',x
              print *,'Ymin is ',yllcorner,'. X given is ',y
              
              STOP
           else
              popDensity(index) = 0.0 !ignoring out of bounds values and setting them to zero

           end if
        else 

           if ((keyPop(nRow,nColumn)<=0).or.(keyArea(nRow,nColumn)<=0)) then
              popDensity(index) = 0.0 !ignoring negative values for population
              !print *, 'setting neg vals to zero'
           else
              
              popDensity(index) = keyPop(nRow,nColumn)/keyArea(nRow,nColumn)


           end if
        end if

     end do



end subroutine agsGetDensitySingle


subroutine agsGetDensitySingleCont(keyPop,keyArea,keycols,keyrows,xllcorner,yllcorner,&
     cellsize,xMax,yMax,xyLocations,popDensity,ignoreBounds,samples)
  implicit none

  integer , intent(in) :: keycols,keyrows,ignoreBounds,samples
  double precision , dimension(keyrows,keycols) , intent(in) :: keyPop,keyArea
  double precision , dimension(samples,2) , intent(in) :: xyLocations ! essentially lattitude and longitude locations

  double precision , intent(in) :: xllcorner, yllcorner, cellsize , xMax , yMax
  double precision , dimension(samples),  intent(out) :: popDensity
  integer :: nRow, nColumn,iloc,jloc,index
  double precision :: x,y,popdensitytemp
  double precision :: x0,y0,x1,y1,popX,popY
  double precision :: fR1,fR2,fQ11,fQ21,fQ22,fQ12
  double precision :: px1,px2,py1,py2

  do index = 1,samples


     x = xyLocations(index,1)
     y = xyLocations(index,2)

     
        if (x>180.) then !making sure there are some upper bounds on longitude

           x = x - 360.
          
        elseif (x<-180.) then
           x = x + 360.
        end if

        if (y > 90.) then
           y = 180. - y
        elseif (y < -90.) then
           y = -180 - y
        end if


        
        !nRow = keyrows - floor((y - yllcorner)/cellsize)
        nRow = 1 + floor((yMax - y)/cellsize)

        nColumn = 1 + floor((x-xllcorner)/cellsize)

        x0 = (nColumn-1)*cellsize+xllcorner + .5*cellsize
        y0 = yMax - (nRow-1)*cellsize - .5*cellsize
        x1 = x-x0
        y1 = y-y0

        if (x==xMax) then
           nColumn = nColumn - 1
        end if
       

        if ((x>xMax).or.(y>yMax).or.(x<xllcorner).or.(y<yllcorner)) then

           if (ignoreBounds.NE.1) then

              
              print *, 'Error in agsGetVals.f90, check X Y bounds in ASC file \n'
              print *,'Xmax is ',Xmax,'. X given is ',x
              print *,'Ymax is ',Ymax,'. X given is ',y

              print *,'Xmin is ',xllcorner,'. X given is ',x
              print *,'Ymin is ',yllcorner,'. X given is ',y
              
              STOP
           else
              popDensity(index) = 0.0 !ignoring out of bounds values and setting them to zero

           end if
        else 
          


           if ((keyPop(nRow,nColumn)<=0).or.(keyArea(nRow,nColumn)<=0)) then
              popDensity(index) = 0.0 !ignoring negative values for population
              !print *, 'setting neg vals to zero'
           else

              if (x1.GE.0) then ! check the ordering
                 iloc = nColumn+1
                 px2 = x0 + cellsize
              else
                 iloc = nColumn-1
                 px2 = x0 - cellsize
              end if

              if (y1.GE.0) then ! check the ordering
                 jloc = nRow-1
                 py2 = y0 + cellsize
              else
                 jloc = nRow+1
                 py2 = y0 - cellsize
              end if

              if (jloc.LE.0) then
                 jloc = 1
              elseif (jloc.GT.keyRows) then
                 jloc = keyRows
              endif


              if (iloc.LE.0) then
                 iloc = 1
              elseif (iloc.GT.keyCols) then
                 iloc = keyCols
              endif

              if ((abs(x1).GT..5*cellsize).or.(abs(y1).GT..5*cellsize)) then

                 print *, 'x1 or y1 > cellsize. Error in safetyMetricFortran.f90'
                 print *,'x1 ',x1,'y1 ',y1,'cellsize ',cellsize
                 print *,'nColumn',nColumn,nRow,xllcorner,yllcorner,x,y
                 STOP 'check x1 or y1'
              end if
              ! bilinear interpolation suggested in http://en.wikipedia.org/wiki/Bilinear_interpolation
              px1 = x0
              py1 = y0
              fQ11 = keyPop(nRow,nColumn)/keyArea(nRow,nColumn)
              fQ22 = keyPop(jloc,iloc)/keyArea(jloc,iloc)
              fQ12 = keyPop(jloc,nColumn)/keyArea(jloc,nColumn)
              fQ21 = keyPop(nRow,iloc)/keyArea(nRow,iloc)
              
              if ((keyArea(nRow,nColumn).LE.0).or.(fQ11.LT.0))  then
                 fQ11 = 0.0
              end if
              if ((keyArea(jloc,iloc).LE.0).or.(fQ22.LT.0)) then
                 fQ22 = 0.0
              end if
              if ((keyArea(jloc,nColumn).LE.0).or.(fQ12.LT.0)) then
                 fQ12 = 0.0
              end if
              if ((keyArea(nRow,iloc).LE.0).or.(fQ21.LT.0)) then
                 fQ21 = 0.0
              end if
              fR1 = (px2-x)/(px2-px1)*fQ11 + (x-px1)/(px2-px1)*fQ21
              fR2 = (px2-x)/(px2-px1)*fQ12 + (x-px1)/(px2-px1)*fQ22
              popDensity(index)= (py2-y)/(py2-py1)*fR1 + (y-py1)/(py2-py1)*fR2

  
           end if
        end if
     end do


      end subroutine agsGetDensitySingleCont


SUBROUTINE locpt (x0, y0, x, y, n, l, m)
!code obtained from http://jblevins.org/mirror/amiller/locpt.f90
!All code written by Alan Miller is released into the public domain
!Slightly modified by Francisco Capristan
!-----------------------------------------------------------------------
! GIVEN A POLYGONAL LINE CONNECTING THE VERTICES (X(I),Y(I)) (I = 1,...,N)
! TAKEN IN THIS ORDER.  IT IS ASSUMED THAT THE POLYGONAL PATH IS A LOOP,
! WHERE (X(N),Y(N)) = (X(1),Y(1)) OR THERE IS AN ARC FROM (X(N),Y(N)) TO
! (X(1),Y(1)).  N.B. The polygon may cross itself any number of times.

! (X0,Y0) IS AN ARBITRARY POINT AND L AND M ARE VARIABLES.
! On output, L AND M ARE ASSIGNED THE FOLLOWING VALUES ...

!    L = -1   IF (X0,Y0) IS OUTSIDE THE POLYGONAL PATH
!    L =  0   IF (X0,Y0) LIES ON THE POLYGONAL PATH
!    L =  1   IF (X0,Y0) IS INSIDE THE POLYGONAL PATH

! M = 0 IF (X0,Y0) IS ON OR OUTSIDE THE PATH.  IF (X0,Y0) IS INSIDE THE
! PATH THEN M IS THE WINDING NUMBER OF THE PATH AROUND THE POINT (X0,Y0).

! Fortran 66 version by A.H. Morris
! Converted to ELF90 compatibility by Alan Miller, 15 February 1997

!-----------------------

IMPLICIT NONE
double precision, INTENT(IN)     :: x0, y0
INTEGER, INTENT(IN)  :: n
double precision , dimension(n),intent(in)::x,y
INTEGER, INTENT(OUT) :: l,m


!     Local variables
INTEGER :: i, n0
double precision    :: angle, eps, pi, pi2, sum, theta, theta1, thetai, tol, u, v

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1.0 + EPS > 1.0

eps = EPSILON(1.0)

!-----------------------------------------------------------------------
! NOTE THAT m is an integer but then real operations are applied to it.
! After some testing this gives the right behavior (m as integer). I tried as a double
! and the code gave false positives (false inside polygon). This is probably because m needs
! to be rounded. To keep changes to a minimum I decided to keep m as an integer. Note that a 
! compiler warning may happen because of the real operations to determine m
n0 = n
IF (x(1) == x(n) .AND. y(1) == y(n)) n0 = n - 1
pi = ATAN2(0.0, -1.0)
pi2 = 2.0*pi
tol = 4.0*eps*pi
l = -1
m = 0

u = x(1) - x0
v = y(1) - y0
IF (u == 0.0 .AND. v == 0.0) GO TO 20
IF (n0 < 2) RETURN
theta1 = ATAN2(v, u)

sum = 0.0
theta = theta1
DO i = 2, n0
  u = x(i) - x0
  v = y(i) - y0
  IF (u == 0.0 .AND. v == 0.0) GO TO 20
  thetai = ATAN2(v, u)
  
  angle = ABS(thetai - theta)
  IF (ABS(angle - pi) < tol) GO TO 20
  IF (angle > pi) angle = angle - pi2
  IF (theta > thetai) angle = -angle
  sum = sum + angle
  theta = thetai
END DO

angle = ABS(theta1 - theta)
IF (ABS(angle - pi) < tol) GO TO 20
IF (angle > pi) angle = angle - pi2
IF (theta > theta1) angle = -angle
sum = sum + angle

!     SUM = 2*PI*M WHERE M IS THE WINDING NUMBER

m = ABS(sum)/pi2 + 0.2
IF (m == 0) RETURN
l = 1
IF (sum < 0.0) m = -m
RETURN

!     (X0, Y0) IS ON THE BOUNDARY OF THE PATH

20 l = 0
RETURN
END SUBROUTINE locpt
