
subroutine debrisPropagation(finalConditions,initialState,DebrisVel,&
     mass,Sref,&
     nCD,MinfCD,CD,& !Drag inputs
     CLOption,nCL,MinfCL,CL,& !Lift inputs
     LoverD,& !LoverD inputs
     atmosOption,altitudeList, densityList,UList,VList,WList,& !! atmospheric parameter inputs
     GEOptions,filename,nList,PlanetModel,dtInterval,thetag0)




  use planetEarthCharacteristics
  use constants
  !use newpos_e10_I
  implicit none

  ! Debris propagation routine
  ! Created by Francisco Capristan for the FAA COE-CST effort
  ! July 2011

  ! Modified by Francisco Capristan
  ! June 2012
  ! Notes : Added CD and CL as a function of Minf

  !#### Inputs#################################################################
  ! InitialState     : [x,y,z,Vx,Vy,Vz] : inputs in ECEF coordinate frame SI UNITS

  ! DebrisVel          : [Vx_debris,Vy_debris,Vz_debris] : debris velocities ECEF coordinate frame
  ! mass                : debris piece mass [kg]
  ! Sref                : debris Sref [m2]
  ! nCD                 : length of CD array. if nCd = 1 then Constant Cd
  ! MinfCD              : array containing Mach numbers for CD interpolation. Must be length nCD 
  ! CD                  : debris piece drag coefficient function of Minf. Must be length nCD
  ! CLOption            : if 1 then CL is used. if 0 then constant LoverD used instead
  ! nCL                 : length of CL array. if nCL = 1, then constant CL
  ! MinfCL              : array containing Mach numbers for CL interpolation. Must be length nCL 
  ! CL                  : debris piece lift coefficient as a function of Minf.
  ! LoverD              : used if ClOption =0. Specifies constant L over D
  ! atmosOption         : if 0 Exponential density used NO WIND. if 1 given density profile used NO WIND. If 2 given density and wind profile used. -1 Vacuum
  ! altitudeList        : altitude vector in decreasing order (to be used in atmospheric profile)
  ! densityList         : density profile corresponding to altitudeList
  ! U/V/W  list         : wind profile corresponding to altitudeList. Usually obtained from GRAM
  ! GEOptions           : if 1, then the trajectory will be writen to filename 
  ! filename            : filename to write trajectory solution suitable for google earth visualization. Only used if GEOptions =1
  ! nList               : number of data points in atmospheric profile    
  ! PlanetModel         : 0 => spherical , 1=> oblate. J2 effects ignored (not calculating orbits, Drag and Lift are the dominating forces)         
  ! dtInterval          : approx value to do rk45 stepping...will be used to set upper limit for dt
  !### Outputs:

  ! FinalConditions :
  !                     finalConditions = (/altitudeFinal,latitudeFinal,longitudeFinal,VrelMag,flag/) 

  !                      altitudeFinal       : Final altitude [m], if all successful, then it must be <0
  !                      latitudeFinal       : final latitude [degrees]
  !                      longitudeFinal      : final longitude [degrees]
  !                      VrelMag             : velocity (relative to rotating planet) magnitude
  !                      flag                : returns potential issues (0 -> OK RUN , 1 -> debris landing, but time step appears too big,
  !                                            2 -> debris appears to be in orbit, or more iterations are needed...If not suppose to be in orbit =>try increasing number of max RK4 iterations...or decreasing dt
  !                                            3-> debris is in orbit
  !#############################################################################


  ! Input Parameters
  !##############################################
  double precision , dimension(3) , intent(in) :: DebrisVel
  double precision , dimension(6) , intent(in) :: initialState
  double precision , intent(in) :: mass,Sref, dtInterval
  character*16 ,intent(in) :: filename
  ! Aerodynamic inputs
  integer , intent(in) :: nCd,nCL
  double precision , dimension(nCD),intent(in) :: MinfCD,CD
  double precision , dimension(nCL),intent(in) :: MinfCL,CL
  double precision , intent(in) :: LoverD,thetag0



  ! Input Options
  integer , intent(in) :: GEOptions,CLOption,atmosOption,nList,PlanetModel
  double precision , dimension(nList),intent(in)::altitudeList,densityList,UList,VList,WList
  !##############################################

  !Output Parameters
  !##############################################
  double precision ::   altitudeFinal, VrelMag
  double precision :: flag
  double precision , dimension(7),intent(out) :: finalConditions

  !##############################################



  !##############################################

  !For Calculations
  !##############################################
  double precision :: longitude,latitude

  double precision :: maxLat,minLat,maxLon,minLon

  double precision, dimension(3) :: V,r,Vrk4 !r and V in cartesian coordinates centered at the center of the earth
  double precision :: t,tref! for calculations
  double precision , dimension(3) ::   Vrotfinal
  double precision , dimension(3) ::  Rorig
  double precision , dimension(3) :: latlonalt, VrelFinal
  !##############################################


  ! RK suite parameters
  integer ( kind = 4 ), parameter :: neqn = 6

  double precision abserr
  !external f01
  integer  index1
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iwork(5)
  !real ( kind = 8 ) :: pi = 3.141592653589793D+00
  double precision relerr
  integer ( kind = 4 ), parameter :: step_num = 12
  !  double precision t
  double precision tout
  double precision work(100+21*neqn)
  double precision y(neqn),yold(neqn),dtintervalLocal,oldAltitudeLocal





  maxLat = -99999.9 ! initialing min and max vals
  maxLon = -99999.9
  minLon =  99999.9
  minLat =  99999.9


  finalConditions = -9999



  dtintervalLocal = dtinterval
  t=0 ! initializing time 
  tref = 0!
  r(1) = initialState(1)
  r(2) = initialState(2)
  r(3) = initialState(3)

  V(1) = initialState(4) + DebrisVel(1) ! debris velocity due to explosion added as an impulse
  V(2) = initialState(5) + DebrisVel(2)
  V(3) = initialState(6) + DebrisVel(3)
  Vrk4 = V
  Y(1) = r(1)
  Y(2) = r(2)
  Y(3) = r(3)
  Y(4) = V(1)
  Y(5) = V(2)
  Y(6) = V(3)
  tout = 0.0
  latlonalt = latLonAltCalculation(y(1:3),tout,thetag0,PlanetModel)

  latitude = latlonalt(1)
  longitude = latlonalt(2)
  altitudeFinal = latlonalt(3)
  !print *,'CD is ',CD

  if(GEOptions.EQ.1)then 

     open(unit=7,file=filename,status='unknown')
     write(7,*) 'Longitude(deg) Latitude(deg) Altitude(m) Time(sec)'
     write(7,'(4F25.13)')  longitude,latitude,altitudeFinal,tref

  endif

  abserr = 0.001D+00
  relerr = 0.001D+00

  iflag = 1

  do while ((latlonalt(3)>=0.0).and.(tout <=5.0*3600)) ! stop when it reaches the ground, or 5 hours and still not on the ground
     tout = tout + dtintervalLocal
     index1 = index1 + 1
     finalConditions(6) = latlonalt(3)
     yold = y
     oldAltitudeLocal = latlonalt(3)
     call ode (F, neqn, y, t, tout, relerr, abserr, iflag, work, iwork )

     if ( iflag /= 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST01 - Fatal error in debrisProp.f90!'
        write ( *, '(a,i8)' ) '  ODE returned IFLAG = ', iflag
        exit
     endif
     latlonalt = latLonAltCalculation(y(1:3),tout,thetag0,PlanetModel)

     if(GEOptions.EQ.1)then 
        open(unit=7,file=filename,status='unknown')
        write(7,'(4F25.13)') latlonalt(2),latlonalt(1),latlonalt(3),tout
     endif
     
    !if ((latlonalt(3)<=-10.0).or.((oldAltitudeLocal>=30.0).and.(latlonalt(3)<0))) then
        !print *,latlonalt(3)
    !    tout = tout - dtintervalLocal
    !    y = yold
    !    latlonalt = latLonAltCalculation(y(1:3),tout,thetag0,PlanetModel)
    !    print *,'adjusting tout'
    !    if (dtintervalLocal<=0.00001)then
    !       print *, 'dt too small. Current dt = ',dtintervalLocal
    !       exit
    !    else
    !       dtintervalLocal = .5*dtintervalLocal
               
    !    endif
    ! else
    !    if(GEOptions.EQ.1)then
    !      open(unit=7,file=filename,status='unknown')
    !      write(7,'(4F25.13)') latlonalt(2),latlonalt(1),latlonalt(3),tout
    !    endif

    ! endif




  end do



  flag = 0

  if (latlonalt(3)>0)then
     flag = 2
  end if


  Vrotfinal(1) = -y(2)*omega
  Vrotfinal(2) = y(1)*omega
  Vrotfinal(3) = 0
  VrelFinal = (/y(4),y(5),y(6)/) - Vrotfinal
  VrelMag = norm(VrelFinal) 

  !finalConditions = (/altitudeFinal,latitudeFinal,longitudeFinal,VrelMag,flag/) 
  finalConditions = (/latlonalt(3),latlonalt(1),latlonalt(2),VrelMag,flag,finalConditions(6),dtintervalLocal/)
 
  !print *,'Accel ',norm(y(4:6) - (/-(y(5)*omega),y(4)*omega,0/))
  !print *,'Altitude ',latlonalt(3)
  !print *, 'Vrot ',norm(Vrotfinal)
  !print *, 'VrelMag ', VrelMag
  !print *, 'VFinal ',norm(y(1:3))
  !print *, 'Yvec ',y(1:3)

  close(unit=7)


  !*********

contains 


  subroutine F(ccT,Y,YP)
    implicit none
    double precision :: ccT
    double precision :: Y(6),YP(6)
    double precision , dimension(3) :: resdmdvdt

    resdmdvdt = dvdt(Y(1:3),Y(4:6),ccT)
    YP(1) = Y(4)
    YP(2) = Y(5)
    YP(3) = Y(6)
    YP(4) = resdmdvdt(1)
    YP(5) = resdmdvdt(2)
    YP(6) = resdmdvdt(3)

    return
  END subroutine F



  function dvdt(r,V,ct)
    use constants
    !use GRAMAtmos
    double precision , intent(in) :: ct
    double precision , dimension(3) , intent(in) :: r, V
    double precision , dimension(3) :: weight, drag, dvdt, Vinf, Vrot, interVec, liftDir,lift, r2, latlonalt,dragDir
    double precision , dimension(3) :: VwindLocal,Vwind
    double precision :: altitude , rho, g, newden,newU,newV,newW!,rotAngle
    double precision :: clat,clon,Minf,speedOfSound,CDLocal,CLLocal,VinfMag
    double precision , dimension(3,3) :: localRotation
    ! print *, 'Solving dvdt'   



    latlonalt = latLonAltCalculation(r,ct,thetag0,PlanetModel)!
    clat = latlonalt(1)
    clon = latlonalt(2)
    altitude = latlonalt(3)


    !rotAngle = ct*angRate
    !altitude = norm(r) - planetRadius
    Vrot(1) = -r(2)*omega
    Vrot(2) = r(1)*omega
    Vrot(3) = 0
    !rho = density(altitude)
    Vwind(1)=0
    Vwind(2)=0
    Vwind(3)=0


    if((atmosOption==2).OR.(atmosOption==1)) then

       call atmosInterpolation(newden,newU,newV,newW,altitude,altitudeList,densityList,UList,VList,WList,nList)


       if (altitude>altitudeList(1)) then
          if (altitude<200000) then !above 150 km set density to zero
             rho = (densityList(1)/density(altitudeList(1)))*density(altitude)
             !rho = 0
             !altitude = altitudeList(1)
             !   print *, 'Using exponential model for density'
          else
             rho = 0.0
          end if
       else
          rho = newden !using density given by user
       end if
       !      rho = density(altitude)

       if (altitude>200000.) then
          rho = 0.0
       end if




       if (atmosOption==2) then

          VwindLocal = (/newW,newU,newV/)
          localRotation = transpose(matmul(rRotation(2,-clat*Pi/180),rRotation(3,thetag0+omega*ct+clon*Pi/180)))
          Vwind = matmul(localRotation,VwindLocal)

       end if
    elseif (atmosOption==0) then

       rho = density(altitude)

    elseif (atmosOption ==-1) then
       rho = 0

    else
       print *, 'Error in debrisPropagation: Unknown atmosOption Flag'
       STOP

    end if



    !print *,rho-newden
    g = gravity(altitude)
    weight = mass*(-g*r/norm(r))
    Vinf = V- Vwind-Vrot
    VinfMag = norm(Vinf)
    call getspeedofSound(altitude,speedOfSound)
    Minf = VinfMag/speedOfSound
    !print *,Minf

    if (atmosOption.NE.-1) then

       ! Drag Calculation
       if (nCD>1) then !using CD as a function of Minf

          call LinearCDInterpolation(MinfCD,CD,nCD,Minf,CDLocal)
          ! print *,CDlocal
       else
          CDlocal = CD(1)
       end if

       dragdir = -Vinf !not a unit vector
       drag = 0.5*CDLocal*rho*Sref*VinfMag*dragdir


       !  drag = drag*(1-0.3*sin(rotAngle))

       !! Lift calculation
       !! getting lift direction
       interVec = cross(r,dragdir)
       if (norm(interVec)> 0.00001) then
          interVec = intervec/norm(intervec)
          liftDir = cross(dragdir,interVec)
          liftDir = liftDir/norm(liftDir)
       else
          r2(1) = r(1) - Rorig(1)
          r2(2) = r(2) - Rorig(2)
          r2(3) = r(3) - Rorig(3)
          interVec = cross(r2,dragdir)
          liftDir = cross(dragdir,interVec)
          liftDir = liftDir/norm(liftDir)
          !print *, 'Switching to offset r for LIFT \n'

       end if

       if (CLOption==1)then

          if (nCL>1)then
             call LinearCLInterpolation(MinfCL,CL,nCL,Minf,CLLocal)



          else
             CLLocal = CL(1)
          end if
          lift = 0.5D0*CLLocal*rho*Sref*VinfMag*VinfMag*liftDir
       else ! using constant L over D
          !lift = (LoverD*norm(drag))*liftDir
          lift = (LoverD*CDLocal*0.5D0*rho*VinfMag**2*Sref)*liftDir
          !lift = -norm(drag)*LoverD*liftDir
       end if

    else
       drag = (/0,0,0/)
       lift = (/0,0,0/)
    end if


    ! liftRot = lift*cos(rotAngle) - interVec*sin(rotAngle)*norm(lift) ! this line if a rotating lift vector is desired

    !liftRot = sin(rotAngle)*lift ! this line if a sinusoidal lift vector wrt Time is desired. Set rotAngle to zero for regular lift
    !dvdt = (drag + liftRot + weight)/mass

    dvdt = (drag + lift + weight)/mass
    ! print *,norm(weight)
    !print *,Cdlocal,speedofSound,Minf

  end function dvdt

  function density(altitude)
    use constants
    implicit none

    double precision, intent(in) :: altitude
    double precision :: density

    !  g = gravity(altitude)
    density = rho0*exp(-altitude/Hscale)


  end function density


  function gravity(altitude)
    use constants
    implicit none

    double precision , intent(in) :: altitude
    double precision :: gravity

    gravity = muPlanet/((altitude+planetRadius)**2.D0)
    ! print *, gravity
  end function gravity


  function norm(x)
    use constants
    double precision, dimension(3), intent(in) :: x
    double precision :: norm

    norm = sqrt(x(1)**2.D0 + x(2)**2.D0 + x(3)**2.D0)

  end function norm



  function latlonaltCalculation(r,t,thetag0,model)
    use constants
    use planetEarthCharacteristics
    implicit none
    ! calculates spherical lat lon 
    ! uses ellipse for altitude if model = 1

    double precision , dimension(3) , intent(in) :: r
    integer , intent(in) :: model
    double precision , intent(in) :: t,thetag0
    double precision :: normr, rdelta,phigd,Nphi,phigd_np1,tol,height
    double precision , dimension(3):: rlocal, latlonaltCalculation

    tol = 999.
    rlocal =  matmul(Rrotation(3,thetag0+omega*t),r) ! accounting for Rotating planet
    normr = norm(rlocal)
    height = 0.0

    latlonaltCalculation(2) = atan2(rlocal(2),rlocal(1))*180.0D0/Pi  !Longitude in degrees

    phigd = asin(rlocal(3)/normr) !Latitude in radians

    if (model==1)then
       rdelta = sqrt(rlocal(1)**2+rlocal(2)**2)
       do while (tol>1e-7)
          Nphi = Requator/sqrt(1-ecc**2*(sin(phigd))**2)
          phigd_np1 = atan((rlocal(3) + Nphi*ecc**2*sin(phigd))/rdelta)
          tol = abs(phigd - phigd_np1)
          phigd = phigd_np1
       end do
       height = rdelta/cos(phigd) - Nphi
    elseif (model==0)then
       height = normr - planetRadius
    end if
    latlonaltCalculation(1) = phigd*180.0/pi ! latitude in degrees
    latlonaltCalculation(3) = height

  end function latlonaltCalculation










  function Rrotation(axis,angle) ! rotation from N frame to axis desired cRn
    use constants
    implicit none
    double precision , intent(in) :: angle
    integer , intent(in) :: axis
    double precision , dimension(3,3) :: Rrotation

    if (axis==1) then
       Rrotation(1,1) = 1
       Rrotation(1,2) = 0
       Rrotation(1,3) = 0
       Rrotation(2,1) = 0
       Rrotation(2,2) = cos(angle)
       Rrotation(2,3) = sin(angle)
       Rrotation(3,1) = 0
       Rrotation(3,2) = -sin(angle)
       Rrotation(3,3) = cos(angle)
    elseif (axis ==2) then
       Rrotation(1,1) = cos(angle)
       Rrotation(1,2) = 0
       Rrotation(1,3) = -sin(angle)
       Rrotation(2,1) = 0
       Rrotation(2,2) = 1
       Rrotation(2,3) = 0
       Rrotation(3,1) = sin(angle)
       Rrotation(3,2) = 0
       Rrotation(3,3) = cos(angle)
    elseif (axis ==3) then
       Rrotation(1,1) = cos(angle)
       Rrotation(1,2) = sin(angle)
       Rrotation(1,3) = 0
       Rrotation(2,1) = -sin(angle)
       Rrotation(2,2) = cos(angle)
       Rrotation(2,3) = 0
       Rrotation(3,1) = 0
       Rrotation(3,2) = 0
       Rrotation(3,3) = 1
    end if
  end function Rrotation

  function cross(a,b) 
    use constants
    implicit none

    double precision , dimension (3), intent(in) :: a,b
    double precision, dimension (3)::cross

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  end function cross





end subroutine debrisPropagation


