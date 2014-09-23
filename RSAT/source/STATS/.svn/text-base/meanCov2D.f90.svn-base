subroutine meanCovMatrix(mean,cov,xy,nelements)
  implicit none


  integer , intent(in) :: nelements
  double precision , dimension(nelements,2),intent(in) :: xy
  double precision , dimension(2), intent(out) :: mean
  double precision , dimension(2,2), intent(out) :: cov

  call getMean(xy,nelements,mean)
  call getCovar(mean,xy,nelements,cov)

end subroutine meanCovMatrix

subroutine getMean(xy,nelements,mean)
  implicit none
  integer , intent(in) :: nelements
  integer :: i
  double precision , dimension(nelements,2),intent(in) :: xy
  double precision , dimension(2) , intent(out) :: mean

  mean = (/0.0,0.0/)
  do i = 1,nelements
     mean(1) = mean(1) + xy(i,1)
     Mean(2) = mean(2) + xy(i,2)
  end do
  mean = mean/(nelements)

end subroutine getMean

subroutine getCovar(mean,xy,nelements,covar)
  implicit none
  integer , intent(in) :: nelements
  integer :: i
  double precision , dimension(nelements,2),intent(in) :: xy
  double precision , dimension(2) , intent(in):: mean
  double precision , dimension(2,2),intent(out) :: covar

  covar(1,1) = 0
  covar(1,2) = 0
  covar(2,1) = 0
  covar(2,2) = 0

  do i = 1,nelements
     covar(1,1) = covar(1,1) + (xy(i,1) - mean(1))**2
     covar(2,2) = covar(2,2) + (xy(i,2) - mean(2))**2
     covar(1,2) = covar(1,2) + (xy(i,1) - mean(1))*(xy(i,2)-mean(2))

  end do

  covar = covar/real(nelements - 1)
  covar(2,1) = covar(1,2)
end subroutine getCovar



subroutine getEigenVector(covar,eigenVector)
  implicit none
  double precision , dimension(2,2) , intent(in):: covar
  double precision , dimension(2) , intent(out) :: eigenVector
  double precision , dimension(2) :: eigenvals
  double precision :: v1,v2,normval,rho,sigma1,sigma2
  sigma1 = sqrt(covar(1,1))
  sigma2 = sqrt(covar(2,2))
  rho = covar(1,2)/(sigma1*sigma2)
  !print *,abs(rho)
  if (abs(rho)<1e-14) then ! this means there is no correlation, so already in principal axis 
     eigenVector(1) = 1.0
     eigenVector(2) = 0.0
  else
     call getEigenvals(covar,eigenvals)

     v1 = 1
     v2 = -(covar(1,1)-eigenvals(1))/(covar(1,2))
     normval = sqrt(v1**2+v2**2)
     eigenVector(1) = v1/normval
     eigenVector(2) = v2/normval

     !print *,eigenvector
     !print *,eigenvals
  end if
end subroutine getEigenVector



subroutine getEigenVals(covar,eigenvals)
implicit none
double precision , dimension(2,2) , intent(in):: covar
double precision , dimension(2) , intent(out) :: eigenVals
double precision :: a,b,c

a = 1
b = -(covar(1,1)+covar(2,2))
c = covar(1,1)*covar(2,2)-covar(1,2)*covar(1,2)

eigenvals(1) =  (-b+sqrt(b**2-4*a*c))/(2*a)
eigenvals(2) =  (-b-sqrt(b**2-4*a*c))/(2*a)
end subroutine getEigenVals

subroutine getSigma(vector,nelements,sigma)
! routine to get standard deviation
  implicit none
  integer , intent(in) :: nelements
  double precision , dimension(nelements) , intent(in) :: vector
  double precision , intent(out) :: sigma
  integer :: i
  double precision :: xbar

  xbar = 0
  sigma = 0
  
  do i = 1,nelements
     xbar = vector(i) + xbar
  end do

  xbar = xbar/nelements
  do i = 1,nelements
     sigma = sigma + (vector(i)-xbar)**2
  end do
  sigma = sigma/(nelements - 1)
  sigma = sigma**.5
end subroutine getSigma

subroutine getPQ(xy,nelements,U,P,Q)
  implicit none
  integer , intent(in) :: nelements
  double precision , dimension(nelements,2) , intent(in) :: xy
  double precision , dimension(2) , intent(in) :: U
  double precision , dimension(nelements) , intent(out) :: P,Q
  integer :: i

  if (abs(U(2))>=1e-15) then
  do i=1,nelements
     P(i) = xy(i,1)*U(1) + xy(i,2)*U(2)
     Q(i) = -xy(i,1)*U(2) + xy(i,2)*U(1)
  end do
  else
     P = xy(:,1)
     Q = xy(:,2)
  end if

end subroutine getPQ

subroutine getIQR(Pin,nelements,val)
  implicit none
  integer , intent(in) :: nelements
  double precision , dimension(nelements) , intent(in) :: Pin
  double precision , dimension(nelements) :: P
  double precision , intent(out):: val
  integer :: mid,low,high

  P = Pin
  !call Qsort(P)
  call shell(nelements,P)
!SUBROUTINE SHELL(N,ARR)

  !print *,P,nelements
  if (mod(nelements,2)==0) then
     mid = nelements/2
     low = mid/2
     high = nelements - low
     val = .5*(P(high+1)+P(high+2)) - .5*(P(low+1)+P(low+2))

  else
     mid = (nelements + 1) /2
     low = mid/2
     high = nelements -low + 1
     val = P(high) - P(low)
  end if
  
end subroutine getIQR

  


subroutine H2Parameters(xy,nelements,H2inv,detH2)
  implicit none
  integer , intent(in) :: nelements
  double precision , dimension(nelements,2) , intent(in) :: xy
  double precision , dimension(3) , intent(out) :: H2inv
  double precision , intent(out) :: detH2
  double precision , dimension(2) :: U , average !eigenvector info
  double precision :: sigma1,sigma2,h1,h2,IQR1,IQR2,multVal,tempvalsN
  double precision , dimension(2,2) :: Covar
  double precision , dimension(nelements) :: P,Q


  call getMean(xy,nelements,average)
  call getCovar(average,xy,nelements,Covar)
  call getEigenVector(Covar,U)
  call getPQ(xy,nelements,U,P,Q)
  call getSigma(P,nelements,sigma1)
  call getSigma(Q,nelements,sigma2)

  call getIQR(P,nelements,IQR1)
  call getIQR(Q,nelements,IQR2)

  !print *,"average",average
  !print *,"variance",Covar
  !print *,"U",U
  !print *,U
  tempvalsN = nelements
  multVal = tempvalsN**(-2.0D-1)
  h1 = 1.06D0*min(sigma1,(IQR1/1.34D0))*multVal
  h2 = 1.06D0*min(sigma2,(IQR2/1.34D0))*multVal

  !print *, "h vals", h1,h2
  !print *, "sigmas", sigma1, sigma2
  !print *, "mult",nelements**(-.2),multVal,nelements
  !print *, 'IQR', IQR1, IQR2
  detH2 = (h1**2.0D0)*(h2**2.0D0)
 ! print *,detH2
  H2inv(1) = U(1)**2.0D0/(h1**2.0D0) + U(2)**2.0D0/(h2**2.0D0)
  H2inv(2) = U(2)**2.0D0/(h1**2.0D0) + U(1)**2.0D0/(h2**2.0D0)
  H2inv(3) = U(1)*U(2)/(h1**2.0D0) - U(1)*U(2)/(h2**2.0D0)

  !print *,h1,h2


end subroutine H2Parameters



subroutine normalBivariate(mean,cov,X,Y,Z,cols,rows)
      implicit none
      integer , intent(in) :: rows,cols
      double precision , parameter :: Pi = 3.141592653589793238462643D0

      double precision , dimension(2) ,intent(in):: mean
      double precision , dimension(2,2) , intent(in) :: cov

      double precision , dimension(rows,cols) , intent(in) :: X,Y
      double precision , dimension(rows,cols) , intent(out) :: Z
      double precision :: rho,xval,yval,sigmax,sigmay,sigmaxinv,sigmayinv,sigmaxyinv
      integer :: i,j

      sigmax = (cov(1,1))**.5
      sigmay = (cov(2,2))**.5
      sigmaxinv = 1./(sigmax)
      sigmayinv = 1./(sigmay)
      sigmaxyinv = 1./(sigmax*sigmay)
      rho = cov(1,2)/(sigmax*sigmay)

      do i = 1,rows
         do j = 1,cols
            xval = X(i,j)
            yval = Y(i,j)
            Z(i,j) = exp((-0.5D0)/(1.0D0-rho**2.0D0)*((sigmaxinv*(xval-mean(1)))**2.0D0 + &
                 (sigmayinv*(yval - mean(2)))**2.0D0 - 2.0D0*sigmaxyinv*rho*(xval-mean(1))*(yval-mean(2))))
         end do
      end do
      Z = Z/(2.0D0*Pi*sigmax*sigmay*(1.0D0-rho**2.0D0)**.5D0)
      
    end subroutine normalBivariate

subroutine KDE(xy,nelements,X,Y,Z,cols,rows)
  implicit none
  integer , intent(in) :: nelements , rows,cols
  double precision , dimension(nelements,2) , intent(in) :: xy 
  double precision , dimension(rows,cols) , intent(in) :: X,Y
  double precision , dimension(rows,cols) , intent(out) :: Z
  double precision , parameter :: Pi = 3.141592653589793238462643D0
  double precision :: xHx
  double precision , dimension(3) :: H2inv
  double precision :: detH2 , tempval, denom,tempElements
  integer :: i,j,k
  call H2Parameters(xy,nelements,H2inv,detH2)
  tempElements = nelements
  denom = 2.0D0*Pi*(detH2**.5D0)*tempElements
  !print *,"Denominator Fortran",denom
  do i=1,rows
     do j = 1,cols
        tempval = 0.0
        do k = 1,nelements
           xHx = ((X(i,j)-xy(k,1))**2.0D0)*H2inv(1) + 2.0D0*(X(i,j)-xy(k,1))*(Y(i,j)-xy(k,2))*H2inv(3) +&
                ((Y(i,j)-xy(k,2))**2.0D0)*H2inv(2)
           tempval = exp(-0.50D0*xHx) + tempval
        end do
        Z(i,j) = tempval/denom
     end do
  end do




end subroutine KDE

recursive subroutine Qsort(A)
  double precision, intent(in out), dimension(:) :: A
  integer :: iq,lenA
  lenA = size(A)
  if(size(A) > 1) then
    print *,'here1',A

     call Partition(A,lenA, iq)

     call Qsort(A(:iq-1))
     call Qsort(A(iq:))
  endif
end subroutine Qsort

subroutine Partition(A,lenA, marker)
  integer , intent(in) :: lenA
  double precision, intent(in out), dimension(lenA) :: A
  integer, intent(out) :: marker
  integer :: i, j
  double precision :: temp
  double precision :: x      ! pivot point

  x = A(1)
  i= 0

  j= size(A) + 1
  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition



SUBROUTINE SHELL(N,ARR)
  parameter(ALN2I=1./0.69314718,TINY=1.E-5)
  double precision :: ARR(N),t
  LOGNB2=INT(ALOG(FLOAT(N))*ALN2I+TINY)
  m=n
  do nn=1,LOGNB2
     m=m/2; k=n-m
     do j=1,k
        i=j
10      continue
        l=i+m
        if(ARR(l).LT.ARR(i)) then
           t=ARR(i)
           ARR(i)=ARR(l)
           ARR(l)=t
           i=i-m
           if(i.GE.1) GOTO 10
        end if
     end do
  end do
  
END SUBROUTINE SHELL







!!$subroutine areaOfInterest(mean,cov,nsigma)
!!$  implicit none
!!$  double precision , dimension(2) ,intent(in):: mean
!!$  double precision , dimension(2,2) , intent(in) :: cov
!!$  double precision , intent(in) :: nsigma
!!$  double precision :: sigmax, sigmay
!!$  double precision :: 
!!$  integer :: nvals
!!$
!!$  sigmax = (cov(1,1))**.5
!!$  sigmay = (cov(2,2))**.5
!!$
!!$
!!$
!!$
!!$
!!$end subroutine areaOfInterest
!!$
!!$subroutine myRange(minval,dval,nvals,rangeout)
!!$  double precision , intent(in) :: minval,maxval,dval
!!$  integer , intent(in):: nvals
!!$  double precision , dimension(nvals) , intent(out) :: rangeout
!!$  integer ::i
!!$  !nvals = ceiling((maxval-minval)/dval)
!!$  rangeout(1) = minval
!!$  do i = 2,nvals
!!$     rangeout(i) = rangeout(i-1) + dval
!!$  end do
!!$
!!$end subroutine myRange

  
