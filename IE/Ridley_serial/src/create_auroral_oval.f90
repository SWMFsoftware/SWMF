!=============================================================================
subroutine create_auroral_oval(currentIn, thetaIn, psiIn, &
     ColatOut, WidthOut, StrenOut)

  ! This subroutine replaces the original oval algorithm.  From the FAC
  ! pattern, it creates an auroral oval that is roughly fit to the upward
  ! FACs.  The strength and day/night offset are set by using MLT-binned
  ! regional averages.  This version of the subroutine is meant to follow
  ! the methodology of the original closely except to prevent aliasing of
  ! the oval position as a function of time.  The main changes are to
  ! include all longitudes in the calculation of the oval position and to
  ! weight the day-night shift by the strength of the binned FAC values.

  use ModNumConst,    ONLY: cDegToRad, cRadToDeg, cHalfPi, cPi
  use ModIonosphere,  ONLY: IONO_nTheta, IONO_nPsi
  use IE_ModIo,       ONLY: NameIonoDir
  use ModConductance, ONLY: DoOvalShift, DoFitCircle
  use ModIoUnit,      ONLY: UnitTMP_
  use IE_ModMain,     ONLY: Time_Array, nSolve
  
  implicit none

  ! Input/Output variables:  
  real, intent(in),  dimension(IONO_nTheta,IONO_nPsi):: CurrentIn,thetaIn,psiIn
  real, intent(out), dimension(IONO_nPsi) :: ColatOut, widthOut, strenOut

  ! Working variables:
  logical :: IsNorth
  integer :: i, j, jMax(1), dJ, iStart, iLonStart, iLonEnd
  real    :: dLat, FacSum, FacMidn, FacNoon, FacDawn
  real    :: WidthMean, WidthDay, ColatDev, ColatMean, ColatNoon, ColatMidn
  real, dimension(IONO_nTheta,IONO_nPsi):: current, theta, psi
  real, dimension(IONO_nPsi) :: FacMax=0, ColatMax=0, Width=0
  real, dimension(8)         :: FacRegional=0, ColatRegional=0, ColatScaled=0

  ! Variables for the gradient descent fit
  integer, parameter:: MaxIter = 20
  integer:: nIter
  real, dimension(IONO_nPsi):: x_I, xShift_I, y2_I, r_I, w_I
  real:: dStep, Tolerance2, ColatMeanPrev, ColatDevPrev
  real:: GradMean, GradDev, GradMeanPrev, GradDevPrev

  ! Testing variables:
  logical       :: DoTest, DoTestMe
  logical, save :: IsFirstWrite = .true.
  character(len=100), save    :: NameFile, StringFormat
  character(len=*), parameter :: NameSub = 'create_auroral_oval'
  !--------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  ! Initialize arrays to zero:
  FacMax=0;  ColatMax=0; Width=0
  
  ! Reverse arrays if in Southern hemisphere.
  if (ThetaIn(1,1) < cHalfPi) then
     Current = CurrentIn
     Theta   = ThetaIn
     Psi     = PsiIn
     IsNorth = .true.
  else
     do i = 1, IONO_nTheta
        do j = 1, IONO_nPsi
           Current(IONO_nTheta-(i-1), j) = CurrentIn(i,j)
           Theta(  IONO_nTheta-(i-1), j) = cPi - ThetaIn(i,j)
           Psi(    IONO_nTheta-(i-1), j) = PsiIn(i,j)
        enddo
     enddo
     IsNorth = .false.
  endif

  ! Calculate some useful variables:
  dLat = theta(2,1) - theta(1,1)
  
  ! Find the max upward FAC and the corresponding colatitude:
  do j=1, IONO_nPsi        ! Loop through longitude
     do i=4, IONO_nTheta   ! Loop through colat, start away from pole
        if(current(i,j)>FacMax(j)) then
           FacMax(j)   = current(i,j)
           ColatMax(j) = theta(i,j)
           iStart      = i
        end if
     end do

     ! At current longitude, find width of oval:
     do i = iStart, IONO_nTheta
        if(current(i,j) <= FacMax(j)/4.0) then
           Width(j) = theta(i,j) - ColatMax(j)
           exit
        endif
     end do

     ! Limit width of oval: 
     if (Width(j) < dLat) Width(j) = ColatMax(j)/5.0
  end do

  ! Sums and means across longitudes, set floors & ceilings:
  FacSum    = sum(FacMax)
  FacDawn   = FacSum
  ColatMean = max(sum(ColatMax*FacMax)/FacSum, 15*cDegToRad)
  WidthMean = min(sum(Width*FacMax)   /FacSum,  6*cDegToRad)
  WidthDay  = max(WidthMean/2.0, cDegToRad)

  if(DoTestMe)then
     write(*,*) NameSub//" DEBUG!-----------------------------"
     write(*,*) "FacSum =    ", FacSum
     write(*,*) "ColatMean = ", ColatMean * cRadToDeg
     write(*,*) "WidthMean = ", WidthMean * cRadToDeg
     write(*,*) "WidthDay  = ", WidthDay  * cRadToDeg
  end if


  ! Calculate a day-night shift in oval location.
  ! Set floor for oval size.
  if (ColatMean < 15.*cDegToRad) then
     ! Small oval: set floor and no shift.
     ColatMean = 15.*cDegToRad
     ColatDev  = 0.0
     ColatOut = ColatMean
  elseif (.not.DoOvalShift) then
     ! No shift if DoOvalShift is false. 
     ColatDev = 0.0
     ColatOut = ColatMean
  elseif(DoFitCircle)then
     ! Fit a circle to the points defined by Psi(j), ColatMax(j) with weight FaxMax(j)/FacSum
     ! Use steepest gradient search method to minimize the radial distance
     ! of the points to the circle.

     ! Precalculate some useful quantities. We flatten out the sphere, and
     ! Colat is taken as the radial coordinate, while
     ! Psi is taken as the angular coordinate of a polar coordinate system.
     x_I  = cos(Psi(1,:))*ColatMax*cRadToDeg       ! x coordinates of points
     y2_I = (sin(Psi(1,:))*ColatMax*cRadToDeg)**2  ! y coordinates squared
     w_I  = FacMax / FacSum                   ! weights normalized to 1

     dStep = 5.0       ! initial step size (in degrees)
     Tolerance2 = 0.01 ! square of step size limit for exiting (in degrees)

     ColatMean = ColatMean*cRadToDeg         ! initial guess in degrees
     ColatDev  = 0.0                         ! initial guess

     ! Barzilai-Borwein method: https://en.wikipedia.org/wiki/Gradient_descent
     ! works for convex functions
     do nIter = 1, MaxIter
        ! Calculate derivatives of the distance function
        ! f = sum( weight * (r - ColatMean)^2 )
        ! where r = sqrt( (x - s)^2 + y^2)
        ! with respect to the parameters ColatMean and ColatDev
        xShift_I = x_I - ColatDev
        r_I = max(sqrt(xShift_I**2 + y2_I), 1e-20)
        GradMean = -2*sum(w_I*(r_I - ColatMean))              ! df/dColatMean
        GradDev  = -2*sum(w_I*(r_I - ColatMean)*xShift_I/r_I) ! df/dColatDev

        if(nIter == 1)then
           dStep = dStep / sqrt(GradMean**2 + GradDev**2)
        else
           dStep = ((ColatMean - ColatMeanPrev)*(GradMean - GradMeanPrev)    &
                +   (ColatDev  - ColatDevPrev )*(GradDev  - GradDevPrev ))   &
                /  ((GradMean - GradMeanPrev)**2 + (GradDev - GradDevPrev)**2)
        end if
        ! Store previous values
        ColatMeanPrev = ColatMean
        ColatDevPrev  = ColatDev

        ! Step in the direction of -gradient
        ColatMean = ColatMean - dStep*GradMean
        ColatDev  = ColatDev  - dStep*GradDev
     
        ! Store previous gradients
        GradMeanPrev = GradMean
        GradDevPrev  = GradDev

        if(DoTestMe) &
             write(*,*) NameSub,' nIter, ColatMeanPrev, ColatDevPrev, ColatMean, ColatDev=', &
             nIter, ColatMeanPrev, ColatDevPrev, ColatMean, ColatDev

        if(Tolerance2 > dStep**2 * (GradMean**2 + GradDev**2)) EXIT

     end do
     if(nIter > MaxIter)then
        write(*,*) NameSub,' failed to converge!!!'
        write(*,*) NameSub,' ColatMeanPrev, ColatMean=', ColatMeanPrev, ColatMean
        write(*,*) NameSub,' ColatDevPrev,  ColatDev =', ColatDevPrev, ColatDev
     end if

     if(DoTestMe)then
        write(*,*) NameSub,' converged in ',nIter,' iterations'
        write(*,*) NameSub,' ColatMean, ColatDev =', ColatMean, ColatDev
     end if

     ! ColatDev cannot be too large, limit it
     if(abs(ColatDev) > ColatMean/2)then
        ColatDev = sign(ColatMean/2, ColatDev)
        if(DoTestMe) write(*,*) NameSub,' limited ColatDev=', ColatDev
     end if

     ! Transform back to radians
     ColatMean = ColatMean*cDegToRad
     ColatDev  = ColatDev*cDegToRad

     ! Calculate distance from the center (ColatOut) for a given Psi value
     ! This requires a bit of geometry...
     ColatOut = ColatDev*cos(Psi(1,:)) + sqrt(ColatMean**2 - ColatDev**2*sin(Psi(1,:)))

  else
     ! Generate regional values (dusk, dawn, etc.):
     dJ = int(IONO_nPsi/8.)
     do i = 1, 8
        ! Get indices of longitude window:
        iLonStart = max(1, (i-1)*dJ)
        iLonEnd   = min(i*dJ, IONO_nPsi)
        ! Find maximum value and associated width inside window:
        jMax = maxloc( FacMax(iLonStart:iLonEnd) )
        FacRegional(i)   = FacMax(  jMax(1)+iLonStart)
        ColatRegional(i) = ColatMax(jMax(1)+iLonStart)
     end do

     ! Weighted regional colat values:
     ColatScaled = FacRegional * ColatRegional
     ColatNoon  = (ColatScaled(1)+ColatScaled(8))/(FacRegional(1)+FacRegional(8))
     ColatMidn  = (ColatScaled(4)+ColatScaled(5))/(FacRegional(4)+FacRegional(5))

     ! Regional FAC strengths:
     FacNoon = ( FacRegional(1)+FacRegional(8)  )/2.0
     FacMidn = ( FacRegional(4)+FacRegional(5)  )/2.0
     FacDawn = ( FacRegional(6)+FacRegional(7)  )/2.0

     ! Oval shift is average of day and night oval limits.
     ColatDev=( FacNoon*(ColatNoon-ColatMean) - FacMidn*(ColatMidn-ColatMean) )&
          / (FacNoon+FacMidn)
     ! Restrict shift to reasonable range:
     if (abs(ColatDev) > ColatMean/2.0) &
          ColatDev = ColatDev*ColatMean/2.0/abs(ColatDev)

     ! This is not a circle
     ColatOut = ColatMean - ColatDev*cos(Psi(1,:))  ! SIGN IS WRONG.

  end if
  
  ! Use above values to craft auroral oval:
  WidthOut = WidthDay  + (WidthMean-WidthDay)*sin(Psi(1,:)/2.0)
  StrenOut = FacDawn   +     (FacSum-FacDawn)*sin(Psi(1,:)/2.0)

  ! For testing purposes, write oval info to file.
  if(.not.DoTestMe .or. .not. IsNorth) return

  write(NameFile,'(a,i8.8,a)')trim(NameIonoDir)//'aurora_n',nSolve,'.dat'
  open(unit=UnitTmp_, file=NameFile, status='replace')
  write(UnitTmp_, '(a,2f8.2)')'Auroral oval fit: ColatMean, ColatDev=', &
       ColatMean*cRadToDeg, ColatDev*cRadToDeg
  write(UnitTmp_,'(a)') 'j Phi ThetaMax ThetaFit Weight Width'
  do j=1, IONO_nPsi
     write(UnitTmp_,'(i4, 5es13.4)') j, Psi(1,j)*cRadToDeg, ColatMax(j)*cRadToDeg, &
          ColatOut(j)*cRadToDeg, facMax(j)/FacSum, Width(j)*cRadToDeg
  end do
  close(UnitTmp_)
  
!!!  if(IsFirstWrite)then
!!!     ! Open file:
!!!     write(NameFile,'(a,i8.8,a)')trim(NameIonoDir)//'aurora_n',nSolve,'.txt'
!!!     open(unit=UnitTmp_, file=NameFile, status='replace')
!!!     ! Write header w/ longitudes:
!!!     write(UnitTmp_, '(a)', advance='NO')'Auroral oval location at Lon='
!!!     do j=1, IONO_nPsi
!!!        write(UnitTmp_,'(1x,f6.2)', advance='NO') Psi(1,j)*cRadToDeg
!!!     end do
!!!     write(UnitTmp_,'(a)')'','YYYY MM DD HH MN SS msc oval_colat'
!!!
!!!     ! Set format codes for writing output:
!!!     write(StringFormat,'("(i4,5i3,i4,", i4,"(1x,f6.2))")') IONO_nPsi
!!!     IsFirstWrite=.false.
!!!
!!!  else
!!!     ! Open file in append mode:
!!!     open(unit=UnitTmp_, file=NameFile, status='old', position='append')
!!!  end if
!!!
!!!  ! Write record & close file:
!!!  write(UnitTmp_, StringFormat) Time_Array(1:7), ColatOut*cRadToDeg
!!!  close(UnitTmp_)
  
end subroutine create_auroral_oval
!=============================================================================
