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
  use ModConductance, ONLY: DoOvalShift
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

     if(DoTestMe) &
          write(*,'(a,i3.3,E12.3,2(f8.2))')'J, FacMax, ColatMax, Width= ', &
          j, facMax(j), ColatMax(j)*cRadToDeg, Width(j)*cRadToDeg
  end do
  
  ! Sums & means across longitudes, set floors & ceilings:
  FacSum    = sum(FacMax)
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

  ! Calculate a day-night shift in oval location.
  ! Set floor for oval size.
  if (ColatMean < 15.*cDegToRad) then
     ! Small oval: set floor and no shift.
     ColatMean = 15.*cDegToRad
     ColatDev  = 0.0
  else if (.not.DoOvalShift) then
     ! No shift if DoOvalShift is false. 
     ColatDev = 0.0
  else
     ! Oval shift is average of day and night oval limits.
     ColatDev=( FacNoon*(ColatNoon-ColatMean) - FacMidn*(ColatMidn-ColatMean) )&
          / (FacNoon+FacMidn)
     ! Restrict shift to reasonable range:
     if (abs(ColatDev) > ColatMean/2.0) &
          ColatDev = ColatDev*ColatMean/2.0/abs(ColatDev)
  end if
  
  ! Use above values to craft auroral oval:
  ColatOut = ColatMean -             ColatDev*cos(Psi(1,:))  ! SIGN IS WRONG.
  WidthOut = WidthDay  + (WidthMean-WidthDay)*sin(Psi(1,:)/2.0)
  StrenOut = FacDawn   +     (FacSum-FacDawn)*sin(Psi(1,:)/2.0)

  ! For testing purposes, write oval info to file.
  if(.not.DoTestMe .or. .not. IsNorth) return
  
  if(IsFirstWrite)then
     ! Open file:
     write(NameFile,'(a,i8.8,a)')trim(NameIonoDir)//'aurora_n',nSolve,'.txt'
     open(unit=UnitTmp_, file=NameFile, status='replace')
     ! Write header w/ longitudes:
     write(UnitTmp_, '(a)', advance='NO')'Auroral oval location at Lon='
     do j=1, IONO_nPsi
        write(UnitTmp_,'(1x,f6.2)', advance='NO') Psi(1,j)*cRadToDeg
     end do
     write(UnitTmp_,'(a)')'','YYYY MM DD HH MN SS msc oval_colat'

     ! Set format codes for writing output:
     write(StringFormat,'("(i4,5i3,i4,", i4,"(1x,f6.2))")') IONO_nPsi
     IsFirstWrite=.false.
  else
     ! Open file in append mode:
     open(unit=UnitTmp_, file=NameFile, status='old', position='append')
  end if

  ! Write record & close file:
  write(UnitTmp_, StringFormat) Time_Array(1:7), ColatOut*cRadToDeg
  close(UnitTmp_)
  
end subroutine create_auroral_oval
!=============================================================================
