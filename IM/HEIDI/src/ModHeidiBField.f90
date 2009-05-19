module ModHeidiBField

  implicit none

contains

  subroutine initialize_b_field (L, nStep, bField_I, RadialDistance_I, Length_I, dLength_I)

    integer, intent(in) :: nStep                   ! Number of points along the field line
    real, intent(in)    :: L                       ! L shell value
    real, intent(out)   :: bField_I(nStep)         ! Magnetic field values
    real, intent(out)   :: RadialDistance_I(nStep) ! Radial distance
    real, intent(out)   :: Length_I(nStep)         ! Length of the dipole field line
    real, intent(out)   :: dLength_I(nStep-1)      ! Length interval between i and i+1
    real, parameter     :: DipoleStrength =  0.32  ! nTm^-3
    real                :: dLat                    ! Latitude cell size
    real                :: LatMin                  ! Minimum Latitude 
    real                :: LatMax                  ! Maximum Latitude 
    real                :: Lat                     ! Latitude
    real                :: SinLat2, CosLat2
    real                :: LatMin2, LatMax2, LatMinMax
    real                :: Beta,dLatNew,beta1,beta2
    real, parameter     :: Epsilon = 1.e-10
    integer             :: i
    !----------------------------------------------------------------------------------

    LatMax = acos(sqrt(1./L))
    LatMin = -LatMax
    dLat   = (LatMax-LatMin)/(nStep-1)
    !uniform grid along the field line
    Lat = LatMin
    
    do i =1, nStep
       SinLat2 = sin(Lat)**2
       CosLat2 = 1.0 - SinLat2
       bField_I(i) = DipoleStrength*sqrt(1.0+3.0*SinLat2)/(L*CosLat2)**3
       RadialDistance_I(i) = L*CosLat2
       Length_I(i) = dipole_length(L,LatMin,Lat) 
       Lat = Lat + dLat
    end do

    !non-uniform grid along the field line
!!$    LatMin2 = LatMin**2
!!$    LatMax2 = LatMax**2
!!$    LatMinMax = LatMin*LatMax
!!$  
!!$    beta1 = 6.0*(-nStep+Epsilon*nStep+1.0)*(nStep-1.0) 
!!$    beta2 = nStep*(2.0*LatMax2*nStep + 2.0*LatMinMax*nStep + 2.0*LatMin2*nStep &
!!$         - 4.0*LatMinMax - LatMin2 - LatMax2)
!!$    Beta = - beta1/beta2  
!!$
!!$    do i =1, nStep
!!$       Lat = LatMin+(i-1)*dLat
!!$       SinLat2 = sin(Lat)**2
!!$       CosLat2 = 1.0 - SinLat2
!!$       bField_I(i) = DipoleStrength*sqrt(1.0+3.0*SinLat2)/(L*CosLat2)**3
!!$       RadialDistance_I(i) = L*CosLat2
!!$       Length_I(i) = dipole_length(L,LatMin,Lat) 
!!$       dLatNew =(Beta*(LatMin+(i-1)*dLat)**2+Epsilon)*dLat
!!$       Lat = Lat+dLatNew
!!$    end do


    do i = 1, nStep-1
       dLength_I(i) = Length_I(i+1) - Length_I(i)
    end do

        
  end subroutine initialize_b_field
  !==================================================================================

  subroutine find_mirror_points (nStep, nPitch, PitchAngle_I, bField_I, bMirror_I,iMirror_II)

    integer, intent(in) :: nStep                ! Number of points along the field line
    integer, intent(in) :: nPitch               ! Pitch angle grid
    real, intent(in)    :: PitchAngle_I(nPitch) ! Pitch angle values
    real, intent(in)    :: bField_I(nStep)      ! Magnetic field values
    real, intent(out)   :: bMirror_I(nPitch)    ! B magnitude at mirror points    
    real                :: bMin                 ! Minimum value of magnetic field along a field line
    integer             :: iMinB                ! Location of minimum B
    integer,intent(out) :: iMirror_II(2,nPitch) ! Location of each mirror point for all pitch angles
    integer             :: i_I(1)
    integer             :: iStep, iPitch
    !----------------------------------------------------------------------------------

    i_I  = minloc(bField_I)
    iMinB= i_I(1)
    bMin = bField_I(iMinB)

    do iPitch = 1, nPitch

       bMirror_I(iPitch) = bMin/(sin(PitchAngle_I(iPitch)))**2

       do iStep = iMinB,1, -1
          if (bField_I(iStep) >= bMirror_I(iPitch)) then 
             iMirror_II(1,iPitch) = (iStep + 1)
             EXIT
          end if
       end do

       do iStep = iMinB, nStep
          if (bField_I(iStep) >= bMirror_I(iPitch)) then 
             iMirror_II(2,iPitch) = (iStep - 1)
             EXIT
          end if
       end do

    enddo
  end subroutine find_mirror_points

  !==================================================================================

  subroutine second_adiabatic_invariant(nStep, iMirror_I, bMirror, bField_I, dLength_I,L, SecondAdiabInv)

    ! Calculate integral of sqrt((B-Bm)/Bm) between the mirror points using the
    ! trapezoidal rule

    integer             :: nStep
    integer, intent(in) :: iMirror_I(2)
    real, intent(in)    :: bMirror  
    real, intent(in)    :: dLength_I(nStep-1)
    real, intent(in)    :: bField_I(nStep) 
    real, intent(out)   :: SecondAdiabInv
    real, intent(in)    :: L
    real                :: InvL
    integer             :: iStep, iFirst, iLast
    real                :: DeltaS1, DeltaS2, b1, b2, Coeff
    !----------------------------------------------------------------------------------
    iFirst = iMirror_I(1)
    iLast  = iMirror_I(2)
    SecondAdiabInv = 0.0

    if (iFirst > iLast) RETURN

    InvL = 1.0/L
    Coeff = InvL/sqrt(bMirror)

    DeltaS1 = abs((bMirror-bField_I(iFirst))*(dLength_I(iFirst-1))/(bField_I(iFirst-1)-bField_I(iFirst)))


    SecondAdiabInv= SecondAdiabInv + Coeff*(2./3.)*DeltaS1*sqrt(bMirror-bField_I(iFirst))

    do iStep = iFirst, iLast-1
       b1 = bField_I(iStep)
       b2 =  bField_I(iStep+1)
       SecondAdiabInv = SecondAdiabInv + Coeff*(2./3.)*dLength_I(iStep)/(b1 - b2) &
            *( sqrt(bMirror  - b2)**3 - sqrt(bMirror  - b1)**3 )  

    end do

    DeltaS2 = abs((bMirror-bField_I(iLast))*(dLength_I(iLast))/(bField_I(iLast+1)-bField_I(iLast)))
    SecondAdiabInv= SecondAdiabInv + Coeff*(2./3.)*DeltaS2*(sqrt(bMirror-bField_I(iLast)))

  end subroutine second_adiabatic_invariant

  !==================================================================================
  subroutine half_bounce_path_length(nStep, iMirror_I, bMirror, bField_I, dLength_I,L, HalfPathLength)

    ! Calculate integral of ds/sqrt((B-Bm)/Bm) between the mirror points using the
    ! trapezoidal rule

    integer             :: nStep
    integer, intent(in) :: iMirror_I(2)
    real, intent(in)    :: bMirror  
    real, intent(in)    :: dLength_I(nStep-1)
    real, intent(in)    :: bField_I(nStep) 
    real, intent(out)   :: HalfPathLength
    real, intent(in)    :: L
    real                :: Inv2L
    integer             :: iStep, iFirst, iLast
    real                :: DeltaS1, DeltaS2,b1,b2, Coeff
    !----------------------------------------------------------------------------------
    iFirst = iMirror_I(1)
    iLast  = iMirror_I(2)
    HalfPathLength = 0.0

    if (iFirst > iLast) RETURN

    Inv2L = 1.0/(2.*L)

    Coeff = Inv2L*sqrt(bMirror)
    
    DeltaS1 = abs((bMirror-bField_I(iFirst))*(dLength_I(iFirst-1))/(bField_I(iFirst-1)-bField_I(iFirst)))
    HalfPathLength= HalfPathLength + Coeff*2.*DeltaS1/(sqrt(bMirror-bField_I(iFirst)))

    do iStep = iFirst, iLast-1
       b1 = bField_I(iStep)
       b2 = bField_I(iStep+1)
       HalfPathLength = HalfPathLength + Coeff*2.*dLength_I(iStep)/(b1 - b2) &
            *( sqrt(bMirror  - b2) - sqrt(bMirror  - b1) )

    end do

    DeltaS2 = abs((bMirror-bField_I(iLast))*(dLength_I(iLast))/(bField_I(iLast+1)-bField_I(iLast)))
    HalfPathLength= HalfPathLength + Coeff*2.*DeltaS2/(sqrt(bMirror-bField_I(iLast)))

  end subroutine half_bounce_path_length

  !================================================================================== 

  real function dipole_length(L, LatMin, LatMax)

    implicit none

    real               :: L                       ! L shell value
    real               :: LatMin                  ! Minimum Latitude
    real               :: LatMax                  ! Maximum Latitude
    real               :: x, y
    !----------------------------------------------------------------------------------

    x = sin(LatMax)
    y = sin(LatMin)

    dipole_length = abs(L*(0.5*Sqrt(3.0*x**2.0+1.0)*x + asinh(sqrt(3.0)*x)/2.0/sqrt(3.0)) &
         - L*(0.5*Sqrt(3.0*y**2.0+1.0)*y + asinh(sqrt(3.0)*y)/2.0/sqrt(3.0))) 
    
  end function dipole_length

  !==================================================================================

  real function asinh(x)              

    implicit none

    real, intent(in) :: x
    !----------------------------------------------------------------------------------

    asinh = log(x+sqrt(x**2+1.0))

  end function asinh

  !==================================================================================

  real function second_adiab_invariant(x)		

    ! This is an approximate formula for second adiabatic invariant 
    ! = the integral of sqrt(1 - B/Bm)ds for a dipole B field. 
    !	function I(mu) taken from Ejiri, JGR,1978

    real, intent(in)     :: x ! cosine of the equatorial pitch angle
    real                 :: y, alpha, beta, a1, a2, a3, a4
    real, parameter      :: Pi = 3.141592654
    !----------------------------------------------------------------------------------

    y=sqrt(1.-x*x)
    alpha=1.+alog(2.+sqrt(3.))/2./sqrt(3.)
    beta=alpha/2.-Pi*sqrt(2.)/12.
    a1=0.055
    a2=-0.037
    a3=-0.074
    a4=0.056
    second_adiab_invariant = &
         2.*alpha*(1.-y)+2.*beta*y*alog(y)+4.*beta*(y-sqrt(y))+  &
         3.*a1*(y**(1./3.)-y)+6.*a2*(y**(2./3.)-y)+6.*a4*(y-y**(4./3.))  &
         -2.*a3*y*alog(y)

  end function second_adiab_invariant

  !==================================================================================

  real function  analytic_h(x)		

    ! This is an approximate formula for second adiabatic invariant 
    ! = the integral of ds/sqrt(1 - B/Bm) for a dipole B field. 
    !	function f(y) taken from Ejiri, JGR,1978

    real, intent(in)     :: x ! cosine of the equatorial pitch angle 
    real                 :: y, alpha, beta, a1, a2, a3, a4
    real, parameter      :: Pi = 3.141592654
    !----------------------------------------------------------------------------------

    y=sqrt(1-X*X)
    alpha=1.+alog(2.+sqrt(3.))/2./sqrt(3.)
    beta=alpha/2.-Pi*sqrt(2.)/12.
    a1=0.055
    a2=-0.037
    a3=-0.074
    a4=0.056
    analytic_h = alpha-beta*(y+sqrt(y))+a1*y**(1./3.)+a2*y**(2./3.)+  &
         a3*y+a4*y**(4./3.)


  end function analytic_h
  !==================================================================================

  subroutine test_general_b

    real                 :: L = 10.0
    integer, parameter   :: nStep = 10001
    integer, parameter   :: nPitch = 1
    real                 :: PitchAngle_I(nPitch) 
    real                 :: bField_I(nstep) 
    real                 :: Length_I(nStep) 
    real                 :: RadialDistance_I(nStep) 
    real                 :: Ds_I(nStep)
    real                 :: SecondAdiabInv, IntegralBAnalytic
    real                 :: HalfPathLength, IntegralHAnalytic
    real                 :: bMirror_I(nPitch),bMirror
    integer              :: iMirror_II(2,nPitch)
    real, parameter      :: Pi = 3.141592654
    real                 :: Percent1, Percent2
    integer              :: iStep
    !----------------------------------------------------------------------------------
    open (unit = 2, file = "Convergence_nonuniform.dat")
    write (2,*)'Numerical values for the second adiabatic invariant integration btw a mirror point and eq'
    write (2,*)'nStep  IntegralBAnalytic   2nd_orderB  IntegralHAnalytic  2nd_orderH Percent1, Percent2'

    do iStep = 101, nStep,100
       call initialize_b_field(L, iStep, bField_I, RadialDistance_I, Length_I, Ds_I)
       PitchAngle_I(1) = Pi/10.
       IntegralBAnalytic = second_adiab_invariant(cos(PitchAngle_I(1)))
       IntegralHAnalytic = analytic_h(cos(PitchAngle_I(1)))

       call find_mirror_points (iStep, nPitch, PitchAngle_I, bField_I, bMirror_I,iMirror_II)

       bMirror = bMirror_I(1)
!!$       
!!$       write(*,*) 'iStep', iStep
!!$       write(*,*) '============================================='
!!$       write(*,*) 'bMirror', bMirror
!!$       write(*,*) '============================================='
!!$       write(*,*) 'iMirror_II(:,1)',iMirror_II(:,1)
!!$       write(*,*) '============================================='
!!$       write(*,*) 'bField_I',bField_I
!!$       write(*,*) '============================================='
!!$       write(*,*) 'Ds_I',Ds_I
!!$       write(*,*) '============================================='
!!$       write(*,*) 'L',L
!!$       write(*,*) '============================================='
              
       call second_adiabatic_invariant(iStep, iMirror_II(:,1), bMirror, bField_I, Ds_I,L, SecondAdiabInv)
!!$       
!!$       write(*,*) '============================================='
!!$       write(*,*) 'SecondAdiabInv',SecondAdiabInv
!!$       write(*,*) '============================================='      
!!$       write(*,*) 'GOT HERE'
       
       Percent1 = 100*abs(IntegralBAnalytic - SecondAdiabInv)/IntegralBAnalytic

       call half_bounce_path_length(iStep, iMirror_II(:,1), bMirror, bField_I, Ds_I,L, HalfPathLength)


       Percent2 = 200*abs(IntegralHAnalytic - HalfPathLength)/(IntegralHAnalytic+HalfPathLength)

       write (2,*) iStep, IntegralBAnalytic, SecondAdiabInv,IntegralHAnalytic,&
            HalfPathLength ,Percent1, Percent2

!       write (*,*) iStep, IntegralBAnalytic, SecondAdiabInv,IntegralHAnalytic,&
!            HalfPathLength !,Percent1, Percent2

    end do

    close(2)

  end subroutine test_general_b
  !==================================================================================

end module ModHeidiBField

