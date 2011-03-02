module NeutralHydrogenModel

  implicit none
contains

  !===========================================================================================
  subroutine get_ostgaard_density(RadialDistance,  HDensity)
    !\
    ! Calculates the geocoronal density as provided by the Ostgaard et al. (2003) model
    !/
    use ModHeidiInput, ONLY: Ct, SolarZenithAngle

    real, intent(in)  :: RadialDistance
    real, intent(out) :: HDensity
    real              :: n1, n2, alpha1, alpha2
    real              :: Exponent1, Exponent2
    !-------------------------------------------------------------------------------------------    

    call get_ostgaard_density_parameters(n1, n2, alpha1, alpha2)

    Exponent1 = -RadialDistance/alpha1
    Exponent2 = -RadialDistance/alpha2

    HDensity = Ct * (n1 * exp(Exponent1) + n2 * exp(Exponent2))

  end subroutine get_ostgaard_density
  !===========================================================================================
  subroutine get_ostgaard_density_parameters(n1, n2, alpha1, alpha2)

    use ModNumConst,   ONLY: cDegToRad
    use ModHeidiInput, ONLY: SolarZenithAngle

    real, intent(out) :: n1, n2, alpha1, alpha2
    !-------------------------------------------------------------------------------------------

    if (SolarZenithAngle == 90.0) then
       n1 = 10000.0  ! in cm ^-3
       n2 = 70.0     ! in cm ^-3
       alpha1 = 1.02 ! in Re
       alpha2 = 8.2  ! in Re
    end if

    if (SolarZenithAngle == 100.0) then
       n1 = 10100.0  ! in cm ^-3
       n2 = 80.0     ! in cm ^-3
       alpha1 = 1.01 ! in Re
       alpha2 = 7.9  ! in Re
    end if

    if (SolarZenithAngle == 110.0) then
       n1 = 10300.0  ! in cm ^-3
       n2 = 100.0    ! in cm ^-3
       alpha1 = 0.99 ! in Re
       alpha2 = 7.1  ! in Re
    end if

    if (SolarZenithAngle == 120.0) then
       n1 = 10600.0  ! in cm ^-3
       n2 = 130.0    ! in cm ^-3
       alpha1 = 0.96 ! in Re
       alpha2 = 6.3  ! in Re
    end if

    if (SolarZenithAngle == 130.0) then
       n1 = 10900.0  ! in cm ^-3
       n2 = 180.0    ! in cm ^-3
       alpha1 = 0.93 ! in Re
       alpha2 = 5.7  ! in Re
    end if

    if (SolarZenithAngle == 140.0) then
       n1 = 11300.0  ! in cm ^-3
       n2 = 220.0    ! in cm ^-3
       alpha1 = 0.90 ! in Re
       alpha2 = 5.2  ! in Re
    end if

    if (SolarZenithAngle == 150.0) then
       n1 = 11600.0  ! in cm ^-3
       n2 = 250.0    ! in cm ^-3
       alpha1 = 0.88 ! in Re
       alpha2 = 4.9  ! in Re
    end if

    if (SolarZenithAngle == 160.0) then
       n1 = 11800.0  ! in cm ^-3
       n2 = 280.0    ! in cm ^-3
       alpha1 = 0.86 ! in Re
       alpha2 = 4.8  ! in Re
    end if

    if (SolarZenithAngle == 170.0) then
       n1 = 12000.0  ! in cm ^-3
       n2 = 300.0    ! in cm ^-3
       alpha1 = 0.85 ! in Re
       alpha2 = 4.7  ! in Re
    end if

    if (SolarZenithAngle == 180.0) then
       n1 = 12000.0  ! in cm ^-3
       n2 = 310.0    ! in cm ^-3
       alpha1 = 0.85 ! in Re
       alpha2 = 4.6  ! in Re
    end if


  end subroutine get_ostgaard_density_parameters
  !===========================================================================================
  subroutine get_hodges_density(Theta, Phi, HDensity, R_I)
    !\
    ! Calculates the geocoronal density as provided by the Hodges et al. (1994) model
    !/
    use ModNumConst,  ONLY: cPi

    integer, parameter   :: nR=40 ! Number of radial distances in the Hodges model  
    real, intent(in)     :: Phi, Theta
    real,    intent(out) :: HDensity(nR), R_I(nR)
    real                 :: N_I(nR), A_lm(nR), B_lm(nR), Y_lm
    real                 :: Z
    integer              :: l, m, iR
    !-------------------------------------------------------------------------------------------

    HDensity = 0.0
    B_lm = 0.0

    do iR = 1, nR
       Z = 0.0
       do l = 0, 3
          do m = 0, l
             call get_Ylm(Theta, Phi, l, m, Y_lm)
             call get_AlmBlm(l, m, A_lm, B_lm, R_I, N_I)
             Z = Z + (A_lm(iR) * cos(m * Phi) + B_lm(iR) * sin(m* Phi)) * Y_lm
          end do
       end do
       HDensity(iR) = N_I(iR) * sqrt(4. * cPi) * Z
    end do
  
  end subroutine get_hodges_density
  !===========================================================================================
  subroutine get_Ylm(Theta, Phi, l, m, Y_lm)
    
    use ModNumConst, ONLY: cPi

    real, intent(in)    :: Theta, Phi
    integer, intent(in) :: l, m
    real, intent(out)   :: Y_lm

    !Local variables
    real :: Y_00, Y_10, Y_11, Y_20, Y_21, Y_22, Y_30, Y_31, Y_32, Y_33
    real :: cos2Theta, cos3Theta, sin2Theta, sin3Theta
    !-------------------------------------------------------------------------------------------

    cos2Theta = (cos(Theta))**2
    sin2Theta = (sin(Theta))**2
    cos3Theta = (cos(Theta))**3  
    sin3Theta = (sin(Theta))**3

    Y_00 =   sqrt(1./(4.*cPi))

    Y_10 =   sqrt(3./(4.*cPi)) * cos(Theta)
    Y_11 = - sqrt(3./(8.*cPi)) * sin(Theta)
    
    Y_20 =   sqrt(5./(4.*cPi)) * (3./2. * cos2Theta - 0.5)
    Y_21 = - sqrt(15./(8.*cPi)) * sin(Theta) * cos(Theta)
    Y_22 = 1./4. * sqrt(15./(2.*cPi)) * sin2Theta
    
    Y_30 =           sqrt(7./(4.*cPi)) * (5./2. * cos3Theta - 3./2. * cos(Theta))
    Y_31 = - 1./4. * sqrt(21./(4.*cPi)) * sin(Theta) * (5. * cos2Theta -1.)
    Y_32 =   1./4. * sqrt(105./(2.*cPi)) * sin2Theta * cos(Theta)
    Y_33 = - 1./4. * sqrt(35./(4.*cPi)) * sin3Theta
    
    if ((l==0) .and. (m==0)) Y_lm = Y_00
    if ((l==1) .and. (m==0)) Y_lm = Y_10
    if ((l==1) .and. (m==1)) Y_lm = Y_11

    if ((l==2) .and. (m==0)) Y_lm = Y_20
    if ((l==2) .and. (m==1)) Y_lm = Y_21
    if ((l==2) .and. (m==2)) Y_lm = Y_22

    if ((l==3) .and. (m==0)) Y_lm = Y_30
    if ((l==3) .and. (m==1)) Y_lm = Y_31
    if ((l==3) .and. (m==2)) Y_lm = Y_32
    if ((l==3) .and. (m==3)) Y_lm = Y_33

  end subroutine get_Ylm
  !===========================================================================================
  subroutine get_AlmBlm(l, m, A_lm, B_lm, R_I, N_I)

    use ModIoUnit,     ONLY : UnitTmp_
    use ModHeidiIO,    ONLY : NameInputDirectory
    use ModHeidiInput, ONLY : TypeHModel, TypeSeason, WhichF107

    
    integer,            intent(in)  :: l, m
    real, dimension(n), intent(out) :: A_lm, B_lm 
    real, dimension(n), intent(out) :: R_I, N_I  
    real, dimension(n)              :: A_00, A_10, A_11, A_20, A_21, &
         A_22, A_30, A_31, A_32, A_33
    real, dimension(n)              :: B_00, B_10, B_11, B_20, B_21, &
         B_22, B_30, B_31, B_32, B_33
    integer, parameter   :: n=40 ! Number of radial distances in the Hodges model  
    character (len=100)  :: StringHeader
    integer              :: i
    !-------------------------------------------------------------------------------------------

     open(UNITTMP_,FILE=NameInputDirectory//trim(TypeHModel)//'_'//trim(TypeSeason)//'F107_'//trim(WhichF107)//'.dat',status='old')
     !\
     ! According to Hodges 1994, the A and B coefficients have been scaled by 10^4 in the provided tables. 
     !/

     A_00 = 1.0e4
     
     do i = 1, 2
        read (UnitTmp_,*) StringHeader
     end do
     
     do i = 1, n
        read(UnitTmp_,*) R_I(i), N_I(i), A_10(i), A_11(i), B_11(i), A_20(i), A_21(i), B_21(i),&
             A_22(i), B_22(i), A_30(i), A_31(i), B_31(i),  A_32(i), B_32(i), A_33(i), B_33(i)
     end do

    close(UnitTmp_)

    if ((l==0) .and. (m==0)) then
       A_lm = A_00
       B_lm = B_00
    end if

    if ((l==1) .and. (m==0)) then
       A_lm = A_10
       B_lm = B_10
    end if
    if ((l==1) .and. (m==1)) then
       A_lm = A_11
       B_lm = B_11
    end if

    if ((l==2) .and. (m==0)) then
       A_lm = A_20
       B_lm = B_20
    end if
    if ((l==2) .and. (m==1)) then
       A_lm = A_21
       B_lm = B_21
    end if
    if ((l==2) .and. (m==2)) then
       A_lm = A_22
       B_lm = B_22
    end if

    if ((l==3) .and. (m==0)) then
       A_lm = A_30
       B_lm = B_30
    end if
    if ((l==3) .and. (m==1)) then
       A_lm = A_31
       B_lm = B_31
    end if
    if ((l==3) .and. (m==2)) then
       A_lm = A_32
       B_lm = B_32
    end if
    if ((l==3) .and. (m==3)) then
       A_lm = A_33
       B_lm = B_33
    end if

    A_lm = A_lm * 1.e-4
    B_lm = B_lm * 1.e-4

  end subroutine get_AlmBlm
  !===========================================================================================
end module NeutralHydrogenModel
