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
    use ModHeidiSize,  ONLY : nRHodges

    
    integer,            intent(in)  :: l, m
    real, dimension(nRHodges), intent(out) :: A_lm, B_lm 
    real, dimension(nRHodges), intent(out) :: R_I, N_I  
    real, dimension(nRHodges)              :: A_00, A_10, A_11, A_20, A_21, &
         A_22, A_30, A_31, A_32, A_33
    real, dimension(nRHodges)              :: B_00, B_10, B_11, B_20, B_21, &
         B_22, B_30, B_31, B_32, B_33
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
     
     do i = 1, nRHodges
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
  subroutine get_interpolated_hodge_density(RhoH_III)

    use ModHeidiInput,  ONLY: StretchingFactorA,StretchingFactorB, &
         TypeHModel, TypeSeason, WhichF107
    use ModHeidiSize,   ONLY: nRHodges, nPoint, nT, nR
    use ModHeidiMain,   ONLY: Phi, LZ, Re
    use ModNumConst,    ONLY: cPi
    use ModHeidiBField, ONLY: get_cubic_root
    use ModIoUnit,      ONLY: UnitTmp_
    use ModHeidiIO,     ONLY: NameInputDirectory
    use ModHeidiBField

    integer, parameter :: nUniform = 100
    real, intent(out) :: RhoH_III(nPoint,nR,nT)
    real, dimension(nRHodges)              :: A_00, A_10, A_11, A_20, A_21, &
         A_22, A_30, A_31, A_32, A_33
    real, dimension(nRHodges)              :: B_00, B_10, B_11, B_20, B_21, &
         B_22, B_30, B_31, B_32, B_33
    real :: HDensity(nRHodges), R_I(nRHodges), N_I(nRHodges)
    real :: HodgesDensity_III(nPoint, nRHodges,nT)
    real :: Theta_III(nPoint, nRHodges, nT)

    !Local Variables
    real    :: LatMin, LatMax, dLat, Lat, a, b    
    real    :: ThetaMax, ThetaMin, Theta, dTheta
    real    :: aa, bb, cc, dd, gamma,sigma, alpha
    real    :: cos2Lat1, CosPhi2, SinPhi2
    real    :: A_lm(nRHodges), B_lm(nRHodges)
    real    :: xc, yc, zc
    real    :: Z, Y_lm
    complex :: root(3)
    integer :: i, iPoint, iR, iPhi, nroot
    integer :: l,m
    real    :: rMin, rMax, dR
    real    :: LnHodgesDensity_III(nPoint, nRHodges, nT)
    real    :: RhoHUniform_III(nPoint, nUniform, nT)
    real    :: LnRhoH
    real    :: Rad_I(nUniform), Weight, r
    real    :: cos2Lat,sin2Lat,cos2Phi,sin2Phi
    character (len=100)  :: StringHeader
    integer :: iRHodges
    real                 :: bFieldMagnitude_III(nPoint,nR,nT) ! Magnitude of magnetic field 
    real                 :: RadialDistance_III(nPoint,nR,nT)
    real                 :: GradBCrossB_VIII(3,nPoint,nR,nT)
    real                 :: GradB_VIII(3,nPoint,nR,nT)
    real                 :: dLength_III(nPoint-1,nR,nT)      ! Length interval between i and i+1  
    real                 :: Length_III(nPoint,nR,nT) 
    real                 :: dBdt_III(nPoint,nR,nT)

    !-------------------------------------------------------------------------------------------
    
    !\
    ! Read the file containing the spherical harmonics expansion coefficients from Hodges 1994 tables.
    !/ 
    open(UNITTMP_,FILE=NameInputDirectory//trim(TypeHModel)//'_'//trim(TypeSeason)//'F107_'//trim(WhichF107)//'.dat',status='old')
    !\
    ! According to Hodges 1994, the A and B coefficients have been scaled by 10^4 in the provided tables. 
    !/
    
    A_00 = 1.0e4
    RhoHUniform_III = 0.0
    RhoH_III = 0.0

    do i = 1, 2
       read (UnitTmp_,*) StringHeader
    end do
    
    do i = 1, nRHodges
       read(UnitTmp_,*) R_I(i), N_I(i), A_10(i), A_11(i), B_11(i), A_20(i), A_21(i), B_21(i),&
            A_22(i), B_22(i), A_30(i), A_31(i), B_31(i),  A_32(i), B_32(i), A_33(i), B_33(i)
    end do
    
    close(UnitTmp_)
    
    
    !\
    ! Calculate the Latitude and colatitude needed for the calculation 
    ! of the neutral hydrogen densities.
    !/
    a = StretchingFactorA
    b = StretchingFactorB
    

    !\
    ! Need the Radial Distance to any point along a field line
    !/

    call initialize_b_field(LZ, Phi, nPoint, nR, nT, bFieldMagnitude_III, &
       RadialDistance_III,Length_III, dLength_III,GradBCrossB_VIII,GradB_VIII,dBdt_III)

    dd = 0.0
    do iPhi = 1, nT
       sigma   = cos(Phi(iPhi))
       CosPhi2 = (cos(Phi(iPhi)))**2
       SinPhi2 = (sin(Phi(iPhi)))**2
       alpha = a + b * cos(Phi(iPhi)) 
       gamma = alpha**2 * CosPhi2 + SinPhi2 
       aa = gamma - 1.0
       bb = 1.0
       
       do iR = 1, nRHodges
         ! cc = -1./(LZ(iR)*LZ(iR))
          cc = -1./(R_I(iR) * 1.e3/Re * R_I(iR) * 1.e3/Re)

          call get_cubic_root(aa,bb,dd,cc,root,nroot)
          do i =1, nroot
             if ((aimag(root(i)) <= 1.e-5).and. (real(root(i))<=1.0) .and. (real(root(i))>=0.0)) &
                  cos2Lat1 = real(root(i))
          end do
          
          LatMax = acos(sqrt(cos2Lat1))
          ThetaMin = cPi/2. - LatMax
          ThetaMax = cPi - LatMax
          dTheta = (ThetaMax-ThetaMin)/(nPoint-1)
          Theta = ThetaMax
          
          do iPoint = 1, nPoint
             
             cos2Lat = (cos(Lat))**2
             sin2Lat = (sin(Lat))**2
             cos2Phi = (cos(Phi(iPhi)))**2
             sin2Phi = (sin(Phi(iPhi)))**2   
             
             Theta_III(iPoint,iR,iPhi) =  Theta
             Theta = Theta - dTheta

          end do
       end do
    end do
  
  
    !\
    ! Start calculation the geocoronal hydrogen density according to Hodges model
    !/
    
    HDensity = 0.0
    B_lm = 0.0
    
    ! need to do more work.......check array sizes... nR vs nRHodges
    do iPhi = 1, nT
       do iPoint = 1, nPoint
          do iR = 1, nRHodges
             Z = 0.0
             do l = 0, 3
                do m = 0, l
                   call get_Ylm(Theta_III(iPoint,iR,iPhi), Phi(iPhi), l, m, Y_lm)
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
                    
                    
                    Z = Z + (A_lm(iR) * cos(m * Phi(iPhi)) + B_lm(iR) * sin(m* Phi(iPhi))) * Y_lm
                 end do
              end do
              
              HDensity(iR) = N_I(iR) * sqrt(4. * cPi) * Z
              HodgesDensity_III(iPoint,iR, iPhi) = HDensity(iR)
             
              LnHodgesDensity_III(iPoint,iR, iPhi) = log(HDensity(iR)*10.**6)  ! need the density in m^-3 NOT cm^-3!!!!

           end do
        end do
     end do
     

     !\
     !Interpolate the hydrogen density to a new, well refined grid.
     !/
     
     rMax = R_I(nRHodges)/Re * 1.e3   ! in Re
     rMin = R_I(1)/Re * 1.e3   
     dR = (rMax-rMin)/(nUniform-1)
     
     do i = 1, nUniform
        Rad_I(i) = rMin + dR*(i-1)
        do iPhi = 1, nT
           do iPoint = 1, nPoint
              do iRHodges = 1, nRHodges
                 if (R_I(iRHodges)/Re * 1.e3>=Rad_I(i)) then
                    iR = iRHodges-1
                    
                    LnRhoH = LnHodgesDensity_III(iPoint,iR,iPhi)+ (Rad_I(i)-R_I(iR)/Re * 1.e3)*&
                         (LnHodgesDensity_III(iPoint,iR+1, iPhi) - LnHodgesDensity_III(iPoint,iR,iPhi))/&
                         (R_I(iR+1)/Re * 1.e3 - R_I(iR)/Re * 1.e3)
                    
                    EXIT
                 end if
                 RhoHUniform_III(iPoint,i,iPhi) = exp(LnRhoH)
                 
              end do
           end do
        end do
     end do
          
     do iPhi =1, nT
       do iR = 1, nR
          do iPoint =1, nPoint
             
             r = RadialDistance_III(iPoint, iR, iPhi)
             if (r>=rmax) then
                RhoH_III(iPoint,iR, iPhi) = RhoHUniform_III(iPoint,nUniform,iPhi)
             else
                i = 1 + abs((r - rMin))/dR
                Weight = (r-Rad_I(i))/dR
                RhoH_III(iPoint,iR,iPhi) = Weight*RhoHUniform_III(iPoint,i+1,iPhi) +&
                     (1. - Weight) * RhoHUniform_III(iPoint,i,iPhi)

                if (RhoH_III(iPoint,iR,iPhi) <=0.0) RhoH_III(iPoint,iR,iPhi) = 0.0
                
             end if
          end do
          
       end do
    end do
 


       write(*,*) 'Rho', RhoH_III(48,13,12)
       STOP

  end subroutine get_interpolated_hodge_density

  !===========================================================================================
 subroutine get_bailey_density(RhoH_III)

    use ModHeidiInput,  ONLY: StretchingFactorA,StretchingFactorB, &
         TypeHModel, TypeSeason, WhichF107
    use ModHeidiSize,   ONLY: nPoint, nT, nR
    use ModHeidiMain,   ONLY: Phi, LZ, Re
    use ModNumConst,    ONLY: cPi
    use ModHeidiBField, ONLY: get_cubic_root
    use ModHeidiBField
    
  
    real, intent(out) :: RhoH_III(nPoint,nR,nT)
    
    real :: A_lm, B_lm, Y_lm, Z
    real :: Theta_III(nPoint, nR, nT)
    real :: LatMin, LatMax, dLat, Lat, a, b    
    real :: ThetaMax, ThetaMin, Theta, dTheta
    real :: aa, bb, cc, dd, gamma,sigma, alpha
    real :: cos2Lat1, CosPhi2, SinPhi2, cos2Lat, sin2Lat, cos2Phi, sin2Phi
    real :: bFieldMagnitude_III(nPoint,nR,nT) ! Magnitude of magnetic field 
    real :: RadialDistance_III(nPoint,nR,nT)
    real :: GradBCrossB_VIII(3,nPoint,nR,nT)
    real :: GradB_VIII(3,nPoint,nR,nT)
    real :: dLength_III(nPoint-1,nR,nT)      ! Length interval between i and i+1  
    real :: Length_III(nPoint,nR,nT) 
    real :: dBdt_III(nPoint,nR,nT)
    real :: n_III(nPoint,nR,nT)
    complex :: root(3)
    integer :: i, iPoint, iR, iPhi, nroot, l, m
    character (len=100)  :: StringHeader
    
    
    real, parameter :: A_00 = 1.00
    real, parameter :: B_00 = 0.00
    real, parameter :: C_00 = 0.00
    real, parameter :: D_00 = 0.00
    real, parameter :: D_10 = 0.00
    real, parameter :: D_20 = 0.00
    real, parameter :: C_10 = 0.00
    real, parameter :: C_20 = 0.00
    
    real, parameter :: A_10 = -4.8992e-2
    real, parameter :: B_10 = -1.8720e-6
    real, parameter :: A_11 = -3.8248e-1
    real, parameter :: B_11 =  9.0636e-6
    real, parameter :: C_11 = -4.8547e-2
    real, parameter :: D_11 = -2.1587e-6
    real, parameter :: A_20 =  1.5739e-1
    real, parameter :: B_20 = -6.1959e-6
    real, parameter :: A_21 = -6.9198e-2 
    real, parameter :: B_21 =  4.5477e-6
    real, parameter :: C_21 =  2.1922e-1
    real, parameter :: D_21 = -7.0881e-6
    real, parameter :: A_22 = -1.0148e-1
    real, parameter :: B_22 =  1.4873e-6 
    real, parameter :: C_22 = -8.8242e-2
    real, parameter :: D_22 =  4.2384e-6

    real, parameter :: p    = 4.1118e13
    real, parameter :: k    = -2.5446
    !-------------------------------------------------------------------------------------------
    
    RhoH_III = 0.0

    !\
    ! Calculate the Latitude and colatitude needed for the calculation 
    ! of the neutral hydrogen densities.
    !/
    a = StretchingFactorA
    b = StretchingFactorB
    
    !\
    ! Need the Radial Distance to any point along a field line
    !/

    call initialize_b_field(LZ, Phi, nPoint, nR, nT, bFieldMagnitude_III, &
       RadialDistance_III,Length_III, dLength_III,GradBCrossB_VIII,GradB_VIII,dBdt_III)

    dd = 0.0
    do iPhi = 1, nT
       sigma   = cos(Phi(iPhi))
       CosPhi2 = (cos(Phi(iPhi)))**2
       SinPhi2 = (sin(Phi(iPhi)))**2
       alpha = a + b * cos(Phi(iPhi)) 
       gamma = alpha**2 * CosPhi2 + SinPhi2 
       aa = gamma - 1.0
       bb = 1.0
       do iR = 1, nR
          cc = -1./(LZ(iR)*LZ(iR))
          call get_cubic_root(aa,bb,dd,cc,root,nroot)
          do i =1, nroot
             if ((aimag(root(i)) <= 1.e-5).and. (real(root(i))<=1.0) .and. (real(root(i))>=0.0)) &
                  cos2Lat1 = real(root(i))
          end do
          
          LatMax = acos(sqrt(cos2Lat1))
          ThetaMin = cPi/2. - LatMax
          ThetaMax = cPi - LatMax
          dTheta = (ThetaMax-ThetaMin)/(nPoint-1)
          Theta = ThetaMax
          
          do iPoint = 1, nPoint
             cos2Lat = (cos(Lat))**2
             sin2Lat = (sin(Lat))**2
             cos2Phi = (cos(Phi(iPhi)))**2
             sin2Phi = (sin(Phi(iPhi)))**2   
             Theta_III(iPoint,iR,iPhi) =  Theta
             Theta = Theta - dTheta
          end do
       end do
    end do
    
    !\
    ! Start calculation the geocoronal hydrogen density according to Bailey 2011 model.
    ! Valid only for the June 11, 2008 event.
    !/
    
    RadialDistance_III = RadialDistance_III * Re/1000. ! need radial distance in km

    do iPhi = 1, nT
       do iR = 1, nR
          do iPoint = 1, nPoint
            

             Z = 0.0
             do l = 0, 2
                do m = 0, l
                   call get_Ylm(Theta_III(iPoint,iR,iPhi), Phi(iPhi), l, m, Y_lm)
                   
                   if ((l==0) .and. (m==0)) then
                       A_lm = A_00 + B_00 * RadialDistance_III(iPoint, iR, iPhi)
                       B_lm = C_00 + D_00 * RadialDistance_III(iPoint, iR, iPhi)
                    end if
                    
                    if ((l==1) .and. (m==0)) then
                       A_lm = A_10 + B_10 * RadialDistance_III(iPoint, iR, iPhi)
                       B_lm = C_10 + D_10 * RadialDistance_III(iPoint, iR, iPhi)
                    end if
                    if ((l==1) .and. (m==1)) then
                       A_lm = A_11 + B_11 * RadialDistance_III(iPoint, iR, iPhi)
                       B_lm = C_11 + D_11 * RadialDistance_III(iPoint, iR, iPhi)
                    end if
                    
                    if ((l==2) .and. (m==0)) then
                       A_lm = A_20 + B_20 * RadialDistance_III(iPoint, iR, iPhi)
                       B_lm = C_20 + D_20 * RadialDistance_III(iPoint, iR, iPhi)
                    end if
                    if ((l==2) .and. (m==1)) then
                       A_lm = A_21 + B_21 * RadialDistance_III(iPoint, iR, iPhi)
                       B_lm = C_21 + D_21 * RadialDistance_III(iPoint, iR, iPhi)
                    end if
                    if ((l==2) .and. (m==2)) then
                       A_lm = A_22 + B_22 * RadialDistance_III(iPoint, iR, iPhi)
                       B_lm = C_22 + D_22 * RadialDistance_III(iPoint, iR, iPhi)
                    end if
                    
              
                    Z = Z + (A_lm * cos(m * Phi(iPhi)) + B_lm * sin(m* Phi(iPhi))) * Y_lm
                    
                 end do
              end do
              
              n_III(iPoint, iR, iPhi) = p * (RadialDistance_III(iPoint, iR, iPhi)) ** k * 1.e4  ! need the density in m^-3 NOT cm^-3!!!!
              RhoH_III(iPoint,iR,iPhi) = sqrt(4. * cPi) * n_III(iPoint, iR, iPhi) * Z  

           end do
        end do
     end do
     
     


   end subroutine get_bailey_density

  !===========================================================================================

end module NeutralHydrogenModel
