! Copyright (C) 2002 Regents of the University of Michigan, 
! portions used with permission 
! For more information, see http://csem.engin.umich.edu/tools/swmf
!======================================
!                                     |
!    Module for Ionosphere Model      |
!                                     |
!======================================

module ModIonosphere
  use ModNumConst
  use IE_ModSize

  implicit none
  save

  !\
  ! Ionosphere solution parameters
  !/
  real, parameter ::                            &
       IONO_TOLER = 5.0e-05,                    &
       IONO_MU = 1.256637e-06,                  &
       IONO_Theta_0 = 0.0001,                   &
       IONO_Min_EFlux = 0.1e-16,                &     ! W/m2
       IONO_Min_Ave_E = 0.5,                    &     ! keV
       Polar_Rain = 0.1e-2                            ! W/m2

  integer, parameter :: IONO_Model_No_Hall = 1, &
       IONO_Model_With_Hall = 2,                &
       IONO_Model_With_Simple_Aurora = 3,       &
       IONO_Model_With_Complex_Aurora = 4

  real :: IONO_Bdp,   &
       IONO_Radius=1.0, IONO_Height=1.0, Radius

  real :: cpcp_north=0.0, cpcp_south=0.0

  ! Variables for empirical conductance (iModels 4 & 5):
  ! File names for coefficients:
  character(len=100) :: &
       NameHalFile = 'cond_hal_coeffs.dat', &
       NamePedFile = 'cond_ped_coeffs.dat'

  logical :: UseCMEEFitting  = .false.
  real :: LatNoConductanceSI = 45.0
  real :: FactorHallCMEE = 7.5, FactorPedCMEE = 5.0

  ! Coefficients for conductance based on FAC:
  real, allocatable, dimension(:,:) ::  &
       hal_a0_up,ped_a0_up,      &
       hal_a0_do,ped_a0_do,      &
       hal_a1_up,ped_a1_up,      &
       hal_a1_do,ped_a1_do,      &
       hal_a2_up,ped_a2_up,      &
       hal_a2_do,ped_a2_do

  ! Grid for conductance coefficients:
  integer :: i_cond_nmlts=-1, i_cond_nlats=-1
  real, allocatable :: cond_mlts(:)
  real, allocatable :: cond_lats(:)


  !\
  ! Ionosphere Solution on the whole grid
  !/
  real, allocatable :: IONO_Phi(:,:)
  real, allocatable :: IONO_IonNumFlux(:,:)
  real, allocatable :: IONO_Joule(:,:)
  real, allocatable :: IONO_Jr(:,:)
  real, allocatable :: IONO_Ave_E(:,:)
  real, allocatable :: IONO_Eflux(:,:)
  real, allocatable :: IONO_SigmaP(:,:)
  real, allocatable :: IONO_SigmaH(:,:)

  !\
  ! Ionosphere solution array definitions
  !/
  real, allocatable :: IONO_NORTH_Phi(:,:)
  real, allocatable :: IONO_SOUTH_Phi(:,:)
  real, allocatable :: IONO_NORTH_X(:,:)
  real, allocatable :: IONO_NORTH_Y(:,:)
  real, allocatable :: IONO_NORTH_Z(:,:)
  real, allocatable :: IONO_NORTH_Theta(:,:)
  real, allocatable :: IONO_NORTH_Psi(:,:)
  real, allocatable :: IONO_SOUTH_X(:,:)
  real, allocatable :: IONO_SOUTH_Y(:,:)
  real, allocatable :: IONO_SOUTH_Z(:,:)
  real, allocatable :: IONO_SOUTH_Theta(:,:)
  real, allocatable :: IONO_SOUTH_Psi(:,:)
  real, allocatable :: IONO_NORTH_Ex(:,:)
  real, allocatable :: IONO_NORTH_Ey(:,:)
  real, allocatable :: IONO_NORTH_Ez(:,:)
  real, allocatable :: IONO_NORTH_ETh(:,:)
  real, allocatable :: IONO_NORTH_EPs(:,:)
  real, allocatable :: IONO_SOUTH_Ex(:,:)
  real, allocatable :: IONO_SOUTH_Ey(:,:)
  real, allocatable :: IONO_SOUTH_Ez(:,:)
  real, allocatable :: IONO_SOUTH_ETh(:,:)
  real, allocatable :: IONO_SOUTH_EPs(:,:)
  real, allocatable :: IONO_NORTH_Ux(:,:)
  real, allocatable :: IONO_NORTH_Uy(:,:)
  real, allocatable :: IONO_NORTH_Uz(:,:)
  real, allocatable :: IONO_NORTH_UTh(:,:)
  real, allocatable :: IONO_NORTH_UPs(:,:)
  real, allocatable :: IONO_SOUTH_Ux(:,:)
  real, allocatable :: IONO_SOUTH_Uy(:,:)
  real, allocatable :: IONO_SOUTH_Uz(:,:)
  real, allocatable :: IONO_SOUTH_UTh(:,:)
  real, allocatable :: IONO_SOUTH_UPs(:,:)
  real, allocatable :: IONO_NORTH_EFlux(:,:)
  real, allocatable :: IONO_NORTH_Ave_E(:,:)
  real, allocatable :: IONO_SOUTH_EFlux(:,:)
  real, allocatable :: IONO_SOUTH_Ave_E(:,:)
  real, allocatable :: IONO_NORTH_Sigma0(:,:)
  real, allocatable :: IONO_NORTH_SigmaH(:,:)
  real, allocatable :: IONO_NORTH_SigmaP(:,:)
  real, allocatable :: IONO_NORTH_SigmaThTh(:,:)
  real, allocatable :: IONO_NORTH_SigmaThPs(:,:)
  real, allocatable :: IONO_NORTH_SigmaPsPs(:,:)
  real, allocatable :: IONO_SOUTH_Sigma0(:,:)
  real, allocatable :: IONO_SOUTH_SigmaH(:,:)
  real, allocatable :: IONO_SOUTH_SigmaP(:,:)
  real, allocatable :: IONO_SOUTH_SigmaThTh(:,:)
  real, allocatable :: IONO_SOUTH_SigmaThPs(:,:)
  real, allocatable :: IONO_SOUTH_SigmaPsPs(:,:)
  real, allocatable :: IONO_NORTH_dSigmaThTh_dTheta(:,:)
  real, allocatable :: IONO_NORTH_dSigmaThPs_dTheta(:,:)
  real, allocatable :: IONO_NORTH_dSigmaPsPs_dTheta(:,:)
  real, allocatable :: IONO_NORTH_dSigmaThTh_dPsi(:,:)
  real, allocatable :: IONO_NORTH_dSigmaThPs_dPsi(:,:)
  real, allocatable :: IONO_NORTH_dSigmaPsPs_dPsi(:,:)
  real, allocatable :: IONO_SOUTH_dSigmaThTh_dTheta(:,:)
  real, allocatable :: IONO_SOUTH_dSigmaThPs_dTheta(:,:)
  real, allocatable :: IONO_SOUTH_dSigmaPsPs_dTheta(:,:)
  real, allocatable :: IONO_SOUTH_dSigmaThTh_dPsi(:,:)
  real, allocatable :: IONO_SOUTH_dSigmaThPs_dPsi(:,:)
  real, allocatable :: IONO_SOUTH_dSigmaPsPs_dPsi(:,:)
  real, allocatable :: SAVE_NORTH_SigmaH(:,:)
  real, allocatable :: SAVE_NORTH_SigmaP(:,:)
  real, allocatable :: SAVE_SOUTH_SigmaH(:,:)
  real, allocatable :: SAVE_SOUTH_SigmaP(:,:)
  real, allocatable :: IONO_NORTH_Joule(:,:)
  real, allocatable :: IONO_SOUTH_Joule(:,:)
  real, allocatable :: IONO_NORTH_IonNumFlux(:,:)
  real, allocatable :: IONO_SOUTH_IonNumFlux(:,:)
  real, allocatable :: IONO_NORTH_JR(:,:)
  real, allocatable :: IONO_NORTH_JTh(:,:)
  real, allocatable :: IONO_NORTH_JPs(:,:)
  real, allocatable :: IONO_NORTH_Jx(:,:)
  real, allocatable :: IONO_NORTH_Jy(:,:)
  real, allocatable :: IONO_NORTH_Jz(:,:)
  real, allocatable :: IONO_SOUTH_JR(:,:)
  real, allocatable :: IONO_SOUTH_JTh(:,:)
  real, allocatable :: IONO_SOUTH_JPs(:,:)
  real, allocatable :: IONO_SOUTH_Jx(:,:)
  real, allocatable :: IONO_SOUTH_Jy(:,:)
  real, allocatable :: IONO_SOUTH_Jz(:,:)
  real, allocatable :: IONO_NORTH_TGCM_JR(:,:)
  real, allocatable :: IONO_SOUTH_TGCM_JR(:,:)
  real, allocatable :: IONO_NORTH_Fake_JR(:,:)
  real, allocatable :: IONO_SOUTH_Fake_JR(:,:)
  real, allocatable :: iono_north_im_jr(:,:)
  real, allocatable :: iono_south_im_jr(:,:)
  real, allocatable :: iono_north_im_avee(:,:)
  real, allocatable :: iono_south_im_avee(:,:)
  real, allocatable :: iono_north_im_eflux(:,:)
  real, allocatable :: iono_south_im_eflux(:,:)

  logical, allocatable :: IsFilledWithIm(:,:)

  real, allocatable :: IONO_NORTH_invB(:,:)
  real, allocatable :: IONO_SOUTH_invB(:,:)
  real, allocatable :: IONO_NORTH_rho(:,:)
  real, allocatable :: IONO_SOUTH_rho(:,:)
  real, allocatable :: IONO_NORTH_p(:,:)
  real, allocatable :: IONO_SOUTH_p(:,:)
  real, allocatable :: IONO_NORTH_t(:,:)
  real, allocatable :: IONO_SOUTH_t(:,:)
  real, allocatable :: IONO_NORTH_dLat(:,:)
  real, allocatable :: IONO_SOUTH_dLat(:,:)
  real, allocatable :: IONO_NORTH_dLon(:,:)
  real, allocatable :: IONO_SOUTH_dLon(:,:)

  ! Pentadiagonal matrix for the Poisson equation
  real, allocatable :: C_A(:,:)
  real, allocatable :: C_B(:,:)
  real, allocatable :: C_C(:,:)
  real, allocatable :: C_D(:,:)
  real, allocatable :: C_E(:,:)
  real, dimension(:), allocatable :: d_I, e_I, f_I, e1_I, f1_I
  logical :: north, DoPrecond
  integer :: nThetaUsed, nX

  real, dimension(IONO_nTheta) :: dTheta_North, dTheta_South
  real, dimension(IONO_nPsi)   :: dPsi_North, dPsi_South

contains
  !===========================================================================
  subroutine load_conductances()

    use ModIoUnit, ONLY: UnitTmp_
    use ModUtilities, ONLY: open_file, close_file

    ! Local variables:
    character (len=100) :: Line
    integer :: i, j, iError, nMltTemp=-1, nLatTemp=-1
    
    ! Testing variables:
    character(len=*), parameter :: NameSub='load_conductances'
    logical :: DoTest, DoTestMe

    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if (DoTest)then
       write(*,*)'IE DEBUG: reading conductance files at'
       write(*,*) NameHalFile
       write(*,*) NamePedFile
    end if

    ! Start with Hall Conductance:
    if(DoTest) write(*,*)NameSub//': Opening Hall cond. file '//NameHalFile
    call open_file(file='IE/'//NameHalFile, status="old")

    ! Skip until DIMENSIONS are found:
    do
       ! Read line; break at EOF.
       read(UnitTmp_, *, iostat=iError) Line
       if (iError /= 0) EXIT
       ! Parse dimensions of arrays:
       if(index(Line,'#DIMENSIONS')>0) then
          read(UnitTmp_, *, iostat=iError) i_cond_nmlts
          read(UnitTmp_, *, iostat=iError) i_cond_nlats
          exit
       end if
    end do

    ! Check if dimensions found.  If not, stop program.
    if( (i_cond_nmlts==-1) .or. (i_cond_nlats==-1) ) call CON_stop(&
         NameSub//' Cannot find #DIMENSION in Hall conductance file.')

    if(DoTest)write(*,*) NameSub//': Size of conductance files (mlt, lat): ', &
         i_cond_nmlts, i_cond_nlats
    
    ! Allocate conductance arrays.  Include MLT ghost cell.
    ! Hall coefficients:
    allocate( hal_a0_up(i_cond_nmlts+1, i_cond_nlats) )
    allocate( hal_a1_up(i_cond_nmlts+1, i_cond_nlats) )
    allocate( hal_a2_up(i_cond_nmlts+1, i_cond_nlats) )
    allocate( hal_a0_do(i_cond_nmlts+1, i_cond_nlats) )
    allocate( hal_a1_do(i_cond_nmlts+1, i_cond_nlats) )
    allocate( hal_a2_do(i_cond_nmlts+1, i_cond_nlats) )
    ! Pedersen coefficients:
    allocate( ped_a0_up(i_cond_nmlts+1, i_cond_nlats) )
    allocate( ped_a1_up(i_cond_nmlts+1, i_cond_nlats) )
    allocate( ped_a2_up(i_cond_nmlts+1, i_cond_nlats) )
    allocate( ped_a0_do(i_cond_nmlts+1, i_cond_nlats) )
    allocate( ped_a1_do(i_cond_nmlts+1, i_cond_nlats) )
    allocate( ped_a2_do(i_cond_nmlts+1, i_cond_nlats) )
    ! Coefficient grids:
    allocate(cond_mlts(i_cond_nmlts+1), cond_lats(i_cond_nlats) )

    ! Read lines until #START:
    do
       read(UnitTmp_, *, iostat=iError) Line
       if (iError /= 0)           EXIT
       if(index(Line,'#START')>0) EXIT
    end do
    
    ! Read and load conductance coefficients & grid:
    do i=1, i_cond_nlats
       do j=1, i_cond_nmlts
          ! Parse a single line:
          read(UnitTmp_,*,iostat=iError) cond_lats(i), cond_mlts(j), &
               hal_a0_up(j,i), hal_a0_do(j,i), &
               hal_a1_up(j,i), hal_a1_do(j,i), &
               hal_a2_up(j,i), hal_a2_do(j,i)
          ! Stop code if IO error:
          if(iError/=0) then
             write(*,*)NameSub//': FILE ERROR at i,j = ', i, j
             call CON_stop(NameSub//': FILE ERROR for Hall input')
          end if
       end do
    end do
    
    ! Close Hall conductance file:
    call close_file
    
    if(DoTest)then
       ! Write out first and last conductances to screen for visual checking:
       write(*,*)NameSub//': Visual check for conductance coeff reading'
       j = 1
       i = 1
       write(*,*) 'At lat, lon = ', cond_lats(i), cond_mlts(j)
       write(*,'(a, 3(1x,E12.4))') '     HALL_UP A0, A1, A2 = ', &
            hal_a0_up(j,i), hal_a1_up(j,i), hal_a2_up(j,i)
       write(*,'(a, 3(1x,E12.4))') '     HALL_DO A0, A1, A2 = ', &
            hal_a0_do(j,i), hal_a1_do(j,i), hal_a2_do(j,i)

       j = i_cond_nmlts
       i = i_cond_nlats
       write(*,*) 'At lat, lon = ', cond_lats(i), cond_mlts(j)
       write(*,'(a, 3(1x,E12.4))') '     HALL_UP A0, A1, A2 = ', &
            hal_a0_up(j,i), hal_a1_up(j,i), hal_a2_up(j,i)
       write(*,'(a, 3(1x,E12.4))') '     HALL_DO A0, A1, A2 = ', &
            hal_a0_do(j,i), hal_a1_do(j,i), hal_a2_do(j,i)
    end if


    
    ! Load Pedersen Conductance:
    if(DoTest) write(*,*)NameSub//': Opening Pedersen cond. file '//NamePedFile
    call open_file(file='IE/'//NamePedFile, status="old")
    
    ! Skip until DIMENSIONS are found
    do
       ! Read line; break at EOF.
       read(UnitTmp_, *, iostat=iError) Line
       if (iError /= 0) EXIT
       ! Parse dimensions of arrays:
       if(index(Line,'#DIMENSIONS')>0) then
          read(UnitTmp_, *, iostat=iError) nMltTemp
          read(UnitTmp_, *, iostat=iError) nLatTemp
          exit
       end if
    end do
    ! Check if dimensions found.  If not, stop program.
    if( (nLatTemp==-1) .or. (nMltTemp==-1) ) call CON_stop(&
         NameSub//' Cannot find #DIMENSION in Pedersen conductance file.')
    ! Check if match Hall file.  
    if( (nLatTemp/=i_cond_nlats) .or. (nMltTemp/=i_cond_nmlts) ) call CON_stop(&
         NameSub//' Hall & Pedersen input file dimensions do not match.')

    ! Read lines until #START:
    do
       read(UnitTmp_, *, iostat=iError) Line
       if(iError /= 0)            EXIT
       if(index(Line,'#START')>0) EXIT
    end do
    
    ! Read and load conductance coefficients & grid:
    do i=1, i_cond_nlats
       do j=1, i_cond_nmlts
          ! Parse a single line:
          read(UnitTmp_,*,iostat=iError) cond_lats(i), cond_mlts(j), &
               ped_a0_up(j,i), ped_a0_do(j,i), &
               ped_a1_up(j,i), ped_a1_do(j,i), &
               ped_a2_up(j,i), ped_a2_do(j,i)
          ! Stop code if IO error:
          if(iError/=0) then
             write(*,*)NameSub//'FILE ERROR at i,j = ', i, j
             call CON_stop(NameSub//' FILE ERROR for Hall input')
          end if
       end do
    end do
    
    ! Close Pedersen conductance file:
    call close_file
    
    ! Wrap values around MLT 00 == 24:
    cond_mlts(i_cond_nmlts+1)   = cond_mlts(1)+24.0
    hal_a0_up(i_cond_nmlts+1,:) = hal_a0_up(1,:)
    ped_a0_up(i_cond_nmlts+1,:) = ped_a0_up(1,:)
    hal_a0_do(i_cond_nmlts+1,:) = hal_a0_do(1,:)
    ped_a0_do(i_cond_nmlts+1,:) = ped_a0_do(1,:)
    hal_a1_up(i_cond_nmlts+1,:) = hal_a1_up(1,:)
    ped_a1_up(i_cond_nmlts+1,:) = ped_a1_up(1,:)
    hal_a1_do(i_cond_nmlts+1,:) = hal_a1_do(1,:)
    ped_a1_do(i_cond_nmlts+1,:) = ped_a1_do(1,:)
    hal_a2_up(i_cond_nmlts+1,:) = hal_a2_up(1,:)
    ped_a2_up(i_cond_nmlts+1,:) = ped_a2_up(1,:)
    hal_a2_do(i_cond_nmlts+1,:) = hal_a2_do(1,:)
    ped_a2_do(i_cond_nmlts+1,:) = ped_a2_do(1,:)
    
    if(DoTest)then
       ! Write out first and last conductances to screen for visual checking:
       write(*,*)NameSub//': Visual check for conductance coeff reading'
       j = 1
       i = 1
       write(*,*) 'At lat, lon = ', cond_lats(i), cond_mlts(j)
       write(*,'(a, 3(1x,E12.4))') '     PEDER_UP A0, A1, A2 = ', &
            ped_a0_up(j,i), ped_a1_up(j,i), ped_a2_up(j,i)
       write(*,'(a, 3(1x,E12.4))') '     PEDER_DO A0, A1, A2 = ', &
            ped_a0_do(j,i), ped_a1_do(j,i), ped_a2_do(j,i)

       j = i_cond_nmlts
       i = i_cond_nlats
       write(*,*) 'At lat, lon = ', cond_lats(i), cond_mlts(j)
       write(*,'(a, 3(1x,E12.4))') '     PEDER_UP A0, A1, A2 = ', &
            ped_a0_up(j,i), ped_a1_up(j,i), ped_a2_up(j,i)
       write(*,'(a, 3(1x,E12.4))') '     PEDER_DO A0, A1, A2 = ', &
            ped_a0_do(j,i), ped_a1_do(j,i), ped_a2_do(j,i)
    end if
    
  end subroutine load_conductances
  !===========================================================================  
  
  subroutine init_mod_ionosphere

    if(allocated(IONO_Phi)) RETURN

    ! Initialize these global grid arrays to 0 (for output before solve)
    allocate(IONO_Phi(2*IONO_nTheta-1,IONO_nPsi));        IONO_Phi = 0
    allocate(IONO_IonNumFlux(2*IONO_nTheta-1,IONO_nPsi)); IONO_IonNumFlux = 0
    allocate(IONO_Joule(2*IONO_nTheta-1,IONO_nPsi));      IONO_Joule = 0
    allocate(IONO_Jr(2*IONO_nTheta-1,IONO_nPsi));         IONO_Jr = 0
    allocate(IONO_Ave_E(2*IONO_nTheta-1,IONO_nPsi));      IONO_Ave_E = 0
    allocate(IONO_Eflux(2*IONO_nTheta-1,IONO_nPsi));      IONO_Eflux = 0
    allocate(IONO_SigmaP(2*IONO_nTheta-1,IONO_nPsi));     IONO_SigmaP = 0
    allocate(IONO_SigmaH(2*IONO_nTheta-1,IONO_nPsi));     IONO_SigmaH = 0

    allocate(IONO_NORTH_PHI(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_PHI(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_X(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Y(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Z(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Theta(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Psi(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_X(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Y(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Z(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Theta(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Psi(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Ex(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Ey(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Ez(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_ETh(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_EPs(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Ex(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Ey(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Ez(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_ETh(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_EPs(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Ux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Uy(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Uz(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_UTh(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_UPs(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Ux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Uy(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Uz(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_UTh(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_UPs(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_EFlux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Ave_E(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_EFlux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Ave_E(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Sigma0(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_SigmaH(IONO_nTheta,IONO_nPsi)); IONO_NORTH_SigmaH = 0
    allocate(IONO_NORTH_SigmaP(IONO_nTheta,IONO_nPsi)); IONO_NORTH_SigmaP = 0
    allocate(IONO_NORTH_SigmaThTh(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_SigmaThPs(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_SigmaPsPs(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Sigma0(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_SigmaH(IONO_nTheta,IONO_nPsi)); IONO_SOUTH_SigmaH = 0
    allocate(IONO_SOUTH_SigmaP(IONO_nTheta,IONO_nPsi)); IONO_SOUTH_SigmaP = 0
    allocate(IONO_SOUTH_SigmaThTh(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_SigmaThPs(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_SigmaPsPs(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_dSigmaThTh_dTheta(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_dSigmaThPs_dTheta(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_dSigmaPsPs_dTheta(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_dSigmaThTh_dPsi(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_dSigmaThPs_dPsi(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_dSigmaPsPs_dPsi(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_dSigmaThTh_dTheta(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_dSigmaThPs_dTheta(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_dSigmaPsPs_dTheta(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_dSigmaThTh_dPsi(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_dSigmaThPs_dPsi(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_dSigmaPsPs_dPsi(IONO_nTheta,IONO_nPsi))
    allocate(SAVE_NORTH_SigmaH(IONO_nTheta,IONO_nPsi))
    allocate(SAVE_NORTH_SigmaP(IONO_nTheta,IONO_nPsi))
    allocate(SAVE_SOUTH_SigmaH(IONO_nTheta,IONO_nPsi))
    allocate(SAVE_SOUTH_SigmaP(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Joule(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Joule(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_IonNumFlux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_IonNumFlux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_JR(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_JTh(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_JPs(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Jx(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Jy(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Jz(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_JR(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_JTh(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_JPs(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Jx(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Jy(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Jz(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_TGCM_JR(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_TGCM_JR(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_Fake_JR(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_Fake_JR(IONO_nTheta,IONO_nPsi))
    allocate(IONO_north_im_jr(IONO_nTheta,IONO_nPsi))
    allocate(IONO_south_im_jr(IONO_nTheta,IONO_nPsi))
    allocate(IONO_north_im_avee(IONO_nTheta,IONO_nPsi))
    allocate(IONO_south_im_avee(IONO_nTheta,IONO_nPsi))
    allocate(IONO_north_im_eflux(IONO_nTheta,IONO_nPsi))
    allocate(IONO_south_im_eflux(IONO_nTheta,IONO_nPsi))
    allocate(IsFilledWithIm(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_invB(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_invB(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_rho(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_rho(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_p(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_p(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_t(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_t(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_dLat(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_dLat(IONO_nTheta,IONO_nPsi))
    allocate(IONO_NORTH_dLon(IONO_nTheta,IONO_nPsi))
    allocate(IONO_SOUTH_dLon(IONO_nTheta,IONO_nPsi))
    allocate(C_A(IONO_nTheta,IONO_nPsi))
    allocate(C_B(IONO_nTheta,IONO_nPsi))
    allocate(C_C(IONO_nTheta,IONO_nPsi))
    allocate(C_D(IONO_nTheta,IONO_nPsi))
    allocate(C_E(IONO_nTheta,IONO_nPsi))

    ! Read empirical conductance values from files:
    call load_conductances()
    
  end subroutine init_mod_ionosphere
  !===========================================================================
  subroutine clean_mod_ionosphere

    if(.not.allocated(IONO_Phi)) RETURN

    deallocate(IONO_Phi)
    deallocate(IONO_IonNumFlux)
    deallocate(IONO_Joule)
    deallocate(IONO_Jr)
    deallocate(IONO_Ave_E)
    deallocate(IONO_Eflux)
    deallocate(IONO_SigmaP)
    deallocate(IONO_SigmaH)
    deallocate(IONO_NORTH_PHI)
    deallocate(IONO_SOUTH_PHI)
    deallocate(IONO_NORTH_X)
    deallocate(IONO_NORTH_Y)
    deallocate(IONO_NORTH_Z)
    deallocate(IONO_NORTH_Theta)
    deallocate(IONO_NORTH_Psi)
    deallocate(IONO_SOUTH_X)
    deallocate(IONO_SOUTH_Y)
    deallocate(IONO_SOUTH_Z)
    deallocate(IONO_SOUTH_Theta)
    deallocate(IONO_SOUTH_Psi)
    deallocate(IONO_NORTH_Ex)
    deallocate(IONO_NORTH_Ey)
    deallocate(IONO_NORTH_Ez)
    deallocate(IONO_NORTH_ETh)
    deallocate(IONO_NORTH_EPs)
    deallocate(IONO_SOUTH_Ex)
    deallocate(IONO_SOUTH_Ey)
    deallocate(IONO_SOUTH_Ez)
    deallocate(IONO_SOUTH_ETh)
    deallocate(IONO_SOUTH_EPs)
    deallocate(IONO_NORTH_Ux)
    deallocate(IONO_NORTH_Uy)
    deallocate(IONO_NORTH_Uz)
    deallocate(IONO_NORTH_UTh)
    deallocate(IONO_NORTH_UPs)
    deallocate(IONO_SOUTH_Ux)
    deallocate(IONO_SOUTH_Uy)
    deallocate(IONO_SOUTH_Uz)
    deallocate(IONO_SOUTH_UTh)
    deallocate(IONO_SOUTH_UPs)
    deallocate(IONO_NORTH_EFlux)
    deallocate(IONO_NORTH_Ave_E)
    deallocate(IONO_SOUTH_EFlux)
    deallocate(IONO_SOUTH_Ave_E)
    deallocate(IONO_NORTH_Sigma0)
    deallocate(IONO_NORTH_SigmaH)
    deallocate(IONO_NORTH_SigmaP)
    deallocate(IONO_NORTH_SigmaThTh)
    deallocate(IONO_NORTH_SigmaThPs)
    deallocate(IONO_NORTH_SigmaPsPs)
    deallocate(IONO_SOUTH_Sigma0)
    deallocate(IONO_SOUTH_SigmaH)
    deallocate(IONO_SOUTH_SigmaP)
    deallocate(IONO_SOUTH_SigmaThTh)
    deallocate(IONO_SOUTH_SigmaThPs)
    deallocate(IONO_SOUTH_SigmaPsPs)
    deallocate(IONO_NORTH_dSigmaThTh_dTheta)
    deallocate(IONO_NORTH_dSigmaThPs_dTheta)
    deallocate(IONO_NORTH_dSigmaPsPs_dTheta)
    deallocate(IONO_NORTH_dSigmaThTh_dPsi)
    deallocate(IONO_NORTH_dSigmaThPs_dPsi)
    deallocate(IONO_NORTH_dSigmaPsPs_dPsi)
    deallocate(IONO_SOUTH_dSigmaThTh_dTheta)
    deallocate(IONO_SOUTH_dSigmaThPs_dTheta)
    deallocate(IONO_SOUTH_dSigmaPsPs_dTheta)
    deallocate(IONO_SOUTH_dSigmaThTh_dPsi)
    deallocate(IONO_SOUTH_dSigmaThPs_dPsi)
    deallocate(IONO_SOUTH_dSigmaPsPs_dPsi)
    deallocate(SAVE_NORTH_SigmaH)
    deallocate(SAVE_NORTH_SigmaP)
    deallocate(SAVE_SOUTH_SigmaH)
    deallocate(SAVE_SOUTH_SigmaP)
    deallocate(IONO_NORTH_Joule)
    deallocate(IONO_SOUTH_Joule)
    deallocate(IONO_NORTH_IonNumFlux)
    deallocate(IONO_SOUTH_IonNumFlux)
    deallocate(IONO_NORTH_JR)
    deallocate(IONO_NORTH_JTh)
    deallocate(IONO_NORTH_JPs)
    deallocate(IONO_NORTH_Jx)
    deallocate(IONO_NORTH_Jy)
    deallocate(IONO_NORTH_Jz)
    deallocate(IONO_SOUTH_JR)
    deallocate(IONO_SOUTH_JTh)
    deallocate(IONO_SOUTH_JPs)
    deallocate(IONO_SOUTH_Jx)
    deallocate(IONO_SOUTH_Jy)
    deallocate(IONO_SOUTH_Jz)
    deallocate(IONO_NORTH_TGCM_JR)
    deallocate(IONO_SOUTH_TGCM_JR)
    deallocate(IONO_NORTH_Fake_JR)
    deallocate(IONO_SOUTH_Fake_JR)
    deallocate(IONO_north_im_jr)
    deallocate(IONO_south_im_jr)
    deallocate(IONO_north_im_avee)
    deallocate(IONO_south_im_avee)
    deallocate(IONO_north_im_eflux)
    deallocate(IONO_south_im_eflux)
    deallocate(IsFilledWithIm)
    deallocate(IONO_NORTH_invB)
    deallocate(IONO_SOUTH_invB)
    deallocate(IONO_NORTH_rho)
    deallocate(IONO_SOUTH_rho)
    deallocate(IONO_NORTH_p)
    deallocate(IONO_SOUTH_p)
    deallocate(IONO_NORTH_t)
    deallocate(IONO_SOUTH_t)
    deallocate(IONO_NORTH_dLat)
    deallocate(IONO_SOUTH_dLat)
    deallocate(IONO_NORTH_dLon)
    deallocate(IONO_SOUTH_dLon)
    deallocate(C_A)
    deallocate(C_B)
    deallocate(C_C)
    deallocate(C_D)
    deallocate(C_E)

    ! Clean up all conductance arrays:
    ! Hall coefficients:
    deallocate( hal_a0_up )
    deallocate( hal_a1_up )
    deallocate( hal_a2_up )
    deallocate( hal_a0_do )
    deallocate( hal_a1_do )
    deallocate( hal_a2_do )
    ! Pedersen coefficients:
    deallocate( ped_a0_up )
    deallocate( ped_a1_up )
    deallocate( ped_a2_up )
    deallocate( ped_a0_do )
    deallocate( ped_a1_do )
    deallocate( ped_a2_do )
    ! Coefficient grids:
    deallocate(cond_mlts, cond_lats)
    
  end subroutine clean_mod_ionosphere

end module ModIonosphere
