program unit_test_update_iono_plas_pot
  use ModSeGrid
  use ModSeMpi
  use ModMPI
  use CON_planet, ONLY: init_planet_const, set_planet_defaults

  integer :: iError
  
  !-----------------------------------------------------------------------------

  !****************************************************************************
  ! Initiallize MPI and get number of processors and rank of given processor
  !****************************************************************************

  write(*,*) 'Initiallizing MPI'

  !---------------------------------------------------------------------------
  call MPI_INIT(iError)
  iComm = MPI_COMM_WORLD

  call MPI_COMM_RANK(iComm,iProc,iError)
  call MPI_COMM_SIZE(iComm,nProc,iError)

  !\
  ! Initialize the planetary constant library and set Earth
  ! as the default planet.
  !/
  write(*,*) 'Initiallizing Planet'

  call init_planet_const
  call set_planet_defaults


  call se_update_state_iono_plas_pot_test
end program unit_test_update_iono_plas_pot

  !============================================================================
  ! UNIT test for update_se_state_iono. Tests update in case with no 
  ! plasmasphere included. Ionospheres are kept decoupled
  subroutine se_update_state_iono_plas_pot_test
    use ModSeBackground
    use ModSeGrid, only:create_se_test_grid,set_grid_pot,nLine,nPoint,nIono,&
       nPlas,nTop,nEnergy,Efield_IC,FieldLineGrid_IC,&
       nThetaAlt_IIC,nThetaAlt_II,Lshell_I
    use ModSeCross,only: SIGS,SIGI,SIGA
    use ModSePlot, only: plot_state,plot_state_pot,plot_omni_iono_pot,&
         plot_omni_pot, plot_along_field
    use ModSeState,only: allocate_state_arrays, calc_integrated_output,&
         iphiup,iphidn,phiup,phidn,&
         initiono_pot,initplas_pot,check_time_pot,delt,epsilon,&
         update_se_state_iono_pot,update_se_state_pot,&
         liphiup,liphidn,lphiup,lphidn,specup, specdn, &
         HeatingRate_IC,NumberDens_IC,NumberFlux_IC
    use ModNumConst,    ONLY: cRadToDeg
    implicit none
    integer :: iLine=1, flag=1
    logical :: DoSavePreviousAndReset = .true.
    integer :: nStep
    logical :: IsIono1, IsOpen
    real    :: Ap(7), F107=80, F107A=80, time=0
    real,parameter :: cCmToM = 1.0e-2  
    !--------------------------------------------------------------------------

    ! First set up the grid that we will update the state in (this is the same 
    ! as the unit test for the grid).
    write(*,*) 'creating grid'
    call create_se_test_grid

    ! input electric field (assume constant in plasmapshere and zero in iono)
    ! choose values to have a 5 V drop from top of iono to equator
    Efield_IC(iLine,:)=0.0
    Efield_IC(iLine,nIono+1:nTop-1)=2.0&
         /((FieldLineGrid_IC(iLine,nTop)-FieldLineGrid_IC(iLine,nIono))*cCmToM)
    Efield_IC(iLine,nTop+1:nPoint-nIono)=-2.0&
         /((FieldLineGrid_IC(iLine,nTop)-FieldLineGrid_IC(iLine,nIono))*cCmToM)
    
    ! set the region of exisitence for the defined Efield
    call set_grid_pot(iLine)

    ! Allocate the state arrays
    write(*,*) 'allocating state arrays'
    call allocate_state_arrays
    
    ! Allocate the background right
    write(*,*) 'allocating background arrays'
    call allocate_background_arrays
    
    !align dipole and rotation
    DoAlignDipoleRot = .true.
    
    ! set location of field line
    mLat_I(iLine)=acos(sqrt(1.0/Lshell_I(iLine)))*cRadToDeg
    mLon_I(iLine)=0.0

    !set glat and glon coords
    call set_footpoint_locations(iLine)

    ! Get the neutral atmosphere and photo e production spectrum
    AP(:)=4.0
    call get_neutrals_and_pe_spectrum(iLine,F107,F107A,AP)

    ! Fill the background arrays
    write(*,*) 'filling background arrays'
    call fill_thermal_plasma_empirical(iLine,F107,F107A,time)
    
    ! Define the initial state in the ionosphere
    iphiup(iLine,:,:,:)=0.0
    iphidn(iLine,:,:,:)=0.0
    
    liphiup(iLine,:,:,:)=0.0
    liphidn(iLine,:,:,:)=0.0
    
    phiup(iLine,:,:,:)=0.0
    phidn(iLine,:,:,:)=0.0
    
    iphiup(iLine,:,:,:)=0.0
    iphidn(iLine,:,:,:)=0.0
 
    
    ! Set the timestep and convergence criteria
    delt=1.0e5
    epsilon = 0.1

    write(*,*) 'Starting Time loop'
    nStep = 0

    !set IsIono1
    IsIono1 = .true.

    !do test for closed line
    IsOpen = .false.
    
    TIME_LOOP: do while (flag == 1)
       ! Initialize the plasmasphere
       write(*,*) 'Initializing ionosphere'
       call initiono_pot(1,DoSavePreviousAndReset)
       
       ! update the SE state for iono1
       write(*,*) 'update se state for iono1'
       call update_se_state_iono_pot(iLine,IsIono1,eThermalDensity_IC(iLine,:),&
            eThermalTemp_IC(iLine,:),nNeutralSpecies,&
            NeutralDens1_IIC(iLine,:,:),SIGS,SIGI,SIGA,&
            ePhotoProdSpec1_IIC(iLine,:,:))

       ! update the SE state for iono2
       write(*,*) 'update se state for iono2'
       call update_se_state_iono_pot(iLine,.not.IsIono1,&
            eThermalDensity_IC(iLine,:),&
            eThermalTemp_IC(iLine,:),nNeutralSpecies,&
            NeutralDens2_IIC(iLine,:,:),SIGS,SIGI,SIGA,&
            ePhotoProdSpec2_IIC(iLine,:,:))

       ! Initialize the plasmasphere
       write(*,*) 'Initializing plasmasphere'
       call initplas_pot(1,DoSavePreviousAndReset)
       
       ! update the SE state
       write(*,*) 'update se state_pot'
       call update_se_state_pot(1, eThermalDensity_IC(iLine,:),&
            eThermalTemp_IC(iLine,:),IsOpen)     
       
       ! check convergence
       write(*,*) 'check for convergence'
       call check_time_pot(1,flag)
       
       ! increment step
       nStep=nStep+1
       
       ! Find heating rate, Se number density and flux
       write(*,*) 'Getting Integrals for output'
       call calc_integrated_output(iLine,eThermalDensity_IC(iLine,:))
       write(*,*)'done'
       ! plot output
       call plot_state_pot(1,90,nStep,time,iphiup,iphidn,phiup,phidn)
       
       call plot_state_pot(1,1,nStep,time,iphiup,iphidn,phiup,phidn)
       call plot_state_pot(1,2,nStep,time,iphiup,iphidn,phiup,phidn)
       call plot_state_pot(1,3,nStep,time,iphiup,iphidn,phiup,phidn)
       call plot_state_pot(1,4,nStep,time,iphiup,iphidn,phiup,phidn)
       call plot_state_pot(1,5,nStep,time,iphiup,iphidn,phiup,phidn)
       call plot_state_pot(1,6,nStep,time,iphiup,iphidn,phiup,phidn)
       
       
       call plot_omni_iono_pot(iLine,nStep,time,specup,specdn,.true.)
       call plot_omni_iono_pot(iLine,nStep,time,specup,specdn,.false.)
       call plot_omni_pot     (iline,nStep,time,specup,specdn)

    end do TIME_LOOP
    
  end subroutine se_update_state_iono_plas_pot_test


!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  use ModSeMpi, ONLY : iProc,iComm
  use ModMpi
  implicit none
  character (len=*), intent(in) :: StringError

  ! Local variables:
  integer :: iError,nError
  !----------------------------------------------------------------------------

  write(*,*)'Stopping execution! me=',iProc,&
       ' with msg:'
  write(*,*)StringError
  call MPI_abort(iComm, nError, iError)
  stop

end subroutine CON_stop

subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe

  DoTest = .false.; DoTestMe = .false.

end subroutine CON_set_do_test

subroutine CON_io_unit_new(iUnit)

  use ModIoUnit, ONLY: io_unit_new
  implicit none
  integer, intent(out) :: iUnit

  iUnit = io_unit_new()

end subroutine CON_io_unit_new



