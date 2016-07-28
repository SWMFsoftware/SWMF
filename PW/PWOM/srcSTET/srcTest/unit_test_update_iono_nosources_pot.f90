program unit_test_update_iono_nosources
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

  call se_update_state_iono_nosources_pot

end program unit_test_update_iono_nosources

!============================================================================
! UNIT test for update_se_state_iono_pot. Tests update in case with no 
! plasmasphere included and no sources, just a plasmasphere BC and a potential. 
! Ionospheres are kept decoupled. only testing transport in iono
subroutine se_update_state_iono_nosources_pot
  use ModSeGrid, only:create_se_test_grid,set_grid_pot,nLine,nPoint,nIono,&
       nPlas,nTop,nEnergy,Efield_IC,FieldLineGrid_IC,nThetaAlt_IIC,nThetaAlt_II
  use ModSePlot, only: plot_state,plot_state_pot,plot_omni_iono_pot
  use ModSeState,only: allocate_state_arrays, iphiup,iphidn,phiup,phidn,&
       InitIono_pot,check_time_pot,delt,epsilon, update_se_state_iono_pot, &
       liphiup,liphidn,lphiup,lphidn,specup, specdn    
  integer :: iLine=1, flag=1
  logical :: DoSavePreviousAndReset = .true.
  integer :: nStep
  logical :: IsIono1
  real    :: Ap(7), F107=80, F107A=80, t=0
  
  ! The values below are normally taken from a module but are delclared here 
  ! since we are doing an simple idealized test and do not want to use those
  ! modules in this unit test.
  
  integer, parameter :: nNeutralSpecies = 3
  
  ! sig arrays
  real, allocatable :: SIGS(:,:),SIGI(:,:,:),SIGA(:,:,:)
  
  ! Neutrals
  real, allocatable :: NeutralDens1_IIC(:,:,:),NeutralDens2_IIC(:,:,:)
  
  ! Arrays that hold the photo electron production spectrum in iono 1 or 2
  real, allocatable :: ePhotoProdSpec1_IIC(:,:,:),ePhotoProdSpec2_IIC(:,:,:)
  
  ! thermal background  
  real, allocatable :: eThermalDensity_C(:),eThermalTemp_C(:)
  
  integer :: iIono
  real,parameter :: cCmToM = 1.0e-2  
!--------------------------------------------------------------------------
  
  ! First set up the grid that we will update the state in (this is the same 
  ! as the unit test for the grid).
  write(*,*) 'creating grid'
  call create_se_test_grid
  
  ! input electric field (assume constant in plasmapshere and zero in iono)
  ! choose values to have a 5 V drop from top of iono to equator
  Efield_IC(iLine,:)=0.0
  Efield_IC(iLine,nIono+1:nTop-1)=5.0&
       /((FieldLineGrid_IC(iLine,nTop)-FieldLineGrid_IC(iLine,nIono))*cCmToM)
  Efield_IC(iLine,nTop+1:nPoint-nIono)=-5.0&
       /((FieldLineGrid_IC(iLine,nTop)-FieldLineGrid_IC(iLine,nIono))*cCmToM)
  
  ! set the region of exisitence for the defined Efield
  call set_grid_pot(iLine)
  

  !allocate background and production arrays
  if(.not.allocated(NeutralDens1_IIC)) &
       allocate(NeutralDens1_IIC(nLine,nNeutralSpecies,nIono))
  if(.not.allocated(NeutralDens2_IIC)) &
       allocate(NeutralDens2_IIC(nLine,nNeutralSpecies,nIono))
  
  if(.not.allocated(ePhotoProdSpec1_IIC)) &
       allocate(ePhotoProdSpec1_IIC(nLine,nEnergy,nIono))
  
  if(.not.allocated(ePhotoProdSpec2_IIC)) &
       allocate(ePhotoProdSpec2_IIC(nLine,nEnergy,nIono))
  
  if(.not.allocated(eThermalDensity_C)) &
       allocate(eThermalDensity_C(nPoint))
  if(.not.allocated(eThermalTemp_C)) &
       allocate(eThermalTemp_C(nPoint))
  
  !allocate sig arrays if not already done
  if (.not.allocated(SIGS)) allocate(SIGS(nNeutralSpecies,nEnergy)) 
  if (.not.allocated(SIGI)) allocate(SIGI(nNeutralSpecies,nEnergy,nEnergy)) 
  if (.not.allocated(SIGA)) allocate(SIGA(nNeutralSpecies,nEnergy,nEnergy)) 
  
  
  ! Allocate the state arrays
  write(*,*) 'allocating state arrays'
  call allocate_state_arrays
  
  !Set collision cross section to zero
  SIGS(:,:)   = 0.0
  SIGI(:,:,:) = 0.0
  SIGA(:,:,:) = 0.0 
  
  ! Set the photoelectron production to zero
  ePhotoProdSpec1_IIC = 0.0
  ePhotoProdSpec2_IIC(:,:,:) = 0.0
  
  ! Set the neutral density to zero
  NeutralDens1_IIC(:,:,:)=0.0
  NeutralDens2_IIC(:,:,:)=0.0
  
  ! Fill the background arrays
  write(*,*) 'filling background arrays'
  
  ! Define the initial state in the ionosphere
  iphiup(iLine,:,:,:)=0.0
  iphidn(iLine,:,:,:)=0.0
  
  liphiup(iLine,:,:,:)=0.0
  liphidn(iLine,:,:,:)=0.0
  
  phiup(iLine,:,:,:)=1.0e5
  phidn(iLine,:,:,:)=1.0e5
  
  lphiup(iLine,:,:,:)=1.0e5
  lphidn(iLine,:,:,:)=1.0e5
  
  ! Fill the background arrays
  eThermalDensity_C(:) = 0.00001
  eThermalTemp_C(:) = 1.0
  
  
  ! Set the timestep and convergence criteria
  delt=1.0e5
  epsilon = 0.1
  
  write(*,*) 'Starting Time loop'
  nStep = 0
  
  !do test for iono1
  IsIono1 = .true.
  TIME_LOOP: do while (flag == 1)
     ! Initialize the plasmasphere
     write(*,*) 'Initializing ionosphere'
     call initiono_pot(1,DoSavePreviousAndReset)
     
     ! update the SE state for iono1
     write(*,*) 'update se state for iono1'
     write(*,*) 'max(ePhotoProdSpec1_IIC)',maxval(ePhotoProdSpec1_IIC)
     
     ! update state for iono1
     call update_se_state_iono_pot(iLine,IsIono1,eThermalDensity_C,&
          eThermalTemp_C,nNeutralSpecies,&
          NeutralDens1_IIC(iLine,:,:),SIGS,SIGI,SIGA,&
          ePhotoProdSpec1_IIC(iLine,:,:))
     
     ! update the SE state for iono2
     !       write(*,*) 'update se state for iono2'
     call update_se_state_iono_pot(iLine,.not.IsIono1,eThermalDensity_C,&
          eThermalTemp_C,nNeutralSpecies,&
          NeutralDens2_IIC(iLine,:,:),SIGS,SIGI,SIGA,&
          ePhotoProdSpec2_IIC(iLine,:,:))
     
     ! check convergence
     write(*,*) 'check for convergence'
     call check_time_pot(1,flag)
     
     ! increment step
     nStep=nStep+1
     
     !       !debug stuff
     !       do iIono=1,2*nIono
     !          write(*,*) 'iIono,iphiup(iLine,0:5,iIono,5)',iIono,iphiup(iLine,0:5,iIono,1)
     !       enddo
     
     !       do iIono=1,nPlas
     !          write(*,*) 'iPlas,phiup(iLine,0:5,iPlas,5)',iIono,phiup(iLine,0:5,iIono,1)
     !       enddo
     
     ! plot output
     call plot_state_pot(1,90,nStep,time,iphiup,iphidn,phiup,phidn)
     
     call plot_state_pot(1,1,nStep,time,iphiup,iphidn,phiup,phidn)
     call plot_state_pot(1,2,nStep,time,iphiup,iphidn,phiup,phidn)
     call plot_state_pot(1,3,nStep,time,iphiup,iphidn,phiup,phidn)
     call plot_state_pot(1,4,nStep,time,iphiup,iphidn,phiup,phidn)
     call plot_state_pot(1,5,nStep,time,iphiup,iphidn,phiup,phidn)
     call plot_state_pot(1,6,nStep,time,iphiup,iphidn,phiup,phidn)


     !      call plot_omni_iono(iLine,nStep,time,specup,specdn,.true.)
     call plot_omni_iono_pot(iLine,nStep,time,specup,specdn,.false.)
     
  end do TIME_LOOP
  
end subroutine se_update_state_iono_nosources_pot


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

