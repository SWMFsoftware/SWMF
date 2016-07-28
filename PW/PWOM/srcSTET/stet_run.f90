!============================================================================
! run STET to get a new steady state or to advance some amount of time for 
! a particular line, iLine.
subroutine stet_run(iLine,IsOpen,DoCouplePWOM)
  use ModSeGrid, only:create_se_test_grid,set_grid_pot,nLine,nPoint,nIono,&
       nPlas,Efield_IC,DoIncludePotential,IsVerbose
  use ModSeBackground,only: plot_background,plot_ephoto_prod,&
       eThermalDensity_IC,eThermalTemp_IC,nNeutralSpecies,&
       NeutralDens1_IIC,ePhotoProdSpec1_IIC,NeutralDens2_IIC,ePhotoProdSpec2_IIC
  use ModSeState,only: check_time,check_time_pot,&
       update_se_state_iono,update_se_state_iono_pot,update_se_state, &
       update_se_state_pot,initplas, initplas_pot, initiono,initiono_pot,&
       calc_integrated_output,HeatingRate_IC,NumberDens_IC,NumberFlux_IC,&
       TotalIonizationRate_IC,&
       iphiup,iphidn,phiup,phidn,liphiup,liphidn,lphiup,lphidn,&
       specup, specdn, epsilon,delt, Time
  use ModSeCross,only: SIGS,SIGI,SIGA    
  use ModSePlot  
  implicit none
  
  integer, intent(in) :: iLine
  logical, intent(in) :: IsOpen,DoCouplePWOM

  integer :: flag=1, nStep=0
  logical :: DoSavePreviousAndReset = .true.

  logical,parameter :: IsIono1=.true.
  real, parameter :: DtCouplePWOM=120.0 !this should come from PWOM in future
  !--------------------------------------------------------------------------
  
  
  !\
  ! The main update
  !/

  ! Set the timestep and convergence criteria
  delt=1.0e5
  epsilon = 0.05

  ! Define the initial state
  iphiup(iLine,:,:,:)=0.0
  iphidn(iLine,:,:,:)=0.0
  
  liphiup(iLine,:,:,:)=0.0
  liphidn(iLine,:,:,:)=0.0
  
  phiup(iLine,:,:,:)=0.0
  phidn(iLine,:,:,:)=0.0
  
  lphiup(iLine,:,:,:)=0.0
  lphidn(iLine,:,:,:)=0.0
  
  nStep = 0
  
  !Start Timeloop
  flag = 1
  TIME_LOOP: do while (flag == 1)
     if(.not.DoIncludePotential) then
        ! Initialize the ionosphere
        if(IsVerbose) write(*,*) 'Initializing iono'
        call initiono(iLine,DoSavePreviousAndReset)
        
        ! update the SE state for iono1
        call update_se_state_iono(iLine,IsIono1,eThermalDensity_IC(iLine,:),&
             eThermalTemp_IC(iLine,:),nNeutralSpecies,&
             NeutralDens1_IIC(iLine,:,:),SIGS,SIGI,SIGA,&
             ePhotoProdSpec1_IIC(iLine,:,:))
        
        ! update the SE state for iono2 if field line is closed
        if (.not.IsOpen) then
           call update_se_state_iono(iLine,.not.IsIono1,&
                eThermalDensity_IC(iLine,:),&
                eThermalTemp_IC(iLine,:),nNeutralSpecies,&
                NeutralDens2_IIC(iLine,:,:),SIGS,SIGI,SIGA,&
                ePhotoProdSpec2_IIC(iLine,:,:))
        endif
        
        ! Initialize the plasmasphere
        if(IsVerbose) write(*,*) 'Initializing plasmasphere'
        call initplas(iLine,DoSavePreviousAndReset)
        
        ! update the SE state
        if(IsVerbose) write(*,*) 'update se state'
        call update_se_state(iLine, eThermalDensity_IC(iLine,:),&
             eThermalTemp_IC(iLine,:),IsOpen)
        
        ! check convergence
        if(IsVerbose) write(*,*) 'check for convergence'
        call check_time(iLine,flag)
     else
                ! Initialize the ionosphere
        if(IsVerbose) write(*,*) 'Initializing iono1'
        call initiono_pot(iLine,DoSavePreviousAndReset)
        
        ! update the SE state for iono1

        call update_se_state_iono_pot(iLine,IsIono1,eThermalDensity_IC(iLine,:),&
             eThermalTemp_IC(iLine,:),nNeutralSpecies,&
             NeutralDens1_IIC(iLine,:,:),SIGS,SIGI,SIGA,&
             ePhotoProdSpec1_IIC(iLine,:,:))
        
        ! update the SE state for iono2 if field line is closed
        if (.not.IsOpen) then
           call update_se_state_iono_pot(iLine,.not.IsIono1,&
                eThermalDensity_IC(iLine,:),&
                eThermalTemp_IC(iLine,:),nNeutralSpecies,&
                NeutralDens2_IIC(iLine,:,:),SIGS,SIGI,SIGA,&
                ePhotoProdSpec2_IIC(iLine,:,:))
        endif

        ! Initialize the plasmasphere
        if(IsVerbose) write(*,*) 'Initializing plasmasphere'
        call initplas_pot(iLine,DoSavePreviousAndReset)
        
        ! update the SE state
        if(IsVerbose) write(*,*) 'update se state'
        call update_se_state_pot(iLine, eThermalDensity_IC(iLine,:),&
             eThermalTemp_IC(iLine,:),IsOpen)
        
        ! check convergence
        if(IsVerbose) write(*,*) 'check for convergence'
        call check_time_pot(iLine,flag)
     endif
     ! increment step
     nStep=nStep+1
     
  end do TIME_LOOP

  ! Find heating rate, Se number density and flux
  if(IsVerbose) write(*,*) 'Getting Integrals for output'
  call calc_integrated_output(iLine,nNeutralSpecies,eThermalDensity_IC(iLine,:))
  
  !\
  ! plot model output
  !/
  if ((Time-real(floor((Time+1.0e-5)/DtSavePlot))*Dtsaveplot<DtCouplePWOM) &
       .or. .not.DoCouplePWOM)then 
     if(DoIncludePotential) then
        
        ! plot background
        call plot_background(iLine,nStep,time)
        
        ! plot ephoto production
        call plot_ephoto_prod(iLine,nStep,time)
        
        !     call plot_state_pot(iLine,90,nStep,time,iphiup,iphidn,phiup,phidn)     
        !     call plot_state_pot(iLine,1,nStep,time,iphiup,iphidn,phiup,phidn)
        !     call plot_state_pot(iLine,2,nStep,time,iphiup,iphidn,phiup,phidn)
        !     call plot_state_pot(iLine,3,nStep,time,iphiup,iphidn,phiup,phidn)
        !     call plot_state_pot(iLine,4,nStep,time,iphiup,iphidn,phiup,phidn)
        !     call plot_state_pot(iLine,5,nStep,time,iphiup,iphidn,phiup,phidn)
        !     call plot_state_pot(iLine,6,nStep,time,iphiup,iphidn,phiup,phidn)
        
        call plot_omni_iono_pot(iLine,nStep,time,specup,specdn,.true.)
        call plot_omni_iono_pot(iLine,nStep,time,specup,specdn,.false.)
        
        call plot_omni_pot(iLine,nStep,time,specup,specdn)
        call plot_along_field(iLine,time,HeatingRate_IC,NumberDens_IC,NumberFlux_IC,TotalIonizationRate_IC)
        
        
     else
        ! plot background
        call plot_background(iLine,nStep,time)
        
        ! plot ephoto production
        call plot_ephoto_prod(iLine,nStep,time)
        
        call plot_state(iLine,nStep,time,iphiup,iphidn,phiup,phidn)
        call plot_omni_iono(iLine,nStep,time,specup,specdn,.true.)
        call plot_omni_iono(iLine,nStep,time,specup,specdn,.false.)
        call plot_omni_line(iLine,nStep,time,specup,specdn)
        call plot_along_field(iLine,time,HeatingRate_IC,NumberDens_IC,NumberFlux_IC,TotalIonizationRate_IC)
     end if
  end if
end subroutine stet_run
