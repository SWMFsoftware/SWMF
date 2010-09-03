!^CFG COPYRIGHT UM
! ModUser for flux emergence problems on Sun
! We calculate the thin radiative cooling, vertical damping, 
! and call coronal heating from PARAM.in to model the solar 
! atmosphere, following Abbett (2007). The implementation of
! tabular equation of state takes into account of the 
! ionization energy. 
!
! 2010-03-01: committed in
!=========================================================================

module ModUser
  use ModMain,        ONLY: GlobalBlk
  use ModSize,        ONLY: nI,nJ,nK,MinI,MaxI,MinJ,MaxJ,MinK,MaxK,MaxBlock
  use ModUserEmpty ,                                   &
       IMPLEMENTED1 => user_read_inputs,               &
       IMPLEMENTED2 => user_init_session,              &
       IMPLEMENTED3 => user_set_ICs,                   &
       IMPLEMENTED4 => user_initial_perturbation,      &
       IMPLEMENTED5 => user_set_outerbcs,              &
       IMPLEMENTED6 => user_calc_sources,              &
       IMPLEMENTED7 => user_update_states,             &
       IMPLEMENTED8 => user_set_plot_var,              &
       IMPLEMENTED9 => user_material_properties
  
  include 'user_module.h'

  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'Flux Emergence, Fang 2010-03-01'
  !\
  ! UseVerticalDamping: adds damping to vertical velocity
  ! UseThinRadiation:   adds thin radiative cooling
  ! UseUniformInitialState: initialize the domain with uniform 
  !                      density and temperature
  ! InitialDensity, InitialTemperature: the initial values if
  !                      UseUniformInitialState
  ! InitialBx, InitialBy, InitalBz: initial magnetic field, added through
  !                      user_initial_perturbation
  ! RhoThinCutoff: the cutoff density for thin radiation( thin radiation =0 
  !                      if rho> RhoThinCutoff)
  ! NumberDensFloor: Minimum density value to prevent negative coronal density 
  ! TimeVerticalDamping: timescale over which the vertical velocity damps
  ! z_photo: height of photosphere
  ! TemperatureGradient: gradient of temperature at the bottom of domain
  !/
  logical :: UseVerticalDamping = .false.
  logical :: UseThinRadiation = .false.
  logical :: UseUniformInitialState = .false.
  real    :: InitialDensity, InitialTemperature, InitialBx, &
       InitialBy, InitialBz
  real    :: RhoThinCutoff, NumberDensFloor, TimeVerticalDamping
  real    :: z_photo, TemperatureGradient

  !rstari = 0.594354e-3/8.31, mu = 0.594354 set in init_session
  real ::  mu, rstari

  ! Flux Rope Variables
  ! x2c_rope: y coord of rope axis
  ! x3c_rope: z coord of rope axis
  ! ra_rope : radius of gaussian decay of the magnitude of B field
  ! qfac_rope: twisting factor
  ! lamb_rope: length of buoyant section
  ! b0_rope: magnetic field strength at the rope center
  ! buoyancy_rope : amount of buoyancy in the central section
  !/
  logical :: UseRope = .false.
  real    :: x2c_rope, x3c_rope,ra_rope,qfac_rope,lamb_rope,b0_rope, &
       buoyancy_rope
  !\
  ! Indexes of lookup EOS tables, CHIANTI table, initial relaxed reference state
  !/
  integer :: iTablePPerE = -1, iTableCvTe = -1, iTablePERhoT = -1, &
       iTableChianti = -1, iTableInitialState = -1, iTableOpacity
  ! Temperature at cell center
  real, allocatable:: Temperature_GB(:,:,:,:)
  
contains
  
  !=========================================================================
  subroutine user_read_inputs
    use ModMain,      ONLY: lverbose, UseGravity
    use ModProcMH,    ONLY: iProc
    use ModReadParam, ONLY: read_var, read_line, read_command
    use ModIO,        ONLY: write_prefix, write_myname, iUnitOut
    use ModCoronalHeating, ONLY: DtUpdateFlux, UnsignedFluxHeight

    implicit none
    
    integer:: i
    character (len=100) :: NameCommand
    !-----------------------------------------------------------------------
    if(iProc==0.and.lVerbose > 0)then
       call write_prefix; write(iUnitOut,*)'user_read_input RADMHD starts'
    endif
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#RADMHD")
          call read_var('z_photo',z_photo)
          call read_var('UnsignedFluxHeight',UnsignedFluxHeight)
          call read_Var('UseThinRadiation', UseThinRadiation)
          call read_var('RhoThinCutoff',RhoThinCutoff)
          call read_Var('UseVerticalDamping',UseVerticalDamping)
          call read_var('TimeVerticalDamping', TimeVerticalDamping)
          call read_var('TemperatureGradient',TemperatureGradient)
          call read_var('DtUptateFlux',DtUpdateFlux)
          call read_var('UseUniformInitalState', UseUniformInitialState)
          call read_var('InitialDensity', InitialDensity)
          call read_var('InitialTemperature',InitialTemperature)
          call read_var('InitialBx',InitialBx)
          call read_var('InitialBy',InitialBy)
          call read_var('InitialBz',InitialBz)
          call read_var('NumberDensFloor',NumberDensFloor)
       case('#ROPE')
          call read_var('UseRope',UseRope)
          call read_var('x2c_rope',x2c_rope)
          call read_var('x3c_rope',x3c_rope)
          call read_var('ra_rope',ra_rope)
          call read_var('qfac_rope',qfac_rope)
          call read_var('lamb_rope',lamb_rope)
          call read_var('buoyancy_rope',buoyancy_rope)
          call read_var('b0_rope',b0_rope)
       case('#USERINPUTEND')
          if(iProc==0.and.lVerbose > 0)then
             call write_prefix; write(iUnitOut,*)'user_read_input RADMHD ends'
          endif
          EXIT
       case default
          if(iProc==0) then
             call write_myname; write(*,*) &
                  'ERROR: Invalid user defined #COMMAND in user_read_inputs. '
             write(*,*) '--Check user_read_inputs for errors'
             write(*,*) '--Check to make sure a #USERINPUTEND command was used'
             write(*,*) '  *Unrecognized command was: '//NameCommand
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if
       end select
    end do
  end subroutine user_read_inputs
  !==========================================================================
  subroutine user_init_session
    use ModProcMH,      ONLY: iProc
    use ModLookupTable, ONLY: i_lookup_table
    use ModPhysics,     ONLY: AverageIonCharge
    use ModMultiFluid,  ONLY: MassIon_I
    use ModConst,       ONLY: cProtonMass, cBoltzmann
    implicit none
    
    character (len=*), parameter :: NameSub = 'user_init_session'
    !------------------------------------------------------------------------
    mu    =  MassIon_I(1)/(1 + AverageIonCharge)
    rstari = mu/(cBoltzmann/cProtonMass)  

    ! initialize the indexes for lookup tables
    iTableInitialState = i_lookup_table('RhoUzExtraEP(Z,Const)')
    iTablePPerE        = i_lookup_table('pPerE(rho,e/rho)')
    iTableCvTe         = i_lookup_table('CvTe(rho,p/rho)')
    iTablePERhoT       = i_lookup_table('pe(rho,T)')
    iTableChianti      = i_lookup_table('prl(T,Const)')
    iTableOpacity      = i_lookup_table('Opacity(rho,T)')
    if(iProc==0) write(*,*) NameSub, &
         'iTableInitialState, PPerE, CvTe, PressureEn , Chianti, opacity = ', &
         iTableInitialState, iTablePPerE, iTableCvTe, iTablePERhoT, &
         iTableChianti, iTableOpacity

    if(.not.allocated(Temperature_GB)) &
         allocate(Temperature_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))

  end subroutine user_init_session
  !==========================================================================
  subroutine user_set_ICs
    use ModMain,       ONLY: ProcTest, unusedBLK
    use ModProcMH,     ONLY: iProc
    use ModAdvance,    ONLY: State_VGB
    use ModVarIndexes
    implicit none

    logical :: oktest, oktest_me
    !------------------------------------------------------------------------
    
    if((iProc==PROCtest) .and. (GlobalBLK==1))then
       write(*,*)'Initializing Flux Emergence problem'
       write(*,*)'Parameters:'
       write(*,*)'iProc = ',iProc
       write(*,*) 'InitialBx = ',InitialBx,'InitialBy = ',InitialBy, &
            'InitialBz = ',InitialBz
       write(*,*)'InitialDensity = ',InitialDensity,'InitialTemperature = ', &
            InitialTemperature, ' in CGS units'
       call set_oktest('user_set_ICs',oktest,oktest_me)
    else
       oktest=.false.; oktest_me=.false.
    end if
    
    if (unusedBLK(GlobalBLK)) return
    if(UseUniformInitialState)then
       call set_uniform_ICs
       if((iProc==0).and.(GlobalBLK==1))then
          write(*,*) '------------------------------------------------------'
          write(*,*) 'Set up an uniform initial atmosphere with parameters:'
          write(*,*) 'InitialDensity = ',InitialDensity, &
               ', InitialTemperature = ',InitialTemperature, ' in CGS units'
          write(*,*) '------------------------------------------------------'
       end if
    else
       call set_perturbed_ICs
       if((iProc==0).and.(GlobalBLK==1))then
          write(*,*) '------------------------------------------------------'
          write(*,*) '            start with a relaxed atmosphere           '
          write(*,*) '------------------------------------------------------'
       end if
    end if
  contains
    !========================================================================
    subroutine set_uniform_ICs
      use ModLookupTable, ONLY: interpolate_lookup_table
      use ModPhysics,     ONLY: Si2No_V, UnitRho_, UnitP_, UnitN_, inv_gm1, &
           UnitEnergyDens_, UnitTemperature_, cRadiationNo
      implicit none
      
      real    :: pNo, Value_V(1:2), ExtraEintNo, InitialDensitySi
      integer :: i,j,k
      !---------------------------------------------------------------------
      ! calculate the Extra Internal Energy and Pressure (SI)
      ! from given initial conditions: InitialDensity, InitialTemperature(CGS)
      InitialDensitySi = InitialDensity*1e3
      if(iTablePERhoT>0)then
         call interpolate_lookup_table(iTablePERhoT, InitialDensitySi, &
              InitialTemperature, Value_V, DoExtrapolate = .false.)
         pNo         = Value_V(1)*Si2No_V(UnitP_)
         ExtraEintNo = (Value_V(2)-Value_V(1)*inv_gm1)*Si2No_V(UnitEnergyDens_)
      else
         pNo         = InitialDensitySi*InitialTemperature/rstari&
              *Si2No_V(UnitP_)
         ExtraEintNo = 0.
      end if

      ! Set the initial condition for a uniform atmosphere
      do i=1, nI; do j=1, nJ; do k=1, nK
         State_VGB(rhoUx_ ,i,j,k,GlobalBLK) = 0.
         State_VGB(rhoUy_ ,i,j,k,GlobalBLK) = 0.
         State_VGB(rhoUz_ ,i,j,k,GlobalBLK) = 0.
         State_VGB(Bx_    ,i,j,k,GlobalBLK) = 0.
         State_VGB(By_    ,i,j,k,GlobalBLK) = 0.
         State_VGB(Bz_    ,i,j,k,GlobalBLK) = 0.
         State_VGB(rho_   ,i,j,k,GlobalBLK) = &
              InitialDensitySi*Si2No_V(UnitRho_)
         State_VGB(Erad_  ,i,j,k,GlobalBLK) = &
              cRadiationNo*(InitialTemperature*Si2No_V(UnitTemperature_))**4
         State_VGB(P_     ,i,j,k,GlobalBlk) = pNo
         State_VGB(ExtraEInt_,i,j,k,GlobalBlk) = &
              ExtraEintNo
         Temperature_GB(i,j,k,GlobalBLK)       = InitialTemperature
      end do; end do ; end do
      
    end subroutine set_uniform_ICs
    !========================================================================
    subroutine set_perturbed_ICs
      use ModPhysics,     ONLY: Si2No_V, UnitRho_, UnitU_, UnitEnergyDens_,&
           UnitP_, UnitX_, No2Si_V
      use ModGeometry,    ONLY: z_BLK
      use ModLookupTable, ONLY: interpolate_lookup_table
      implicit none

      real    :: RandomChange, InitialState(1:4), TeSi, Const, &
           InitialRho, InitialUz, InitialExtraE, InitialP
      integer :: i, j, k
      !---------------------------------------------------------------------
      Const = 2.5
      do k =1, nK    
         ! interpolate the tabular data of reference initial state and get
         ! an relaxed initial state
         call interpolate_lookup_table(iTableInitialState, &
              z_BLK(1,1,k,GlobalBLK)*No2Si_V(UnitX_), Const, InitialState, &
              DoExtrapolate = .false.)
         InitialRho    = InitialState(1)
         InitialUz     = InitialState(2)
         InitialExtraE = InitialState(3)
         InitialP      = InitialState(4)
         State_VGB(rho_,:,:,k,GlobalBLK)          = &
              InitialRho*Si2No_V(UnitRho_)
         State_VGB(rhoUx_:rhoUy_,:,:,k,GlobalBLK) = 0.
         State_VGB(rhoUz_,:,:,k,GlobalBLK)        = &
              InitialRho*InitialUz*Si2No_V(UnitRho_)*Si2No_V(UnitU_)
         State_VGB(Bx_:Erad_ ,:,:,k,GlobalBLK)    = 0.
         State_VGB(ExtraEint_,:,:,k,GlobalBLK)    = &
              InitialExtraE*Si2No_V(UnitEnergyDens_)
         State_VGB(p_,:,:,k,GlobalBLK)            = &
              InitialP*Si2No_V(UnitP_)
         call user_material_properties(State_VGB(:,1,1,k,GlobalBLK), &
              TeOut = TeSi)
         Temperature_GB(:,:,k,GlobalBLK) = TeSi
         ! Add random perturbation to energy and pressure values of 
         ! cells below the photosphere height
         if(z_BLK(4,4,k,GlobalBlk) .lt. -1000.)then
            do i=1, nI; do j=1,nJ
               call random_number(RandomChange)
               RandomChange = (RandomChange-0.5)*2
               State_VGB(ExtraEint_,i,j,k,GlobalBLK) =          &
                    State_VGB(ExtraEint_,i,j,k,GlobalBLK) +     &
                    1.e3*RandomChange*Si2No_V(UnitEnergyDens_)
               State_VGB(p_,i,j,k,GlobalBLK) =                  &
                    State_VGB(p_,i,j,k,GlobalBLK) +             &
                    1.e3*RandomChange*Si2No_V(UnitP_)
               call user_material_properties(State_VGB(:,i,j,k,GlobalBLK), &
                    TeOut = Temperature_GB(i,j,k,GlobalBLK))
            end do; end do
         end if
      end do
    end subroutine set_perturbed_ICs
    
  end subroutine user_set_ICs
  
  !=========================================================================
  
  subroutine user_initial_perturbation
    use ModEnergy,   only: calc_energy
    use ModMain,     ONLY: unusedBlk, nBlockMax
    use ModGeometry, ONLY: z_BLK,y_BLK,x_BLK
    use ModAdvance,  ONLY: State_VGB
    use ModPhysics
    use ModVarIndexes
    implicit none
    
    integer:: iBlock, i,j,k
    real :: dp_ratio,prof,rsq,rasq,EInternal, PressureSi
    
    character (len=*), parameter :: NameSub = 'user_initial_perturbation'
    !-----------------------------------------------------------------------   
    rasq = ra_rope*ra_rope
    do iBlock = 1, nBlockMax
       if( UnusedBlk(iBlock) ) CYCLE
       do k=1,nK
          ! Add initial magnetic field 
          if(z_BLK(4,4,k,iBlock) .lt. z_photo)then
             prof = (1 - tanh(z_BLK(4,4,k,iBlock) - z_photo))/2.
             State_VGB(Bx_,:,:,k,iBlock) = &
                  State_VGB(Bx_,:,:,k,iBlock) + &
                  InitialBx*1.e-4*Si2No_V(UnitB_)*prof
             State_VGB(By_,:,:,k,iBlock) = &
                  State_VGB(By_,:,:,k,iBlock) + &
                  InitialBy*1.e-4*Si2No_V(UnitB_)*prof
          end if
          State_VGB(Bz_,:,:,k,iBlock) = &
               State_VGB(Bz_,:,:,k,iBlock) + InitialBz*1.e-4*Si2No_V(UnitB_)
          !\
          ! If UseRope, Add flux rope in
          ! Set negative pressure to 1.e-10
          !/
          if(UseRope)then
             do i=1,nI ; do j = 1, nJ
                rsq = (y_BLK(i,j,k,iBlock) - x2c_rope)**2 + &
                     (z_BLK(i,j,k,iBlock) - x3c_rope)**2
                if(rsq < rasq**1.5e1)then
                   prof = b0_rope*exp(-rsq/rasq)
                else
                   prof = 0.
                end if
                dp_ratio =  0.5*prof*prof*(-1.+0.5*qfac_rope*qfac_rope* &
                     (1. - 2.*rsq/rasq))/State_VGB(p_,i,j,k,iBlock)
                State_VGB(Bx_,i,j,k,iBlock) = &
                     State_VGB(Bx_,i,j,k,iBlock) + prof
                State_VGB(By_,i,j,k,iBlock) = &
                     State_VGB(By_,i,j,k,iBlock) &
                     - prof*qfac_rope*(z_BLK(i,j,k,iBlock) - x3c_rope)/ra_rope
                State_VGB(Bz_,i,j,k,iBlock) = &
                     State_VGB(Bz_,i,j,k,iBlock) + &
                     prof*qfac_rope*(y_BLK(i,j,k,iBlock) - x2c_rope)/ra_rope
                State_VGB(P_,i,j,k,iBlock) = &
                     State_VGB(p_,i,j,k,iBlock)*(1.+dp_ratio)
                State_VGB(rho_,i,j,k,iBlock) = &
                     State_VGB(rho_,i,j,k,iBlock)*(1. + &
                     exp(-(x_BLK(i,j,k,iBlock)/lamb_rope)**2)*&
                     dp_ratio*buoyancy_rope)
                if(State_VGB(p_,i,j,k,iBlock).le.0.)&
                     State_VGB(p_,i,j,k,iBlock)=1.e-10
             end do; end do
          end if
       end do
       call calc_energy(-1,nI+2,-1,nJ+2,-1,nK+2,iBlock,1,1)
    end do
  end subroutine user_initial_perturbation
  
  !==========================================================================
  subroutine user_set_outerbcs(iBlock,iSide, TypeBc, IsFound)
    use ModVarIndexes, ONLY: rho_, rhoUz_, Bz_, p_, Erad_, ExtraEInt_
    use ModPhysics,    ONLY: No2Si_V, UnitX_, UnitP_, UnitEnergyDens_, &
         inv_gm1, UnitTemperature_, unitRho_, cRadiationNo, cLight, Si2No_V
    use ModGeometry,   ONLY: dz_BLK
    use ModAdvance,    ONLY: State_VGB
    use ModImplicit,   ONLY: iEradImpl, StateSemi_VGB

    integer,          intent(in)  :: iBlock, iSide
    character(len=20),intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    integer :: i,j
    real :: EInternalSI1p, EInternalSI2g,PressureSi, EinternalSI1g, DeltaTemp, &
         Trad, Fincident, Coef, OpacityPlanckSi_W(1), OpacityPlanck

    character (len=*), parameter :: NameSub = 'user_set_outerbcs'
    !------------------------------------------------------------------------

    select case(iSide)
    case(5)
       select case(TypeBc)
       case('TempGrad')
          ! set temperature gradient at the bottom boundary
          State_VGB(rho_:p_,:,:,0,iBlock) = &
               State_VGB(rho_:p_,:,:,1,iBlock)
          State_VGB(rho_:p_,:,:,-1,iBlock) = &
               State_VGB(rho_:p_,:,:,1,iBlock)
          State_VGB(rhoUz_,:,:,0,iBlock) = &
               -State_VGB(rhoUz_,:,:,1,iBlock)
          State_VGB(rhoUz_,:,:,-1,iBlock) = &
               -State_VGB(rhoUz_,:,:,1,iBlock)
          State_VGB(Bz_,:,:,0,iBlock) = &
               -State_VGB(Bz_,:,:,1,iBlock)
          State_VGB(Bz_,:,:,-1,iBlock) = &
               -State_VGB(Bz_,:,:,1,iBlock)

          DeltaTemp = TemperatureGradient*dz_BLK(iBlock)*No2Si_V(UnitX_)/&
               No2Si_V(UnitTemperature_)
          do i=-1,nI+2; do j=-1,nJ+2
             EInternalSI1p = No2Si_V(UnitEnergyDens_)*    &
                  (inv_gm1*State_VGB(p_,i,j,1,iBlock) + &
                  State_VGB(ExtraEInt_,i,j,1,iBlock))
             EinternalSI2g = EInternalSI1p + State_VGB(rho_,i,j,0,iBlock)* &
                  DeltaTemp*inv_gm1*No2Si_V(UnitRho_)/rstari
             call user_material_properties(State_VGB(:,i,j,0,iBlock), &
                  EInternalIn=EInternalSI2g, PressureOut=PressureSi)
             State_VGB(p_,i,j,0,iBlock) = PressureSi/No2Si_V(UnitP_)
             State_VGB(ExtraEInt_,i,j,0,iBlock) = (EInternalSi2g - PressureSi&
                  *inv_gm1)/No2Si_V(UnitEnergyDens_)
             EinternalSI1g = EInternalSI2g + State_VGB(rho_,i,j,-1,iBlock)* &
                  DeltaTemp*inv_gm1*No2Si_V(UnitRho_)/rstari
             call user_material_properties(State_VGB(:,i,j,-1,iBlock), &
                  EInternalIn=EInternalSI1g, PressureOut=PressureSi)
             State_VGB(p_,i,j,-1,iBlock) = PressureSi/No2Si_V(UnitP_)
             State_VGB(ExtraEInt_,i,j,-1,iBlock) = (EInternalSI1g - PressureSi&
                  *inv_gm1)/No2Si_V(UnitEnergyDens_)
             Trad = (State_VGB(Erad_,i,j,1,iBlock)/cRadiationNo)**0.25
             State_VGB(Erad_,i,j,0,iBlock) = cRadiationNo*(Trad + DeltaTemp)**4
             State_VGB(Erad_,i,j,-1,iBlock) = cRadiationNo*(Trad + 2*DeltaTemp)**4
          end do; end do
          IsFound = .true.
       case('usersemi')
          DeltaTemp = TemperatureGradient*dz_BLK(iBlock)*No2Si_V(UnitX_)/&
               No2Si_V(UnitTemperature_)
          do i = -1, nI+1; do j = -1, nJ+1
             call user_material_properties(State_VGB(:,i,j,0,iBlock), &
                  i, j, 0, iBlock, OpacityPlanckOut_W = OpacityPlanckSi_W)
             OpacityPlanck = OpacityPlanckSi_W(1)/Si2No_V(UnitX_)
             Coef = 2.0/(3.0*OpacityPlanck*dz_BLK(iBlock))
             Trad = (State_VGB(Erad_,i,j,1,iBlock)/cRadiationNo)**0.25 + 0.5*DeltaTemp
             Fincident = 16./3.*cRadiationNo*TemperatureGradient*Si2No_V(UnitTemperature_)&
                  /Si2No_V(UnitX_)*Trad**3/OpacityPlanck
             StateSemi_VGB(iEradImpl, i, j, nK+1, iBlock) = (Fincident*4./cLight + &
                  StateSemi_VGB(iEradImpl, i, j, nK, iBlock)*(Coef - 0.5))/(Coef + 0.5)
          end do; end do
          IsFound = .true.
       case('usersemilinear')
          do i = -1, nI+1; do j = -1, nJ+1
             call user_material_properties(State_VGB(:,i,j,0,iBlock), &
                  i, j, 0, iBlock, OpacityPlanckOut_W = OpacityPlanckSi_W)
             OpacityPlanck = OpacityPlanckSi_W(1)/Si2No_V(UnitX_)
             Coef = 2.0/(3.0*OpacityPlanck*dz_BLK(iBlock))
             StateSemi_VGB(iEradImpl, i, j, nK+1, iBlock) = &
                  StateSemi_VGB(iEradImpl, i, j, nK, iBlock)*(Coef - 0.5)/(Coef + 0.5)
          end do; end do
          IsFound = .true.
       end select
    case(6)
       select case(TypeBc)
          ! open upper boundary condition                         
       case('open')
          State_VGB(:,:,:,nK+1,iBlock) = State_VGB(:,:,:,nK,iBlock)
          State_VGB(:,:,:,nK+2,iBlock) = State_VGB(:,:,:,nK,iBlock)
          IsFound = .true.
          ! closed upper boundary condition, not allowing upward motions
       case('no_inflow')
          State_VGB(:,:,:,nK+1,iBlock) = State_VGB(:,:,:,nK,iBlock)
          State_VGB(:,:,:,nK+2,iBlock) = State_VGB(:,:,:,nK,iBlock)
          State_VGB(rhoUz_,:,:,nK+1,iBlock) = &
               abs(State_VGB(rhoUz_,:,:,nK,iBlock))
          State_VGB(rhoUz_,:,:,nK+2,iBlock) = &
               abs(State_VGB(rhoUz_,:,:,nK,iBlock))
          IsFound = .true.
       case('usersemi','usersemilinear')
          do i = -1, nI+1; do j = -1, nJ+1
             call user_material_properties(State_VGB(:,i,j,nK+1,iBlock), &
                  i, j, nK+1, iBlock, OpacityPlanckOut_W = OpacityPlanckSi_W)
             OpacityPlanck = OpacityPlanckSi_W(1)/Si2No_V(UnitX_)
             Coef = 2.0/(3.0*OpacityPlanck*dz_BLK(iBlock))
             StateSemi_VGB(iEradImpl, i, j, nK+1, iBlock) = &
                  StateSemi_VGB(iEradImpl, i, j, nK, iBlock)*(Coef - 0.5)/(Coef + 0.5)
          end do; end do
          IsFound = .true.
       end select
    end select

  end subroutine user_set_outerbcs
  
  !==========================================================================

  subroutine user_calc_sources
    use ModAdvance,     ONLY: Source_VC, State_VGB
    use ModGeometry,    ONLY: z_BLK
    use ModVarIndexes,  ONLY: Energy_, rhoUz_
    implicit none

    integer :: i, j, k
    real    :: RadiativeCooling, EInternalSource, rhoUzSource,&
         DampingRhoUz, DampingEnergy
    
    character (len=*), parameter :: NameSub = 'user_calc_sources'
    !------------------------------------------------------------------------
    RadiativeCooling = 0.
    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       call user_material_properties(State_VGB(:,i,j,k,GlobalBLK), &
            TeOut = Temperature_GB(i,j,k,GlobalBLK))
       if(UseVerticalDamping) call get_vertical_damping( &
            State_VGB(:,i,j,k,GlobalBLK), &
            z_BLK(i,j,k,GlobalBLK),DampingRhoUz, DampingEnergy)
       if(UseThinRadiation) call get_radiative_cooling(&
            State_VGB(:,i,j,k,GlobalBLK), Temperature_GB(i,j,k,GlobalBLK), &
            z_BLK(i,j,k,GlobalBLK), RadiativeCooling)

       EInternalSource = RadiativeCooling
       rhoUzSource     = DampingRhoUz

       Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + &
            EinternalSource + DampingEnergy
       Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k) + rhoUzSource
    end do; end do; end do
  end subroutine user_calc_sources
  
  !==========================================================================
  !----------------------------------------------!
  !  subroutines containing energy source terms  !            
  !----------------------------------------------!
  subroutine get_vertical_damping(State_V, Z_V, DampingRhoUz, DampingEnergy)
    use ModGeometry,    ONLY: z_BLK
    use ModPhysics,     ONLY: Si2No_V, UnitT_
    use ModVarIndexes,  ONLY: rhoUz_, rho_, nVar
    implicit none
    
    real, intent(in) :: State_V(nVar), Z_V
    real, intent(out):: DampingRhoUz, DampingEnergy
    !-----------------------------------------------------------------------  
    if(UseUniformInitialState)then
       DampingRhoUz  = -State_V(rhoUz_)/ &
            (TimeVerticalDamping*Si2No_V(UnitT_))
       DampingEnergy = -State_V(rhoUz_)**2/State_V(rho_)/&
            (TimeVerticalDamping*Si2No_V(UnitT_))
    elseif(Z_V .gt. z_photo)then
       DampingRhoUz  = -State_V(rhoUz_)/ &
            (TimeVerticalDamping*Si2No_V(UnitT_))
       DampingEnergy = -State_V(rhoUz_)**2/State_V(rho_)/&
            (TimeVerticalDamping*Si2No_V(UnitT_))
    else
       DampingRhoUz  = 0.
       DampingEnergy = 0.
    end if
  end subroutine get_vertical_damping
  !=========================================================================
  subroutine get_radiative_cooling(State_V, TeSi, Z_V, RadiativeCooling)
    use ModVarIndexes, ONLY: rho_, p_, nVar
    use ModPhysics,    ONLY: UnitRho_, UnitEnergyDens_, Si2No_V, UnitT_
    use ModLookupTable,ONLY: interpolate_lookup_table
    use ModConst,      ONLY: cProtonMass
    implicit none
    
    real, intent(in) :: State_V(1:nVar)
    real, intent(in) :: TeSi
    real, intent(in) :: Z_V
    real, intent(out):: RadiativeCooling

    real, parameter :: RadiationCutoff = - 5e8, atrl = 1e4, Const = 2.5
    real :: Fraction = 0.0, CoolingFunctionCgs = 0.0, &
         csw, MassDensCgs, NumberDensCgs
    real :: Chianti(1:1)
    !------------------------------------------------------------------------
    csw = atrl/RhoThinCutoff
    MassDensCgs     = State_V(rho_)/Si2No_V(UnitRho_)*1.e-3
    NumberDensCgs   = MassDensCgs/(mu*cProtonMass*1.e3)

    ! calculate the thin radiative loss above the z_phto height
    if (z_V > z_photo) then
       ! Smoothing function is 1 if rho<RhoThinCutoff , 0 if not
       Fraction = 0.5 - 0.5*&
            tanh(atrl*(MassDensCgs/RhoThinCutoff - 1.))
       ! Calculate the cooling function culve dependent on temperature
       if (TeSi <= 3.e5 ) Fraction = Fraction*(TeSi/3.e5)**2.5
       if (TeSi <= 8e+03)then
          CoolingFunctionCgs = (1.0606e-06*TeSi)**11.7
       elseif(TeSi <= 1e4)then
          CoolingFunctionCgs = (1.397e-08*TeSi)**6.15
       else
          call interpolate_lookup_table(iTableChianti, TeSi, Const, &
               Chianti, DoExtrapolate = .false.)
          CoolingFunctionCgs = Chianti(1)
       end if
       ! thin radiative cooling = -\CoolinFunction * n_{e}*n_{p}
       RadiativeCooling = -Fraction*NumberDensCgs**2*CoolingFunctionCgs
       if(RadiativeCooling/MassDensCgs <= RadiationCutoff) &
            RadiativeCooling = RadiationCutoff*MassDensCgs
       RadiativeCooling = RadiativeCooling*0.1*Si2No_V(UnitEnergyDens_)/&
            Si2No_V(UnitT_)
    else 
       RadiativeCooling = 0.
    end if
    
  end subroutine get_radiative_cooling

  !==========================================================================
  
  subroutine user_update_states(iStage,iBlock)
    use ModVarIndexes, ONLY: rho_, p_, ExtraEInt_
    use ModAdvance,    ONLY: State_VGB
    use ModPhysics,    ONLY: Si2No_V, No2Si_V, UnitRho_, UnitP_, &
         UnitEnergyDens_, cProtonMass, inv_gm1
    use ModEnergy,     ONLY: calc_energy_cell
    implicit none
    
    integer, intent(in) :: iStage,iBlock
    integer:: i,j,k
    real   :: MassDensFloor, EnergyFloor, EinternalSi, PressureSi
    character (len=*), parameter :: NameSub = 'user_update_states'
    !--------------------------------------------------------------------
    ! Define the minimum mass density and energy value
    ! Corresponding to ~NumberDensFloor and 200 K plasma
    MassDensFloor  = NumberDensFloor*1e6*cProtonMass*Si2No_V(UnitRho_)
    EnergyFloor    = MassDensFloor*No2Si_V(UnitRho_)/rstari*200.

    call update_states_MHD(iStage,iBlock)
    do k = -1,nK+2; do j=-1,nJ+2; do i=-1,nI+2
       ! check if the density or pressure is below the minimum value
       if(State_VGB(rho_,i,j,k,iBlock) < MassDensFloor) &
            State_VGB(rho_,i,j,k,iBlock) = MassDensFloor
       ! Total internal energy, ExtraEInt + P/(\gamma -1),
       EInternalSi = (inv_gm1*State_VGB(P_,i,j,k,iBlock) + &
            State_VGB(ExtraEInt_,i,j,k,iBlock))*No2Si_V(UnitEnergyDens_)
       if(EInternalSi < EnergyFloor) EInternalSi = EnergyFloor
       ! get pressure and extra internal energy from the EOS table
       call user_material_properties(State_VGB(:,i,j,k,iBlock), &
            EInternalIn=EInternalSi, PressureOut=PressureSi)
       State_VGB(P_,i,j,k,iBlock) = PressureSi*Si2No_V(UnitP_)
       State_VGB(ExtraEInt_,i,j,k,iBlock) = Si2No_V(UnitEnergyDens_)*&
            (EInternalSI - PressureSi*inv_gm1)
    end do;end do; end do
    ! calculate the total energy
    call calc_energy_cell(iBlock)
  end subroutine user_update_states

  !===========================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G,PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)
    implicit none
    
    integer,          intent(in) :: iBlock
    character(len=*), intent(in) :: NameVar
    logical,          intent(in) :: IsDimensional
    real,             intent(out):: PlotVar_G(-1:nI+2, -1:nJ+2, -1:nK+2)
    real,             intent(out):: PlotVarBody
    logical,          intent(out):: UsePlotVarBody
    character(len=*), intent(out):: NameTecVar
    character(len=*), intent(out):: NameTecUnit
    character(len=*), intent(out):: NameIdlUnit
    logical,          intent(out):: IsFound

    character (len=*), parameter :: Name='user_set_plot_var'
    !-----------------------------------------------------------------
    UsePlotVarBody = .true.
    PlotVarBody    = 0.0
    IsFound        = .true.
    select case(NameVar)
    case('tempe')
       NameTecVar  = 'T'
       NameTecUnit = '[K]'
       NameIdlUnit = 'K'
       PlotVar_G   = Temperature_GB(:,:,:,iBlock)
    case default 
       IsFound = .false.
       call stop_mpi(Name//': unknown plot variables = '//NameVar)
    end select
  end subroutine user_set_plot_var
  
  !========================================================================
  
  subroutine user_material_properties(State_V, i,j,k,iBlock,iDir, &
       EinternalIn, TeIn, NatomicOut, AverageIonChargeOut, &
       EinternalOut, TeOut, PressureOut,   &
       CvOut, GammaOut, HeatCondOut, IonHeatCondOut, TeTiRelaxOut, &
       OpacityPlanckOut_W, OpacityRosselandOut_W, PlanckOut_W)

    use ModLookupTable,ONLY: interpolate_lookup_table
    use ModPhysics,    ONLY: No2Si_V, UnitEnergyDens_, UnitP_, &
         UnitRho_, inv_gm1, gm1, gamma0
    use ModVarIndexes, ONLY: nVar, Rho_, p_, ExtraEInt_
    use ModAdvance,    ONLY: nWave

    !--------------------------------------------------------------------------
    ! The State_V vector is in normalized units
    real, intent(in) :: State_V(nVar)
    integer, optional, intent(in) :: i, j, k, iBlock, iDir
    real, optional, intent(in)  :: EinternalIn             ! [J/m^3]
    real, optional, intent(in)  :: TeIn                    ! [K]
    real, optional, intent(out) :: NatomicOut              ! [1/m^3]
    real, optional, intent(out) :: AverageIonChargeOut     ! dimensionless
    real, optional, intent(out) :: EinternalOut            ! [J/m^3]
    real, optional, intent(out) :: TeOut                   ! [K]
    real, optional, intent(out) :: PressureOut             ! [Pa]   
    real, optional, intent(out) :: CvOut                   ! [J/(K*m^3)]  
    real, optional, intent(out) :: GammaOut
    real, optional, intent(out) :: HeatCondOut             ! [Jm^2/(Ks)]   
    real, optional, intent(out) :: IonHeatCondOut          ! [J/(m*K*s)]
    real, optional, intent(out) :: TeTiRelaxOut            ! [1/s]  
    real, optional, intent(out) :: OpacityPlanckOut_W(nWave)      ! [1/m] 
    real, optional, intent(out) :: OpacityRosselandOut_W(nWave)   ! [1/m] 
    real, optional, intent(out) :: PlanckOut_W(nWave)      ! [J/m^3] 

    real    :: pSi, enSi,RhoSi, TeSi, pPerE(1:1), PressureEnSi(2)
    real    :: Value_V(2)
    character (len=*), parameter :: NameSub = 'user_material_properties'
    !-------------------------------------------------------------------------
    ! Density, transformed to SI
    RhoSi = No2Si_V(UnitRho_)*State_V(Rho_)
    
    ! Calculate Pressure and Internal Energy from EOS table or Ideal EOS
    if(present(EinternalIn))then
       if(iTablePPerE > 0)then
          call interpolate_lookup_table(iTablePPerE, RhoSi, &
               EinternalIn/RhoSi, pPerE, DoExtrapolate = .false.)
          pSi = pPerE(1)*EinternalIn
       else
          pSi = EinternalIn*gm1
       end if
       enSi = EinternalIn
    elseif(present(TeIn))then
       TeSi = TeIn
       call interpolate_lookup_table(iTablePERhoT,RhoSi,TeSi, &
            Value_V, DoExtrapolate = .false.)
       pSi = Value_V(1)
       enSi= Value_V(2)
    else
       pSi = State_V(p_)*No2Si_V(UnitP_)
       enSi= pSi*inv_gm1 + State_V(ExtraEInt_)*No2Si_V(UnitEnergyDens_)
    end if
    if(present(PressureOut)) PressureOut = pSi
    if(present(EInternalOut)) EInternalOut = enSi

    ! Calculate Temperature from EOS table or Ideal EOS
    if(present(TeOut) .or. present(CvOut) .or. present(HeatCondOut).or.&
         present(OpacityPlanckOut_W).or.present(OpacityRosselandOut_W))then
       if(iTableCvTe > 0)then
          call interpolate_lookup_table(iTableCvTe, RhoSi, pSi/RhoSi, &
               Value_V, DoExtrapolate = .false.)
          TeSi = Value_V(2)     
          if(present(CvOut)) CvOut =  Value_V(1)
       else
          TeSi = pSi/RhoSi*rstari
       endif
    endif
    if(present(TeOut)) TeOut = TeSi

    if(present(OpacityPlanckOut_W).or.present(OpacityRosselandOut_W))then
       call interpolate_lookup_table(iTableOpacity, RhoSi, TeSi, &
            Value_V, DoExtrapolate = .false.)
    end if
    if(present(OpacityPlanckOut_W))OpacityPlanckOut_W = Value_V(1)*&
         (0.5 + 0.5*tanh((RhoSi - 1e-4)*1e5))
    if(present(OpacityRosselandOut_W))OpacityRosselandOut_W = Value_V(2) + &
         (0.5 - 0.5*tanh((RhoSi - 1e-4)*1e5))*1e3

    if(present(HeatCondOut))&
         HeatCondOut = 1.e-11*(max(TeSi,3e5))**2.5

  end subroutine user_material_properties
  
end module ModUser
