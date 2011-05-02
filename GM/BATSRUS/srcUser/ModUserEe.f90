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
  logical :: UseCoronalField = .false.
  logical :: UseUniformInitialState = .false.
  logical :: UseUniformT = .false.
  logical :: UseEnergyPert = .false.
  real    :: InitialDensity, InitialTemperature, InitialBx, &
       InitialBy, InitialBz
  real    :: RhoThinCutoff, NumberDensFloor, TimeVerticalDamping
  real    :: z_photo, TemperatureGradient
  real    :: BotDensity, BotPressure, BotExtraE
  !rstari = 0.594354e-3/8.31, mu = 0.594354 set in init_session
  real,parameter ::  mu = 0.594354, rstari = 0.594354e-3/8.31

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
  ! Indexes of EOS tables, CHIANTI table, initial relaxed reference state
  !/
  integer :: iTableEOS = -1, &
       iTableChianti = -1, iTableInitialState = -1, iTableGamma = -1
  real, allocatable:: srcthin_GB(:,:,:,:)

contains

  !=========================================================================
  subroutine user_read_inputs
    use ModMain,      ONLY: lverbose
    use ModProcMH,    ONLY: iProc
    use ModReadParam, ONLY: read_var, read_line, read_command
    use ModIO,        ONLY: write_prefix, write_myname, iUnitOut
    use ModCoronalHeating, ONLY: DtUpdateFlux, UnsignedFluxHeight

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
          call read_var('UseCoronalField',UseCoronalField)
          call read_var('RhoThinCutoff',RhoThinCutoff)
          call read_Var('UseVerticalDamping',UseVerticalDamping)
          call read_var('TimeVerticalDamping', TimeVerticalDamping)
          call read_var('TemperatureGradient',TemperatureGradient)
          call read_var('DtUptateFlux',DtUpdateFlux)
          call read_var('UseUniformInitalState', UseUniformInitialState)
          call read_var('UseUniformT', UseUniformT)
          call read_var('UseEnergyPert', UseEnergyPert)
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
    use ModLookupTable, ONLY: i_lookup_table, interpolate_lookup_table
    use ModPhysics,     ONLY: No2Si_V, UnitX_, Si2No_V, UnitRho_, UnitP_, &
         UnitEnergyDens_
    use ModMultiFluid,  ONLY: MassIon_I
    use ModConst,       ONLY: cProtonMass, cBoltzmann
    use ModGeometry,    ONLY: z1

    real :: InitialState(1:4)

    character (len=*), parameter :: NameSub = 'user_init_session'
    !------------------------------------------------------------------------
    !Gbody = -1 
    !mu    =  MassIon_I(1)/(1 + AverageIonCharge)
    !rstari = mu/(cBoltzmann/cProtonMass)  

    ! initialize the indexes for lookup tables
    iTableInitialState = i_lookup_table('RhoUzExtraEP(Z,Const)')
    iTableEOS          = i_lookup_table('peent(T,rho)')
    iTableChianti      = i_lookup_table('prl(T,Const)')
    iTableGamma        = i_lookup_table('Gamma(T,rho)')

    if(iProc==0) write(*,*) NameSub, &
         'iTableInitialState, EOS , Chianti, gamma = ', &
         iTableInitialState, iTableEOS, iTableChianti, iTableGamma

    call interpolate_lookup_table(iTableInitialState, &
         z1*No2Si_V(UnitX_), 2.5, InitialState, &
         DoExtrapolate = .false.)
    BotDensity  = InitialState(1)*Si2No_V(UnitRho_)
    BotExtraE   = InitialState(3)*Si2No_V(UnitEnergyDens_)
    BotPressure = InitialState(4)*Si2No_V(UnitP_)

    if(.not.allocated(srcthin_GB)) &
         allocate(srcthin_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))

  end subroutine user_init_session
  !==========================================================================
  subroutine user_set_ICs
    use ModMain,       ONLY: ProcTest, unusedBLK
    use ModProcMH,     ONLY: iProc
    use ModAdvance,    ONLY: State_VGB
    use ModVarIndexes

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

    if (unusedBLK(GlobalBLK)) RETURN

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
      use ModPhysics,     ONLY: Si2No_V, UnitRho_, UnitP_, inv_gm1, &
           UnitEnergyDens_, UnitX_, No2Si_V
      use ModGeometry,    ONLY: z_BLK

      real :: p, Value_V(1:3), ExtraEint, InitialDensitySi, Rho
      real :: g1, inv_g1, inv_g1m1, z0

      integer :: i,j,k
      !---------------------------------------------------------------------
      ! calculate the Extra Internal Energy and Pressure (SI)
      ! from given initial conditions: InitialDensity, InitialTemperature(CGS)
      InitialDensitySi = InitialDensity*1e3
      if(iTableEOS>0)then
         call interpolate_lookup_table(iTableEOS, InitialTemperature, &
              InitialDensitySi, Value_V, DoExtrapolate = .false.)
         p         = Value_V(1)*Si2No_V(UnitP_)
         ExtraEint = (Value_V(2)-Value_V(1)*inv_gm1)*Si2No_V(UnitEnergyDens_)
      else
         p         = InitialDensitySi*InitialTemperature/rstari*Si2No_V(UnitP_)
         ExtraEint = 0.
      end if
      Rho = InitialDensitySi*Si2No_V(UnitRho_)

      ! Set the initial condition for a uniform atmosphere
      do k=1,nK; do j=1,nJ; do i=1,nI
         State_VGB(rhoUx_ ,i,j,k,GlobalBLK) = 0.
         State_VGB(rhoUy_ ,i,j,k,GlobalBLK) = 0.
         State_VGB(rhoUz_ ,i,j,k,GlobalBLK) = 0.
         State_VGB(Bx_    ,i,j,k,GlobalBLK) = 0.
         State_VGB(By_    ,i,j,k,GlobalBLK) = 0.
         State_VGB(Bz_    ,i,j,k,GlobalBLK) = 0.
         State_VGB(rho_   ,i,j,k,GlobalBLK) = Rho              
         State_VGB(Erad_  ,i,j,k,GlobalBLK) = 0.
         !cRadiationNo*(InitialTemperature*Si2No_V(UnitTemperature_))**4
         if(UseUniformT)then
            State_VGB(P_,i,j,k,GlobalBlk) = p
            State_VGB(ExtraEInt_,i,j,k,GlobalBlk) = ExtraEint
         else
            if(iTableEOS < -2)then
               g1       = 1.7
               inv_g1   = 1.0/g1
               inv_g1m1 = 1.0/(g1-1.0)
               z0 = -20.0
               State_VGB(rho_,i,j,k,GlobalBLK) = Rho
               State_VGB(rho_,i,j,k,GlobalBLK) = Rho*(1 - &
                    (1-inv_g1)*2.73e2*1.67e-27/(1.38e-23*InitialTemperature)* &
                    (min(z_BLK(i,j,k,GlobalBLK), z0) - z0)&
                    *No2Si_V(UnitX_))**inv_g1m1
               State_VGB(p_,i,j,k,GlobalBLK) = p*(1 - &
                    (1-inv_g1)*2.73e2*1.67e-27/(1.38e-23*InitialTemperature)* &
                    (min(z_BLK(i,j,k,GlobalBLK), z0) - z0)*No2Si_V(UnitX_))** &
                    (1+inv_g1m1)
               State_VGB(ExtraEInt_,i,j,k,GlobalBlk) = 0.0
            end if
         end if
      end do; end do ; end do

    end subroutine set_uniform_ICs
    !========================================================================
    subroutine set_perturbed_ICs

      use ModPhysics,     ONLY: Si2No_V, UnitRho_, UnitU_, UnitEnergyDens_,&
           UnitP_, UnitX_, No2Si_V
      use ModGeometry,    ONLY: z_BLK
      use ModLookupTable, ONLY: interpolate_lookup_table

      real    :: InitialState(1:4), Const, &
           InitialRho, InitialUz, InitialExtraE, InitialP
      integer :: k
      !---------------------------------------------------------------------
      Const = 2.5
      do k = 1, nK
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
      end do
    end subroutine set_perturbed_ICs

  end subroutine user_set_ICs

  !=========================================================================

  subroutine user_initial_perturbation
    use ModEnergy,   only: calc_energy_cell
    use ModMain,     ONLY: unusedBlk, nBlockMax
    use ModGeometry, ONLY: z_BLK,y_BLK,x_BLK
    use ModAdvance,  ONLY: State_VGB
    use ModPhysics
    use ModVarIndexes

    integer:: iBlock, i,j,k
    real :: dp_ratio,prof,rsq,rasq,EInternal, RandomChange

    character (len=*), parameter :: NameSub = 'user_initial_perturbation'
    !-----------------------------------------------------------------------   
    rasq = ra_rope*ra_rope
    do iBlock = 1, nBlockMax
       if( UnusedBlk(iBlock) ) CYCLE
       do k=1,nK
          ! Add initial magnetic field 
          if(UseCoronalField)then
             if(z_BLK(4,4,k,iBlock) < z_photo)then
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
          end if
          ! Add random perturbation to energy and pressure values of 
          ! cells below the photosphere height                  
          if(UseEnergyPert.and.(z_BLK(4,4,k,iBlock) < z_photo ))then
             do j=1, nJ; do i=1,nI
                call random_number(RandomChange)
                RandomChange = (RandomChange-0.5)*2
                if(iTableEOS > 0)then
                   State_VGB(ExtraEint_,i,j,k,iBlock) =          &
                        State_VGB(ExtraEint_,i,j,k,iBlock)*(1.0 +  &
                        1.e-3*RandomChange)
                end if
                State_VGB(p_,i,j,k,iBlock) =                  &
                     State_VGB(p_,i,j,k,iBlock)*(1.0 +        &
                     1.e-3*RandomChange) 
                EInternal = (inv_gm1*State_VGB(p_,i,j,k,iBlock) + &
                     State_VGB(ExtraEint_,i,j,k,iBlock))* &
                     No2Si_V(UnitEnergyDens_)
             end do; end do
          end if
          !\
          ! If UseRope, Add flux rope in
          ! Set negative pressure to 1.e-10
          !/
          if(UseRope)then
             do j=1,nJ ; do i = 1, nI
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
                if(State_VGB(p_,i,j,k,iBlock)<=0.)&
                     State_VGB(p_,i,j,k,iBlock)=1.e-10
             end do; end do
          end if
       end do

       call calc_energy_cell(iBlock)
    end do

  end subroutine user_initial_perturbation

  !==========================================================================
  subroutine user_set_outerbcs(iBlock,iSide, TypeBc, IsFound)
    use ModVarIndexes!, ONLY: rho_, rhoUz_, Bz_, p_, Erad_, ExtraEInt_
    use ModPhysics,    ONLY: UnitEnergyDens_, inv_gm1, Si2No_V
    use ModGeometry,   ONLY: dz_BLK
    use ModAdvance,    ONLY: State_VGB

    integer,          intent(in)  :: iBlock, iSide
    character(len=20),intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    integer :: i, j, k
    real :: EinternalSi

    character (len=*), parameter :: NameSub = 'user_set_outerbcs'
    !------------------------------------------------------------------------

    select case(iSide)
    case(5)
       select case(TypeBc)
       case('fixvalue')
          State_VGB(rho_,:,:,-1:0,iBlock) = BotDensity
          State_VGB(rhoUx_:rhoUy_,:,:,-1:0,iBlock) = 0.0
          State_VGB(rhoUz_,:,:,0,iBlock) = -State_VGB(rhoUz_,:,:,1,iBlock)
          State_VGB(rhoUz_,:,:,-1,iBlock) = -State_VGB(rhoUz_,:,:,2,iBlock)
          State_VGB(Bx_:By_,:,:,-1:0,iBlock) = 0.0
          State_VGB(Bz_,:,:,0,iBlock) = State_VGB(Bz_,:,:,1,iBlock)
          State_VGB(Bz_,:,:,-1,iBlock) = State_VGB(Bz_,:,:,1,iBlock)
          State_VGB(p_,:,:,-1:0,iBlock) = BotPressure
          State_VGB(ExtraEint_,:,:,-1:0,iBlock) = BotExtraE
          IsFound = .true.
       case('forcebalance')
          ! the density at the boundary is fixed.
          ! pressure gradient at the boundary is set to balance the 
          ! gravitational force.
          State_VGB(rho_,:,:,-1:0,iBlock) = BotDensity
          State_VGB(rhoUx_:rhoUy_,:,:,-1:0,iBlock) = 0.0
          State_VGB(rhoUz_,:,:,0,iBlock) = -State_VGB(rhoUz_,:,:,1,iBlock)
          State_VGB(rhoUz_,:,:,-1,iBlock) = -State_VGB(rhoUz_,:,:,1,iBlock)
          State_VGB(Bx_:By_,:,:,-1:0,iBlock) = 0.0
          State_VGB(Bz_,:,:,0,iBlock) = State_VGB(Bz_,:,:,1,iBlock)
          State_VGB(Bz_,:,:,-1,iBlock) = State_VGB(Bz_,:,:,1,iBlock)
          State_VGB(p_,:,:,0,iBlock) =  State_VGB(p_,:,:,2,iBlock) + &
               State_VGB(rho_,:,:,1,iBlock)*2.0*dZ_BLK(iBlock)
          State_VGB(p_,:,:,-1,iBlock) = State_VGB(p_,:,:,1,iBlock) + &
               State_VGB(rho_,:,:,0,iBlock)*2.0*dZ_BLK(iBlock)
          do k = -1, 0; do j = -1, nJ+2; do i = -1, nI+2 
             call user_material_properties(State_VGB(:,i,j,k,iBlock), &
                  EinternalOut = EinternalSi)
             State_VGB(ExtraEint_,i,j,k,iBlock) = EinternalSi* &
                  Si2No_V(UnitEnergyDens_) - State_VGB(p_,i,j,k,iBlock)*inv_gm1
          end do; end do; end do
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
       end select
    end select

  end subroutine user_set_outerbcs
  !==========================================================================
  subroutine user_calc_sources
    use ModAdvance,     ONLY: Source_VC, State_VGB
    use ModGeometry,    ONLY: z_BLK
    use ModPhysics,     ONLY: No2Si_V,UnitEnergyDens_,UnitT_
    use ModVarIndexes,  ONLY: Energy_, rhoUz_

    integer :: i, j, k
    real    :: RadiativeCooling = 0.0, EInternalSource =0.0, rhoUzSource =0.0,&
         DampingRhoUz = 0.0, DampingEnergy = 0.0, TeSi

    character (len=*), parameter :: NameSub = 'user_calc_sources'
    !------------------------------------------------------------------------
    RadiativeCooling = 0.

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(UseVerticalDamping) call get_vertical_damping( &
            State_VGB(:,i,j,k,GlobalBLK), &
            z_BLK(i,j,k,GlobalBLK),DampingRhoUz, DampingEnergy)
       if(UseThinRadiation)then
          call user_material_properties(State_VGB(:,i,j,k,GlobalBLK), &
               TeOut=TeSi)
          call get_radiative_cooling(&
               State_VGB(:,i,j,k,GlobalBLK), TeSi, &
               z_BLK(i,j,k,GlobalBLK), RadiativeCooling)
       end if

       EInternalSource = RadiativeCooling
       rhoUzSource     = DampingRhoUz
       srcthin_GB(i,j,k,GlobalBLK) = RadiativeCooling &
            *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitT_)*10.  ! CGS unit
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
    use ModPhysics,     ONLY: Si2No_V, UnitT_
    use ModVarIndexes,  ONLY: rhoUz_, rho_, nVar

    real, intent(in) :: State_V(nVar), Z_V
    real, intent(out):: DampingRhoUz, DampingEnergy
    !-----------------------------------------------------------------------
    !if(UseUniformInitialState)then
    !   DampingRhoUz  = -State_V(rhoUz_)/ &
    !        (TimeVerticalDamping*Si2No_V(UnitT_))
    !   DampingEnergy = -State_V(rhoUz_)**2/State_V(rho_)/&
    !        (TimeVerticalDamping*Si2No_V(UnitT_))
    !else
    if(Z_V > z_photo)then
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
    use ModVarIndexes, ONLY: rho_, nVar
    use ModPhysics,    ONLY: UnitRho_, UnitEnergyDens_, Si2No_V, UnitT_
    use ModLookupTable,ONLY: interpolate_lookup_table
    use ModConst,      ONLY: cProtonMass

    real, intent(in) :: State_V(1:nVar)
    real, intent(in) :: TeSi
    real, intent(in) :: Z_V
    real, intent(out):: RadiativeCooling

    real, parameter :: RadiationCutoff = - 1.0e9, atrl = 1.0e3, Const = 2.5
    real :: Fraction = 0.0, CoolingFunctionCgs = 0.0, &
         MassDensCgs, NumberDensCgs
    real :: Chianti(1:1)
    !------------------------------------------------------------------------

    MassDensCgs     = State_V(rho_)/Si2No_V(UnitRho_)*1.e-3
    NumberDensCgs   = MassDensCgs/(mu*cProtonMass*1.e3)

    ! calculate the thin radiative loss above the z_phto height
    if (z_V > -30.0) then
       ! Smoothing function is 1 if rho<RhoThinCutoff , 0 if not
       Fraction = 0.5 - 0.5*&
            tanh(atrl*(MassDensCgs/RhoThinCutoff - 1.))
       ! Calculate the cooling function culve dependent on temperature
       if (TeSi <= 8.0e+03)&
            CoolingFunctionCgs = (1.0606e-06*TeSi)**11.7
       if ((TeSi <= 1.0e+04).and.(TeSi > 8.0e+03) ) &
            CoolingFunctionCgs = (1.397e-08*TeSi)**6.15
       if(TeSi > 1.0e+04)then
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

    integer, intent(in) :: iStage,iBlock

    integer:: i,j,k
    real   :: MassDensFloor, EnergyFloor, EinternalSi, PressureSi

    character (len=*), parameter :: NameSub = 'user_update_states'
    !--------------------------------------------------------------------
    ! Define the minimum mass density and energy value
    ! Corresponding to ~NumberDensFloor and 200 K plasma
    MassDensFloor  = NumberDensFloor*1e6*cProtonMass*Si2No_V(UnitRho_)
    EnergyFloor    = MassDensFloor*No2Si_V(UnitRho_)/rstari*2000.

    call update_states_MHD(iStage,iBlock)

    do k = 1,nK; do j=1,nJ; do i=1,nI
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
            (EInternalSi - PressureSi*inv_gm1)
    end do;end do; end do

    ! calculate the total energy
    call calc_energy_cell(iBlock)

  end subroutine user_update_states

  !==========================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G,PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)
    use ModAdvance,  ONLY:State_VGB

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

    integer :: i,j,k

    character (len=*), parameter :: Name='user_set_plot_var'
    !-----------------------------------------------------------------
    UsePlotVarBody = .true.
    PlotVarBody    = 0.0
    IsFound        = .true.
    select case(NameVar)
    case('srcthin')
       NameTecVar = 'srcthin'
       NameTecUnit = '[erg/s/cm^3]'
       NameIdlUnit = 'erg/s'
       PlotVar_G   = srcthin_GB(:,:,:,iBlock)
    case('tempe')
       NameTecVar  = 'T'
       NameTecUnit = '[K]'
       NameIdlUnit = 'K'
       do k = 0, nK+1; do j = 0, nJ+1; do i =0, nI+1
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               TeOut = PlotVar_G(i,j,k))
       end do; end do; end do
    case('gamma')
       NameTecVar  = 'gamma'
       NameTecUnit = ''
       NameIdlUnit = ''
       do k = 0, nK+1; do j = 0, nJ+1; do i =0, nI+1
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               GammaOut=PlotVar_G(i,j,k))
       end do; end do; end do
    case('entropy')
       NameTecVar  = 'entropy'
       NameTecUnit = ''
       NameIdlUnit = ''
       do k = 0, nK+1; do j = 0, nJ+1; do i =0, nI+1
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               EntropyOut=PlotVar_G(i,j,k))
       end do; end do; end do
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
       OpacityPlanckOut_W, OpacityRosselandOut_W, PlanckOut_W, &
       EntropyOut)

    use ModLookupTable,ONLY: interpolate_lookup_table
    use ModPhysics,    ONLY: No2Si_V, UnitP_, UnitRho_, inv_gm1
    use ModVarIndexes, ONLY: nVar, Rho_, p_, ExtraEInt_
    use ModAdvance,    ONLY: nWave
    !------------------------------------------------------------------------
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
    real, optional, intent(out) :: EntropyOut

    real :: pSi, EinternalSi, RhoSi, TeSi
    real :: Value_V(1:3), Value_I(1)

    character (len=*), parameter :: NameSub = 'user_material_properties'
    !------------------------------------------------------------------------
    ! Density, transformed to SI
    RhoSi = No2Si_V(UnitRho_)*State_V(Rho_)

    ! Find the thermodynamic variables from an equation-of-state (EOS).
    ! The columns of the EOS table represent in SI units:
    !   10log(temperature), 10log(density),
    !   pressure, internal energy, entropy, specific heat, speed-of-sound gamma
    if(present(TeIn))then
       TeSi = TeIn
       call interpolate_lookup_table(iTableEOS, TeSi, RhoSi, &
            Value_V, DoExtrapolate = .false.)
    elseif(present(EinternalOut))then
       pSi = No2Si_V(UnitP_)*State_V(p_)
       ! Find temperature from density and pressure
       call interpolate_lookup_table(iTableEOS, 1, pSi, RhoSi, &
            Value_V, Arg1Out = TeSi, DoExtrapolate=.false.)
    else
       if(present(EinternalIn)) then 
          EinternalSi = EinternalIn
       else
          EinternalSi = &
               (State_V(p_)*inv_gm1 + State_V(ExtraEint_))*No2Si_V(UnitP_)
       end if
       ! Find temperature from density and internal energy
       call interpolate_lookup_table(iTableEOS, 2, EinternalSi, RhoSi, &
            Value_V, Arg1Out = TeSi, DoExtrapolate=.false.)
    end if

    if(present(TeOut)) TeOut = TeSi
    if(present(PressureOut)) PressureOut = Value_V(1)
    if(present(EInternalOut)) EInternalOut = Value_V(2)
    if(present(EntropyOut)) EntropyOut = Value_V(3)

    ! the specific heat needs to be replaced by an EOS version
    !if(present(CvOut)) CvOut = Value_V(4)
    if(present(CvOut)) CvOut = 1./rstari*inv_gm1*RhoSi !CvSi

    ! speed-of-sound gamma ???
    !if(present(GammaOut)) GammaOut = Value_V(5)
    if(present(GammaOut)) then
       call interpolate_lookup_table(iTableGamma, TeSi, RhoSi, &
            Value_I, DoExtrapolate = .false.)
       GammaOut = Value_I(1)
    end if

  end subroutine user_material_properties

end module ModUser
