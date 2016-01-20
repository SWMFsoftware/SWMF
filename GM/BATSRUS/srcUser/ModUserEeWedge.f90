!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!
! ModUser for flux emergence problems on Sun
! We calculate the thin radiative cooling, vertical damping, 
! and call coronal heating from PARAM.in to model the solar 
! atmosphere, following Abbett (2007). The implementation of
! tabular equation of state takes into account of the 
! ionization energy. 
!=========================================================================

module ModUser
  use ModSize, ONLY: x_, y_, z_, &
       nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, MaxBlock
  use ModUserEmpty ,                                   &
       IMPLEMENTED1 => user_read_inputs,               &
       IMPLEMENTED2 => user_init_session,              &
       IMPLEMENTED3 => user_set_ICs,                   &
       IMPLEMENTED4 => user_initial_perturbation,      &
       IMPLEMENTED5 => user_set_cell_boundary,              &
       IMPLEMENTED6 => user_calc_sources,              &
       IMPLEMENTED7 => user_update_states,             &
       IMPLEMENTED8 => user_set_plot_var,              &
       IMPLEMENTED9 => user_material_properties

  include 'user_module.h'

  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'Flux emergence in a wedge'
  !\
  ! UseVerticalDamping: adds damping to vertical velocity
  ! UseThinRadiation:   adds thin radiative cooling
  ! InitialBr, InitialBt: initial magnetic field, added through
  !                      user_initial_perturbation
  ! RhoThinCutoff: the cutoff density for thin radiation( thin radiation =0 
  !                      if rho> RhoThinCutoff)
  ! NumberDensFloor: Minimum density value to prevent negative coronal density 
  ! TimeVerticalDamping: timescale over which the vertical velocity damps
  ! r_photo: height of photosphere
  ! UseAtmReset: reset atmosphere when it starts to collaps (this can happen
  !              before activating the coronal heating)
  !/
  logical :: UseVerticalDamping = .false.
  logical :: UseThinRadiation = .false.
  logical :: UseCoronalField = .false.
  logical :: UseEnergyPert = .false.
  logical :: UseAtmReset = .false.
  real    :: InitialBr, InitialBt
  real    :: RhoThinCutoff, NumberDensFloor, TimeVerticalDamping
  real    :: r_photo
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
  integer :: iTableEOS = -1, iTableRadCool = -1, iTableInitialState = -1
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
          call read_var('r_photo',r_photo)
          call read_var('UnsignedFluxHeight',UnsignedFluxHeight)
          call read_Var('UseThinRadiation', UseThinRadiation)
          call read_var('UseCoronalField',UseCoronalField)
          call read_var('RhoThinCutoff',RhoThinCutoff)
          call read_Var('UseVerticalDamping',UseVerticalDamping)
          call read_var('TimeVerticalDamping', TimeVerticalDamping)
          call read_var('UseAtmReset', UseAtmReset)
          call read_var('DtUptateFlux',DtUpdateFlux)
          call read_var('UseEnergyPert', UseEnergyPert)
          call read_var('InitialBr',InitialBr)
          call read_var('InitialBphi',InitialBphi)
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

    character (len=*), parameter :: NameSub = 'user_init_session'
    !------------------------------------------------------------------------
    !Gbody = -1 
    !mu    =  MassIon_I(1)/(1 + AverageIonCharge)
    !rstari = mu/(cBoltzmann/cProtonMass)  

    ! initialize the indexes for lookup tables
    iTableInitialState = i_lookup_table('RhoExtraEP(R)')
    iTableEOS          = i_lookup_table('eos(T,rho)')
    iTableRadCool      = i_lookup_table('radcool')

    if(iProc==0) write(*,*) NameSub, &
         'iTableInitialState, EOS , Chianti = ', &
         iTableInitialState, iTableEOS, iTableRadCool

    if(.not.allocated(srcthin_GB)) &
         allocate(srcthin_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))

  end subroutine user_init_session
  !==========================================================================
  subroutine user_set_ICs(iBlock)

    use ModAdvance,     ONLY: State_VGB
    use ModGeometry,    ONLY: r_BLK
    use ModVarIndexes
    use ModLookupTable, ONLY: interpolate_lookup_table
    use ModPhysics,     ONLY: Si2No_V, UnitRho_, UnitEnergyDens_, UnitP_, &
         UnitX_, No2Si_V

    integer, intent(in) :: iBlock

    integer :: i
    real :: InitialState_V(3), InitialRho, InitialExtraE, InitialP
    !------------------------------------------------------------------------

    do i = 1, nI
       ! interpolate the tabular data of reference initial state and get
       ! a relaxed initial state
       call interpolate_lookup_table(iTableInitialState, &
            r_BLK(i,1,1,iBlock)*No2Si_V(UnitX_), InitialState_V, &
            DoExtrapolate = .false.)

       InitialRho    = InitialState_V(1)*Si2No_V(UnitRho_)
       InitialExtraE = InitialState_V(2)*Si2No_V(UnitEnergyDens_)
       InitialP      = InitialState_V(3)*Si2No_V(UnitP_)

       State_VGB(rho_,i,:,:,iBlock) = InitialRho
       State_VGB(rhoUx_:rhoUz_,i,:,:,iBlock) = 0.0
       State_VGB(Bx_:Bz_,i,:,:,iBlock) = 0.
       if(Erad_ > 1) State_VGB(Erad_,i,:,:,iBlock) = 0.
       State_VGB(ExtraEint_,i,:,:,iBlock) = InitialExtraE
       State_VGB(p_,i,:,:,iBlock) = InitialP
    end do

  end subroutine user_set_ICs

  !=========================================================================

  subroutine user_initial_perturbation

    use ModEnergy,   only: calc_energy_cell
    use ModMain,     ONLY: Unused_B, nBlockMax
    use ModGeometry, ONLY: Xyz_DGB, r_BLK
    use ModAdvance,  ONLY: State_VGB
    use ModPhysics,  ONLY: Si2No_V, UnitB_, UnitP_, UnitRho_, UnitX_, &
         UnitEnergyDens_, No2Si_V, InvGammaMinus1
    use ModVarIndexes
    use ModConst,    ONLY: mSun, rSun, cProtonMass, cGravitation, cBoltzmann
    use ModLookupTable, ONLY: interpolate_lookup_table

    integer :: iBlock, i, j, k
    real :: dp_ratio, Prof, rsq, rasq, RandomChange
    real :: Runit_D(3), PhiUnit_D(3)

    ! atmosphere parameters 
    real :: Value_V(3)
    real, parameter ::  mu = 0.594354
    real :: x2c, x3c, qfac, wdth
    real :: r, r_ph, r_ch, r_tr, r_cr, g_ch, g_tr, g_cr 
    real :: rho_ch, p_ch, rho_tr, p_tr, rho_cr, p_cr 
    real :: T_ch, T_tr, T_cr, H_ch, H_tr, H_cr
    real :: p, ExtraEint

    character (len=*), parameter :: NameSub = 'user_initial_perturbation'
    !-----------------------------------------------------------------------   
    rasq = ra_rope*ra_rope
    wdth = 2.0*1.5e5*Si2No_V(UnitX_)

    ! atmosphere parameters
    r_ph = 1.0015
    r_ch = 1.0033
    r_tr = 1.012 
    r_cr = 1.120
    T_ch = 1.5e4
    T_tr = 1.0e5
    T_cr = 1.0e6
    p_ch = 2.00e-03*1.0e-1  !mks units
    rho_ch = p_ch*(mu*cProtonMass)/(cBoltzmann*T_ch) !1.434e-14

    g_ch  = cGravitation*mSun/(r_ch*rSun)**2
    g_tr  = cGravitation*mSun/(r_tr*rSun)**2
    g_cr  = cGravitation*mSun/(r_cr*rSun)**2

    H_ch  = cBoltzmann*T_ch/(g_ch*mu*cProtonMass)*Si2No_V(UnitX_)
    H_tr  = cBoltzmann*T_tr/(g_tr*mu*cProtonMass)*Si2No_V(UnitX_)
    H_cr  = cBoltzmann*T_cr/(g_cr*mu*cProtonMass)*Si2No_V(UnitX_)

    p_tr   =  p_ch*exp(-(r_tr - r_ch)/H_ch) 
    rho_tr =  p_tr*(mu*cProtonMass)/(cBoltzmann*T_tr)

    p_cr   =  p_tr*exp(-(r_cr - r_tr)/H_tr)
    rho_cr =  p_cr*(mu*cProtonMass)/(cBoltzmann*1.0e8)

    rho_ch = rho_ch*Si2No_V(UnitRho_)
    rho_tr = rho_tr*Si2No_V(UnitRho_)
    rho_cr = rho_cr*Si2No_V(UnitRho_)
    p_ch   =   p_ch*Si2No_V(UnitP_)
    p_tr   =   p_tr*Si2No_V(UnitP_)
    p_cr   =   p_cr*Si2No_V(UnitP_)

    Phiunit_D(3) =  0.0

    do iBlock = 1, nBlockMax
       if(Unused_B(iBlock)) CYCLE

       do i = 1, nI
          ! Add initial magnetic field 
          if(UseCoronalField)then
             do k = 1, nK; do j = 1, nJ
                Runit_D = Xyz_DGB(:,i,j,k,iBlock)/r_BLK(i,j,k,iBlock)
                Phiunit_D(1) =  Xyz_DGB(y_,i,j,k,iBlock)/r_BLK(i,j,k,iBlock)
                Phiunit_D(2) = -Xyz_DGB(x_,i,j,k,iBlock)/r_BLK(i,j,k,iBlock)

                if(r_BLK(i,j,k,iBlock) < r_photo)then
                   Prof = 0.5*(1 - tanh(r_BLK(i,j,k,iBlock) - r_photo))
                else
                   Prof = 0.0
                end if
                State_VGB(Bx_:Bz_,i,j,k,iBlock) = &
                     State_VGB(Bx_:Bz_,i,j,k,iBlock) &
                     + (PhiUnit_D*InitialBphi*Prof + Runit_D*InitialBr) &
                     *1e-4*Si2No_V(UnitB_)
             end do; end do
          end if

          ! Add random perturbation to energy and pressure values of 
          ! cells below the photosphere height                  
          if(UseEnergyPert .and. (r_BLK(i,1,1,iBlock) < r_photo))then
             do k = 1, nK; do j = 1, nJ
                call random_number(RandomChange)
                RandomChange = (RandomChange - 0.5)*2
                if(iTableEOS > 0) State_VGB(ExtraEint_,i,j,k,iBlock) = &
                     State_VGB(ExtraEint_,i,j,k,iBlock) &
                     *(1.0 + 1e-3*RandomChange)
                State_VGB(p_,i,j,k,iBlock) = State_VGB(p_,i,j,k,iBlock) &
                     *(1.0 + 1e-3*RandomChange)
             end do; end do
          end if

          ! If UseRope, Add flux rope in
          ! Set negative pressure to 1.e-10
          if(UseRope)then
             ! Flux rope does not yet have the toroidal curvature
             do k = 1, nK ; do j = 1, nJ
                rsq = (Xyz_DGB(y_,i,j,k,iBlock) - x2c_rope)**2 + &
                     (Xyz_DGB(z_,i,j,k,iBlock) - x3c_rope)**2
                if(rsq < rasq**1.5e1)then
                   prof = b0_rope*exp(-rsq/rasq)
                else
                   prof = 0.
                end if
                dp_ratio =  0.5*prof*prof*(-1.+0.5*qfac_rope*qfac_rope* &
                     (1. - 2.*rsq/rasq))/State_VGB(p_,i,j,k,iBlock)
                !B: field needs to be rotated
                State_VGB(Bx_,i,j,k,iBlock) = &
                     State_VGB(Bx_,i,j,k,iBlock) + prof
                State_VGB(By_,i,j,k,iBlock) = &
                     State_VGB(By_,i,j,k,iBlock) &
                     - prof*qfac_rope*(Xyz_DGB(z_,i,j,k,iBlock) - x3c_rope) &
                     /ra_rope
                State_VGB(Bz_,i,j,k,iBlock) = &
                     State_VGB(Bz_,i,j,k,iBlock) + &
                     prof*qfac_rope*(Xyz_DGB(y_,i,j,k,iBlock) - x2c_rope) &
                     /ra_rope
                State_VGB(P_,i,j,k,iBlock) = &
                     State_VGB(p_,i,j,k,iBlock)*(1.+dp_ratio)
                State_VGB(rho_,i,j,k,iBlock) = State_VGB(rho_,i,j,k,iBlock) &
                     *(1. +  exp(-(Xyz_DGB(x_,i,j,k,iBlock)/lamb_rope)**2) &
                     *dp_ratio*buoyancy_rope)
                if(State_VGB(p_,i,j,k,iBlock) <= 0.) &
                     State_VGB(p_,i,j,k,iBlock) = 1.e-10
             end do; end do
          end if

          ! If UseAtmReset, reset upper atmosphere 
          if(UseAtmReset)then  
             r = r_BLK(i,1,1,iBlock)
             do k = 1, nK; do j = 1, nJ
                !CH
                if( (r > r_ch) .and. (r < r_tr) )then
                   State_VGB(rho_,i,j,k,iBlock) = rho_ch*exp(-(r-r_ch)/H_ch)
                   State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
                   if(iTableEOS > 0)then
                      call interpolate_lookup_table(iTableEOS, T_ch, &
                           State_VGB(rho_,i,j,k,iBlock)*No2Si_V(UnitRho_), &
                           Value_V, DoExtrapolate = .false.)
                      p         = Value_V(1)*Si2No_V(UnitP_)
                      ExtraEint = (Value_V(2)-Value_V(1)*InvGammaMinus1) &
                           *Si2No_V(UnitEnergyDens_)
                   else
                      p         = p_ch*exp(-(r-r_ch)/H_ch)
                      ExtraEint = 0.
                   end if
                   State_VGB(ExtraEint_,i,j,k,iBlock) = ExtraEint
                   State_VGB(p_,i,j,k,iBlock) = p
                end if
                !TR
                if(r > r_tr) then
                   !if( (r > r_tr) .and. (r < r_cr) )then
                   State_VGB(rho_,i,j,k,iBlock) = rho_tr*exp(-(r-r_tr)/H_tr)
                   State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
                   if(iTableEOS > 0)then
                      call interpolate_lookup_table(iTableEOS, T_tr, &
                           State_VGB(rho_,i,j,k,iBlock)*No2Si_V(UnitRho_), &
                           Value_V, DoExtrapolate = .false.)
                      p         = Value_V(1)*Si2No_V(UnitP_)
                      ExtraEint = (Value_V(2)-Value_V(1)*InvGammaMinus1) &
                           *Si2No_V(UnitEnergyDens_)
                   else
                      p         = p_tr*exp(-(r-r_tr)/H_tr)
                      ExtraEint = 0.
                   end if
                   State_VGB(ExtraEint_,i,j,k,iBlock) = ExtraEint
                   State_VGB(p_,i,j,k,iBlock) = p
                end if
                !CR
                if((r > r_cr) .and. (r < r_tr))then
                   State_VGB(rho_,i,j,k,iBlock) = rho_cr*exp(-(r-r_cr)/H_cr)
                   State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
                   if(iTableEOS > 0)then
                      call interpolate_lookup_table(iTableEOS, T_cr, &
                           State_VGB(rho_,i,j,k,iBlock)*No2Si_V(UnitRho_), &
                           Value_V, DoExtrapolate = .false.)
                      p         = Value_V(1)*Si2No_V(UnitP_)
                      ExtraEint = (Value_V(2)-Value_V(1)*InvGammaMinus1) &
                           *Si2No_V(UnitEnergyDens_)
                   else
                      p         = p_cr*exp(-(r-r_cr)/H_cr)
                      ExtraEint = 0.
                   end if
                   State_VGB(ExtraEint_,i,j,k,iBlock) = ExtraEint
                   State_VGB(p_,i,j,k,iBlock) = p
                end if
             end do; end do
          end if
       end do

       call calc_energy_cell(iBlock)
    end do

  end subroutine user_initial_perturbation

  !==========================================================================
  subroutine user_set_cell_boundary(iBlock,iSide, TypeBc, IsFound)

    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUz_, Bx_, Bz_, p_, ExtraEInt_
    use ModAdvance,    ONLY: State_VGB
    use ModGeometry,   ONLY: Xyz_DGB, r_BLK
    use ModLookupTable, ONLY: interpolate_lookup_table
    use ModPhysics,    ONLY: Si2No_V, No2Si_V, UnitRho_, UnitP_, &
         UnitEnergyDens_, UnitX_

    integer,          intent(in)  :: iBlock, iSide
    character(len=*), intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    integer :: i, j, k
    real :: Rho, p, ExtraEint, Runit_D(3), RhoUr, InitialState_V(3)

    character (len=*), parameter :: NameSub = 'user_set_cell_boundary'
    !------------------------------------------------------------------------

    select case(iSide)
    case(1)
       select case(TypeBc)
       case('fixvalue')

          do i = MinI, 0
             call interpolate_lookup_table(iTableInitialState, &
                  r_BLK(i,1,1,iBlock)*No2Si_V(UnitX_), InitialState_V, &
                  DoExtrapolate = .false.)

             Rho       = InitialState_V(1)*Si2No_V(UnitRho_)
             ExtraEint = InitialState_V(2)*Si2No_V(UnitEnergyDens_)
             p         = InitialState_V(3)*Si2No_V(UnitP_)

             do k = MinK, MaxK; do j = MinJ, MaxJ

                State_VGB(Rho_,i,j,k,iBlock) = Rho
                State_VGB(p_,i,j,k,iBlock) = p
                State_VGB(ExtraEint_,i,j,k,iBlock) = ExtraEint

                Runit_D = Xyz_DGB(:,1,j,k,iBlock)/r_BLK(1,j,k,iBlock)

                ! Reflect total vector to reduce surface strong flows at
                ! bottom boundary
                State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = &
                     -State_VGB(RhoUx_:RhoUz_,1-i,j,k,iBlock) 
                ! Float Br
                State_VGB(Bx_:Bz_,i,j,k,iBlock) = &
                     sum(State_VGB(Bx_:Bz_,1,j,k,iBlock)*Runit_D)*Runit_D

             end do; end do
          end do
          IsFound = .true.
       end select
    case(2)
       select case(TypeBc)
       case('no_inflow')
          ! closed upper boundary condition, not allowing upward motions
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = nI+1, MaxI
             ! Float
             State_VGB(:,i,j,k,iBlock) = State_VGB(:,nI,j,k,iBlock)

             Runit_D = Xyz_DGB(:,nI,j,k,iBlock)/r_BLK(nI,j,k,iBlock)

             RhoUr = sum(Runit_D*State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock))

             State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = &
                  State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) &
                  + (-RhoUr + abs(RhoUr))*Runit_D

          end do; end do; end do
          IsFound = .true.
       end select
    end select

  end subroutine user_set_cell_boundary
  !==========================================================================
  subroutine user_calc_sources(iBlock)

    use ModAdvance,     ONLY: Source_VC, State_VGB
    use ModPhysics,     ONLY: No2Si_V,UnitEnergyDens_,UnitT_
    use ModVarIndexes,  ONLY: Energy_, RhoUx_, RhoUz_
    use ModGeometry,    ONLY: Xyz_DGB, r_BLK

    integer, intent(in) :: iBlock

    integer :: i, j, k
    real :: RadiativeCooling = 0.0, EInternalSource = 0.0, rhoUrSource = 0.0, &
         DampingRhoUr = 0.0, DampingEnergy = 0.0, TeSi
    real :: Runit_D(3)

    character (len=*), parameter :: NameSub = 'user_calc_sources'
    !------------------------------------------------------------------------
    RadiativeCooling = 0.

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       Runit_D = Xyz_DGB(:,i,j,k,iBlock)/r_BLK(i,j,k,iBlock)

       if(UseVerticalDamping) call get_vertical_damping( &
            State_VGB(:,i,j,k,iBlock), &
            r_BLK(i,j,k,iBlock), Runit_D, DampingRhoUr, DampingEnergy)

       if(UseThinRadiation)then
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               TeOut=TeSi)
          call get_radiative_cooling(&
               State_VGB(:,i,j,k,iBlock), TeSi, &
               r_BLK(i,j,k,iBlock), RadiativeCooling)
       end if

       EInternalSource = RadiativeCooling
       rhoUrSource     = DampingRhoUr
       srcthin_GB(i,j,k,iBlock) = RadiativeCooling &
            *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitT_)*10.  ! CGS unit
       Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + &
            EinternalSource + DampingEnergy

       Source_VC(RhoUx_:RhoUz_,i,j,k) = Source_VC(RhoUx_:RhoUz_,i,j,k) &
            + rhoUrSource*Runit_D
    end do; end do; end do

  end subroutine user_calc_sources

  !==========================================================================
  !----------------------------------------------!
  !  subroutines containing energy source terms  !            
  !----------------------------------------------!
  subroutine get_vertical_damping(State_V, r, Runit_D, DampingRhoUr, &
       DampingEnergy)
    use ModPhysics,     ONLY: Si2No_V, UnitT_
    use ModVarIndexes,  ONLY: RhoUx_, rhoUz_, rho_, nVar

    real, intent(in) :: State_V(nVar), r, Runit_D(3)
    real, intent(out):: DampingRhoUr, DampingEnergy

    real :: RhoUr
    !-----------------------------------------------------------------------
    if(r > r_photo)then
       RhoUr = sum(State_V(RhoUx_:RhoUz_)*Runit_D)
       DampingRhoUr  = -RhoUr/ &
            (TimeVerticalDamping*Si2No_V(UnitT_))
       DampingEnergy = -RhoUr**2/State_V(rho_)/&
            (TimeVerticalDamping*Si2No_V(UnitT_))
    else
       DampingRhoUr  = 0.
       DampingEnergy = 0.
    end if
  end subroutine get_vertical_damping
  !=========================================================================
  subroutine get_radiative_cooling(State_V, TeSi, r, RadiativeCooling)

    use ModVarIndexes, ONLY: rho_, nVar
    use ModPhysics,    ONLY: UnitRho_, UnitEnergyDens_, Si2No_V, UnitT_
    use ModLookupTable,ONLY: interpolate_lookup_table
    use ModConst,      ONLY: cProtonMass

    real, intent(in) :: State_V(1:nVar)
    real, intent(in) :: TeSi
    real, intent(in) :: r
    real, intent(out):: RadiativeCooling

    real, parameter :: RadiationCutoff = - 3.0e9, atrl = 1.0e3
    real, parameter :: mf=1.0

    real :: Fraction = 0.0, CoolingFunctionCgs = 0.0, &
         MassDensCgs, NumberDensCgs
    real :: CoolingTableOut_I(1)
    real, parameter :: RadNorm = 1e22
    !------------------------------------------------------------------------

    MassDensCgs     = State_V(rho_)/Si2No_V(UnitRho_)*1.e-3
    NumberDensCgs   = MassDensCgs/(mu*cProtonMass*1.e3)

    ! calculate the thin radiative loss above the r_photo height
    !   if (r > -30.0) then  
    !    coordinate transform 1.0 + tmp(0)/6.96E+8
    if (r > 1.0 -30.0*1.50E5/6.96e8) then
       ! Smoothing function is 1 if rho<RhoThinCutoff , 0 if not
       Fraction = mf*(0.5 - 0.5*tanh(atrl*(MassDensCgs/RhoThinCutoff - 1.)))
       ! Calculate the cooling function culve dependent on temperature
       if(TeSi <= 8e3)then
          CoolingFunctionCgs = (1.0606e-6*TeSi)**11.7
       elseif(TeSi <= 1e4)then
          CoolingFunctionCgs = (1.397e-8*TeSi)**6.15
       else
          ! Table variable should be normalized to radloss_cgs * 1e22
          ! since we don't want to deal with such tiny numbers
          call interpolate_lookup_table(iTableRadCool, TeSi, &
               CoolingTableOut_I, DoExtrapolate = .false.)
          CoolingFunctionCgs = CoolingTableOut_I(1)/RadNorm
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
         UnitEnergyDens_, cProtonMass, InvGammaMinus1
    use ModEnergy,     ONLY: calc_energy_cell

    integer, intent(in) :: iStage,iBlock

    integer:: i, j, k
    real   :: MassDensFloor, EnergyFloor, EinternalSi, PressureSi

    character (len=*), parameter :: NameSub = 'user_update_states'
    !--------------------------------------------------------------------
    ! Define the minimum mass density and energy value
    ! Corresponding to ~NumberDensFloor and 200 K plasma
    MassDensFloor  = NumberDensFloor*1e6*cProtonMass*Si2No_V(UnitRho_)
    EnergyFloor    = MassDensFloor*No2Si_V(UnitRho_)/rstari*2000.

    call update_states_MHD(iStage,iBlock)

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       ! check if the density or pressure is below the minimum value
       if(State_VGB(rho_,i,j,k,iBlock) < MassDensFloor) &
            State_VGB(rho_,i,j,k,iBlock) = MassDensFloor
       ! Total internal energy, ExtraEInt + P/(\gamma -1),
       EInternalSi = (InvGammaMinus1*State_VGB(P_,i,j,k,iBlock) + &
            State_VGB(ExtraEInt_,i,j,k,iBlock))*No2Si_V(UnitEnergyDens_)
       if(EInternalSi < EnergyFloor) EInternalSi = EnergyFloor
       ! get pressure and extra internal energy from the EOS table
       call user_material_properties(State_VGB(:,i,j,k,iBlock), &
            EInternalIn=EInternalSi, PressureOut=PressureSi)
       State_VGB(P_,i,j,k,iBlock) = PressureSi*Si2No_V(UnitP_)
       State_VGB(ExtraEInt_,i,j,k,iBlock) = Si2No_V(UnitEnergyDens_)*&
            (EInternalSi - PressureSi*InvGammaMinus1)
    end do; end do; end do

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
    real,             intent(out):: PlotVar_G(MinI:MaxI, MinJ:MaxJ, MinK:MaxK)
    real,             intent(out):: PlotVarBody
    logical,          intent(out):: UsePlotVarBody
    character(len=*), intent(out):: NameTecVar
    character(len=*), intent(out):: NameTecUnit
    character(len=*), intent(out):: NameIdlUnit
    logical,          intent(out):: IsFound

    integer :: i, j, k

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
       OpacityPlanckOut_W, OpacityEmissionOut_W, OpacityRosselandOut_W, &
       PlanckOut_W, EntropyOut)

    use ModLookupTable,ONLY: interpolate_lookup_table
    use ModPhysics,    ONLY: No2Si_V, UnitP_, UnitRho_, InvGammaMinus1
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
    real, optional, intent(out) :: OpacityEmissionOut_W(nWave)    ! [1/m] 
    real, optional, intent(out) :: OpacityRosselandOut_W(nWave)   ! [1/m] 
    real, optional, intent(out) :: PlanckOut_W(nWave)      ! [J/m^3] 
    real, optional, intent(out) :: EntropyOut

    real :: pSi, EinternalSi, RhoSi, TeSi
    real :: Value_V(1:5)

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
               (State_V(p_)*InvGammaMinus1+State_V(ExtraEint_))*No2Si_V(UnitP_)
       end if
       ! Find temperature from density and internal energy
       call interpolate_lookup_table(iTableEOS, 2, EinternalSi, RhoSi, &
            Value_V, Arg1Out = TeSi, DoExtrapolate=.false.)
    end if

    if(present(TeOut)) TeOut = TeSi
    if(present(PressureOut)) PressureOut = Value_V(1)
    if(present(EInternalOut)) EInternalOut = Value_V(2)
    if(present(EntropyOut)) EntropyOut = Value_V(3)
    if(present(CvOut)) CvOut = Value_V(4)
    if(present(GammaOut)) GammaOut = Value_V(5)

  end subroutine user_material_properties

end module ModUser
