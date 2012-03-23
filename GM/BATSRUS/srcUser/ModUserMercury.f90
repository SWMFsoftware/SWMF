!^CFG COPYRIGHT UM
!========================================================================
module ModUser
  ! This is the default user module which contains empty methods defined
  ! in ModUserEmpty.f90

  use ModUserEmpty,                          &
       IMPLEMENTED1 => user_init_session,    &
       IMPLEMENTED2 => user_set_ics,         &
       IMPLEMENTED3 => user_update_states,   &
       IMPLEMENTED4 => user_set_resistivity, &
       IMPLEMENTED5 => user_read_inputs

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = 'Mercury, Lars Daldorff'

  real :: PlanetDensity=-1., PlanetPressure=-1., PlanetMaxResistivity=-1., &
       PlanetRadius=-1., PlanetResistiveR_I(2)= (/-1.,-1./)

  real :: PlanetDensitySi=-1., PlanetPressureSi=-1., PlanetMaxResistivitySi=-1., &
       PlanetRadiusSi=-1., PlanetResistiveRSi_I(2)=(/-1.,-1./)

  real :: ResistivityRate

contains

  !===========================================================================

  subroutine user_init_session

    use CON_planet,     ONLY: RadiusPlanet, MassPlanet
    use ModNumConst,    ONLY: cPi
    use ModPhysics,     ONLY: Si2No_V,No2Si_V,UnitRho_, &
                              UnitP_, UnitX_
    use ModIO,          ONLY: write_prefix, write_myname, iUnitOut
    use ModProcMH,      ONLY: iProc
    use ModResistivity, ONLY: Si2NoEta

    CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(A20,E10.3)"
    !-------------------------------------------------------------------

    !!! Introduce PlanetDensitySi etc., read those and convert 
    !!! from that. Initialize these to -1. PlanetRadius should be
    !!! in dimensional units.

    if (index(TypeGeometry,'spherical')>0) &
         call stop_mpi('ERROR: Correct PARAM.in, need spherical grid.')


    if(PlanetDensitySi < 0.0) &
         PlanetDensitySI  = 3.0*MassPlanet/(4.0*cPi*RadiusPlanet**3)

    if(PlanetRadius < 0.0) &
         PlanetRadius = RadiusPlanet*Si2No_V(UnitX_)

    if(PlanetPressureSi < 0.0) &
         PlanetPressureSi = 1.0e-8*No2Si_V(UnitP_)

    if(PlanetMaxResistivitySi < 0.0) &
         PlanetMaxResistivitySi = 1.0e13

    if(sum(PlanetResistiveR_I) < 0.0 ) &
         PlanetResistiveR_I(:) = (/1.0, 0.0/)
    
 
    PlanetDensity           = PlanetDensitySi*Si2No_V(UnitRho_)
    PlanetPressure          = PlanetPressureSi*Si2No_V(UnitP_)
    PlanetMaxResistivity    = PlanetMaxResistivitySi*Si2NoEta

    PlanetRadiusSi          = PlanetRadius*No2Si_V(UnitX_)
    PlanetResistiveRSi_I(:) = PlanetResistiveR_I(:)*No2Si_V(UnitX_)

    ResistivityRate = &
         PlanetMaxResistivity/(PlanetResistiveR_I(1)-PlanetResistiveR_I(2))

    if(iProc==0) then
       call write_myname
       write(*,*) ''
       write(*,FMT1) '  Planet density  = ',PlanetDensitySi
       write(*,FMT1) '  Planet pressure = ',PlanetPressureSi
       write(*,FMT1) '  Planet radius   = ',PlanetRadiusSi
       write(*,*) ''
       write(*,"(A45,E10.3)") '  Planet resistivity goes linear from radius',&
            PlanetResistiveRSi_I(1)
       write(*,"(A6,E10.3,A27,E10.3)") '  to ', &
            PlanetResistiveRSi_I(2), &
            ' where it has maximum value :', &
            PlanetMaxResistivitySi 
       write(*,*) ''
       write(*,*) ''
    end if

  end subroutine user_init_session
  
  !===========================================================================

  subroutine user_read_inputs

    use ModMain
    use ModProcMH,      ONLY: iProc
    use ModReadParam
    character(len=100) :: NameCommand
    !-------------------------------------------------------------------
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#RESISTIVEPLANET")
          call read_var('PlanetDensitySi'       , PlanetDensitySi)
          call read_var('PlanetPressureSi'      , PlanetPressureSi)
          call read_var('PlanetMaxResistivitySi', PlanetMaxResistivitySi) 
          call read_var('PlanetRadiusDim'        , PlanetRadius)
          call read_var('OuterResistiveRadiusDim', PlanetResistiveR_I(1))
          call read_var('InerResistiveRadiusDim' , PlanetResistiveR_I(2))
          if(sum(PlanetResistiveR_I) >= 0.0)then
             if(PlanetResistiveR_I(1) <= PlanetResistiveR_I(2)) then
                write(*,*) 'ERROR: OuterResistivRadius has to be '//&
                     'larger than InnerResistiveRadius'
                call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
             end if
          end if
       case('#USERINPUTEND')
          EXIT
       case default
          if(iProc==0) then
             write(*,*) &
                  'ERROR: Invalid user defined #COMMAND in user_read_inputs. '
             write(*,*) '--Check user_read_inputs for errors'
             write(*,*) '--Check to make sure a #USERINPUTEND command was used'
             write(*,*) ' *Unrecognized command was: '//NameCommand
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if
       end select
    end do

  end subroutine user_read_inputs

  !===========================================================================
  
  subroutine user_set_ics

    use ModMain,ONLY: globalBLK
    use ModAdvance,    ONLY: State_VGB
    use ModMain,       ONLY: UnUsedBlk, nBlockMax
    use ModGeometry,   ONLY: R_BLK, Rmin_BLK
    use ModSize
    use ModVarIndexes, ONLY: Rho_, rhoUx_, rhoU_, P_,Bx_,Bz_
    use ModNumConst,   ONLY:cPi

    use CON_planet,    ONLY: RadiusPlanet, MassPlanet
    use ModNumConst,   ONLY:cPi
    use ModPhysics,    ONLY: Si2No_V,UnitRho_ 

    use ModMultiFluid, ONLY: select_fluid, iFluid, nFluid, iP,iEnergy, &
         iRho, iRhoUx, iRhoUy, iRhoUz

    integer :: iBlock,i,j,k
    character (len=*), parameter :: NameSub = 'user_set_ics'
    !-------------------------------------------------------------------
    !The routine is called for each block, the number of block should
    !be passed to the routine using globalBLK

    iBlock = globalBLK

    if(Rmin_BLK(iBlock) > PlanetRadius) return

    do iFluid = 1, nFluid
       call select_fluid
       do k=1,nK; do j=1,nJ; do i=1,nI
          if(R_BLK(i+2,j,k,iBlock) > PlanetRadius) CYCLE
          State_VGB(iRho,i,j,k,iBlock) = PlanetDensity
          State_VGB(iP,i,j,k,iBlock) = PlanetPressure
          State_VGB(irhoUx:irhoUz,i,j,k,iBlock) = 0.0
       end do; end do; end do
    end do

    do k=1,nK; do j=1,nJ; do i=1,nI
       State_VGB(Bx_:Bz_,i,j,k,iBlock) = 0.0
    end do; end do; end do

  end subroutine user_set_ics

  !===========================================================================

  subroutine user_update_states(iStage,iBlock)

    use CON_planet,    ONLY: RadiusPlanet, MassPlanet
    use ModAdvance,    ONLY: State_VGB
    use ModMain,       ONLY: UnUsedBlk, nBlockMax
    use ModGeometry,   ONLY: X_BLK, Y_BLK, Z_BLK, R_BLK, Rmin_BLK
    use ModSize
    use ModPhysics,    ONLY: Si2No_V,UnitRho_ 
    use ModVarIndexes, ONLY: nVar, Rho_, RhoU_, RhoUx_, RhoUz_, Bx_, Bz_, p_
    use ModEnergy,     ONLY: calc_energy_cell

    integer,intent(in) :: iStage, iBlock

    real :: r_D(3), dRhoUr_D(3), RhoUr
    integer :: i, j, k, iVar
    character (len=*), parameter :: NameSub = 'user_set_ics'
    !----------------------------------------------------------------------

    call update_states_MHD(iStage,iBlock)

    if(Rmin_BLK(iBlock) <= PlanetRadius) then
       do k=1,nK; do j=1,nJ; do i=1,nI
          if(R_BLK(i+2,j,k,iBlock) >= PlanetRadius ) CYCLE
          State_VGB(Rho_,i,j,k,iBlock) = PlanetDensity
          State_VGB(P_,i,j,k,iBlock) = PlanetPressure
          State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock) = 0.0
       end do; end do; end do



       do k=MinK,MaxK; do j=MinJ,MaxJ; do i=1,nI
          if(R_BLK(i-1,j,k,iBlock) > PlanetRadius ) CYCLE
          if(R_BLK(i,j,k,iBlock) <= PlanetRadius ) CYCLE

          r_D  = (/ x_BLK(i,j,k,iBlock), y_BLK(i,j,k,iBlock), &
               z_BLK(i,j,k,iBlock) /) / r_BLK(i,j,k,iBlock)

          RhoUr = dot_product(State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock),r_D)
          if(RhoUr > 0.0) then
             ! If flow is out of the planet, remove the radial componet 
             ! of the momentum so that the flow is tangential
             dRhoUr_D = -r_D*RhoUr
          else
             ! If flow is into the planet do nothing so the flow is absorbed
             dRhoUr_D = (/0.0,0.0,0.0/)
          end if

          ! Two cell boundary layer inside the planet
          !! -1
          State_VGB(Rho_,i-1,j,k,iBlock) = State_VGB(Rho_,i,j,k,iBlock)
          State_VGB(RhoUx_:RhoUz_,i-1,j,k,iBlock) = &
               State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) + dRhoUr_D
          State_VGB(Bz_+1:nVar,i-1,j,k,iBlock) = &
               State_VGB(Bz_+1:nVar,i,j,k,iBlock)
          !! -2
          State_VGB(Rho_,i-2,j,k,iBlock) = State_VGB(Rho_,i,j,k,iBlock)
          State_VGB(RhoUx_:RhoUz_,i-2,j,k,iBlock) = &
               State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) + dRhoUr_D
          State_VGB(Bz_+1:nVar,i-2,j,k,iBlock) = &
               State_VGB(Bz_+1:nVar,i,j,k,iBlock)

       end do; end do; end do

    end if

    call calc_energy_cell(iBlock)

  end subroutine user_update_states

  !===========================================================================

  subroutine user_set_resistivity(iBlock, Eta_G)
    use ModEnergy
    use BATL_size, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModGeometry,   ONLY: R_BLK, Rmin_BLK
    use ModResistivity, ONLY: Eta0
    integer, intent(in) :: iBlock
    real, intent(out) :: Eta_G(MinI:MaxI, MinJ:MaxJ, MinK:MaxK) 

    integer ::i,j,k

    Eta_G = Eta0

    if(Rmin_BLK(iBlock) > PlanetRadius) return

    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
       if(R_BLK(i,j,k,iBlock) > PlanetResistiveR_I(1) ) CYCLE
       if(R_BLK(i,j,k,iBlock) < PlanetResistiveR_I(2) ) CYCLE
       ! to avoid eta jumps adding Eta_G
       Eta_G(i,j,k) = Eta_G(i,j,k) + (PlanetResistiveR_I(1) - R_BLK(i,j,k,iBlock))* &
            ResistivityRate
    end do; end do; end do

  end subroutine user_set_resistivity

end module ModUser
