!^CFG COPYRIGHT UM
!========================================================================
module ModUser
  ! This is the default user module which contains empty methods defined
  ! in ModUserEmpty.f90

  use ModUserEmpty,                       &
       IMPLEMENTED1 => user_init_session, &
       IMPLEMENTED2 => user_set_ics,      &
       IMPLEMENTED3 => user_update_states,&
       IMPLEMENTED4 => user_set_resistivity,&
       IMPLEMENTED5 => user_read_inputs,&
       IMPLEMENTED6 => user_set_plot_var
  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = 'Moon Module, Rubin 2011'

  real :: SW_Pe, SW_Pi
  real :: planetDensity, planetPressure, planetMaxResistivity,&
       planetR, planetResistivR_I(2), ResistivetyRate

contains

  subroutine user_init_session

    use CON_planet,    ONLY: RadiusPlanet, MassPlanet
    use ModNumConst,   ONLY:cPi
    use ModPhysics,    ONLY: Si2No_V,No2Si_V,UnitRho_, UnitP_, UnitX_
    use ModIO, ONLY: write_prefix, write_myname, iUnitOut
    use ModProcMH,     ONLY:iProc

    CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(A20,E10.3)"
    !-------------------------------------------------------------------

    planetDensity  = planetDensity*Si2No_V(UnitRho_)
    planetPressure = planetPressure*Si2No_V(UnitP_)
    planetR        = planetR*Si2No_V(UnitX_)
    planetResistivR_I = planetResistivR_I*Si2No_V(UnitX_)

    if(planetDensity < 0.0) then
       print *, "WHY here!!! ", planetDensity
       planetDensity  = SI2No_V(UnitRho_)*3.0*MassPlanet/(4.0*cPi*RadiusPlanet**3)
    end if

    if(planetR < 0.0) &
         planetR = SI2No_V(UnitX_)*RadiusPlanet

    if(planetPressure < 0.0) &
         planetPressure = 1.0e-8

    if(planetMaxResistivity < 0.0) &
         planetMaxResistivity = 1.0e13

    if(sum(planetResistivR_I) < 0.0 ) &
         planetResistivR_I(:) = (/1.0, 0.0/)

    ResistivetyRate = &
         planetMaxResistivity/(planetResistivR_I(1)-planetResistivR_I(2))

    if(iProc==0) then
       call write_myname
       write(*,*) ''
       write(*,FMT1) '  Planet density  = ',planetDensity*No2Si_V(UnitRho_)
       write(*,FMT1) '  Planet pressure = ',planetPressure*No2Si_V(UnitP_)
       write(*,FMT1) '  Planet radius   = ',planetR*No2Si_V(UnitX_)
       write(*,*) ''
       write(*,"(A45,E10.3)") '  Planet resistivety goes linear from radius',&
            planetResistivR_I(1)*No2Si_V(UnitX_) 
       write(*,"(A6,E10.3,A27,E10.3)") '  to ',planetResistivR_I(2)*No2Si_V(UnitX_), &
            ' where it have max value :', &
            planetMaxResistivity 
       write(*,*) ''
       write(*,*) ''
    end if


  end subroutine user_init_session


  subroutine user_read_inputs

    use ModMain
    use ModProcMH,     ONLY:iProc
    use ModReadParam
    use ModIO, ONLY: write_prefix, write_myname, iUnitOut

    character (len=100) :: NameCommand
    !-------------------------------------------------------------------

    if(iProc==0.and.lVerbose > 0)then
       call write_prefix; write(iUnitOut,*)
    endif

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#RESISTIVPLANET")
          call read_var('Density' , planetDensity)
          call read_var('Pressure' , planetPressure)
          call read_var('MaxResistivity', planetMaxResistivity) 
          call read_var('PlanetRadius' , planetR)
          call read_var('OuterResistivRadius' , planetResistivR_I(1))
          call read_var('InerResistivRadius' , planetResistivR_I(2))
          if(sum(planetResistivR_I) >= 0.0)then 
             if(planetResistivR_I(1) <= planetResistivR_I(2)) then
                write(*,*) ''
                write(*,*) 'Error'
                write(*,*) 'OuterResistivRadius has to be lager then InerResistivRadius'
                call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
             end if
          end if
       case('#USERINPUTEND')
          if(iProc==0.and.lVerbose > 0)then
             call write_prefix;
             write(iUnitOut,*)'User read_input Comet ends'
          endif
          EXIT
       case default
          if(iProc==0) then
             call write_myname; write(*,*) &
                  'ERROR: Invalid user defined #COMMAND in user_read_inputs. '
             write(*,*) '--Check user_read_inputs for errors'
             write(*,*) '--Check to make sure a #USERINPUTEND command was used'
             write(*,*) ' *Unrecognized command was: '//NameCommand
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if
       end select
    end do

  end subroutine user_read_inputs

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

    if(Rmin_BLK(iBlock) > planetR) return

    do iFluid = 1, nFluid
       call select_fluid
       do k=1,nK; do j=1,nJ; do i=1,nI
          if(R_BLK(i+2,j,k,iBlock) > planetR) CYCLE
          State_VGB(iRho,i,j,k,iBlock) = planetDensity
          State_VGB(iP,i,j,k,iBlock) = planetPressure
          State_VGB(irhoUx:irhoUz,i,j,k,iBlock) = 0.0
       end do; end do; end do
    end do

    do k=1,nK; do j=1,nJ; do i=1,nI
       State_VGB(Bx_:Bz_,i,j,k,iBlock) = 0.0
    end do; end do; end do

  end subroutine user_set_ics

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

    real :: r_D(3),rhovr_D(3), rhovr
    integer :: i, j, k, iVar
    character (len=*), parameter :: NameSub = 'user_set_ics'
    !----------------------------------------------------------------------

    call update_states_MHD(iStage,iBlock)

    if(Rmin_BLK(iBlock) <= planetR) then
       do k=1,nK; do j=1,nJ; do i=1,nI
          if(R_BLK(i+2,j,k,iBlock) >= planetR ) CYCLE
          State_VGB(Rho_,i,j,k,iBlock) = planetDensity
          State_VGB(P_,i,j,k,iBlock) = planetPressure
          State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock) = 0.0
       end do; end do; end do



       do k=MinK,MaxK; do j=MinJ,MaxJ; do i=1,nI
          if(R_BLK(i-1,j,k,iBlock) > planetR ) CYCLE
          if(R_BLK(i,j,k,iBlock) <= planetR ) CYCLE

          r_D  = (/ X_BLK(i,j,k,iBlock),Y_BLK(i,j,k,iBlock),Z_BLK(i,j,k,iBlock)/)/R_BLK(i,j,k,iBlock)

          rhovr = dot_product(State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock),r_D)
          if(rhovr > 0.0) then
             ! If low out of the planet boundary we set the radial componet to zero for the momentum inside
             rhovr_D = r_D*rhovr
          else
             ! If flow into the planet it will be absourbed
             rhovr_D = (/0.0,0.0,0.0/)
          end if

          ! Two cell boundary layer inside the planet
          !! -1
          State_VGB(Rho_,i-1,j,k,iBlock) = State_VGB(Rho_,i,j,k,iBlock)
          State_VGB(RhoUx_:RhoUz_,i-1,j,k,iBlock) = State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) - rhovr_D
          State_VGB(Bz_+1:nVar,i-1,j,k,iBlock) = State_VGB(Bz_+1:nVar,i,j,k,iBlock)
          !! -2
          State_VGB(Rho_,i-2,j,k,iBlock) = State_VGB(Rho_,i,j,k,iBlock)
          State_VGB(RhoUx_:RhoUz_,i-2,j,k,iBlock) = State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) - rhovr_D
          State_VGB(Bz_+1:nVar,i-2,j,k,iBlock) = State_VGB(Bz_+1:nVar,i,j,k,iBlock)

       end do; end do; end do

    end if

    call calc_energy_cell(iBlock)

  end subroutine user_update_states

  subroutine user_set_resistivity(iBlock, Eta_G)
    use ModEnergy
    use BATL_size, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModGeometry,   ONLY: R_BLK, Rmin_BLK
    use ModResistivity, ONLY: Eta0
    integer, intent(in) :: iBlock
    real, intent(out) :: Eta_G(MinI:MaxI, MinJ:MaxJ, MinK:MaxK) 

    integer ::i,j,k

    Eta_G = Eta0

    if(Rmin_BLK(iBlock) > planetR) return

    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
       if(R_BLK(i,j,k,iBlock) > planetResistivR_I(1) ) CYCLE
       if(R_BLK(i,j,k,iBlock) < planetResistivR_I(2) ) CYCLE
       ! to avoid eta jumps adding Eta_G
       Eta_G(i,j,k) = Eta_G(i,j,k) + (planetResistivR_I(1) - R_BLK(i,j,k,iBlock))* &
            ResistivetyRate
    end do; end do; end do

  end subroutine user_set_resistivity

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModSize, ONLY: nI, nJ, nK

    integer,          intent(in)   :: iBlock
    character(len=*), intent(in)   :: NameVar
    logical,          intent(in)   :: IsDimensional
    real,             intent(out)  :: PlotVar_G(-1:nI+2, -1:nJ+2, -1:nK+2)
    real,             intent(out)  :: PlotVarBody
    logical,          intent(out)  :: UsePlotVarBody
    character(len=*), intent(inout):: NameTecVar
    character(len=*), intent(inout):: NameTecUnit
    character(len=*), intent(inout):: NameIdlUnit
    logical,          intent(out)  :: IsFound

    character (len=*), parameter :: NameSub = 'user_set_plot_var'
    !-------------------------------------------------------------------

    !  The planet 2-cell wide layer is included in the computational domain,
    !  but are not true values. This boundary layer shod be masked or the 
    !  planet body stretched out to cover this layer


  end subroutine user_set_plot_var

end module ModUser
