!#NOTPUBLIC  email:xzjia@umich.edu  expires:12/31/2099
!This code is a copyright protected software (c) 2002- University of Michigan
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
  character (len=*), parameter :: NameUserModule = 'Mercury, Lars Daldorff & Xianzhe Jia'

  real :: PlanetDensity=-1., PlanetPressure=-1., PlanetRadius=-1.

  real :: PlanetDensitySi=-1., PlanetPressureSi=-1., PlanetRadiusSi=-1.

  integer :: nLayer =0 ! Number of points in planet resistivity profile
  real, allocatable :: PlanetRadiusSi_I(:),PlanetRadius_I(:),&
       ResistivetySi_I(:),Resistivety_I(:)
  real, allocatable :: ResistivityRate(:)

  integer :: iLayer

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
    use ModGeometry,    ONLY: TypeGeometry
    CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(A22,E10.3)"
    !-------------------------------------------------------------------

!!! Introduce PlanetDensitySi etc., read those and convert 
!!! from that. Initialize these to -1. PlanetRadius should be
!!! in dimensional units.

    if (TypeGeometry /= 'spherical_lnr') &
         call stop_mpi('ERROR: Correct PARAM.in, need spherical grid.')


    if(PlanetDensitySi < 0.0) &
         PlanetDensitySI  = 3.0*MassPlanet/(4.0*cPi*RadiusPlanet**3)

    if(PlanetRadius < 0.0) &
         PlanetRadius = RadiusPlanet*Si2No_V(UnitX_)

    if(PlanetPressureSi < 0.0) &
         PlanetPressureSi = 1.0e-8*No2Si_V(UnitP_)

 
    PlanetDensity           = PlanetDensitySi*Si2No_V(UnitRho_)
    PlanetPressure          = PlanetPressureSi*Si2No_V(UnitP_)

    PlanetRadiusSi          = PlanetRadius*No2Si_V(UnitX_)

    if(nLayer > 1) then
       PlanetRadiusSi_I = PlanetRadius_I*No2Si_V(UnitX_)
       Resistivety_I = ResistivetySi_I*Si2NoEta
       do iLayer=2,nLayer
          ResistivityRate(iLayer-1) = &
               (Resistivety_I(iLayer) -Resistivety_I(iLayer-1))/&
               (PlanetRadius_I(iLayer) - PlanetRadius_I(iLayer-1))
       end do
    end if

    print *," Rate eta", ResistivityRate


     if(iProc==0) then
       call write_myname
       write(*,*) ''
       write(*,*) '   Resistiv Planet Model'       
       write(*,*) '   ---------------------'       
       write(*,*) ''
       write(*,FMT1) '  Planet density  = ',PlanetDensitySi
       write(*,FMT1) '  Planet pressure = ',PlanetPressureSi
       write(*,FMT1) '  Planet radius   = ',PlanetRadiusSi
       if(nLayer > 0 ) then
          write(*,*) ''
          write(*,*) '   |-------- Planet Resistivety Profile -----|'
          write(*,*) '       Radius(SI)            Resistivety(SI)'
          do iLayer =1,nLayer 
             write(*,"(A7,E10.3,A15,E10.3)") " ",PlanetRadiusSi_I(iLayer)," ",&
                  ResistivetySi_I(iLayer)
          end do
       else
          write(*,*) 'Conducting Planet (eta =0)'
       end if
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
          call read_var('PlanetRadius'        , PlanetRadius)

          call read_var('nResistivPoints', nLayer)
          if(nLayer == 1) then
             write(*,*) ' We need minimum 2 points for including resistivety profile'
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if

          if(nLayer > 1) then
             allocate(ResistivityRate(nLayer-1),&
                  PlanetRadiusSi_I(nLayer),&
                  PlanetRadius_I(nLayer), &
                  ResistivetySi_I(nLayer),&
                  Resistivety_I(nLayer))

             do iLayer=1,nLayer
                call read_var('Radius',PlanetRadius_I(iLayer))
                call read_var('Resistivety', ResistivetySi_I(iLayer))
             end do

             


             !Check values
             do iLayer=2,nLayer
                if(PlanetRadius_I(iLayer-1) < &
                     PlanetRadius_I(iLayer)) then
                   write(*,*) 'ERROR: Shoud be decreasing Radius.'
                   call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
                end if
             end do
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
  
  subroutine user_set_ics(iBlock)

    use ModAdvance,    ONLY: State_VGB
    use ModMain,       ONLY: Unused_B, nBlockMax
    use ModGeometry,   ONLY: R_BLK
    use ModSize,       ONLY: nI, nJ, nK, nG
    use ModVarIndexes, ONLY: Rho_, rhoUx_, rhoU_, P_,Bx_,Bz_
    use ModNumConst,   ONLY:cPi
    use CON_planet,    ONLY: RadiusPlanet, MassPlanet
    use ModPhysics,    ONLY: Si2No_V,UnitRho_ 
    use ModMultiFluid, ONLY: select_fluid, iFluid, nFluid, iP, iEnergy, &
         iRho, iRhoUx, iRhoUy, iRhoUz

    integer, intent(in) :: iBlock

    integer :: i,j,k
    character (len=*), parameter :: NameSub = 'user_set_ics'
    !-------------------------------------------------------------------

    if(R_BLK(1,1,1,iBlock) > PlanetRadius) RETURN

    do iFluid = 1, nFluid
       call select_fluid
       do k=1,nK; do j=1,nJ; do i=1,nI
          if(R_BLK(i+nG,j,k,iBlock) > PlanetRadius) CYCLE
          State_VGB(iRho,i,j,k,iBlock) = PlanetDensity
          State_VGB(iP,i,j,k,iBlock)   = PlanetPressure
          State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = 0.0
       end do; end do; end do
    end do

    do k=1,nK; do j=1,nJ; do i=1,nI
       if(R_BLK(i,j,k,iBlock) > PlanetRadius) CYCLE
       State_VGB(Bx_:Bz_,i,j,k,iBlock) = 0.0
    end do; end do; end do

  end subroutine user_set_ics

  !===========================================================================

  subroutine user_update_states(iStage,iBlock)

    use CON_planet,    ONLY: RadiusPlanet, MassPlanet
    use ModAdvance,    ONLY: State_VGB
    use ModMain,       ONLY: Unused_B, nBlockMax
    use ModGeometry,   ONLY: Xyz_DGB, R_BLK, Rmin_BLK
    use ModSize,       ONLY: nI, nJ, nK, nG
    use ModPhysics,    ONLY: Si2No_V,UnitRho_ 
    use ModVarIndexes, ONLY: nVar, Rho_, RhoU_, RhoUx_, RhoUz_, Bx_, Bz_, p_
    use ModEnergy,     ONLY: calc_energy_cell

    integer,intent(in) :: iStage, iBlock

    real :: r_D(3), dRhoUr_D(3), RhoUr
    integer :: i, iG, j, k, iVar
    character (len=*), parameter :: NameSub = 'user_set_ics'
    !----------------------------------------------------------------------

    call update_states_MHD(iStage,iBlock)

    if(Rmin_BLK(iBlock) <= PlanetRadius) then
       do i = 1, nI
          if(R_BLK(i+nG,1,1,iBlock) >= PlanetRadius) CYCLE
          do k = 1, nK; do j = 1, nJ; 
             ! Set density, pressure and momentum inside the planet
             ! and the nG ghost cells to fixed values.
             State_VGB(Rho_,i,j,k,iBlock) = PlanetDensity
             State_VGB(P_,i,j,k,iBlock)   = PlanetPressure
             State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
          end do; end do
       end do

       if(r_BLK(MaxI,1,1,iBlock) >= PlanetRadius)then
          do i = MaxI, 1, -1
            ! Find the i index just outside the planet radius
            if(R_BLK(i-1,1,1,iBlock) > PlanetRadius ) CYCLE
            EXIT
          end do
          do k = MinK, MaxK; do j = MinJ, MaxJ

             ! Get radial momentum
             r_D = Xyz_DGB(x_:z_,i,j,k,iBlock)/ r_BLK(i,j,k,iBlock)
             RhoUr = dot_product(State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock), r_D)

             if(RhoUr > 0.0) then
                ! If flow is out of the planet, remove the radial component 
                ! of the momentum so that the flow is tangential
                dRhoUr_D = -r_D*RhoUr
             else
                ! If flow is into the planet do nothing so the flow is absorbed
                dRhoUr_D = 0.0
             end if

             ! Set nG cells inside the planet as a boundary condition
             do iG = i-nG, i-1
                State_VGB(Rho_,iG,j,k,iBlock) = State_VGB(Rho_,i,j,k,iBlock)
                State_VGB(P_,iG,j,k,iBlock)   = State_VGB(P_,i,j,k,iBlock)
		State_VGB(RhoUx_:RhoUz_,iG,j,k,iBlock) = &
                     State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) + dRhoUr_D
             end do
          end do; end do

       end if

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

    integer ::i,j,k,iLayer
    !-------------------------------------------------------------------
    Eta_G = Eta0

    if(nLayer <2 ) RETURN
    if(Rmin_BLK(iBlock) > PlanetRadius_I(1)) RETURN

    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
       do iLayer=nLayer-1,1,-1
          if(R_BLK(i,j,k,iBlock) < PlanetRadius_I(iLayer+1) ) CYCLE
          if(R_BLK(i,j,k,iBlock) > PlanetRadius_I(iLayer) ) CYCLE
          ! to avoid eta jumps adding Eta_G
          Eta_G(i,j,k) = Eta_G(i,j,k) +Resistivety_I(iLayer+1)+(R_BLK(i,j,k,iBlock)-PlanetRadius_I(iLayer+1))* &
               ResistivityRate(iLayer)
       end do
    end do; end do; end do

  end subroutine user_set_resistivity

end module ModUser
