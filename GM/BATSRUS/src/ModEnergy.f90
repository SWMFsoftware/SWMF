!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModEnergy

  use BATL_lib, ONLY: &
       test_start, test_stop, iTest, jTest, kTest, iBlockTest

  use ModProcMH,     ONLY: iProc
  use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUz_, Bx_, Bz_, Hyp_, p_, Pe_, IsMhd
  use ModMultiFluid, ONLY: nFluid, iFluid, IonLast_, &
       iRho, iRhoUx, iRhoUz, iP, iP_I, DoConserveNeutrals, &
       select_fluid, MassFluid_I, iRho_I, iRhoIon_I, MassIon_I, ChargeIon_I
  use ModSize,       ONLY: nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK
  use ModAdvance,    ONLY: State_VGB, Energy_GBI, StateOld_VGB, EnergyOld_CBI,&
       UseNonConservative, nConservCrit, IsConserv_CB, UseElectronPressure
  use ModPhysics,    ONLY: GammaMinus1_I, InvGammaMinus1_I, InvGammaMinus1, &
       pMin_I, PeMin, Tmin_I, TeMin

  implicit none

  private ! except

  public:: calc_energy_or_pressure  ! p -> e or e -> p
  public:: calc_pressure            ! e -> p in a range of cells/fluids
  public:: calc_pressure_point      ! e -> p in 1 cell for all fluids
  public:: calc_pressure_cell       ! e -> p in physical cells
  public:: calc_energy              ! p -> e in a range of cells/fluids
  public:: calc_energy_point        ! p -> e in one cell for all fluids
  public:: calc_energy_cell         ! p -> e in physical cells
  public:: calc_energy_ghost        ! p -> e in ghost cells
  public:: calc_old_pressure        ! eold -> pold (ModPartImplicit)
  public:: calc_old_energy          ! pold -> eold (ModPartImplicit)
  public:: correctP                 ! obsolete method for CT and projection

  ! It is possible to add the hyperbolic scalar energy to the total energy
  ! (Tricco & Price 2012, J. Comput. Phys. 231, 7214).
  ! Experiments so far do not show any benefit, but the code is preserved.
  logical, parameter:: UseHypEnergy = .false.

contains
  !============================================================================

  subroutine calc_energy_or_pressure(iBlock)

    ! Calculate pressure from energy or energy from pressure depending on
    ! the value of UseNonConservative and IsConserv_CB

    integer, intent(in) :: iBlock
    integer::i,j,k
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'calc_energy_or_pressure'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    if(DoTest)write(*,*)NameSub, &
         ': UseNonConservative, DoConserveNeutrals, nConservCrit=', &
         UseNonConservative, DoConserveNeutrals, nConservCrit

    if(.not. UseNonConservative)then
       if(DoConserveNeutrals) then
          ! All cells are conservative
          call calc_pressure_cell(iBlock)
       else
          ! Ions are conservative, neutrals are non-conservative
          call calc_pressure(1, nI, 1, nJ, 1, nK, iBlock, 1         , IonLast_)
          call calc_energy(  1, nI, 1, nJ, 1, nK, iBlock, IonLast_+1, nFluid)
       end if
       RETURN
    end if

    if(UseNonConservative .and. nConservCrit <= 0)then
       ! All cells are non-conservative
       call calc_energy_cell(iBlock)
       RETURN
    end if

    ! Make sure pressure used for energy is larger than floor value
    call limit_pressure(1, nI, 1, nJ, 1, nK, iBlock, 1, nFluid)

    ! A mix of conservative and non-conservative scheme (at least for the ions)
    FLUIDLOOP: do iFluid = 1, nFluid

       if(iFluid > IonLast_ .and. .not. DoConserveNeutrals) then
          ! Do all neutrals non-conservative and exit from the loop
          call calc_energy(1, nI, 1, nJ, 1, nK, iBlock, iFluid, nFluid)
          EXIT FLUIDLOOP
       end if

       call select_fluid
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(IsConserv_CB(i,j,k,iBlock)) then
             State_VGB(iP,i,j,k,iBlock) =                                 &
                  GammaMinus1_I(iFluid)*( Energy_GBI(i,j,k,iBlock,iFluid) &
                  - 0.5*sum(State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock)**2)     &
                  /State_VGB(iRho,i,j,k,iBlock) )
          else
             Energy_GBI(i,j,k,iBlock,iFluid) =                        &
                  InvGammaMinus1_I(iFluid)*State_VGB(iP,i,j,k,iBlock) &
                  + 0.5*sum(State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock)**2) &
                  /State_VGB(iRho,i,j,k,iBlock)
          end if
       end do; end do; end do

       if(iFluid > 1 .or. .not. IsMhd) CYCLE FLUIDLOOP

       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(IsConserv_CB(i,j,k,iBlock)) then
             State_VGB(iP,i,j,k,iBlock) = State_VGB(iP,i,j,k,iBlock) &
                  - GammaMinus1_I(iFluid)*0.5 * &
                  sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)**2)
          else
             Energy_GBI(i,j,k,iBlock,iFluid) = &
                  Energy_GBI(i,j,k,iBlock,iFluid) + &
                  0.5*sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)**2)
          end if
       end do; end do; end do

       if(Hyp_ > 1 .and. UseHypEnergy)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(IsConserv_CB(i,j,k,iBlock)) then
                State_VGB(iP,i,j,k,iBlock) = State_VGB(iP,i,j,k,iBlock) &
                     - GammaMinus1_I(iFluid)*0.5 * &
                     State_VGB(Hyp_,i,j,k,iBlock)**2
             else
                Energy_GBI(i,j,k,iBlock,iFluid) = &
                     Energy_GBI(i,j,k,iBlock,iFluid) + &
                     0.5*State_VGB(Hyp_,i,j,k,iBlock)**2
             end if
          end do; end do; end do
       endif

    end do FLUIDLOOP

    ! Make sure pressure used for energy is larger than floor value
    call limit_pressure(1, nI, 1, nJ, 1, nK, iBlock, 1, nFluid, &
         DoUpdateEnergy = .true.)

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine calc_energy_or_pressure
  !============================================================================

  subroutine calc_old_pressure(iBlock)

    ! Calculate pressure from energy for the old state
    !
    !   P = (gamma-1)*(E - 0.5*rho*u^2 - 0.5*b1^2)
    !

    integer, intent(in) :: iBlock
    integer :: i, j, k

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'calc_old_pressure'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    do iFluid = 1, nFluid
       call select_fluid
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          StateOld_VGB(iP, i, j, k,iBlock) = &
               GammaMinus1_I(iFluid)*                           &
               ( EnergyOld_CBI(i,j,k,iBlock,iFluid) - 0.5 *     &
               sum(StateOld_VGB(iRhoUx:iRhoUz,i,j,k,iBlock)**2) &
               /StateOld_VGB(iRho,i,j,k,iBlock) )
       end do; end do; end do

       if(iFluid > 1 .or. .not. IsMhd) CYCLE

       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          StateOld_VGB(iP,i,j,k,iBlock) = StateOld_VGB(iP,i,j,k,iBlock) &
               - GammaMinus1_I(iFluid)*0.5* &
               sum(StateOld_VGB(Bx_:Bz_,i,j,k,iBlock)**2)
       end do; end do; end do

       if(Hyp_ > 1 .and. UseHypEnergy)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             StateOld_VGB(iP,i,j,k,iBlock) = StateOld_VGB(iP,i,j,k,iBlock) &
                  - GammaMinus1_I(iFluid)*0.5* &
                  StateOld_VGB(Hyp_,i,j,k,iBlock)**2
          end do; end do; end do
       end if

    end do

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine calc_old_pressure
  !============================================================================

  subroutine calc_old_energy(iBlock)

    ! Calculate energy from pressure for the old state
    !
    !   E = P/(gamma-1) + 0.5*rho*u^2 + 0.5*b1^2

    integer, intent(in) :: iBlock
    integer :: i, j, k

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'calc_old_energy'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    do iFluid = 1, nFluid
       call select_fluid
       ! Calculate thermal plus kinetic energy
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if (StateOld_VGB(iRho,i,j,k,iBlock) <= 0.0)then
             EnergyOld_CBI(i,j,k,iBlock,iFluid) = 0.0
          else
             EnergyOld_CBI(i,j,k,iBlock,iFluid) =                 &
                  InvGammaMinus1_I(iFluid)*StateOld_VGB(iP,i,j,k,iBlock) &
                  + 0.5*(sum(StateOld_VGB(iRhoUx:iRhoUz,i,j,k,iBlock)**2)/&
                  StateOld_VGB(iRho,i,j,k,iBlock))
          end if
       end do; end do; end do

       if(iFluid > 1 .or. .not. IsMhd) CYCLE

       ! Add magnetic energy for ion fluid
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          EnergyOld_CBI(i,j,k,iBlock,iFluid) = &
               EnergyOld_CBI(i,j,k,iBlock,iFluid) + &
               0.5*sum(StateOld_VGB(Bx_:Bz_,i,j,k,iBlock)**2)
       end do; end do; end do

       if(Hyp_ > 1 .and. UseHypEnergy)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             EnergyOld_CBI(i,j,k,iBlock,iFluid) = &
                  EnergyOld_CBI(i,j,k,iBlock,iFluid) + &
                  0.5*StateOld_VGB(Hyp_,i,j,k,iBlock)**2
          end do; end do; end do
       end if

    end do

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine calc_old_energy
  !============================================================================

  subroutine calc_pressure(iMin, iMax, jMin, jMax, kMin, kMax, iBlock, &
       iFluidMin, iFluidMax)

    ! Calculate pressure from energy
    !
    !   P = (gamma-1)*(E - 0.5*rho*u^2 - 0.5*b1^2)
    !

    integer, intent(in) :: iMin, iMax, jMin, jMax, kMin, kMax, iBlock
    integer, intent(in) :: iFluidMin, iFluidMax
    integer :: i, j, k
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'calc_pressure'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    if(DoTest)write(*,*)NameSub, &
         ': iMin,iMax,jMin,jMax,kMin,kMax,iFluidMin,iFluidMax=', &
         iMin,iMax,jMin,jMax,kMin,kMax,iFluidMin,iFluidMax

    do iFluid = iFluidMin, iFluidMax
       call select_fluid
       do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
          State_VGB(iP,i,j,k,iBlock) = &
               GammaMinus1_I(iFluid)*                        &
               (Energy_GBI(i,j,k,iBlock,iFluid) - 0.5*       &
               sum(State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock)**2) &
               /State_VGB(iRho,i,j,k,iBlock) )
       end do; end do; end do

       if(iFluid > 1 .or. .not. IsMhd) CYCLE

       ! Subtract magnetic energy
       do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
          State_VGB(iP,i,j,k,iBlock) = State_VGB(iP,i,j,k,iBlock) &
               - GammaMinus1_I(iFluid)*0.5* &
               sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)**2)
       end do; end do; end do

       if(Hyp_ > 1 .and. UseHypEnergy)then
          ! Subtract energy associated with the Hyp scalar
          do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
             State_VGB(iP,i,j,k,iBlock) = State_VGB(iP,i,j,k,iBlock) &
                  - GammaMinus1_I(iFluid)*0.5*State_VGB(Hyp_,i,j,k,iBlock)**2
          end do; end do; end do
       end if
    end do

    call limit_pressure(iMin, iMax, jMin, jMax, kMin, kMax, &
         iBlock, iFluidMin, iFluidMax, DoUpdateEnergy = .true.)

    if(DoTest)then
       write(*,*)NameSub,':Energy_GBI=',Energy_GBI(iTest,jTest,kTest,iBlock,:)
       write(*,*)NameSub,':State_VGB=',State_VGB(:,iTest,jTest,kTest,iBlock)
    end if

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine calc_pressure
  !============================================================================

  subroutine calc_energy(iMin, iMax, jMin, jMax, kMin, kMax, iBlock, &
       iFluidMin, iFluidMax)

    ! Calculate total energy (excluding B0):
    !
    !   E = p/(gamma-1) + 0.5*rho*u^2 + 0.5*b1^2
    !

    integer, intent(in) :: iMin, iMax, jMin, jMax, kMin, kMax, iBlock
    integer, intent(in) :: iFluidMin, iFluidMax
    integer::i,j,k
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'calc_energy'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    call limit_pressure(iMin, iMax, jMin, jMax, kMin, kMax, &
         iBlock, iFluidMin, iFluidMax)

    do iFluid = iFluidMin, iFluidMax
       call select_fluid
       ! Calculate thermal plus kinetic energy
       do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
          if (State_VGB(iRho,i,j,k,iBlock) <= 0.0)then
             Energy_GBI(i,j,k,iBlock,iFluid) = 0.0
          else
             Energy_GBI(i,j,k,iBlock,iFluid) = &
                  InvGammaMinus1_I(iFluid)*State_VGB(iP,i,j,k,iBlock) &
                  +0.5*(sum(State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock)**2)/&
                  State_VGB(iRho,i,j,k,iBlock))
          end if
       end do; end do; end do

       if(iFluid > 1 .or. .not. IsMhd) CYCLE

       ! Add magnetic energy for ion fluid
       do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
          Energy_GBI(i,j,k,iBlock,iFluid) = &
               Energy_GBI(i,j,k,iBlock,iFluid) + &
               0.5*sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)**2)
       end do; end do; end do

       if(Hyp_ > 1 .and. UseHypEnergy)then
          do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
             Energy_GBI(i,j,k,iBlock,iFluid) = &
                  Energy_GBI(i,j,k,iBlock,iFluid) + &
                  0.5*State_VGB(Hyp_,i,j,k,iBlock)**2
          end do; end do; end do
       end if

    end do

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine calc_energy
  !============================================================================

  subroutine calc_pressure_cell(iBlock)
    integer, intent(in) :: iBlock
    !--------------------------------------------------------------------------
    call calc_pressure(1,nI,1,nJ,1,nK,iBlock,1,nFluid)
  end subroutine calc_pressure_cell
  !============================================================================

  subroutine calc_pressure_ghost(iBlock)
    integer, intent(in) :: iBlock
    !--------------------------------------------------------------------------
    call calc_pressure(MinI,MaxI,MinJ,MaxJ,MinK,MaxK,iBlock,1,nFluid)
  end subroutine calc_pressure_ghost
  !============================================================================

  subroutine calc_pressure_point(i, j, k, iBlock)
    integer, intent(in) :: i, j, k, iBlock
    !--------------------------------------------------------------------------
    call calc_pressure(i,i,j,j,k,k,iBlock,1,nFluid)
  end subroutine calc_pressure_point
  !============================================================================

  subroutine calc_energy_cell(iBlock)
    integer, intent(in) :: iBlock
    !--------------------------------------------------------------------------
    call calc_energy(1,nI,1,nJ,1,nK,iBlock,1,nFluid)
  end subroutine calc_energy_cell
  !============================================================================

  subroutine calc_energy_ghost(iBlock, DoResChangeOnlyIn)

    use BATL_lib, ONLY: DiLevelNei_IIIB

    integer, intent(in) :: iBlock
    logical, optional, intent(in) :: DoResChangeOnlyIn

    integer :: i, j, k
    logical :: DoResChangeOnly

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'calc_energy_ghost'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    DoResChangeOnly =.false.
    if(present(DoResChangeOnlyIn)) DoResChangeOnly = DoResChangeOnlyIn

    if( DoResChangeOnly ) then
       if( .not.any(abs(DiLevelNei_IIIB(-1:1,-1:1,-1:1,iBlock)) == 1) ) RETURN
    end if

    call limit_pressure(MinI, MaxI, MinJ, MaxJ, MinK, MaxK, iBlock, 1, nFluid)

    do iFluid = 1, nFluid

       call select_fluid

       if(IsMhd .and. iFluid == 1) then
          ! MHD energy
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             if(State_VGB(Rho_,i,j,k,iBlock) <= 0.0)then
                Energy_GBI(i,j,k,iBlock,iFluid) = 0.0
             else
                Energy_GBI(i,j,k,iBlock,iFluid) = &
                     InvGammaMinus1*State_VGB(p_,i,j,k,iBlock) &
                     + 0.5*sum(State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock)**2) &
                     /State_VGB(Rho_,i,j,k,iBlock) + &
                     0.5*sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)**2)

                if(Hyp_ > 1 .and. UseHypEnergy) &
                     Energy_GBI(i,j,k,iBlock,iFluid) = &
                     Energy_GBI(i,j,k,iBlock,iFluid) &
                     + 0.5*State_VGB(Hyp_,i,j,k,iBlock)**2

             end if
          end do; end do; end do
       else
          ! HD energy
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             if(State_VGB(iRho,i,j,k,iBlock) <= 0.0)then
                Energy_GBI(i,j,k,iBlock,iFluid) = 0.0
             else
                Energy_GBI(i,j,k,iBlock,iFluid) = &
                     InvGammaMinus1_I(iFluid)*State_VGB(iP,i,j,k,iBlock) &
                     + 0.5*sum(State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock)**2) &
                     /State_VGB(iRho,i,j,k,iBlock)
             end if
          end do; end do; end do
       end if
    end do

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine calc_energy_ghost
  !============================================================================

  subroutine calc_energy_point(i, j, k, iBlock)
    integer, intent(in) :: i, j, k, iBlock
    !--------------------------------------------------------------------------
    call calc_energy(i,i,j,j,k,k,iBlock,1,nFluid)
  end subroutine calc_energy_point
  !============================================================================

  subroutine limit_pressure(iMin, iMax, jMin, jMax, kMin, kMax, iBlock, &
       iFluidMin, iFluidMax, DoUpdateEnergy)

    ! Keep pressure(s) in State_VGB above pMin_I limit
    ! If DoUpdateEnergy is present, also modify energy to remain consistent

    use ModAdvance, ONLY: UseAnisoPressure
    use ModMultiFluid, ONLY: IsIon_I, iPparIon_I

    integer, intent(in) :: iMin, iMax, jMin, jMax, kMin, kMax, iBlock
    integer, intent(in) :: iFluidMin, iFluidMax
    logical, intent(in), optional:: DoUpdateEnergy ! if present should be true

    integer:: i, j, k
    real :: NumDens, p, pMin, Ne

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'limit_pressure'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    do iFluid = iFluidMin, iFluidMax
       if(pMin_I(iFluid) < 0.0) CYCLE
       pMin = pMin_I(iFluid)
       iP = iP_I(iFluid)
       if(present(DoUpdateEnergy))then
          do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
             p = State_VGB(iP,i,j,k,iBlock)
             if(p < pMin)then
                State_VGB(iP,i,j,k,iBlock) = pMin
                Energy_GBI(i,j,k,iBlock,iFluid) =      &
                     Energy_GBI(i,j,k,iBlock,iFluid) + &
                     InvGammaMinus1_I(iFluid)*(pMin - p)
             end if
             if(UseAnisoPressure .and. IsIon_I(iFluid)) &
                  State_VGB(iPparIon_I(iFluid),i,j,k,iBlock) = &
                  max(pMin, State_VGB(iPparIon_I(iFluid),i,j,k,iBlock))
          end do; end do; end do
       else
          do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
             State_VGB(iP,i,j,k,iBlock) = max(pMin, State_VGB(iP,i,j,k,iBlock))
             if(UseAnisoPressure .and. IsIon_I(iFluid))&
                  State_VGB(iPparIon_I(iFluid),i,j,k,iBlock) = &
                  max(pMin, State_VGB(iPparIon_I(iFluid),i,j,k,iBlock))
          end do; end do; end do
       end if
    end do

    do iFluid = iFluidMin, iFluidMax
       if(Tmin_I(iFluid) < 0.0) CYCLE
       iP = iP_I(iFluid)
       if(present(DoUpdateEnergy))then
          do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
             p = State_VGB(iP,i,j,k,iBlock)
             NumDens=State_VGB(iRho_I(iFluid),i,j,k,iBlock)/MassFluid_I(iFluid)
             pMin = NumDens*Tmin_I(iFluid)
             if(p < pMin)then
                State_VGB(iP,i,j,k,iBlock) = pMin
                Energy_GBI(i,j,k,iBlock,iFluid) =      &
                     Energy_GBI(i,j,k,iBlock,iFluid) + &
                     InvGammaMinus1_I(iFluid)*(pMin - p)
             end if
             if(UseAnisoPressure .and. IsIon_I(iFluid))&
                  State_VGB(iPparIon_I(iFluid),i,j,k,iBlock) = &
                  max(pMin, State_VGB(iPparIon_I(iFluid),i,j,k,iBlock))
          end do; end do; end do
       else
          do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
             NumDens=State_VGB(iRho_I(iFluid),i,j,k,iBlock)/MassFluid_I(iFluid)
             pMin = NumDens*Tmin_I(iFluid)
             State_VGB(iP,i,j,k,iBlock) = max(pMin, State_VGB(iP,i,j,k,iBlock))
             if(UseAnisoPressure .and. IsIon_I(iFluid))&
                  State_VGB(iPparIon_I(iFluid),i,j,k,iBlock) = &
                  max(pMin, State_VGB(iPparIon_I(iFluid),i,j,k,iBlock))
          end do; end do; end do
       end if
    end do

    if(UseElectronPressure .and. iFluidMin == 1)then
       if(PeMin > 0.0)then
          do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
             State_VGB(Pe_,i,j,k,iBlock) = &
                  max(PeMin, State_VGB(Pe_,i,j,k,iBlock))
          end do; end do; end do
       end if
       if(TeMin > 0.0)then
          do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
             Ne = sum(ChargeIon_I*State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I)
             State_VGB(Pe_,i,j,k,iBlock) = &
                  max(Ne*TeMin, State_VGB(Pe_,i,j,k,iBlock))
          end do; end do; end do
       end if
    end if

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine limit_pressure
  !============================================================================

  ! moved from file exchange_messages.f90
  subroutine correctP(iBlock)

    ! Make pressure and energy consistent and maintain thermal energy ratio
    ! at a reasonable value (this also excludes negative pressure)

    use ModProcMH
    use ModMain,       ONLY: nI,nJ,nK
    use ModVarIndexes, ONLY: rho_, rhoUx_, rhoUy_, rhoUz_, Bx_, By_, Bz_, P_
    use ModAdvance,    ONLY: State_VGB, Energy_GBI
    use ModPhysics,    ONLY: GammaMinus1, InvGammaMinus1, Pratio_hi, Pratio_lo
    use ModGeometry,   ONLY: true_cell
    use BATL_lib,      ONLY: Xyz_DGB

    integer, intent(in) :: iBlock

    integer :: i,j,k
    real :: inv_dratio, qp, qe, qth, qratio, qd, qde, qpmin, &
         qdesum, qdesumabs, qderelmax

    real, dimension(1:nI,1:nJ,1:nK) :: P_old

    integer :: ierror1=-1, ierror2=-1, ierror3=-1, loc(3)

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'correctP'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    if(iBlock==iBlockTest)then
    else
       DoTest=.false.; DoTest=.false.
    end if

    qpmin=1.
    qdesum=0.
    qdesumabs=0.
    qderelmax=0.

    P_old=State_VGB(P_,1:nI,1:nJ,1:nK,iBlock)

    inv_dratio=1./(Pratio_hi-Pratio_lo)

    do k=1,nK; do j=1,nJ; do i=1,nI

       if(.not.true_cell(i,j,k,iBlock))CYCLE

       ! Pressure and total energy
       qp=P_old(i,j,k)
       qe=Energy_GBI(i,j,k,iBlock,1)

       if(DoTest.and.i==iTest.and.J==jTest.and.K==kTest)&
            write(*,*)'CorrectP at me,BLK,i,j,k=',&
            iProc,iBlockTest,iTest,jTest,kTest, &
            ', initial P,E=',qp,qe

       ! Memorize smallest pressure
       qpmin=min(qp,qpmin)

       ! Thermal energy
       qth=InvGammaMinus1*qp

       ! Deviation=extra total energy=qe-inv_gm1*qp-(rhoU**2/rho+B**2)/2
       qd=qE-qth                                                         &
            -0.5*(State_VGB(rhoUx_,i,j,k,iBlock)**2+                         &
            State_VGB(rhoUy_,i,j,k,iBlock)**2+                               &
            State_VGB(rhoUz_,i,j,k,iBlock)**2)/State_VGB(rho_,i,j,k,iBlock)  &
            -0.5*(State_VGB(Bx_,i,j,k,iBlock)**2+                            &
            State_VGB(By_,i,j,k,iBlock)**2+                                  &
            State_VGB(Bz_,i,j,k,iBlock)**2)

       ! Limited thermal/total energy ratio for correction
       qratio=min(Pratio_hi,max(Pratio_lo,min(qth,qth+qd)/qe))

       ! Total energy is modified by qde (=0 if qratio==Pratio_hi)
       qde=qd*(Pratio_hi-qratio)*inv_dratio

       ! Collect total energy change
       qdesum   =qdesum   -qde
       qdesumabs=qdesumabs+abs(qde)
       qderelmax=max(qderelmax,qde/qe)

       ! Pressure is modified
       State_VGB(P_,i,j,k,iBlock)=GammaMinus1*(qth+qd-qde)

       ! We should now have E=inv_gm1*P+(rhoU**2/rho+B**2)/2:
       !
       ! qp*inv_gm1+qd-qde + (rhoU**2/rho+B**2)/2 = qe-qde = E
       !
       ! Correct!

       if(DoTest.and.i==iTest.and.J==jTest.and.K==kTest)then
          write(*,*)'qp,qth,qe,qd,qratio,qde=',qp,qth,qe,qd,qratio,qde
          write(*,*)'CorrectP, final P=',State_VGB(P_,i,j,k,iBlock)
       end if

    end do; end do; end do

    if(qpmin<0.)then
       if(ierror1==-1)then
          loc=minloc(P_old)
          write(*,*)'Negative P at me,iBlock,I,J,K,x,y,z,val',&
               iProc,iBlock,loc,&
               Xyz_DGB(:,loc(1),loc(2),loc(3),iBlock), &
               P_old(loc(1),loc(2),loc(3))
       end if
       call error_report('Negative P in exchange msgs, min(P)', &
            qpmin,ierror1,.true.)
    end if
    if(qderelmax>1.0E-3)then
       call error_report('E change in exchange_msgs, dE',qdesum,ierror2,.false.)
       call error_report('|E| change in exchange_msgs, d|E|',qdesumabs,ierror3,&
            .false.)
    end if

    if(DoTest)write(*,*)'CorrectP qpmin=',qpmin

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine correctP
  !============================================================================

end module ModEnergy
!==============================================================================

