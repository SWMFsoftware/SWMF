module ModEnergy

  use ModProcMH,     ONLY: iProc
  use ModMain,       ONLY: BlkTest,iTest,jTest,kTest,ProcTest
  use ModMultiFluid
  use ModSize,       ONLY: nI, nJ, nK, gcn, MaxBlock
  use ModAdvance,    ONLY: State_VGB, Energy_GBI, StateOld_VCB, EnergyOld_CBI,&
       UseNonConservative, nConservCrit, IsConserv_CB, UseElectronPressure
  use ModPhysics,    ONLY: Gm1, Inv_Gm1, pMin_I
  use ModVarIndexes, ONLY: Pe_, WaveFirst_, WaveLast_
  use ModWaves,      ONLY: UseWavePressure

  implicit none

contains

  !============================================================================
  subroutine calc_energy_or_pressure(iBlock)

    ! Calculate pressure from energy or energy from pressure depending on 
    ! the value of UseNonConservative and IsConserv_CB

    integer, intent(in) :: iBlock
    integer::i,j,k
    logical:: DoTest,DoTestMe
    character(len=*),parameter:: NameSub='calc_energy_or_pressure'
    !--------------------------------------------------------------------------
    if(iBlock==BlkTest .and. iProc==ProcTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    if(DoTestMe)write(*,*)NameSub, &
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

    call limit_pressure(1, nI, 1, nJ, 1, nK, iBlock, iFluid, nFluid)

    ! A mix of conservative and non-conservative scheme (at least for the ions)
    FLUIDLOOP: do iFluid = 1, nFluid

       if(iFluid > IonLast_ .and. .not. DoConserveNeutrals) then
          ! Do all neutrals non-conservative and exit from the loop
          call calc_energy(1, nI, 1, nJ, 1, nK, iBlock, iFluid, nFluid)
          EXIT FLUIDLOOP
       end if

       call select_fluid
       do k=1, nK; do j=1, nJ; do i=1, nI
          if(IsConserv_CB(i,j,k,iBlock)) then
             State_VGB(iP,i,j,k,iBlock) =                             &
                  gm1*( Energy_GBI(i,j,k,iBlock,iFluid)               &
                  - 0.5*sum(State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock)**2) &
                  /State_VGB(iRho,i,j,k,iBlock) )
          else
             Energy_GBI(i,j,k,iBlock,iFluid) =                        &
                  inv_gm1*State_VGB(iP,i,j,k,iBlock)                  &
                  + 0.5*sum(State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock)**2) &
                  /State_VGB(iRho,i,j,k,iBlock)
          end if
       end do; end do; end do

       if(nIonFluid == 1 .and. iFluid == 1)then
          if(UseElectronPressure)then
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                if(IsConserv_CB(i,j,k,iBlock))then
                   State_VGB(iP,i,j,k,iBlock) = State_VGB(iP,i,j,k,iBlock) &
                        - State_VGB(Pe_,i,j,k,iBlock)
                else
                   Energy_GBI(i,j,k,iBlock,iFluid) = &
                        Energy_GBI(i,j,k,iBlock,iFluid) &
                        + inv_gm1*State_VGB(Pe_,i,j,k,iBlock)
                end if
             end do; end do; end do
          end if
          if(UseWavePressure)then
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                if(IsConserv_CB(i,j,k,iBlock))then
                   State_VGB(iP,i,j,k,iBlock) = State_VGB(iP,i,j,k,iBlock) &
                        - gm1*sum(State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock))
                else
                   Energy_GBI(i,j,k,iBlock,iFluid) = &
                        Energy_GBI(i,j,k,iBlock,iFluid) &
                        + sum(State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock))
                end if
             end do; end do; end do
          end if
       end if

       if(iFluid > 1 .or. .not. IsMhd) CYCLE FLUIDLOOP

       do k=1, nK; do j=1, nJ; do i=1, nI
          if(IsConserv_CB(i,j,k,iBlock)) then
             State_VGB(iP, i, j, k,iBlock) = State_VGB(iP, i, j, k,iBlock) &
                  - gm1*0.5*sum(State_VGB(Bx_:Bz_,i, j, k,iBlock)**2)
          else
             Energy_GBI(i, j, k, iBlock, iFluid) = &
                  Energy_GBI(i, j, k, iBlock, iFluid) + &
                  0.5*sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)**2)
          end if
       end do; end do; end do
    end do FLUIDLOOP

    call limit_pressure(1, nI, 1, nJ, 1, nK, iBlock, iFluid, nFluid)

  end subroutine calc_energy_or_pressure

  !ADJOINT SPECIFIC BEGIN
  !============================================================================
  subroutine calc_energy_or_pressure_adjoint(iBlock)

    ! Adjoint version of above

    integer, intent(in) :: iBlock
!!$    integer::i,j,k
!!$    logical:: DoTest,DoTestMe
!!$    character(len=*),parameter:: NameSub='calc_energy_or_pressure_adjoint'
!!$    !--------------------------------------------------------------------------
!!$
!!$    if(.not. UseNonConservative)then
!!$       if(DoConserveNeutrals) then
!!$          ! All cells are conservative
!!$          call calc_pressure_cell_adjoint(iBlock)
!!$       else
!!$          ! Ions are conservative, neutrals are non-conservative
!!$          call calc_pressure_adjoint(1, nI, 1, nJ, 1, nK, iBlock, 1         , IonLast_)
!!$          call calc_energy_adjoint(  1, nI, 1, nJ, 1, nK, iBlock, IonLast_+1, nFluid)
!!$       end if
!!$       RETURN
!!$    end if
!!$
!!$    if(UseNonConservative .and. nConservCrit <= 0)then
!!$       ! All cells are non-conservative
!!$       call calc_energy_cell_adjoint(iBlock)
!!$       RETURN
!!$    end if
!!$
!!$    call stop_mpi(" in calc_energy_or_pressure_adjoint: not yet implemented.")

  end subroutine calc_energy_or_pressure_adjoint
  !ADJOINT SPECIFIC END
  
  !============================================================================

  subroutine calc_old_pressure(iBlock)

    ! Calculate pressure from energy for the old state
    !
    !   P = (gamma-1)*(E - 0.5*rho*u^2 - 0.5*b1^2)
    !
    
    integer, intent(in) :: iBlock
    integer :: i, j, k
    !--------------------------------------------------------------------------
    do iFluid = 1, nFluid
       call select_fluid
       do k=1, nK; do j=1, nJ; do i=1, nI
          StateOld_VCB(iP, i, j, k,iBlock) = &
               gm1*(EnergyOld_CBI(i,j,k,iBlock,iFluid) - 0.5*   &
               sum(StateOld_VCB(iRhoUx:iRhoUz,i,j,k,iBlock)**2)  &
               /StateOld_VCB(iRho,i,j,k,iBlock) )
       end do; end do; end do

       if(nIonFluid == 1 .and. iFluid == 1)then
          if(UseElectronPressure)then
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                StateOld_VCB(iP,i,j,k,iBlock) = StateOld_VCB(iP,i,j,k,iBlock) &
                     - StateOld_VCB(Pe_,i,j,k,iBlock)
             end do; end do; end do
          end if
          if(UseWavePressure)then
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                StateOld_VCB(iP,i,j,k,iBlock) = StateOld_VCB(iP,i,j,k,iBlock) &
                     - gm1*sum(StateOld_VCB(WaveFirst_:WaveLast_,i,j,k,iBlock))
             end do; end do; end do
          end if
       end if

       if(iFluid > 1 .or. .not. IsMhd) CYCLE
       
       do k=1, nK; do j=1, nJ; do i=1, nI
          StateOld_VCB(iP,i,j,k,iBlock) = StateOld_VCB(iP,i,j,k,iBlock) &
               - gm1*0.5*sum(StateOld_VCB(Bx_:Bz_,i,j,k,iBlock)**2)
       end do; end do; end do
    end do

    call limit_old_pressure(iBlock)

  end subroutine calc_old_pressure

  !============================================================================

  subroutine calc_old_energy(iBlock)

    ! Calculate energy from pressure for the old state
    !
    !   E = P/(gamma-1) + 0.5*rho*u^2 + 0.5*b1^2
    
    integer, intent(in) :: iBlock
    integer :: i, j, k
    !--------------------------------------------------------------------------

    call limit_old_pressure(iBlock)

    do iFluid = 1, nFluid
       call select_fluid
       ! Calculate thermal plus kinetic energy
       do k=1, nK; do j=1, nJ; do i=1, nI
          if (StateOld_VCB(iRho,i,j,k,iBlock) <= 0.0)then
             EnergyOld_CBI(i,j,k,iBlock,iFluid) = 0.0
          else
             EnergyOld_CBI(i, j, k, iBlock, iFluid) = &
                  inv_gm1*StateOld_VCB(iP,i,j,k,iBlock) &
                  + 0.5*(sum(StateOld_VCB(iRhoUx:iRhoUz,i,j,k,iBlock)**2)/&
                  StateOld_VCB(iRho,i,j,k,iBlock))
          end if
       end do; end do; end do

       if(nIonFluid == 1 .and. iFluid == 1)then
          if(UseElectronPressure)then
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                EnergyOld_CBI(i,j,k,iBlock,iFluid) = &
                     EnergyOld_CBI(i,j,k,iBlock,iFluid) &
                     + inv_gm1*StateOld_VCB(Pe_,i,j,k,iBlock)
             end do; end do; end do
          end if
          if(UseWavePressure)then
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                EnergyOld_CBI(i,j,k,iBlock,iFluid) = &
                     EnergyOld_CBI(i,j,k,iBlock,iFluid) &
                     + sum(StateOld_VCB(WaveFirst_:WaveLast_,i,j,k,iBlock))
             end do; end do; end do
          end if
       end if
       
       if(iFluid > 1 .or. .not. IsMhd) CYCLE

       ! Add magnetic energy for ion fluid
       do k=1, nK; do j=1, nJ; do i=1, nI
          EnergyOld_CBI(i,j,k,iBlock, iFluid) = &
               EnergyOld_CBI(i,j,k,iBlock, iFluid) + &
               0.5*sum(StateOld_VCB(Bx_:Bz_,i,j,k,iBlock)**2)
       end do; end do; end do

    end do

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
    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub='calc_pressure'
    !-------------------------------------------------------------------------
    if(iBlock==BlkTest .and. iProc==ProcTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if

    if(DoTestMe)write(*,*)NameSub, &
         ': iMin,iMax,jMin,jMax,kMin,kMax,iFluidMin,iFluidMax=', &
         iMin,iMax,jMin,jMax,kMin,kMax,iFluidMin,iFluidMax

    do iFluid = iFluidMin, iFluidMax
       call select_fluid
       do k=kMin, kMax; do j=jMin, jMax; do i=iMin, iMax
          State_VGB(iP, i, j, k, iBlock) = &
               gm1*(Energy_GBI(i, j, k, iBlock, iFluid) - 0.5*   &
               sum(State_VGB(iRhoUx:iRhoUz,i, j, k, iBlock)**2)  &
               /State_VGB(iRho,i, j, k, iBlock) )
       end do; end do; end do

       if(nIonFluid == 1 .and. iFluid == 1)then
          if(UseElectronPressure)then
             do k=kMin, kMax; do j=jMin, jMax; do i=iMin, iMax
                State_VGB(iP,i,j,k,iBlock) = State_VGB(iP,i,j,k,iBlock) &
                     - State_VGB(Pe_,i,j,k,iBlock)
             end do; end do; end do
          end if
          if(UseWavePressure)then
             do k=kMin, kMax; do j=jMin, jMax; do i=iMin, iMax
                State_VGB(iP,i,j,k,iBlock) = State_VGB(iP,i,j,k,iBlock) &
                     - gm1*sum(State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock))
             end do; end do; end do
          end if
       end if

       if(iFluid > 1 .or. .not. IsMhd) CYCLE

       do k=kMin, kMax; do j=jMin, jMax; do i=iMin, iMax
          State_VGB(iP, i, j, k,iBlock) = State_VGB(iP, i, j, k,iBlock) &
               - gm1*0.5*sum(State_VGB(Bx_:Bz_,i, j, k,iBlock)**2)
       end do; end do; end do
    end do

    call limit_pressure(iMin, iMax, jMin, jMax, kMin, kMax, &
         iBlock, iFluidMin, iFluidMax)
    
    if(DoTestMe)then
       write(*,*)NameSub,':Energy_GBI=',Energy_GBI(iTest,jTest,kTest,iBlock,:)
       write(*,*)NameSub,':State_VGB=',State_VGB(:,iTest,jTest,kTest,iBlock)
    end if

  end subroutine calc_pressure

  !ADJOINT SPECIFIC BEGIN
  !===========================================================================
  
  subroutine calc_pressure_adjoint(iMin, iMax, jMin, jMax, kMin, kMax, iBlock,&
       iFluidMin, iFluidMax)

    ! Calculate pressure from energy, adjoint version
    !
    !   P = (gamma-1)*(E - 0.5*rho*u^2 - 0.5*b1^2)
    !

    use ModAdjoint, ONLY: Adjoint_VGB, AdjEnergy_GBI
   
    integer, intent(in) :: iMin, iMax, jMin, jMax, kMin, kMax, iBlock
    integer, intent(in) :: iFluidMin, iFluidMax
!!$    integer :: i, j, k, d  
!!$    logical:: DoTest, DoTestMe
!!$    character(len=*), parameter:: NameSub='calc_pressure_adjoint'
!!$    !--------------------------------------------------------------------------
!!$
!!$    do iFluid = iFluidMin, iFluidMax
!!$       call select_fluid
!!$       do k=kMin, kMax; do j=jMin, jMax; do i=iMin, iMax
!!$          !State_VGB(iP, i, j, k, iBlock) = &
!!$          !     gm1*(Energy_GBI(i, j, k, iBlock, iFluid) - 0.5*   &
!!$          !     sum(State_VGB(iRhoUx:iRhoUz,i, j, k, iBlock)**2)  &
!!$          !     /State_VGB(iRho,i, j, k, iBlock) )
!!$
!!$          AdjEnergy_GBI(i,j,k,iBlock,iFluid) = AdjEnergy_GBI(i,j,k,iBlock,iFluid) + &
!!$               gm1*Adjoint_VGB(iP,i,j,k,iBlock)
!!$          Adjoint_VGB(iRho,i,j,k,iBlock) = Adjoint_VGB(iRho,i,j,k,iBlock) + &
!!$               0.5*(sum(State_VGB(iRhoUx:iRhoUz, i, j, k, iBlock)**2)/&
!!$               State_VGB(iRho, i, j, k, iBlock)**2)*&
!!$               Adjoint_VGB(iP,i,j,k,iBlock)
!!$          do d=iRhoUx,iRhoUz
!!$             Adjoint_VGB(d,i,j,k,iBlock) = Adjoint_VGB(d,i,j,k,iBlock) - &
!!$                  State_VGB(d, i, j, k, iBlock)/State_VGB(iRho, i, j, k, iBlock)*&
!!$                  Adjoint_VGB(iP,i,j,k,iBlock)
!!$          end do
!!$          Adjoint_VGB(iP,i,j,k,iBlock) = 0.0
!!$          
!!$       end do; end do; end do
!!$
!!$       if(iFluid > 1 .or. .not. IsMhd) CYCLE
!!$
!!$       call stop_mpi(" in calc_pressure_adjoint: not yet implemented.")
!!$    end do
  
  end subroutine calc_pressure_adjoint
  !ADJOINT SPECIFIC END

  !===========================================================================

  subroutine calc_energy(iMin, iMax, jMin, jMax, kMin, kMax, iBlock, &
       iFluidMin, iFluidMax)

    ! Calculate total energy (excluding B0):
    !
    !   E = p/(gamma-1) + 0.5*rho*u^2 + 0.5*b1^2
    !

    integer, intent(in) :: iMin, iMax, jMin, jMax, kMin, kMax, iBlock
    integer, intent(in) :: iFluidMin, iFluidMax
    integer::i,j,k
    !--------------------------------------------------------------------------

    call limit_pressure(iMin, iMax, jMin, jMax, kMin, kMax, &
         iBlock, iFluidMin, iFluidMax)

    do iFluid = iFluidMin, iFluidMax
       call select_fluid
       ! Calculate thermal plus kinetic energy
       do k=kMin, kMax; do j=jMin, jMax; do i=iMin, iMax
          if (State_VGB(iRho,i,j,k,iBlock) <= 0.0)then
             Energy_GBI(i, j, k, iBlock, iFluid) = 0.0
          else
             Energy_GBI(i, j, k, iBlock, iFluid) = &
                  inv_gm1*State_VGB(iP,i,j,k,iBlock) &
                  +0.5*(sum(State_VGB(iRhoUx:iRhoUz, i, j, k, iBlock)**2)/&
                  State_VGB(iRho, i, j, k, iBlock))
          end if
       end do; end do; end do

       if(nIonFluid == 1 .and. iFluid == 1)then
          if(UseElectronPressure)then
             do k=kMin, kMax; do j=jMin, jMax; do i=iMin, iMax
                Energy_GBI(i,j,k,iBlock,iFluid) = &
                     Energy_GBI(i,j,k,iBlock,iFluid) &
                     + inv_gm1*State_VGB(Pe_,i,j,k,iBlock)
             end do; end do; end do
          end if
          if(UseWavePressure)then
             do k=kMin, kMax; do j=jMin, jMax; do i=iMin, iMax
                Energy_GBI(i,j,k,iBlock,iFluid) = &
                     Energy_GBI(i,j,k,iBlock,iFluid) &
                     + sum(State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock))
             end do; end do; end do
          end if
       end if
       
       if(iFluid > 1 .or. .not. IsMhd) CYCLE

       ! Add magnetic energy for ion fluid
       do k=kMin, kMax; do j=jMin, jMax; do i=iMin, iMax
          Energy_GBI(i, j, k, iBlock, iFluid) = &
               Energy_GBI(i, j, k, iBlock, iFluid) + &
               0.5*sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)**2)
       end do; end do; end do

    end do

  end subroutine calc_energy

  !ADJOINT SPECIFIC BEGIN
  !===========================================================================
  subroutine calc_energy_adjoint(iMin, iMax, jMin, jMax, kMin, kMax, iBlock, &
       iFluidMin, iFluidMax)

!!$    ! Calculate total energy (excluding B0):
!!$    !
!!$    !   E = p/(gamma-1) + 0.5*rho*u^2 + 0.5*b1^2
!!$    !
!!$
!!$    use ModAdjoint, ONLY: Adjoint_VGB, AdjEnergy_GBI
!!$
    integer, intent(in) :: iMin, iMax, jMin, jMax, kMin, kMax, iBlock
    integer, intent(in) :: iFluidMin, iFluidMax
!!$    integer::i,j,k,d
!!$    !--------------------------------------------------------------------------
!!$
!!$    do iFluid = iFluidMin, iFluidMax
!!$       call select_fluid
!!$       ! Calculate thermal plus kinetic energy
!!$       do k=kMin, kMax; do j=jMin, jMax; do i=iMin, iMax
!!$          if (State_VGB(iRho,i,j,k,iBlock) <= 0.0)then
!!$             Energy_GBI(i, j, k, iBlock, iFluid) = 0.0
!!$             AdjEnergy_GBI(i,j,k,iBlock, iFluid) = 0.0
!!$          else
!!$             !Energy_GBI(i, j, k, iBlock, iFluid) = &
!!$             !     inv_gm1*State_VGB(iP,i,j,k,iBlock) &
!!$             !     +0.5*(sum(State_VGB(iRhoUx:iRhoUz, i, j, k, iBlock)**2)/&
!!$             !     State_VGB(iRho, i, j, k, iBlock))
!!$             Adjoint_VGB(iP,i,j,k,iBlock) = Adjoint_VGB(iP,i,j,k,iBlock) + &
!!$                  inv_gm1*AdjEnergy_GBI(i,j,k,iBlock,iFluid)
!!$             Adjoint_VGB(iRho,i,j,k,iBlock) = Adjoint_VGB(iRho,i,j,k,iBlock) - &
!!$                  0.5*(sum(State_VGB(iRhoUx:iRhoUz, i, j, k, iBlock)**2)/&
!!$                  State_VGB(iRho, i, j, k, iBlock)**2)*&
!!$                  AdjEnergy_GBI(i,j,k,iBlock,iFluid)
!!$             do d=iRhoUx,iRhoUz
!!$                Adjoint_VGB(d,i,j,k,iBlock) = Adjoint_VGB(d,i,j,k,iBlock) + &
!!$                     State_VGB(d, i, j, k, iBlock)/State_VGB(iRho, i, j, k, iBlock)*&
!!$                     AdjEnergy_GBI(i,j,k,iBlock,iFluid)
!!$             end do
!!$             AdjEnergy_GBI(i,j,k,iBlock, iFluid) = 0.0
!!$          end if
!!$       end do; end do; end do
!!$       
!!$       if(iFluid > 1 .or. .not. IsMhd) CYCLE
!!$
!!$       call stop_mpi(" in calc_energy_adjoint: not yet implemented.")
!!$    end do

  end subroutine calc_energy_adjoint
  !ADJOINT SPECIFIC END


  !===========================================================================

  subroutine calc_pressure_cell(iBlock)
    integer, intent(in) :: iBlock
    call calc_pressure(1,nI,1,nJ,1,nK,iBlock,1,nFluid)
  end subroutine calc_pressure_cell
  
  !ADJOINT SPECIFIC BEGIN
  !===========================================================================
  subroutine calc_pressure_cell_adjoint(iBlock)
    integer, intent(in) :: iBlock
    call calc_pressure_adjoint(1,nI,1,nJ,1,nK,iBlock,1,nFluid)
  end subroutine calc_pressure_cell_adjoint
  !ADJOINT SPECIFIC END

  !===========================================================================

  subroutine calc_pressure_ghost(iBlock)
    integer, intent(in) :: iBlock
    call calc_pressure(-1,nI+2,-1,nJ+2,-1,nK+2,iBlock,1,nFluid)
  end subroutine calc_pressure_ghost

  !===========================================================================

  subroutine calc_pressure_point(i, j, k, iBlock)
    integer, intent(in) :: i, j, k, iBlock
    call calc_pressure(i,i,j,j,k,k,iBlock,1,nFluid)
  end subroutine calc_pressure_point

  !===========================================================================

  subroutine calc_pressure1_cell(iBlock)
    integer, intent(in) :: iBlock
    call calc_pressure(1,nI,1,nJ,1,nK,iBlock,1,1)
  end subroutine calc_pressure1_cell

  !===========================================================================

  subroutine calc_pressure1_ghost(iBlock)
    integer, intent(in) :: iBlock
    call calc_pressure(-1,nI+2,-1,nJ+2,-1,nK+2,iBlock,1,1)
  end subroutine calc_pressure1_ghost

  !===========================================================================

  subroutine calc_pressure1_point(i, j, k, iBlock)
    integer, intent(in) :: i, j, k, iBlock
    call calc_pressure(i,i,j,j,k,k,iBlock,1,1)
  end subroutine calc_pressure1_point

  !===========================================================================

  subroutine calc_energy_cell(iBlock)
    integer, intent(in) :: iBlock
    call calc_energy(1,nI,1,nJ,1,nK,iBlock,1,nFluid)
  end subroutine calc_energy_cell

  !ADJOINT SPECIFIC BEGIN
  !===========================================================================
  subroutine calc_energy_cell_adjoint(iBlock)
    integer, intent(in) :: iBlock
    call calc_energy_adjoint(1,nI,1,nJ,1,nK,iBlock,1,nFluid)
  end subroutine calc_energy_cell_adjoint
  !ADJOINT SPECIFIC END

  !===========================================================================

  subroutine calc_energy_ghost(iBlock)
    integer, intent(in) :: iBlock
    call calc_energy(-1,nI+2,-1,nJ+2,-1,nK+2,iBlock,1,nFluid)
  end subroutine calc_energy_ghost

  !ADJOINT SPECIFIC BEGIN
  !===========================================================================
  subroutine calc_energy_ghost_adjoint(iBlock)
    integer, intent(in) :: iBlock
    call calc_energy_adjoint(-1,nI+2,-1,nJ+2,-1,nK+2,iBlock,1,nFluid)
  end subroutine calc_energy_ghost_adjoint
  !ADJOINT SPECIFIC END

  !===========================================================================

  subroutine calc_energy_point(i, j, k, iBlock)
    integer, intent(in) :: i, j, k, iBlock
    call calc_energy(i,i,j,j,k,k,iBlock,1,nFluid)
  end subroutine calc_energy_point

  !===========================================================================

  subroutine calc_energy1_cell(iBlock)
    integer, intent(in) :: iBlock
    call calc_energy(1,nI,1,nJ,1,nK,iBlock,1,1)
  end subroutine calc_energy1_cell

  !===========================================================================

  subroutine calc_energy1_ghost(iBlock)
    integer, intent(in) :: iBlock
    call calc_energy(-1,nI+2,-1,nJ+2,-1,nK+2,iBlock,1,1)
  end subroutine calc_energy1_ghost

  !===========================================================================

  subroutine calc_energy1_point(i, j, k, iBlock)
    integer, intent(in) :: i, j, k, iBlock
    call calc_energy(i,i,j,j,k,k,iBlock,1,1)
  end subroutine calc_energy1_point

  !===========================================================================
  subroutine limit_pressure(iMin, iMax, jMin, jMax, kMin, kMax, &
       iBlock, iFluidMin, iFluidMax)

    ! Keep pressure(s) in State_VGB above pMin_I limit

    integer, intent(in) :: iMin, iMax, jMin, jMax, kMin, kMax, iBlock
    integer, intent(in) :: iFluidMin, iFluidMax

    integer:: i, j, k
    !------------------------------------------------------------------------
    do iFluid = iFluidMin, iFluidMax
       if(pMin_I(iFluid) < 0.0) CYCLE
       iP = iP_I(iFluid)
       do k=kMin, kMax; do j=jMin, jMax; do i=iMin, iMax
          State_VGB(iP, i, j, k, iBlock) = max(pMin_I(iFluid), &
               State_VGB(iP, i, j, k, iBlock))
       end do; end do; end do
    end do

  end subroutine limit_pressure

  !===========================================================================
  subroutine limit_old_pressure(iBlock)

    ! Keep pressure(s) in StateOld_VCB above pMin_I limit

    integer, intent(in) :: iBlock

    integer:: i, j, k
    !------------------------------------------------------------------------
    do iFluid = 1, nFluid
       if(pMin_I(iFluid) < 0.0) CYCLE
       iP = iP_I(iFluid)
       do k=1, nK; do j=1, nJ; do i=1, nI
          StateOld_VCB(iP, i, j, k, iBlock) = max(pMin_I(iFluid), &
               StateOld_VCB(iP, i, j, k, iBlock))
       end do; end do; end do
    end do

  end subroutine limit_old_pressure

end module ModEnergy

