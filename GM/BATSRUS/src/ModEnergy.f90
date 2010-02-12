module ModEnergy

  use ModProcMH,  ONLY: iProc
  use ModMultiFluid
  use ModSize,    ONLY: nI, nJ, nK, gcn, MaxBlock
  use ModAdvance, ONLY: State_VGB, Energy_GBI, StateOld_VCB, EnergyOld_CBI,&
       UseNonConservative, nConservCrit, IsConserv_CB, UseAnisoPressure
  use ModPhysics, ONLY: Gm1, Inv_Gm1

  implicit none

  ! Variables used for converting energy to anisotropic pressure 
  real :: Ppar, Pperp, pOld, pNew 

contains

  !============================================================================
  subroutine calc_energy_or_pressure(iBlock)

    ! Calculate pressure from energy or energy from pressure depending on 
    ! the value of UseNonConservative and IsConserv_CB

    integer, intent(in) :: iBlock
    integer::i,j,k
    !--------------------------------------------------------------------------
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

    ! A mix of conservative and non-conservative scheme for aniso pressure
    ! Only work for single fluid now, nfluid = 1
    if(UseAnisoPressure)then
       do k=1, nK; do j=1,nJ; do i=1,nI
          if(IsConserv_CB(i,j,k,iBlock))then
             Pperp = State_VGB(Pperp_,i,j,k,iBlock)
             Ppar  = State_VGB(Ppar_,i,j,k,iBlock)
             pOld  = 2/3.*Pperp + 1/3.*Ppar
             pNew  = gm1*(Energy_GBI(i,j,k,iBlock,1) - 0.5* &
                  sum(State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock)**2)   &
                  /State_VGB(iRho,i,j,k,iBlock) )                 &
                  - gm1*0.5*sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)**2)
             
             State_VGB(Ppar_,i,j,k,iBlock)  = pNew/pOld*Ppar
             State_VGB(Pperp_,i,j,k,iBlock) = pNew/pOld*Pperp
          else
             Energy_GBI(i,j,k,iBlock,1) = &
                  inv_gm1*(2/3.*State_VGB(iP,i,j,k,iBlock) &
                  + 1/3.*State_VGB(iP-1,i,j,k,iBlock)) &
                  +0.5*(sum(State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock)**2)/&
                  State_VGB(iRho,i,j,k,iBlock))
          end if
       end do; end do; end do;
    end if

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

  end subroutine calc_energy_or_pressure

  !============================================================================

  subroutine calc_old_pressure(iBlock)

    ! Calculate pressure from energy for the old state
    !
    !   P = (gamma-1)*(E - 0.5*rho*u^2 - 0.5*b1^2)
    !
    ! For anisotropic pressure
    !   pNew = 2/3*(E - 0.5*rho*u^2 - 0.5*b1^2)
    !   pperp = 3*pNew/(2*pperp + ppar)*pperp = pNew/pOld*pperp
    !   ppar  = 3*pNew/(2*pperp + ppar)*ppar  = pNew/pOld*ppar

    integer, intent(in) :: iBlock
    integer :: i,j,k
    !--------------------------------------------------------------------------
    if(UseAnisoPressure)then
       ! Only works for single fluid now
       do k=1, nK; do j=1, nJ; do i=1, nI
          Pperp = StateOld_VCB(Pperp_,i,j,k,iBlock)
          Ppar  = StateOld_VCB(Ppar_,i,j,k,iBlock)
          pOld  = 2/3.*Pperp + 1/3.*Ppar
          pNew  = gm1*(EnergyOld_CBI(i,j,k,iBlock,1) - 0.5* &
               sum(StateOld_VCB(iRhoUx:iRhoUz,i,j,k,iBlock)**2)   &
               /StateOld_VCB(iRho,i,j,k,iBlock) )                 &
               - gm1*0.5*sum(StateOld_VCB(Bx_:Bz_,i,j,k,iBlock)**2)

          StateOld_VCB(Ppar_,i,j,k,iBlock)  = pNew/pOld*Ppar
          StateOld_VCB(Pperp_,i,j,k,iBlock) = pNew/pOld*Pperp
       end do; end do; end do
       RETURN
    end if

    do iFluid = 1, nFluid
       call select_fluid
       do k=1, nK; do j=1, nJ; do i=1, nI
          StateOld_VCB(iP, i, j, k,iBlock) = &
               gm1*(EnergyOld_CBI(i, j, k, iBlock, iFluid) - 0.5*   &
               sum(StateOld_VCB(iRhoUx:iRhoUz,i, j, k, iBlock)**2)  &
               /StateOld_VCB(iRho,i, j, k, iBlock) )
       end do; end do; end do

       if(iFluid > 1 .or. .not. IsMhd) CYCLE
       
       do k=1, nK; do j=1, nJ; do i=1, nI
          StateOld_VCB(iP, i, j, k,iBlock) = StateOld_VCB(iP, i, j, k,iBlock) &
               - gm1*0.5*sum(StateOld_VCB(Bx_:Bz_,i,j,k,iBlock)**2)
       end do; end do; end do
    end do
  end subroutine calc_old_pressure

  !============================================================================

  subroutine calc_pressure(iMin, iMax, jMin, jMax, kMin, kMax, iBlock, &
       iFluidMin, iFluidMax)

    ! Calculate pressure from energy
    !
    !   P = (gamma-1)*(E - 0.5*rho*u^2 - 0.5*b1^2)
    !
    ! For anisotropic pressure, correction as
    !   pNew = 2/3*(E - 0.5*rho*u^2 - 0.5*b1^2)
    !   pperp = 3*pNew/(2*pperp + ppar)*pperp = pNew/pOld*pperp
    !   ppar  = 3*pOld/(2*pperp + ppar)*ppar  = pNew/pOld*ppar

    integer, intent(in) :: iMin, iMax, jMin, jMax, kMin, kMax, iBlock
    integer, intent(in) :: iFluidMin, iFluidMax
    integer :: i,j,k    
    !--------------------------------------------------------------------------
    if(UseAnisoPressure)then
       ! Only works for single fluid now
       do k=1, nK; do j=1, nJ; do i=1, nI
          Pperp = State_VGB(Pperp_,i,j,k,iBlock)
          Ppar  = State_VGB(Ppar_,i,j,k,iBlock)
          pOld  = 2/3.*Pperp + 1/3.*Ppar
          pNew  = gm1*(Energy_GBI(i,j,k,iBlock,1) - 0.5* &
               sum(State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock)**2)   &
               /State_VGB(iRho,i,j,k,iBlock)                  &
               - 0.5*sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)**2))

          State_VGB(Ppar_,i,j,k,iBlock)  = pNew/pOld*Ppar
          State_VGB(Pperp_,i,j,k,iBlock) = pNew/pOld*Pperp
       end do; end do; end do
       RETURN
    end if

    do iFluid = iFluidMin, iFluidMax
       call select_fluid
       do k=kMin, kMax; do j=jMin, jMax; do i=iMin, iMax
          State_VGB(iP, i, j, k, iBlock) = &
               gm1*(Energy_GBI(i, j, k, iBlock, iFluid) - 0.5*   &
               sum(State_VGB(iRhoUx:iRhoUz,i, j, k, iBlock)**2)  &
               /State_VGB(iRho,i, j, k, iBlock) )
       end do; end do; end do

       if(iFluid > 1 .or. .not. IsMhd) CYCLE
       do k=kMin, kMax; do j=jMin, jMax; do i=iMin, iMax
          State_VGB(iP, i, j, k,iBlock) = State_VGB(iP, i, j, k,iBlock) &
               - gm1*0.5*sum(State_VGB(Bx_:Bz_,i, j, k,iBlock)**2)
       end do; end do; end do          
    end do
    
  end subroutine calc_pressure

  !===========================================================================

  subroutine calc_energy(iMin, iMax, jMin, jMax, kMin, kMax, iBlock, &
       iFluidMin, iFluidMax)

    ! Calculate total energy (excluding B0):
    !
    !   E = p/(gamma-1) + 0.5*rho*u^2 + 0.5*b1^2
    !
    ! For anisotropic pressure:
    !
    !   E = (2/3*pperp + 1/3*ppar)/(5/3-1) + 0.5*rho*u^2 + 0.5*b1^2
    

    integer, intent(in) :: iMin, iMax, jMin, jMax, kMin, kMax, iBlock
    integer, intent(in) :: iFluidMin, iFluidMax
    integer::i,j,k
    !--------------------------------------------------------------------------

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
             if(UseAnisoPressure) &
                  Energy_GBI(i, j, k, iBlock, iFluid) = &
                  Energy_GBI(i, j, k, iBlock, iFluid) &
                  +0.5*(State_VGB(iP-1,i,j,k,iBlock)-State_VGB(iP,i,j,k,iBlock))
          end if
       end do; end do; end do
       
       if(iFluid > 1 .or. .not. IsMhd) CYCLE

       ! Add magnetic energy for ion fluid
       do k=kMin, kMax; do j=jMin, jMax; do i=iMin, iMax
          Energy_GBI(i, j, k, iBlock, iFluid) = &
               Energy_GBI(i, j, k, iBlock, iFluid) + &
               0.5*sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)**2)
       end do; end do; end do

    end do

  end subroutine calc_energy

  !===========================================================================

  subroutine calc_pressure_cell(iBlock)
    integer, intent(in) :: iBlock
    call calc_pressure(1,nI,1,nJ,1,nK,iBlock,1,nFluid)
  end subroutine calc_pressure_cell

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

  !===========================================================================

  subroutine calc_energy_ghost(iBlock)
    integer, intent(in) :: iBlock
    call calc_energy(-1,nI+2,-1,nJ+2,-1,nK+2,iBlock,1,nFluid)
  end subroutine calc_energy_ghost

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

end module ModEnergy

