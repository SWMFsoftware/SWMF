!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModEnergy

  use ModProcMH,     ONLY: iProc
  use ModMain,       ONLY: BlkTest,iTest,jTest,kTest,ProcTest
  use ModMultiFluid, ONLY: nFluid, iFluid, IonLast_, &
       iRho, iRhoUx, iRhoUy, iRhoUz, iP, iP_I, DoConserveNeutrals, &
       select_fluid
  use ModSize,       ONLY: nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK
  use ModAdvance,    ONLY: State_VGB, Bx_, By_, Bz_, IsMhd, &
       Energy_GBI, StateOld_VCB, EnergyOld_CBI,&
       UseNonConservative, nConservCrit, IsConserv_CB
  use ModPhysics,    ONLY: Gm1, Inv_Gm1, pMin_I

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
             Energy_GBI(i,j,k,iBlock,iFluid) = 0.0
          else
             Energy_GBI(i, j, k, iBlock, iFluid) = &
                  inv_gm1*State_VGB(iP,i,j,k,iBlock) &
                  +0.5*(sum(State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock)**2)/&
                  State_VGB(iRho,i,j,k,iBlock))
          end if
       end do; end do; end do

       if(iFluid > 1 .or. .not. IsMhd) CYCLE

       ! Add magnetic energy for ion fluid
       do k=kMin, kMax; do j=jMin, jMax; do i=iMin, iMax
          Energy_GBI(i,j,k,iBlock, iFluid) = &
               Energy_GBI(i,j,k,iBlock, iFluid) + &
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
    call calc_pressure(MinI,MaxI,MinJ,MaxJ,MinK,MaxK,iBlock,1,nFluid)
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
    call calc_pressure(MinI,MaxI,MinJ,MaxJ,MinK,MaxK,iBlock,1,1)
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

  subroutine calc_energy_ghost(iBlock, DoResChangeOnlyIn)

    use BATL_lib, ONLY: DiLevelNei_IIIB

    integer, intent(in) :: iBlock
    logical, optional, intent(in) :: DoResChangeOnlyIn

    integer :: i,j,k
    logical :: DoResChangeOnly
    !--------------------------------------------------------------------------

    DoResChangeOnly =.false.
    if(present(DoResChangeOnlyIn)) DoResChangeOnly = DoResChangeOnlyIn

    if( DoResChangeOnly ) then
       if( .not.any(abs(DiLevelNei_IIIB(-1:1,-1:1,-1:1,iBlock)) == 1) ) RETURN
    end if

!!$    !------------------- for calculation on ghost cells only ----------------
!!$    if( DoResChangeOnly ) then
!!$       
!!$       if( any(abs(DiLevelNei_IIIB(-1,:,:,iBlock)) == 1) ) then
!!$          !call limit_pressure(-1,0,MinJ,MaxJ,MinK,MaxK,iBlock,1,nFluid)
!!$          call calc_energy(-1,0,MinJ,MaxJ,MinK,MaxK,iBlock,1,nFluid)
!!$       end if
!!$
!!$       if( any(abs(DiLevelNei_IIIB(1,:,:,iBlock)) == 1) ) then
!!$          !call limit_pressure(nI+1,nI+2,MinJ,MaxJ,MinK,MaxK,iBlock,1,nFluid)
!!$          call calc_energy(nI+1,nI+2,MinJ,MaxJ,MinK,MaxK,iBlock,1,nFluid)
!!$       end if
!!$
!!$       if( any(abs(DiLevelNei_IIIB(:,-1,:,iBlock)) == 1) ) then
!!$          !call limit_pressure(MinI,MaxI,-1,0,MinK,MaxK,iBlock,1,nFluid)
!!$          call calc_energy(MinI,MaxI,-1,0,MinK,MaxK,iBlock,1,nFluid)
!!$       end if
!!$
!!$       if( any(abs(DiLevelNei_IIIB(:,1,:,iBlock)) == 1) ) then
!!$          !call limit_pressure(MinI,MaxI,nJ+1,nJ+2,MinK,MaxK,iBlock,1,nFluid)
!!$          call calc_energy(MinI,MaxI,nJ+1,nJ+2,MinK,MaxK,iBlock,1,nFluid)
!!$       end if
!!$       if( any(abs(DiLevelNei_IIIB(:,:,-1,iBlock)) == 1) ) then
!!$          !call limit_pressure(MinI,MaxI,MinJ,MaxJ,-1,0,iBlock,1,nFluid)
!!$          call calc_energy(MinI,MaxI,MinJ,MaxJ,-1,0,iBlock,1,nFluid)
!!$       end if
!!$       if( any(abs(DiLevelNei_IIIB(:,:,1,iBlock)) == 1) ) then
!!$          !call limit_pressure(MinI,MaxI,MinJ,MaxJ,nK+1,nK+2,iBlock,1,nFluid)
!!$          call calc_energy(MinI,MaxI,MinJ,MaxJ,nK+1,nK+2,iBlock,1,nFluid)
!!$       end if
!!$
!!$
!!$    end if

    call limit_pressure(MinI,MaxI,MinJ,MaxJ,MinK,MaxK,iBlock,1,nFluid)

    do iFluid = 1,nFluid

       call select_fluid

       if(IsMhd .and. iFluid == 1) then
          ! MHD energy
          where(State_VGB(iRho, MinI:MaxI, MinJ:MaxJ, MinK:MaxK, iBlock) <= 0.0)
             Energy_GBI(MinI:MaxI, MinJ:MaxJ, MinK:MaxK, iBlock, iFluid) = 0.0
          elsewhere
             Energy_GBI(MinI:MaxI, MinJ:MaxJ, MinK:MaxK, iBlock, iFluid) = &                  
                  inv_gm1*State_VGB(iP,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,iBlock) &                
                  +0.5*((State_VGB(iRhoUx, MinI:MaxI, MinJ:MaxJ, MinK:MaxK, iBlock)**2 + &    
                  State_VGB(iRhoUy, MinI:MaxI, MinJ:MaxJ, MinK:MaxK, iBlock)**2 + &           
                  State_VGB(iRhoUz, MinI:MaxI, MinJ:MaxJ, MinK:MaxK, iBlock)**2)/&            
                  State_VGB(iRho, MinI:MaxI, MinJ:MaxJ, MinK:MaxK, iBlock)) + &               
                  0.5*(State_VGB(Bx_,MinI:MaxI, MinJ:MaxJ, MinK:MaxK, iBlock)**2 +&           
                  State_VGB(By_,MinI:MaxI, MinJ:MaxJ, MinK:MaxK, iBlock)**2 +&                
                  State_VGB(Bz_,MinI:MaxI, MinJ:MaxJ, MinK:MaxK, iBlock)**2 )                 
          end where
       else
          ! HD energy
          where(State_VGB(iRho, MinI:MaxI, MinJ:MaxJ, MinK:MaxK, iBlock) <= 0.0)
             Energy_GBI(MinI:MaxI, MinJ:MaxJ, MinK:MaxK, iBlock, iFluid) = 0.0
          elsewhere
             Energy_GBI(MinI:MaxI, MinJ:MaxJ, MinK:MaxK, iBlock, iFluid) = &                  
                  inv_gm1*State_VGB(iP,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,iBlock) &                
                  +0.5*((State_VGB(iRhoUx, MinI:MaxI, MinJ:MaxJ, MinK:MaxK, iBlock)**2 + &    
                  State_VGB(iRhoUy, MinI:MaxI, MinJ:MaxJ, MinK:MaxK, iBlock)**2 + &           
                  State_VGB(iRhoUz, MinI:MaxI, MinJ:MaxJ, MinK:MaxK, iBlock)**2)/&            
                  State_VGB(iRho, MinI:MaxI, MinJ:MaxJ, MinK:MaxK, iBlock))               
          end where
       end if

    end do

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
    call calc_energy(MinI,MaxI,MinJ,MaxJ,MinK,MaxK,iBlock,1,1)
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

  !============================================================================
  ! moved form file exchange_messages.f90 
  subroutine correctP(iBlock)

    ! Make pressure and energy consistent and maintain thermal energy ratio 
    ! at a reasonable value (this also excludes negative pressure)

    use ModProcMH
    use ModMain,       ONLY: nI,nJ,nK,Itest,Jtest,Ktest,BLKtest
    use ModVarIndexes, ONLY: rho_, rhoUx_, rhoUy_, rhoUz_, Bx_, By_, Bz_, P_
    use ModAdvance,    ONLY: State_VGB, Energy_GBI
    use ModPhysics,    ONLY: gm1, inv_gm1, Pratio_hi, Pratio_lo
    use ModGeometry,   ONLY: true_cell
    use BATL_lib,      ONLY: Xyz_DGB
    implicit none

    integer, intent(in) :: iBlock

    integer :: i,j,k
    real :: inv_dratio, qp, qe, qth, qratio, qd, qde, qpmin, &
         qdesum, qdesumabs, qderelmax

    real, dimension(1:nI,1:nJ,1:nK) :: P_old

    integer :: ierror1=-1, ierror2=-1, ierror3=-1, loc(3)

    logical :: oktest, oktest_me
    !--------------------------------------------------------------------------

    if(iBlock==BLKtest)then
       call set_oktest('correctP',oktest,oktest_me)
    else
       oktest=.false.; oktest_me=.false.
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

       if(oktest_me.and.i==Itest.and.J==Jtest.and.K==Ktest)&
            write(*,*)'CorrectP at me,BLK,i,j,k=',&
            iProc,BLKtest,Itest,Jtest,Ktest, &
            ', initial P,E=',qp,qe

       ! Memorize smallest pressure
       qpmin=min(qp,qpmin)

       ! Thermal energy
       qth=inv_gm1*qp

       ! Deviation=extra total energy=qe-inv_gm1*qp-(rhoU**2/rho+B**2)/2
       qd=qE-qth                                                         &
            -0.5*(State_VGB(rhoUx_,i,j,k,iBlock)**2+                         &
            State_VGB(rhoUy_,i,j,k,iBlock)**2+                               &
            State_VGB(rhoUz_,i,j,k,iBlock)**2)/State_VGB(rho_,i,j,k,iBlock)      &
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
       State_VGB(P_,i,j,k,iBlock)=gm1*(qth+qd-qde)

       ! We should now have E=inv_gm1*P+(rhoU**2/rho+B**2)/2:
       !
       ! qp*inv_gm1+qd-qde + (rhoU**2/rho+B**2)/2 = qe-qde = E
       !
       ! Correct!

       if(oktest_me.and.i==Itest.and.J==Jtest.and.K==Ktest)then
          write(*,*)'qp,qth,qe,qd,qratio,qde=',qp,qth,qe,qd,qratio,qde
          write(*,*)'CorrectP, final P=',State_VGB(P_,i,j,k,iBlock)
       end if

    end do; end do; end do

    if(qpmin<0.)then
       if(ierror1==-1)then
          loc=minloc(P_old)
          write(*,*)'Negative P at me,iBLK,I,J,K,x,y,z,val',&
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

    if(oktest_me)write(*,*)'CorrectP qpmin=',qpmin

  end subroutine correctP


end module ModEnergy

