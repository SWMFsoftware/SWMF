!^CFG COPYRIGHT UM
!^CFG FILE IMPLICIT
!==============================================================================
subroutine implicit_init

  ! Set number of implicit blocks and variables, 
  ! and conversion array between explicit and implicit block indices
  ! The implicit blocks are contiguous (all used) from 1 ... nImplBLK

  use ModMain
  use ModImplicit
  use ModAdvance, ONLY: iTypeAdvance_B, ImplBlock_
  implicit none

  logical :: IsInitialized=.false.
  integer :: iBLK, iBlockImpl
  !---------------------------------------------------------------------------

  nImplBLK=count(iTypeAdvance_B(1:nBlock) == ImplBlock_)

  ! Check for too many implicit blocks
  if(nImplBLK>MaxImplBLK)then
     write(*,*)'ERROR: Too many implicit blocks!'
     write(*,*)'MaxImplBLK < nImplBLK :',MaxImplBLK,nImplBLK
     call stop_mpi( &
          'Change number of processors, reduce number of implicit blocks,'// &
          ' or increase MaxImplBLK in ModSize.f90 !')
  end if

  ! Number of implicit variables
  nImpl = nImplBLK*nwIJK
  ! Create conversion array and find the test block
  implBLKtest=1
  iBlockImpl=0
  do iBLK=1,nBlock
     if (iTypeAdvance_B(iBLK) == ImplBlock_) then
        iBlockImpl = iBlockImpl + 1
        impl2iBLK(iBlockImpl)=iBLK
        if(iBLK==BLKtest)implBLKtest=iBlockImpl
     endif
  end do

  ! The index of the test variable in the linear array
  implVARtest=VARtest+nw*(Itest-1+nI*(Jtest-1+nJ*(Ktest-1+nK*(implBLKtest-1))))

  if(.not.IsInitialized)then
     residual = bigdouble
     IsInitialized=.true.
  end if

end subroutine implicit_init
!==============================================================================

subroutine explicit2implicit(imin,imax,jmin,jmax,kmin,kmax,Var_VGB)

  ! Convert data structure Var_VGB of the implicit code to the explicit code

  use ModMain
  use ModAdvance, ONLY : State_VGB, Energy_GBI, nVar
  use ModMultiFluid, ONLY: select_fluid, iFluid, nFluid, iP
  use ModImplicit
  use ModGrayDiffusion, ONLY: get_impl_gray_diff_state
  implicit none

  integer,intent(in) :: imin,imax,jmin,jmax,kmin,kmax
  real, intent(out)  :: Var_VGB(nw,imin:imax,jmin:jmax,kmin:kmax,MaxImplBLK)

  integer :: implBLK, iBLK, iVar
  logical :: DoTest, DoTestMe

  character(len=*), parameter:: NameSub = 'explicit2implicit'
  !---------------------------------------------------------------------------
  call set_oktest(NameSub,DoTest,DoTestMe)
  if(DoTestMe)write(*,*)'Starting explicit2implicit: ',&
       'imin,imax,jmin,jmax,kmin,kmax=',imin,imax,jmin,jmax,kmin,kmax

  if(DoTestMe)write(*,*)'E=',Energy_GBI(Itest,Jtest,Ktest,BLKtest,:)

  call timing_start('expl2impl')

  if(UseSemiImplicit)then
     select case(TypeSemiImplicit)
     case("radiation")
        call get_impl_gray_diff_state(Var_VGB)
     case default
        call stop_mpi(NameSub//': no get_impl_state implemented for' &
             //TypeSemiImplicit)
     end select
  else
     do implBLK=1,nImplBLK
        iBLK = impl2iBLK(implBLK)
        Var_VGB(:,:,:,:,implBLK) = &
             State_VGB(:,imin:imax,jmin:jmax,kmin:kmax,iBLK)
        do iFluid = 1, nFluid
           call select_fluid
           Var_VGB(iP,:,:,:,implBLK) = &
                Energy_GBI(imin:imax,jmin:jmax,kmin:kmax,iBLK,iFluid)
        end do
     end do
  end if

  call timing_stop('expl2impl')

  if(DoTestMe.and.nImplBLK>0)write(*,*)'Finished explicit2implicit: Var_VGB=',&
       Var_VGB(:,iTest,jTest,kTest,implBLKtest)

end subroutine explicit2implicit

!==============================================================================

subroutine impl2expl(Var_VC, iBLK)

  ! Convert the implicit block Var_VC to block iBLK of the explicit code

  use ModSize,     ONLY : nI, nJ, nK
  use ModAdvance,  ONLY : nVar, State_VGB, Energy_GBI
  use ModEnergy,   ONLY : calc_pressure_cell
  use ModMultiFluid, ONLY: iFluid, nFluid, iP_I, iP
  implicit none

  real, intent(in)    :: Var_VC(nVar,nI,nJ,nK)
  integer, intent(in) :: iBLK
  integer :: iVar
  !---------------------------------------------------------------------------

  call timing_start('impl2expl')

  State_VGB(1:nVar,1:nI,1:nJ,1:nK,iBLK) = Var_VC

  do iFluid = 1, nFluid
     iP = iP_I(iFluid)
     Energy_GBI(1:nI,1:nJ,1:nK,iBLK,iFluid) = Var_VC(iP,:,:,:)
  end do
  call calc_pressure_cell(iBLK)

  call timing_stop('impl2expl')

end subroutine impl2expl

!==============================================================================

subroutine implicit2explicit(Var_VCB)

  use ModMain, ONLY: nI,nJ,nK,MaxImplBLK
  use ModImplicit, ONLY: nw, nImplBLK, impl2iBLK, &
       UseSemiImplicit, TypeSemiImplicit
  use ModGrayDiffusion, ONLY: update_impl_gray_diff
  implicit none

  real :: Var_VCB(nw,nI,nJ,nK,MaxImplBLK)
  integer :: implBLK, iBLK

  character(len=*), parameter:: NameSub = 'implicit2explicit'
  !---------------------------------------------------------------------------

  do implBLK=1,nImplBLK
     iBLK=impl2iBLK(implBLK)
     if(UseSemiImplicit)then
        select case(TypeSemiImplicit)
        case('radiation')
           call update_impl_gray_diff(iBLK, Var_VCB(:,:,:,:,implBLK))
        case default
           call stop_mpi(NameSub//': no update_impl implemented for' &
                //TypeSemiImplicit)
        end select
     else
        call impl2expl(Var_VCB(:,:,:,:,implBLK),iBLK)
     end if
  end do

end subroutine implicit2explicit

!=============================================================================
subroutine get_residual(IsLowOrder, DoCalcTimestep, DoSubtract, Var_VCB, &
     Res_VCB)

  ! If IsLowOrder is true apply low  order scheme
  ! otherwise             apply high order scheme
  !
  ! If DoCalcTimestep is true calculate time step based on CFL condition
  !
  ! If DoSubtract is true return  Res_VCB = Var_VCB(t+dtexpl)-Var_VCB(t) 
  ! otherwise return              Res_VCB = Var_VCB(t+dtexpl)

  use ModMain
  use ModAdvance, ONLY : FluxType,time_BLK
  use ModGeometry, ONLY : true_cell
  use ModImplicit

  use ModMpi
  implicit none

  logical, intent(in) :: IsLowOrder, DoCalcTimestep, DoSubtract
  real, intent(in)    :: Var_VCB(nVar,nI,nJ,nK,MaxImplBLK)
  ! The actual Var_VCB and Res_VCB arguments may be the same array: 
  ! intent(inout)
  real, intent(inout) :: Res_VCB(nVar,nI,nJ,nK,MaxImplBLK)

  real    :: CflTmp
  integer :: nOrderTmp, nStageTmp, implBLK, iBLK
  character (len=10) :: FluxTypeTmp

  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------

  call set_oktest('get_residual',DoTest,DoTestMe)

  call timing_start('get_residual')

  if(DoTestMe.and.nImplBLK>0)&
       write(*,*)'get_residual DoSubtract,IsLowOrder,Var_VCB=',&
       DoSubtract,IsLowOrder,Var_VCB(Ktest,VARtest,Itest,Jtest,implBLKtest)

  nStageTmp       = nStage
  nStage          = 1
  if(IsLowOrder)then
     nOrderTmp    = nOrder
     nOrder       = nOrder_impl
     FluxTypeTmp  = FluxType
     FluxType     = FluxTypeImpl
  endif
  if(UseDtFixed)then
     do implBLK=1,nimplBLK
        iBLK=impl2iBLK(implBLK)
        time_BLK(:,:,:,iBLK)=0.0
        where(true_cell(1:nI,1:nJ,1:nK,iBLK)) &
             time_BLK(1:nI,1:nJ,1:nK,iBLK) = DtExpl
     end do
  else
     CflTmp = Cfl
     Cfl    = 0.5
  end if

  ! Res_VCB = Var_VCB(t+dt)
  call implicit2explicit(Var_VCB)
  call exchange_messages
  call advance_expl(DoCalcTimestep)
  call explicit2implicit(1,nI,1,nJ,1,nK,Res_VCB)

  if(DoSubtract) Res_VCB(:,:,:,:,1:nImplBLK) = Res_VCB(:,:,:,:,1:nImplBLK) &
       - Var_VCB(:,:,:,:,1:nImplBLK)

  if(DoTestMe.and.nImplBLK>0)write(*,*)'get_residual Res_VCB:',&
       Res_VCB(VARtest,Itest,Jtest,Ktest,implBLKtest)

  ! Restore global variables
  nStage      = nStageTmp
  if(IsLowOrder)then
     nOrder   = nOrderTmp
     FluxType = FluxTypeTmp 
  end if
  if (.not.UseDtFixed) Cfl = CflTmp

  call timing_stop('get_residual')

end subroutine get_residual
!==============================================================================
subroutine get_semi_impl_rhs(StateImpl_VGB, Rhs_VCB)

  use ModImplicit, ONLY: StateSemi_VGB, nw, nImplBlk, impl2iblk, &
       TypeSemiImplicit
  use ModMain, ONLY: dt
  use ModSize, ONLY: nI, nJ, nK, MaxImplBlk
  use ModGrayDiffusion, ONLY: get_gray_diffusion_rhs
  implicit none

  real, intent(in)  :: StateImpl_VGB(nw,0:nI+1,0:nJ+1,0:nK+1,MaxImplBlk)
  real, intent(out) :: Rhs_VCB(nw,nI,nJ,nK,MaxImplBlk)

  integer :: iImplBlock, iBlock, i, j, k, iVar

  character(len=*), parameter:: NameSub = 'get_semi_impl_rhs'
  !------------------------------------------------------------------------
  ! Fill in StateSemi so it can be message passed
  do iImplBlock = 1, nImplBLK
     iBlock = impl2iBLK(iImplBlock)
     do k = 1, nK; do j = 1, nJ; do i = 1, nI; do iVar = 1, nw
        StateSemi_VGB(iVar,i,j,k,iBlock) = StateImpl_VGB(iVar,i,j,k,iImplBlock)
     end do; end do; end do; end do
  end do
  !                       DoOneLayer DoFacesOnly No UseMonoteRestrict
  call message_pass_cells8(.true.,    .true.,     .false., nw, StateSemi_VGB)

  do iImplBlock = 1, nImplBLK
     iBlock = impl2iBLK(iImplBlock)
     select case(TypeSemiImplicit)
     case('radiation')
        Rhs_VCB(:,:,:,:,iImplBlock) = 0.0
        call get_gray_diffusion_rhs(iBlock, StateSemi_VGB(:,:,:,:,iBlock), &
             Rhs_VCB(:,:,:,:,iImplBlock), IsLinear=.false.)
     case default
        call stop_mpi(NameSub//': no get_rhs implemented for' &
             //TypeSemiImplicit)
     end select

     Rhs_VCB(:,:,:,:,iImplBlock) = dt*Rhs_VCB(:,:,:,:,iImplBlock)
  end do

end subroutine get_semi_impl_rhs
!==============================================================================
subroutine get_semi_impl_residual(StateImpl_VCB)

  use ModImplicit, ONLY: StateSemi_VGB, nw, nImplBlk, impl2iblk, &
       TypeSemiImplicit
  use ModMain, ONLY: dt
  use ModSize, ONLY: nI, nJ, nK, MaxImplBlk
  use ModGrayDiffusion, ONLY: get_gray_diffusion_rhs
  implicit none

  real, intent(inout) :: StateImpl_VCB(nw,nI,nJ,nK,MaxImplBlk)

  integer :: iImplBlock, iBlock, i, j, k, iVar
  real, allocatable, save :: Rhs_VC(:,:,:,:)

  character(len=*), parameter:: NameSub = 'get_semi_impl_residual'
  !------------------------------------------------------------------------

  if(.not.allocated(Rhs_VC)) allocate(Rhs_VC(nW,nI,nJ,nK))

  ! Fill in StateSemi so it can be message passed
  do iImplBlock = 1, nImplBLK
     iBlock = impl2iBLK(iImplBlock)
     do k = 1, nK; do j = 1, nJ; do i = 1, nI; do iVar = 1, nw
        StateSemi_VGB(iVar,i,j,k,iBlock) = StateImpl_VCB(iVar,i,j,k,iImplBlock)
     end do; end do; end do; end do
  end do
  !                       DoOneLayer DoFacesOnly No UseMonoteRestrict
  call message_pass_cells8(.true.,    .true.,     .false., nw, StateSemi_VGB)

  do iImplBlock = 1, nImplBLK
     iBlock = impl2iBLK(iImplBlock)

     select case(TypeSemiImplicit)
     case('radiation')
        call get_gray_diffusion_rhs(iBlock, &
             StateSemi_VGB(:,:,:,:,iBlock), Rhs_VC, IsLinear = .true.)
     case default
        call stop_mpi(NameSub//': no get_rhs implemented for' &
             //TypeSemiImplicit)
     end select

     StateImpl_VCB(:,:,:,:,iImplBlock) = dt*Rhs_VC

  end do

end subroutine get_semi_impl_residual
!==============================================================================
subroutine get_semi_impl_jacobian

  use ModImplicit, ONLY: nw, nImplBlk, impl2iblk, TypeSemiImplicit, &
       nStencil, MAT, ImplCoeff, wnrm
  use ModGrayDiffusion, ONLY: get_gray_diff_jacobian
  use ModMain, ONLY: nI, nJ, nK, nDim, Dt

  implicit none

  integer :: iImplBlock, iBlock, i, j, k, iStencil, iVar, jVar
  character(len=*), parameter:: NameSub = 'get_semi_impl_jacobian'
  !---------------------------------------------------------------------------
  do iImplBlock = 1, nImplBLK
     iBlock = impl2iBLK(iImplBlock)

     ! Get dR/dU
     select case(TypeSemiImplicit)
     case('radiation')
        call get_gray_diff_jacobian(iBlock, nw, MAT(:,:,:,:,:,:,iImplBlock))
     case default
        call stop_mpi(NameSub//': no get_rhs implemented for' &
             //TypeSemiImplicit)
     end select

     ! Form A = 1 - ImplCoeff*Dt*dR/dU matrix with normalization
     do iStencil = 1, nStencil; do k = 1, nK; do j = 1, nJ; do i = 1, nI
        do jVar = 1, nw; do iVar = 1, nw
           if(MAT(iVar, jVar, i, j, k, iStencil, iImplBlock) == 0.0) CYCLE 
           MAT(iVar, jVar, i, j, k, iStencil, iImplBlock) = &
                - MAT(iVar, jVar, i, j, k, iStencil, iImplBlock) &
                * dt * ImplCoeff * wnrm(jVar) / wnrm(iVar)
        end do; end do
     end do; end do; end do; end do
     do k = 1, nK; do j = 1, nJ; do i = 1, nI; do iVar = 1, nw
        MAT(iVar, iVar, i, j, k, 1, iImplBlock) = &
             MAT(iVar, iVar, i, j, k, 1, iImplBlock) + 1.0
     end do; end do; end do; end do
  end do

end subroutine get_semi_impl_jacobian
!==============================================================================
subroutine getsource(iBLK,Var_VCB,SourceImpl_VC)

  ! Get sources for block iBLK using implicit data Var_VCB

  use ModMain
  use ModVarIndexes
  use ModAdvance, ONLY : Source_VC  ! To communicate to calc_sources
  use ModImplicit, ONLY : nw,E_
  implicit none

  integer, intent(in) :: iBLK
  real, intent(in)    :: Var_VCB(nI,nJ,nK,nw)
  real, intent(out)   :: SourceImpl_VC(nw,nI,nJ,nK)

  logical :: qUseDivbSource
  integer::iVar, iFluid
  !--------------------------------------------------------------------------

  call timing_start('getsource')

  qUseDivbSource   =UseDivbSource
  UseDivbSource    =.false.
  
  call impl2expl(Var_VCB,iBLK)
  globalBLK = iBLK

  !!! Explicit time dependence  t+ImplCoeff*dt !!!
  !call calc_point_sources(t+ImplCoeff*dt)
  call calc_sources

  SourceImpl_VC = Source_VC(1:nVar,:,:,:)
  ! Overwrite pressure source terms with energy source term
  SourceImpl_VC(iP_I,:,:,:) = Source_VC(Energy_:Energy_+nFluid-1,:,:,:)

  UseDivbSource   =qUseDivbSource
  call timing_stop('getsource')

end subroutine getsource

!==============================================================================
subroutine get_face_flux(StateCons_VC,B0_DC,nI,nJ,nK,iDim,iBlock,Flux_VC)

  ! We need the cell centered physical flux function, but to keep
  ! the implicit scheme general for all equations, we reuse
  ! subroutine get_physical_flux from ModFaceFlux.

  use ModVarIndexes,ONLY: nFluid, nVar, Energy_
  use ModProcMH,   ONLY: iProc
  use ModMain,     ONLY: nDim, x_, y_, z_, &
       ProcTest, BlkTest,iTest,jTest,kTest
  use ModFaceFlux, ONLY: nFlux, iFace, jFace, kFace, Area, &
       set_block_values, set_cell_values, get_physical_flux, &
       HallJx, HallJy, HallJz, DoTestCell
  use ModHallResist, ONLY: UseHallResist, HallJ_CD
  use ModMultiFluid, ONLY: iFluid, nFluid, iP_I, iP

  implicit none

  integer, intent(in):: nI,nJ,nK,idim,iBlock
  real, intent(in)   :: StateCons_VC(nVar,nI,nJ,nK)
  real, intent(in)   :: B0_DC(nDim,nI,nJ,nK)
  real, intent(out)  :: Flux_VC(nVar,nI,nJ,nK)

  real :: Primitive_V(nVar), Conservative_V(nFlux), Flux_V(nFlux)

  real :: Un_I(nFluid+1), En
  integer :: i, j, k

  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------

  if(iBlock==BLKtest .and. iProc==PROCtest)then
     call set_oktest('get_face_flux', DoTest, DoTestMe)
  else
     DoTest=.false.; DoTestMe=.false.
  end if

  call set_block_values(iBlock, iDim)
  ! Set iFace=i, jFace=j, kFace=k so that 
  ! call set_cell_values and call get_physical_flux work
  ! This is not quite right but good enough for the preconditioner
  do k = 1, nK; kFace=k; do j = 1, nJ; jFace=j; do i = 1, nI; iFace=i

     DoTestCell = DoTestMe .and. &
          i==iTest .and. j==jTest .and. k==kTest

     Primitive_V = StateCons_VC( :,i, j, k)
     call conservative_to_primitive(Primitive_V)

     Conservative_V(1:nVar) = StateCons_VC( :,i, j, k)
     do iFluid=1, nFluid
        iP = iP_I(iFluid)
        Conservative_V(iP) = Primitive_V(iP)
        Conservative_V(nVar+iFluid) = StateCons_VC( iP,i, j, k)
     end do

     if(UseHallResist)then
        HallJx = HallJ_CD(i, j, k, x_)
        HallJy = HallJ_CD(i, j, k, y_)
        HallJz = HallJ_CD(i, j, k, z_)
     end if

     call set_cell_values
     call get_physical_flux(Primitive_V, &
          B0_DC(x_, i, j, k), &
          B0_DC(y_, i, j, k), &
          B0_DC(z_, i, j, k), &
          Conservative_V, Flux_V, Un_I, En)

     Flux_VC( 1:nVar,i, j, k)= Flux_V(1:nVar)*Area

     ! Replace pressure flux with energy flux
     Flux_VC(iP_I,i,j,k) = Flux_V(Energy_:Energy_+nFluid-1)*Area

  end do; end do; end do

end subroutine get_face_flux

!==============================================================================
subroutine get_cmax_face(Var_VF,B0_DF,qnI,qnJ,qnK,iDim,iBlock,Cmax)

  use ModProcMH,   ONLY: iProc
  use ModMain,     ONLY: nDim, x_, y_, z_, ProcTest, BlkTest,iTest,jTest,kTest
  use ModImplicit, ONLY: nw
  use ModFaceFlux, ONLY: DoTestCell, iFace, jFace, kFace, Area, &
       set_block_values, set_cell_values, get_speed_max, nFluid, &
       DoLf, DoAw, DoRoe, DoHll, DoHlld, UnLeft_I, UnRight_I
  use ModAdvance,  ONLY: eFluid_

  implicit none

  integer, intent(in):: qnI,qnJ,qnK,idim,iBlock
  real, intent(in)   :: Var_VF(nw,qnI,qnJ,qnK)
  real, intent(in)   :: B0_DF(ndim,qnI,qnJ,qnK)
  real, intent(out)  :: Cmax(qnI,qnJ,qnK)

  real :: Primitive_V(nw), Cmax_I(nFluid)

  character(len=*), parameter:: NameSub = 'get_cmax_face'
  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------

  if(iBlock==BLKtest .and. iProc==PROCtest)then
     call set_oktest('get_cmax_face', DoTest, DoTestMe)
  else
     DoTest=.false.; DoTestMe=.false.
  end if

  DoLf  = .true.
  DoAw  = .false.
  DoRoe = .false.
  DoHll = .false.
  DoHlld= .false.

  ! The electron speed is set to zero (I can't remember why)
  UnLeft_I(eFluid_)  = 0.0
  UnRight_I(eFluid_) = 0.0

  call set_block_values(iBlock, iDim)
  do kFace = 1, qnK; do jFace = 1, qnJ; do iFace = 1, qnI

     DoTestCell = DoTestMe .and. &
          iFace==iTest .and. jFace==jTest .and. kFace==kTest

     Primitive_V = Var_VF(:,iFace, jFace, kFace)

     call conservative_to_primitive(Primitive_V)

     call set_cell_values

     call get_speed_max(Primitive_V, &
          B0_DF( x_,iFace, jFace, kFace), &
          B0_DF( y_,iFace, jFace, kFace), &
          B0_DF( z_,iFace, jFace, kFace), &
          cmax_I = Cmax_I)

     cmax(iFace, jFace, kFace) = maxval(Cmax_I)*Area

  end do; end do; end do

  if(DoTestMe)write(*,*) NameSub,': Area, cmax=', &
       Area, cmax(iTest, jTest, kTest)

end subroutine get_cmax_face

!==============================================================================
subroutine conservative_to_primitive(State_V)

  use ModImplicit, ONLY: nw
  use ModVarIndexes, ONLY: Bx_, Bz_, IsMhd, nFluid
  use ModMultiFluid, ONLY: select_fluid,           &
       iFluid, iRho, iRhoUx, iUx, iRhoUz, iUz, iP
  use ModPhysics, ONLY: gm1
  implicit none
  real, intent(inout):: State_V(nw)
  real :: InvRho
  !---------------------------------------------------------------------------
  do iFluid = 1, nFluid
     call select_fluid
     InvRho = 1.0/State_V(iRho)
     State_V(iP) = gm1*(State_V(iP) - &
          0.5*sum(State_V(iRhoUx:iRhoUz)**2)*InvRho)
     if(iFluid == 1 .and. IsMhd) &
          State_V(iP) = State_V(iP) - 0.5*gm1*sum(State_V(Bx_:Bz_)**2)
     State_V(iUx:iUz) = InvRho*State_V(iRhoUx:iRhoUz)
  end do

end subroutine conservative_to_primitive

!==============================================================================
subroutine getdt_courant(qdt)

  use ModProcMH
  use ModMain
  use ModAdvance, ONLY : B0_DGB
  use ModGeometry, ONLY : dx_BLK,dy_BLK,dz_BLK,dxyz,true_cell,true_BLK,vInv_CB
  use ModImplicit
  use ModMpi
  implicit none

  real, intent(out) :: qdt

  real :: cmax(nI,nJ,nK), B0_DC(nDim,nI,nJ,nK), qdt_local
  integer :: idim, implBLK, iBLK, iError

  logical :: DoTest, DoTestMe
  !-------------------------------------------------------------------------
  call set_oktest('getdt_courant',DoTest,DoTestMe)

  ! First calculate max(cmax/dx) for each cell and dimension
  qdt_local=0.0
  do implBLK=1,nImplBLK; 
     iBLK=impl2iBLK(implBLK); 
     dxyz(x_)=dx_BLK(iBLK); dxyz(y_)=dy_BLK(iBLK); dxyz(z_)=dz_BLK(iBLK)
     if(UseB0)then
        B0_DC = B0_DGB(:,1:nI,1:nJ,1:nK,iBLK)
     else
        B0_DC = 0.0
     end if

     do idim=1,ndim

        call get_cmax_face(Impl_VGB(1:nw,1:nI,1:nJ,1:nK,implBLK),B0_DC,&
             nI, nJ, nK, iDim, iBlk, Cmax)

        if(.not.true_BLK(iBLK))then
           where(.not.true_cell(1:nI,1:nJ,1:nK,iBLK))cmax=0.0
        end if

        qdt_local=max(qdt_local,maxval(cmax*vInv_CB(:,:,:,iBlk)))

        if(DoTestMe)write(*,*)'getdt_courant idim,dx,cmax,1/qdt=',&
             idim,cmax(Itest,Jtest,Ktest),qdt_local
     end do
  end do

  ! Take global maximum
  call MPI_allreduce(qdt_local,qdt,1,MPI_REAL,MPI_MAX,iComm,iError)

  if(DoTestMe)write(*,*)'1/dt_local,1/dt=',qdt_local,qdt

  ! Take inverse, and reduce so it is OK for 3D calculation
  qdt=0.3/qdt

  if(DoTestMe)write(*,*)'getdt_courant final dt=',qdt

end subroutine getdt_courant

