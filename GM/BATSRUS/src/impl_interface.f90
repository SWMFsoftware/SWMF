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

subroutine explicit2implicit(imin,imax,jmin,jmax,kmin,kmax,w)

  ! Convert data structure w of the implicit code to the explicit code

  use ModMain
  use ModAdvance, ONLY : State_VGB, E_BLK,nVar
  use ModImplicit
  implicit none

  integer,intent(in) :: imin,imax,jmin,jmax,kmin,kmax
  real, intent(out)  :: w(imin:imax,jmin:jmax,kmin:kmax,nw,MaxImplBLK)
  integer :: implBLK, iBLK, iVar
  logical :: oktest, oktest_me
  !---------------------------------------------------------------------------
  call set_oktest('explicit2implicit',oktest,oktest_me)
  if(oktest_me)write(*,*)'Starting explicit2implicit: ',&
       'imin,imax,jmin,jmax,kmin,kmax=',imin,imax,jmin,jmax,kmin,kmax

  if(oktest_me)write(*,*)'E=',E_BLK(Itest,Jtest,Ktest,BLKtest)

  call timing_start('expl2impl')

  do implBLK=1,nImplBLK
     iBLK = impl2iBLK(implBLK)
     do iVar=1, nVar
        if(iVar==E_)then
           w(:,:,:,E_,implBLK) = E_BLK(imin:imax,jmin:jmax,kmin:kmax,iBLK)
        else
           w(:,:,:,iVar,implBLK) = &
                State_VGB(iVar,imin:imax,jmin:jmax,kmin:kmax,iBLK)
        end if
     end do
  end do

  call timing_stop('expl2impl')

  if(oktest_me.and.nImplBLK>0)write(*,*)'Finished explicit2implicit: w=',&
       w(Itest,Jtest,Ktest,VARtest,implBLKtest)

end subroutine explicit2implicit

!==============================================================================

subroutine impl2expl(w,iBLK)

  ! Convert the implicit block w to block iBLK of the explicit code

  use ModMain
  use ModVarIndexes
  use ModAdvance, ONLY : &
       State_VGB, E_BLK
  use ModImplicit, ONLY : nw,E_
  use ModPhysics, ONLY : gm1
  implicit none

  real, intent(in)    :: w(nI,nJ,nK,nw)
  integer, intent(in) :: iBLK
  integer :: iVar
  !---------------------------------------------------------------------------

  call timing_start('impl2expl')

  do iVar=1,nVar
     if(iVar==E_)then
        E_BLK(1:nI,1:nJ,1:nK,iBLK) = w(:,:,:,E_)
     else
        State_VGB(iVar,1:nI,1:nJ,1:nK,iBLK) = w(:,:,:,iVar)
     end if
  end do

  State_VGB(P_, 1:nI,1:nJ,1:nK,iBLK) = gm1*(E_BLK(1:nI,1:nJ,1:nK,iBLK)-0.5*( &
       (State_VGB(rhoUx_,1:nI,1:nJ,1:nK,iBLK)**2                             &
       +State_VGB(rhoUy_,1:nI,1:nJ,1:nK,iBLK)**2                             &
       +State_VGB(rhoUz_,1:nI,1:nJ,1:nK,iBLK)**2                             &
       )/State_VGB(rho_,1:nI,1:nJ,1:nK,iBLK)                                 &
       +State_VGB(Bx_,1:nI,1:nJ,1:nK,iBLK)**2                                &
       +State_VGB(By_,1:nI,1:nJ,1:nK,iBLK)**2                                &
       +State_VGB(Bz_,1:nI,1:nJ,1:nK,iBLK)**2)                             )

  call timing_stop('impl2expl')

end subroutine impl2expl

!==============================================================================

subroutine implicit2explicit(w)

  use ModMain, ONLY : nI,nJ,nK,MaxImplBLK
  use ModImplicit, ONLY : nw,nImplBLK,impl2iBLK
  implicit none

  real :: w(nI,nJ,nK,nw,MaxImplBLK)
  integer :: implBLK, iBLK
  !---------------------------------------------------------------------------

  do implBLK=1,nImplBLK
     iBLK=impl2iBLK(implBLK)
     call impl2expl(w(:,:,:,:,implBLK),iBLK)
  end do

end subroutine implicit2explicit

!=============================================================================
subroutine get_residual(low_order, do_calc_timestep, do_subtract, w, RES)

  ! If low_order is true apply low  order scheme
  ! otherwise            apply high order scheme
  !
  ! If do_calc_timestep is true calculate time step based on CFL condition
  !
  ! If do_subtract is true return RES = w(t+dtexpl)-w(t) 
  ! otherwise return              RES = w(t+dtexpl)

  use ModMain
  use ModAdvance, ONLY : FluxType,time_BLK
  use ModGeometry, ONLY : true_cell
  use ModImplicit

  use ModMpi
  implicit none

  logical, intent(in) :: low_order, do_calc_timestep, do_subtract
  real, intent(in)    :: w(nI,nJ,nK,nw,MaxImplBLK)
  real, intent(out)   :: RES(nI,nJ,nK,nw,MaxImplBLK)

  real    :: qcfl
  integer :: qnorder, nStageTmp, implBLK, iBLK
  character (len=10) :: FluxTypeTmp

  logical:: sourceunsplit0

  real*8  :: time_before
  logical :: oktest, oktest_me
  !--------------------------------------------------------------------------

  call set_oktest('get_residual',oktest,oktest_me)

  call timing_start('get_residual')

  if(oktest_me.and.nImplBLK>0)&
       write(*,*)'get_residual do_subtract,low_order,w=',&
       do_subtract,low_order,w(Itest,Jtest,Ktest,VARtest,implBLKtest)

  nStageTmp    =nStage
  nStage       =1
  if(low_order)then
     qnorder      =norder
     norder       =norder_impl
     FluxTypeTmp  =FluxType
     FluxType     =FluxTypeImpl

     ! If the source is not to be treated implicitly pretend split sources
     if(.not.implsource)then
        sourceunsplit0=sourceunsplit
        sourceunsplit=.false.
     endif
  endif
  if(UseDtFixed)then
     do implBLK=1,nimplBLK
        iBLK=impl2iBLK(implBLK)
        time_BLK(:,:,:,iBLK)=0.0
        where(true_cell(1:nI,1:nJ,1:nK,iBLK)) &
             time_BLK(1:nI,1:nJ,1:nK,iBLK)=dtexpl
     end do
  else
     qcfl=cfl
     cfl=0.5
  end if

  ! RES = w(t+dt)
  call implicit2explicit(w)
  call exchange_messages
  call advance_expl(do_calc_timestep)
  call explicit2implicit(1,nI,1,nJ,1,nK,RES)

  if(do_subtract) &
       RES(1:nI,1:nJ,1:nK,1:nw,1:nImplBLK) = &
       RES(1:nI,1:nJ,1:nK,1:nw,1:nImplBLK) &
       - w(1:nI,1:nJ,1:nK,1:nw,1:nImplBLK)

  if(oktest_me.and.nImplBLK>0)write(*,*)'get_residual RES:',&
       RES(Itest,Jtest,Ktest,VARtest,implBLKtest)

  ! Restore global variables
  nStage=nStageTmp
  if(low_order)then
     norder=qnorder
     FluxType=FluxTypeTmp 
     if(.not.implsource)sourceunsplit=sourceunsplit0
  end if
  if(.not.UseDtFixed)cfl=qcfl

  call timing_stop('get_residual')

end subroutine get_residual

!==============================================================================
subroutine getsource(iBLK,w,s)

  ! Get source s for block iBLK using implicit data w

  use ModMain
  use ModVarIndexes
  use ModAdvance, ONLY : Source_VC  ! To communicate to calc_sources
  use ModImplicit, ONLY : nw,E_
  implicit none

  integer, intent(in) :: iBLK
  real, intent(in)    :: w(nI,nJ,nK,nw)
  real, intent(out)   :: s(nI,nJ,nK,nw)

  logical :: qUseDivbSource
    integer::iVar
  !--------------------------------------------------------------------------

  call timing_start('getsource')

  qUseDivbSource   =UseDivbSource
  UseDivbSource    =.false.
  

  call impl2expl(w,iBLK)
  globalBLK = iBLK

  !!! Explicit time dependence  t+ImplCoeff*dt !!!
  !call calc_point_sources(t+ImplCoeff*dt)
  call calc_sources

  do iVar=1,nVar-1
     s(:,:,:,iVar)  =Source_VC(iVar,1:nI,1:nJ,1:nK)
  end do

  s(:,:,:,E_)    =Source_VC(Energy_,1:nI,1:nJ,1:nK)

  UseDivbSource   =qUseDivbSource
  call timing_stop('getsource')

end subroutine getsource

!==============================================================================
subroutine get_face_flux(StateCons_CV,B0_CD,nI,nJ,nK,iDim,iBlock,Flux_CV)

  use ModProcMH,   ONLY: iProc
  use ModMain,     ONLY: nDim, x_, y_, z_, &
       ProcTest, BlkTest,iTest,jTest,kTest
  use ModImplicit, ONLY: nw
  use ModFaceFlux, ONLY: iBlockFace, iFace, jFace, kFace, &
       set_block_values, &
       set_cell_values_x, set_cell_values_y, set_cell_values_z, &
       get_physical_flux, &
       HallJx, HallJy, HallJz, DoTestCell, Area
  use ModHallResist, ONLY: UseHallResist, HallJ_CD
  use ModImplicit, ONLY: p_, e_

  implicit none

  integer, intent(in):: nI,nJ,nK,idim,iBlock
  real, intent(in)   :: StateCons_CV(nI,nJ,nK,nw)
  real, intent(in)   :: B0_CD(nI,nJ,nK,nDim)
  real, intent(out)  :: Flux_CV(nI,nJ,nK,nw)

  real :: Primitive_V(nw), Conservative_V(nw+1), Flux_V(nw+1)

  real :: Un, En, HallUn

  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------

  if(iBlock==BLKtest .and. iProc==PROCtest)then
     call set_oktest('get_face_flux', DoTest, DoTestMe)
  else
     DoTest=.false.; DoTestMe=.false.
  end if

  iBlockFace = iBlock

  call set_block_values(iDim)
  do kFace = 1, nK; do jFace = 1, nJ; do iFace = 1, nI

     DoTestCell = DoTestMe .and. &
          iFace==iTest .and. jFace==jTest .and. kFace==kTest

     Primitive_V = StateCons_CV(iFace, jFace, kFace, :)
     call conservative_to_primitive(Primitive_V)

     Conservative_V(1:nw) = StateCons_CV(iFace, jFace, kFace, :)
     Conservative_V(p_)     = Primitive_V(p_)
     Conservative_V(nw+1)   = StateCons_CV(iFace, jFace, kFace, E_)

     if(UseHallResist)then
        HallJx = HallJ_CD(iFace, jFace, kFace, x_)
        HallJy = HallJ_CD(iFace, jFace, kFace, y_)
        HallJz = HallJ_CD(iFace, jFace, kFace, z_)
     end if

     select case(iDim)
     case(x_)
        call set_cell_values_x
     case(y_)
        call set_cell_values_y
     case(z_)
        call set_cell_values_z
     end select
     call get_physical_flux(iDim, Primitive_V, &
          B0_CD(iFace, jFace, kFace, x_), &
          B0_CD(iFace, jFace, kFace, y_), &
          B0_CD(iFace, jFace, kFace, z_), &
          Conservative_V, Flux_V, Un, En, HallUn)

     Flux_CV(iFace, jFace, kFace, 1:E_-1)= Flux_V(1:E_-1)
     Flux_CV(iFace, jFace, kFace, E_)    = Flux_V(nw+1)

  end do; end do; end do

end subroutine get_face_flux

!==============================================================================
subroutine get_cmax_face(w,B0,qnI,qnJ,qnK,iDim,iBlock,cmax)

  use ModProcMH,   ONLY: iProc
  use ModMain,     ONLY: nDim, x_, y_, z_, ProcTest, BlkTest,iTest,jTest,kTest
  use ModImplicit, ONLY: nw
  use ModFaceFlux, ONLY: iBlockFace, iFace, jFace, kFace, &
       set_block_values, &
       set_cell_values_x, set_cell_values_y, set_cell_values_z, get_speed_max,&
       Area, DoLf, DoAw, DoRoe, DoHll, HallUnLeft, HallUnRight, DoTestCell

  implicit none

  integer, intent(in):: qnI,qnJ,qnK,idim,iBlock
  real, intent(in)   :: w(qnI,qnJ,qnK,nw)
  real, intent(in)   :: B0(qnI,qnJ,qnK,ndim)
  real, intent(out)  :: Cmax(qnI,qnJ,qnK)

  real :: Primitive_V(nw)

  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------

  if(iBlock==BLKtest .and. iProc==PROCtest)then
     call set_oktest('get_cmax_face', DoTest, DoTestMe)
  else
     DoTest=.false.; DoTestMe=.false.
  end if

  iBlockFace = iBlock
  DoLf  = .true.
  DoAw  = .false.
  DoRoe = .false.
  DoHll = .false.
  HallUnLeft = 0.0
  HallUnRight = 0.0

  call set_block_values(iDim)
  do kFace = 1, qnK; do jFace = 1, qnJ; do iFace = 1, qnI

     DoTestCell = DoTestMe .and. &
          iFace==iTest .and. jFace==jTest .and. kFace==kTest

     Primitive_V = w(iFace, jFace, kFace, :)

     call conservative_to_primitive(Primitive_V)

     select case(iDim)
     case(x_)
        call set_cell_values_x
     case(y_)
        call set_cell_values_y
     case(z_)
        call set_cell_values_z
     end select
     call get_speed_max(iDim, Primitive_V, &
          B0(iFace, jFace, kFace, x_), &
          B0(iFace, jFace, kFace, y_), &
          B0(iFace, jFace, kFace, z_), &
          cmax = cmax(iFace, jFace, kFace))
     cmax(iFace, jFace, kFace) = cmax(iFace, jFace, kFace)
  end do; end do; end do

end subroutine get_cmax_face

!==============================================================================
subroutine conservative_to_primitive(State_V)

  use ModImplicit, ONLY: nw, p_, e_
  use ModVarIndexes, ONLY: Rho_, Ux_, Uz_, RhoUx_, RhoUz_, Bx_, Bz_
  use ModPhysics, ONLY: gm1
  implicit none
  real, intent(inout):: State_V(nw)
  real :: InvRho
  !---------------------------------------------------------------------------
  InvRho = 1.0/State_V(Rho_)
  State_V(p_)      = gm1*(State_V(e_) - &
       0.5*(sum(State_V(RhoUx_:RhoUz_)**2)*InvRho + sum(State_V(Bx_:Bz_)**2)))
  State_V(Ux_:Uz_) = InvRho*State_V(RhoUx_:RhoUz_)

end subroutine conservative_to_primitive

!==============================================================================
subroutine getdt_courant(qdt)

  use ModProcMH
  use ModMain
  use ModAdvance, ONLY : B0xCell_BLK,B0yCell_BLK,B0zCell_BLK
  use ModGeometry, ONLY : dx_BLK,dy_BLK,dz_BLK,dxyz,true_cell,true_BLK,vInv_CB
  use ModImplicit
  use ModMpi
  implicit none

  real, intent(out) :: qdt

  real :: cmax(nI,nJ,nK), B0cell(nI,nJ,nK,ndim), qdt_local
  integer :: idim, implBLK, iBLK, iError

  logical :: oktest, oktest_me
  !-------------------------------------------------------------------------
  call set_oktest('getdt_courant',oktest,oktest_me)

  ! First calculate max(cmax/dx) for each cell and dimension
  qdt_local=0.0
  do implBLK=1,nImplBLK; 
     iBLK=impl2iBLK(implBLK); 
     dxyz(x_)=dx_BLK(iBLK); dxyz(y_)=dy_BLK(iBLK); dxyz(z_)=dz_BLK(iBLK)

     B0cell(:,:,:,x_)=B0xCell_BLK(1:nI,1:nJ,1:nK,iBLK)
     B0cell(:,:,:,y_)=B0yCell_BLK(1:nI,1:nJ,1:nK,iBLK)
     B0cell(:,:,:,z_)=B0zCell_BLK(1:nI,1:nJ,1:nK,iBLK)

     do idim=1,ndim

        call get_cmax_face(w_k(1:nI,1:nJ,1:nK,1:nw,implBLK),B0cell,&
             nI, nJ, nK, iDim, iBlk, Cmax)

        if(.not.true_BLK(iBLK))then
           where(.not.true_cell(1:nI,1:nJ,1:nK,iBLK))cmax=0.0
        end if

        qdt_local=max(qdt_local,maxval(cmax*vInv_CB(:,:,:,iBlk)))

        if(oktest_me)write(*,*)'getdt_courant idim,dx,cmax,1/qdt=',&
             idim,cmax(Itest,Jtest,Ktest),qdt_local
     end do
  end do

  ! Take global maximum
  call MPI_allreduce(qdt_local,qdt,1,MPI_REAL,MPI_MAX,iComm,iError)

  if(oktest_me)write(*,*)'1/dt_local,1/dt=',qdt_local,qdt

  ! Take inverse, and reduce so it is OK for 3D calculation
  qdt=0.3/qdt

  if(oktest_me)write(*,*)'getdt_courant final dt=',qdt

end subroutine getdt_courant

