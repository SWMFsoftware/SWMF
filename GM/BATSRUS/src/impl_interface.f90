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
subroutine getpthermal(w,qnI,qnJ,qnK,p)

  use ModVarIndexes, ONLY : rho_,rhoUx_,rhoUy_,rhoUz_,Bx_,By_,Bz_
  use ModImplicit, ONLY : nw,E_
  use ModPhysics, ONLY : gm1
  implicit none

  integer, intent(in):: qnI,qnJ,qnK
  real, intent(in)   :: w(qnI,qnJ,qnK,nw)
  real, intent(out)  :: p(qnI,qnJ,qnK)

  !--------------------------------------------------------------------------

  p=gm1*(w(:,:,:,e_)-0.5*(                                          &
       (w(:,:,:,rhoUx_)**2+w(:,:,:,rhoUy_)**2+w(:,:,:,rhoUz_)**2)   &
       /w(:,:,:,rho_)                                               &
        +w(:,:,:,Bx_)**2+w(:,:,:,By_)**2+w(:,:,:,Bz_)**2            &
       ))

end subroutine getpthermal
!==============================================================================

subroutine getptotal(w,qnI,qnJ,qnK,p)

  use ModVarIndexes, ONLY : Bx_,By_,Bz_
  use ModImplicit, ONLY : nw
  implicit none

  integer, intent(in):: qnI,qnJ,qnK
  real, intent(in)   :: w(qnI,qnJ,qnK,nw)
  real, intent(out)  :: p(qnI,qnJ,qnK)

  !--------------------------------------------------------------------------

  call getpthermal(w,qnI,qnJ,qnK,p)
  p=p+0.5*(w(:,:,:,Bx_)**2+w(:,:,:,By_)**2+w(:,:,:,Bz_)**2)

end subroutine getptotal

!==============================================================================
subroutine getcmax(w,B0,Dxyz_D,qnI,qnJ,qnK,idim,implBLK,cmax)

  use ModProcMH
  use ModMain
  use ModVarIndexes
  use ModImplicit
  use ModPhysics, ONLY : g
  use ModHallResist, ONLY: IonMassPerCharge, HallCmaxFactor
  use ModNumConst,   ONLY: cPi
  implicit none

  integer, intent(in):: qnI,qnJ,qnK,idim,implBLK
  real, intent(in)   :: w(qnI,qnJ,qnK,nw)
  real, intent(in)   :: B0(qnI,qnJ,qnK,ndim)
  real, intent(in)   :: Dxyz_D(ndim)
  real, intent(out)  :: cmax(qnI,qnJ,qnK)

  ! used to be automatic arrays
  real, dimension(:,:,:), allocatable :: csound2,cfast2
  integer :: iError

  logical :: oktest, oktest_me
  !--------------------------------------------------------------------------

  if(implBLK==implBLKtest.and.iProc==PROCtest)then
     call set_oktest('getcmax',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if

  ! Allocate arrays that used to be automatic
  allocate(csound2(qnI,qnJ,qnK), cfast2(qnI,qnJ,qnK), stat=iError)
  call alloc_check(iError,"getcmax arrays")

  if(oktest_me)then
     if(okdebug)then
        write(*,*)'Starting getcmax: ',&
             'idim,rho,E,Bx,By,Bz,B0z,B0y,B0z=', idim, &
             w(Itest,Jtest,Ktest,rho_),&
             w(Itest,Jtest,Ktest,E_),&
             w(Itest,Jtest,Ktest,Bx_),&
             w(Itest,Jtest,Ktest,By_),&
             w(Itest,Jtest,Ktest,Bz_),&
             B0(Itest,Jtest,Ktest,x_),&
             B0(Itest,Jtest,Ktest,y_),&
             B0(Itest,Jtest,Ktest,z_)
     else
        write(*,*)'Starting getcmax: idim=',idim
     end if
  end if

  ! csound^2 = g*p/rho
  call getpthermal(w,qnI,qnJ,qnK,csound2)

  if(oktest_me)write(*,*)'getcmax: p=',csound2(Itest,Jtest,Ktest)

  csound2=g*csound2/w(:,:,:,rho_)

  ! cfast^2 = csound^2 + B^2/rho  (cfast orthogonal to B)
  cfast2= csound2 + &
       ((w(:,:,:,Bx_)+B0(:,:,:,x_))**2                &
       +(w(:,:,:,By_)+B0(:,:,:,y_))**2                &
       +(w(:,:,:,Bz_)+B0(:,:,:,z_))**2)/w(:,:,:,rho_)

  if(oktest_me)then
     write(*,*)'getcmax: csound2, cfast2=',csound2(Itest,Jtest,Ktest),&
          cfast2(Itest,Jtest,Ktest)
     write(*,*)'getcmax: discr=',cfast2(Itest,Jtest,Ktest)**2-&
          4*csound2(Itest,Jtest,Ktest)*&
          (w(Itest,Jtest,Ktest,B_+idim)+B0(Itest,Jtest,Ktest,idim))**2/&
          w(Itest,Jtest,Ktest,rho_)
     write(*,*)'getcmax: min(discr)=',&
          minval(cfast2**2-4*csound2*(w(:,:,:,B_+idim)+B0(:,:,:,idim))**2/&
          w(:,:,:,rho_)),&
          minloc(cfast2**2-4*csound2*(w(:,:,:,B_+idim)+B0(:,:,:,idim))**2/&
          w(:,:,:,rho_))
  end if

  ! cmax = |ux| + sqrt(0.5*(cfast^2 + sqrt(cfast^4 - 4*csound^2*Bx^2/rho)))
  cmax= abs(w(:,:,:,rhoU_+idim)/w(:,:,:,rho_)) + &
       sqrt(0.5*(cfast2+ sqrt(max(0.0,cfast2**2          &
       -4*csound2*(w(:,:,:,B_+idim)+B0(:,:,:,idim))**2/w(:,:,:,rho_)))))

  if(HallCmaxFactor>0.0) cmax = cmax + HallCmaxFactor* &
          cPi*abs(w(:,:,:,B_+idim)+B0(:,:,:,idim))*IonMassPerCharge &
          /(w(:,:,:,rho_)*Dxyz_D(iDim))

  if(oktest_me)write(*,*)'Finished getcmax: cmax=',cmax(Itest,Jtest,Ktest)

  ! Deallocate arrays that used to be automatic
  deallocate(csound2, cfast2)

end subroutine getcmax

!==============================================================================
subroutine getdt_courant(qdt)

  use ModProcMH
  use ModMain
  use ModAdvance, ONLY : B0xCell_BLK,B0yCell_BLK,B0zCell_BLK
  use ModGeometry, ONLY : dx_BLK,dy_BLK,dz_BLK,dxyz,true_cell,true_BLK
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

        call getcmax(w_k(1:nI,1:nJ,1:nK,1:nw,implBLK),B0cell,dxyz,&
             nI,nJ,nK,idim,implBLK,cmax)

        if(.not.true_BLK(iBLK))then
           where(.not.true_cell(1:nI,1:nJ,1:nK,iBLK))cmax=0.0
        end if

        qdt_local=max(qdt_local,maxval(cmax)/dxyz(idim))

        if(oktest_me)write(*,*)'getdt_courant idim,dx,cmax,1/qdt=',&
             idim,dxyz(idim),cmax(Itest,Jtest,Ktest),qdt_local
     end do
  end do

  ! Take global maximum
  call MPI_allreduce(qdt_local,qdt,1,MPI_REAL,MPI_MAX,iComm,iError)

  if(oktest_me)write(*,*)'1/dt_local,1/dt=',qdt_local,qdt

  ! Take inverse, and reduce so it is OK for 3D calculation
  qdt=0.3/qdt

  if(oktest_me)write(*,*)'getdt_courant final dt=',qdt

end subroutine getdt_courant

!==============================================================================
subroutine getflux(w,B0,qnI,qnJ,qnK,iw,idim,implBLK,f)

  ! Calculate flux f=F_idim[iw] from w and B0 

  use ModProcMH
  use ModMain
  use ModVarIndexes
  use ModImplicit
  use ModHallResist, ONLY: UseHallResist, HallJ_CD
  implicit none

  integer, intent(in) :: qnI,qnJ,qnK,iw,idim,implBLK
  real,    intent(in) :: w(qnI,qnJ,qnK,nw)
  real,    intent(in) :: B0(qnI,qnJ,qnK,ndim)
  real,    intent(out):: f(qnI,qnJ,qnK)

  integer :: kdim, iBLK
  logical :: oktest, oktest_me
  !--------------------------------------------------------------------------
  if(iProc==PROCtest.and.implBLK==implBLKtest.and.iw==VARtest &
       .and.iDim==DimTest)then
     call set_oktest('getflux',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  endif
  if(oktest_me)write(*,*)'Getflux idim,w:',idim,w(Itest,Jtest,Ktest,VARtest)

  iBLK = impl2iBLK(implBLK)

  select case(iw)
  case(rho_)
     ! f_i[rho]=m_i
     f = w(:,:,:,rhoU_+idim)

  case(rhoUx_, rhoUy_, rhoUz_)
     ! f_i[m_k]=m_i*m_k/rho - b_k*b_i [+ptotal if i==k]
     !         -B0_k*b_i - B0_i*b_k [+B0_j*b_j if i==k]
     kdim = iw-rhoU_
     if(idim==kdim)then
        ! f_i[m_i]=ptotal + m_i**2/rho - b_i**2 - 2*B0_i*b_i + B0_j*b_j
        call getptotal(w,qnI,qnJ,qnK,f)
        f = f + w(:,:,:,iw)**2/w(:,:,:,rho_) - w(:,:,:,B_+idim)**2 &
             -2*B0(:,:,:,idim)*w( :,:,:,B_+idim) &
             +B0(:,:,:,x_)*w(:,:,:,Bx_)          &
             +B0(:,:,:,y_)*w(:,:,:,By_)          &
             +B0(:,:,:,z_)*w(:,:,:,Bz_)
     else
        ! f_i[m_k]=m_i*m_k/rho - b_k*b_i -B0_k*b_i - B0_i*b_k
        f = w(:,:,:,rhoU_+idim)*w(:,:,:,iw)/w(:,:,:,rho_) &
             -w( :,:,:,B_+kdim)*w(:,:,:,B_+idim)          &
             -B0(:,:,:,kdim)   *w(:,:,:,B_+idim)          &
             -B0(:,:,:,idim)   *w(:,:,:,B_+kdim)
     endif

  case(E_)
     ! f_i[e]=(m_i*(ptotal+e+(b_k*B0_k))-(b_i+B0_i)*(b_k*m_k))/rho 
     call getptotal(w,qnI,qnJ,qnK,f)
     f = (w(:,:,:,rhoU_+idim)                 & ! (m_i
          *(f + w( :,:,:,E_)                  & ! *(ptotal + e
          +B0(:,:,:,x_)*w(:,:,:,Bx_)          & !   +B0.b)
          +B0(:,:,:,y_)*w(:,:,:,By_)          &
          +B0(:,:,:,z_)*w(:,:,:,Bz_))         & 
          -(w(:,:,:,B_+idim)+B0(:,:,:,idim))  & ! -(b_i+B0_i)
          *(w(:,:,:,Bx_)*w(:,:,:,rhoUx_)      & !  *(b.m)
          +w(:,:,:,By_)*w(:,:,:,rhoUy_)       &
          +w(:,:,:,Bz_)*w(:,:,:,rhoUz_))      & ! 
          )/w(:,:,:,rho_)                       ! )/rho

  case(Bx_,By_,Bz_)
     kdim = iw-B_
     if(idim==kdim) then
        ! f_i[b_i] should be exactly 0
        f = 0.0
     else
        if(UseHallResist)then
           ! Take Hall effect into account. The off-diagonal part of 
           ! Eta.J = H*(J x B)/rho where H = HallFactor*IonMassPerCharge
           ! so the electric field (ignoring the diagonal part of J) becomes 
           ! E = -u x B + eta.J = -(rho U - H*J)xB/rho
           ! In effect the momentum rhoU is preplaced with (rhoU - H*J)
           !
           ! f_i[b_k]=( (m_i-H*J_i)*(b_k+B0_k) - (m_k-H*J_k)*(b_i+B0_i))/rho
           f =  ( (w(:,:,:,rhoU_+idim)-HallJ_CD(:,:,:,idim)) &
                *(w(:,:,:,iw)+ B0(:,:,:,kdim)) &
                -(w(:,:,:,rhoU_+kdim)-HallJ_CD(:,:,:,kdim)) &
                *(w(:,:,:,B_+idim) + B0(:,:,:,idim)) &
                )/w(:,:,:,rho_)
        else
           ! f_i[b_k]=(m_i*(b_k+B0_k) - m_k*(b_i+B0_i))/rho
           f =  (w(:,:,:,rhoU_+idim)*(w(:,:,:,iw)      + B0(:,:,:,kdim)) &
                -w(:,:,:,rhoU_+kdim)*(w(:,:,:,B_+idim) + B0(:,:,:,idim)) &
                )/w(:,:,:,rho_)
        end if
     endif
  case default
     ! We assume that all other variables behave like advected scalars !!!
     ! f_i[scalar]=m_i/rho*scalar
     f = w(:,:,:,rhoU_+idim)/w(:,:,:,rho_)*w(:,:,:,iw)
  end select

  if(oktest_me)write(*,*)'getflux final f=',f(Itest,Jtest,Ktest)

end subroutine getflux
