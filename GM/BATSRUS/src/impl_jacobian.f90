!^CFG COPYRIGHT UM
!^CFG FILE IMPLICIT
subroutine impl_jacobian(implBLK,JAC)

  ! Calculate Jacobian matrix for block implBLK:
  !
  !    JAC = I - dt*beta*dR/dw
  !
  ! using 1st order Rusanov scheme for the residual RES:
  !
  !    R_i = 1./Dx*[                 -(Fx_i+1/2 - Fx_i-1/2 )
  !                + 0.5*cmax_i+1/2 * (W_i+1    - W_i      )
  !                + 0.5*cmax_i-1/2 * (W_i-1    - W_i      )
  !                + 0.5*Q_i        * (Bx_i+1   - Bx_i-1   ) ]
  !         +1./Dy*[...]
  !         +1./Dz*[...]
  !         +S
  !
  ! where W contains the conservative variables, 
  ! Fx, Fy, Fz are the fluxes, cmax is the maximum speed,
  ! Q are the coefficients B, U and U.B in the Powell source terms, and
  ! S are the local source terms.
  !
  !    Fx_i+1/2 = 0.5*(FxL_i+1/2 + FxR_i+1/2)
  !
  ! For first order scheme 
  !
  !    FxL_i+1/2 = Fx[ W_i  , B0_i+1/2 ]
  !    FxR_i+1/2 = Fx[ W_i+1, B0_i+1/2 ]
  !
  ! We neglect terms containing d(cmax)/dW, and obtain a generalized eq.18:
  !
  ! Main diagonal stencil==1:
  !    dR_i/dW_i   = 0.5/Dx*[ (dFxR/dW-cmax)_i-1/2 - (dFxL/dW+cmax)_i+1/2 ]
  !                + 0.5/Dy*[ (dFyR/dW-cmax)_j-1/2 - (dFyL/dW+cmax)_j+1/2 ]
  !                + 0.5/Dz*[ (dFzR/dW-cmax)_k-1/2 - (dFzL/dW+cmax)_k+1/2 ]
  !                + dQ/dW_i*divB
  !                + dS/dW
  !
  ! Subdiagonal stencil==2:
  !    dR_i/dW_i-1 = 0.5/Dx* [ (dFxL/dW+cmax)_i-1/2
  !                           - Q_i*dBx_i-1/dW_i-1 ]
  ! Superdiagonal stencil==3:
  !    dR_i/dW_i+1 = 0.5/Dx* [-(dFxR/dW-cmax)_i+1/2
  !                           + Q_i*dBx_i+1/dW_i+1 ]
  !
  ! and similar terms for stencil=4,5,6,7.
  ! 
  ! The partial derivatives are calculated numerically 
  ! (except for the trivial dQ/dW and dB/dW terms):
  !
  !  dF_iw/dW_jw = [F_iw(W + eps*W_jw) - F_iw(W)] / eps
  !  dS_iw/dW_jw = [S_iw(W + eps*W_jw) - S_iw(W)] / eps

  use ModProcMH
  use ModMain
  use ModvarIndexes
  use ModAdvance, ONLY : B0_DX, B0_DY,B0_DZ,set_b0_face,&
       time_BLK
  use ModGeometry, ONLY : dx_BLK,dy_BLK,dz_BLK,dxyz,true_cell, &
       FaceAreaI_DFB, FaceAreaJ_DFB, FaceAreaK_DFB
  use ModImplicit
  use ModHallResist, ONLY: UseHallResist, hall_factor
  use ModHeatConduction, ONLY: UseParallelConduction, impl_jac_heat_conduction
  use ModGeometry, ONLY: vInv_CB, UseCovariant
  implicit none

  integer, intent(in) :: implBLK
  real,    intent(out):: JAC(nw,nw,nI,nJ,nK,nstencil)

  integer :: iBLK
  real, dimension(nw,nI,nJ,nK)                :: Impl_VC      , ImplEps_VC
  real, dimension(nI+1,nJ+1,nK+1)             :: dfdwLface, dfdwRface
  real, dimension(nw,nI,nJ,nK)                :: s_VC, sEps_VC, sPowell_VC
  real :: DivB(nI,nJ,nK)
  real :: B0_DFD(nDim,nI+1,nJ+1,nK+1,nDim), Cmax_DF(nDim,nI+1,nJ+1,nK+1)

  real   :: qeps, coeff
  logical:: divbsrc, UseDivbSource0
  integer:: i,j,k,i1,i2,i3,j1,j2,j3,k1,k2,k3,istencil,iw,jw,idim,qj

  logical :: oktest, oktest_me

  real :: FluxLeft_VFD(nW,nI+1,nJ+1,nK+1,nDim) ! Unperturbed left flux
  real :: FluxRight_VFD(nW,nI,nJ,nK,nDim)      ! Unperturbed right flux
  real :: FluxEpsLeft_VF(nW,nI+1,nJ+1,nK+1)    ! Perturbed left flux
  real :: FluxEpsRight_VF(nW,nI,nJ,nK)         ! Perturbed right flux
  real :: FaceArea_F(nI, nJ, nK)               ! Only the inner faces

  real :: HallFactor_G(0:nI+1,0:nJ+1,0:nK+1)
  !----------------------------------------------------------------------------

  if(iProc==PROCtest.and.implBLK==implBLKtest)then
     call set_oktest('impl_jacobian',oktest, oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if

  qeps=sqrt(JacobianEps)

  divbsrc= UseDivbSource

  ! Extract state for this block
  Impl_VC = Impl_VGB(1:nw,1:nI,1:nJ,1:nK,implBLK)
  iBLK = impl2iBLK(implBLK)
  dxyz(1)=dx_BLK(iBLK); dxyz(2)=dy_BLK(iBLK); dxyz(3)=dz_BLK(iBLK)
  if(UseB0)then
     call set_b0_face(iBLK)
     B0_DFD(:,1:nI+1,1:nJ  ,1:nK  ,x_)=B0_DX(:,1:nI+1,1:nJ,1:nK)
     B0_DFD(:,1:nI  ,1:nJ+1,1:nK  ,y_)=B0_DY(:,1:nI,1:nJ+1,1:nK)
     B0_DFD(:,1:nI  ,1:nJ  ,1:nK+1,z_)=B0_DZ(:,1:nI,1:nJ,1:nK+1)
  end if

  if(UseHallResist)call impl_init_hall

  ! Initialize matrix to zero (to be safe)
  JAC=0.0

  ! Initialize reference flux and the cmax array
  do idim=1,ndim

     i1 = 1+kr(1,idim); j1= 1+kr(2,idim); k1= 1+kr(3,idim)
     i2 =nI+kr(1,idim); j2=nJ+kr(2,idim); k2=nK+kr(3,idim);

     call get_face_flux(Impl_VC,B0_DFD(:,i1:i2,j1:j2,k1:k2,idim),&
          nI,nJ,nK,iDim,iBLK,FluxLeft_VFD(:,i1:i2,j1:j2,k1:k2,iDim))

     call get_face_flux(Impl_VC,B0_DFD(:,1:nI,1:nJ,1:nK,iDim),&
          nI,nJ,nK,iDim,iBLK,FluxRight_VFD(:,:,:,:,iDim))

     ! Average w for each cell interface into ImplEps_VC
     i1 = 1-kr(1,idim); j1= 1-kr(2,idim); k1= 1-kr(3,idim)

     ! Calculate orthogonal cmax for each interface in advance
     call get_cmax_face(                            &
          0.5*(Impl_VGB(1:nw, 1:i2, 1:j2, 1:k2,implBLK)+ &
          Impl_VGB(1:nw,i1:nI,j1:nJ,k1:nK,implBLK)),     &
          B0_DFD(1:ndim,1:i2,1:j2,1:k2,idim),       &
          i2,j2,k2,idim,iBlk,Cmax_DF(idim,1:i2,1:j2,1:k2))

     ! cmax always occurs as -ImplCoeff*0.5/dx*cmax
     coeff = -ImplCoeff*0.5 
     Cmax_DF(idim,1:i2,1:j2,1:k2)=coeff*Cmax_DF(idim,1:i2,1:j2,1:k2)
  enddo

  ! Initialize divB and sPowell_VC arrays 
  if(UseB .and. divbsrc)call impl_divbsrc_init

  ! Set s_VC=S(Impl_VC)
  if(implsource)call getsource(iBLK,Impl_VC,s_VC)

  ! The w to be perturbed and jw is the index for the perturbed variable
  ImplEps_VC=Impl_VC

  !DEBUG
  !write(*,*)'Initial ImplEps_VC=Impl_VC'
  !write(*,'(a,8(f10.6))')'Impl_VC(nK)   =', Impl_VC(:,Itest,Jtest,Ktest)
  !write(*,'(a,8(f10.6))')'Impl_VC(nK+1) =', Impl_VC(:,Itest,Jtest,Ktest+1)
  !write(*,'(a,8(f10.6))')'ImplEps_VC(nK)  =',ImplEps_VC(:,Itest,Jtest,Ktest)
  !write(*,'(a,8(f10.6))')'ImplEps_VC(nK+1)=',ImplEps_VC(:,Itest,Jtest,Ktest+1)

  do jw=1,nw; 
     ! Remove perturbation from previous jw if there was a previous one
     if(jw>1)ImplEps_VC(jw-1,:,:,:)=Impl_VC(jw-1,:,:,:)

     ! Perturb new jw variable
     coeff=qeps*wnrm(jw)
     ImplEps_VC(jw,:,:,:)=Impl_VC(jw,:,:,:) + coeff

     do idim=1,ndim
        ! Index limits for faces and shifted centers
        i1 = 1+kr(1,idim); j1= 1+kr(2,idim); k1= 1+kr(3,idim)
        i2 =nI+kr(1,idim); j2=nJ+kr(2,idim); k2=nK+kr(3,idim);
        i3 =nI-kr(1,idim); j3=nJ-kr(2,idim); k3=nK-kr(3,idim);

        call get_face_flux(ImplEps_VC,B0_DFD(:,i1:i2,j1:j2,k1:k2,idim),&
             nI,nJ,nK,iDim,iBLK,FluxEpsLeft_VF(:,i1:i2,j1:j2,k1:k2))

        call get_face_flux(ImplEps_VC,B0_DFD(:,1:nI,1:nJ,1:nK,iDim),&
             nI,nJ,nK,iDim,iBLK,FluxEpsRight_VF)

        ! Calculate dfdw=(feps-f0)/eps for each iw variable and both
        ! left and right sides
        do iw=1,nw

           !call getflux(ImplEps_VC,B0_DFD(:,i1:i2,j1:j2,k1:k2,idim),&
           !     nI,nJ,nK,iw,idim,implBLK,fepsLface(i1:i2,j1:j2,k1:k2))

           !call getflux(ImplEps_VC,B0_DFD(:,1:nI,1:nJ,1:nK,idim),&
           !     nI,nJ,nK,iw,idim,implBLK,fepsRface(1:nI,1:nJ,1:nK))

           ! dfdw = F_iw(W + eps*W_jw) - F_iw(W)] / eps is multiplied by 
           ! -ImplCoeff/2/dx*wnrm(jw)/wnrm(iw) in all formulae
           coeff=-ImplCoeff*0.5/qeps/wnrm(iw) 

           dfdwLface(i1:i2,j1:j2,k1:k2)=coeff*&
                (FluxEpsLeft_VF(iw,i1:i2,j1:j2,k1:k2) &
                -FluxLeft_VFD(iw,   i1:i2,j1:j2,k1:k2,iDim)) 
           dfdwRface( 1:nI, 1:nJ, 1:nK)=coeff*&
                (FluxEpsRight_VF( iW, 1:nI, 1:nJ, 1:nK) &
                -FluxRight_VFD( iW,   1:nI, 1:nJ, 1:nK, iDim))

           if(oktest_me)write(*,'(a,i1,i2,6(f15.8))') &
                'iw,jw,f0L,fepsL,dfdwL,R:', &
                iw,jw,&
                FluxLeft_VFD(iw,Itest,Jtest,Ktest,idim),&
                FluxEpsLeft_VF(iW,Itest,Jtest,Ktest),&
                dfdwLface(Itest,Jtest,Ktest),&
                FluxRight_VFD(iw,Itest,Jtest,Ktest,idim),&
                FluxEpsRight_VF(iW,Itest,Jtest,Ktest),&
                dfdwRface(Itest,Jtest,Ktest)

           !DEBUG
           !if(idim==3.and.iw==4.and.jw==2)&
           !write(*,*)'BEFORE addcmax dfdw(iih)=',&
           !          dfdwLface(Itest,Jtest,Ktest-1)

           ! Add contribution of cmax to dfdwL and dfdwR
           if(iw==jw)then
              ! FxL_i-1/2 <-- (FxL + cmax)_i-1/2
              dfdwLface(i1:i2,j1:j2,k1:k2)=dfdwLface(i1:i2,j1:j2,k1:k2)&
                   +Cmax_DF(idim,i1:i2,j1:j2,k1:k2)
              ! FxR_i+1/2 <-- (FxR - cmax)_i+1/2
              dfdwRface( 1:nI, 1:nJ, 1:nK)=dfdwRface( 1:nI, 1:nJ, 1:nK)&
                   -Cmax_DF(idim, 1:nI, 1:nJ, 1:nK)
           endif

           ! Divide flux*area by volume
           dfdwLface(i1:i2,j1:j2,k1:k2) = dfdwLface(i1:i2,j1:j2,k1:k2) &
                *vInv_CB(1:nI, 1:nJ, 1:nK, iBlk)
           dfdwRface( 1:nI, 1:nJ, 1:nK) = dfdwRface( 1:nI, 1:nJ, 1:nK) &
                *vInv_CB(1:nI, 1:nJ, 1:nK, iBlk)

           !DEBUG
           !if(idim==3.and.iw==4.and.jw==2)&
           !write(*,*)'AFTER  addcmax dfdwL(iih)=',&
           !           dfdwLface(Itest,Jtest,Ktest-1)

           ! Contribution of fluxes to main diagonal (middle cell)
           ! dR_i/dW_i = 0.5/Dx*[ (dFxR/dW-cmax)_i-1/2 - (dFxL/dW+cmax)_i+1/2 ]

           JAC(iw,jw,:,:,:,1)=JAC(iw,jw,:,:,:,1) &
                +dfdwRface( 1:nI, 1:nJ, 1:nK)                      &
                -dfdwLface(i1:i2,j1:j2,k1:k2)


           ! Add Q*dB/dw to dfdwL and dfdwR for upper and lower diagonals
           ! These diagonals are non-zero for the inside interfaces only
           ! which corresponds to the range i1:nI,j1:nJ,k1:nK.
           if(UseB    .and. divbsrc .and. &
                (     (iw >= RhoUx_ .and. iw <= RhoUz_) &
                .or.  (iW >= Bx_    .and. iW <= Bz_   ) &
                .or.   iw == E_                         &
                ) )then
              if(UseCovariant .and. jw>=Bx_ .and. jw<=Bz_)then
                 ! The source terms are always multiplied by coeff
                 coeff=-ImplCoeff*0.5*wnrm(jw)/wnrm(iw)
                 ! Get the corresponding face area
                 select case(iDim)
                 case(x_)
                    FaceArea_F(i1:nI,j1:nJ,k1:nK) = &
                         FaceAreaI_DFB(jw-B_,i1:nI,j1:nJ,k1:nK,iBLK)
                 case(y_)
                    FaceArea_F(i1:nI,j1:nJ,k1:nK) = &
                         FaceAreaJ_DFB(jw-B_,i1:nI,j1:nJ,k1:nK,iBLK)
                 case(z_)
                    FaceArea_F(i1:nI,j1:nJ,k1:nK) = &
                         FaceAreaK_DFB(jw-B_,i1:nI,j1:nJ,k1:nK,iBlk)
                 end select

                 ! Relative to the right face flux Q is shifted to the left
                 dfdwLface(i1:nI,j1:nJ,k1:nK)=dfdwLface(i1:nI,j1:nJ,k1:nK)+ &
                      coeff*sPowell_VC(iw,i1:nI,j1:nJ,k1:nK) &
                      *FaceArea_F(i1:nI,j1:nJ,k1:nK) &
                      *vInv_CB(i1:nI,j1:nJ,k1:nK,iBlk)

                 dfdwRface(i1:nI,j1:nJ,k1:nK)=dfdwRface(i1:nI,j1:nJ,k1:nK)+ &
                      coeff*sPowell_VC(iw, 1:i3, 1:j3, 1:k3) &
                      *FaceArea_F(i1:nI,j1:nJ,k1:nK) &
                      *vInv_CB(1:i3,1:j3,1:k3,iBlk)

              elseif(jw==B_+idim)then
                 ! The source terms are always multiplied by coeff
                 coeff=-ImplCoeff*0.5/dxyz(idim)*wnrm(jw)/wnrm(iw)

                 ! Relative to the right face flux Q is shifted to the left
                 dfdwLface(i1:nI,j1:nJ,k1:nK)=dfdwLface(i1:nI,j1:nJ,k1:nK)+ &
                      coeff*sPowell_VC(iw,i1:nI,j1:nJ,k1:nK)

                 dfdwRface(i1:nI,j1:nJ,k1:nK)=dfdwRface(i1:nI,j1:nJ,k1:nK)+ &
                      coeff*sPowell_VC(iw, 1:i3, 1:j3, 1:k3)
              end if
           end if
           JAC(iw,jw,i1:nI,j1:nJ,k1:nK,2*idim  )=  dfdwLface(i1:nI,j1:nJ,k1:nK)
           JAC(iw,jw, 1:i3, 1:j3, 1:k3,2*idim+1)= -dfdwRface(i1:nI,j1:nJ,k1:nK)
        enddo ! iw
     enddo ! idim
     if(oktest_me)then
        write(*,*)'After fluxes jw=',jw,' stencil, row, JAC'
        do istencil=1,nstencil
           do qj=1,nw
              write(*,'(i1,a,i1,a,20(f9.5))')istencil,',',qj,':',&
                   JAC(:,qj,Itest,Jtest,Ktest,istencil)
           end do
        enddo
     endif

     !Derivatives of local source terms 
     if(implsource)then
        if(oktest_me)write(*,*)'Adding dS/dw'

        ! w2=S(Impl_VC+eps*W_jw)
        call getsource(iBLK,ImplEps_VC,sEps_VC)
        do iw=1,nw
           ! JAC(..1) += dS/dW_jw
           coeff=-ImplCoeff/qeps/wnrm(iw)
           JAC(iw,jw,:,:,:,1)=JAC(iw,jw,:,:,:,1)&
                +coeff*(sEps_VC(iw,:,:,:)-s_VC(iw,:,:,:))
        enddo
     endif
  enddo

  if(oktest_me)write(*,*)'After fluxes and sources:  JAC(...,1):', &
       JAC(1:nw,1:nw,Itest,Jtest,Ktest,1)

  ! Contribution of middle to Powell's source terms
  if(UseB .and. divbsrc)then
     ! JAC(...1) += d(Q/divB)/dW*divB
     call impl_divbsrc_middle

     if(oktest_me)then
        write(*,*)'After divb sources: row, JAC(...,1):'
        do qj=1,nw
           write(*,'(i1,a,20(f9.5))')qj,':',&
                JAC(:,qj,Itest,Jtest,Ktest,1)
        end do
     end if
  end if

  ! Add extra terms for Hall resistivity
  if(UseHallResist .and. .not. UseCovariant)call impl_hall_resist
  if(UseHallResist .and.       UseCovariant)call impl_hall_resist_general

  if(UseGrayDiffusion) call impl_gray_diffusion
  if(UseParallelConduction) call impl_jac_heat_conduction(iBLK, nw, JAC)

  ! Multiply JAC by the implicit timestep dt
  if(time_accurate)then
     do k=1,nK; do j=1,nJ; do i=1,nI
        if(true_cell(i,j,k,iBLK))then
           JAC(:,:,i,j,k,:) = JAC(:,:,i,j,k,:)*dt
        else
           ! Set JAC = 0.0 inside body
           JAC(:,:,i,j,k,:) = 0.0
        end if
     end do; end do; end do
  else
     ! Local time stepping has time_BLK=0.0 inside the body
     do istencil=1,nstencil; do jw=1,nw; do iw=1,nw
        JAC(iw,jw,:,:,:,istencil)=JAC(iw,jw,:,:,:,istencil) &
             *time_BLK(1:nI,1:nJ,1:nK,iBLK)*implCFL
     end do; end do; end do
  endif

  if(oktest_me)then
     write(*,*)'After boundary correction and *dt: row, JAC(...,1):'
     do qj=1,nw
        write(*,'(i1,a,20(f9.5))')qj,':',&
             JAC(:,qj,Itest,Jtest,Ktest,1)
     end do
  end if

  ! Add unit matrix to main diagonal
  do iw=1,nw
     JAC(iw,iw,:,:,:,1)=JAC(iw,iw,:,:,:,1)+1.0
  end do

  if(oktest_me)then
     write(*,*)'After adding I: row, JAC(...,1):'
     do qj=1,nw
        write(*,'(i1,a,20(f9.5))')qj,':',&
             JAC(:,qj,Itest,Jtest,Ktest,1)
     end do
  end if

  ! Restore UseDivbSource
  if(UseB .and. divbsrc)UseDivbSource=UseDivbSource0

contains

  !===========================================================================
  subroutine impl_divbsrc_init

    ! Switch off UseDivbSource for addsource
    UseDivbSource0=UseDivbSource
    UseDivbSource=.false.

    ! Calculate div B for middle cell contribution to Powell's source terms
    if(UseCovariant)then
       do k=1,nK; do j=1,nJ; do i=1,nI
          divb(i,j,k) = 0.5*vInv_CB(i,j,k,iBlk)*( &
               sum(Impl_VGB(Bx_:Bz_,i+1,j,k,implBLK)*FaceAreaI_DFB(:,i+1,j,k,iBlk))&
               -sum(Impl_VGB(Bx_:Bz_,i-1,j,k,implBLK)*FaceAreaI_DFB(:,i,j,k,iBlk))+&
               sum(Impl_VGB(Bx_:Bz_,i,j+1,k,implBLK)*FaceAreaJ_DFB(:,i,j+1,k,iBlk))&
               -sum(Impl_VGB(Bx_:Bz_,i,j-1,k,implBLK)*FaceAreaJ_DFB(:,i,j,k,iBlk))+&
               sum(Impl_VGB(Bx_:Bz_,i,j,k+1,implBLK)*FaceAreaK_DFB(:,i,j,k+1,iBlk))&
               -sum(Impl_VGB(Bx_:Bz_,i,j,k-1,implBLK)*FaceAreaK_DFB(:,i,j,k,iBlk))) 
       end do; end do; end do
   else
       divb=0.5*&
            ((Impl_VGB(Bx_,2:nI+1,1:nJ,1:nK,implBLK)           &
            -Impl_VGB(Bx_,0:nI-1,1:nJ,1:nK,implBLK))/dxyz(x_) &
            +(Impl_VGB(By_,1:nI,2:nJ+1,1:nK,implBLK)           &
            -Impl_VGB(By_,1:nI,0:nJ-1,1:nK,implBLK))/dxyz(y_) &
            +(Impl_VGB(Bz_,1:nI,1:nJ,2:nK+1,implBLK)           &
            -Impl_VGB(Bz_,1:nI,1:nJ,0:nK-1,implBLK))/dxyz(z_))
    end if

    ! Make sure that sPowell_VC is defined for all indexes
    sPowell_VC = 0.0
    do k=1,nK; do j=1,nJ; do i=1,nI
       ! Calculate coefficients Q that multiply div B in Powell source terms
       ! Q(rhoU)= B
       sPowell_VC(RhoUx_:RhoUz_,i,j,k)=Impl_VC(Bx_:Bz_,i,j,k)

       ! Q(B)   = U
       sPowell_VC(Bx_:Bz_,i,j,k) = Impl_VC(RhoUx_:RhoUz_,i,j,k) &
            /Impl_VC(Rho_,i,j,k) 

       ! Q(E)   = U.B
       sPowell_VC(E_,i,j,k) = &
            sum(Impl_VC(Bx_:Bz_,i,j,k)*Impl_VC(RhoUx_:RhoUz_,i,j,k)) &
            /Impl_VC(Rho_,i,j,k)
    end do; end do; end do

  end subroutine impl_divbsrc_init

  !===========================================================================
  subroutine impl_divbsrc_middle

    integer:: i,j,k

    ! JAC(...1) += dQ/dW_i*divB

    ! Q(rhoU)= -divB*B
    ! dQ(rhoU)/dB = -divB
    do k=1,nK; do j=1,nJ; do i=1,nI
       JAC(rhoUx_,Bx_,i,j,k,1)=JAC(rhoUx_,Bx_,i,j,k,1)&
            +ImplCoeff*wnrm(Bx_)/wnrm(rhoUx_)*divb(i,j,k) 
       JAC(rhoUy_,By_,i,j,k,1)=JAC(rhoUy_,By_,i,j,k,1)&
            +ImplCoeff*wnrm(By_)/wnrm(rhoUy_)*divb(i,j,k) 
       JAC(rhoUz_,Bz_,i,j,k,1)=JAC(rhoUz_,Bz_,i,j,k,1)&
            +ImplCoeff*wnrm(Bz_)/wnrm(rhoUz_)*divb(i,j,k) 

       ! Q(B)= -divB*rhoU/rho
       ! dQ(B)/drho = +divB*rhoU/rho**2
       JAC(Bx_,rho_,i,j,k,1)=JAC(Bx_,rho_,i,j,k,1) &
            -ImplCoeff*wnrm(rho_)/wnrm(Bx_)*divb(i,j,k) &
            *Impl_VC(rhoUx_,i,j,k)/Impl_VC(rho_,i,j,k)**2
       JAC(By_,rho_,i,j,k,1)=JAC(By_,rho_,i,j,k,1) &
            -ImplCoeff*wnrm(rho_)/wnrm(By_)*divb(i,j,k) &
            *Impl_VC(rhoUy_,i,j,k)/Impl_VC(rho_,i,j,k)**2
       JAC(Bz_,rho_,i,j,k,1)=JAC(Bz_,rho_,i,j,k,1) &
            -ImplCoeff*wnrm(rho_)/wnrm(Bz_)*divb(i,j,k) &
            *Impl_VC(rhoUz_,i,j,k)/Impl_VC(rho_,i,j,k)**2

       ! dQ(B)/drhoU= -divB/rho
       JAC(Bx_,rhoUx_,i,j,k,1)=JAC(Bx_,rhoUx_,i,j,k,1)&
            +ImplCoeff*wnrm(rhoUx_)/wnrm(Bx_)*divb(i,j,k)/Impl_VC(rho_,i,j,k) 
       JAC(By_,rhoUy_,i,j,k,1)=JAC(By_,rhoUy_,i,j,k,1)&
            +ImplCoeff*wnrm(rhoUy_)/wnrm(By_)*divb(i,j,k)/Impl_VC(rho_,i,j,k) 
       JAC(Bz_,rhoUz_,i,j,k,1)=JAC(Bz_,rhoUz_,i,j,k,1)&
            +ImplCoeff*wnrm(rhoUz_)/wnrm(Bz_)*divb(i,j,k)/Impl_VC(rho_,i,j,k) 

       ! Q(E)= -divB*rhoU.B/rho
       ! dQ(E)/drho = +divB*rhoU.B/rho**2
       JAC(E_,rho_,i,j,k,1)=JAC(E_,rho_,i,j,k,1)&
            -ImplCoeff*wnrm(rho_)/wnrm(E_)*divb(i,j,k)*&
            (Impl_VC(rhoUx_,i,j,k)*Impl_VC(Bx_,i,j,k)&
            +Impl_VC(rhoUy_,i,j,k)*Impl_VC(By_,i,j,k)&
            +Impl_VC(rhoUz_,i,j,k)*Impl_VC(Bz_,i,j,k))&
            /Impl_VC(rho_,i,j,k)**2

       ! dQ(E)/drhoU = -divB*B/rho
       JAC(E_,rhoUx_,i,j,k,1)=JAC(E_,rhoUx_,i,j,k,1) &
            +ImplCoeff*wnrm(rhoUx_)/wnrm(E_)*divb(i,j,k) &
            *Impl_VC(Bx_,i,j,k)/Impl_VC(rho_,i,j,k) 
       JAC(E_,rhoUy_,i,j,k,1)=JAC(E_,rhoUy_,i,j,k,1) &
            +ImplCoeff*wnrm(rhoUy_)/wnrm(E_)*divb(i,j,k) &
            *Impl_VC(By_,i,j,k)/Impl_VC(rho_,i,j,k) 
       JAC(E_,rhoUz_,i,j,k,1)=JAC(E_,rhoUz_,i,j,k,1) &
            +ImplCoeff*wnrm(rhoUz_)/wnrm(E_)*divb(i,j,k) &
            *Impl_VC(Bz_,i,j,k)/Impl_VC(rho_,i,j,k) 

       ! dQ(E)/dB = -divB*rhoU/rho
       JAC(E_,Bx_,i,j,k,1)=JAC(E_,Bx_,i,j,k,1) &
            +ImplCoeff*wnrm(Bx_)/wnrm(E_)*divb(i,j,k) &
            *Impl_VC(rhoUx_,i,j,k)/Impl_VC(rho_,i,j,k) 
       JAC(E_,By_,i,j,k,1)=JAC(E_,By_,i,j,k,1) &
            +ImplCoeff*wnrm(By_)/wnrm(E_)*divb(i,j,k) &
            *Impl_VC(rhoUy_,i,j,k)/Impl_VC(rho_,i,j,k) 
       JAC(E_,Bz_,i,j,k,1)=JAC(E_,Bz_,i,j,k,1) &
            +ImplCoeff*wnrm(Bz_)/wnrm(E_)*divb(i,j,k) &
            *Impl_VC(rhoUz_,i,j,k)/Impl_VC(rho_,i,j,k) 
    end do; end do; end do

  end subroutine impl_divbsrc_middle

  !===========================================================================
  subroutine impl_init_hall

    use ModHallResist, ONLY: HallJ_CD, IonMassPerCharge_G, &
         BxPerN_G, ByPerN_G, BzPerN_G, set_ion_mass_per_charge

    use ModAdvance, ONLY: B0_DGB
    
    use ModHallResist, ONLY: &                 
         DgenDxyz_DDC, set_block_jacobian_cell 
    real :: DbDgen_DD(3,3)                     

    real :: InvDx2, InvDy2, InvDz2, InvN_G(0:nI+1,0:nJ+1,0:nK+1)

    !----------------------------------------------------------------------
    ! Calculate cell centered currents to be used by getflux
    
    call set_ion_mass_per_charge(iBlk)

    InvDx2 = 0.5/Dxyz(x_); InvDy2 = 0.5/Dxyz(y_); InvDz2 = 0.5/Dxyz(z_)

    if(UseCovariant)then                      

       call set_block_jacobian_cell(iBlk)

       do k=1,nK; do j=1,nJ; do i=1,nI
          DbDgen_DD(:,1) = InvDx2* &
               (Impl_VGB(Bx_:Bz_,i+1,j,k,implBLK)-Impl_VGB(Bx_:Bz_,i-1,j,k,implBLK))
          DbDgen_DD(:,2) = InvDy2* &
               (Impl_VGB(Bx_:Bz_,i,j+1,k,implBLK)-Impl_VGB(Bx_:Bz_,i,j-1,k,implBLK))
          DbDgen_DD(:,3) = InvDz2* &
               (Impl_VGB(Bx_:Bz_,i,j,k+1,implBLK)-Impl_VGB(Bx_:Bz_,i,j,k-1,implBLK))

          ! Jx = Dbz/Dy - Dby/Dz
          HallJ_CD(i,j,k,x_) = &
               sum(DbDgen_DD(z_,:)*DgenDxyz_DDC(:,y_,i,j,k)) - &
               sum(DbDgen_DD(y_,:)*DgenDxyz_DDC(:,z_,i,j,k))

          ! Jy = Dbx/Dz - Dbz/Dx
          HallJ_CD(i,j,k,y_) = &
               sum(DbDgen_DD(x_,:)*DgenDxyz_DDC(:,z_,i,j,k)) - &
               sum(DbDgen_DD(z_,:)*DgenDxyz_DDC(:,x_,i,j,k))

          ! Jz = Dby/Dx - Dbx/Dy
          HallJ_CD(i,j,k,z_) = &
               sum(DbDgen_DD(y_,:)*DgenDxyz_DDC(:,x_,i,j,k)) - &
               sum(DbDgen_DD(x_,:)*DgenDxyz_DDC(:,y_,i,j,k))

       end do; end do; end do

    else                                        

       do k=1,nK; do j=1,nJ; do i=1,nI
          HallJ_CD(i,j,k,x_) = &
               +InvDy2*(Impl_VGB(Bz_,i,j+1,k,implBLK)-Impl_VGB(Bz_,i,j-1,k,implBLK)) &
               -InvDz2*(Impl_VGB(By_,i,j,k+1,implBLK)-Impl_VGB(By_,i,j,k-1,implBLK))
       end do; end do; end do

       do k=1,nK; do j=1,nJ; do i=1,nI
          HallJ_CD(i,j,k,y_) = &
               +InvDz2*(Impl_VGB(Bx_,i,j,k+1,implBLK)-Impl_VGB(Bx_,i,j,k-1,implBLK)) &
               -InvDx2*(Impl_VGB(Bz_,i+1,j,k,implBLK)-Impl_VGB(Bz_,i-1,j,k,implBLK))
       end do; end do; end do

       do k=1,nK; do j=1,nJ; do i=1,nI
          HallJ_CD(i,j,k,z_) = &
               +InvDx2*(Impl_VGB(By_,i+1,j,k,implBLK)-Impl_VGB(By_,i-1,j,k,implBLK)) &
               -InvDy2*(Impl_VGB(Bx_,i,j+1,k,implBLK)-Impl_VGB(Bx_,i,j-1,k,implBLK))
       end do; end do; end do
    end if                                    

    do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
       HallFactor_G(i,j,k) = hall_factor(0,i,j,k,iBlk)
    end do; end do; end do

    !write(*,*)'iBlock, max(IonMassPerCh), max(HallFactor), max(HallJ) =', &
    !     implBlk,&
    !     maxval(IonMassPerCharge_G(1:nI,1:nJ,1:nK)),&
    !     maxval(abs(HallJ_CD(:,:,:,:)))

    do k=1,nK; do j=1,nJ; do i=1,nI
       HallJ_CD(i,j,k,:) = IonMassPerCharge_G(i,j,k) &
            *HallFactor_G(i,j,k)*HallJ_CD(i,j,k,:)
    end do; end do; end do

    !write(*,*)'iBlock, max(HallJ)=',implBlk,maxval(abs(HallJ_CD(:,:,:,:)))

    InvN_G = HallFactor_G * IonMassPerCharge_G / &
         Impl_VGB(Rho_,0:nI+1,0:nJ+1,0:nK+1,implBLK)
    if(UseB0)then
       BxPerN_G = (B0_DGB(x_,0:nI+1,0:nJ+1,0:nK+1,iBlk) + &
            Impl_VGB(Bx_,0:nI+1,0:nJ+1,0:nK+1,implBLK))*InvN_G
       ByPerN_G = (B0_DGB(y_,0:nI+1,0:nJ+1,0:nK+1,iBlk) + &
            Impl_VGB(By_,0:nI+1,0:nJ+1,0:nK+1,implBLK))*InvN_G
       BzPerN_G = (B0_DGB(z_,0:nI+1,0:nJ+1,0:nK+1,iBlk) + &
            Impl_VGB(Bz_,0:nI+1,0:nJ+1,0:nK+1,implBLK))*InvN_G
    else
       BxPerN_G = Impl_VGB(Bx_,0:nI+1,0:nJ+1,0:nK+1,implBLK)*InvN_G
       ByPerN_G = Impl_VGB(By_,0:nI+1,0:nJ+1,0:nK+1,implBLK)*InvN_G
       BzPerN_G = Impl_VGB(Bz_,0:nI+1,0:nJ+1,0:nK+1,implBLK)*InvN_G
    end if
  end subroutine impl_init_hall
  !===========================================================================
  subroutine impl_hall_resist

    ! Add partial derivatives of the Hall term to the Jacobian that 
    ! are not calculated by the general algorithm

    ! Cartesian code

    use ModHallResist, ONLY: BxPerN_G, ByPerN_G, BzPerN_G
    !real :: ResistDiag should be Eta0Resist_ND from calc_resistive_flux

    real :: Coeff
    !integer:: iVar
    !-----------------------------------------------------------------------

    ! For diffusion term
    ! R(Bi) = d^2(Bi)/dx^2 + d^2(Bi)/dy^2 + d^2(Bi)/dz^2 
    ! dR(Bx)/dBx, dR(By)/dBy, dR(Bz)/dBz 
    
    !do iDim = 1, nDim
    !   Coeff= - ImplCoeff*ResistDiag/Dxyz(iDim)**2
    !   do k=1,nK; do j=1,nJ; do i=1,nI
    !      do iVar=Bx_, Bz_
    !         JAC(iVar,iVar,i,j,k,1)= &
    !              JAC(iVar,iVar,i,j,k,1)         - 2.0*Coeff
    !         JAC(iVar,iVar,i,j,k,2*nDim)= &
    !              JAC(iVar,iVar,i,j,k,2*nDim)    + Coeff
    !         JAC(iVar,iVar,i,j,k,2*nDim+1)= & 
    !              JAC(iVar,iVar,i,j,k,2*nDim+1)  + Coeff
    !      end do
    !   end do; end do; end do
    !end do

    ! dR(Bx)/dBy
    Coeff = ImplCoeff*wnrm(By_)/wnrm(Bx_)/Dxyz(z_)**2
    ! Main diagonal
    do k=1,nK; do j=1,nJ; do i=1,nI
       JAC(Bx_,By_,i,j,k,1)= JAC(Bx_,By_,i,j,k,1) &
            + Coeff*(BzPerN_G(i,j,k-1)+BzPerN_G(i,j,k+1))
    end do; end do; end do
    !J+1
    do k=1,nK-1; do j=1,nJ; do i=1,nI
       JAC(Bx_,By_,i,j,k,7)=JAC(Bx_,By_,i,j,k,7) - Coeff*BzPerN_G(i,j,k+1)
    end do; end do; end do
    !J-1
    do k=2,nK; do j=1,nJ; do i=1,nI
       JAC(Bx_,By_,i,j,k,6)=JAC(Bx_,By_,i,j,k,6) - Coeff*BzPerN_G(i,j,k-1)
    end do; end do; end do

    ! dR(Bx)/dBz
    Coeff = ImplCoeff*wnrm(Bz_)/wnrm(Bx_)/Dxyz(y_)**2
    ! Main diagonal
    do k=1,nK; do j=1,nJ; do i=1,nI
       JAC(Bx_,Bz_,i,j,k,1)= JAC(Bx_,Bz_,i,j,k,1) &
            - Coeff*(ByPerN_G(i,j-1,k)+ByPerN_G(i,j+1,k))
    end do; end do; end do
    !J+1
    do k=1,nK; do j=1,nJ-1; do i=1,nI
       JAC(Bx_,Bz_,i,j,k,5)= JAC(Bx_,Bz_,i,j,k,5) + Coeff*ByPerN_G(i,j+1,k)
    end do; end do; end do
    !J-1
    do k=1,nK; do j=2,nJ; do i=1,nI
       JAC(Bx_,Bz_,i,j,k,4)= JAC(Bx_,Bz_,i,j,k,4) + Coeff*ByPerN_G(i,j-1,k)
    end do; end do; end do

    ! dR(By)/dBz
    Coeff = ImplCoeff*wnrm(Bz_)/wnrm(By_)/Dxyz(x_)**2
    ! Main diagonal
    do k=1,nK; do j=1,nJ; do i=1,nI
       JAC(By_,Bz_,i,j,k,1)= JAC(By_,Bz_,i,j,k,1) &
            + Coeff*(BxPerN_G(i-1,j,k)+BxPerN_G(i+1,j,k))
    end do; end do; end do
    !I+1
    do k=1,nK; do j=1,nJ; do i=1,nI-1
       JAC(By_,Bz_,i,j,k,3)= JAC(By_,Bz_,i,j,k,3) - Coeff*BxPerN_G(i+1,j,k)
    end do; end do; end do
    !I-1
    do k=1,nK; do j=1,nJ; do i=2,nI
       JAC(By_,Bz_,i,j,k,2)= JAC(By_,Bz_,i,j,k,2) - Coeff*BxPerN_G(i-1,j,k)
    end do; end do; end do

    ! dR(By)/dBx
    Coeff = ImplCoeff*wnrm(Bx_)/wnrm(By_)/Dxyz(z_)**2
    ! Main diagonal
    do k=1,nK; do j=1,nJ; do i=1,nI
       JAC(By_,Bx_,i,j,k,1)= JAC(By_,Bx_,i,j,k,1) &
            - Coeff*(BzPerN_G(i,j,k-1)+BzPerN_G(i,j,k+1))
    end do; end do; end do
    !K+1
    do k=1,nK-1; do j=1,nJ; do i=1,nI
       JAC(By_,Bx_,i,j,k,7)=JAC(By_,Bx_,i,j,k,7) + Coeff*BzPerN_G(i,j,k+1)
    end do; end do; end do
    !K-1
    do k=2,nK; do j=1,nJ; do i=1,nI
       JAC(By_,Bx_,i,j,k,6)=JAC(By_,Bx_,i,j,k,6) + Coeff*BzPerN_G(i,j,k-1)
    end do; end do; end do

    ! dR(Bz)/dBx
    Coeff = ImplCoeff*wnrm(Bx_)/wnrm(Bz_)/Dxyz(y_)**2
    ! Main diagonal
    do k=1,nK; do j=1,nJ; do i=1,nI
       JAC(Bz_,Bx_,i,j,k,1)= JAC(Bz_,Bx_,i,j,k,1) &
            + Coeff*(ByPerN_G(i,j-1,k)+ByPerN_G(i,j+1,k))
    end do; end do; end do
    !J+1
    do k=1,nK; do j=1,nJ-1; do i=1,nI
       JAC(Bz_,Bx_,i,j,k,5)=JAC(Bz_,Bx_,i,j,k,5) - Coeff*ByPerN_G(i,j+1,k)
    end do; end do; end do
    !J-1
    do k=1,nK; do j=2,nJ; do i=1,nI
       JAC(Bz_,Bx_,i,j,k,4)=JAC(Bz_,Bx_,i,j,k,4) - Coeff*ByPerN_G(i,j-1,k)
    end do; end do; end do

    ! dR(Bz)/dBy
    Coeff = ImplCoeff*wnrm(By_)/wnrm(Bz_)/Dxyz(x_)**2
    ! Main diagonal
    do k=1,nK; do j=1,nJ; do i=1,nI
       JAC(Bz_,By_,i,j,k,1)= JAC(Bz_,By_,i,j,k,1) &
            - Coeff*(BxPerN_G(i-1,j,k)+BxPerN_G(i+1,j,k))
    end do; end do; end do
    !I+1
    do k=1,nK; do j=1,nJ; do i=1,nI-1
       JAC(Bz_,By_,i,j,k,3)= JAC(Bz_,By_,i,j,k,3) + Coeff*BxPerN_G(i+1,j,k)
    end do; end do; end do
    !I-1
    do k=1,nK; do j=1,nJ; do i=2,nI
       JAC(Bz_,By_,i,j,k,2)= JAC(Bz_,By_,i,j,k,2) + Coeff*BxPerN_G(i-1,j,k)
    end do; end do; end do

  end subroutine impl_hall_resist

  !============================================================================

  subroutine impl_hall_resist_general

    use ModHallResist, ONLY: set_block_jacobian_cell, DgenDxyz_DDC, &
         BxPerN_G, ByPerN_G, BzPerN_G

    use ModNumConst,   ONLY: iLeviCivita_III

    ! Notation follows the Toth et. al. Hall MHD on Block Adaptive Grids paper

    integer :: iDim, jDim, kDim, lDim, iFace, i, j, k, i2, j2, k2

    integer :: iSub, iSup                  ! indexes of sub and super diagonals
    integer :: jB, lB                      ! index of Bj and Bl
    integer :: iklEpsilon, jklEpsilon      ! Levi-Civita symbols
    
    real :: Coeff, Term, TermSub, TermSup

    real :: BPerN_CD(nI,nJ,nK,3), Area_DDF(3, 3, nI+1, nJ+1, nK+1)
    !-------------------------------------------------------------------------

    call set_block_jacobian_cell(iBlk)

    BPerN_CD(:,:,:,1) = BxPerN_G(1:nI,1:nJ,1:nK)
    BPerN_CD(:,:,:,2) = ByPerN_G(1:nI,1:nJ,1:nK)
    BPerN_CD(:,:,:,3) = BzPerN_G(1:nI,1:nJ,1:nK)

    Area_DDF(:,1,:,1:nJ,1:nK) = FaceAreaI_DFB(:,:,:,:,iBlk)
    Area_DDF(:,2,1:nI,:,1:nK) = FaceAreaJ_DFB(:,:,:,:,iBlk)
    Area_DDF(:,3,1:nI,1:nJ,:) = FaceAreaK_DFB(:,:,:,:,iBlk)

    do kDim = 1,3; do lDim = 1,3
       if(kDim == lDim) CYCLE

       do jDim = 1, 3
          ! Normalization for dR(B_j)/dB_l
          jB = B_ + jDim; 
          lB = B_ + lDim; 
          Coeff = ImplCoeff*wnrm(lB)/wnrm(jB)

          jklEpsilon = iLeviCivita_III(jDim, kDim, lDim)

          do iDim = 1,3
             if(iDim == jDim) CYCLE  ! Terms cancel out

             iklEpsilon = iLeviCivita_III(iDim, kDim, lDim)
             if(iklEpsilon == 0 .and. jklEpsilon == 0) CYCLE

             do iFace = 1,3      ! 3 directions of gen. coordinates

                iSub = 2*iFace   ! stencil index of subdiagonal elements
                iSup = iSub + 1  ! stencil index of superdiagonal elements

                do k=1,nK; do j=1,nJ; do i=1,nI

                   ! Index for the 'right' face of the cell
                   i2 = i + kr(1,iFace)
                   j2 = j + kr(2,iFace)
                   k2 = k + kr(3,iFace)

                   ! Area(iFace)/V*T_ks*(Bi/n*jklEpsilon - Bj/n*iklEpsilon)

                   Term = Coeff *vInv_CB(i,j,k,iBlk) &
                        *DgenDxyz_DDC(kDim,iFace,i,j,k) / dxyz(iFace) &
                        *(BPerN_CD(i,j,k,iDim)*jklEpsilon &
                        - BPerN_CD(i,j,k,jDim)*iklEpsilon)

                   TermSub = Area_DDF(iDim,iFace,i ,j ,k )*Term
                   TermSup = Area_DDF(iDim,iFace,i2,j2,k2)*Term

                   !if(iBlk==1.and.jB==Bz_.and.lB==Bx_.and. &
                   !     i==4.and.j==4.and.k==1.and.iSub==4)then
                   !
                   !if(iBlk==1.and.jB==By_.and.lB==Bz_.and.iFace==1.and.&
                   !     iDim==1.and.kDim==1.and.i==2.and.j==2.and.k==2)then
                   !
                   !   write(*,*)'i2,j2,k2=',i2,j2,k2
                   !   write(*,*)'vInv,Area1,Area2=',vInv_CB(i,j,k,iBlk),&
                   !        Area_DDF(iDim,iFace,i ,j ,k ), &
                   !        Area_DDF(iDim,iFace,i2,j2,k2)
                   !
                   !   write(*,*)'Tks=',DgenDxyz_DDC(kDim,iFace,i,j,k)
                   !   write(*,*)'Coeff,BPerN(iDim),BPerN(jDim)=',&
                   !        Coeff,BPerN_CD(i,j,k,iDim),BPerN_CD(i,j,k,jDim)
                   !
                   !   write(*,*)'jklEpsilon, iklEpsilon=',&
                   !        jklEpsilon, iklEpsilon
                   !
                   !   write(*,*)'TermSub=',TermSub
                   !end if

                   ! Exclude boundaries
                   if(  iFace==1.and.i>1 .or. &
                        iFace==2.and.j>1 .or. &
                        iFace==3.and.k>1)&
                        JAC(jB,lB,i,j,k,iSub) = JAC(jB,lB,i,j,k,iSub) + TermSub

                   if(  iFace==1.and.i<nI .or. &
                        iFace==2.and.j<nJ .or. &
                        iFace==3.and.k<nK)&
                        JAC(jB,lB,i,j,k,iSup) = JAC(jB,lB,i,j,k,iSup) + TermSup

                   ! The main diagonal is -1 * the sum of off-diagonal terms
                   JAC(jB,lB,i,j,k,1) = JAC(jB,lB,i,j,k,1) - TermSub - TermSup

                end do; end do; end do
             end do
          end do
       end do
    end do; end do

  end subroutine impl_hall_resist_general
  !===========================================================================
  subroutine impl_gray_diffusion

    ! Add partial derivatives of the gray diffusion term to the Jacobian that 
    ! are not calculated by the general algorithm

    use ModAdvance,       ONLY: Eradiation_
    use ModGrayDiffusion, ONLY: DiffCoef_VFDB

    integer :: iVar, i, j, k, iDim, Di, Dj, Dk
    real :: Coeff, DiffLeft, DiffRight
    !-----------------------------------------------------------------------
    iVar = Eradiation_

    do iDim = 1, nDim
       Coeff = -ImplCoeff/Dxyz(iDim)**2
       Di = kr(iDim,1)
       Dj = kr(iDim,2)
       Dk = kr(iDim,3)
       do k=1,nK; do j=1,nJ; do i=1,nI
          DiffLeft  = DiffCoef_VFDB(1,i,j,k,iDim,iBlk)
          DiffRight = DiffCoef_VFDB(1,i+Di,j+Dj,k+Dk,iDim,iBlk)
          if(iDim==1.and.i==1 .or. iDim==2.and.j==1 .or. iDim==3.and.k==1)&
               DiffLeft = 0.0
          if(iDim==1.and.i==nI .or. iDim==2.and.j==nJ .or. iDim==3.and.k==nK)&
               DiffRight = 0.0
          JAC(iVar,iVar,i,j,k,1) = JAC(iVar,iVar,i,j,k,1) &
               - Coeff*(DiffLeft + DiffRight)
          JAC(iVar,iVar,i,j,k,2*iDim)   = JAC(iVar,iVar,i,j,k,2*iDim) &
               + Coeff*DiffLeft
          JAC(iVar,iVar,i,j,k,2*iDim+1) = JAC(iVar,iVar,i,j,k,2*iDim+1) &
               + Coeff*DiffRight
       end do; end do; end do
    end do

  end subroutine impl_gray_diffusion

end subroutine impl_jacobian
