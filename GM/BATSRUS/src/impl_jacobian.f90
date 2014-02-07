!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

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
  use ModNumConst, ONLY: i_DD
  use ModvarIndexes
  use ModAdvance, ONLY: time_BLK
  use ModB0, ONLY: B0_DX, B0_DY, B0_DZ, set_b0_face
  use ModImplicit
  use ModHallResist, ONLY: UseHallResist, hall_factor
  use ModRadDiffusion, ONLY: add_jacobian_rad_diff
  use ModResistivity, ONLY: UseResistivity, add_jacobian_resistivity
  use ModGeometry, ONLY: true_cell
  use BATL_lib, ONLY: IsCartesianGrid, IsRzGeometry, &
       FaceNormal_DDFB, CellSize_DB, CellVolume_GB

  implicit none

  integer, intent(in) :: implBLK
  real,    intent(out):: JAC(nw,nw,nI,nJ,nK,nstencil)

  integer :: iBLK
  real, dimension(nw,nI,nJ,nK)                :: Impl_VC      , ImplEps_VC
  real, dimension(nI+1,nJ+1,nK+1)             :: dfdwLface, dfdwRface
  real, dimension(nw,nI,nJ,nK)                :: s_VC, sEps_VC, sPowell_VC
  real :: DivB(nI,nJ,nK)
  real :: B0_DFD(MaxDim,nI+1,nJ+1,nK+1,MaxDim), Cmax_DF(MaxDim,nI+1,nJ+1,nK+1)

  real   :: qeps, coeff
  logical:: divbsrc, UseDivbSource0
  integer:: i,j,k,i1,i2,i3,j1,j2,j3,k1,k2,k3,istencil,iw,jw,idim,qj

  logical :: oktest, oktest_me

  real :: Dxyz(MaxDim)
  real :: FluxLeft_VFD(nW,nI+1,nJ+1,nK+1,MaxDim) ! Unperturbed left flux
  real :: FluxRight_VFD(nW,nI,nJ,nK,MaxDim)      ! Unperturbed right flux
  real :: FluxEpsLeft_VF(nW,nI+1,nJ+1,nK+1)    ! Perturbed left flux
  real :: FluxEpsRight_VF(nW,nI,nJ,nK)         ! Perturbed right flux
  real :: FaceArea_F(nI, nJ, nK)               ! Only the inner faces

  real :: HallFactor_G(0:nI+1,j0_:nJp1_,k0_:nKp1_)
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
  Dxyz = CellSize_DB(:,iBlk)
  if(UseB0)then
     call set_b0_face(iBLK)
     B0_DFD(:,1:nI+1,1:nJ  ,1:nK  ,x_)=B0_DX(:,1:nI+1,1:nJ,1:nK)
     B0_DFD(:,1:nI  ,1:nJ+1,1:nK  ,y_)=B0_DY(:,1:nI,1:nJ+1,1:nK)
     B0_DFD(:,1:nI  ,1:nJ  ,1:nK+1,z_)=B0_DZ(:,1:nI,1:nJ,1:nK+1)
  else
     B0_DFD =0.0
  end if

  if(UseHallResist)call impl_init_hall

  ! Initialize matrix to zero (to be safe)
  JAC=0.0

  ! Initialize reference flux and the cmax array
  do iDim = 1, nDim

     i1 = 1+i_DD(1,idim); j1= 1+i_DD(2,idim); k1= 1+i_DD(3,idim)
     i2 =nI+i_DD(1,idim); j2=nJ+i_DD(2,idim); k2=nK+i_DD(3,idim);

     call get_face_flux(Impl_VC,B0_DFD(:,i1:i2,j1:j2,k1:k2,idim),&
          nI,nJ,nK,iDim,iBLK,FluxLeft_VFD(:,i1:i2,j1:j2,k1:k2,iDim))

     call get_face_flux(Impl_VC,B0_DFD(:,1:nI,1:nJ,1:nK,iDim),&
          nI,nJ,nK,iDim,iBLK,FluxRight_VFD(:,:,:,:,iDim))

     ! Average w for each cell interface into ImplEps_VC
     i1 = 1-i_DD(1,idim); j1= 1-i_DD(2,idim); k1= 1-i_DD(3,idim)

     ! Calculate orthogonal cmax for each interface in advance
     call get_cmax_face(                            &
          0.5*(Impl_VGB(1:nw, 1:i2, 1:j2, 1:k2,implBLK)+ &
          Impl_VGB(1:nw,i1:nI,j1:nJ,k1:nK,implBLK)),     &
          B0_DFD(:,1:i2,1:j2,1:k2,idim),       &
          i2,j2,k2,idim,iBlk,Cmax_DF(idim,1:i2,1:j2,1:k2))

     ! cmax always occurs as -ImplCoeff*0.5/dx*cmax
     coeff = -0.5 
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

     do iDim = 1, nDim
        ! Index limits for faces and shifted centers
        i1 = 1+i_DD(1,idim); j1= 1+i_DD(2,idim); k1= 1+i_DD(3,idim)
        i2 =nI+i_DD(1,idim); j2=nJ+i_DD(2,idim); k2=nK+i_DD(3,idim);
        i3 =nI-i_DD(1,idim); j3=nJ-i_DD(2,idim); k3=nK-i_DD(3,idim);

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
           coeff=-0.5/qeps/wnrm(iw) 

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
                /CellVolume_GB(1:nI, 1:nJ, 1:nK, iBlk)
           dfdwRface( 1:nI, 1:nJ, 1:nK) = dfdwRface( 1:nI, 1:nJ, 1:nK) &
                /CellVolume_GB(1:nI, 1:nJ, 1:nK, iBlk)

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
                .or.  (iw == E_ .and. UseImplicitEnergy) &
                ) )then
              if(.not.IsCartesianGrid .and. jw>=Bx_ .and. jw<=B_+nDim)then
                 ! The source terms are always multiplied by coeff
                 coeff=-0.5*wnrm(jw)/wnrm(iw)
                 ! Get the corresponding face area
                 FaceArea_F(i1:nI,j1:nJ,k1:nK) = &
                      FaceNormal_DDFB(jw-B_,iDim,i1:nI,j1:nJ,k1:nK,iBLK)

                 ! Relative to the right face flux Q is shifted to the left
                 dfdwLface(i1:nI,j1:nJ,k1:nK)=dfdwLface(i1:nI,j1:nJ,k1:nK)+ &
                      coeff*sPowell_VC(iw,i1:nI,j1:nJ,k1:nK) &
                      *FaceArea_F(i1:nI,j1:nJ,k1:nK) &
                      /CellVolume_GB(i1:nI,j1:nJ,k1:nK,iBlk)

                 dfdwRface(i1:nI,j1:nJ,k1:nK)=dfdwRface(i1:nI,j1:nJ,k1:nK)+ &
                      coeff*sPowell_VC(iw, 1:i3, 1:j3, 1:k3) &
                      *FaceArea_F(i1:nI,j1:nJ,k1:nK) &
                      /CellVolume_GB(1:i3,1:j3,1:k3,iBlk)

              elseif(jw==B_+idim)then
                 ! The source terms are always multiplied by coeff
                 coeff=-0.5/dxyz(idim)*wnrm(jw)/wnrm(iw)

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
           coeff=-1.0/qeps/wnrm(iw)
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

  ! Add extra terms for (Hall) resistivity
  if(UseResistivity .or. UseHallResist) &
       call add_jacobian_resistivity(iBlk, JAC)

  ! Add extra terms for radiative diffusion
  if(UseRadDiffusion) &
       call add_jacobian_rad_diff(iBLK, JAC)

  ! Multiply JAC by the implicit timestep dt and ImplCoeff
  if(time_accurate)then
     do k=1,nK; do j=1,nJ; do i=1,nI
        if(true_cell(i,j,k,iBLK))then
           JAC(:,:,i,j,k,:) = JAC(:,:,i,j,k,:)*dt*ImplCoeff
        else
           ! Set JAC = 0.0 inside body
           JAC(:,:,i,j,k,:) = 0.0
        end if
     end do; end do; end do
  else
     ! Local time stepping has time_BLK=0.0 inside the body
     do istencil=1,nstencil; do jw=1,nw; do iw=1,nw
        JAC(iw,jw,:,:,:,istencil)=JAC(iw,jw,:,:,:,istencil) &
             *time_BLK(1:nI,1:nJ,1:nK,iBLK)*implCFL*ImplCoeff
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
    if(IsCartesianGrid)then
       if(IsRzGeometry)call CON_stop('impl_divbsrc_init not working for RZ')

       divb = &
            ( Impl_VGB(Bx_,2:nI+1,1:nJ,1:nK,implBLK)           &
            - Impl_VGB(Bx_,0:nI-1,1:nJ,1:nK,implBLK))/dxyz(x_)
       if(nJ>1) divb = divb &
            +(Impl_VGB(By_,1:nI,2:nJ+1,1:nK,implBLK)           &
            - Impl_VGB(By_,1:nI,0:nJ-1,1:nK,implBLK))/dxyz(y_)
       if(nK>1) divb = divb &
            +(Impl_VGB(Bz_,1:nI,1:nJ,2:nK+1,implBLK)           &
            - Impl_VGB(Bz_,1:nI,1:nJ,0:nK-1,implBLK))/dxyz(z_)
       divb=0.5*divb

    else
       do k=1,nK; do j=1,nJ; do i=1,nI
          divb(i,j,k) = &
               sum (Impl_VGB(Bx_:B_+nDim,i+1,j,k,implBLK) &
               *    FaceNormal_DDFB(:,1,i+1,j,k,iBlk))&
               -sum(Impl_VGB(Bx_:B_+nDim,i-1,j,k,implBLK) &
               *    FaceNormal_DDFB(:,1,i,j,k,iBlk))  &
               +sum(Impl_VGB(Bx_:B_+nDim,i,j+1,k,implBLK) &
               *    FaceNormal_DDFB(:,2,i,j+1,k,iBlk))&
               -sum(Impl_VGB(Bx_:B_+nDim,i,j-1,k,implBLK) &
               *    FaceNormal_DDFB(:,2,i,j,k,iBlk))

          if(nK>1) divb(i,j,k) = divb(i,j,k) &
               +sum(Impl_VGB(Bx_:B_+nDim,i,j,k+1,implBLK) &
               *    FaceNormal_DDFB(:,3,i,j,k+1,iBlk))&
               -sum(Impl_VGB(Bx_:B_+nDim,i,j,k-1,implBLK) &
               *    FaceNormal_DDFB(:,3,i,j,k,iBlk))

          divb(i,j,k) = 0.5/CellVolume_GB(i,j,k,iBlk)*divb(i,j,k)

       end do; end do; end do
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

       if(.not. UseImplicitEnergy) CYCLE
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
            +wnrm(Bx_)/wnrm(rhoUx_)*divb(i,j,k) 
       JAC(rhoUy_,By_,i,j,k,1)=JAC(rhoUy_,By_,i,j,k,1)&
            +wnrm(By_)/wnrm(rhoUy_)*divb(i,j,k) 
       JAC(rhoUz_,Bz_,i,j,k,1)=JAC(rhoUz_,Bz_,i,j,k,1)&
            +wnrm(Bz_)/wnrm(rhoUz_)*divb(i,j,k) 

       ! Q(B)= -divB*rhoU/rho
       ! dQ(B)/drho = +divB*rhoU/rho**2
       JAC(Bx_,rho_,i,j,k,1)=JAC(Bx_,rho_,i,j,k,1) &
            -wnrm(rho_)/wnrm(Bx_)*divb(i,j,k) &
            *Impl_VC(rhoUx_,i,j,k)/Impl_VC(rho_,i,j,k)**2
       JAC(By_,rho_,i,j,k,1)=JAC(By_,rho_,i,j,k,1) &
            -wnrm(rho_)/wnrm(By_)*divb(i,j,k) &
            *Impl_VC(rhoUy_,i,j,k)/Impl_VC(rho_,i,j,k)**2
       JAC(Bz_,rho_,i,j,k,1)=JAC(Bz_,rho_,i,j,k,1) &
            -wnrm(rho_)/wnrm(Bz_)*divb(i,j,k) &
            *Impl_VC(rhoUz_,i,j,k)/Impl_VC(rho_,i,j,k)**2

       ! dQ(B)/drhoU= -divB/rho
       JAC(Bx_,rhoUx_,i,j,k,1)=JAC(Bx_,rhoUx_,i,j,k,1)&
            +wnrm(rhoUx_)/wnrm(Bx_)*divb(i,j,k)/Impl_VC(rho_,i,j,k) 
       JAC(By_,rhoUy_,i,j,k,1)=JAC(By_,rhoUy_,i,j,k,1)&
            +wnrm(rhoUy_)/wnrm(By_)*divb(i,j,k)/Impl_VC(rho_,i,j,k) 
       JAC(Bz_,rhoUz_,i,j,k,1)=JAC(Bz_,rhoUz_,i,j,k,1)&
            +wnrm(rhoUz_)/wnrm(Bz_)*divb(i,j,k)/Impl_VC(rho_,i,j,k) 

       if(.not.UseImplicitEnergy) CYCLE

       ! Q(E)= -divB*rhoU.B/rho
       ! dQ(E)/drho = +divB*rhoU.B/rho**2
       JAC(E_,rho_,i,j,k,1)=JAC(E_,rho_,i,j,k,1)&
            -wnrm(rho_)/wnrm(E_)*divb(i,j,k)*&
            (Impl_VC(rhoUx_,i,j,k)*Impl_VC(Bx_,i,j,k)&
            +Impl_VC(rhoUy_,i,j,k)*Impl_VC(By_,i,j,k)&
            +Impl_VC(rhoUz_,i,j,k)*Impl_VC(Bz_,i,j,k))&
            /Impl_VC(rho_,i,j,k)**2

       ! dQ(E)/drhoU = -divB*B/rho
       JAC(E_,rhoUx_,i,j,k,1)=JAC(E_,rhoUx_,i,j,k,1) &
            +wnrm(rhoUx_)/wnrm(E_)*divb(i,j,k) &
            *Impl_VC(Bx_,i,j,k)/Impl_VC(rho_,i,j,k) 
       JAC(E_,rhoUy_,i,j,k,1)=JAC(E_,rhoUy_,i,j,k,1) &
            +wnrm(rhoUy_)/wnrm(E_)*divb(i,j,k) &
            *Impl_VC(By_,i,j,k)/Impl_VC(rho_,i,j,k) 
       JAC(E_,rhoUz_,i,j,k,1)=JAC(E_,rhoUz_,i,j,k,1) &
            +wnrm(rhoUz_)/wnrm(E_)*divb(i,j,k) &
            *Impl_VC(Bz_,i,j,k)/Impl_VC(rho_,i,j,k) 

       ! dQ(E)/dB = -divB*rhoU/rho
       JAC(E_,Bx_,i,j,k,1)=JAC(E_,Bx_,i,j,k,1) &
            +wnrm(Bx_)/wnrm(E_)*divb(i,j,k) &
            *Impl_VC(rhoUx_,i,j,k)/Impl_VC(rho_,i,j,k) 
       JAC(E_,By_,i,j,k,1)=JAC(E_,By_,i,j,k,1) &
            +wnrm(By_)/wnrm(E_)*divb(i,j,k) &
            *Impl_VC(rhoUy_,i,j,k)/Impl_VC(rho_,i,j,k) 
       JAC(E_,Bz_,i,j,k,1)=JAC(E_,Bz_,i,j,k,1) &
            +wnrm(Bz_)/wnrm(E_)*divb(i,j,k) &
            *Impl_VC(rhoUz_,i,j,k)/Impl_VC(rho_,i,j,k) 
    end do; end do; end do

  end subroutine impl_divbsrc_middle

  !===========================================================================
  subroutine impl_init_hall

    ! Calculate cell centered currents to be used by getflux

    use ModHallResist, ONLY: HallJ_CD, IonMassPerCharge_G, &
         set_ion_mass_per_charge

    use ModGeometry, ONLY: DgenDxyz_DDC, set_block_jacobian_cell 

    real :: DbDgen_DD(3,3)                     

    real :: InvDx2, InvDy2, InvDz2

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub='impl_init_hall'
    !----------------------------------------------------------------------
    if(iProc == PROCtest.and.implBLK==implBLKtest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if

    call set_ion_mass_per_charge(iBlk)

    InvDx2 = 0.5/Dxyz(x_); InvDy2 = 0.5/Dxyz(y_); InvDz2 = 0.5/Dxyz(z_)

    if(IsCartesianGrid)then

       do k=1,nK; do j=1,nJ; do i=1,nI
          ! Jx = dBz/dy - dBy/dz
          if(nJ>1) HallJ_CD(i,j,k,x_) =                    &
               +InvDy2*(Impl_VGB(Bz_,i,j+1,k,implBLK)      &
               -        Impl_VGB(Bz_,i,j-1,k,implBLK))
          if(nK>1) HallJ_CD(i,j,k,x_) = HallJ_CD(i,j,k,x_) &
               -InvDz2*(Impl_VGB(By_,i,j,k+1,implBLK)      &
               -        Impl_VGB(By_,i,j,k-1,implBLK))
       end do; end do; end do

       do k=1,nK; do j=1,nJ; do i=1,nI
          ! Jy = dBx/dz - dBz/dx
          HallJ_CD(i,j,k,y_) = &
               -InvDx2*(Impl_VGB(Bz_,i+1,j,k,implBLK)      &
               -        Impl_VGB(Bz_,i-1,j,k,implBLK))
          if(nK>1) HallJ_CD(i,j,k,y_) = HallJ_CD(i,j,k,y_) &
               +InvDz2*(Impl_VGB(Bx_,i,j,k+1,implBLK)      &
               -        Impl_VGB(Bx_,i,j,k-1,implBLK))
       end do; end do; end do

       do k=1,nK; do j=1,nJ; do i=1,nI
          ! Jz = dBy/dx - dBx/dy
          HallJ_CD(i,j,k,z_) = &
               +InvDx2*(Impl_VGB(By_,i+1,j,k,implBLK)      &
               -        Impl_VGB(By_,i-1,j,k,implBLK))
          if(nJ>1) HallJ_CD(i,j,k,z_) = HallJ_CD(i,j,k,z_) &
               -InvDy2*(Impl_VGB(Bx_,i,j+1,k,implBLK)      &
               -        Impl_VGB(Bx_,i,j-1,k,implBLK))
       end do; end do; end do

    else                                        

       call set_block_jacobian_cell(iBlk)

       DbDgen_DD = 0.0 !!! make it MaxDim*nDim and use Dim1_, Dim2_, Dim3_

       do k=1,nK; do j=1,nJ; do i=1,nI
          DbDgen_DD(:,1) = InvDx2*&
               (Impl_VGB(Bx_:Bz_,i+1,j,k,implBLK) &
               -Impl_VGB(Bx_:Bz_,i-1,j,k,implBLK))
          if(nJ>1) DbDgen_DD(:,2) = InvDy2* &
               (Impl_VGB(Bx_:Bz_,i,j+1,k,implBLK) &
               -Impl_VGB(Bx_:Bz_,i,j-1,k,implBLK))
          if(nK>1) DbDgen_DD(:,3) = InvDz2* &
               (Impl_VGB(Bx_:Bz_,i,j,k+1,implBLK) &
               -Impl_VGB(Bx_:Bz_,i,j,k-1,implBLK))

          ! Jx = dBz/dy - dBy/dz
          if(nJ>1) HallJ_CD(i,j,k,x_) = &
               + sum(DbDgen_DD(z_,:)*DgenDxyz_DDC(:,y_,i,j,k)) 
          if(nK>1) HallJ_CD(i,j,k,x_) = HallJ_CD(i,j,k,x_) &
               - sum(DbDgen_DD(y_,:)*DgenDxyz_DDC(:,z_,i,j,k))

          ! Jy = dBx/dz - dBz/dx
          HallJ_CD(i,j,k,y_) = &
               - sum(DbDgen_DD(z_,:)*DgenDxyz_DDC(:,x_,i,j,k))
          if(nK>1)HallJ_CD(i,j,k,y_) = HallJ_CD(i,j,k,y_) &
               + sum(DbDgen_DD(x_,:)*DgenDxyz_DDC(:,z_,i,j,k))

          ! Jz = dBy/dx - dBx/dy
          HallJ_CD(i,j,k,z_) = &
               + sum(DbDgen_DD(y_,:)*DgenDxyz_DDC(:,x_,i,j,k)) 
          if(nJ>1) HallJ_CD(i,j,k,z_) = HallJ_CD(i,j,k,z_) &
               - sum(DbDgen_DD(x_,:)*DgenDxyz_DDC(:,y_,i,j,k))

       end do; end do; end do

    end if                                    
    if(DoTestMe) write(*,*) NameSub,' HallJ_CD=',HallJ_CD(iTest,jTest,kTest,:)

    do k = k0_,nKp1_; do j=j0_,nJp1_; do i=0,nI+1
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

  end subroutine impl_init_hall

end subroutine impl_jacobian
