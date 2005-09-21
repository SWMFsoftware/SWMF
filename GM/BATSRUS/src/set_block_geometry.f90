!^CFG COPYRIGHT UM
subroutine set_root_block_geometry
  use ModProcMH
  use ModSize
  use ModMain, ONLY : TypeBc_I,unusedBLK,nBlock,nBlockMax
  use ModAMR, ONLY : availableBLKs
  use ModGeometry,ONLY: xyzStart_BLK,dx_BLK,dy_BLK,dz_BLK,&
       R2_BLK,&                !^CFG IF SECONDBODY
       TypeGeometry,&          !^CFG IF NOT CARTESIAN
       R_BLK,x_BLK,y_BLK,z_BLK,Dxyz,XyzMin_D,XyzMax_D
  use ModNumConst
  use ModParallel, ONLY: proc_dims,periodic3D
  use ModMpi
  implicit none

  integer :: i, j, k, iBLK
  real :: dx, dy, dz                         

  ! set the array periodic3D for the periodic boundary

  select case(TypeGeometry)       !^CFG IF NOT CARTESIAN
     case('cartesian')            !^CFG IF NOT CARTESIAN
        periodic3D(1)=any(TypeBc_I(east_:west_)=='periodic')
        periodic3D(2)=any(TypeBc_I(south_:north_)=='periodic')
        periodic3D(3)=any(TypeBc_I(bot_:top_)=='periodic')
     case('spherical','spherical_lnr')        !^CFG IF NOT CARTESIAN BEGIN
        periodic3D(1)=.false. 
        periodic3D(2)=.true.  
        periodic3D(3)=.false.
     case('cylindrical')
        periodic3D(1)=.false. 
        periodic3D(2)=.true.  
        periodic3D(3)= any(TypeBc_I(bot_:top_)=='periodic')
     case default
        call stop_mpi('Unknown TypeGeometry = '//TypeGeometry)
     end select                   !^CFG END CARTESIAN

  xyzStart_BLK = -777777.
  dx_BLK = -777777.
  dy_BLK = -777777.
  dz_BLK = -777777.
  x_BLK  = -777777.
  y_BLK  = -777777.
  z_BLK  = -777777.
  R_BLK  = -777777.  
  R2_BLK = -777777.           !^CFG IF SECONDBODY

  ! Set the block geometry for all root blocks
  if (iProc == 0) then
     iBLK = 0
     Dxyz(:)=(XyzMax_D(:)-XyzMin_D(:))/(nCells(:)*proc_dims(:))

     do k = 1, proc_dims(3)
        do j = 1, proc_dims(2)
           do i = 1, proc_dims(1)
              iBLK = iBLK+1

              unusedBLK(iBLK) = .false.

              dx_BLK(iBLK) = Dxyz(1)  
              dy_BLK(iBLK) = Dxyz(2)
              dz_BLK(iBLK) = Dxyz(3)
              xyzStart_BLK(1,iBLK) = Dxyz(1) * (cHalf + nI*(i-1))+XyzMin_D(1)
              xyzStart_BLK(2,iBLK) = Dxyz(2) * (cHalf + nJ*(j-1))+XyzMin_D(2)
              xyzStart_BLK(3,iBLK) = Dxyz(3) * (cHalf + nK*(k-1))+XyzMin_D(3)
              call fix_block_geometry(iBLK)
           end do
        end do
     end do
  end if

  ! Let every PE know the available blocks on PE 0
  nBlock    = product(proc_dims)
  nBlockMax = nBlock
  availableBLKs(0,0)=nBlockMax+1

end subroutine set_root_block_geometry
!==============================================================================
subroutine fix_block_geometry(iBLK)

  use ModMain, ONLY : body1,body1_,body2_,ExtraBc_,&
       UseExtraBoundary,DoFixExtraBoundary,unusedBLK                       
  use ModMain, ONLY : UseBody2                       !^CFG IF SECONDBODY
  use ModNodes
  use ModGeometry
  use ModNumConst
  use ModPhysics, ONLY : Rbody
  use ModPhysics, ONLY : xBody2,yBody2,zBody2 !^CFG IF SECONDBODY
  use ModParallel, ONLY : periodic3D
  implicit none

  integer, intent(in) :: iBLK

  integer :: i,j,k, iBoundary
  real :: dx,dy,dz,fAx,fAy,fAz,cV,VInv
  !---------------------------------------------------------------------------
  dx = dx_BLK(iBLK)
  dy = dy_BLK(iBLK)
  dz = dz_BLK(iBLK)

  select case(TypeGeometry)                        !^CFG IF NOT CARTESIAN
  case('cartesian')                                !^CFG IF NOT CARTESIAN
     !\
     ! Assign cell centered coordinates.
     !/
     ! This case is the normal one using a spherical Radius.
     do i = 1-gcn, nI+gcn
        do j = 1-gcn, nJ+gcn
           do k = 1-gcn, nK+gcn
              x_BLK(i,j,k,iBLK) =  (i-1)*dx + xyzStart_BLK(1,iBLK)
              y_BLK(i,j,k,iBLK) =  (j-1)*dy + xyzStart_BLK(2,iBLK)
              z_BLK(i,j,k,iBLK) =  (k-1)*dz + xyzStart_BLK(3,iBLK)     
              R_BLK(i,j,k,iBLK) = sqrt( &
                   x_BLK(i,j,k,iBLK)**2+   &
                   y_BLK(i,j,k,iBLK)**2+   &
                   z_BLK(i,j,k,iBLK)**2)
           end do
        end do
     end do
     cV = dx*dy*dz                                  !^CFG IF CARTESIAN BEGIN

     fAx = dy*dz       !Cell face areas
     fAy = dz*dx 
     fAz = dx*dy 


     fAx_BLK(iBLK) = fAx                                       
     fAy_BLK(iBLK) = fAy
     fAz_BLK(iBLK) = fAz
     cV_BLK(iBLK) = cV                                            
     VInv_CB(:,:,:,iBLK) = cOne/cV                !^CFG END CARTESIAN
!     call set_covar_cartesian_geometry(iBLK)     !^CFG IF NOT CARTESIAN 
  case('spherical')                               !^CFG IF NOT CARTESIAN 
!     call set_covar_spherical_geometry(iBLK)     !^CFG IF NOT CARTESIAN 
  case('spherical_lnr')                           !^CFG IF NOT CARTESIAN 
!     call set_covar_spherical2_geometry(iBLK)    !^CFG IF NOT CARTESIAN 
  case('cylindrical')                             !^CFG IF NOT CARTESIAN 
!     call set_covar_cylindrical_geometry(iBLK)   !^CFG IF NOT CARTESIAN
  case default                                    !^CFG IF NOT CARTESIAN BEGIN
     call stop_mpi('Unknown TypeGeometry = '//TypeGeometry)
  end select                                      !^CFG END CARTESIAN 

  if (.not. UseBody2) then                        !^CFG IF SECONDBODY BEGIN
     ! if no second body set R2 to zero
     R2_BLK(:,:,:,iBLK) = cZero
  else                    
     ! calculate the radius as measured from the second body
     R2_BLK(:,:,:,iBLK) = sqrt( &
          (x_BLK(:,:,:,iBLK)-xBody2)**2 + &
          (y_BLK(:,:,:,iBLK)-yBody2)**2 + &
          (z_BLK(:,:,:,iBLK)-zBody2)**2)
  end if                                          !^CFG END SECONDBODY

  Rmin_BLK(iBLK)  = minval(R_BLK(:,:,:,iBLK))
  Rmin2_BLK(iBLK) = minval(R2_BLK(:,:,:,iBLK))    !^CFG IF SECONDBODY

  ! Compute node geometry
  !^CFG IF CARTESIAN BEGIN
  do i=0,nI; do j=0,nJ; do k=0,nK
     NodeX_IIIB(i,j,k,iBLK) = cEighth*sum(x_BLK(i:i+1,j:j+1,k:k+1,iBLK))
     NodeY_IIIB(i,j,k,iBLK) = cEighth*sum(y_BLK(i:i+1,j:j+1,k:k+1,iBLK))
     NodeZ_IIIB(i,j,k,iBLK) = cEighth*sum(z_BLK(i:i+1,j:j+1,k:k+1,iBLK))
  end do; end do; end do
  !^CFG END CARTESIAN
  !call calc_node_coords_covar            !^CFG IF NOT CARTESIAN

  far_field_BCs_BLK(iBLK) = &
       (((xyzStart_BLK(1,iBLK)-dx_BLK(iBLK))<XyzMin_D(1).or.&
       (xyzStart_BLK(1,iBLK)+nI*dx_BLK(iBLK))>XyzMax_D(1)) &
       .and. .not.periodic3D(1)) .or. &
       (((xyzStart_BLK(2,iBLK)-dy_BLK(iBLK))<XyzMin_D(2).or.&
       (xyzStart_BLK(2,iBLK)+nJ*dy_BLK(iBLK))>XyzMax_D(2)) &
       .and. .not.periodic3D(2)).or. &
       (((xyzStart_BLK(3,iBLK)-dz_BLK(iBLK))<XyzMin_D(3).or.&
       (xyzStart_BLK(3,iBLK)+nK*dz_BLK(iBLK))>XyzMax_D(3)) &
       .and. .not.periodic3D(3))  

  !\
  ! TRUE_CELL: if not inside a body or outside the outer face boundary
  !/
  true_cell(:,:,:,iBLK)=.true.
  IsBoundaryBlock_IB(:,iBLK)=.false.
  do iBoundary = MinBoundary, MaxBoundary
     IsBoundaryBlock_IB(iBoundary,iBLK) = .true.
  end do
  IsBoundaryBlock_IB(ExtraBc_,iBLK) = &
       DoFixExtraBoundary &                         !^CFG IF FACEOUTERBC BEGIN
       .and.TypeGeometry=='cartesian' &             !^CFG IF NOT CARTESIAN
       .or. &                                       !^CFG END FACEOUTERBC
       UseExtraBoundary

  ! set true_cell array
  call set_boundary_cells(iBLK)

  do iBoundary = MinBoundary, min(MaxBoundary,Body1_)
     IsBoundaryBlock_IB(iBoundary,iBLK)=any(IsBoundaryCell_GI(:,:,:,iBoundary))
     true_cell(:,:,:,iBLK) = &
          true_cell(:,:,:,iBLK) .and. .not.IsBoundaryCell_GI(:,:,:,iBoundary)
  end do

  BodyFlg_B(iBLK) = .not. all(true_cell(:,:,:,iBLK))
                                                    !^CFG IF FACEOUTERBC BEGIN
  if(TypeGeometry=='cartesian')     &               !^CFG IF NOT CARTESIAN
       DoFixExtraBoundary_B(iBLK) = DoFixExtraBoundary &
       .and.any(IsBoundaryCell_GI(:,:,:,ExtraBc_))
                                                    !^CFG END FACEOUTERBC
  IsBoundaryCell_GI(:,:,:,ExtraBc_) = &
       UseExtraBoundary .and. IsBoundaryCell_GI(:,:,:,ExtraBc_)

  do iBoundary=ExtraBc_,MaxBoundary
     IsBoundaryBlock_IB(iBoundary,iBLK)=any(IsBoundaryCell_GI(:,:,:,iBoundary))
     true_cell(:,:,:,iBLK) = true_cell(:,:,:,iBLK) &
          .and. .not.IsBoundaryCell_GI(:,:,:,iBoundary)
  end do
  
  BodyFlg_B(iBLK)= BodyFlg_B(iBLK) .and. any(true_cell(:,:,:,iBLK))

  !\
  ! body_BLK: if any cell INCLUDING ghost cells is inside body(ies)
  !/
  body_BLK(iBLK) = .not. all(true_cell(:,:,:,iBLK))
  !\
  ! TRUE_BLK: if all cells EXCLUDING ghost cells are outside body(ies)
  !/
  true_BLK(iBLK) = all(true_cell(1:nI,1:nJ,1:nK,iBLK))

contains
  !===========================================================================
  subroutine calc_node_coords_covar
    implicit none
    real,dimension(3,3):: A_DD, A1_DD
    real,dimension(3)  :: B_D
    real               :: DetInv
    !------------------------------------------------------------------------
    do i=0,nI; do j=0,nJ; do k=0,nK

       A_DD(1,1)=x_BLK(i+1,j,k,iBLK)-x_BLK(i,j,k,iBLK)
       A_DD(1,2)=y_BLK(i+1,j,k,iBLK)-y_BLK(i,j,k,iBLK)
       A_DD(1,3)=z_BLK(i+1,j,k,iBLK)-z_BLK(i,j,k,iBLK)

       A_DD(2,1)=x_BLK(i,j+1,k,iBLK)-x_BLK(i,j,k,iBLK)
       A_DD(2,2)=y_BLK(i,j+1,k,iBLK)-y_BLK(i,j,k,iBLK)
       A_DD(2,3)=z_BLK(i,j+1,k,iBLK)-z_BLK(i,j,k,iBLK)

       A_DD(3,1)=x_BLK(i,j,k+1,iBLK)-x_BLK(i,j,k,iBLK)
       A_DD(3,2)=y_BLK(i,j,k+1,iBLK)-y_BLK(i,j,k,iBLK)
       A_DD(3,3)=z_BLK(i,j,k+1,iBLK)-z_BLK(i,j,k,iBLK)

       DetInv=cOne/det(A_DD)

       B_D(1)=cHalf*(R_BLK(i+1,j,k,iBLK)-R_BLK(i,j,k,iBLK))*&
            (R_BLK(i+1,j,k,iBLK)+R_BLK(i,j,k,iBLK))
       B_D(2)=cHalf*(R_BLK(i,j+1,k,iBLK)-R_BLK(i,j,k,iBLK))*&
            (R_BLK(i,j+1,k,iBLK)+R_BLK(i,j,k,iBLK))
       B_D(3)=cHalf*(R_BLK(i,j,k+1,iBLK)-R_BLK(i,j,k,iBLK))*&
            (R_BLK(i,j,k+1,iBLK)+R_BLK(i,j,k,iBLK))

       A1_DD(:,2:3)=A_DD(:,2:3)
       A1_DD(:,1)=B_D

       NodeX_IIIB(i,j,k,iBLK) = det(A1_DD)*DetInv

       A1_DD(:,1)=A_DD(:,1)
       A1_DD(:,3)=A_DD(:,3)
       A1_DD(:,2)=B_D

       NodeY_IIIB(i,j,k,iBLK) = det(A1_DD)*DetInv
       
       A1_DD(:,1:2)=A_DD(:,1:2)
       A1_DD(:,3)=B_D

       NodeZ_IIIB(i,j,k,iBLK) = det(A1_DD)*DetInv
    end do; end do; end do
    
  end subroutine calc_node_coords_covar
  !===========================================================================
  real function det(A_DD)
    implicit none
    real,dimension(3,3),intent(in)::A_DD
    det=A_DD(1,1)*(A_DD(2,2)*A_DD(3,3)-&
         A_DD(3,2)*A_DD(2,3))-&
         A_DD(1,2)*(A_DD(2,1)*A_DD(3,3)-&
         A_DD(2,3)*A_DD(3,1))+&
         A_DD(1,3)*(A_DD(2,1)*A_DD(3,2)-&
         A_DD(2,2)*A_DD(3,1))
  end function det

end subroutine fix_block_geometry
!=============================================================================
subroutine set_xyzminmax_cart
  use ModGeometry, ONLY:x1,x2,y1,y2,z1,z2,XyzMin_D,XyzMax_D
  use ModMain,ONLY:x_,y_,z_
  implicit none

  XyzMin_D(x_) = x1
  XyzMin_D(y_) = y1
  XyzMin_D(z_) = z1
  XyzMax_D(x_) = x2
  XyzMax_D(y_) = y2
  XyzMax_D(z_) = z2
end subroutine set_xyzminmax_cart

!^CFG IF NOT CARTESIAN BEGIN
!=============================================================================
subroutine set_xyzminmax_sph
  use ModGeometry, ONLY: x1,x2,y1,y2,z1,z2,XyzMin_D,XyzMax_D
  use ModNumConst, ONLY: cZero,cTwoPI,cPi
  use ModMain, ONLY: R_,Phi_,Theta_
  implicit none

  XyzMin_D(R_)     = cZero
  XyzMin_D(Phi_)   = cZero
  XyzMin_D(Theta_) = cZero
  XyzMax_D(R_)     = sqrt(max(x1*x1,x2*x2)+max(y1*y1,y2*y2)+max(z1*z1,z2*z2))
  XyzMax_D(Phi_)   = cTwoPi
  XyzMax_D(Theta_) = cPi
end subroutine set_xyzminmax_sph
!=============================================================================
subroutine set_xyzminmax_sph2
  use ModGeometry, ONLY: x1,x2,y1,y2,z1,z2,XyzMin_D,XyzMax_D
  use ModNumConst, ONLY: cZero,cTwoPI,cPi,cHalf
  use ModMain, ONLY: R_,Phi_,Theta_
  implicit none

  XyzMin_D(R_)     = cZero
  XyzMin_D(Phi_)   = cZero
  XyzMin_D(Theta_) = cZero
  XyzMax_D(R_)     = cHalf*alog(&
       max(x1*x1,x2*x2)+max(y1*y1,y2*y2)+max(z1*z1,z2*z2))
  XyzMax_D(Phi_)   = cTwoPi
  XyzMax_D(Theta_) = cPi
end subroutine set_xyzminmax_sph2
!=============================================================================
subroutine set_xyzminmax_cyl
  use ModGeometry, ONLY:x1,x2,y1,y2,z1,z2,XyzMin_D,XyzMax_D
  use ModNumConst, ONLY:cZero,cTwoPI
  use ModMain, ONLY:R_,Phi_,z_
  implicit none

  XyzMin_D(R_)     = cZero
  XyzMin_D(Phi_)   = cZero
  XyzMin_D(z_)     = z1
  XyzMax_D(R_)     = sqrt(max(x1*x1,x2*x2)+max(y1*y1,y2*y2))
  XyzMax_D(Phi_)   = cTwoPi
  XyzMax_D(z_)     = z2

end subroutine set_xyzminmax_cyl
!^CFG END CARTESIAN
!==============================================================================
subroutine set_boundary_cells(iBLK)

  use ModProcMH
  use ModMain
  use ModPhysics,  ONLY: Rbody
  use ModGeometry, ONLY: R_BLK, IsBoundaryBlock_IB, IsBoundaryCell_GI, &
       MinBoundary,MaxBoundary
  use ModPhysics,  ONLY: Rbody2                            !^CFG IF SECONDBODY
  use ModGeometry, ONLY: R2_BLK                            !^CFG IF SECONDBODY
  use ModGeometry,ONLY:x1,x2,y1,y2,z1,z2,x_BLK,y_BLK,z_BLK !^CFG IF FACEOUTERBC

  implicit none
  integer,intent(in)::iBLK
  !----------------------------------------------------------------------------
 
  IsBoundaryCell_GI=.false.  
  !^CFG IF SECONDBODY BEGIN
  if(IsBoundaryBlock_IB(Body2_,iBLK))&               
       IsBoundaryCell_GI(:,:,:,Body2_) = &
       UseBody2 .and. R2_BLK(:,:,:,iBLK) < RBody2  
  !^CFG END SECONDBODY

  if(IsBoundaryBlock_IB(Body1_,iBLK)) &
       IsBoundaryCell_GI(:,:,:,Body1_) = &
       body1    .and. R_BLK(:,:,:,iBLK) < Rbody

  !^CFG IF USERFILES BEGIN
  if(IsBoundaryBlock_IB(ExtraBc_,iBLK))&             
       call set_extra_boundary_cells(iBLK)
  !^CFG END USERFILES

  !^CFG IF FACEOUTERBC BEGIN
  if(IsBoundaryBlock_IB(East_,iBLK)) &
       IsBoundaryCell_GI(:,:,:,East_)=x_BLK(:,:,:,iBLK)<x1  
  if(IsBoundaryBlock_IB(West_,iBLK)) &
       IsBoundaryCell_GI(:,:,:,West_)=x_BLK(:,:,:,iBLK)>x2 
  if(IsBoundaryBlock_IB(South_,iBLK)) &
       IsBoundaryCell_GI(:,:,:,South_)= y_BLK(:,:,:,iBLK)<y1
  if(IsBoundaryBlock_IB(North_,iBLK)) &
       IsBoundaryCell_GI(:,:,:,North_)= y_BLK(:,:,:,iBLK)>y2 
  if(IsBoundaryBlock_IB(Bot_,iBLK)) &
       IsBoundaryCell_GI(:,:,:,Bot_)= z_BLK(:,:,:,iBLK)<z1
  if(IsBoundaryBlock_IB(Top_,iBLK)) &
       IsBoundaryCell_GI(:,:,:,Top_)= z_BLK(:,:,:,iBLK)>z2  
  !^CFG END FACEOUTERBC

end subroutine set_boundary_cells
