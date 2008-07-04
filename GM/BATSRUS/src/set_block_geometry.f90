!^CFG COPYRIGHT UM
subroutine set_root_block_geometry
  use ModProcMH
  use ModSize
  use ModMain, ONLY : TypeBc_I,unusedBLK,nBlock,nBlockMax
  use ModMain,ONLY  : Phi_,Theta_,z_             
  use ModAMR, ONLY : availableBLKs
  use ModGeometry,ONLY: xyzStart_BLK,dx_BLK,dy_BLK,dz_BLK,&
       R2_BLK,&                !^CFG IF SECONDBODY
       TypeGeometry,is_axial_geometry,&          
       R_BLK,x_BLK,y_BLK,z_BLK,Dxyz,XyzMin_D,XyzMax_D
  use ModNumConst
  use ModParallel, ONLY: proc_dims,periodic3D
  use ModMpi
  implicit none

  integer :: i, j, k, iBLK
  real :: dx, dy, dz                         

  ! set the array periodic3D for the periodic boundary

  if(is_axial_geometry())then                     
     Periodic3D(Phi_)=.true.
     if(index(TypeGeometry,'spherical')>0)then
        Periodic3D(Theta_)=.false.
     else
        periodic3D(z_)=any(TypeBc_I(bot_:top_)=='periodic')
     end if
  else                                          
     periodic3D(1)=any(TypeBc_I(east_:west_)=='periodic')
     periodic3D(2)=any(TypeBc_I(south_:north_)=='periodic')
     periodic3D(3)=any(TypeBc_I(bot_:top_)=='periodic')
  end if                                                
 
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
       UseExtraBoundary,DoFixExtraBoundaryOrPole,unusedBLK,ProcTest,BlkTest   
  use ModMain, ONLY : UseBody2                       !^CFG IF SECONDBODY
  use ModNodes
  use ModGeometry
  use ModNumConst
  use ModPhysics, ONLY : Rbody
  use ModPhysics, ONLY : xBody2,yBody2,zBody2 !^CFG IF SECONDBODY
  use ModParallel, ONLY : periodic3D
  use ModBoundaryCells
  implicit none

  integer, intent(in) :: iBLK

  integer :: i,j,k, iBoundary
  real :: dx,dy,dz,fAx,fAy,fAz,cV,VInv
  real,dimension(nDim)::XyzOfNode111_D               

  logical:: DoTest, DoTestMe
  character(len=*), parameter :: NameSub='fix_block_geometry'
  !---------------------------------------------------------------------------

  if(iBlk==BlkTest .and. iProc==ProcTest)then
     call set_oktest(NameSub, DoTest, DoTestMe)
  else
     DoTest = .false.; DoTestMe = .false.
  end if

  !\--------------------------------
  ! Assign cell centered coordinates.
  !/--------------------------------


  if(.not.UseCovariant)then                  
     dx = dx_BLK(iBLK)
     dy = dy_BLK(iBLK)
     dz = dz_BLK(iBLK)

     ! This case is the normal one using a spherical Radius.
     do k = 1-gcn, nK+gcn
        do j = 1-gcn, nJ+gcn
           do i = 1-gcn, nI+gcn
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
     cV = dx*dy*dz                                  

     fAx = dy*dz       !Cell face areas
     fAy = dz*dx 
     fAz = dx*dy 

     fAx_BLK(iBLK) = fAx                           
     fAy_BLK(iBLK) = fAy
     fAz_BLK(iBLK) = fAz
     cV_BLK(iBLK) = cV                               
     vInv_CB(:,:,:,iBLK) = cOne/cV                

     !Calculate the node coordinates
     do i=1,1+nI; do j=1,1+nJ; do k=1,1+nK
        NodeX_NB(i,j,k,iBLK) = 0.125*sum(x_BLK(i-1:i,j-1:j,k-1:k,iBLK))
        NodeY_NB(i,j,k,iBLK) = 0.125*sum(y_BLK(i-1:i,j-1:j,k-1:k,iBLK))
        NodeZ_NB(i,j,k,iBLK) = 0.125*sum(z_BLK(i-1:i,j-1:j,k-1:k,iBLK))
     end do; end do; end do                  
  else                                      
     !Cell center coordinates are calculated directly as the
     !transformed generalized coordinates
     call gen_to_xyz_arr(XyzStart_BLK(:,iBLK),&
                         dx_BLK(iBLK),dy_BLK(iBLK),dz_BLK(iBLK),&
                         1-gcn,nI+gcn,1-gcn,nJ+gcn,1-gcn,nK+gcn,&
                         x_BLK(:,:,:,iBLK),&
                         y_BLK(:,:,:,iBLK),&     
                         z_BLK(:,:,:,iBLK))

     R_BLK(:,:,:,iBLK)=sqrt(&
          x_BLK(:,:,:,iBLK)**2+y_BLK(:,:,:,iBLK)**2+z_BLK(:,:,:,iBLK)**2)

     if(UseVertexBasedGrid)then
        
        !Node coordinates are calculated directly, not derived from
        !the cell centered x,y,z_BLK. In general case there is no close 
        !relation between the node coordinates and the "cell center" 
        !coordinates

        
        XyzOfNode111_D(1)=XyzStart_BLK(1,iBLK)-dx_BLK(iBLK)*cHalf
        XyzOfNode111_D(2)=XyzStart_BLK(2,iBLK)-dy_BLK(iBLK)*cHalf
        XyzOfNode111_D(3)=XyzStart_BLK(3,iBLK)-dz_BLK(iBLK)*cHalf

        
        call gen_to_xyz_arr(XyzOfNode111_D,&
                            dx_BLK(iBLK),dy_BLK(iBLK),dz_BLK(iBLK),&
                            1,1+nI,1,1+nJ,1,1+nK,&
                            NodeX_NB(:,:,:,iBLK),&
                            NodeY_NB(:,:,:,iBLK),&
                            NodeZ_NB(:,:,:,iBLK))


        !
        !Mark the block for fixing the covariant geometry afterwards
        !(when the refinement level of the neighboring blocks are
        !all known). The condition for fixing the block is
        !any(OldLevel_IIIB(:,:,:,iBLK)/=NeiLev_BLK(:,:,:,iBLK).and.
        !(OldLevel_IIIB(:,:,:,iBLK)==-1.or.NeiLev_BLK(:,:,:,iBLK).  
        !To force the fix procedure for the given block set
        
        OldLevel_IIIB(:,:,:,iBLK)=0

     else

        !Node coordinates are calculated as the for the point
        !of the faces intersection
 
        call calc_node_coords_covar

     end if
     !Face area vectors and cell volumes are expressed 
     !in terms of the node coordinates
     
     call fix_covariant_geometry(iBLK)
   
  end if                                          

  Rmin_BLK(iBLK)  = minval(R_BLK(:,:,:,iBLK))

  if (.not. UseBody2) then                        !^CFG IF SECONDBODY BEGIN
     ! if no second body set R2 to zero
     R2_BLK(:,:,:,iBLK) = cZero
  else                    
     ! calculate the radius as measured from the second body
     R2_BLK(:,:,:,iBLK) = sqrt( &
          (x_BLK(:,:,:,iBLK)-xBody2)**2 + &
          (y_BLK(:,:,:,iBLK)-yBody2)**2 + &
          (z_BLK(:,:,:,iBLK)-zBody2)**2)
  end if                                          

 
  Rmin2_BLK(iBLK) = minval(R2_BLK(:,:,:,iBLK))    !^CFG END SECONDBODY


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


  if(DoTestMe)then
     write(*,*)NameSub,': far_field_bcs_blk=',far_field_bcs_BLK(iBlk)
     write(*,*)NameSub,': xyzStart_BLK=',xyzStart_BLK(:,BlkTest)
     write(*,*)NameSub,': dx,dy,dz= ',dx_BLK(BlkTest),dy_BLK(iBLK),dz_BLK(iBLK)
     write(*,*)NameSub,': xyzmin=',XyzMin_D(:)
     write(*,*)NameSub,': xyzmax=',XyzMax_D(:)
  end if

  !\
  ! TRUE_CELL: if not inside a body or outside the outer face boundary
  !/
  true_cell(:,:,:,iBLK)=.true.
  IsBoundaryBlock_IB(:,iBLK)=.false.
  do iBoundary = MinBoundary, MaxBoundary
     IsBoundaryBlock_IB(iBoundary,iBLK) = .true.
  end do
  IsBoundaryBlock_IB(ExtraBc_,iBLK) = &
       (DoFixExtraBoundaryOrPole.and..not.is_axial_geometry()) &             
       .or.UseExtraBoundary

  ! set true_cell array
  call set_boundary_cells(iBLK)

  do iBoundary = MinBoundary, min(MaxBoundary,Body1_)
     IsBoundaryBlock_IB(iBoundary,iBLK)=any(IsBoundaryCell_GI(:,:,:,iBoundary))
     true_cell(:,:,:,iBLK) = &
          true_cell(:,:,:,iBLK) .and. .not.IsBoundaryCell_GI(:,:,:,iBoundary)
  end do

  BodyFlg_B(iBLK) = .not. all(true_cell(:,:,:,iBLK))
                                                    
  if(.not.is_axial_geometry())     &              
       DoFixExtraBoundary_B(iBLK) = DoFixExtraBoundaryOrPole &
       .and.any(IsBoundaryCell_GI(:,:,:,ExtraBc_))
                                                    
  IsBoundaryCell_GI(:,:,:,ExtraBc_) = &
       UseExtraBoundary .and. IsBoundaryCell_GI(:,:,:,ExtraBc_)

  do iBoundary=ExtraBc_,MaxBoundary
     if(SaveBoundaryCells.and.iBoundary>=MinBoundarySaved)then
        IsBoundaryCell_IGB(iBoundary,:,:,:,iBLK)=&
             IsBoundaryCell_GI(:,:,:,iBoundary)
        true_cell(1:nI,1:nJ,1:nK,iBLK) = true_cell(1:nI,1:nJ,1:nK,iBLK) &
          .and. .not.IsBoundaryCell_GI(1:nI,1:nJ,1:nK,iBoundary)
     else
        true_cell(:,:,:,iBLK) = true_cell(:,:,:,iBLK) &
             .and. .not.IsBoundaryCell_GI(:,:,:,iBoundary)
        IsBoundaryBlock_IB(iBoundary,iBLK)=&
             any(IsBoundaryCell_GI(:,:,:,iBoundary))
     end if
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

       NodeX_NB(i+1,j+1,k+1,iBLK) = det(A1_DD)*DetInv

       A1_DD(:,1)=A_DD(:,1)
       A1_DD(:,3)=A_DD(:,3)
       A1_DD(:,2)=B_D

       NodeY_NB(i+1,j+1,k+1,iBLK) = det(A1_DD)*DetInv
       
       A1_DD(:,1:2)=A_DD(:,1:2)
       A1_DD(:,3)=B_D

       NodeZ_NB(i+1,j+1,k+1,iBLK) = det(A1_DD)*DetInv
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

subroutine set_boundary_cells(iBLK)

  use ModProcMH
  use ModMain
  use ModPhysics,  ONLY: Rbody
  use ModGeometry, ONLY: R_BLK, IsBoundaryBlock_IB, IsBoundaryCell_GI, &
       MinBoundary,MaxBoundary
  use ModPhysics,  ONLY: Rbody2                            !^CFG IF SECONDBODY
  use ModGeometry, ONLY: R2_BLK                            !^CFG IF SECONDBODY
  use ModGeometry,ONLY:x1,x2,y1,y2,z1,z2,x_BLK,y_BLK,z_BLK 
  use ModUser, ONLY: user_set_boundary_cells

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


  if(IsBoundaryBlock_IB(ExtraBc_,iBLK))&             
       call user_set_boundary_cells(iBLK)



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
 

end subroutine set_boundary_cells
