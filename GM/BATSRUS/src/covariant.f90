!^CFG COPYRIGHT UM
!^CFG FILE COVARIANT
!subroutine gen_to_xyz_arr maps the piece of an equally spaced grid
!in the space of GENERALIZED COORDINATES to the cartesian space xyz
subroutine gen_to_xyz_arr(&
     GenCoord111_D,    &!(in)Gen.coords.of the point (1,1,1)
     dGen1,dGen2,dGen3,&!(in)Mesh sizes in Gen coords.
     iStart,iMax,      &!(in)The first and the last value of i index
     jStart,jMax,      &!(in)The first and the last value of j index
     kStart,kMax,      &!(in)The first and the last value of k index
     X_C,              &!(out)Cartesian x-coord. of the mapped points
     Y_C,              &!(out)Cartesian y-coord. of the mapped points    
     Z_C)               !(out)Cartesian z-coord. of the mapped points
  use ModNumConst
  use ModCovariant,ONLY:TypeGeometry
  use ModMain,     ONLY:nDim,R_,Phi_,Theta_,x_,y_,z_
  use ModUser
  implicit none
  real,intent(in):: GenCoord111_D(nDim),dGen1,dGen2,dGen3
  integer,intent(in)::  iStart,iMax,jStart,jMax,kStart,kMax
  real,dimension(iStart:iMax,jStart:jMax,kStart:kMax),&
       intent(out)::X_C,Y_C,Z_C
!----------------------------------------------------------
  integer::i,j,k
  real::R,Theta,Phi,sinTheta,cosTheta,sinPhi,cosPhi
!----------------------------------------------------------
  
  
  select case(TypeGeometry)
  case('cartesian')
     !Gen1=x , Gen2=y, Gen3=z
     do k = kStart, kMax
        do j = jStart, jMax
           do i = iStart, iMax
              X_C(i,j,k) =  (i-1)*dGen1 + GenCoord111_D(x_)
              Y_C(i,j,k) =  (j-1)*dGen2 + GenCoord111_D(y_)
              Z_C(i,j,k) =  (k-1)*dGen3 + GenCoord111_D(z_)
           end do
        end do
     end do
  case('spherical')
     ! Gen1=R, Gen2=Phi, Gen3=Theta
     do k = kStart, kMax
        Theta      =  (k-1)*dGen3 + GenCoord111_D(Theta_)
        sinTheta   =  sin(Theta)
        cosTheta   =  cos(Theta)
        do j = jStart, jMax
           Phi     =  (j-1)*dGen2 + GenCoord111_D(Phi_)
           sinPhi  =  sin(Phi)
           cosPhi  =  cos(Phi)
           do i = iStart, iMax
              R    =  (i-1)*dGen1 + GenCoord111_D(R_)
             
                        
              X_C(i,j,k) = R*sinTheta*cosPhi                      
              Y_C(i,j,k) = R*sinTheta*sinPhi                      
              Z_C(i,j,k) = R*cosTheta
           end do
        end do
     end do
  case('spherical_lnr')
     ! Gen1=log(R), Gen2=Phi, Gen3=Theta
     do k = kStart, kMax
        Theta      =      (k-1)*dGen3 + GenCoord111_D(Theta_)
        sinTheta   =  sin(Theta)
        cosTheta   =  cos(Theta)
        do j = jStart, jMax
           Phi     =      (j-1)*dGen2 + GenCoord111_D(Phi_)
           sinPhi  =  sin(Phi)
           cosPhi  =  cos(Phi)
           do i = iStart, iMax
              R    =   exp((i-1)*dGen1 + GenCoord111_D(R_))
             
              X_C(i,j,k) = R*sinTheta*cosPhi                      
              Y_C(i,j,k) = R*sinTheta*sinPhi                      
              Z_C(i,j,k) = R*cosTheta
           end do
        end do
     end do
  case('cylindrical')
    ! Gen1=r, Gen2=Phi, Gen3=z
     do k = kStart, kMax
        Z_C(:,:,k) = (k-1)*dGen3 + GenCoord111_D(z_)
        do j = jStart, jMax
           Phi     =  (j-1)*dGen2 + GenCoord111_D(Phi_)
           sinPhi  =  sin(Phi)
           cosPhi  =  cos(Phi)
           do i = iStart, iMax
              R    =       (i-1)*dGen1 + GenCoord111_D(R_)
             
              X_C(i,j,k) = R*cosPhi                      
              Y_C(i,j,k) = R*sinPhi                      
           end do
        end do
     end do     
  case default
     call stop_mpi('Unknown geometry: '//TypeGeometry)
  end select
end subroutine gen_to_xyz_arr
!-------------------------------------------------------------!
subroutine fix_covariant_geometry(iBLK)
  use ModCovariant
  use ModNodes,ONLY:NodeX_NB,NodeY_NB,NodeZ_NB
  use ModGeometry,ONLY: vInv_CB
  use ModMain,ONLY:x_,y_,z_
  implicit none
  integer,intent(in)::iBLK
  real,dimension(nDim,1:nI+1,1:nJ+1,1:nK+1)::XyzNode_DN
  !It is easy to see that volume can be represented as follows:
  !\int{dV}=\int{{\bf r}\cdot d{\bf S}/nDim. The following array store
  !the dot products of face area vectors by the raduis vector of
  !the "face center"
  real,dimension(1:nI+1,nJ,nK)::RDotFaceAreaI_F 
  real,dimension(1:nI,nJ+1,nK)::RDotFaceAreaJ_F
  real,dimension(1:nI,nJ,nK+1)::RDotFaceAreaK_F
  
  integer::i,j,k

  do k=1,nK+1; do j=1,nJ+1; do i=1,nI+1
     XyzNode_DN(x_,i,j,k)=NodeX_NB(i,j,k,iBLK)
     XyzNode_DN(y_,i,j,k)=NodeY_NB(i,j,k,iBLK)
     XyzNode_DN(z_,i,j,k)=NodeZ_NB(i,j,k,iBLK)
  end do; end do; end do
  !\
  ! Face area vector and its dot product by the face 
  ! center radius-vector. FACE I
  !/
  call get_face_area_i(XyzNode_DN,&
                       1,nI+1,1,nJ,1,nK,&
                       FaceAreaI_DFB(:,:,:,:,iBLK))

  do k=1,nK; do j=1,nJ; do i=1,nI+1
     RDotFaceAreaI_F(i,j,k)=cQuarter*dot_product(&
          XyzNode_DN(:,i,j  ,k  )+ &
          XyzNode_DN(:,i,j+1,k  )+ &
          XyzNode_DN(:,i,j+1,k+1)+ &
          XyzNode_DN(:,i,j  ,k+1), &
          FaceAreaI_DFB(:,i,j,k,iBLK))
  end do; end do; end do

  !\
  ! Face area vector and its dot product by the face 
  ! center radius-vector. FACE J
  !/
  call get_face_area_j(XyzNode_DN,&
                       1,nI,1,nJ+1,1,nK,&
                       FaceAreaJ_DFB(:,:,:,:,iBLK))

  do k=1,nK; do j=1,nJ+1; do i=1,nI
     RDotFaceAreaJ_F(i,j,k)=cQuarter*dot_product(&
          XyzNode_DN(:,i  ,j,k  )+ &
          XyzNode_DN(:,i  ,j,k+1)+ &
          XyzNode_DN(:,i+1,j,k+1)+ &
          XyzNode_DN(:,i+1,j,k  ), &
          FaceAreaJ_DFB(:,i,j,k,iBLK))
  end do; end do; end do 

  !\
  ! Face area vector and its dot product by the face 
  ! center radius-vector. FACE K
  !/
  call get_face_area_k(XyzNode_DN,&
                       1,nI,1,nJ,1,nK+1,&
                       FaceAreaK_DFB(:,:,:,:,iBLK))

  do k=1,nK+1; do j=1,nJ; do i=1,nI
     RDotFaceAreaK_F(i,j,k)=cQuarter*dot_product(&
          XyzNode_DN(:,i  ,j  ,k)+ &
          XyzNode_DN(:,i+1,j  ,k)+ &
          XyzNode_DN(:,i+1,j+1,k)+ &
          XyzNode_DN(:,i  ,j+1,k), &
          FaceAreaK_DFB(:,i,j,k,iBLK))
  end do; end do; end do 
  
  !Calculate Volume (inverse)
  vInv_CB(:,:,:,iBLK)=nDim/(&
       RDotFaceAreaI_F(2:nI+1,:,:)-RDotFaceAreaI_F(1:nI,:,:)+&
       RDotFaceAreaJ_F(:,2:nJ+1,:)-RDotFaceAreaJ_F(:,1:nJ,:)+&
       RDotFaceAreaK_F(:,:,2:nK+1)-RDotFaceAreaK_F(:,:,1:nK) )

  select case(TypeGeometry)
  case('cartesian')                               
     call fix_cartesian_geometry(iBLK)      
  case('spherical')                               
     call fix_spherical_geometry(iBLK)      
  case('spherical_lnr')                           
     call fix_spherical2_geometry(iBLK)     
  case('cylindrical')                             
     call fix_cylindrical_geometry(iBLK)   
  case default                                    
     call stop_mpi('Unknown TypeGeometry = '//TypeGeometry)
  end select
end subroutine fix_covariant_geometry
!-------------------------------------------------------------------
subroutine fix_cartesian_geometry(iBLK)
  use ModCovariant
  use ModGeometry
  use ModNumConst
  implicit none
  integer, intent(in) :: iBLK

  FaceArea2MinI_B(iBLK)=cZero
  FaceArea2MinJ_B(iBLK)=cZero
  FaceArea2MinK_B(iBLK)=cZero
end subroutine fix_cartesian_geometry
!---------------------------------------------------------------------
subroutine fix_spherical_geometry(iBLK)
  use ModMain
  use ModGeometry,ONLY:x_BLK,y_BLK,z_BLK,dx_BLK,dy_BLK,dz_BLK,&
       DoFixExtraBoundary_B,XyzStart_BLK,XyzMin_D,XyzMax_D,vInv_CB
  use ModCovariant
  use ModNumConst
  implicit none

  integer, intent(in) :: iBLK
  
  real::dR,dPhi,dTheta
 
  dR=dx_BLK(iBLK)
  dPhi=dy_BLK(iBLK)
  dTheta=dz_BLK(iBLK)



  DoFixExtraBoundary_B(iBLK) = XyzStart_BLK(Theta_,iBLK)-dz_BLK(iBLK)&
       <XyzMin_D(Theta_).or.&
       XyzStart_BLK(Theta_,iBLK)+nK*dz_BLK(iBLK)>XyzMax_D(Theta_)

  FaceArea2MinI_B(iBLK)=cZero
  FaceArea2MinJ_B(iBLK)=cZero
  if(DoFixExtraBoundary_B(iBLK)) then
     FaceArea2MinK_B(iBLK)=dR*XyzStart_BLK(1,iBLK)*&
          (cTwo*tan(cHalf*dPhi))*&       
          (cTwo*tan(cHalf*dTheta))
     FaceArea2MinK_B(iBLK)=FaceArea2MinK_B(iBLK)**2
  else
     FaceArea2MinK_B(iBLK)=cZero
  end if
end subroutine fix_spherical_geometry

!-----------------------------------------------------------------
!real function SpherFaceAreaK(i,k,RStart,ThetaStart,dR,dPhi,dTheta)
!  use ModNumConst
!  implicit none
!  integer, intent(in)::i,k
!  real,intent(in)::RStart,ThetaStart,dR,dPhi,dTheta
!  SpherFaceAreaK= &                              
!       dR*sin(ThetaStart+(k-1-cHalf)*dTheta)*&               
!       (cTwo*tan(cHalf*dPhi))/(cos(cHalf*dTheta)*sin(dTheta))
!end function SpherFaceAreaK

!---------------------------------------------------------------------
subroutine fix_spherical2_geometry(iBLK)
  use ModMain
  use ModGeometry,ONLY:x_BLK,y_BLK,z_BLK,dx_BLK,dy_BLK,dz_BLK,&
       DoFixExtraBoundary_B,XyzStart_BLK,XyzMin_D,XyzMax_D,vInv_CB
  use ModCovariant
  use ModNumConst
  implicit none

  integer, intent(in) :: iBLK
  real::dR2,dPhi,dTheta

  dR2=dx_BLK(iBLK)
  dPhi=dy_BLK(iBLK)
  dTheta=dz_BLK(iBLK)

 

  DoFixExtraBoundary_B(iBLK) = XyzStart_BLK(Theta_,iBLK)-dz_BLK(iBLK)&
       <XyzMin_D(Theta_).or.&
       XyzStart_BLK(Theta_,iBLK)+nK*dz_BLK(iBLK)>XyzMax_D(Theta_)
  
  FaceArea2MinI_B(iBLK)=cZero
  FaceArea2MinJ_B(iBLK)=cZero
  if(DoFixExtraBoundary_B(iBLK)) then
     FaceArea2MinK_B(iBLK)=exp(cTwo*XyzStart_BLK(1,iBLK))*&
          sinh(dR2)*(cOne+cosh(dR2))*cHalf*&
          (cTwo*tan(cHalf*dPhi))*&       
          (cTwo*tan(cHalf*dTheta))
     FaceArea2MinK_B(iBLK)=FaceArea2MinK_B(iBLK)**2
  else
     FaceArea2MinK_B(iBLK)=cZero
  end if

end subroutine fix_spherical2_geometry


!-----------------------------------------------------------------
subroutine fix_cylindrical_geometry(iBLK)
  use ModMain
  use ModGeometry
  implicit none
  integer,intent(in)::iBLK
  DoFixExtraBoundary_B(iBLK)= XyzStart_BLK(R_,iBLK)-dx_BLK(iBLK)<XyzMin_D(R_)
  FaceArea2MinI_B(iBLK)=cZero
  FaceArea2MinJ_B(iBLK)=cZero
  FaceArea2MinK_B(iBLK)=cZero
end subroutine fix_cylindrical_geometry
!-------------------------------------------------------------------

subroutine calc_b0source_covar(iBlock)  
  use ModProcMH  
  use ModMain,ONLY:UseB0Source,x_,y_,z_,R_,Theta_,Phi_
  use ModSize
  use ModParallel, ONLY :&
       neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth
  use ModCovariant
  use ModAdvance, ONLY : &
       B0xFace_x_BLK,B0yFace_x_BLK,B0zFace_x_BLK, &
       B0xFace_y_BLK,B0yFace_y_BLK,B0zFace_y_BLK, &
       B0xFace_z_BLK,B0yFace_z_BLK,B0zFace_z_BLK, &
       B0SourceMatrix_DDCB
  use ModGeometry,ONLY: dx_BLK,dy_BLK,dz_BLK,XyzStart_BLK,TypeGeometry,&
       vInv_CB
  use ModNumConst
  implicit none

  integer, intent(in) :: iBlock
  integer::i,j,k,iFace,jFace,kFace,iDirB0,iDirFA,iSide
  integer::i2,j2
  real::divB0,x,y,z,R,Phi,Theta
  real,dimension(3)::GenCoord111_D,dGenFine_D
  real,dimension(ndim,0:1,0:1,East_:Top_)::FaceArea_DIIS, B0_DIIS
  real,dimension(nDim,0:1,0:1,0:1)   :: RefB0_DIII
  real,dimension(nDim,-1:2,-1:2,-1:2):: RefXyz_DIII
  real,dimension(nDim,0:2,0:2,0:2)   :: RefXyzNodes_DIII
  real,dimension(nDim,0:1,0:1,0:1)   :: RefFaceArea_DIII

  dGenFine_D(1)=cHalf*dx_BLK(iBlock)
  dGenFine_D(2)=cHalf*dy_BLK(iBlock)
  dGenFine_D(3)=cHalf*dz_BLK(iBlock)

  do k=1,nK
     do j=1,nJ
        do i=1,nI
           FaceArea_DIIS=cZero
           B0_DIIS=cZero

           if(i==nI.and.neiLWest(iBlock)==-1)then
              call correct_b0_face(West_)
              B0xFace_x_BLK(i+1,j,k,iBlock)=&
                   sum(B0_DIIS(x_,:,:,West_))*cQuarter
              B0yFace_x_BLK(i+1,j,k,iBlock)=&
                   sum(B0_DIIS(y_,:,:,West_))*cQuarter
              B0zFace_x_BLK(i+1,j,k,iBlock)=&
                   sum(B0_DIIS(z_,:,:,West_))*cQuarter
           else
              FaceArea_DIIS(:,0,0,West_)=FaceAreaI_DFB(:,i+1,j,k,iBlock)
              B0_DIIS(x_,0,0,West_)= B0xFace_x_BLK(i+1,j,k,iBlock)
              B0_DIIS(y_,0,0,West_)= B0yFace_x_BLK(i+1,j,k,iBlock)
              B0_DIIS(z_,0,0,West_)= B0zFace_x_BLK(i+1,j,k,iBlock)
           end if

           if(i==1.and.neiLEast(iBlock)==-1)then
              call correct_b0_face(East_)
              B0xFace_x_BLK(i,j,k,iBlock)=sum(B0_DIIS(x_,:,:,East_))*cQuarter
              B0yFace_x_BLK(i,j,k,iBlock)=sum(B0_DIIS(y_,:,:,East_))*cQuarter
              B0zFace_x_BLK(i,j,k,iBlock)=sum(B0_DIIS(z_,:,:,East_))*cQuarter
           else
              FaceArea_DIIS(:,0,0,East_)=FaceAreaI_DFB(:,i,j,k,iBlock)
              B0_DIIS(x_,0,0,East_)= B0xFace_x_BLK(i,j,k,iBlock)
              B0_DIIS(y_,0,0,East_)= B0yFace_x_BLK(i,j,k,iBlock)
              B0_DIIS(z_,0,0,East_)= B0zFace_x_BLK(i,j,k,iBlock)
           end if

           If(j==nJ.and.neiLNorth(iBlock)==-1)then
              call correct_b0_face(North_)
              B0xFace_y_BLK(i,j+1,k,iBlock)=&
                   sum(B0_DIIS(x_,:,:,North_))*cQuarter
              B0yFace_y_BLK(i,j+1,k,iBlock)=&
                   sum(B0_DIIS(y_,:,:,North_))*cQuarter
              B0zFace_y_BLK(i,j+1,k,iBlock)=&
                   sum(B0_DIIS(z_,:,:,North_))*cQuarter
           else
              FaceArea_DIIS(:,0,0,North_)=FaceAreaJ_DFB(:,i,j+1,k,iBlock) 
              B0_DIIS(x_,0,0,North_)= B0xFace_y_BLK(i,j+1,k,iBlock)
              B0_DIIS(y_,0,0,North_)= B0yFace_y_BLK(i,j+1,k,iBlock)
              B0_DIIS(z_,0,0,North_)= B0zFace_y_BLK(i,j+1,k,iBlock)
           end If

           if(j==1.and.neiLsouth(iBlock)==-1)then
              call correct_b0_face(South_)
              B0xFace_y_BLK(i,j,k,iBlock)=sum(B0_DIIS(x_,:,:,South_))*cQuarter
              B0yFace_y_BLK(i,j,k,iBlock)=sum(B0_DIIS(y_,:,:,South_))*cQuarter
              B0zFace_y_BLK(i,j,k,iBlock)=sum(B0_DIIS(z_,:,:,South_))*cQuarter
           else
              FaceArea_DIIS(:,0,0,South_)=FaceAreaJ_DFB(:,i,j,k,iBlock)   
              B0_DIIS(x_,0,0,South_)= B0xFace_y_BLK(i,j,k,iBlock)
              B0_DIIS(y_,0,0,South_)= B0yFace_y_BLK(i,j,k,iBlock)
              B0_DIIS(z_,0,0,South_)= B0zFace_y_BLK(i,j,k,iBlock)
           end if

           if(k==nK.and.neiLTop(iBlock)==-1)then
              call correct_b0_face(Top_)
              B0xFace_z_BLK(i,j,k+1,iBlock)=&
                   sum(B0_DIIS(x_,:,:,Top_))*cQuarter
              B0yFace_z_BLK(i,j,k+1,iBlock)=&
                   sum(B0_DIIS(y_,:,:,Top_))*cQuarter
              B0zFace_z_BLK(i,j,k+1,iBlock)=&
                   sum(B0_DIIS(z_,:,:,Top_))*cQuarter
           else
              FaceArea_DIIS(:,0,0,Top_)=FaceAreaK_DFB(:,i,j,k+1,iBlock)   
              B0_DIIS(x_,0,0,Top_)= B0xFace_z_BLK(i,j,k+1,iBlock)
              B0_DIIS(y_,0,0,Top_)= B0yFace_z_BLK(i,j,k+1,iBlock)
              B0_DIIS(z_,0,0,Top_)= B0zFace_z_BLK(i,j,k+1,iBlock)
           end if

           if(k==1.and.neiLBot(iBlock)==-1)then
              call correct_b0_face(Bot_)
              B0xFace_z_BLK(i,j,k,iBlock)=&
                   sum(B0_DIIS(x_,:,:,Bot_))*cQuarter
              B0yFace_z_BLK(i,j,k,iBlock)=&
                   sum(B0_DIIS(y_,:,:,Bot_))*cQuarter
              B0zFace_z_BLK(i,j,k,iBlock)=&
                   sum(B0_DIIS(z_,:,:,Bot_))*cQuarter
           else
              FaceArea_DIIS(:,0,0,Bot_)=FaceAreaK_DFB(:,i,j,k,iBlock)  
              B0_DIIS(x_,0,0,Bot_)= B0xFace_z_BLK(i,j,k,iBlock)
              B0_DIIS(y_,0,0,Bot_)= B0yFace_z_BLK(i,j,k,iBlock)
              B0_DIIS(z_,0,0,Bot_)= B0zFace_z_BLK(i,j,k,iBlock)
           end if

           if(UseB0Source)then
              DivB0 = cZero
              B0SourceMatrix_DDCB(:,:,i,j,k,iBlock)=cZero

              do iSide=East_,Bot_,2
                 FaceArea_DIIS(:,:,:,iSide)=-FaceArea_DIIS(:,:,:,iSide)
              end do

              do iSide=East_,Top_
                 do j2=0,1
                    do i2=0,1
                       DivB0= DivB0+&
                            dot_product(B0_DIIS(:,i2,j2,iSide),&
                            FaceArea_DIIS(:,i2,j2,iSide))    
                       do iDirB0=x_,z_
                          do iDirFA=x_,z_
                             B0SourceMatrix_DDCB(iDirFA,iDirB0,i,j,k,iBlock)= & 
                                  B0SourceMatrix_DDCB(iDirFA,iDirB0,i,j,k,iBlock)&
                                  +FaceArea_DIIS(iDirB0,i2,j2,iSide)*&
                                  B0_DIIS(iDirFA,i2,j2,iSide)&  
                                  -FaceArea_DIIS(iDirFA,i2,j2,iSide)*&
                                  B0_DIIS(iDirB0,i2,j2,iSide)
                          end do
                       end do
                    end do
                 end do
              end do

              B0SourceMatrix_DDCB(1,1,i,j,k,iBlock) = -DivB0           
              B0SourceMatrix_DDCB(2,2,i,j,k,iBlock) = -DivB0
              B0SourceMatrix_DDCB(3,3,i,j,k,iBlock) = -DivB0    

              B0SourceMatrix_DDCB(:,:,i,j,k,iBlock) = &
                   B0SourceMatrix_DDCB(:,:,i,j,k,iBlock)*&
                   vInv_CB(i,j,k,iBlock)
           end if
        end do
     end do
  end do
contains
  subroutine correct_b0_face(iSide)
    implicit none
    integer,intent(in)::iSide
    select case(iSide)
    case(East_,West_)
       iFace=1+nI*(iSide-East_)
       !Arguments of the get_nodes_and_ref_b0 are:
       !the face center generalized coordinates 
       !minus
       !the generalized coordinates of the (1,1,1) cell center,
       !divided by dGenRef_D
       call get_nodes_and_ref_b0(2*iFace-3,2*j-2,2*k-2)
       call get_face_area_i(RefXyzNodes_DIII(:,1:1,:,:),&
                            1,1,0,1,0,1,&
                            RefFaceArea_DIII(:,1:1,:,:))
       do j2=0,1
          do i2=0,1
             B0_DIIS(:,i2,j2,iSide)=(RefB0_DIII(:,0,i2,j2)+&
                  RefB0_DIII(:,1,i2,j2))*cHalf
             FaceArea_DIIS(:,i2,j2,iSide)=RefFaceArea_DIII(:,1,i2,j2)
          end do
       end do
    case(South_,North_)
       jFace=1+nJ*(iSide-South_)
       !Arguments of the get_nodes_and_ref_b0 are:
       !the face center generalized coordinates 
       !minus
       !the generalized coordinates of the (1,1,1) cell center,
       !divided by dXyzRef_D
       call get_nodes_and_ref_b0(2*i-2,2*jFace-3,2*k-2)
       call get_face_area_j(RefXyzNodes_DIII(:,:,1:1,:),&
                            0,1,1,1,0,1,&
                            RefFaceArea_DIII(:,:,1:1,:))
       do j2=0,1
          do i2=0,1
             B0_DIIS(:,i2,j2,iSide)=(RefB0_DIII(:,i2,0,j2)+&
                  RefB0_DIII(:,i2,1,j2))*cHalf
             FaceArea_DIIS(:,i2,j2,iSide)=RefFaceArea_DIII(:,i2,1,j2)
          end do
       end do
    case(Bot_,Top_)
       kFace=1+nK*(iSide-Bot_)
       !Arguments of the get_nodes_and_ref_b0 are:
       !the face center generalized coordinates 
       !minus
       !the generalized coordinates of the (1,1,1) cell center,
       !divided by dXyzRef_D
       call get_nodes_and_ref_b0(2*i-2,2*j-2,2*kFace-3)
       call get_face_area_k(RefXyzNodes_DIII(:,:,:,1:1),&
                            0,1,0,1,1,1,&
                            RefFaceArea_DIII(:,:,:,1:1))
       do j2=0,1
          do i2=0,1
             B0_DIIS(:,i2,j2,iSide)=(RefB0_DIII(:,i2,j2,0)+&
                  RefB0_DIII(:,i2,j2,1))*cHalf
             FaceArea_DIIS(:,i2,j2,iSide)=RefFaceArea_DIII(:,i2,j2,1)
          end do
       end do
    end select
  end subroutine correct_b0_face

  subroutine get_nodes_and_ref_b0(iRef,jRef,kRef)
    implicit none
    integer,intent(in)::iRef,jRef,kRef
    integer::k2         
    !Get the face center generalized coordinates, which is also
    !the generalized coordinates of the refined grid node (111),
    !if UseVertexBasedGrid=.true.
 
    GenCoord111_D=XyzStart_BLK(:,iBlock)+(/iRef,jRef,kRef/)*dGenFine_D
    
    !Get cell center coordinates of the refined grid
    !and nodes coordinates
    if(UseVertexBasedGrid)then
       call gen_to_xyz_arr(GenCoord111_D+cHalf*dGenFine_D,&
            DGenFine_D(1),DGenFine_D(2),dGenFine_D(3),&
            0,1,0,1,0,1,&
            RefXyz_DIII(x_,0:1,0:1,0:1),&
            RefXyz_DIII(y_,0:1,0:1,0:1),&
            RefXyz_DIII(z_,0:1,0:1,0:1))       
       call gen_to_xyz_arr(GenCoord111_D,&
            DGenFine_D(1),DGenFine_D(2),dGenFine_D(3),&
            0,2,0,2,0,2,&
            RefXyzNodes_DIII(x_,0:2,0:2,0:2),&
            RefXyzNodes_DIII(y_,0:2,0:2,0:2),&
            RefXyzNodes_DIII(z_,0:2,0:2,0:2))
    else
       !We need a wider stencil to construct the nodes
       call gen_to_xyz_arr(GenCoord111_D+cHalf*dGenFine_D,&
            DGenFine_D(1),DGenFine_D(2),dGenFine_D(3),&
            -1,2,-1,2,-1,2,&
            RefXyz_DIII(x_,-1:2,-1:2,-1:2),&
            RefXyz_DIII(y_,-1:2,-1:2,-1:2),&
            RefXyz_DIII(z_,-1:2,-1:2,-1:2))
!REDO!!!!!!!!!
       call cell_centers_to_nodes
    !  call gen_to_xyz_arr(GenCoord111_D,&
    !        DGenFine_D(1),DGenFine_D(2),dGenFine_D(3),&
    !        0,2,0,2,0,2,&
    !        RefXyzNodes_DIII(x_,0:2,0:2,0:2),&
    !        RefXyzNodes_DIII(y_,0:2,0:2,0:2),&
    !        RefXyzNodes_DIII(z_,0:2,0:2,0:2))
   end if
     
    do k2=0,1; do j2=0,1; do i2=0,1
       call  get_B0(&
            RefXyz_DIII(x_,i2,j2,k2),&
            RefXyz_DIII(y_,i2,j2,k2),&
            RefXyz_DIII(z_,i2,j2,k2),&
            RefB0_DIII(:,i2,j2,k2))
    end do; end do; end do	
  end subroutine get_nodes_and_ref_b0
  subroutine cell_centers_to_nodes
    real,dimension(3,3):: A_DD, A1_DD
    real,dimension(3)  :: B_D
    real               :: DetInv
    integer::k2
    real,dimension(-1:2,-1:2,-1:2)::R2_III
    !------------------------------------------------------------------------
    do k2=-1,2; do j2=-1,2; do i2=-1,2
       R2_III(i2,j2,k2)=sum(RefXyz_DIII(:,i2,j2,k2)**2)
    end do; end do; end do

    do k2=0,2; do j2=0,2; do i2=0,2
       
       A_DD(1,:)=RefXyz_DIII(:,i2,j2,k2)-RefXyz_DIII(:,i2-1,j2,k2)
       A_DD(2,:)=RefXyz_DIII(:,i2,j2,k2)-RefXyz_DIII(:,i2,j2-1,k2)
       A_DD(3,:)=RefXyz_DIII(:,i2,j2,k2)-RefXyz_DIII(:,i2,j2,k2-1)

       DetInv=cOne/det(A_DD)

       B_D(1)=cHalf*(R2_III(i2,j2,k2)-R2_III(i2-1,j2,k2))
       B_D(2)=cHalf*(R2_III(i2,j2,k2)-R2_III(i2,j2-1,k2))
       B_D(3)=cHalf*(R2_III(i2,j2,k2)-R2_III(i2,j2,k2-1))

       A1_DD(:,2:3)=A_DD(:,2:3)
       A1_DD(:,1)=B_D

       RefXyzNodes_DIII(x_,i2,j2,k2) = det(A1_DD)*DetInv

       A1_DD(:,1)=A_DD(:,1)
       A1_DD(:,3)=A_DD(:,3)
       A1_DD(:,2)=B_D

       RefXyzNodes_DIII(y_,i2,j2,k2) = det(A1_DD)*DetInv
       
       A1_DD(:,1:2)=A_DD(:,1:2)
       A1_DD(:,3)=B_D

       RefXyzNodes_DIII(z_,i2,j2,k2) = det(A1_DD)*DetInv
    end do; end do; end do
    
  end subroutine cell_centers_to_nodes
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
end subroutine calc_b0source_covar

subroutine covariant_force_integral(i,j,k,iBLK,Fai_S)
  use ModSize
  use ModCovariant
  use ModAdvance,ONLY: fbody_x_BLK,fbody_y_BLK,fbody_z_BLK
  use ModGeometry,ONLY: vInv_CB
  implicit none

  integer,intent(in)::i,j,k,iBLK
  real,dimension(East_:Top_),intent(in)::Fai_S
  real:: FaceArea_DS(3,east_:top_),VInv


  VInv=vInv_CB(i,j,k,iBLK)

  FaceArea_DS(:,East_ :West_ )=FaceAreaI_DFB(:,i:i+1,j,k,iBLK)
  FaceArea_DS(:,South_:North_)=FaceAreaJ_DFB(:,i,j:j+1,k,iBLK)
  FaceArea_DS(:,Bot_  :Top_  )=FaceAreaK_DFB(:,i,j,k:k+1,iBLK)
 
  fbody_x_BLK(i,j,k,iBLK) = VInv*&              
       dot_product(FaceArea_DS(1,:),Fai_S(:))                     
  fbody_y_BLK(i,j,k,iBLK) = VInv*&              
       dot_product(FaceArea_DS(2,:),Fai_S(:))                     
  fbody_z_BLK(i,j,k,iBLK) = VInv*&              
       dot_product(FaceArea_DS(3,:),Fai_S(:))        

end subroutine covariant_force_integral
!-----------------------------------------------------------------------------
subroutine covariant_gradient(iBlock, Var_G,&     
     GradientX_G, GradientY_G, GradientZ_G)
  use ModSize
  use ModMain, ONLY: x_, y_, z_
  use ModCovariant
  use ModGeometry,ONLY:body_blk, true_cell, &
       vInv_CB
  use ModNumConst
  implicit none

  integer,intent(in) :: iBlock

  real, dimension(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn),&
       intent(in) :: Var_G
  real, dimension(0:nI+1, 0:nJ+1, 0:nK+1),&
       intent(out) ::  GradientX_G, GradientY_G, GradientZ_G

  real, dimension(0:nI+1, 0:nJ+1, 0:nK+1) :: OneTrue_G

  integer :: i, j, k

  real,dimension(3,east_:top_) :: FaceArea_DS
  real,dimension(east_:top_) :: Difference_S

  real::VInvHalf

  !To fill in the ghostcells
  GradientX_G = cZero
  GradientY_G = cZero
  GradientZ_G = cZero

  if(.not.body_BLK(iBlock)) then
     do k=1,nK; do j=1,nJ; do i=1,nI
        VInvHalf=chalf*vInv_CB(i,j,k,iBlock)

        FaceArea_DS(:,East_ :West_ )=FaceAreaI_DFB(:,i:i+1,j,k,iBlock)
        FaceArea_DS(:,South_:North_)=FaceAreaJ_DFB(:,i,j:j+1,k,iBlock)
        FaceArea_DS(:,Bot_  :Top_  )=FaceAreaK_DFB(:,i,j,k:k+1,iBlock)

        Difference_S(East_) = -Var_G(i-1,j,k)
        Difference_S(West_) = +Var_G(i+1,j,k)
        Difference_S(South_)= -Var_G(i,j-1,k)
        Difference_S(North_)= +Var_G(i,j+1,k)
        Difference_S(Bot_)  = -Var_G(i,j,k-1)
        Difference_S(Top_)  = +Var_G(i,j,k+1)
  
        GradientX_G(i,j,k) = &
             dot_product(FaceArea_DS(x_,:),Difference_S)*VInvHalf
        GradientY_G(i,j,k) = &
             dot_product(FaceArea_DS(y_,:),Difference_S)*VInvHalf
        GradientZ_G(i,j,k) = &
             dot_product(FaceArea_DS(z_,:),Difference_S)*VInvHalf
     end do; end do; end do
  else
     where(true_cell(0:nI+1, 0:nJ+1, 0:nK+1,iBlock)) 
        OneTrue_G=cOne
     elsewhere
        OneTrue_G=cZero
     end where
     do k=1,nK;  do j=1,nJ;  do i=1,nI
        VInvHalf=&
             chalf*vInv_CB(i,j,k,iBlock)*OneTrue_G(i,j,k)
        
        FaceArea_DS(:,East_ :West_ )=FaceAreaI_DFB(:,i:i+1,j,k,iBlock)
        FaceArea_DS(:,South_:North_)=FaceAreaJ_DFB(:,i,j:j+1,k,iBlock)
        FaceArea_DS(:,Bot_  :Top_  )=FaceAreaK_DFB(:,i,j,k:k+1,iBlock)        

        Difference_S(East_) =  OneTrue_G(i-1,j,k)*&
             (Var_G(i,j,k)-Var_G(i-1,j,k))+&
             (cOne- OneTrue_G(i-1,j,k))*&
              OneTrue_G(i+1,j,k)*&
              (Var_G(i+1,j,k)-Var_G(i,j,k))

        Difference_S(West_) =  OneTrue_G(i+1,j,k)*&
             (Var_G(i+1,j,k)-Var_G(i,j,k))+&
             (cOne- OneTrue_G(i+1,j,k))*&
             OneTrue_G(i-1,j,k)*&
             (Var_G(i,j,k)-Var_G(i-1,j,k))

        Difference_S(South_)=  OneTrue_G(i,j-1,k)*&
             (Var_G(i,j,k)-Var_G(i,j-1,k))+&
             (cOne-OneTrue_G(i,j-1,k))*&
             OneTrue_G(i,j+1,k)*&
             (Var_G(i,j+1,k)-Var_G(i,j,k))

        Difference_S(North_)=  OneTrue_G(i,j+1,k)*&
             (Var_G(i,j+1,k)-Var_G(i,j,k))+&
             (cOne-OneTrue_G(i,j+1,k))*&
             OneTrue_G(i,j-1,k)*&
             (Var_G(i,j,k)-Var_G(i,j-1,k))

        Difference_S(Bot_)  =  OneTrue_G(i,j,k-1)*&
             (Var_G(i,j,k)-Var_G(i,j,k-1))+&
             (cOne-OneTrue_G(i,j,k-1))*&
             OneTrue_G(i,j,k+1)*&
             (Var_G(i,j,k+1)-Var_G(i,j,k))

        Difference_S(Top_)  =  OneTrue_G(i,j,k+1)*&
             (Var_G(i,j,k+1)-Var_G(i,j,k))+&
             (cOne-OneTrue_G(i,j,k+1))*&
             OneTrue_G(i,j,k-1)*&
             (Var_G(i,j,k)-Var_G(i,j,k-1))

        GradientX_G(i,j,k) = &
             dot_product(FaceArea_DS(x_,:),Difference_S)*VInvHalf
        GradientY_G(i,j,k) = &
             dot_product(FaceArea_DS(y_,:),Difference_S)*VInvHalf
        GradientZ_G(i,j,k) = &
             dot_product(FaceArea_DS(z_,:),Difference_S)*VInvHalf
     end do; end do; end do
  end if

end subroutine covariant_gradient
!-----------------------------------------------------------------------------
subroutine covariant_curlb(i,j,k,iBLK,CurlB_D)
  use ModSize
  use ModVarIndexes,ONLY: Bx_,By_,Bz_
  use ModCovariant
  use ModNumConst
  use Modgeometry,ONLY:vInv_CB
  use ModAdvance,ONLY:State_VGB
  implicit none
  integer,intent(in)::i,j,k,iBLK
  real,dimension(3),intent(out)::CurlB_D

  real,dimension(3,east_:top_)::MagneticField_DS, FaceArea_DS
  real::VInvHalf

  VInvHalf=chalf*vInv_CB(i,j,k,iBLK)

  FaceArea_DS(:,East_ :West_ )=FaceAreaI_DFB(:,i:i+1,j,k,iBLK)
  FaceArea_DS(:,South_:North_)=FaceAreaJ_DFB(:,i,j:j+1,k,iBLK)
  FaceArea_DS(:,Bot_  :Top_  )=FaceAreaK_DFB(:,i,j,k:k+1,iBLK)   

  MagneticField_DS(:,East_ )=-State_VGB(Bx_:Bz_,i-1,j,k,iBLK)
  MagneticField_DS(:,West_ )= State_VGB(Bx_:Bz_,i+1,j,k,iBLK)
  MagneticField_DS(:,South_)=-State_VGB(Bx_:Bz_,i,j-1,k,iBLK)
  MagneticField_DS(:,North_)= State_VGB(Bx_:Bz_,i,j+1,k,iBLK)
  MagneticField_DS(:,Bot_)  =-State_VGB(Bx_:Bz_,i,j,k-1,iBLK)
  MagneticField_DS(:,Top_)  = State_VGB(Bx_:Bz_,i,j,k+1,iBLK)
 
  CurlB_D(1)=dot_product(FaceArea_DS(2,:),MagneticField_DS(3,:))-&
       dot_product(FaceArea_DS(3,:),MagneticField_DS(2,:))
  CurlB_D(2)=dot_product(FaceArea_DS(3,:),MagneticField_DS(1,:))-&
       dot_product(FaceArea_DS(1,:),MagneticField_DS(3,:))
  CurlB_D(3)=dot_product(FaceArea_DS(1,:),MagneticField_DS(2,:))-&
       dot_product(FaceArea_DS(2,:),MagneticField_DS(1,:))
  CurlB_D=VInvHalf*CurlB_D
end subroutine covariant_curlb
!===========================================================
!===========================================================
subroutine covar_curlb_plotvar(iDir,iBLK,PlotVar_G)
  use ModSize
  use ModNumConst
  implicit none

  integer,intent(in):: iDir
  integer,intent(in):: iBLK
  real,dimension(-1:nI+2,-1:nJ+2,-1:nK+2),intent(out):: PlotVar_G
  real,dimension(1:nDim):: CurlB_D
  integer::i,j,k
  PlotVar_G=cZero
  do k=1,nK
     do j=1,nJ
        do i=1,nI
           call covariant_curlb(i,j,k,iBLK,CurlB_D)
           PlotVar_G(i,j,k)=CurlB_D(iDir)
        end do
     end do
  end do

end subroutine covar_curlb_plotvar
!===========================================================
subroutine covar_curlbr_plotvar(iBLK,PlotVar_G)
  use ModSize
  use ModMain,ONLY:x_,y_,z_
  use ModGeometry,ONLY:x_BLK,y_BLK,z_BLK,R_BLK
  use ModNumConst
  implicit none

  integer,intent(in):: iBLK
  real,dimension(-1:nI+2,-1:nJ+2,-1:nK+2),intent(out):: PlotVar_G

  real,dimension(1:nDim):: CurlB_D
  integer::i,j,k

  PlotVar_G=cZero
  do k=1,nK
     do j=1,nJ
        do i=1,nI
           call covariant_curlb(i,j,k,iBLK,CurlB_D)
           PlotVar_G(i,j,k)=(CurlB_D(x_)*x_BLK(i,j,k,iBLK)+&
                CurlB_D(y_)*y_BLK(i,j,k,iBLK)+&
                CurlB_D(z_)*z_BLK(i,j,k,iBLK))/&
                R_BLK(i,j,k,iBLK)
        end do
     end do
  end do

end subroutine covar_curlbr_plotvar


!=======================================================================
subroutine save_bn_faceI(iFaceOut,iFaceIn,iBlock)
  use ModMain,ONLY: nJ,nK, nDim,x_,y_,z_
  use ModVarIndexes,ONLY:Bx_, By_, Bz_
  use ModAdvance, ONLY:BnL_,BnR_,CorrectedFlux_VXB,&
       LeftState_VX,RightState_VX
  use ModCovariant
  implicit none

  integer,intent(in) :: iFaceOut,iFaceIn,iBlock

  integer :: j,k



  do k=1,nK; do j=1,nJ
     CorrectedFlux_VXB(BnL_,j,k,iFaceOut,iBlock) = &
          dot_product(LeftState_VX(Bx_:Bz_,iFaceIn,j,k),&
          FaceAreaI_DFB(:,iFaceIn,j,k,iBlock))
     CorrectedFlux_VXB(BnR_,j,k,iFaceOut,iBlock) = &
          dot_product(RightState_VX(Bx_:Bz_,iFaceIn,j,k),&
          FaceAreaI_DFB(:,iFaceIn,j,k,iBlock))
  end do; end do

end subroutine save_bn_faceI

!-------------------------------------------------------------------------
subroutine save_bn_faceJ(jFaceOut,jFaceIn,iBlock)
  use ModMain,ONLY: nI,nK, nDim,x_,y_,z_
  use ModVarIndexes,ONLY: Bx_, By_, Bz_
  use ModAdvance, ONLY: BnL_,BnR_,CorrectedFlux_VYB,&
       LeftState_VY,RightState_VY
  use ModCovariant
  implicit none

  integer,intent(in) :: jFaceOut,jFaceIn,iBlock

  integer :: i,k

  do k=1,nK; do i=1,nI
     CorrectedFlux_VYB(BnL_,i,k,jFaceOut,iBlock) = &
          dot_product(LeftState_VY(Bx_:Bz_,i,jFaceIn,k),&
          FaceAreaJ_DFB(:,i,jFaceIn,k,iBlock))
     CorrectedFlux_VYB(BnR_,i,k,jFaceOut,iBlock) = &
          dot_product(RightState_VY(Bx_:Bz_,i,jFaceIn,k),&
          FaceAreaJ_DFB(:,i,jFaceIn,k,iBlock))
  end do; end do

end subroutine save_bn_faceJ

!---------------------------------------------------------------------------
subroutine save_bn_faceK(kFaceOut,kFaceIn,iBlock)
  use ModMain,ONLY: nI, nJ, nDim, x_, y_, z_ 
  use ModVarIndexes,ONLY:Bx_, By_, Bz_
  use ModAdvance, ONLY:BnL_,BnR_,CorrectedFlux_VZB,&
       LeftState_VZ,RightState_VZ
  use ModCovariant
  implicit none

  integer,intent(in) :: kFaceOut,kFaceIn,iBlock

  integer :: i,j
 


  do j=1,nJ; do i=1,nI
     CorrectedFlux_VZB(BnL_,i,j,kFaceOut,iBlock) = &
          dot_product(LeftState_VZ(Bx_:Bz_,i,j,kFaceIn),&
          FaceAreaK_DFB(:,i,j,kFaceIn,iBlock))
     CorrectedFlux_VZB(BnR_,i,j,kFaceOut,iBlock) =&
          dot_product(RightState_VZ(Bx_:Bz_,i,j,kFaceIn),&
          FaceAreaK_DFB(:,i,j,kFaceIn,iBlock))
  end do; end do

end subroutine save_bn_faceK
!---------------------------------------------------------------------------
subroutine apply_bn_faceI(iFaceIn,iFaceOut,iBlock)
  use ModMain,ONLY: nJ,nK, nDim,x_,y_,z_
  use ModVarIndexes,ONLY:Bx_, By_, Bz_
  use ModAdvance, ONLY: BnL_,BnR_,CorrectedFlux_VXB,&
       LeftState_VX,RightState_VX
  use ModCovariant
  use ModGeometry,ONLY:true_cell
  implicit none

  integer,intent(in) :: iFaceOut,iFaceIn,iBlock

  integer :: j,k
  real,dimension(nDim) :: B_D,FaceArea_D
  
  real:: FaceArea2,DeltaBDotFA

  do k=1,nK; do j=1,nJ
     if(.not.all(true_cell(iFaceOut-1:iFaceOut,j,k,iBlock)))CYCLE
     FaceArea_D=FaceAreaI_DFB(:,iFaceOut,j,k,iBlock)
     FaceArea2=dot_product(FaceArea_D,FaceArea_D)

     B_D=LeftState_VX(Bx_:Bz_,iFaceOut,j,k)

     DeltaBDotFA = (CorrectedFlux_VXB(BnL_,j,k,iFaceIn,iBlock) -&
          dot_product(B_D,FaceArea_D))/FaceArea2

     LeftState_VX(Bx_:Bz_,iFaceOut,j,k)=B_D(x_)+DeltaBDotFA*FaceArea_D
  
     B_D=RightState_VX(Bx_:Bz_,iFaceOut,j,k)

     DeltaBDotFA = (CorrectedFlux_VXB(BnR_,j,k,iFaceIn,iBlock) -&
          dot_product(B_D,FaceArea_D))/FaceArea2

     RightState_VX(Bx_:Bz_,iFaceOut,j,k)=B_D+DeltaBDotFA*FaceArea_D
  end do; end do

end subroutine apply_bn_faceI

!-------------------------------------------------------------------------
subroutine apply_bn_faceJ(jFaceIn,jFaceOut,iBlock)
  use ModMain,ONLY: nI,nK, nDim,x_,y_,z_
  use ModVarIndexes,ONLY:Bx_, By_, Bz_
  use ModAdvance, ONLY: BnL_,BnR_,CorrectedFlux_VYB,&
       LeftState_VY,RightState_VY
  use ModCovariant
  use ModGeometry,ONLY:true_cell
  implicit none

  integer,intent(in) :: jFaceOut,jFaceIn,iBlock

  integer :: i,k
  real,dimension(nDim) :: B_D,FaceArea_D

  real:: FaceArea2,DeltaBDotFA

  do k=1,nK; do i=1,nI
     if(.not.all(true_cell(i,jFaceOut-1:jFaceOut,k,iBlock)))CYCLE
     FaceArea_D=FaceAreaJ_DFB(:,i,jFaceOut,k,iBlock)
     FaceArea2=dot_product(FaceArea_D,FaceArea_D)

     B_D=LeftState_VY(Bx_:Bz_,i,jFaceOut,k)

     DeltaBDotFA = (CorrectedFlux_VYB(BnL_,i,k,jFaceIn,iBlock)-&
          dot_product(B_D,FaceArea_D))/FaceArea2

     LeftState_VY(Bx_:Bz_,i,jFaceOut,k)=B_D+DeltaBDotFA*FaceArea_D

     B_D=RightState_VY(Bx_:Bz_,i,jFaceOut,k)

     DeltaBDotFA = (CorrectedFlux_VYB(BnR_,i,k,jFaceIn,iBlock)-&
          dot_product(B_D,FaceArea_D))/FaceArea2

     RightState_VY(Bx_:Bz_,i,jFaceOut,k)=B_D+DeltaBDotFA*FaceArea_D

  end do; end do

end subroutine apply_bn_faceJ

!---------------------------------------------------------------------------
subroutine apply_bn_faceK(kFaceIn,kFaceOut,iBlock)
  use ModMain,ONLY: nI, nJ, nDim, x_, y_, z_
  use ModVarIndexes,ONLY:Bx_, By_, Bz_
  use ModAdvance, ONLY: BnL_,BnR_,CorrectedFlux_VZB,&
       LeftState_VZ,RightState_VZ
  use ModCovariant
  use ModGeometry,ONLY:true_cell
  implicit none

  integer,intent(in) :: kFaceOut,kFaceIn,iBlock

  integer :: i,j
  real,dimension(nDim) :: B_D, FaceArea_D

  real:: FaceArea2,DeltaBDotFA

  do j=1,nJ; do i=1,nI
     if(.not.all(true_cell(i,j,kFaceOut-1:kFaceOut,iBlock)))CYCLE
     FaceArea_D=FaceAreaK_DFB(:,i,j,kFaceOut,iBlock)
     FaceArea2=dot_product(FaceArea_D,FaceArea_D)

     B_D=LeftState_VZ(Bx_:Bz_,i,j,kFaceOut)

     DeltaBDotFA = ( CorrectedFlux_VZB(BnL_,i,j,kFaceIn,iBlock) -&
          dot_product(B_D,FaceArea_D))/FaceArea2

     LeftState_VZ(Bx_:Bz_,i,j,kFaceOut) = B_D+DeltaBDotFA*FaceArea_D  

     B_D=RightState_VZ(Bx_:Bz_,i,j,kFaceOut)
     
     DeltaBDotFA = (CorrectedFlux_VZB(BnR_,i,j,kFaceIn,iBlock) -&
          dot_product(B_D,FaceArea_D))/FaceArea2

     RightState_VZ(Bx_:Bz_,i,j,kFaceOut) = B_D+DeltaBDotFA*FaceArea_D  
  end do; end do

end subroutine apply_bn_faceK
!---------------------------------------------------------------------------


!=========================================End covariant.f90=================



!============To be moved to src/covariant_facefluxes.f90=====================



!=========End src/covariant_facefluxes.f90========================

real function integrate_BLK_covar(qnum,qa)             

  ! Return the volume integral of qa, ie. sum(qa*cV_BLK) 
  ! for all used blocks and true cells
  ! Do for each processor separately if qnum=1, otherwise add them all

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,nBLK,nIJK,nBlockMax,unusedBLK
  use ModGeometry, ONLY :&
                          vInv_CB,&                   
                          true_BLK,true_cell
  use ModMpi
  implicit none 

  ! Arguments

  integer, intent(in) :: qnum
  real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(in) :: qa

  ! Local variables:
  real    :: qsum, qsum_all
  integer :: iBLK, iError,i,j,k

  logical :: oktest, oktest_me

  !---------------------------------------------------------------------------

  call set_oktest('integrate_BLK',oktest, oktest_me)

  qsum=0.0
                                                     
  do iBLK=1,nBlockMax
     if(.not.unusedBLK(iBLK)) then
        if(true_BLK(iBLK)) then
           do k=1,nK
              do j=1,nJ
                 do i=1,nI
                    qsum=qsum + qa(i,j,k,iBLK)/&
                         vInv_CB(i,j,k,iBLK)
                 end do
              end do
           end do
        else
           do k=1,nK
              do j=1,nJ
                 do i=1,nI
                    if(true_cell(i,j,k,iBLK))&
                    qsum=qsum + qa(i,j,k,iBLK)/&
                         vInv_CB(i,j,k,iBLK)
                 end do
              end do
           end do
        end if
     end if
  end do
                                                    
  if(qnum>1)then
     call MPI_allreduce(qsum, qsum_all, 1,  MPI_REAL, MPI_SUM, &
          iComm, iError)
     integrate_BLK_covar=qsum_all
     if(oktest)write(*,*)'me,sum,sum_all:',iProc,qsum,qsum_all
  else
     integrate_BLK_covar=qsum
     if(oktest)write(*,*)'me,qsum:',iProc,qsum
  end if
  
end function integrate_BLK_covar

subroutine integrate_domain_covar(Sum_V, Pressure_GB)
  use ModAdvance,   ONLY: State_VGB, tmp2_BLK, nVar
  use ModMain,      ONLY: nI, nJ, nK, nIJK, nBlock, MaxBlock, UnusedBLK
  use ModGeometry,  ONLY: vInv_CB, true_BLK, true_cell
  use ModVarIndexes,ONLY: P_
  use ModNumConst,  ONLY: cZero, cOne
  implicit none 

  ! Arguments
  real, intent(out) :: Sum_V(nVar)
  real, intent(out) :: Pressure_GB(-1:nI+2,-1:nJ+2,-1:nK+2,MaxBlock)

  ! Local variables:
  real    :: CellVolume
  integer :: iBlock, iVar, i, j, k
  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------

  call set_oktest('integrate_domain',DoTest, DoTestMe)

  Sum_V=cZero
                                                     
  do iBlock = 1, nBlock
     if(unusedBLK(iBlock)) CYCLE
     if(true_BLK(iBlock)) then
        do k=1,nK; do j=1,nJ; do i=1,nI
           CellVolume=cOne/vInv_CB(i,j,k,iBlock)
           do iVar=1,nVar
              Sum_V(iVar)=Sum_V(iVar) + State_VGB(iVar,i,j,k,iBlock)*CellVolume
           end do
        end do; end do; end do
     else
        do k=1,nK; do j=1,nJ; do i=1,nI
           if(.not.true_cell(i,j,k,iBlock))CYCLE
           CellVolume=cOne/vInv_CB(i,j,k,iBlock)
           do iVar = 1, nVar
              Sum_V(iVar)=Sum_V(iVar) + State_VGB(iVar,i,j,k,iBlock)*CellVolume
           end do
        end do; end do; end do
     end if
     Pressure_GB(1:nI,1:nJ,1:nK,iBlock) = State_VGB(P_,1:nI,1:nJ,1:nK,iBlock)
  end do

end subroutine integrate_domain_covar
!==============End of library programms=================================
