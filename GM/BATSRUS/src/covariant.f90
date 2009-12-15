!^CFG COPYRIGHT UM
subroutine gen_to_xyz_arr(&
     GenCoord111_D,    &!(in)Gen.coords.of the point (1,1,1)
     dGen1,dGen2,dGen3,&!(in)Mesh sizes in Gen coords.
     iStart,iMax,      &!(in)The first and the last value of i index
     jStart,jMax,      &!(in)The first and the last value of j index
     kStart,kMax,      &!(in)The first and the last value of k index
     X_C,              &!(out)Cartesian x-coord. of the mapped points
     Y_C,              &!(out)Cartesian y-coord. of the mapped points    
     Z_C)               !(out)Cartesian z-coord. of the mapped points
 
  !subroutine gen_to_xyz_arr maps the piece of an equally spaced grid
  !in the space of GENERALIZED COORDINATES to the cartesian space xyz

  use ModNumConst
  use ModCovariant,ONLY:TypeGeometry,rTorusLarge,rTorusSmall,&
                        gen_to_r, r_to_gen
  use ModMain,     ONLY:nDim,R_,Phi_,Theta_,x_,y_,z_
  implicit none
  real,intent(in):: GenCoord111_D(nDim),dGen1,dGen2,dGen3
  integer,intent(in)::  iStart,iMax,jStart,jMax,kStart,kMax
  real,dimension(iStart:iMax,jStart:jMax,kStart:kMax),&
       intent(out)::X_C,Y_C,Z_C

  integer::i,j,k
  real::R,Theta,Phi,sinTheta,cosTheta,sinPhi,cosPhi,Gen1
  real::PoloidalAngle,Z,StretchCoef
  real,external::wall_radius
  !----------------------------------------------------------

  select case(TypeGeometry)
  case('cartesian','rz')
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
     ! Gen1=R, Gen2=Phi, Gen3=Latitude
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

              X_C(i,j,k) = R*cosTheta*cosPhi                      
              Y_C(i,j,k) = R*cosTheta*sinPhi                      
              Z_C(i,j,k) = R*sinTheta
           end do
        end do
     end do
  case('spherical_lnr')
     ! Gen1=log(R), Gen2=Phi, Gen3=Latitude
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

              X_C(i,j,k) = R*cosTheta*cosPhi                      
              Y_C(i,j,k) = R*cosTheta*sinPhi                      
              Z_C(i,j,k) = R*sinTheta
           end do
        end do
     end do
  case('spherical_genr')
     ! Gen1=general variable from 0 to 1, Gen2=Phi, Gen3=Latitude
     do k = kStart, kMax
        Theta      =      (k-1)*dGen3 + GenCoord111_D(Theta_)
        sinTheta   =  sin(Theta)
        cosTheta   =  cos(Theta)
        do j = jStart, jMax
           Phi     =      (j-1)*dGen2 + GenCoord111_D(Phi_)
           sinPhi  =  sin(Phi)
           cosPhi  =  cos(Phi)
           do i = iStart, iMax
              Gen1 = (i-1)*dGen1 + GenCoord111_D(R_)
              R    =   gen_to_r(Gen1)

              X_C(i,j,k) = R*cosTheta*cosPhi                      
              Y_C(i,j,k) = R*cosTheta*sinPhi                      
              Z_C(i,j,k) = R*sinTheta
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
  case('axial_torus')
     ! Gen1=r, Gen2=Phi, Gen3=z
     do k = kStart, kMax
        do j = jStart, jMax
           Phi     =  (j-1)*dGen2 + GenCoord111_D(Phi_)
           sinPhi  =  sin(Phi)
           cosPhi  =  cos(Phi)
           do i = iStart, iMax
              Z = (k-1)*dGen3 + GenCoord111_D(z_)
              R    =  (i-1)*dGen1 + GenCoord111_D(R_)-rTorusLarge
              if(.not.(R==cZero.and.Z==cZero))then
                 PoloidalAngle=atan2(Z,R)
                 if(PoloidalAngle<cZero)&
                      PoloidalAngle=PoloidalAngle+cTwoPi
                 StretchCoef = wall_radius(PoloidalAngle)/rTorusSmall*&
                      max(abs(cos(PoloidalAngle)),&
                      abs(sin(PoloidalAngle)))
                 Z = Z*StretchCoef
                 R=R*StretchCoef
              end if
              Z_C(i,j,k)=Z
              R=R+rTorusLarge
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
subroutine xyz_to_gen(XyzIn_D,GenOut_D)
  use ModNumConst
  use ModCovariant
  use ModMain,     ONLY:nDim,R_,Phi_,Theta_,x_,y_,z_
  implicit none
  real,dimension(nDim),intent(in) ::XyzIn_D
  real,dimension(nDim),intent(out)::GenOut_D
  real::PoloidalAngle,R,Z,StretchCoef
  real,external::wall_radius

  select case(TypeGeometry)           
  case('cartesian','rz')                   
     GenOut_D=XyzIn_D
  case('spherical','spherical_lnr','spherical_genr')   
     call xyz_to_spherical(XyzIn_D(x_),XyzIn_D(y_),XyzIn_D(z_),&
          GenOut_D(R_),GenOut_D(Phi_),GenOut_D(Theta_))
     !From colatitude to latitude:
     GenOut_D(Theta_)=cHalfPi-GenOut_D(Theta_)
     if(TypeGeometry=='spherical_lnr')&
          GenOut_D(R_)=log(max(GenOut_D(R_),cTiny))
     if(TypeGeometry=='spherical_genr')&
          GenOut_D(R_)=r_to_gen(GenOut_D(R_))

  case('cylindrical')
     GenOut_D(R_)   = sqrt(sum(XyzIn_D(x_:y_)**2))
     GenOut_D(Phi_) = modulo(atan2(XyzIn_D(y_), XyzIn_D(x_)), cTwoPi)
     GenOut_D(z_)   = XyzIn_D(z_)

  case('axial_torus')
     if(all(XyzIn_D(x_:y_)==cZero))&
          call stop_mpi(&
          'axil_torus geometry does not work for points at the pole')
     GenOut_D(Phi_)=atan2(XyzIn_D(y_),XyzIn_D(x_))
     if(GenOut_D(Phi_)<cZero)GenOut_D(Phi_)=GenOut_D(Phi_)+cTwoPi
     R=sqrt(sum(XyzIn_D(x_:y_)**2))-rTorusLarge
     Z=XyzIn_D(z_)
     if(.not.(Z==cZero.and.R==cZero))then
        PoloidalAngle=atan2(Z,R)
        if(PoloidalAngle<cZero)&
             PoloidalAngle= PoloidalAngle+cTwoPi
        StretchCoef=rTorusSmall/&
             (max(abs(cos(PoloidalAngle)),abs(sin(PoloidalAngle)))*&
             wall_radius(PoloidalAngle))
        R=R*StretchCoef
        Z=Z*StretchCoef
     end if
     GenOut_D(z_)=Z
     GenOut_D(R_)=R+rTorusLarge
  case default
     call stop_mpi('Unknown TypeGeometry='//TypeGeometry)
  end select
end subroutine xyz_to_gen
!=============================================================================
real function wall_radius(PoloidalAngle)
  !This functions calculates the poloidal radius of the vacuum chamber wall
  !as a function of the poloidal angle
  use ModCovariant
  implicit none
  real::Residual,PoloidalAngle
  integer::iPoint
  real,parameter::dAngle=cTwoPi/nToroidalBoundaryPoints
  if(.not.IsInitializedTorusGeometry)then
     wall_radius=rTorusSmall
  else
     iPoint=int(PoloidalAngle/dAngle)
     if(iPoint==nToroidalBoundaryPoints)iPoint=0
     !Use linear interpolation:
     Residual=PoloidalAngle-iPoint*dAngle
     wall_radius=TorusSurface_I(iPoint)*(cOne-Residual)+&
          TorusSurface_I(iPoint+1)*Residual
  end if
end function wall_radius

!=============================================================================
subroutine fix_pole(iBlock)
  use ModSize
  use ModMain,ONLY:R_,Theta_
  use ModGeometry,ONLY:dx_BLK,dz_BLK,XyzStart_BLK
  use ModGeometry,ONLY:IsPole_B=>DoFixExtraBoundary_B
  use ModCovariant
  implicit none
  integer, intent(in)::iBlock
  real::dR,dTheta,Area2
  integer::k

  !-------------------------------------------------------------!
  FaceArea2MinI_B(iBlock)=cZero
  FaceArea2MinJ_B(iBlock)=cZero
  FaceArea2MinK_B(iBlock)=cZero

  if(.not.is_axial_geometry())return
  if(index(TypeGeometry,'spherical')>0)then
     !calculate a minimal non-zero face area squared along theta direction
     dTheta=dz_BLK(iBlock)
     IsPole_B(iBlock)= XyzStart_BLK(Theta_,iBlock)-dTheta < -cHalfPi &
          .or.XyzStart_BLK(Theta_,iBlock)+nK*dTheta > cHalfPi
     if(IsPole_B(iBlock))then
        Area2=sum(FaceAreaK_DFB(:,1,1,2,iBlock)**2)
        !minimal Area2 is calculated over cell faces non-intersecting the pole
        !no dependence on Phi_

        Area2=min(Area2,sum(FaceAreaK_DFB(:,1,1,nK,iBlock)**2))

        FaceArea2MinK_B(iBlock)=Area2
     end if
  elseif(index(TypeGeometry,'cylindrical')>0.or.&
       index(TypeGeometry,'axial')>0)then
     dR=dx_BLK(iBlock)
     IsPole_B(iBlock)= XyzStart_BLK(R_,iBlock)-dR < cZero 
     if(IsPole_B(iBlock))then
        Area2=sum(FaceAreaI_DFB(:,2,1,1,iBlock)**2)
        do k=1,nK 
           Area2=min(Area2,sum(FaceAreaI_DFB(:,2,1,k,iBlock)**2))
        end do
        FaceArea2MinI_B(iBlock)=Area2
     end if
  end if
end subroutine fix_pole
!=======================================================================
subroutine fix_geometry_at_reschange(iBlock)
  use ModCovariant
  use ModNodes,ONLY:NodeX_NB,NodeY_NB,NodeZ_NB
  use ModGeometry,ONLY: vInv_CB,XyzStart_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModMain,ONLY:x_,y_,z_
  implicit none
  integer,intent(in)::iBlock
  real,dimension(nDim,1:nI+1,1:nJ+1,1:nK+1)::XyzNode_DN
  real,dimension(nDim)::DXyzRef_D
  !It is easy to see that volume can be represented as follows:
  !\int{dV}=\int{{\bf r}\cdot d{\bf S}/nDim. The following array store
  !the dot products of face area vectors by the raduis vector of
  !the "face center"
  real,dimension(1:nI+1,nJ,nK)::RDotFaceAreaI_F 
  real,dimension(1:nI,nJ+1,nK)::RDotFaceAreaJ_F
  real,dimension(1:nI,nJ,nK+1)::RDotFaceAreaK_F

  real,dimension(nDim,1:nI+1,nJ,nK)::FaceCenterI_DF
  real,dimension(nDim,1:nI,nJ+1,nK)::FaceCenterJ_DF
  real,dimension(nDim,1:nI,nJ,nK+1)::FaceCenterK_DF
  real,dimension(nDim)::FaceCenterStart_D,StartNode_D

  character(LEN=*),parameter::NameSub='fix_geometry_at_reschange'

  integer::i,j,k,iDim

  DXyzRef_D(x_)=dx_BLK(iBlock)*cHalf
  DXyzRef_D(y_)=dy_BLK(iBlock)*cHalf
  DXyzRef_D(z_)=dz_BLK(iBlock)*cHalf

  do k=1,nK+1; do j=1,nJ+1; do i=1,nI+1
     XyzNode_DN(x_,i,j,k)=NodeX_NB(i,j,k,iBlock)
     XyzNode_DN(y_,i,j,k)=NodeY_NB(i,j,k,iBlock)
     XyzNode_DN(z_,i,j,k)=NodeZ_NB(i,j,k,iBlock)
  end do; end do; end do

  if(any(OldLevel_IIIB(:,:,:,iBlock)==-1.and.IsNotCorner_III))then
     !Face areas are modified, recover original ones
     call get_face_area_i(XyzNode_DN,&
          1,nI+1,1,nJ,1,nK,&
          FaceAreaI_DFB(:,:,:,:,iBlock))
     call get_face_area_j(XyzNode_DN,&
          1,nI,1,nJ+1,1,nK,&
          FaceAreaJ_DFB(:,:,:,:,iBlock))
     call get_face_area_k(XyzNode_DN,&
          1,nI,1,nJ,1,nK+1,&
          FaceAreaK_DFB(:,:,:,:,iBlock))
  end if


  !\
  ! Face area vector dot product by the face 
  ! center radius-vector. FACE I
  !/
  FaceCenterStart_D=XyzStart_BLK(:,iBlock)
  FaceCenterStart_D(x_)=FaceCenterStart_D(x_)-dx_BLK(iBlock)*cHalf
  call gen_to_xyz_arr(FaceCenterStart_D,&
       dx_BLK(iBlock),dy_BLK(iBlock),dz_BLK(iBlock),&
       1,1+nI,1,nJ,1,nK,&
       FaceCenterI_DF(x_,:,:,:),&
       FaceCenterI_DF(y_,:,:,:),&
       FaceCenterI_DF(z_,:,:,:))
  do k=1,nK; do j=1,nJ; do i=1,nI+1
     RDotFaceAreaI_F(i,j,k)=dot_product(&
          FaceCenterI_DF(:,i,j,k),&
          FaceAreaI_DFB( :,i,j,k,iBlock))
  end do; end do; end do
  !\
  ! Face area vector dot product by the face 
  ! center radius-vector. FACE J
  !/

  FaceCenterStart_D=XyzStart_BLK(:,iBlock)
  FaceCenterStart_D(y_)=FaceCenterStart_D(y_)-dy_BLK(iBlock)*cHalf
  call gen_to_xyz_arr(FaceCenterStart_D,&
       dx_BLK(iBlock),dy_BLK(iBlock),dz_BLK(iBlock),&
       1,nI,1,nJ+1,1,nK,&
       FaceCenterJ_DF(x_,:,:,:),&
       FaceCenterJ_DF(y_,:,:,:),&
       FaceCenterJ_DF(z_,:,:,:))
  do k=1,nK; do j=1,nJ+1; do i=1,nI
     RDotFaceAreaJ_F(i,j,k)=dot_product(&
          FaceCenterJ_DF(:,i,j,k),&
          FaceAreaJ_DFB( :,i,j,k,iBlock))
  end do; end do; end do

  !\
  ! Face area vector dot product by the face 
  ! center radius-vector. FACE K
  !/

  FaceCenterStart_D=XyzStart_BLK(:,iBlock)
  FaceCenterStart_D(z_)=FaceCenterStart_D(z_)-dz_BLK(iBlock)*cHalf
  call gen_to_xyz_arr(FaceCenterStart_D,&
       dx_BLK(iBlock),dy_BLK(iBlock),dz_BLK(iBlock),&
       1,nI,1,nJ,1,nK+1,&
       FaceCenterK_DF(x_,:,:,:),&
       FaceCenterK_DF(y_,:,:,:),&
       FaceCenterK_DF(z_,:,:,:))
  do k=1,nK+1; do j=1,nJ; do i=1,nI
     RDotFaceAreaK_F(i,j,k)=dot_product(&
          FaceCenterK_DF(:,i,j,k),&
          FaceAreaK_DFB( :,i,j,k,iBlock))
  end do; end do; end do




  !Fix faces of i direction
  do i=1,nI+1
     if((i==   1.and.BLKneighborLEV(-1,0,0,iBlock)==-1).or.&
          (i==nI+1.and.BLKneighborLEV(+1,0,0,iBlock)==-1))then
        call refine_face_i(i)
     else
        if((i==1   .and.BLKneighborLEV(-1,-1, 0,iBlock)==-1).or.&
             (i==1+nI.and.BLKneighborLEV(+1,-1, 0,iBlock)==-1).or.&
             BLKneighborLEV( 0,-1, 0,iBlock)==-1)&
             call refine_face_i_edge_j_minus(i)
        if((i==1   .and.BLKneighborLEV(-1,+1, 0,iBlock)==-1).or.&
             (i==1+nI.and.BLKneighborLEV(+1,+1, 0,iBlock)==-1).or.&
             BLKneighborLEV( 0,+1, 0,iBlock)==-1)&
             call refine_face_i_edge_j_plus(i)
        if((i==1   .and.BLKneighborLEV(-1, 0,-1,iBlock)==-1).or.&
             (i==1+nI.and.BLKneighborLEV(+1, 0,-1,iBlock)==-1).or.&
             BLKneighborLEV( 0, 0,-1,iBlock)==-1)&
             call refine_face_i_edge_k_minus(i)
        if((i==1   .and.BLKneighborLEV(-1, 0,+1,iBlock)==-1).or.&
             (i==1+nI.and.BLKneighborLEV(+1, 0,+1,iBlock)==-1).or.&
             BLKneighborLEV( 0, 0,+1,iBlock)==-1)&
             call refine_face_i_edge_k_plus(i)
     end if
  end do

  !Fix faces of J direction
  do j=1,nJ+1
     if((j==   1.and.BLKneighborLEV(0,-1,0,iBlock)==-1).or.&
          (j==nJ+1.and.BLKneighborLEV(0,+1,0,iBlock)==-1))then
        call refine_face_j(j)
     else
        if((j==1   .and.BLKneighborLEV(-1,-1, 0,iBlock)==-1).or.&
             (j==1+nJ.and.BLKneighborLEV(-1,+1, 0,iBlock)==-1).or.&
             BLKneighborLEV(-1, 0, 0,iBlock)==-1)&
             call refine_face_j_edge_i_minus(j)
        if((j==1   .and.BLKneighborLEV(+1,-1, 0,iBlock)==-1).or.&
             (j==1+nJ.and.BLKneighborLEV(+1,+1, 0,iBlock)==-1).or.&
             BLKneighborLEV(+1, 0, 0,iBlock)==-1)&
             call refine_face_j_edge_i_plus(j)
        if((j==1   .and.BLKneighborLEV( 0,-1,-1,iBlock)==-1).or.&
             (j==1+nJ.and.BLKneighborLEV( 0,+1,-1,iBlock)==-1).or.&
             BLKneighborLEV( 0, 0,-1,iBlock)==-1)&
             call refine_face_j_edge_k_minus(j)
        if((j==1   .and.BLKneighborLEV( 0,-1,+1,iBlock)==-1).or.&
             (j==1+nJ.and.BLKneighborLEV( 0,+1,+1,iBlock)==-1).or.&
             BLKneighborLEV( 0, 0,+1,iBlock)==-1)&
             call refine_face_j_edge_k_plus(j)
     end if
  end do
  !Fix faces of k direction
  do k=1,nK+1
     if((k==   1.and.BLKneighborLEV(0,0,-1,iBlock)==-1).or.&
          (k==nK+1.and.BLKneighborLEV(0,0,+1,iBlock)==-1))then
        call refine_face_k(k)
     else
        if((k==1   .and.BLKneighborLEV(-1, 0,-1,iBlock)==-1).or.&
             (k==1+nK.and.BLKneighborLEV(-1, 0,+1,iBlock)==-1).or.&
             BLKneighborLEV(-1, 0, 0,iBlock)==-1)&
             call refine_face_k_edge_i_minus(k)
        if((k==1   .and.BLKneighborLEV(+1, 0,-1,iBlock)==-1).or.&
             (k==1+nK.and.BLKneighborLEV(+1, 0,+1,iBlock)==-1).or.&
             BLKneighborLEV(+1, 0, 0,iBlock)==-1)&
             call refine_face_k_edge_i_plus(k)
        if((k==1   .and.BLKneighborLEV( 0,-1,-1,iBlock)==-1).or.&
             (k==1+nK.and.BLKneighborLEV( 0,-1,+1,iBlock)==-1).or.&
             BLKneighborLEV( 0,-1, 0,iBlock)==-1)&
             call refine_face_k_edge_j_minus(k)
        if((k==1   .and.BLKneighborLEV( 0,+1,-1,iBlock)==-1).or.&
             (k==1+nK.and.BLKneighborLEV( 0,+1,+1,iBlock)==-1).or.&
             BLKneighborLEV( 0,+1, 0,iBlock)==-1)&
             call refine_face_k_edge_j_plus(k)
     end if
  end do


  !Calculate Volume (inverse)
  vInv_CB(:,:,:,iBlock)=nDim/(&
       RDotFaceAreaI_F(2:nI+1,:,:)-RDotFaceAreaI_F(1:nI,:,:)+&
       RDotFaceAreaJ_F(:,2:nJ+1,:)-RDotFaceAreaJ_F(:,1:nJ,:)+&
       RDotFaceAreaK_F(:,:,2:nK+1)-RDotFaceAreaK_F(:,:,1:nK) )

  !Save level of refinement
  OldLevel_IIIB(:,:,:,iBlock)=BLKneighborLEV(:,:,:,iBlock)
  call fix_pole(iBlock)
  call test_block_geometry(iBlock,NameSub)
contains
  !--------------------------------FACE I----------------------------------!
  !Fix face area vectors along I direction
  subroutine refine_face_i(iFace)
    integer,intent(in)::iFace
    real,dimension(nDim,1,2*nJ,2*nK)::RefFaceAreaI_DF,RefFaceCenterI_DF
    real,dimension(nDim,1,2*nJ+1,2*nK+1)::RefNodesI_DN
    StartNode_D=XyzStart_BLK(:,iBlock)-DXyzRef_D
    StartNode_D(x_)=StartNode_D(x_)+(iFace-1)*dx_BLK(iBlock)
    call gen_to_xyz_arr(StartNode_D,&
         DXyzRef_D(x_),DXyzRef_D(y_),DXyzRef_D(z_),&
         1,1,1,2*nJ+1,1,2*nK+1,&
         RefNodesI_DN(x_,:,:,:),&
         RefNodesI_DN(y_,:,:,:),&
         RefNodesI_DN(z_,:,:,:))
    call get_face_area_i(RefNodesI_DN,&
         1,1,1,2*nJ,1,2*nK,&
         RefFaceAreaI_DF)
    FaceCenterStart_D=StartNode_D
    FaceCenterStart_D(y_)=FaceCenterStart_D(y_)+DXyzRef_D(y_)*cHalf
    FaceCenterStart_D(z_)=FaceCenterStart_D(z_)+DXyzRef_D(z_)*cHalf
    call gen_to_xyz_arr(FaceCenterStart_D,&
         DXyzRef_D(x_),DXyzRef_D(y_),DXyzRef_D(z_),&
         1,1,1,2*nJ,1,2*nK,&
         RefFaceCenterI_DF(x_,:,:,:),&
         RefFaceCenterI_DF(y_,:,:,:),&
         RefFaceCenterI_DF(z_,:,:,:))
    do k=1,nK;do j=1,nJ
       do iDim=1,nDim
          FaceAreaI_DFB(iDim,iFace,j,k,iBlock)=sum(&
               RefFaceAreaI_DF(iDim,1,2*j-1:2*j,2*k-1:2*k))
       end do
       RDotFaceAreaI_F(iFace,j,k)=sum(&
            RefFaceAreaI_DF(:,1,2*j-1:2*j,2*k-1:2*k)*&
            RefFaceCenterI_DF(:,1,2*j-1:2*j,2*k-1:2*k))
    end do;end do
  end subroutine refine_face_i
  !--------------------------------EDGES OF FACE I------------------------!
  !          Fix edges. Start from the following formula:
  ! FaceAreaI=cHalf*[&
  ! cross_product(XyzNodes_DIII(:,i,j  ,k),XyzNodes_DIII(:,i,j+1,k  ))+&
  ! cross_product(XyzNodes_DIII(:,i,j+1,k),XyzNodes_DIII(:,i,j+1,k+1))+&
  ! cross_product(XyzNodes_DIII(:,i,j+1,k+1),XyzNodes_DIII(:,i,j  ,k+1))+&
  ! cross_product(XyzNodes_DIII(:,i,j  ,k+1),XyzNodes_DIII(:,i,j  ,k))]
  subroutine refine_face_i_edge_j_minus(iFace)
    integer,intent(in)::iFace
    real,dimension(nDim,1,1,2*nK+1)::RefNodesI_DN
    StartNode_D=XyzStart_BLK(:,iBlock)-DXyzRef_D
    StartNode_D(x_)= StartNode_D(x_)+(iFace-1)*dx_BLK(iBlock)
    call gen_to_xyz_arr(StartNode_D,&
         DXyzRef_D(x_),DXyzRef_D(y_),DXyzRef_D(z_),&
         1,1,1,1,1,2*nK+1,&
         RefNodesI_DN(x_,:,:,:),&
         RefNodesI_DN(y_,:,:,:),&
         RefNodesI_DN(z_,:,:,:))
    do k=1,nK
       FaceAreaI_DFB(:,iFace,1,k,iBlock)=&
            FaceAreaI_DFB(:,iFace,1,k,iBlock)+cHalf*(&
            -cross_product(RefNodesI_DN(:,1,1,2*k+1),RefNodesI_DN(:,1,1,2*k-1))&
            +cross_product(RefNodesI_DN(:,1,1,2*k+1),RefNodesI_DN(:,1,1,2*k  ))&
            +cross_product(RefNodesI_DN(:,1,1,2*k),RefNodesI_DN(:,1,1,2*k-1)))
       RDotFaceAreaI_F(iFace,1,k)=dot_product(&
            FaceCenterI_DF(:,iFace,1,k),&
            FaceAreaI_DFB( :,iFace,1,k,iBlock))
    end do
  end subroutine refine_face_i_edge_j_minus
  subroutine refine_face_i_edge_j_plus(iFace)
    integer,intent(in)::iFace
    real,dimension(nDim,1,1,2*nK+1)::RefNodesI_DN
    StartNode_D=XyzStart_BLK(:,iBlock)-DXyzRef_D
    StartNode_D(x_)= StartNode_D(x_)+(iFace-1)*dx_BLK(iBlock)
    StartNode_D(y_)= StartNode_D(y_)+nJ*dy_BLK(iBlock)
    call gen_to_xyz_arr(StartNode_D,&
         DXyzRef_D(x_),DXyzRef_D(y_),DXyzRef_D(z_),&
         1,1,1,1,1,2*nK+1,&
         RefNodesI_DN(x_,:,:,:),&
         RefNodesI_DN(y_,:,:,:),&
         RefNodesI_DN(z_,:,:,:))
    do k=1,nK
       FaceAreaI_DFB(:,iFace,nJ,k,iBlock)=&
            FaceAreaI_DFB(:,iFace,nJ,k,iBlock)+cHalf*(&
            -cross_product(RefNodesI_DN(:,1,1,2*k-1),RefNodesI_DN(:,1,1,2*k+1))&
            +cross_product(RefNodesI_DN(:,1,1,2*k-1),RefNodesI_DN(:,1,1,2*k  ))&
            +cross_product(RefNodesI_DN(:,1,1,2*k),RefNodesI_DN(:,1,1,2*k+1)))
       RDotFaceAreaI_F(iFace,nJ,k)=dot_product(&
            FaceCenterI_DF(:,iFace,nJ,k),&
            FaceAreaI_DFB( :,iFace,nJ,k,iBlock))
    end do
  end subroutine refine_face_i_edge_j_plus
  subroutine refine_face_i_edge_k_minus(iFace)
    integer,intent(in)::iFace
    real,dimension(nDim,1,2*nJ+1,1)::RefNodesI_DN
    StartNode_D=XyzStart_BLK(:,iBlock)-DXyzRef_D
    StartNode_D(x_)= StartNode_D(x_)+(iFace-1)*dx_BLK(iBlock)
    call gen_to_xyz_arr(StartNode_D,&
         DXyzRef_D(x_),DXyzRef_D(y_),DXyzRef_D(z_),&
         1,1,1,2*nJ+1,1,1,&
         RefNodesI_DN(x_,:,:,:),&
         RefNodesI_DN(y_,:,:,:),&
         RefNodesI_DN(z_,:,:,:))
    do j=1,nJ
       FaceAreaI_DFB(:,iFace,j,1,iBlock)=&
            FaceAreaI_DFB(:,iFace,j,1,iBlock)+cHalf*(&
            -cross_product(RefNodesI_DN(:,1,2*j-1,1),RefNodesI_DN(:,1,2*j+1,1))&
            +cross_product(RefNodesI_DN(:,1,2*j-1,1),RefNodesI_DN(:,1,2*j,1  ))&
            +cross_product(RefNodesI_DN(:,1,2*j  ,1),RefNodesI_DN(:,1,2*j+1,1)))
       RDotFaceAreaI_F(iFace,j,1)=dot_product(&
            FaceCenterI_DF(:,iFace,j,1),&
            FaceAreaI_DFB( :,iFace,j,1,iBlock))
    end do
  end subroutine refine_face_i_edge_k_minus
  subroutine refine_face_i_edge_k_plus(iFace)
    integer,intent(in)::iFace
    real,dimension(nDim,1,2*nJ+1,1)::RefNodesI_DN
    StartNode_D=XyzStart_BLK(:,iBlock)-DXyzRef_D
    StartNode_D(x_)= StartNode_D(x_)+(iFace-1)*dx_BLK(iBlock)
    StartNode_D(z_)= StartNode_D(z_)+nK*dz_BLK(iBlock)
    call gen_to_xyz_arr(StartNode_D,&
         DXyzRef_D(x_),DXyzRef_D(y_),DXyzRef_D(z_),&
         1,1,1,2*nJ+1,1,1,&
         RefNodesI_DN(x_,:,:,:),&
         RefNodesI_DN(y_,:,:,:),&
         RefNodesI_DN(z_,:,:,:))
    do j=1,nJ
       FaceAreaI_DFB(:,iFace,j,nK,iBlock)=&
            FaceAreaI_DFB(:,iFace,j,nK,iBlock)+cHalf*(&
            -cross_product(RefNodesI_DN(:,1,2*j+1,1),RefNodesI_DN(:,1,2*j-1,1))&
            +cross_product(RefNodesI_DN(:,1,2*j+1,1),RefNodesI_DN(:,1,2*j  ,1))&
            +cross_product(RefNodesI_DN(:,1,2*j  ,1),RefNodesI_DN(:,1,2*j-1,1)))
       RDotFaceAreaI_F(iFace,j,nK)=dot_product(&
            FaceCenterI_DF(:,iFace,j,nK),&
            FaceAreaI_DFB( :,iFace,j,nK,iBlock))
    end do
  end subroutine refine_face_i_edge_k_plus
  !--------------------------------FACE J----------------------------------!
  !Fix face area vectors along J direction
  subroutine refine_face_j(jFace)
    integer,intent(in)::jFace
    real,dimension(nDim,2*nI,1,2*nK)::RefFaceAreaJ_DF,RefFaceCenterJ_DF
    real,dimension(nDim,2*nI+1,1,2*nK+1)::RefNodesJ_DN
    StartNode_D=XyzStart_BLK(:,iBlock)-DXyzRef_D
    StartNode_D(y_)=StartNode_D(y_)+(jFace-1)*dy_BLK(iBlock)
    call gen_to_xyz_arr(StartNode_D,&
         DXyzRef_D(x_),DXyzRef_D(y_),DXyzRef_D(z_),&
         1,2*nI+1,1,1,1,2*nK+1,&
         RefNodesJ_DN(x_,:,:,:),&
         RefNodesJ_DN(y_,:,:,:),&
         RefNodesJ_DN(z_,:,:,:))
    call get_face_area_j(RefNodesJ_DN,&
         1,2*nI,1,1,1,2*nK,&
         RefFaceAreaJ_DF)
    FaceCenterStart_D=StartNode_D
    FaceCenterStart_D(x_)=FaceCenterStart_D(x_)+DXyzRef_D(x_)*cHalf
    FaceCenterStart_D(z_)=FaceCenterStart_D(z_)+DXyzRef_D(z_)*cHalf
    call gen_to_xyz_arr(FaceCenterStart_D,&
         DXyzRef_D(x_),DXyzRef_D(y_),DXyzRef_D(z_),&
         1,2*nI,1,1,1,2*nK,&
         RefFaceCenterJ_DF(x_,:,:,:),&
         RefFaceCenterJ_DF(y_,:,:,:),&
         RefFaceCenterJ_DF(z_,:,:,:))
    do k=1,nK;do i=1,nI
       do iDim=1,nDim
          FaceAreaJ_DFB(iDim,i,jFace,k,iBlock)=sum(&
               RefFaceAreaJ_DF(iDim,2*i-1:2*i,1,2*k-1:2*k))
       end do
       RDotFaceAreaJ_F(i,jFace,k)=sum(&
            RefFaceAreaJ_DF(:,2*i-1:2*i,1,2*k-1:2*k)*&
            RefFaceCenterJ_DF(:,2*i-1:2*i,1,2*k-1:2*k))
    end do;end do
  end subroutine refine_face_j
  !--------------------------------EDGES OF FACE J------------------------!
  !          Fix edges. Start from the following formula:
  ! FaceAreaJ=cHalf*[&
  ! cross_product(XyzNodes_DIII(:,i,j  ,k),XyzNodes_DIII(:,i,j,k+1))+&
  ! cross_product(XyzNodes_DIII(:,i,j,k+1),XyzNodes_DIII(:,i+1,j,k+1))+&
  ! cross_product(XyzNodes_DIII(:,i+1,j,k+1),XyzNodes_DIII(:,i+1,j  ,k))+
  ! cross_product(XyzNodes_DIII(:,i+1,j  ,k),XyzNodes_DIII(:,i,j  ,k))]
  subroutine refine_face_j_edge_i_minus(jFace)
    integer,intent(in)::jFace
    real,dimension(nDim,1,1,2*nK+1)::RefNodesJ_DN
    StartNode_D=XyzStart_BLK(:,iBlock)-DXyzRef_D
    StartNode_D(y_)= StartNode_D(y_)+(jFace-1)*dy_BLK(iBlock)
    call gen_to_xyz_arr(StartNode_D,&
         DXyzRef_D(x_),DXyzRef_D(y_),DXyzRef_D(z_),&
         1,1,1,1,1,2*nK+1,&
         RefNodesJ_DN(x_,:,:,:),&
         RefNodesJ_DN(y_,:,:,:),&
         RefNodesJ_DN(z_,:,:,:))
    do k=1,nK
       FaceAreaJ_DFB(:,1,jFace,k,iBlock)=&
            FaceAreaJ_DFB(:,1,jFace,k,iBlock)+cHalf*(&
            -cross_product(RefNodesJ_DN(:,1,1,2*k-1),RefNodesJ_DN(:,1,1,2*k+1))&
            +cross_product(RefNodesJ_DN(:,1,1,2*k-1),RefNodesJ_DN(:,1,1,2*k  ))&
            +cross_product(RefNodesJ_DN(:,1,1,2*k),  RefNodesJ_DN(:,1,1,2*k+1)))
       RDotFaceAreaJ_F(1,jFace,k)=dot_product(&
            FaceCenterJ_DF(:,1,jFace,k),&
            FaceAreaJ_DFB( :,1,jFace,k,iBlock))
    end do
  end subroutine refine_face_j_edge_i_minus
  subroutine refine_face_j_edge_i_plus(jFace)
    integer,intent(in)::jFace
    real,dimension(nDim,1,1,2*nK+1)::RefNodesJ_DN
    StartNode_D=XyzStart_BLK(:,iBlock)-DXyzRef_D
    StartNode_D(y_)= StartNode_D(y_)+(jFace-1)*dy_BLK(iBlock)
    StartNode_D(x_)= StartNode_D(x_)+nI*dx_BLK(iBlock)
    call gen_to_xyz_arr(StartNode_D,&
         DXyzRef_D(x_),DXyzRef_D(y_),DXyzRef_D(z_),&
         1,1,1,1,1,2*nK+1,&
         RefNodesJ_DN(x_,:,:,:),&
         RefNodesJ_DN(y_,:,:,:),&
         RefNodesJ_DN(z_,:,:,:))
    do k=1,nK
       FaceAreaJ_DFB(:,nI,jFace,k,iBlock)=&
            FaceAreaJ_DFB(:,nI,jFace,k,iBlock)+cHalf*(&
            -cross_product(RefNodesJ_DN(:,1,1,2*k+1),RefNodesJ_DN(:,1,1,2*k-1))&
            +cross_product(RefNodesJ_DN(:,1,1,2*k+1),RefNodesJ_DN(:,1,1,2*k  ))&
            +cross_product(RefNodesJ_DN(:,1,1,2*k)  ,RefNodesJ_DN(:,1,1,2*k-1)))
       RDotFaceAreaJ_F(nI,jFace,k)=dot_product(&
            FaceCenterJ_DF(:,nI,jFace,k),&
            FaceAreaJ_DFB( :,nI,jFace,k,iBlock))
    end do
  end subroutine refine_face_j_edge_i_plus
  subroutine refine_face_j_edge_k_minus(jFace)
    integer,intent(in)::jFace
    real,dimension(nDim,2*nI+1,1,1)::RefNodesJ_DN
    StartNode_D=XyzStart_BLK(:,iBlock)-DXyzRef_D
    StartNode_D(y_)= StartNode_D(y_)+(jFace-1)*dy_BLK(iBlock)
    call gen_to_xyz_arr(StartNode_D,&
         DXyzRef_D(x_),DXyzRef_D(y_),DXyzRef_D(z_),&
         1,2*nI+1,1,1,1,1,&
         RefNodesJ_DN(x_,:,:,:),&
         RefNodesJ_DN(y_,:,:,:),&
         RefNodesJ_DN(z_,:,:,:))
    do i=1,nI
       FaceAreaJ_DFB(:,i,jFace,1,iBlock)=&
            FaceAreaJ_DFB(:,i,jFace,1,iBlock)+cHalf*(&
            -cross_product(RefNodesJ_DN(:,2*i+1,1,1),RefNodesJ_DN(:,2*i-1,1,1))&
            +cross_product(RefNodesJ_DN(:,2*i+1,1,1),RefNodesJ_DN(:,2*i  ,1,1  ))&
            +cross_product(RefNodesJ_DN(:,2*i  ,1,1),RefNodesJ_DN(:,2*i-1,1,1)))
       RDotFaceAreaJ_F(i,jFace,1)=dot_product(&
            FaceCenterJ_DF(:,i,jFace,1),&
            FaceAreaJ_DFB( :,i,jFace,1,iBlock))
    end do
  end subroutine refine_face_j_edge_k_minus
  subroutine refine_face_j_edge_k_plus(jFace)
    integer,intent(in)::jFace
    real,dimension(nDim,2*nI+1,1,1)::RefNodesJ_DN
    StartNode_D=XyzStart_BLK(:,iBlock)-DXyzRef_D
    StartNode_D(y_)= StartNode_D(y_)+(jFace-1)*dy_BLK(iBlock)
    StartNode_D(z_)= StartNode_D(z_)+nK*dz_BLK(iBlock)
    call gen_to_xyz_arr(StartNode_D,&
         DXyzRef_D(x_),DXyzRef_D(y_),DXyzRef_D(z_),&
         1,2*nI+1,1,1,1,1,&
         RefNodesJ_DN(x_,:,:,:),&
         RefNodesJ_DN(y_,:,:,:),&
         RefNodesJ_DN(z_,:,:,:))
    do i=1,nI
       FaceAreaJ_DFB(:,i,jFace,nK,iBlock)=&
            FaceAreaJ_DFB(:,i,jFace,nK,iBlock)+cHalf*(&
            -cross_product(RefNodesJ_DN(:,2*i-1,1,1),RefNodesJ_DN(:,2*i+1,1,1))&
            +cross_product(RefNodesJ_DN(:,2*i-1,1,1),RefNodesJ_DN(:,2*i  ,1,1))&
            +cross_product(RefNodesJ_DN(:,2*i  ,1,1),RefNodesJ_DN(:,2*i+1,1,1)))
       RDotFaceAreaJ_F(i,jFace,nK)=dot_product(&
            FaceCenterJ_DF(:,i,jFace,nK),&
            FaceAreaJ_DFB( :,i,jFace,nK,iBlock))
    end do
  end subroutine refine_face_j_edge_k_plus
  !--------------------------------FACE K----------------------------------!
  !Fix face area vectors along K direction
  subroutine refine_face_k(kFace)
    integer,intent(in)::kFace
    real,dimension(nDim,2*nI,2*nJ,1)::RefFaceAreaK_DF,RefFaceCenterK_DF
    real,dimension(nDim,2*nI+1,2*nJ+1,1)::RefNodesK_DN
    StartNode_D=XyzStart_BLK(:,iBlock)-DXyzRef_D
    StartNode_D(z_)=StartNode_D(z_)+(kFace-1)*dz_BLK(iBlock)
    call gen_to_xyz_arr(StartNode_D,&
         DXyzRef_D(x_),DXyzRef_D(y_),DXyzRef_D(z_),&
         1,2*nI+1,1,2*nJ+1,1,1,&
         RefNodesK_DN(x_,:,:,:),&
         RefNodesK_DN(y_,:,:,:),&
         RefNodesK_DN(z_,:,:,:))
    call get_face_area_k(RefNodesK_DN,&
         1,2*nI,1,2*nJ,1,1,&
         RefFaceAreaK_DF)
    FaceCenterStart_D=StartNode_D
    FaceCenterStart_D(x_)=FaceCenterStart_D(x_)+DXyzRef_D(x_)*cHalf
    FaceCenterStart_D(y_)=FaceCenterStart_D(y_)+DXyzRef_D(y_)*cHalf
    call gen_to_xyz_arr(FaceCenterStart_D,&
         DXyzRef_D(x_),DXyzRef_D(y_),DXyzRef_D(z_),&
         1,2*nI,1,2*nJ,1,1,&
         RefFaceCenterK_DF(x_,:,:,:),&
         RefFaceCenterK_DF(y_,:,:,:),&
         RefFaceCenterK_DF(z_,:,:,:))
    do j=1,nJ;do i=1,nI
       do iDim=1,nDim
          FaceAreaK_DFB(iDim,i,j,kFace,iBlock)=sum(&
               RefFaceAreaK_DF(iDim,2*i-1:2*i,2*j-1:2*j,1))
       end do
       RDotFaceAreaK_F(i,j,kFace)=sum(&
            RefFaceAreaK_DF(:,2*i-1:2*i,2*j-1:2*j,1)*&
            RefFaceCenterK_DF(:,2*i-1:2*i,2*j-1:2*j,1))
    end do;end do
  end subroutine refine_face_k
  !--------------------------------EDGES OF FACE K------------------------!
  !          Fix edges. Start from the following formula:
  ! FaceAreaK=cHalf*[&
  ! cross_product(XyzNodes_DIII(:,i,j  ,k),XyzNodes_DIII(:,i+1,j,k))+&
  ! cross_product(XyzNodes_DIII(:,i+1,j,k),XyzNodes_DIII(:,i+1,j+1,k))+&
  ! cross_product(XyzNodes_DIII(:,i+1,j+1,k),XyzNodes_DIII(:,i,j+1  ,k))+&
  ! cross_product(XyzNodes_DIII(:,i,j+1  ,k),XyzNodes_DIII(:,i,j  ,k))]

  subroutine refine_face_k_edge_i_minus(kFace)
    integer,intent(in)::kFace
    real,dimension(nDim,1,2*nJ+1,1)::RefNodesK_DN
    StartNode_D=XyzStart_BLK(:,iBlock)-DXyzRef_D
    StartNode_D(z_)= StartNode_D(z_)+(kFace-1)*dz_BLK(iBlock)
    call gen_to_xyz_arr(StartNode_D,&
         DXyzRef_D(x_),DXyzRef_D(y_),DXyzRef_D(z_),&
         1,1,1,2*nJ+1,1,1,&
         RefNodesK_DN(x_,:,:,:),&
         RefNodesK_DN(y_,:,:,:),&
         RefNodesK_DN(z_,:,:,:))
    do j=1,nJ
       FaceAreaK_DFB(:,1,j,kFace,iBlock)=&
            FaceAreaK_DFB(:,1,j,kFace,iBlock)+cHalf*(&
            -cross_product(RefNodesK_DN(:,1,2*j+1,1),RefNodesK_DN(:,1,2*j-1,1))&
            +cross_product(RefNodesK_DN(:,1,2*j+1,1),RefNodesK_DN(:,1,2*j  ,1))&
            +cross_product(RefNodesK_DN(:,1,2*j  ,1),RefNodesK_DN(:,1,2*j-1,1)))
       RDotFaceAreaK_F(1,j,kFace)=dot_product(&
            FaceCenterK_DF(:,1,j,kFace),&
            FaceAreaK_DFB( :,1,j,kFace,iBlock))
    end do
  end subroutine refine_face_k_edge_i_minus
  subroutine refine_face_k_edge_i_plus(kFace)
    integer,intent(in)::kFace
    real,dimension(nDim,1,2*nJ+1,1)::RefNodesK_DN
    StartNode_D=XyzStart_BLK(:,iBlock)-DXyzRef_D
    StartNode_D(z_)= StartNode_D(z_)+(kFace-1)*dz_BLK(iBlock)
    StartNode_D(x_)= StartNode_D(x_)+nI*dx_BLK(iBlock)
    call gen_to_xyz_arr(StartNode_D,&
         DXyzRef_D(x_),DXyzRef_D(y_),DXyzRef_D(z_),&
         1,1,1,2*nJ+1,1,1,&
         RefNodesK_DN(x_,:,:,:),&
         RefNodesK_DN(y_,:,:,:),&
         RefNodesK_DN(z_,:,:,:))
    do j=1,nJ
       FaceAreaK_DFB(:,nI,j,kFace,iBlock)=&
            FaceAreaK_DFB(:,nI,j,kFace,iBlock)+cHalf*(&
            -cross_product(RefNodesK_DN(:,1,2*j-1,1),RefNodesK_DN(:,1,2*j+1,1))&
            +cross_product(RefNodesK_DN(:,1,2*j-1,1),RefNodesK_DN(:,1,2*j  ,1))&
            +cross_product(RefNodesK_DN(:,1,2*j  ,1),RefNodesK_DN(:,1,2*j+1,1)))
       RDotFaceAreaK_F(nI,j,kFace)=dot_product(&
            FaceCenterK_DF(:,nI,j,kFace),&
            FaceAreaK_DFB( :,nI,j,kFace,iBlock))
    end do
  end subroutine refine_face_k_edge_i_plus
  subroutine refine_face_k_edge_j_minus(kFace)
    integer,intent(in)::kFace
    real,dimension(nDim,2*nI+1,1,1)::RefNodesK_DN
    StartNode_D=XyzStart_BLK(:,iBlock)-DXyzRef_D
    StartNode_D(z_)= StartNode_D(z_)+(kFace-1)*dz_BLK(iBlock)
    call gen_to_xyz_arr(StartNode_D,&
         DXyzRef_D(x_),DXyzRef_D(y_),DXyzRef_D(z_),&
         1,2*nI+1,1,1,1,1,&
         RefNodesK_DN(x_,:,:,:),&
         RefNodesK_DN(y_,:,:,:),&
         RefNodesK_DN(z_,:,:,:))
    do i=1,nI
       FaceAreaK_DFB(:,i,1,kFace,iBlock)=&
            FaceAreaK_DFB(:,i,1,kFace,iBlock)+cHalf*(&
            -cross_product(RefNodesK_DN(:,2*i-1,1,1),RefNodesK_DN(:,2*i+1,1,1))&
            +cross_product(RefNodesK_DN(:,2*i-1,1,1),RefNodesK_DN(:,2*i  ,1,1  ))&
            +cross_product(RefNodesK_DN(:,2*i  ,1,1),RefNodesK_DN(:,2*i+1,1,1)))
       RDotFaceAreaK_F(i,1,kFace)=dot_product(&
            FaceCenterK_DF(:,i,1,kFace),&
            FaceAreaK_DFB( :,i,1,kFace,iBlock))
    end do
  end subroutine refine_face_k_edge_j_minus
  subroutine refine_face_k_edge_j_plus(kFace)
    integer,intent(in)::kFace
    real,dimension(nDim,2*nI+1,1,1)::RefNodesK_DN
    StartNode_D=XyzStart_BLK(:,iBlock)-DXyzRef_D
    StartNode_D(z_)= StartNode_D(z_)+(kFace-1)*dz_BLK(iBlock)
    StartNode_D(y_)= StartNode_D(y_)+nJ*dy_BLK(iBlock)
    call gen_to_xyz_arr(StartNode_D,&
         DXyzRef_D(x_),DXyzRef_D(y_),DXyzRef_D(z_),&
         1,2*nI+1,1,1,1,1,&
         RefNodesK_DN(x_,:,:,:),&
         RefNodesK_DN(y_,:,:,:),&
         RefNodesK_DN(z_,:,:,:))
    do i=1,nI
       FaceAreaK_DFB(:,i,nJ,kFace,iBlock)=&
            FaceAreaK_DFB(:,i,nJ,kFace,iBlock)+cHalf*(&
            -cross_product(RefNodesK_DN(:,2*i+1,1,1),RefNodesK_DN(:,2*i-1,1,1))&
            +cross_product(RefNodesK_DN(:,2*i+1,1,1),RefNodesK_DN(:,2*i  ,1,1))&
            +cross_product(RefNodesK_DN(:,2*i  ,1,1),RefNodesK_DN(:,2*i-1,1,1)))
       RDotFaceAreaK_F(i,nJ,kFace)=dot_product(&
            FaceCenterK_DF(:,i,nJ,kFace),&
            FaceAreaK_DFB( :,i,nJ,kFace,iBlock))
    end do
  end subroutine refine_face_k_edge_j_plus

end subroutine fix_geometry_at_reschange
!=======================================================================
subroutine fix_covariant_geometry(iBLK)
  use ModCovariant
  use ModNodes,ONLY:NodeX_NB,NodeY_NB,NodeZ_NB
  use ModGeometry,ONLY: vInv_CB,XyzStart_BLK,dx_BLK,dy_BLK,dz_BLK, &
       x_BLK,y_BLK,z_BLK
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

  real,dimension(nDim,1:nI+1,nJ,nK)::FaceCenterI_DF
  real,dimension(nDim,1:nI,nJ+1,nK)::FaceCenterJ_DF
  real,dimension(nDim,1:nI,nJ,nK+1)::FaceCenterK_DF
  real,dimension(nDim)::FaceCenterStart_D
  character(LEN=*),parameter::NameSub='fix_covariant_geometry'

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
  if(UseVertexBasedGrid)then
     FaceCenterStart_D=XyzStart_BLK(:,iBLK)
     FaceCenterStart_D(x_)=FaceCenterStart_D(x_)-dx_BLK(iBLK)*cHalf
     call gen_to_xyz_arr(FaceCenterStart_D,&
          dx_BLK(iBLK),dy_BLK(iBLK),dz_BLK(iBLK),&
          1,1+nI,1,nJ,1,nK,&
          FaceCenterI_DF(x_,:,:,:),&
          FaceCenterI_DF(y_,:,:,:),&
          FaceCenterI_DF(z_,:,:,:))
     do k=1,nK; do j=1,nJ; do i=1,nI+1
        RDotFaceAreaI_F(i,j,k)=dot_product(&
             FaceCenterI_DF(:,i,j,k),&
             FaceAreaI_DFB( :,i,j,k,iBLK))
     end do; end do; end do
  else
     do k=1,nK; do j=1,nJ; do i=1,nI+1
        RDotFaceAreaI_F(i,j,k)=0.25*dot_product(&
             XyzNode_DN(:,i,j  ,k  )+ &
             XyzNode_DN(:,i,j+1,k  )+ &
             XyzNode_DN(:,i,j+1,k+1)+ &
             XyzNode_DN(:,i,j  ,k+1), &
             FaceAreaI_DFB(:,i,j,k,iBLK))
     end do; end do; end do
  end if
  !\
  ! Face area vector and its dot product by the face 
  ! center radius-vector. FACE J
  !/
  call get_face_area_j(XyzNode_DN,&
       1,nI,1,nJ+1,1,nK,&
       FaceAreaJ_DFB(:,:,:,:,iBLK))
  if(UseVertexBasedGrid)then
     FaceCenterStart_D=XyzStart_BLK(:,iBLK)
     FaceCenterStart_D(y_)=FaceCenterStart_D(y_)-dy_BLK(iBLK)*cHalf
     call gen_to_xyz_arr(FaceCenterStart_D,&
          dx_BLK(iBLK),dy_BLK(iBLK),dz_BLK(iBLK),&
          1,nI,1,nJ+1,1,nK,&
          FaceCenterJ_DF(x_,:,:,:),&
          FaceCenterJ_DF(y_,:,:,:),&
          FaceCenterJ_DF(z_,:,:,:))
     do k=1,nK; do j=1,nJ+1; do i=1,nI
        RDotFaceAreaJ_F(i,j,k)=dot_product(&
             FaceCenterJ_DF(:,i,j,k),&
             FaceAreaJ_DFB( :,i,j,k,iBLK))
     end do; end do; end do
  else
     do k=1,nK; do j=1,nJ+1; do i=1,nI
        RDotFaceAreaJ_F(i,j,k)=0.25*dot_product(&
             XyzNode_DN(:,i  ,j,k  )+ &
             XyzNode_DN(:,i  ,j,k+1)+ &
             XyzNode_DN(:,i+1,j,k+1)+ &
             XyzNode_DN(:,i+1,j,k  ), &
             FaceAreaJ_DFB(:,i,j,k,iBLK))
     end do; end do; end do 
  end if

  !\
  ! Face area vector and its dot product by the face 
  ! center radius-vector. FACE K
  !/
  call get_face_area_k(XyzNode_DN,&
       1,nI,1,nJ,1,nK+1,&
       FaceAreaK_DFB(:,:,:,:,iBLK))
  if(UseVertexBasedGrid)then
     FaceCenterStart_D=XyzStart_BLK(:,iBLK)
     FaceCenterStart_D(z_)=FaceCenterStart_D(z_)-dz_BLK(iBLK)*cHalf
     call gen_to_xyz_arr(FaceCenterStart_D,&
          dx_BLK(iBLK),dy_BLK(iBLK),dz_BLK(iBLK),&
          1,nI,1,nJ,1,nK+1,&
          FaceCenterK_DF(x_,:,:,:),&
          FaceCenterK_DF(y_,:,:,:),&
          FaceCenterK_DF(z_,:,:,:))
     do k=1,nK+1; do j=1,nJ; do i=1,nI
        RDotFaceAreaK_F(i,j,k)=dot_product(&
             FaceCenterK_DF(:,i,j,k),&
             FaceAreaK_DFB( :,i,j,k,iBLK))
     end do; end do; end do
  else
     do k=1,nK+1; do j=1,nJ; do i=1,nI
        RDotFaceAreaK_F(i,j,k)=0.25*dot_product(&
             XyzNode_DN(:,i  ,j  ,k)+ &
             XyzNode_DN(:,i+1,j  ,k)+ &
             XyzNode_DN(:,i+1,j+1,k)+ &
             XyzNode_DN(:,i  ,j+1,k), &
             FaceAreaK_DFB(:,i,j,k,iBLK))
     end do; end do; end do 
  end if

  !Calculate Volume (inverse)
  vInv_CB(:,:,:,iBLK)=nDim/(&
       RDotFaceAreaI_F(2:nI+1,:,:)-RDotFaceAreaI_F(1:nI,:,:)+&
       RDotFaceAreaJ_F(:,2:nJ+1,:)-RDotFaceAreaJ_F(:,1:nJ,:)+&
       RDotFaceAreaK_F(:,:,2:nK+1)-RDotFaceAreaK_F(:,:,1:nK) )
  call fix_pole(iBLK)
  call test_block_geometry(iBLK,NameSub)
  if(trim(TypeGeometry)=='rz')call fix_rz_geometry(iBLK)
end subroutine fix_covariant_geometry
subroutine test_block_geometry(iBlock,NameSub)
  use ModCovariant
  implicit none
  integer,intent(in)::iBlock
  character(LEN=*),intent(in)::NameSub
  real,dimension(nDim)::FaceArea_D
  integer::i,j,k
  do k=1,nK;do j=1,nJ;do i=1,nI
     FaceArea_D=FaceAreaI_DFB(:,i+1,j,k,iBlock)-&
          FaceAreaI_DFB(:,i  ,j,k,iBlock)+&
          FaceAreaJ_DFB(:,i,j+1,k,iBlock)-&
          FaceAreaJ_DFB(:,i  ,j,k,iBlock)+&
          FaceAreaK_DFB(:,i,j,k+1,iBlock)-&
          FaceAreaK_DFB(:,i  ,j,k,iBlock)
     if(sum(FaceArea_D**2)>cTolerance)then
        write(*,*)'Wrongly defined face areas'
        write(*,*)'i,j,k,iBlock=',i,j,k,iBlock
        write(*,*)'Wrong geometry is created in '//NameSub
        write(*,*)'Face Area Vectors:',&
             FaceAreaI_DFB(:,i+1,j,k,iBlock),&
             FaceAreaI_DFB(:,i  ,j,k,iBlock),&
             FaceAreaJ_DFB(:,i,j+1,k,iBlock),&
             FaceAreaJ_DFB(:,i  ,j,k,iBlock),&
             FaceAreaK_DFB(:,i,j,k+1,iBlock),&
             FaceAreaK_DFB(:,i  ,j,k,iBlock)
        call stop_mpi('Stopped')
     end if
  end do;end do;end do
end subroutine test_block_geometry
!========================================================================
! Formulae for face areas and control volumes are legitimate as the 
! finite ratios Area/d\phi and Volume/d\phi, for the the azimuthal 
! width of the control volume, angle, d\phi, tending to zero. 

subroutine fix_rz_geometry(iBlock)

  use ModCovariant
  use ModNodes,ONLY:NodeY_NB
  use ModGeometry,ONLY: vInv_CB, dz_BLK, y_BLK, dy_BLK, dx_BLK
  use ModSize
  implicit none
  integer,intent(in)::iBlock
  
  integer :: i, j, k
  
  character(LEN=*),parameter::NameSub='fix_rz_geometry'
  !---------------------------------------------------!
  
  !No fluxes along 'k' direction
  FaceAreaK_DFB(:,:,:,:,iBlock) = 0; FaceArea2MinK_B(iBlock)= 1.0
  
  !To transform the face areas and volumes obtained for cartesian geometry
  !we need (1) to eliminate the dependence on meaningless dZ and (2)
  !to add the radius-dependent multipliers (radial coordinate is abs(y)
  do j=1,nJ
     vInv_CB(:,j,:,iBlock) = vInv_CB(:,j,:,iBlock)/abs(y_BLK(1,j,1,iBlock))

     FaceAreaI_DFB(:,:,j,:,iBlock) = FaceAreaI_DFB(:,:,j,:,iBlock) &
           *abs(y_BLK(1,j,1,iBlock))
  end do

  do j=1,nJ+1
     FaceAreaJ_DFB(:,:,j,:,iBlock) = FaceAreaJ_DFB(:,:,j,:,iBlock) &
           *abs(NodeY_NB(1,j,1,iBlock))
  end do

  !At the pole the face area along j direction is zero. Non-zero FaceArea2MinJ 
  !should be set:
  FaceArea2MinJ_B(iBlock) = 1e-30

end subroutine fix_rz_geometry

!========================================================================
subroutine calc_b0source_covar(iBlock)  
  use ModProcMH  
  use ModMain,ONLY:UseB0Source,UseCurlB0,x_,y_,z_,UseB0!,R_,Theta_,Phi_
  use ModSize
  use ModParallel, ONLY :&
       neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth
  use ModCovariant
  use ModAdvance, ONLY : &
       B0_DX, B0_DY, B0_DZ, set_b0_face,&
       B0ResChange_DXSB, B0ResChange_DYSB, B0ResChange_DZSB, &
       CurlB0_DCB,DivB0_CB, NormB0_CB  
  use ModGeometry,ONLY: dx_BLK,dy_BLK,dz_BLK,XyzStart_BLK,&!TypeGeometry,&
       vInv_CB
  use ModNumConst
  implicit none

  integer, intent(in) :: iBlock
  integer::i,j,k,iFace,jFace,kFace,iDirB0,iDirFA,iSide
  integer::i2,j2,iDim
  !real::divB0!,x,y,z,R,Phi,Theta
  real,dimension(3)::GenCoord111_D,dGenFine_D
  real,dimension(ndim,0:1,0:1,East_:Top_)::FaceArea_DIIS, B0_DIIS
  real,dimension(nDim,0:1,0:1,0:1)   :: RefB0_DIII
  real,dimension(nDim,-1:2,-1:2,-1:2):: RefXyz_DIII
  real,dimension(nDim,0:2,0:2,0:2)   :: RefXyzNodes_DIII
  real,dimension(nDim,0:1,0:1,0:1)   :: RefFaceArea_DIII
  real,dimension(nDim,nDim,nI,nJ,nK) :: B0SourceMatrix_DDC
  real::CurlB02,CurlB02_DD(3,3),B0Nabla_DD(3,3)
  real:: eigenvalue_max
  external eigenvalue_max
  if(.not.UseB0)call stop_mpi('Illegal call for calc_B0source_covar')
  dGenFine_D(1)=cHalf*dx_BLK(iBlock)
  dGenFine_D(2)=cHalf*dy_BLK(iBlock)
  dGenFine_D(3)=cHalf*dz_BLK(iBlock)
  call set_b0_face(iBlock)
  do k=1,nK
     do j=1,nJ
        do i=1,nI
           FaceArea_DIIS=cZero
           B0_DIIS=cZero

           if(i==nI.and.neiLWest(iBlock)==-1)then
              call correct_b0_face(West_)
              do iDim=1,3
                 B0ResChange_DXSB(iDim,j,k,West_,iBlock)= &
                      sum(B0_DIIS(iDim,:,:,West_))*0.25
              end do
           else
              FaceArea_DIIS(:,0,0,West_)=FaceAreaI_DFB(:,i+1,j,k,iBlock)
              B0_DIIS(:,0,0,West_)= B0_DX(:,i+1,j,k)
           end if

           if(i==1.and.neiLEast(iBlock)==-1)then
              call correct_b0_face(East_)
              do iDim=1,3
                 B0ResChange_DXSB(iDim,j,k,East_,iBlock)= &
                      sum(B0_DIIS(iDim,:,:,East_))*0.25
              end do
           else
              FaceArea_DIIS(:,0,0,East_)=FaceAreaI_DFB(:,i,j,k,iBlock)
              B0_DIIS(:,0,0,East_)= B0_DX(:,i,j,k)
           end if

           If(j==nJ.and.neiLNorth(iBlock)==-1)then
              call correct_b0_face(North_)
              do iDim=1,3
                 B0ResChange_DYSB(iDim,i,k,North_,iBlock)= &
                      sum(B0_DIIS(iDim,:,:,North_))*0.25
              end do
           else
              FaceArea_DIIS(:,0,0,North_)=FaceAreaJ_DFB(:,i,j+1,k,iBlock) 
              B0_DIIS(:,0,0,North_)= B0_DY(:,i,j+1,k)
           end If

           if(j==1.and.neiLsouth(iBlock)==-1)then
              call correct_b0_face(South_)
              do iDim=1,3
                 B0ResChange_DYSB(iDim,i,k,South_,iBlock)= &
                      sum(B0_DIIS(iDim,:,:,South_))*0.25
              end do
           else
              FaceArea_DIIS(:,0,0,South_)=FaceAreaJ_DFB(:,i,j,k,iBlock)  
              B0_DIIS(:,0,0,South_)= B0_DY(:,i,j,k)
           end if

           if(k==nK.and.neiLTop(iBlock)==-1)then
              call correct_b0_face(Top_)
              do iDim=1,3
                 B0ResChange_DZSB(iDim,i,j,Top_,iBlock)= &
                      sum(B0_DIIS(iDim,:,:,Top_))*0.25
              end do
           else
              FaceArea_DIIS(:,0,0,Top_)=FaceAreaK_DFB(:,i,j,k+1,iBlock)   
              B0_DIIS(:,0,0,Top_)= B0_DZ(:,i,j,k+1)
           end if

           if(k==1.and.neiLBot(iBlock)==-1)then
              call correct_b0_face(Bot_)
              do iDim=1,3
                 B0ResChange_DZSB(iDim,i,j,Bot_,iBlock)= &
                      sum(B0_DIIS(iDim,:,:,Bot_))*0.25
              end do
           else
              FaceArea_DIIS(:,0,0,Bot_)=FaceAreaK_DFB(:,i,j,k,iBlock)  
              B0_DIIS(:,0,0,Bot_)= B0_DZ(:,i,j,k)
           end if

           if((.not.UseB0Source).and.(.not.UseCurlB0))CYCLE
           DivB0_CB(i,j,k,iBlock) = cZero
           B0SourceMatrix_DDC(:,:,i,j,k)=cZero

           do iSide=East_,Bot_,2
              FaceArea_DIIS(:,:,:,iSide)=-FaceArea_DIIS(:,:,:,iSide)
           end do

           do iSide=East_,Top_
              do j2=0,1
                 do i2=0,1
                    DivB0_CB(i,j,k,iBlock)= DivB0_CB(i,j,k,iBlock)+&
                         dot_product(B0_DIIS(:,i2,j2,iSide),&
                         FaceArea_DIIS(:,i2,j2,iSide))    
                    do iDirB0=x_,z_
                       do iDirFA=x_,z_
                          B0SourceMatrix_DDC(iDirFA,iDirB0,i,j,k)= & 
                               B0SourceMatrix_DDC(iDirFA,iDirB0,i,j,k)&
                               +FaceArea_DIIS(iDirB0,i2,j2,iSide)*&
                               B0_DIIS(iDirFA,i2,j2,iSide)&  
                               -FaceArea_DIIS(iDirFA,i2,j2,iSide)*&
                               B0_DIIS(iDirB0,i2,j2,iSide)
                       end do
                    end do
                 end do
              end do
           end do

           DivB0_CB(i,j,k,iBlock)=DivB0_CB(i,j,k,iBlock)*vInv_CB(i,j,k,iBlock)        
           
           
           CurlB0_DCB(z_,i,j,k,iBlock) = &
                B0SourceMatrix_DDC(2,1,i,j,k)*&
                vInv_CB(i,j,k,iBlock)
           CurlB0_DCB(y_,i,j,k,iBlock) = -&
                B0SourceMatrix_DDC(3,1,i,j,k)*&
                vInv_CB(i,j,k,iBlock)
           CurlB0_DCB(x_,i,j,k,iBlock) = &
                B0SourceMatrix_DDC(3,2,i,j,k)*&
                vInv_CB(i,j,k,iBlock)
           CurlB02=sum(CurlB0_DCB(:,i,j,k,iBlock)**2)
           if(.not.UseCurlB0)CYCLE
           if(CurlB02==cZero)then
              NormB0_CB(i,j,k,iBlock)=cZero
           else
              CurlB02_DD=CurlB02*cUnit_DD
              do j2=1,3
                 do i2=1,3
                    CurlB02_DD(i2,j2)=CurlB02_DD(i2,j2)-&
                         CurlB0_DCB(i2,i,j,k,iBlock)*&
                         CurlB0_DCB(j2,i,j,k,iBlock)
                 end do
              end do
              B0Nabla_DD=B0SourceMatrix_DDC(:,:,i,j,k)*&
                   vInv_CB(i,j,k,iBlock)&
                   -DivB0_CB(i,j,k,iBlock)*cUnit_DD
              NormB0_CB(i,j,k,iBlock)= sqrt(sqrt(eigenvalue_max(&
                   matmul(transpose(B0Nabla_DD),matmul(&
                   CurlB02_DD,&
                   B0Nabla_DD))               )&
                   )    )
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

        Difference_S(East_) = -(Var_G(i-1,j,k)+Var_G(i,j,k))
        Difference_S(West_) = +(Var_G(i+1,j,k)+Var_G(i,j,k))
        Difference_S(South_)= -(Var_G(i,j-1,k)+Var_G(i,j,k))
        Difference_S(North_)= +(Var_G(i,j+1,k)+Var_G(i,j,k))
        Difference_S(Bot_)  = -(Var_G(i,j,k-1)+Var_G(i,j,k))
        Difference_S(Top_)  = +(Var_G(i,j,k+1)+Var_G(i,j,k))

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
subroutine covariant_curlb(i,j,k,iBLK,CurlB_D,IsTrueBlock)
  use ModSize
  use ModVarIndexes,ONLY: Bx_,Bz_
  use ModCovariant
  use ModNumConst
  use ModGeometry,ONLY:vInv_CB,true_cell
  use ModAdvance,ONLY:State_VGB
  implicit none
  integer,intent(in)::i,j,k,iBLK
  logical,intent(in)::IsTrueBlock
  real,dimension(3),intent(out)::CurlB_D

  real,dimension(3)::B_D
  real,dimension(3,east_:top_)::MagneticField_DS, FaceArea_DS
  real::VInvHalf


  VInvHalf=chalf*vInv_CB(i,j,k,iBLK)

  FaceArea_DS(:,East_ :West_ )=FaceAreaI_DFB(:,i:i+1,j,k,iBLK)
  FaceArea_DS(:,South_:North_)=FaceAreaJ_DFB(:,i,j:j+1,k,iBLK)
  FaceArea_DS(:,Bot_  :Top_  )=FaceAreaK_DFB(:,i,j,k:k+1,iBLK)   

  if(IsTrueBlock)then
     MagneticField_DS(:,East_ )=-State_VGB(Bx_:Bz_,i-1,j,k,iBLK)
     MagneticField_DS(:,West_ )=+State_VGB(Bx_:Bz_,i+1,j,k,iBLK)
     MagneticField_DS(:,South_)=-State_VGB(Bx_:Bz_,i,j-1,k,iBLK)
     MagneticField_DS(:,North_)=+State_VGB(Bx_:Bz_,i,j+1,k,iBLK)
     MagneticField_DS(:,Bot_)  =-State_VGB(Bx_:Bz_,i,j,k-1,iBLK)
     MagneticField_DS(:,Top_)  =+State_VGB(Bx_:Bz_,i,j,k+1,iBLK)
  else
     if(.not.true_cell(i,j,k,iBLK))then
        CurlB_D=cZero
        return
     end if
     B_D=State_VGB(Bx_:Bz_,i,j,k,iBLK)
     !Input from I faces
     if(.not.true_cell(i-1,j,k,iBLK).and.(.not.true_cell(i+1,j,k,iBLK)))then
        CurlB_D=cZero
        return
     end if
     if(true_cell(i-1,j,k,iBLK))then
        MagneticField_DS(:,East_ )=-(+State_VGB(Bx_:Bz_,i-1,j,k,iBLK)+B_D)
     else
        MagneticField_DS(:,East_ )=-(-State_VGB(Bx_:Bz_,i+1,j,k,iBLK)+2*B_D)
     end if
     if(true_cell(i+1,j,k,iBLK))then
        MagneticField_DS(:,West_ )=+(+State_VGB(Bx_:Bz_,i+1,j,k,iBLK)+B_D)
     else
        MagneticField_DS(:,West_ )=+(-State_VGB(Bx_:Bz_,i-1,j,k,iBLK)+2*B_D)
     end if

     !Input from J faces
     if(.not.true_cell(i,j-1,k,iBLK).and.(.not.true_cell(i,j+1,k,iBLK)))then
        CurlB_D=cZero
        return
     end if

     if(true_cell(i,j-1,k,iBLK))then
        MagneticField_DS(:,South_ )=-(+State_VGB(Bx_:Bz_,i,j-1,k,iBLK)+B_D)
     else
        MagneticField_DS(:,South_ )=-(-State_VGB(Bx_:Bz_,i,j+1,k,iBLK)+2*B_D)
     end if
     if(true_cell(i,j+1,k,iBLK))then
        MagneticField_DS(:,North_ )=+(+State_VGB(Bx_:Bz_,i,j+1,k,iBLK)+B_D)
     else
        MagneticField_DS(:,North_ )=+(-State_VGB(Bx_:Bz_,i,j-1,k,iBLK)+2*B_D)
     end if

     !Input from K faces
     if(.not.true_cell(i,j,k-1,iBLK).and.(.not.true_cell(i,j,k+1,iBLK)))then
        CurlB_D=cZero
        return
     end if

     if(true_cell(i,j,k-1,iBLK))then
        MagneticField_DS(:,Bot_ )=-(+State_VGB(Bx_:Bz_,i,j,k-1,iBLK)+B_D)
     else
        MagneticField_DS(:,Bot_ )=-(-State_VGB(Bx_:Bz_,i,j,k+1,iBLK)+2*B_D)
     end if
     if(true_cell(i,j,k+1,iBLK))then
        MagneticField_DS(:,Top_ )=+(+State_VGB(Bx_:Bz_,i,j,k+1,iBLK)+B_D)
     else
        MagneticField_DS(:,Top_ )=+(-State_VGB(Bx_:Bz_,i,j,k-1,iBLK)+2*B_D)
     end if
  end if
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
  use ModGeometry,ONLY:true_BLK
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
           call covariant_curlb(i,j,k,iBLK,CurlB_D,true_BLK(iBLK))
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
           call covariant_curlb(i,j,k,iBLK,CurlB_D,.true.)
           PlotVar_G(i,j,k)=(CurlB_D(x_)*x_BLK(i,j,k,iBLK)+&
                CurlB_D(y_)*y_BLK(i,j,k,iBLK)+&
                CurlB_D(z_)*z_BLK(i,j,k,iBLK))/&
                R_BLK(i,j,k,iBLK)
        end do
     end do
  end do

end subroutine covar_curlbr_plotvar



!=========================================End covariant.f90=================



!============To be moved to src/covariant_facefluxes.f90=====================



!=========End src/covariant_facefluxes.f90========================

real function integrate_BLK_covar(qnum,qa)             

  ! Return the volume integral of qa, ie. sum(qa*cV_BLK) 
  ! for all used blocks and true cells
  ! Do for each processor separately if qnum=1, otherwise add them all

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,nBLK,nBlockMax,unusedBLK
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
  use ModAdvance,   ONLY: State_VGB,  nVar
  use ModMain,      ONLY: nI, nJ, nK,  nBlock, MaxBlock, UnusedBLK
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

