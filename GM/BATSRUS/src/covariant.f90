!^CFG COPYRIGHT UM
!^CFG FILE NOT CARTESIAN

!=============To be left in src/covariant.f90=============
subroutine set_covar_cartesian_geometry(iBLK)
  use ModCovariant
  use ModGeometry
  use ModNumConst
  implicit none
  integer, intent(in) :: iBLK
  real::dx,dy,dz
  dx=dx_BLK(iBLK)
  dy=dy_BLK(iBLK)
  dz=dz_BLK(iBLK)
  FaceAreaI_FB(0:nI+2,0:nJ+1,0:nK+1,iBLK)=dy*dz/dx                             
  FaceAreaJ_FB(0:nI+1,0:nJ+2,0:nK+1,iBLK)=dz*dx/dy                       
  FaceAreaK_FB(0:nI+1,0:nJ+1,0:nK+2,iBLK)=dx*dy/dz
  vInv_CB(:,:,:,iBLK)=cOne/(dx*dy*dz)
  FaceArea2MinI_B(iBLK)=cZero
  FaceArea2MinJ_B(iBLK)=cZero
  FaceArea2MinK_B(iBLK)=cZero
end subroutine set_covar_cartesian_geometry
!---------------------------------------------------------------------
subroutine set_covar_spherical_geometry(iBLK)
  use ModMain
  use ModGeometry,ONLY:x_BLK,y_BLK,z_BLK,R_BLK,dx_BLK,dy_BLK,dz_BLK,&
       DoFixExtraBoundary_B,XyzStart_BLK,XyzMin_D,XyzMax_D,vInv_CB
  use ModCovariant
  use ModNumConst
  implicit none

  integer, intent(in) :: iBLK
  integer :: i,j,k
  real, dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn) ::&           
       Phi_G,Theta_G
  real::dR,dPhi,dTheta
  real,external::SpherFaceAreaI,SpherFaceAreaJ,SpherFaceAreaK

  dR=dx_BLK(iBLK)
  dPhi=dy_BLK(iBLK)
  dTheta=dz_BLK(iBLK)

  do k = 1-gcn, nK+gcn
     do j = 1-gcn, nJ+gcn
        do i = 1-gcn, nI+gcn
           R_BLK(i,j,k,iBLK) = (i-1)*dR  + xyzStart_BLK(R_,iBLK)           
           Phi_G(i,j,k) = (j-1)*dPhi  + xyzStart_BLK(Phi_,iBLK)            
           Theta_G(i,j,k) = (k-1)*dTheta  + xyzStart_BLK(Theta_,iBLK)      
           x_BLK(i,j,k,iBLK) = R_BLK(i,j,k,iBLK)*&                         
                sin(Theta_G(i,j,k))*cos(Phi_G(i,j,k))                      
           y_BLK(i,j,k,iBLK) = R_BLK(i,j,k,iBLK)*&                         
                sin(Theta_G(i,j,k))*sin(Phi_G(i,j,k))                      
           z_BLK(i,j,k,iBLK) = R_BLK(i,j,k,iBLK)*cos(Theta_G(i,j,k))
        end do
     end do
  end do

  ! face scalar fAs are introduced in such a manner that the face area     
  ! vector is equal to                                                     
  ! ((x,y,z)_BLK(i,j,k,iBLK)-(x,y,z)_BLK(i-1,j,k,iBLK))*fA...(i,j,k,iBLK) 
  ! for the face between i and i-1 cells

  do k=0,nK+1; do j=0,nJ+1; do i=0,nI+2
     FaceAreaI_FB(i,j,k,iBLK) =SpherFaceAreaI(i,k,XyzStart_BLK(R_,iBLK),&
          XyzStart_BLK(Theta_,iBLK),dR,dPhi,dTheta)
  end do; end do;end do

  do k=0,nK+1; do j=0,nJ+2; do i=0,nI+1
     FaceAreaJ_FB(i,j,k,iBLK) =SpherFaceAreaJ(i,k,XyzStart_BLK(R_,iBLK),&
          XyzStart_BLK(Theta_,iBLK),dR,dPhi,dTheta)
  end do; end do;end do

  do k=0,nK+2; do j=0,nJ+1; do i=0,nI+1
     FaceAreaK_FB(i,j,k,iBLK) =SpherFaceAreaK(i,k,XyzStart_BLK(R_,iBLK),&
          XyzStart_BLK(Theta_,iBLK),dR,dPhi,dTheta)
  end do; end do;end do

  do k=1,nK
     do j=1,nJ
        do i=1,nI
           vInv_CB(i,j,k,iBLK) = cThree/(cTwo*tan(cHalf*dPhi)*&  
                cTwo*tan(cHalf*dTheta)*sin(Theta_G(i,j,k))*&                 
                dR*(cThree*R_BLK(i,j,k,iBLK)**2+cQuarter*dR**2))
        end do
     end do
  end do

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
end subroutine set_covar_spherical_geometry
!------------------------------------------------------------------
!------------------------------------------------------------------
real function SpherFaceAreaI(i,k,RStart,ThetaStart,dR,dPhi,dTheta)
  use ModNumConst
  implicit none
  integer, intent(in)::i,k
  real,intent(in)::RStart,ThetaStart,dR,dPhi,dTheta
  SpherFaceAreaI=&
       (RStart+(i-1-cHalf)*dR)**2 &                   
       *sin(ThetaStart+(k-1)*dTheta)*&                              
       (cTwo*tan(cHalf*dPhi))*(cTwo*tan(cHalf*dTheta))/dR                 
end function SpherFaceAreaI
!-----------------------------------------------------------------
real function SpherFaceAreaJ(i,k,RStart,ThetaStart,dR,dPhi,dTheta)
  use ModNumConst
  implicit none
  integer, intent(in)::i,k
  real,intent(in)::RStart,ThetaStart,dR,dPhi,dTheta
  SpherFaceAreaJ=&                                
       dR*(cTwo*tan(cHalf*dTheta))/(sin(dPhi)*&                           
       sin(ThetaStart+(k-1)*dTheta))
end function SpherFaceAreaJ
!-----------------------------------------------------------------
real function SpherFaceAreaK(i,k,RStart,ThetaStart,dR,dPhi,dTheta)
  use ModNumConst
  implicit none
  integer, intent(in)::i,k
  real,intent(in)::RStart,ThetaStart,dR,dPhi,dTheta
  SpherFaceAreaK= &                              
       dR*sin(ThetaStart+(k-1-cHalf)*dTheta)*&               
       (cTwo*tan(cHalf*dPhi))/(cos(cHalf*dTheta)*sin(dTheta))
end function SpherFaceAreaK

!-----------------------------------------------------------------
!---------------------------------------------------------------------
subroutine set_covar_spherical2_geometry(iBLK)
  use ModMain
  use ModGeometry,ONLY:x_BLK,y_BLK,z_BLK,R_BLK,dx_BLK,dy_BLK,dz_BLK,&
       DoFixExtraBoundary_B,XyzStart_BLK,XyzMin_D,XyzMax_D,vInv_CB
  use ModCovariant
  use ModNumConst
  implicit none

  integer, intent(in) :: iBLK
  integer :: i,j,k
  real, dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn) ::&           
       Phi_G,Theta_G
  real::dR2,dPhi,dTheta
  real,external::Spher2FaceAreaI,Spher2FaceAreaJ,Spher2FaceAreaK

  dR2=dx_BLK(iBLK)
  dPhi=dy_BLK(iBLK)
  dTheta=dz_BLK(iBLK)

  do k = 1-gcn, nK+gcn
     do j = 1-gcn, nJ+gcn
        do i = 1-gcn, nI+gcn
           R_BLK(i,j,k,iBLK) = exp((i-1)*dR2  + xyzStart_BLK(R_,iBLK))   
           Phi_G(i,j,k) = (j-1)*dPhi  + xyzStart_BLK(Phi_,iBLK)            
           Theta_G(i,j,k) = (k-1)*dTheta  + xyzStart_BLK(Theta_,iBLK)     
           x_BLK(i,j,k,iBLK) = R_BLK(i,j,k,iBLK)*&                         
                sin(Theta_G(i,j,k))*cos(Phi_G(i,j,k))                      
           y_BLK(i,j,k,iBLK) = R_BLK(i,j,k,iBLK)*&                         
                sin(Theta_G(i,j,k))*sin(Phi_G(i,j,k))                      
           z_BLK(i,j,k,iBLK) = R_BLK(i,j,k,iBLK)*cos(Theta_G(i,j,k))
        end do
     end do
  end do

  ! face scalar fAs are introduced in such a manner that the face area     
  ! vector is equal to                                                     
  ! ((x,y,z)_BLK(i,j,k,iBLK)-(x,y,z)_BLK(i-1,j,k,iBLK))*fA...(i,j,k,iBLK) 
  ! for the face between i and i-1 cells

  do k=0,nK+1; do j=0,nJ+1; do i=0,nI+2
     FaceAreaI_FB(i,j,k,iBLK) =Spher2FaceAreaI(i,k,XyzStart_BLK(R_,iBLK),&
          XyzStart_BLK(Theta_,iBLK),dR2,dPhi,dTheta)
  end do; end do;end do

  do k=0,nK+1; do j=0,nJ+2; do i=0,nI+1
     FaceAreaJ_FB(i,j,k,iBLK) =Spher2FaceAreaJ(i,k,XyzStart_BLK(R_,iBLK),&
          XyzStart_BLK(Theta_,iBLK),dR2,dPhi,dTheta)
  end do; end do;end do

  do k=0,nK+2; do j=0,nJ+1; do i=0,nI+1
     FaceAreaK_FB(i,j,k,iBLK) =Spher2FaceAreaK(i,k,XyzStart_BLK(R_,iBLK),&
          XyzStart_BLK(Theta_,iBLK),dR2,dPhi,dTheta)
  end do; end do;end do


  do k=1,nK
     do j=1,nJ
        do i=1,nI
           vInv_CB(i,j,k,iBLK) = cThree/(cTwo*tan(cHalf*dPhi)*&
                cTwo*tan(cHalf*dTheta)*sin(Theta_G(i,j,k))*&                 
                (R_BLK(i,j,k,iBLK)**3)*(&
                0.75*(sinh(dr2)+sinh(2.0*dr2))+cQuarter*sinh(3.0*dr2)))
        end do
     end do
  end do

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

end subroutine set_covar_spherical2_geometry

real function Spher2FaceAreaI(i,k,RStart,ThetaStart,dR2,dPhi,dTheta)
  use ModNumConst
  implicit none
  integer, intent(in)::i,k
  real,intent(in)::RStart,ThetaStart,dR2,dPhi,dTheta
  Spher2FaceAreaI=&
       cQuarter*exp(RStart+(i-1)*dR2)*((exp(-dR2)+cOne)**2) &
       *sin(ThetaStart+(k-1)*dTheta)*&                              
       (cTwo*tan(cHalf*dPhi))*(cTwo*tan(cHalf*dTheta))/(cOne-exp(-dR2))        
end function Spher2FaceAreaI
!-----------------------------------------------------------------
real function Spher2FaceAreaJ(i,k,RStart,ThetaStart,dR2,dPhi,dTheta)
  use ModNumConst
  implicit none
  integer, intent(in)::i,k
  real,intent(in)::RStart,ThetaStart,dR2,dPhi,dTheta
  Spher2FaceAreaJ=&                                
      sinh(dR2)*exp(RStart+(i-1)*dR2)*&
      (cOne+cosh(dr2))*cHalf*&
      (cTwo*tan(cHalf*dTheta))/(sin(dPhi)*&                           
       sin(ThetaStart+(k-1)*dTheta))
end function Spher2FaceAreaJ
!-----------------------------------------------------------------
real function Spher2FaceAreaK(i,k,RStart,ThetaStart,dR2,dPhi,dTheta)
  use ModNumConst
  implicit none
  integer, intent(in)::i,k
  real,intent(in)::RStart,ThetaStart,dR2,dPhi,dTheta
  Spher2FaceAreaK= &                              
       sinh(dR2)*exp(RStart+(i-1)*dR2)*&
      (cOne+cosh(dr2))*cHalf*&      
       sin(ThetaStart+(k-1-cHalf)*dTheta)*&               
       (cTwo*tan(cHalf*dPhi))/(cos(cHalf*dTheta)*sin(dTheta))
end function Spher2FaceAreaK

!-----------------------------------------------------------------
subroutine set_covar_cylindrical_geometry(iBLK)
  use ModMain
  use ModGeometry
  implicit none
  integer,intent(in)::iBLK
  DoFixExtraBoundary_B(iBLK)= XyzStart_BLK(R_,iBLK)-dx_BLK(iBLK)<XyzMin_D(R_)
end subroutine set_covar_cylindrical_geometry
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
  real,dimension(3)::XyzStartFine_D,dXyzFine_D
  real,dimension(ndim,0:1,0:1,East_:Top_)::FaceArea_DIIS, B0_DIIS
  real,dimension(nDim,0:1,0:1,0:1):: RefB0_DIII,RefXyz_DIII

  dXyzFine_D(1)=cHalf*dx_BLK(iBlock)
  dXyzFine_D(2)=cHalf*dy_BLK(iBlock)
  dXyzFine_D(3)=cHalf*dz_BLK(iBlock)

  XyzStartFine_D(:)=XyzStart_BLK(:,iBlock)-cHalf*dXyzFine_D(:)

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
              call calc_faceareaI(i+1,j,k,iBlock,FaceArea_DIIS(:,0,0,West_))
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
              call calc_faceareaI(i,j,k,iBlock,FaceArea_DIIS(:,0,0,East_))
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
              call calc_faceareaJ(i,j+1,k,iBlock,FaceArea_DIIS(:,0,0,North_)) 
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
              call calc_faceareaJ(i,j,k,iBlock,FaceArea_DIIS(:,0,0,South_))   
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
              call calc_faceareaK(i,j,k+1,iBlock,FaceArea_DIIS(:,0,0,Top_))   
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
              call calc_faceareaK(i,j,k,iBlock,FaceArea_DIIS(:,0,0,Bot_))  
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
    real,external::SpherFaceAreaI,SpherFaceAreaJ,SpherFaceAreaK
    real,external::Spher2FaceAreaI,Spher2FaceAreaJ,Spher2FaceAreaK
    select case(TypeGeometry)
    case('spherical')
       select case(iSide)
       case(East_,West_)
          iFace=1+nI*(iSide-East_)
          call get_spher_refined_b0(2*iFace-3,2*j-2,2*k-2)
          do j2=0,1
             do i2=0,1
                B0_DIIS(:,i2,j2,iSide)=(RefB0_DIII(:,0,i2,j2)+&
                     RefB0_DIII(:,1,i2,j2))*cHalf
                FaceArea_DIIS(:,i2,j2,iSide)=(RefXyz_DIII(:,1,i2,j2)-&
                     RefXyz_DIII(:,0,i2,j2))*&
                     SpherFaceAreaI(1+2*nI*(iSide-East_),2*k-1+j2,&
                     XyzStartFine_D(R_),XyzStartFine_D(Theta_),&
                     dXyzFine_D(R_),dXyzFine_D(Phi_),dXyzFine_D(Theta_))
             end do
          end do
       case(South_,North_)
          jFace=1+nJ*(iSide-South_)
          call get_spher_refined_b0(2*i-2,2*jFace-3,2*k-2)
          do j2=0,1
             do i2=0,1
                B0_DIIS(:,i2,j2,iSide)=(RefB0_DIII(:,i2,0,j2)+&
                     RefB0_DIII(:,i2,1,j2))*cHalf
                FaceArea_DIIS(:,i2,j2,iSide)=(RefXyz_DIII(:,i2,1,j2)-&
                     RefXyz_DIII(:,i2,0,j2))*&
                     SpherFaceAreaJ(2*i-1+i2,2*k-1+j2,&
                     XyzStartFine_D(R_),XyzStartFine_D(Theta_),&
                     dXyzFine_D(R_),dXyzFine_D(Phi_),dXyzFine_D(Theta_))
             end do
          end do
       case(Bot_,Top_)
          kFace=1+nK*(iSide-Bot_)
          call get_spher_refined_b0(2*i-2,2*j-2,2*kFace-3)
          do j2=0,1
             do i2=0,1
                B0_DIIS(:,i2,j2,iSide)=(RefB0_DIII(:,i2,j2,0)+&
                     RefB0_DIII(:,i2,j2,1))*cHalf
                FaceArea_DIIS(:,i2,j2,iSide)=(RefXyz_DIII(:,i2,j2,1)-&
                     RefXyz_DIII(:,i2,j2,0))*&
                     SpherFaceAreaK(2*i-1+i2, 2*nK*(iSide-Bot_)+1, &
                     XyzStartFine_D(R_), XyzStartFine_D(Theta_),   &
                     dXyzFine_D(R_),dXyzFine_D(Phi_),dXyzFine_D(Theta_))
             end do
          end do
       end select
    case('spherical_lnr')
       select case(iSide)
       case(East_,West_)
          iFace=1+nI*(iSide-East_)
          call get_spher2_refined_b0(2*iFace-3,2*j-2,2*k-2)
          do j2=0,1
             do i2=0,1
                B0_DIIS(:,i2,j2,iSide)=(RefB0_DIII(:,0,i2,j2)+&
                     RefB0_DIII(:,1,i2,j2))*cHalf
                FaceArea_DIIS(:,i2,j2,iSide)=(RefXyz_DIII(:,1,i2,j2)-&
                     RefXyz_DIII(:,0,i2,j2))*&
                     Spher2FaceAreaI(1+2*nI*(iSide-East_),2*k-1+j2,&
                     XyzStartFine_D(R_),XyzStartFine_D(Theta_),&
                     dXyzFine_D(R_),dXyzFine_D(Phi_),dXyzFine_D(Theta_))
             end do
          end do
       case(South_,North_)
          jFace=1+nJ*(iSide-South_)
          call get_spher2_refined_b0(2*i-2,2*jFace-3,2*k-2)
          do j2=0,1
             do i2=0,1
                B0_DIIS(:,i2,j2,iSide)=(RefB0_DIII(:,i2,0,j2)+&
                     RefB0_DIII(:,i2,1,j2))*cHalf
                FaceArea_DIIS(:,i2,j2,iSide)=(RefXyz_DIII(:,i2,1,j2)-&
                     RefXyz_DIII(:,i2,0,j2))*&
                     Spher2FaceAreaJ(2*i-1+i2,2*k-1+j2,&
                     XyzStartFine_D(R_),XyzStartFine_D(Theta_),&
                     dXyzFine_D(R_),dXyzFine_D(Phi_),dXyzFine_D(Theta_))
             end do
          end do
       case(Bot_,Top_)
          kFace=1+nK*(iSide-Bot_)
          call get_spher2_refined_b0(2*i-2,2*j-2,2*kFace-3)
          do j2=0,1
             do i2=0,1
                B0_DIIS(:,i2,j2,iSide)=(RefB0_DIII(:,i2,j2,0)+&
                     RefB0_DIII(:,i2,j2,1))*cHalf
                FaceArea_DIIS(:,i2,j2,iSide)=(RefXyz_DIII(:,i2,j2,1)-&
                     RefXyz_DIII(:,i2,j2,0))*&
                     Spher2FaceAreaK(2*i-1+i2, 2*nK*(iSide-Bot_)+1, &
                     XyzStartFine_D(R_), XyzStartFine_D(Theta_),   &
                     dXyzFine_D(R_),dXyzFine_D(Phi_),dXyzFine_D(Theta_))
             end do
          end do
       end select
    case('cartesian')
       select case(iSide)
       case(East_,West_)
          iFace=1+nI*(iSide-East_)
          call get_refined_b0(2*iFace-3,2*j-2,2*k-2)
          do j2=0,1
             do i2=0,1
                B0_DIIS(:,i2,j2,iSide)=(RefB0_DIII(:,0,i2,j2)+&
                     RefB0_DIII(:,1,i2,j2))*cHalf
                FaceArea_DIIS(:,i2,j2,iSide)=(RefXyz_DIII(:,1,i2,j2)-&
                     RefXyz_DIII(:,0,i2,j2))*&
                     dXyzFine_D(y_)* dXyzFine_D(z_)/dXyzFine_D(x_)
             end do
          end do
       case(South_,North_)
          jFace=1+nJ*(iSide-South_)
          call get_refined_b0(2*i-2,2*jFace-3,2*k-2)
          do j2=0,1
             do i2=0,1
                B0_DIIS(:,i2,j2,iSide)=(RefB0_DIII(:,i2,0,j2)+&
                     RefB0_DIII(:,i2,1,j2))*cHalf
                FaceArea_DIIS(:,i2,j2,iSide)=(RefXyz_DIII(:,i2,1,j2)-&
                     RefXyz_DIII(:,i2,0,j2))*&
                     dXyzFine_D(z_)* dXyzFine_D(x_)/dXyzFine_D(y_)
             end do
          end do
       case(Bot_,Top_)
          kFace=1+nK*(iSide-Bot_)
          call get_refined_b0(2*i-2,2*j-2,2*kFace-3)
          do j2=0,1
             do i2=0,1
                B0_DIIS(:,i2,j2,iSide)=(RefB0_DIII(:,i2,j2,0)+&
                     RefB0_DIII(:,i2,j2,1))*cHalf
                FaceArea_DIIS(:,i2,j2,iSide)=(RefXyz_DIII(:,i2,j2,1)-&
                     RefXyz_DIII(:,i2,j2,0))*&
                     dXyzFine_D(x_)* dXyzFine_D(y_)/dXyzFine_D(z_)
             end do
          end do
       end select
    case default
       call stop_mpi('Unknown TypeGeometryB0 = '//TypeGeometry)
    end select
  end subroutine correct_b0_face
  subroutine get_refined_b0(iRef,jRef,kRef)
    implicit none
    integer,intent(in)::iRef,jRef,kRef
    integer::k2         
    do k2=0,1; do j2=0,1; do i2=0,1
       x=XyzStart_BLK(x_,iBlock)+(iRef-cHalf+i2)*cHalf*dx_BLK(iBlock)
       y=XyzStart_BLK(y_,iBlock)+(jRef-cHalf+j2)*cHalf*dy_BLK(iBlock)
       z=XyzStart_BLK(z_,iBlock)+(kRef-cHalf+k2)*cHalf*dz_BLK(iBlock)
       call  get_B0(x,y,z,RefB0_DIII(:,i2,j2,k2))
       RefXyz_DIII(x_,i2,j2,k2)=x
       RefXyz_DIII(y_,i2,j2,k2)=y
       RefXyz_DIII(z_,i2,j2,k2)=z
    end do; end do; end do	
  end subroutine get_refined_b0
  subroutine get_spher_refined_b0(iRef,jRef,kRef)
    implicit none
    integer,intent(in)::iRef,jRef,kRef
    integer::k2         
    do k2=0,1; do j2=0,1; do i2=0,1
       R=XyzStart_BLK(x_,iBlock)+(iRef-cHalf+i2)*cHalf*dx_BLK(iBlock)
       Phi=XyzStart_BLK(y_,iBlock)+(jRef-cHalf+j2)*cHalf*dy_BLK(iBlock)
       Theta=XyzStart_BLK(z_,iBlock)+(kRef-cHalf+k2)*cHalf*dz_BLK(iBlock)
       z=R*cos(Theta)
       r=R*sin(Theta)
       x=r*cos(Phi)
       y=r*sin(Phi)
       call  get_B0(x,y,z,RefB0_DIII(:,i2,j2,k2))
       RefXyz_DIII(x_,i2,j2,k2)=x
       RefXyz_DIII(y_,i2,j2,k2)=y
       RefXyz_DIII(z_,i2,j2,k2)=z
    end do; end do; end do	
  end subroutine get_spher_refined_b0
  subroutine get_spher2_refined_b0(iRef,jRef,kRef)
    implicit none
    integer,intent(in)::iRef,jRef,kRef
    integer::k2         
    do k2=0,1; do j2=0,1; do i2=0,1
       R=exp(XyzStart_BLK(x_,iBlock)+(iRef-cHalf+i2)*cHalf*dx_BLK(iBlock))
       Phi=XyzStart_BLK(y_,iBlock)+(jRef-cHalf+j2)*cHalf*dy_BLK(iBlock)
       Theta=XyzStart_BLK(z_,iBlock)+(kRef-cHalf+k2)*cHalf*dz_BLK(iBlock)
       z=R*cos(Theta)
       r=R*sin(Theta)
       x=r*cos(Phi)
       y=r*sin(Phi)
       call  get_B0(x,y,z,RefB0_DIII(:,i2,j2,k2))
       RefXyz_DIII(x_,i2,j2,k2)=x
       RefXyz_DIII(y_,i2,j2,k2)=y
       RefXyz_DIII(z_,i2,j2,k2)=z
    end do; end do; end do	
  end subroutine get_spher2_refined_b0
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

  call calc_faceareaI(i,j,k,iBLK,FaceArea_DS(:,East_))
  call calc_faceareaI(i+1,j,k,iBLK,FaceArea_DS(:,West_))
  call calc_faceareaJ(i,j,k,iBLK,FaceArea_DS(:,South_))
  call calc_faceareaJ(i,j+1,k,iBLK,FaceArea_DS(:,North_))
  call calc_faceareaK(i,j,k,iBLK,FaceArea_DS(:,Bot_))
  call calc_faceareaK(i,j,k+1,iBLK,FaceArea_DS(:,Top_))

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

        call calc_faceareaI(i,j,k,iBlock,FaceArea_DS(:,East_))
        call calc_faceareaI(i+1,j,k,iBlock,FaceArea_DS(:,West_))
        call calc_faceareaJ(i,j,k,iBlock,FaceArea_DS(:,South_))
        call calc_faceareaJ(i,j+1,k,iBlock,FaceArea_DS(:,North_))
        call calc_faceareaK(i,j,k,iBlock,FaceArea_DS(:,Bot_))
        call calc_faceareaK(i,j,k+1,iBlock,FaceArea_DS(:,Top_))

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
        
        call calc_faceareaI(i,j,k,iBlock,FaceArea_DS(:,East_))
        call calc_faceareaI(i+1,j,k,iBlock,FaceArea_DS(:,West_))
        call calc_faceareaJ(i,j,k,iBlock,FaceArea_DS(:,South_))
        call calc_faceareaJ(i,j+1,k,iBlock,FaceArea_DS(:,North_))
        call calc_faceareaK(i,j,k,iBlock,FaceArea_DS(:,Bot_))
        call calc_faceareaK(i,j,k+1,iBlock,FaceArea_DS(:,Top_))

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

  call calc_faceareaI(i,j,k,iBLK,FaceArea_DS(:,East_))
  call calc_faceareaI(i+1,j,k,iBLK,FaceArea_DS(:,West_))
  call calc_faceareaJ(i,j,k,iBLK,FaceArea_DS(:,South_))
  call calc_faceareaJ(i,j+1,k,iBLK,FaceArea_DS(:,North_))
  call calc_faceareaK(i,j,k,iBLK,FaceArea_DS(:,Bot_))
  call calc_faceareaK(i,j,k+1,iBLK,FaceArea_DS(:,Top_))

  MagneticField_DS(1,East_)=-State_VGB(Bx_,i-1,j,k,iBLK)
  MagneticField_DS(2,East_)=-State_VGB(By_,i-1,j,k,iBLK)
  MagneticField_DS(3,East_)=-State_VGB(Bz_,i-1,j,k,iBLK)

  MagneticField_DS(1,West_)=State_VGB(Bx_,i+1,j,k,iBLK)
  MagneticField_DS(2,West_)=State_VGB(By_,i+1,j,k,iBLK)
  MagneticField_DS(3,West_)=State_VGB(Bz_,i+1,j,k,iBLK)

  MagneticField_DS(1,South_)=-State_VGB(Bx_,i,j-1,k,iBLK)
  MagneticField_DS(2,South_)=-State_VGB(By_,i,j-1,k,iBLK)
  MagneticField_DS(3,South_)=-State_VGB(Bz_,i,j-1,k,iBLK)

  MagneticField_DS(1,North_)=State_VGB(Bx_,i,j+1,k,iBLK)
  MagneticField_DS(2,North_)=State_VGB(By_,i,j+1,k,iBLK)
  MagneticField_DS(3,North_)=State_VGB(Bz_,i,j+1,k,iBLK)

  MagneticField_DS(1,Bot_)=-State_VGB(Bx_,i,j,k-1,iBLK)
  MagneticField_DS(2,Bot_)=-State_VGB(By_,i,j,k-1,iBLK)
  MagneticField_DS(3,Bot_)=-State_VGB(Bz_,i,j,k-1,iBLK)

  MagneticField_DS(1,Top_)=State_VGB(Bx_,i,j,k+1,iBLK)
  MagneticField_DS(2,Top_)=State_VGB(By_,i,j,k+1,iBLK)
  MagneticField_DS(3,Top_)=State_VGB(Bz_,i,j,k+1,iBLK)

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
  real,dimension(nDim) :: B_D,FaceArea_D


  do k=1,nK; do j=1,nJ
     call calc_faceareaI(iFaceIn,j,k,iBlock,FaceArea_D)
     B_D(x_)=LeftState_VX(Bx_,iFaceIn,j,k)
     B_D(y_)=LeftState_VX(By_,iFaceIn,j,k)
     B_D(z_)=LeftState_VX(Bz_,iFaceIn,j,k)
     CorrectedFlux_VXB(BnL_,j,k,iFaceOut,iBlock) = &
          dot_product(B_D,FaceArea_D)
     B_D(x_)=RightState_VX(Bx_,iFaceIn,j,k)
     B_D(y_)=RightState_VX(By_,iFaceIn,j,k)
     B_D(z_)=RightState_VX(Bz_,iFaceIn,j,k)
     CorrectedFlux_VXB(BnR_,j,k,iFaceOut,iBlock) = &
          dot_product(B_D,FaceArea_D)
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
  real,dimension(nDim) :: B_D,FaceArea_D


  do k=1,nK; do i=1,nI
     call calc_faceareaJ(i,jFaceIn,k,iBlock,FaceArea_D)
     B_D(x_)=LeftState_VY(Bx_,i,jFaceIn,k)
     B_D(y_)=LeftState_VY(By_,i,jFaceIn,k)
     B_D(z_)=LeftState_VY(Bz_,i,jFaceIn,k)
     CorrectedFlux_VYB(BnL_,i,k,jFaceOut,iBlock) = &
          dot_product(B_D,FaceArea_D)
     B_D(x_)=RightState_VY(Bx_,i,jFaceIn,k)
     B_D(y_)=RightState_VY(By_,i,jFaceIn,k)
     B_D(z_)=RightState_VY(Bz_,i,jFaceIn,k)
     CorrectedFlux_VYB(BnR_,i,k,jFaceOut,iBlock) = dot_product(B_D,FaceArea_D)
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
  real,dimension(nDim) :: B_D, FaceArea_D


  do j=1,nJ; do i=1,nI
     call calc_faceareaK(i,j,kFaceIn,iBlock,FaceArea_D)
     B_D(x_)=LeftState_VZ(Bx_,i,j,kFaceIn)
     B_D(y_)=LeftState_VZ(By_,i,j,kFaceIn)
     B_D(z_)=LeftState_VZ(Bz_,i,j,kFaceIn)
     CorrectedFlux_VZB(BnL_,i,j,kFaceOut,iBlock) = &
          dot_product(B_D,FaceArea_D)
     B_D(x_)=RightState_VZ(Bx_,i,j,kFaceIn)
     B_D(y_)=RightState_VZ(By_,i,j,kFaceIn)
     B_D(z_)=RightState_VZ(Bz_,i,j,kFaceIn)
     CorrectedFlux_VZB(BnR_,i,j,kFaceOut,iBlock) =&
           dot_product(B_D,FaceArea_D)
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
     call calc_faceareaI(iFaceOut,j,k,iBlock,FaceArea_D)
     FaceArea2=dot_product(FaceArea_D,FaceArea_D)

     B_D(x_)=LeftState_VX(Bx_,iFaceOut,j,k)
     B_D(y_)=LeftState_VX(By_,iFaceOut,j,k)
     B_D(z_)=LeftState_VX(Bz_,iFaceOut,j,k)

     DeltaBDotFA = (CorrectedFlux_VXB(BnL_,j,k,iFaceIn,iBlock) -&
          dot_product(B_D,FaceArea_D))/FaceArea2

     LeftState_VX(Bx_,iFaceOut,j,k)=B_D(x_)+DeltaBDotFA*FaceArea_D(x_)
     LeftState_VX(By_,iFaceOut,j,k)=B_D(y_)+DeltaBDotFA*FaceArea_D(y_)
     LeftState_VX(Bz_,iFaceOut,j,k)=B_D(z_)+DeltaBDotFA*FaceArea_D(z_)

     B_D(x_)=RightState_VX(Bx_,iFaceOut,j,k)
     B_D(y_)=RightState_VX(By_,iFaceOut,j,k)
     B_D(z_)=RightState_VX(Bz_,iFaceOut,j,k)

     DeltaBDotFA = (CorrectedFlux_VXB(BnR_,j,k,iFaceIn,iBlock) -&
          dot_product(B_D,FaceArea_D))/FaceArea2

     RightState_VX(Bx_,iFaceOut,j,k)=B_D(x_)+DeltaBDotFA*FaceArea_D(x_)
     RightState_VX(By_,iFaceOut,j,k)=B_D(y_)+DeltaBDotFA*FaceArea_D(y_)
     RightState_VX(Bz_,iFaceOut,j,k)=B_D(z_)+DeltaBDotFA*FaceArea_D(z_)
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
     call calc_faceareaJ(i,jFaceOut,k,iBlock,FaceArea_D)
     FaceArea2=dot_product(FaceArea_D,FaceArea_D)

     B_D(x_)=LeftState_VY(Bx_,i,jFaceOut,k)
     B_D(y_)=LeftState_VY(By_,i,jFaceOut,k)
     B_D(z_)=LeftState_VY(Bz_,i,jFaceOut,k)

     DeltaBDotFA = (CorrectedFlux_VYB(BnL_,i,k,jFaceIn,iBlock)-&
          dot_product(B_D,FaceArea_D))/FaceArea2

     LeftState_VY(Bx_,i,jFaceOut,k)=B_D(x_)+DeltaBDotFA*FaceArea_D(x_)
     LeftState_VY(By_,i,jFaceOut,k)=B_D(y_)+DeltaBDotFA*FaceArea_D(y_)
     LeftState_VY(Bz_,i,jFaceOut,k)=B_D(z_)+DeltaBDotFA*FaceArea_D(z_)

     B_D(x_)=RightState_VY(Bx_,i,jFaceOut,k)
     B_D(y_)=RightState_VY(By_,i,jFaceOut,k)
     B_D(z_)=RightState_VY(Bz_,i,jFaceOut,k)

     DeltaBDotFA = (CorrectedFlux_VYB(BnR_,i,k,jFaceIn,iBlock)-&
          dot_product(B_D,FaceArea_D))/FaceArea2

     RightState_VY(Bx_,i,jFaceOut,k)=B_D(x_)+DeltaBDotFA*FaceArea_D(x_)
     RightState_VY(By_,i,jFaceOut,k)=B_D(y_)+DeltaBDotFA*FaceArea_D(y_)
     RightState_VY(Bz_,i,jFaceOut,k)=B_D(z_)+DeltaBDotFA*FaceArea_D(z_)

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
     call calc_faceareaK(i,j,kFaceOut,iBlock,FaceArea_D)
     FaceArea2=dot_product(FaceArea_D,FaceArea_D)

     B_D(x_)=LeftState_VZ(Bx_,i,j,kFaceOut)
     B_D(y_)=LeftState_VZ(By_,i,j,kFaceOut)
     B_D(z_)=LeftState_VZ(Bz_,i,j,kFaceOut)

     DeltaBDotFA = ( CorrectedFlux_VZB(BnL_,i,j,kFaceIn,iBlock) -&
          dot_product(B_D,FaceArea_D))/FaceArea2

     LeftState_VZ(Bx_,i,j,kFaceOut) = B_D(x_)+DeltaBDotFA*FaceArea_D(x_) 
     LeftState_VZ(By_,i,j,kFaceOut) = B_D(y_)+DeltaBDotFA*FaceArea_D(y_)
     LeftState_VZ(Bz_,i,j,kFaceOut) = B_D(z_)+DeltaBDotFA*FaceArea_D(z_) 

     B_D(x_)=RightState_VZ(Bx_,i,j,kFaceOut)
     B_D(y_)=RightState_VZ(By_,i,j,kFaceOut)
     B_D(z_)=RightState_VZ(Bz_,i,j,kFaceOut)
     
     DeltaBDotFA = (CorrectedFlux_VZB(BnR_,i,j,kFaceIn,iBlock) -&
          dot_product(B_D,FaceArea_D))/FaceArea2

     RightState_VZ(Bx_,i,j,kFaceOut) = B_D(x_)+DeltaBDotFA*FaceArea_D(x_) 
     RightState_VZ(By_,i,j,kFaceOut) = B_D(y_)+DeltaBDotFA*FaceArea_D(y_)
     RightState_VZ(Bz_,i,j,kFaceOut) = B_D(z_)+DeltaBDotFA*FaceArea_D(z_) 
  end do; end do

end subroutine apply_bn_faceK
!---------------------------------------------------------------------------


!=========================================End covariant.f90=================



!============To be moved to src/covariant_facefluxes.f90=====================



!=========End src/covariant_facefluxes.f90========================

real function integrate_BLK(qnum,qa)             

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
     integrate_BLK=qsum_all
     if(oktest)write(*,*)'me,sum,sum_all:',iProc,qsum,qsum_all
  else
     integrate_BLK=qsum
     if(oktest)write(*,*)'me,qsum:',iProc,qsum
  end if
  
end function integrate_BLK

subroutine integrate_domain(Sum_V, Pressure_GB)
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

end subroutine integrate_domain
!==============End of library programms=================================
