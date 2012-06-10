!^CFG COPYRIGHT UM

subroutine set_potential_force(iBlock)

  ! Calculate and store gravitational and centrifugal forces

  use ModMain
  use ModAdvance, ONLY : fbody_x_BLK,fbody_y_BLK,fbody_z_BLK
  use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, true_cell
  use BATL_lib, ONLY: IsCartesianGrid, CellSize_DB, CellVolume_GB, &
       FaceNormal_DDFB

  implicit none

  integer, intent(in):: iBlock

  integer :: i,j,k

  real:: DxInv, DyInv, DzInv
  real:: Potential_S(6), x, y, z
  real:: FaceArea_DS(3,6), VInv

  character(len=*), parameter:: NameSub = 'set_potential_force'
  !---------------------------------------------------------------------------
  if(.not.UseGravity .and. .not.UseRotatingFrame) RETURN

  !true_cell note: using true_cell to replace an Rbody test does not apply here

  fbody_x_BLK(:,:,:,iBlock) = 0.0
  fbody_y_BLK(:,:,:,iBlock) = 0.0
  fbody_z_BLK(:,:,:,iBlock) = 0.0

  if(IsCartesianGrid)then
     DxInv = 1.0/CellSize_DB(x_,iBlock)                        
     DyInv = 1.0/CellSize_DB(y_,iBlock)
     DzInv = 1.0/CellSize_DB(z_,iBlock)                        
  end if

  do k=1,nK; do j=1,nJ;  do i=1,nI  
     if(.not.true_cell(i,j,k,iBlock))CYCLE

     x = 0.5*(x_BLK(i-1,j,k,iBlock) + x_BLK(i,j,k,iBlock))
     y = 0.5*(y_BLK(i-1,j,k,iBlock) + y_BLK(i,j,k,iBlock))
     z = 0.5*(z_BLK(i-1,j,k,iBlock) + z_BLK(i,j,k,iBlock))
     Potential_S(east_) = GravityForcePotential(x,y,z) &
          + CentrifugalForcePotential(x,y,z)

     x = 0.5*(x_BLK(i+1,j,k,iBlock)+x_BLK(i,j,k,iBlock))
     y = 0.5*(y_BLK(i+1,j,k,iBlock)+y_BLK(i,j,k,iBlock))
     z = 0.5*(z_BLK(i+1,j,k,iBlock)+z_BLK(i,j,k,iBlock))
     Potential_S(west_) = -GravityForcePotential(x,y,z) &
          - CentrifugalForcePotential(x,y,z)

     x = 0.5*(x_BLK(i,j-1,k,iBlock)+x_BLK(i,j,k,iBlock))
     y = 0.5*(y_BLK(i,j-1,k,iBlock)+y_BLK(i,j,k,iBlock))
     z = 0.5*(z_BLK(i,j-1,k,iBlock)+z_BLK(i,j,k,iBlock))
     Potential_S(south_) = GravityForcePotential(x,y,z) &
          + CentrifugalForcePotential(x,y,z)

     x = 0.5*(x_BLK(i,j+1,k,iBlock)+x_BLK(i,j,k,iBlock))
     y = 0.5*(y_BLK(i,j+1,k,iBlock)+y_BLK(i,j,k,iBlock))
     z = 0.5*(z_BLK(i,j+1,k,iBlock)+z_BLK(i,j,k,iBlock))
     Potential_S(north_) = -GravityForcePotential(x,y,z) &
          - CentrifugalForcePotential(x,y,z)

     x = 0.5*(x_BLK(i,j,k-1,iBlock)+x_BLK(i,j,k,iBlock))
     y = 0.5*(y_BLK(i,j,k-1,iBlock)+y_BLK(i,j,k,iBlock))
     z = 0.5*(z_BLK(i,j,k-1,iBlock)+z_BLK(i,j,k,iBlock))
     Potential_S(bot_) = GravityForcePotential(x,y,z) &
          + CentrifugalForcePotential(x,y,z)

     x = 0.5*(x_BLK(i,j,k+1,iBlock)+x_BLK(i,j,k,iBlock))
     y = 0.5*(y_BLK(i,j,k+1,iBlock)+y_BLK(i,j,k,iBlock))
     z = 0.5*(z_BLK(i,j,k+1,iBlock)+z_BLK(i,j,k,iBlock))
     Potential_S(top_) = -GravityForcePotential(x,y,z) &
          - CentrifugalForcePotential(x,y,z)

     if(IsCartesianGrid)then    

        fbody_x_BLK(i,j,k,iBlock) = &
             (Potential_S(east_)  + Potential_S(west_))*DxInv 
        fbody_y_BLK(i,j,k,iBlock) = &
             (Potential_S(south_) + Potential_S(north_))*DyInv   
        fbody_z_BLK(i,j,k,iBlock) = &
             (Potential_S(bot_) + Potential_S(top_))*DzInv

     else                                    
        VInv = 1.0/CellVolume_GB(i,j,k,iBlock)

        FaceArea_DS(:,1:2) = FaceNormal_DDFB(:,1,i:i+1,j,k,iBlock)
        FaceArea_DS(:,3:4) = FaceNormal_DDFB(:,2,i,j:j+1,k,iBlock)
        FaceArea_DS(:,5:6) = FaceNormal_DDFB(:,3,i,j,k:k+1,iBlock)

        fbody_x_BLK(i,j,k,iBlock) = VInv*&
             dot_product(FaceArea_DS(1,:), Potential_S)
        fbody_y_BLK(i,j,k,iBlock) = VInv*&
             dot_product(FaceArea_DS(2,:), Potential_S)
        fbody_z_BLK(i,j,k,iBlock) = VInv*&
             dot_product(FaceArea_DS(3,:), Potential_S)

     end if

  end do; end do ;end do 

contains

  !============================================================================
  real function GravityForcePotential(X0,Y0,Z0)

    use ModMain, ONLY : UseGravity,GravityDir
    use ModMain, ONLY : UseBody2                !^CFG IF SECONDBODY
    use ModPhysics, ONLY : Gbody, cTolerance
    use ModPhysics, ONLY : Gbody2,xBody2,yBody2,zBody2  !^CFG IF SECONDBODY
    implicit none

    real, intent(in) :: X0,Y0,Z0
    real :: R0
    !-----------------------------------------------------------------------
    if(.not.UseGravity)then
       GravityForcePotential=0.0
    else
       select case(GravityDir)
       case(1)
          GravityForcePotential=-Gbody*X0
       case(2)
          GravityForcePotential=-Gbody*Y0
       case(3)
          GravityForcePotential=-Gbody*Z0
       case(0)
          R0 = sqrt(X0*X0 + Y0*Y0 + Z0*Z0+cTolerance**2)
          GravityForcePotential= Gbody/R0
       case default
          write(*,*)'Impossible value for GravityDir=',GravityDir
          call stop_mpi(NameSub)
       end select

       !^CFG IF SECONDBODY BEGIN
       !\
       ! if there is a second body for which gravity is important 
       ! include it's contribution  here in the following if block
       !
       ! Note that the second body IGNORES the gravity direction and
       ! only does gravity radially towards the body.
       !/
       if (UseBody2) then
          r0=sqrt((X0-xBody2)**2 +(Y0-yBody2)**2 + (Z0-zBody2)**2 &
               + cTolerance**2)
          GravityForcePotential = GravityForcePotential + GBody2/r0
       end if
       !^CFG IF SECONDBODY END
    end if

  end function GravityForcePotential
  !============================================================================
  real function  CentrifugalForcePotential(X0,Y0,Z0)

    ! Evaluates the non-dimensional centrifugal force potential at the
    ! specified location (X0,Y0,Z0) for a rotating coordinate frame.

    use ModMain
    use ModPhysics, ONLY : OmegaBody

    implicit none
    real, intent(in) :: X0, Y0, Z0
    !--------------------------------------------------------------------------

    if(UseRotatingFrame)then
       CentrifugalForcePotential = - OMEGAbody**2 * 0.5*(X0**2+Y0**2)
    else
       CentrifugalForcePotential = 0.0
    end if

  end function CentrifugalForcePotential

end subroutine set_potential_force
