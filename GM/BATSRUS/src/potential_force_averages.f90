!^CFG COPYRIGHT UM

subroutine body_force_averages

  ! Calculate and store gravitational and centrifugal forces

  use ModMain
  use ModAdvance, ONLY : fbody_x_BLK,fbody_y_BLK,fbody_z_BLK
  use ModGeometry, ONLY :&
       dx_BLK, dy_BLK, dz_BLK,&                      
       UseCovariant,          &                       
       x_BLK, y_BLK, z_BLK, true_cell

  implicit none

  integer :: i,j,k

  real ::  DxInv,DyInv,DzInv                          
  real ::  Potential_S(6),x,y,z
  !---------------------------------------------------------------------------

  if(.not.UseGravity .and. .not.UseRotatingFrame) RETURN

  !true_cell note: using true_cell to replace an Rbody test does not apply here

  fbody_x_BLK(:,:,:,globalBLK) = 0.0
  fbody_y_BLK(:,:,:,globalBLK) = 0.0
  fbody_z_BLK(:,:,:,globalBLK) = 0.0

  DxInv = 1.0/dx_BLK(globalBLK)                        
  DyInv = 1.0/dy_BLK(globalBLK)
  DzInv = 1.0/dz_BLK(globalBLK)                        

  do k=1,nK; do j=1,nJ;  do i=1,nI  
     if(.not.true_cell(i,j,k,globalBLK))CYCLE

     x = 0.5*(x_BLK(i-1,j,k,globalBLK) + x_BLK(i,j,k,globalBLK))
     y = 0.5*(y_BLK(i-1,j,k,globalBLK) + y_BLK(i,j,k,globalBLK))
     z = 0.5*(z_BLK(i-1,j,k,globalBLK) + z_BLK(i,j,k,globalBLK))
     Potential_S(east_) = GravityForcePotential(x,y,z) &
          + CentrifugalForcePotential(x,y,z)

     x = 0.5*(x_BLK(i+1,j,k,globalBLK)+x_BLK(i,j,k,globalBLK))
     y = 0.5*(y_BLK(i+1,j,k,globalBLK)+y_BLK(i,j,k,globalBLK))
     z = 0.5*(z_BLK(i+1,j,k,globalBLK)+z_BLK(i,j,k,globalBLK))
     Potential_S(west_) = -GravityForcePotential(x,y,z) &
          - CentrifugalForcePotential(x,y,z)

     x = 0.5*(x_BLK(i,j-1,k,globalBLK)+x_BLK(i,j,k,globalBLK))
     y = 0.5*(y_BLK(i,j-1,k,globalBLK)+y_BLK(i,j,k,globalBLK))
     z = 0.5*(z_BLK(i,j-1,k,globalBLK)+z_BLK(i,j,k,globalBLK))
     Potential_S(south_) = GravityForcePotential(x,y,z) &
          + CentrifugalForcePotential(x,y,z)

     x = 0.5*(x_BLK(i,j+1,k,globalBLK)+x_BLK(i,j,k,globalBLK))
     y = 0.5*(y_BLK(i,j+1,k,globalBLK)+y_BLK(i,j,k,globalBLK))
     z = 0.5*(z_BLK(i,j+1,k,globalBLK)+z_BLK(i,j,k,globalBLK))
     Potential_S(north_) = -GravityForcePotential(x,y,z) &
          - CentrifugalForcePotential(x,y,z)

     x = 0.5*(x_BLK(i,j,k-1,globalBLK)+x_BLK(i,j,k,globalBLK))
     y = 0.5*(y_BLK(i,j,k-1,globalBLK)+y_BLK(i,j,k,globalBLK))
     z = 0.5*(z_BLK(i,j,k-1,globalBLK)+z_BLK(i,j,k,globalBLK))
     Potential_S(bot_) = GravityForcePotential(x,y,z) &
          + CentrifugalForcePotential(x,y,z)

     x = 0.5*(x_BLK(i,j,k+1,globalBLK)+x_BLK(i,j,k,globalBLK))
     y = 0.5*(y_BLK(i,j,k+1,globalBLK)+y_BLK(i,j,k,globalBLK))
     z = 0.5*(z_BLK(i,j,k+1,globalBLK)+z_BLK(i,j,k,globalBLK))
     Potential_S(top_) = -GravityForcePotential(x,y,z) &
          - CentrifugalForcePotential(x,y,z)
     if(.not.UseCovariant)then    

        fbody_x_BLK(i,j,k,globalBLK) = &
             (Potential_S(east_)  + Potential_S(west_))*DxInv 
        fbody_y_BLK(i,j,k,globalBLK) = &
             (Potential_S(south_) + Potential_S(north_))*DyInv   
        fbody_z_BLK(i,j,k,globalBLK) = &
             (Potential_S(bot_) + Potential_S(top_))*DzInv

        continue
     else                                    
        call covariant_force_integral(i,j,k,globalBLK,Potential_S)     
     end if

  end do; end do ;end do 
contains

  !==============================================================================
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
          call stop_mpi('ERROR in body_force_averages !!!')
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
          r0=sqrt((X0-xBody2)**2+(Y0-yBody2)**2 + (Z0-zBody2)**2 + cTolerance**2)
          GravityForcePotential = GravityForcePotential + GBody2/R0
       end if
       !^CFG IF SECONDBODY END
    end if

  end function GravityForcePotential
  !=============================================================================
  real function  CentrifugalForcePotential(X0,Y0,Z0)

    ! Evaluates the non-dimensional centrifugal force potential at the
    ! specified location (X0,Y0,Z0) for a rotating coordinate frame.

    use ModMain
    use ModPhysics, ONLY : OMEGAbody

    implicit none
    real, intent(in) :: X0,Y0,Z0
    !--------------------------------------------------------------------------

    if(UseRotatingFrame)then
       CentrifugalForcePotential = - OMEGAbody**2 * 0.5*(X0**2+Y0**2)
    else
       CentrifugalForcePotential = 0.0
    end if

  end function CentrifugalForcePotential

end subroutine body_force_averages
