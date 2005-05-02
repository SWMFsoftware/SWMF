!^CFG COPYRIGHT UM
!---------------------------------------------------------------------------
!                     body_force_averages
!---------------------------------------------------------------------------
!==========================================================================
!==========================================================================
!==========================================================================
module ModGaussQuad
  !\
  ! This module declares the arrays of Gauss-Legendre Quadrature
  ! weights and nodes.  This is done so that they are declared and
  ! stored in one place for the routines that need them.
  ! Three sets of the weights and nodes
  ! are listed below so that if more accuracy is desired the number of
  ! nodes can easily be changed.  It was found that 6 nodes worked well
  ! for this application.  This was a compromise between speed and
  ! accuracy.  See documentation for more information.
  !
  ! Note that the variable NGaussNodes is the number of unique values
  ! of nodes and weights for the integration.  Since the integration
  ! is done on the interval (-1,1) each of the values is used twice.
  ! This means that NGaussNodes is one half the number of integration
  ! points used.
  !/

!!$  !\
!!$  ! This is the block for 4 node integration
!!$  !
!!$  integer, parameter :: nGaussNodes=2
!!$  real, parameter, dimension(1:2) :: GaussNodes = &
!!$       (/0.339981043584856, 0.861136311594053 /)
!!$  real, parameter, dimension(1:2) :: GaussWeights = &
!!$       (/0.652145154862546, 0.347854845137454/)

  !\
  ! This is the block for 6 node integration
  !
  integer, parameter :: nGaussNodes=3
  real, parameter, dimension(1:3) :: GaussNodes = &
       (/0.23861918608319693, 0.66120938646626459, 0.93246951420315205/)
  real, parameter, dimension(1:3) :: GaussWeights = &
       (/0.46791393457269126, 0.36076157304813861, 0.17132449237916234/)

!!$  !\
!!$  ! This is the block for 8 node integration
!!$  !
!!$  integer, parameter :: nGaussNodes=4
!!$  real, parameter, dimension(1:4) :: GaussNodes = &
!!$       (/0.183434642495650, 0.525532409916329, &
!!$       0.796666477413627, 0.960289856497536 /)
!!$  real, parameter, dimension(1:4) :: GaussWeights = &
!!$       (/0.362683783378362, 0.313706645877887,  &
!!$       0.222381034453374, 0.101228536290376 /)

end module ModGaussQuad
!==========================================================================
!==========================================================================
!==========================================================================
!\
! This routine computes the volume-averaged body force for each of the 
! cells on the grid.
!/

subroutine body_force_averages
  use ModMain, ONLY : nI,nJ,nK,globalBLK
  use ModMain, ONLY : UseBody2                !^CFG IF SECONDBODY
  use ModAdvance, ONLY : fbody_x_BLK,fbody_y_BLK,fbody_z_BLK
  use ModGeometry, ONLY : R_BLK
  use ModGeometry, ONLY : R2_BLK        !^CFG IF SECONDBODY
  use ModPhysics, ONLY : Rbody
  use ModPhysics, ONLY : Rbody2               !^CFG IF SECONDBODY
  implicit none

  integer :: i,j,k
  real :: Rcrit, Rcrit2

  !true_cell note: using true_cell to replace an Rbody test does not apply here

  if (.not. UseBody2) then                    !^CFG IF SECONDBODY
  ! if there is only one body do this branch  !^CFG IF SECONDBODY

     Rcrit = 0.50*Rbody

     fbody_x_BLK(:,:,:,globalBLK) = 0.00
     fbody_y_BLK(:,:,:,globalBLK) = 0.00
     fbody_z_BLK(:,:,:,globalBLK) = 0.00

     do k=1,nK
        do j=1,nJ
           do i=1,nI
              if ( R_BLK(i,j,k,globalBLK) > Rcrit .or. &
                   R_BLK(i-1,j,k,globalBLK) > Rcrit .or. &
                   R_BLK(i+1,j,k,globalBLK) > Rcrit .or. &
                   R_BLK(i,j-1,k,globalBLK) > Rcrit .or. &
                   R_BLK(i,j+1,k,globalBLK) > Rcrit .or. &
                   R_BLK(i,j,k-1,globalBLK) > Rcrit .or. &
                   R_BLK(i,j,k+1,globalBLK) > Rcrit ) then 
                 call body_force_cell_average(i,j,k)
              end if
           end do !end i loop
        end do !end j loop
     end do !end k loop

  !^CFG IF SECONDBODY BEGIN
  else               
  ! else do the case for 2 bodies

     Rcrit  = 0.50*RBody
     Rcrit2 = 0.50*RBody2

     fbody_x_BLK(:,:,:,globalBLK) = 0.00
     fbody_y_BLK(:,:,:,globalBLK) = 0.00
     fbody_z_BLK(:,:,:,globalBLK) = 0.00

     do k=1,nK
        do j=1,nJ
           do i=1,nI
              if ( (R_BLK(i,j,k,globalBLK) > Rcrit .or. &
                   R_BLK(i-1,j,k,globalBLK) > Rcrit .or. &
                   R_BLK(i+1,j,k,globalBLK) > Rcrit .or. &
                   R_BLK(i,j-1,k,globalBLK) > Rcrit .or. &
                   R_BLK(i,j+1,k,globalBLK) > Rcrit .or. &
                   R_BLK(i,j,k-1,globalBLK) > Rcrit .or. &
                   R_BLK(i,j,k+1,globalBLK) > Rcrit)     &
                   .and.                                    &
                   (R2_BLK(i,j,k,globalBLK) > Rcrit2 .or. &
                   R2_BLK(i-1,j,k,globalBLK) > Rcrit2 .or. &
                   R2_BLK(i+1,j,k,globalBLK) > Rcrit2 .or. &
                   R2_BLK(i,j-1,k,globalBLK) > Rcrit2 .or. &
                   R2_BLK(i,j+1,k,globalBLK) > Rcrit2 .or. &
                   R2_BLK(i,j,k-1,globalBLK) > Rcrit2 .or. &
                   R2_BLK(i,j,k+1,globalBLK) > Rcrit2)  ) then 
                 call body_force_cell_average(i,j,k)
              end if
           end do !end i loop
        end do !end j loop
     end do !end k loop

  end if
  !^CFG IF SECONDBODY END

end subroutine body_force_averages

!---------------------------------------------------------------------------
!                     body_force_cell_average
!---------------------------------------------------------------------------

!\
! Calculates the cell averages of the gravitational force field
!/

subroutine body_force_cell_average(i,j,k)
  use ModMain, ONLY : globalBLK,UseGravity
  use ModAdvance, ONLY : fbody_x_BLK,fbody_y_BLK,fbody_z_BLK
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,dx_BLK,dy_BLK,dz_BLK,VolumeInverse_I

  implicit none

  integer, intent(in) :: i,j,k

  integer :: ii
  real ::  xmin,xmax,ymin,ymax,zmin,zmax
  real, dimension(1:3) :: Fave

  !------------------------------------------------------------------------

  if(.not.UseGravity) RETURN

  xmin = x_BLK(i,j,k,globalBLK) - 0.50*dx_BLK(globalBLK)
  ymin = y_BLK(i,j,k,globalBLK) - 0.50*dy_BLK(globalBLK)
  zmin = z_BLK(i,j,k,globalBLK) - 0.50*dz_BLK(globalBLK)

  xmax = x_BLK(i,j,k,globalBLK) + 0.50*dx_BLK(globalBLK)
  ymax = y_BLK(i,j,k,globalBLK) + 0.50*dy_BLK(globalBLK)
  zmax = z_BLK(i,j,k,globalBLK) + 0.50*dz_BLK(globalBLK)

  call body_force_GaussQuadXYZ(xmin,xmax,ymin,ymax,zmin,zmax,Fave)

  do ii=1,3
     Fave(ii) = Fave(ii)*VolumeInverse_I(globalBLK)
  end do

  fbody_x_BLK(i,j,k,globalBLK) = fbody_x_BLK(i,j,k,globalBLK) + Fave(1)
  fbody_y_BLK(i,j,k,globalBLK) = fbody_y_BLK(i,j,k,globalBLK) + Fave(2)
  fbody_z_BLK(i,j,k,globalBLK) = fbody_z_BLK(i,j,k,globalBLK) + Fave(3)

end subroutine body_force_cell_average

!----------------------------------------------------------------------
!                          gravity_force
!----------------------------------------------------------------------

!\
! Evaluates the non-dimensional gravitational force at the
! specified location (X0,Y0,Z0).
!/

subroutine gravity_force(X0,Y0,Z0,g0)
  use ModMain, ONLY : UseGravity,GravityDir
  use ModMain, ONLY : UseBody2                !^CFG IF SECONDBODY
  use ModPhysics, ONLY : Gbody
  use ModPhysics, ONLY : Gbody2               !^CFG IF SECONDBODY
  implicit none

  real, intent(in) :: X0,Y0,Z0
  real, intent(out), dimension(1:3) :: g0

  real :: R0
  !-----------------------------------------------------------------------
  if(.not.UseGravity)then
     g0(1:3)=0.0
     RETURN
  end if

  R0 = sqrt(X0*X0 + Y0*Y0 + Z0*Z0)
  R0 = max(R0,1.00E-6)

  select case(GravityDir)
  case(1,2,3)
     g0(1:3)=0.0
     g0(GravityDir) = Gbody
  case(0)
     g0(1) = Gbody*X0/R0**3
     g0(2) = Gbody*Y0/R0**3
     g0(3) = Gbody*Z0/R0**3
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
     g0(1) = g0(1) + GBody2*X0/R0**3
     g0(2) = g0(2) + GBody2*Y0/R0**3
     g0(3) = g0(3) + GBody2*Z0/R0**3
  end if
  !^CFG IF SECONDBODY END


end subroutine gravity_force

!----------------------------------------------------------------------
!                        centripetal_force
!----------------------------------------------------------------------

!\
! Evaluates the non-dimensional centripetal force at the
! specified location (X0,Y0,Z0) for a rotating coordinate frame.
!/

subroutine centripetal_force(X0,Y0,Z0,f0)
  use ModMain
  use ModPhysics, ONLY : OMEGAbody
  implicit none

  real, intent(in) :: X0,Y0,Z0
  real, intent(out), dimension(1:3) :: f0

  select case (problem_type)             !^CFG IF NOT SIMPLE BEGIN
  case (problem_heliosphere)
     if (.not.UseRotatingFrame) then
        f0(1) = 0.00
        f0(2) = 0.00
        f0(3) = 0.00
     else 
        f0(1) = OMEGAbody*OMEGAbody*X0
        f0(2) = OMEGAbody*OMEGAbody*Y0
        f0(3) = 0.00
     endif
  case default                           !^CFG END SIMPLE
     f0(1) = 0.00
     f0(2) = 0.00
     f0(3) = 0.00
  end select                             !^CFG IF NOT SIMPLE

end subroutine centripetal_force

!----------------------------------------------------------------------
!                       body_force_GaussQuadXYZ
!----------------------------------------------------------------------

subroutine body_force_GaussQuadXYZ(xmin,xmax,ymin,ymax,zmin,zmax, &
     integral)
  use ModGaussQuad
  implicit none

  integer :: j
  real :: xr,xm,localdx
  real :: xmin,ymin,zmin
  real :: xmax,ymax,zmax
  real, dimension(1:3) :: tmp1,tmp2,integral

  xm=0.5*(xmax+xmin)
  xr=0.5*(xmax-xmin)
  integral(1)=0.0
  integral(2)=0.0
  integral(3)=0.0

  do j=1,nGaussNodes
     localdx=xr*GaussNodes(j)
     call body_force_GaussQuadYZ(xm+localdx,ymin,ymax,zmin,zmax,tmp1)
     call body_force_GaussQuadYZ(xm-localdx,ymin,ymax,zmin,zmax,tmp2)
     integral(1) = integral(1) + GaussWeights(j)*(tmp1(1)+tmp2(1))
     integral(2) = integral(2) + GaussWeights(j)*(tmp1(2)+tmp2(2))
     integral(3) = integral(3) + GaussWeights(j)*(tmp1(3)+tmp2(3))
  end do

  integral(1) = integral(1)*xr
  integral(2) = integral(2)*xr
  integral(3) = integral(3)*xr

end subroutine body_force_GaussQuadXYZ

!----------------------------------------------------------------------
!                        body_force_GaussQuadYZ
!----------------------------------------------------------------------

subroutine body_force_GaussQuadYZ(x0,ymin,ymax,zmin,zmax, &
     integral)
  use ModGaussQuad
  implicit none

  integer :: j
  real :: yr,ym,localdy
  real :: x0
  real :: ymin,zmin
  real :: ymax,zmax
  real, dimension(1:3) :: tmp1,tmp2,integral

  ym=0.5*(ymax+ymin)
  yr=0.5*(ymax-ymin)
  integral(1)=0.0
  integral(2)=0.0
  integral(3)=0.0

  do j=1,nGaussNodes
     localdy=yr*GaussNodes(j)
     call body_force_GaussQuadZ(x0,ym+localdy,zmin,zmax,tmp1)
     call body_force_GaussQuadZ(x0,ym-localdy,zmin,zmax,tmp2)
     integral(1) = integral(1) + GaussWeights(j)*(tmp1(1)+tmp2(1))
     integral(2) = integral(2) + GaussWeights(j)*(tmp1(2)+tmp2(2))
     integral(3) = integral(3) + GaussWeights(j)*(tmp1(3)+tmp2(3))
  end do

  integral(1) = integral(1)*yr
  integral(2) = integral(2)*yr
  integral(3) = integral(3)*yr

end subroutine body_force_GaussQuadYZ

!----------------------------------------------------------------------
!                      body_force_GaussQuadZ
!----------------------------------------------------------------------

subroutine body_force_GaussQuadZ(x0,y0,zmin,zmax,integral)
  use ModGaussQuad
  implicit none

  integer :: j
  real :: zr,zm,localdz
  real :: x0,y0
  real :: zmin,zmax
  real, dimension(1:3) :: tmp1,tmp2,integral

  zm=0.5*(zmax+zmin)
  zr=0.5*(zmax-zmin)

  integral(1)=0.0
  integral(2)=0.0
  integral(3)=0.0

  do j=1,nGaussNodes
     localdz=zr*GaussNodes(j)
     call gravity_force(x0,y0,zm+localdz,tmp1)
     call gravity_force(x0,y0,zm-localdz,tmp2)
     integral(1) = integral(1) + GaussWeights(j)*(tmp1(1)+tmp2(1))
     integral(2) = integral(2) + GaussWeights(j)*(tmp1(2)+tmp2(2))
     integral(3) = integral(3) + GaussWeights(j)*(tmp1(3)+tmp2(3))
     call centripetal_force(x0,y0,zm+localdz,tmp1)
     call centripetal_force(x0,y0,zm-localdz,tmp2)
     integral(1) = integral(1) + GaussWeights(j)*(tmp1(1)+tmp2(1))
     integral(2) = integral(2) + GaussWeights(j)*(tmp1(2)+tmp2(2))
     integral(3) = integral(3) + GaussWeights(j)*(tmp1(3)+tmp2(3))
  end do

  integral(1) = integral(1)*zr
  integral(2) = integral(2)*zr
  integral(3) = integral(3)*zr

end subroutine body_force_GaussQuadZ

!---------------------------------------------------------------------------
!                     heat_source_averages
!---------------------------------------------------------------------------

!\
! This routine computes the volume-averaged volumetric heating for 
! each of the cells on the grid.
!/

subroutine heat_source_averages
  use ModMain, ONLY : nI,nJ,nK,globalBLK
  use ModMain, ONLY : UseBody2                !^CFG IF SECONDBODY
  use ModAdvance, ONLY : qheat_BLK
  use ModPhysics, ONLY : Rbody
  use ModPhysics, ONLY : Rbody2               !^CFG IF SECONDBODY
  use ModGeometry, ONLY : R_BLK
  use ModGeometry, ONLY: R2_BLK               !^CFG IF SECONDBODY
  implicit none

  integer :: i,j,k
  real :: Rcrit, Rcrit2

  !true_cell note: using true_cell to replace an Rbody test does not apply here

  if (.not. UseBody2) then                 !^CFG IF SECONDBODY

     Rcrit = 0.50*Rbody

     qheat_BLK(:,:,:,globalBLK) = 0.00

     do k=1,nK
        do j=1,nJ
           do i=1,nI
              if ( R_BLK(i,j,k,globalBLK) > Rcrit .or. &
                   R_BLK(i-1,j,k,globalBLK) > Rcrit .or. &
                   R_BLK(i+1,j,k,globalBLK) > Rcrit .or. &
                   R_BLK(i,j-1,k,globalBLK) > Rcrit .or. &
                   R_BLK(i,j+1,k,globalBLK) > Rcrit .or. &
                   R_BLK(i,j,k-1,globalBLK) > Rcrit .or. &
                   R_BLK(i,j,k+1,globalBLK) > Rcrit ) then 
                 call heat_source_cell_average(i,j,k)
              end if
           end do !end i loop
        end do !end j loop
     end do !end k loop

  !^CFG IF SECONDBODY BEGIN
  else                      ! do the block if there are two bodies

     Rcrit  = 0.50*Rbody
     Rcrit2 = 0.50*RBody2

     qheat_BLK(:,:,:,globalBLK) = 0.00

     do k=1,nK
        do j=1,nJ
           do i=1,nI
              if ( (R_BLK(i,j,k,globalBLK) > Rcrit .or. &
                   R_BLK(i-1,j,k,globalBLK) > Rcrit .or. &
                   R_BLK(i+1,j,k,globalBLK) > Rcrit .or. &
                   R_BLK(i,j-1,k,globalBLK) > Rcrit .or. &
                   R_BLK(i,j+1,k,globalBLK) > Rcrit .or. &
                   R_BLK(i,j,k-1,globalBLK) > Rcrit .or. &
                   R_BLK(i,j,k+1,globalBLK) > Rcrit)     &
                   .and.                                   &
                   (R2_BLK(i,j,k,globalBLK) > Rcrit2 .or. &
                   R2_BLK(i-1,j,k,globalBLK) > Rcrit2 .or. &
                   R2_BLK(i+1,j,k,globalBLK) > Rcrit2 .or. &
                   R2_BLK(i,j-1,k,globalBLK) > Rcrit2 .or. &
                   R2_BLK(i,j+1,k,globalBLK) > Rcrit2 .or. &
                   R2_BLK(i,j,k-1,globalBLK) > Rcrit2 .or. &
                   R2_BLK(i,j,k+1,globalBLK) > Rcrit2) ) then 
                 call heat_source_cell_average(i,j,k)
              end if
           end do !end i loop
        end do !end j loop
     end do !end k loop

  end if
  !^CFG IF SECONDBODY END

end subroutine heat_source_averages

!---------------------------------------------------------------------------
!                     heat_source_cell_average
!---------------------------------------------------------------------------

!\
! Calculates the cell averages of the heating source
!/

subroutine heat_source_cell_average(i,j,k)
  use ModMain, ONLY : globalBLK
  use ModAdvance, ONLY : qheat_BLK
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,dx_BLK,dy_BLK,dz_BLK,VolumeInverse_I
  implicit none

  integer, intent(in) :: i,j,k

  integer :: ii
  real ::  xmin,xmax,ymin,ymax,zmin,zmax
  real :: q0ave

  xmin = x_BLK(i,j,k,globalBLK) - 0.50*dx_BLK(globalBLK)
  ymin = y_BLK(i,j,k,globalBLK) - 0.50*dy_BLK(globalBLK)
  zmin = z_BLK(i,j,k,globalBLK) - 0.50*dz_BLK(globalBLK)

  xmax = x_BLK(i,j,k,globalBLK) + 0.50*dx_BLK(globalBLK)
  ymax = y_BLK(i,j,k,globalBLK) + 0.50*dy_BLK(globalBLK)
  zmax = z_BLK(i,j,k,globalBLK) + 0.50*dz_BLK(globalBLK)

  call heat_source_GaussQuadXYZ(xmin,xmax,ymin,ymax,zmin,zmax,q0ave)

  q0ave = q0ave*VolumeInverse_I(globalBLK)

  qheat_BLK(i,j,k,globalBLK) = qheat_BLK(i,j,k,globalBLK) + q0ave

end subroutine heat_source_cell_average

!----------------------------------------------------------------------
!                          heat_source
!----------------------------------------------------------------------

!\
! Evaluates the non-dimensional heating source at the
! specified location (X0,Y0,Z0).
!/

subroutine heat_source(X0,Y0,Z0,q0)
  use ModMain, ONLY : UseBody2                !^CFG IF SECONDBODY
  use ModPhysics
  implicit none

  real, intent(in) :: X0,Y0,Z0
  real :: q0
!^CFG IF NOT SIMPLE BEGIN
  real :: R0,XT,YT,ZT,SIGMAHeat2
  real :: cosTheta, sinTheta, cosPhi, sinPhi, &
       sin2Theta_coronal_hole

  R0 = sqrt(X0*X0 + Y0*Y0 + Z0*Z0)
  R0 = max(R0,1.E-6)

  XT =  cosTHETAtilt*X0 + sinTHETAtilt*Z0
  YT = Y0
  ZT = -sinTHETAtilt*X0 + cosTHETAtilt*Z0

  cosTheta = ZT/R0
  sinTheta = sqrt(XT**2+YT**2)/R0
  cosPhi = XT/sqrt(XT**2+YT**2+cTolerance**2)
  sinPhi = YT/sqrt(XT**2+YT**2+cTolerance**2)

  call coronal_hole_boundary(R0, sin2Theta_coronal_hole)

  if (sinTheta*sinTheta < sin2Theta_coronal_hole) then
     if (R0 <= Rheat) then
        !         SIGMAHeat2 = 0.50*(Rheat - 1.00)
        SIGMAHeat2 = SIGMAheat
        q0 = Qsun*exp(-((R0-Rheat)/SIGMAheat2)**2)
     else
        SIGMAHeat2 = 2.00*SIGMAheat- &
             (sinTheta*sinTheta/sin2Theta_coronal_hole)* &
             (1.00*SIGMAheat)
        q0 = Qsun*exp(-((R0-Rheat)/SIGMAheat2)**2)
     end if
  else
     if (R0 <= Rheat) then
        !         SIGMAHeat2 = 0.50*(Rheat - 1.00)
        SIGMAHeat2 = SIGMAheat
        q0 = Qsun*exp(-((R0-Rheat)/SIGMAheat2)**2)
     else
        q0 = Qsun*exp(-((R0-Rheat)/SIGMAheat)**2)
     end if
  end if

  !^CFG IF SECONDBODY BEGIN
  !\
  ! If there is a second body for which there is a heat source
  ! add its contribution in the if block below
  !/

  if (UseBody2) then
  
  end if
  !^CFG IF SECONDBODY END

!^CFG END SIMPLE
end subroutine heat_source

!----------------------------------------------------------------------
!                       heat_source_GaussQuadXYZ
!----------------------------------------------------------------------

subroutine heat_source_GaussQuadXYZ(xmin,xmax,ymin,ymax,zmin,zmax, &
     integral)
  use ModGaussQuad
  implicit none

  integer :: j
  real :: xr,xm,localdx
  real :: xmin,ymin,zmin
  real :: xmax,ymax,zmax
  real :: tmp1,tmp2,integral

  xm=0.5*(xmax+xmin)
  xr=0.5*(xmax-xmin)
  integral=0.0

  do j=1,nGaussNodes
     localdx=xr*GaussNodes(j)
     call heat_source_GaussQuadYZ(xm+localdx,ymin,ymax,zmin,zmax,tmp1)
     call heat_source_GaussQuadYZ(xm-localdx,ymin,ymax,zmin,zmax,tmp2)
     integral = integral + GaussWeights(j)*(tmp1+tmp2)
  end do

  integral = integral*xr

end subroutine heat_source_GaussQuadXYZ

!----------------------------------------------------------------------
!                        heat_source_GaussQuadYZ
!----------------------------------------------------------------------

subroutine heat_source_GaussQuadYZ(x0,ymin,ymax,zmin,zmax,integral)
  use ModGaussQuad
  implicit none

  integer :: j
  real :: yr,ym,localdy
  real :: x0
  real :: ymin,zmin
  real :: ymax,zmax
  real :: tmp1,tmp2,integral

  ym=0.5*(ymax+ymin)
  yr=0.5*(ymax-ymin)
  integral=0.0

  do j=1,nGaussNodes
     localdy=yr*GaussNodes(j)
     call heat_source_GaussQuadZ(x0,ym+localdy,zmin,zmax,tmp1)
     call heat_source_GaussQuadZ(x0,ym-localdy,zmin,zmax,tmp2)
     integral = integral + GaussWeights(j)*(tmp1+tmp2)
  end do

  integral = integral*yr

end subroutine heat_source_GaussQuadYZ

!----------------------------------------------------------------------
!                      heat_source_GaussQuadZ
!----------------------------------------------------------------------

subroutine heat_source_GaussQuadZ(x0,y0,zmin,zmax,integral)
  use ModGaussQuad
  implicit none

  integer :: j
  real :: zr,zm,localdz
  real :: x0,y0
  real :: zmin,zmax
  real :: tmp1,tmp2,integral

  zm=0.5*(zmax+zmin)
  zr=0.5*(zmax-zmin)

  integral=0.0

  do j=1,nGaussNodes
     localdz=zr*GaussNodes(j)
     call heat_source(x0,y0,zm+localdz,tmp1)
     call heat_source(x0,y0,zm-localdz,tmp2)
     integral = integral + GaussWeights(j)*(tmp1+tmp2)
  end do

  integral = integral*zr

end subroutine heat_source_GaussQuadZ
