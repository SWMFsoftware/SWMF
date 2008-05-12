
!*************************************************************************
!
! MAGNETOSPHERE/IONOSPHERE coupling Routines
!
!*************************************************************************

!-------------------------------------------------------------------------
! ionosphere_currents
!
!
!
!-------------------------------------------------------------------------

subroutine ionosphere_currents(Jx,Jy,Jz,                                     &
                               Ex,Ey,Ez,ETh,EPs,                             &
                               Ux,Uy,Uz,                                     &
                               PHI, SigmaThTh, SigmaThPs, SigmaPsPs,         &
                               X, Y, Z,                                      &
                               Theta, Psi, Radius, nTheta, nPsi,             &
                               dTheta, dPsi)

  !\
  ! For the calculated ionospheric potential solution,
  ! this routine determines the ionospheric currents and
  ! electric fields, as well as convection velocities.
  !/

  use ModIonosphere
  implicit none

  integer :: nTheta, nPsi
  real :: Radius
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
                  PHI, SigmaThTh, SigmaThPs, SigmaPsPs, &
                  Jx,Jy,Jz, &
                  Ex,Ey,Ez,ETh,EPs, &
                  Ux,Uy,Uz, &
                  X, Y, Z, &
                  Theta, Psi

  real, dimension(1:IONO_nTheta) :: dTheta
  real, dimension(1:IONO_nPsi)   :: dPsi

  logical :: north
  integer :: i, j
  real :: cosTheta, sinTheta, cosPhi, sinPhi, &
          cs2, cs3, cs4, &
          ER, JR, JTh, JPs, &
          xx, yy, zz, &
          bx, by, bz, bR, bTh, bPs, BB, BBabs, &
          Vll, Vp_x, Vp_y, Vp_z, VR

  if (Theta(1,1) < 2.00*IONO_Theta_0) then
     north = .true.
  else
     north = .false.
  end if

  ! Compute the ionospheric electric field.

  do j = 1, nPsi

     if (j > 1 .and. j < nPsi ) then 
        do i = 2, nTheta-1
           sinTheta = sin(Theta(i,j))
           ETh(i,j) = -(PHI(i+1,j)-PHI(i-1,j))/                               &
                      (dTheta(i)*Radius)
           EPs(i,j) = -(PHI(i,j+1)-PHI(i,j-1))/                               &
                      (dPsi(j)*Radius*sinTheta)
        end do
        ETh(1,j) = -(PHI(2,j)-PHI(1,j))/                                      &
                   (dTheta(1)*Radius)
        EPs(1,j) = EPs(2,j)
        ETh(nTheta,j) = -(PHI(nTheta,j)-PHI(nTheta-1,j))/                     &
                        (dTheta(nTheta)*Radius)
        EPs(nTheta,j) = EPs(nTheta-1,j)
     else if (j == 1) then
        do i = 2, nTheta-1
           sinTheta = sin(Theta(i,j))
           ETh(i,j) = -(PHI(i+1,j)-PHI(i-1,j))/                               &
                      (dTheta(i)*Radius)
           EPs(i,j) = -(PHI(i,j+1)-PHI(i,nPsi-1))/                            &
                      (dPsi(j)*Radius*sinTheta)
        end do
        ETh(1,j) = -(PHI(2,j)-PHI(1,j))/                                      &
                   (dTheta(1)*Radius)
        EPs(1,j) = EPs(2,j)
        ETh(nTheta,j) = -(PHI(nTheta,j)-PHI(nTheta-1,j))/                     &
                        (dTheta(nTheta)*Radius)
        EPs(nTheta,j) = EPs(nTheta-1,j)
     else
        do i = 2, nTheta-1
           sinTheta = sin(Theta(i,j))
           ETh(i,j) = -(PHI(i+1,j)-PHI(i-1,j))/                               &
                      (dTheta(i)*Radius)
           EPs(i,j) = -(PHI(i,2)-PHI(i,j-1))/                                 &
                      (dPsi(j)*Radius*sinTheta)
        end do
        ETh(1,j) = -(PHI(2,j)-PHI(1,j))/                                      &
                   (dTheta(1)*Radius)
        EPs(1,j) = EPs(2,j)
        ETh(nTheta,j) = -(PHI(nTheta,j)-PHI(nTheta-1,j))/                     &
                        (dTheta(nTheta)*Radius)
        EPs(nTheta,j) = EPs(nTheta-1,j)
     end if
  end do

  ! Compute the ionospheric currents convection velocities.

  do j = 1, nPsi

     do i = 1, nTheta

        cosTheta = cos(Theta(i,j))
        sinTheta = sin(Theta(i,j))
        cosPhi = cos(Psi(i,j))
        sinPhi = sin(Psi(i,j))

        if (north .and. i == nTheta) then
           ER = 0.00
        else if (.not.north .and. i == 1) then
           ER = 0.00
        else
           ER = -0.50*(sinTheta/(cosTheta+IONO_Toler**2))*ETh(i,j)
        end if

        Ex(i,j) = ER*sinTheta*cosPhi + ETh(i,j)*cosTheta*cosPhi - &
                  EPs(i,j)*sinPhi
        Ey(i,j) = ER*sinTheta*sinPhi + ETh(i,j)*cosTheta*sinPhi + &
                  EPs(i,j)*cosPhi
        Ez(i,j) = ER*cosTheta - ETh(i,j)*sinTheta
        
        JR = 0.00
        JTh =  SigmaThTh(i,j)*ETh(i,j) + SigmaThPs(i,j)*EPs(i,j)
        JPs = -SigmaThPs(i,j)*ETh(i,j) + SigmaPsPs(i,j)*EPs(i,j)
        
        if (north) then
           IONO_NORTH_JTh(i,j) = JTh
           IONO_NORTH_JPs(i,j) = JPs
        else
           IONO_SOUTH_JTh(i,j) = JTh
           IONO_SOUTH_JPs(i,j) = JPs
        endif

        Jx(i,j) = JR*sinTheta*cosPhi + JTh*cosTheta*cosPhi - &
                  JPs*sinPhi
        Jy(i,j) = JR*sinTheta*sinPhi + JTh*cosTheta*sinPhi + &
                  JPs*cosPhi
        Jz(i,j) = JR*cosTheta - JTh*sinTheta
        
        cs2 = cosTheta*cosTheta
        cs3 = 1.00 + 3.00*cs2
        cs4 = sqrt(cs3)

        BB = IONO_Bdp*((IONO_Radius/Radius)**3)
        BBabs = abs(BB)
        bR = (BB/BBabs)*2.00*cosTheta/cs4
        bTh = (BB/BBabs)*sinTheta/cs4 
        bPs = 0.00  
        bx = bR*sinTheta*cosPhi + bTh*cosTheta*cosPhi - &
             bPs*sinPhi
        by = bR*sinTheta*sinPhi + bTh*cosTheta*sinPhi + &
             bPs*cosPhi
        bz = bR*cosTheta - bTh*sinTheta
        
        xx = sinTheta*cosPhi
        yy = sinTheta*sinPhi
        zz = cosTheta

        Vp_x = (Ey(i,j)*bz - Ez(i,j)*by)/BBabs
        Vp_y = (Ez(i,j)*bx - Ex(i,j)*bz)/BBabs
        Vp_z = (Ex(i,j)*by - Ey(i,j)*bx)/BBabs

!!$        VR = Vp_x*xx + Vp_y*yy + Vp_z*zz
!!$        if (north .and. i == nTheta) then
!!$           Vll = 0.00
!!$        else if (.not.north .and. i == 1) then
!!$           Vll = 0.00
!!$        else
!!$           Vll = -VR/(bx*xx+by*yy+bz*zz)
!!$        end if

        Vll = 0.00
        
        Ux(i,j) = Vp_x + Vll*bx
        Uy(i,j) = Vp_y + Vll*by
        Uz(i,j) = Vp_z + Vll*bz
     end do  
  end do

end subroutine ionosphere_currents

