!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine ionosphere_currents(iBlock, Jx,Jy,Jz,                             &
                               Ex,Ey,Ez,ETh,EPs,                             &
                               Ux,Uy,Uz,                                     &
                               PHI, SigmaThTh, SigmaThPs, SigmaPsPs,         &
                               X, Y, Z,                                      &
                               Theta, Psi,                                   &
                               dTheta, dPsi)

  !\
  ! For the calculated ionospheric potential solution,
  ! this routine determines the ionospheric currents and
  ! electric fields, as well as convection velocities.
  !/

  use ModIonosphere
  use ModCoordTransform, ONLY: dir_to_xyz, cross_product
  use CON_planet_field,  ONLY: get_planet_field
  use IE_ModMain,        ONLY: Time_Simulation

  implicit none

  integer, parameter :: nTheta = IONO_nTheta, nPsi = IONO_nPsi

  integer, intent(in) :: iBlock
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::  &
                  PHI, SigmaThTh, SigmaThPs, SigmaPsPs, &
                  Jx,Jy,Jz, &
                  Ex,Ey,Ez,ETh,EPs, &
                  Ux,Uy,Uz, &
                  X, Y, Z, &
                  Theta, Psi

  real, dimension(1:IONO_nTheta) :: dTheta
  real, dimension(1:IONO_nPsi)   :: dPsi

  integer :: i, j
  real :: cosTheta, sinTheta, cosPhi, sinPhi, &
          ER, JR, JTh, JPs, &
          Xyz_D(3), NormRadius, b_D(3), Vp_D(3)
  !----------------------------------------------------------------------------
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

        ! The radial gradient of the potential is related to the 
        ! theta gradient as dPhi/dr = dPhi/dTheta*dTheta/dr
        ! where dTheta/dr = sin(theta)/(2cos(theta)) along the 
        ! dipole field line which is an equipotential line.
        ! This formula diverges at the equator.
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

        ! Calculate location in Cartesian coordinates
        call dir_to_xyz(SinTheta,CosTheta,SinPhi,CosPhi,Xyz_D)
        Xyz_D = Xyz_D * (IONO_Radius + IONO_Height) / IONO_Radius
        ! Get magnetic field and normalize it to unity
        call get_planet_field(Time_Simulation,Xyz_D,'SMG NORM',b_D)

        ! Get potential V = E x B/|B^2|
        b_D = b_D/sum(b_D**2)
        Vp_D = cross_product((/Ex(i,j), Ey(i,j), Ez(i,j)/), b_D)

        Ux(i,j) = Vp_D(1)
        Uy(i,j) = Vp_D(2)
        Uz(i,j) = Vp_D(3)

     end do  
  end do

end subroutine ionosphere_currents
