!^CFG COPYRIGHT UM

!*************************************************************************
!
! MAGNETOSPHERE/IONOSPHERE coupling Routines
!
!*************************************************************************

!-------------------------------------------------------------------------
! ionosphere_fac
!
!
!
!-------------------------------------------------------------------------

subroutine ionosphere_fac(iBlock)
  !\
  ! Given the solution for the current from the magnetosphere,
  ! this routine maps the solution for the current down to the
  ! ionospheric boundary and evaluates the incoming and outgoing
  ! field-aligned currents (FAC).
  !/
  use ModIonosphere
  use IE_ModMain, ONLY: UseFullCurrent
  implicit none

  integer, intent(in) :: iBlock ! 1 for north, 2 for south

  integer :: i, j, k

  real, dimension(3) :: MagField_Orientation
  real :: FAC, dist_total, fac_total

  character(len=*), parameter :: NameSub='ionosphere_fac'
  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest,DoTestMe)
  if(DoTest)write(*,*)NameSub,': iBlock=',iBlock

  select case(iBlock)
  case(1)      ! Northern hemisphere

     if(allocated(MAG_NORTH_JR))deallocate(MAG_NORTH_JR)
     allocate(MAG_NORTH_JR(IONO_NORTH_nMagBndPts))

     do i = 1, IONO_NORTH_nMagBndPts

        MagField_Orientation(:) = MAG_NORTH_MagField(i,1:3)

        FAC = MAG_NORTH_Jx(i) * MagField_Orientation(1) + &
             MAG_NORTH_Jy(i) * MagField_Orientation(2) + &
             MAG_NORTH_Jz(i) * MagField_Orientation(3)

        ! -------------------------------------------------------
        ! Do we want to take the FULL current?
        !   If so, do the following
        ! -------------------------------------------------------

        if (UseFullCurrent) then

           ! Determine Sign

           if (FAC.ne.0) then 
              FAC = FAC / abs(FAC)
           else
              FAC = 1.0
           endif

           !
           ! Take total current
           !

           FAC = FAC * sqrt(MAG_NORTH_Jx(i)**2.0 + &
                MAG_NORTH_Jy(i)**2.0 + &
                MAG_NORTH_Jz(i)**2.0)

        else

           !
           ! Multiply times the cos of the dipole tilt angle
           !

           FAC = FAC * MAG_NORTH_MagField(i,5)

        endif

        !
        ! Multiply times the ratio of the magnetic fields

        MAG_NORTH_JR(i) = FAC * MAG_NORTH_MagField(i,4)

        if(DoTest.and.i==1)write(*,*)'ionosphere_fac FAC, MAG_NORTH_JR=',&
             FAC, MAG_NORTH_JR(i)

     end do

     ! Northern hemisphere

     IONO_NORTH_JR = 0.00

     do j = 1, IONO_nPsi
        do i = 1, IONO_nTheta

           k = 1
           dist_total = 1.0e-8
           fac_total = 0.0
           do while (k <= 8)
              if (IONO_NORTH_Mapping_Index(i,j,k) > 0) then
                 fac_total = fac_total +                                  &
                      MAG_NORTH_JR(IONO_NORTH_Mapping_Index(i,j,k)) *     &
                      (1.0 / (IONO_NORTH_Mapping_Distance(i,j,k) + 1.0e-8))
                 dist_total = dist_total + &
                      1.0 / (IONO_NORTH_Mapping_Distance(i,j,k) + 1.0e-8)
                 k = k + 1
              else
                 k = 101
              endif
           enddo

           IONO_NORTH_JR(i,j) = fac_total / dist_total

           !
           ! If we only have 1 mapping point within our field of view, then
           ! let's extrapolate it out as a function of distance
           ! We multiply times dist_total because this is actually 1/distance
           !

           if (k == 101) IONO_NORTH_JR(i,j) = &
                IONO_NORTH_JR(i,j)*(1.0-(1.0/dist_total)/Max_FAC_Distance)

        enddo
     enddo

     ! Ensure that the FAC is periodic.

     IONO_NORTH_JR(1:IONO_nTheta,IONO_nPsi) = &
          0.50*(IONO_NORTH_JR(1:IONO_nTheta,IONO_nPsi)+ &
          IONO_NORTH_JR(1:IONO_nTheta,1)) 
     IONO_NORTH_JR(1:IONO_nTheta,1) = &
          IONO_NORTH_JR(1:IONO_nTheta,IONO_nPsi)


     if(DoTest)then

        write(*,*)NameSub,': sum(abs(MAG_NORTH_Jxyz))=',&
             sum(abs(MAG_NORTH_Jx)), &
             sum(abs(MAG_NORTH_Jy)), &
             sum(abs(MAG_NORTH_Jz))

        write(*,*)NameSub,': sum(abs(IONO_NORTH_Mapping_Index))=',&
             sum(abs(IONO_NORTH_Mapping_Index))

        write(*,*)NameSub,': sum(abs(IONO_NORTH_Mapping_Distance))=',&
             sum(abs(IONO_NORTH_Mapping_Distance))

        write(*,*)NameSub,': sum(abs(MAG_NORTH_MagField))=',&
             sum(abs(MAG_NORTH_MagField))

        write(*,*)NameSub,': sum(abs(IONO_NORTH_JR))=',&
             sum(abs(IONO_NORTH_JR))

     end if

  case(2)      ! Southern hemisphere

     if(allocated(MAG_SOUTH_JR))deallocate(MAG_SOUTH_JR)
     allocate(MAG_SOUTH_JR(IONO_SOUTH_nMagBndPts) )

     do i = 1, IONO_SOUTH_nMagBndPts

        MagField_Orientation(:) = MAG_SOUTH_MagField(i,1:3)

        FAC = MAG_SOUTH_Jx(i) * MagField_Orientation(1) + &
             MAG_SOUTH_Jy(i) * MagField_Orientation(2) + &
             MAG_SOUTH_Jz(i) * MagField_Orientation(3)

        ! -------------------------------------------------------
        ! Do we want to take the FULL current?
        !   If so, do the following
        ! -------------------------------------------------------

        if (UseFullCurrent) then

           ! Determine Sign

           if (FAC.ne.0) then 
              FAC = FAC / abs(FAC)
           else
              FAC = 1.0
           endif

           !
           ! Take total current
           !

           FAC = FAC * sqrt(MAG_SOUTH_Jx(i)**2.0 + &
                MAG_SOUTH_Jy(i)**2.0 + &
                MAG_SOUTH_Jz(i)**2.0)

        else

           !
           ! Multiply times the cos of the dipole tilt angle
           !

           FAC = FAC * MAG_SOUTH_MagField(i,5)

        endif

        !
        ! Multiply times the ratio of the magnetic fields

        MAG_SOUTH_JR(i) = FAC * MAG_SOUTH_MagField(i,4)

     end do

     ! Southern hemisphere

     IONO_SOUTH_JR = 0.00

     do j = 1, IONO_nPsi
        do i = 1, IONO_nTheta

           k = 1
           dist_total = 1.0e-8
           fac_total = 0.0
           do while (k <= 8)
              if (IONO_SOUTH_Mapping_Index(i,j,k) > 0) then
                 fac_total = fac_total +                                     &
                      MAG_SOUTH_JR(IONO_SOUTH_Mapping_Index(i,j,k)) *        &
                      (1.0 / (IONO_SOUTH_Mapping_Distance(i,j,k) + 1.0e-8))
                 dist_total = dist_total + &
                      1.0 / (IONO_SOUTH_Mapping_Distance(i,j,k) + 1.0e-8)
                 k = k + 1
              else
                 k = 101
              endif
           enddo

           IONO_SOUTH_JR(i,j) = fac_total / dist_total
           if (k == 101) IONO_SOUTH_JR(i,j) = &
                IONO_SOUTH_JR(i,j)*(1.0-(1.0/dist_total)/Max_FAC_Distance)

        enddo
     enddo

     ! Ensure that the FAC is periodic.
     IONO_SOUTH_JR(1:IONO_nTheta,IONO_nPsi) = &
          0.50*(IONO_SOUTH_JR(1:IONO_nTheta,IONO_nPsi)+ &
          IONO_SOUTH_JR(1:IONO_nTheta,1)) 
     IONO_SOUTH_JR(1:IONO_nTheta,1) = &
          IONO_SOUTH_JR(1:IONO_nTheta,IONO_nPsi)


     if(DoTest)then

        write(*,*)NameSub,': sum(abs(MAG_SOUTH_Jxyz))=',&
             sum(abs(MAG_SOUTH_Jx)), &
             sum(abs(MAG_SOUTH_Jy)), &
             sum(abs(MAG_SOUTH_Jz))

        write(*,*)NameSub,': sum(abs(IONO_SOUTH_Mapping_Index))=',&
             sum(abs(IONO_SOUTH_Mapping_Index))

        write(*,*)NameSub,': sum(abs(IONO_SOUTH_Mapping_Distance))=',&
             sum(abs(IONO_SOUTH_Mapping_Distance))

        write(*,*)NameSub,': sum(abs(MAG_SOUTH_MagField))=',&
             sum(abs(MAG_SOUTH_MagField))

        write(*,*)NameSub,': sum(abs(IONO_SOUTH_JR))=',&
             sum(abs(IONO_SOUTH_JR))

     end if


  end select

end subroutine ionosphere_fac

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
  use ModCoordTransform, ONLY: dir_to_xyz, cross_product
  use CON_physics, ONLY: get_planet_field
  use IE_ModMain, ONLY: Time_Simulation

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
          ER, JR, JTh, JPs, &
          Xyz_D(3), NormRadius, b_D(3), Vp_D(3), Babs

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

        ! Calculate location in Cartesian coordinates
        call dir_to_xyz(SinTheta,CosTheta,SinPhi,CosPhi,Xyz_D)
        Xyz_D = Xyz_D * (IONO_Radius + IONO_Height) / IONO_Radius
        ! Get magnetic field and normalize it to unity
        call get_planet_field(Time_Simulation,Xyz_D,'SMG NORM',b_D)
        bAbs = sqrt(sum(b_D**2))
        b_D = b_D/bAbs
        ! Get potential V = E x B/|B|
        Vp_D = cross_product((/Ex(i,j), Ey(i,j), Ez(i,j)/), b_D)

        Ux(i,j) = Vp_D(1)
        Uy(i,j) = Vp_D(2)
        Uz(i,j) = Vp_D(3)


!!!!!  This is all wrong at this point. CON_physics needs to fix this.

!!!!!        
!!!!!        xx = NormRadius * sinTheta*cosPhi
!!!!!        yy = NormRadius * sinTheta*sinPhi
!!!!!        zz = NormRadius * cosTheta
!!!!!
!!!!!
!!!!!
!!!!!        BB = IONO_Bdp*((IONO_Radius/Radius)**3)
!!!!!        BBabs = abs(BB)
!!!!!        bR = (BB/BBabs)*2.00*cosTheta/cs4
!!!!!        bTh = (BB/BBabs)*sinTheta/cs4 
!!!!!        bPs = 0.00  
!!!!!        bx = bR*sinTheta*cosPhi + bTh*cosTheta*cosPhi - &
!!!!!             bPs*sinPhi
!!!!!        by = bR*sinTheta*sinPhi + bTh*cosTheta*sinPhi + &
!!!!!             bPs*cosPhi
!!!!!        bz = bR*cosTheta - bTh*sinTheta

!        Vp_x = (Ey(i,j)*bz - Ez(i,j)*by)/BBabs
!        Vp_y = (Ez(i,j)*bx - Ex(i,j)*bz)/BBabs
!        Vp_z = (Ex(i,j)*by - Ey(i,j)*bx)/BBabs

     end do  
  end do

end subroutine ionosphere_currents




