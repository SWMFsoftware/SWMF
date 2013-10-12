!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine Create_Region2_Currents(iBlock)

  use ModIonosphere
  implicit none

  integer,intent(in) :: iBlock ! 1 for north, 2 for south !!! SWMF

  real, dimension(1:IONO_nPsi) :: &
       Loc_of_Oval, Width_of_Oval, Strength_of_Oval

  real :: distance, center, strength, width

  integer :: i,j

  real :: total_up, total_down, factor, area

  !---------------------------------------------------------------------------
  !
  ! We have a slight problem with the currents within the MHD code.  When
  ! the IMF By component gets rather large, the current flowing into each
  ! hemisphere does not balance. If the Earth did not exist, the MHD code
  ! would work quite well; in that the total current flowing into and out
  ! of the ionosphere over the whole Earth balances exactly.  Because we
  ! put an inner boundary on the problem and do not let the current flow
  ! between the hemispheres, we have a problem.  This is true for the real
  ! Earth also - very little current flows between the hemispheres.
  !
  ! We have devised a way to force the current flowing in and out to be
  ! approximately equal.  Namely, we add up the current flowing in and
  ! out, determine which one needs a little more and by how much, then we
  ! multiply the current which needs more by a factor, allowing it to raise
  ! to approximately the same value as the other current.
  !
  ! This is physically unjustified.
  !

  select case(iBlock)

  case(1) ! This is for north:

     IONO_NORTH_Fake_JR = 0.0
     total_up   = 0.0
     total_down = 0.0
     do i=1, IONO_nTheta
        do j=1, IONO_nPsi
           area = dTheta_North(i) * Radius *   &
                dPsi_North(j) * sin(IONO_NORTH_Theta(i,j)) * Radius
           if (IONO_NORTH_JR(i,j) < 0.0) then
              total_down = total_down + abs(IONO_NORTH_JR(i,j)) * area
           else
              total_up = total_up + IONO_NORTH_JR(i,j) * area
           endif
        enddo
     enddo
     if (total_down > total_up) then
        factor = (total_down - total_up)/total_up
        where (IONO_NORTH_JR > 0.0) 
           IONO_NORTH_Fake_JR = IONO_NORTH_JR * factor
        endwhere
     else
        factor = (total_up - total_down)/total_down
        where (IONO_NORTH_JR < 0.0) 
           IONO_NORTH_Fake_JR = IONO_NORTH_JR * factor
        endwhere
     endif

     ! The MHD code does not do well with Region 2 currents. These currents
     ! are created by pressure gradients in the inner magnetosphere which may
     ! be overwhelmed by other features of the code, such as the inner
     ! boundary condition on pressure and density.  If that is the case, then
     ! we may need to create fake region 2 currents.  These currents are
     ! basically just equatorward of the region 1 currents and have the
     ! opposite sense - namely up on the dawn side and down on the dusk side.
     ! We approximate these by a sin function.

     call Determine_Oval_Characteristics(                   &
          IONO_NORTH_JR, IONO_NORTH_Theta, IONO_NORTH_Psi,  &
          Loc_of_Oval, Width_of_Oval, Strength_of_Oval)

     center = Loc_of_Oval(1) + 15.0*cDegToRad
     width  = 2.5*cDegToRad
     strength = -0.5*Strength_of_Oval(1)

     do i=1,IONO_nTheta
        if (IONO_NORTH_JR(i,IONO_nPsi/4).ne.0.0) &
             center = IONO_NORTH_Theta(i,1)
     enddo

     center = center - 5.0*cDegToRad

     do i=1, IONO_nTheta
        do j=1, IONO_nPsi

           distance = IONO_NORTH_Theta(i,j)-center

           if (abs(distance) < width*2.0) then

              IONO_NORTH_Fake_JR(i,j) =                               &
                   strength * sin(IONO_NORTH_Psi(i,j)) *        &
                   exp(-1.0*(distance/width)**2)

           endif

        enddo
     enddo

  case(2) ! This is for south:

     IONO_SOUTH_Fake_JR = 0.0
     total_up   = 0.0
     total_down = 0.0
     do i=1, IONO_nTheta
        do j=1, IONO_nPsi
           area = dTheta_South(i) * Radius *   &
                dPsi_South(j) * sin(IONO_SOUTH_Theta(i,j)) * Radius
           if (IONO_SOUTH_JR(i,j) < 0.0) then
              total_down = total_down + abs(IONO_SOUTH_JR(i,j)) * area
           else
              total_up = total_up + IONO_SOUTH_JR(i,j) * area
           endif
        enddo
     enddo
     if (total_down > total_up) then
        factor = (total_down - total_up)/total_up
        where (IONO_SOUTH_JR > 0.0) 
           IONO_SOUTH_Fake_JR = IONO_SOUTH_JR * factor
        endwhere
     else
        factor = (total_up - total_down)/total_down
        where (IONO_SOUTH_JR < 0.0) 
           IONO_SOUTH_Fake_JR = IONO_SOUTH_JR * factor
        endwhere
     endif

     call Determine_Oval_Characteristics(                   &
          IONO_SOUTH_JR, IONO_SOUTH_Theta, IONO_SOUTH_Psi,  &
          Loc_of_Oval, Width_of_Oval, Strength_of_Oval)

     center = cPi - (Loc_of_Oval(1) + 15.0*cDegToRad)
     strength = -0.5*Strength_of_Oval(1)

     do i=IONO_nTheta,1,-1
        if (IONO_SOUTH_JR(i,IONO_nPsi/4).ne.0.0) &
             center = IONO_SOUTH_Theta(i,1)
     enddo

     center = center + 5.0*cDegToRad

     do i=1, IONO_nTheta
        do j=1, IONO_nPsi

           distance = IONO_SOUTH_Theta(i,j)-center

           if (abs(distance) < width*2.0) then

              IONO_SOUTH_Fake_JR(i,j) =                               &
                   strength * sin(IONO_SOUTH_Psi(i,j)) *        &
                   exp(-1.0*(distance/width)**2)

           endif

        enddo
     enddo
  end select
end subroutine Create_Region2_Currents
