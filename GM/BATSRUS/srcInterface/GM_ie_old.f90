! This file contains the obsolete algorithms for GM-IE coupling

!==========================================================================
module ModMapPotential
  use ModInnerBC
  use ModNumConst
  use ModUtilities, ONLY: check_allocate
  implicit none
  save
  real, dimension(:,:), allocatable ::     &
       IONO_NORTH_PHI,IONO_SOUTH_PHI,               & !Ionospheric potential
       IONO_NORTH_X,IONO_NORTH_Y,IONO_NORTH_Z,      & !Ionospheric coordinates
       IONO_SOUTH_X,IONO_SOUTH_Y,IONO_SOUTH_Z         !Ionospheric coordinates
 
  integer, dimension(:,:), allocatable ::             &
       IONOtoMAG_NORTH_hemisphere,IONOtoMAG_SOUTH_hemisphere, &  
       IONOtoMAG_NORTH_lat_index, IONOtoMAG_SOUTH_lat_index,   &
       IONOtoMAG_NORTH_lon_index, IONOtoMAG_SOUTH_lon_index
  
  real, dimension(:,:), allocatable ::                &
       IONOtoMAG_NORTH_lat_factor, IONOtoMAG_SOUTH_lat_factor, &
       IONOtoMAG_NORTH_lon_factor, IONOtoMAG_SOUTH_lon_factor
contains
  subroutine init_map_potential_arrays
    integer::iError
    allocate( IONO_NORTH_PHI(IONO_nTheta,IONO_nPsi),&
              IONO_SOUTH_PHI(IONO_nTheta,IONO_nPsi),& 
              IONO_NORTH_X(IONO_nTheta,IONO_nPsi),  &
              IONO_NORTH_Y(IONO_nTheta,IONO_nPsi),  &
              IONO_NORTH_Z(IONO_nTheta,IONO_nPsi),  & 
              IONO_SOUTH_X(IONO_nTheta,IONO_nPsi),  &
              IONO_SOUTH_Y(IONO_nTheta,IONO_nPsi),  &
              IONO_SOUTH_Z(IONO_nTheta,IONO_nPsi),  & 
              IONOtoMAG_NORTH_hemisphere(IONO_nTheta,IONO_nPsi),&
              IONOtoMAG_SOUTH_hemisphere(IONO_nTheta,IONO_nPsi),&
              IONOtoMAG_NORTH_lat_index(IONO_nTheta,IONO_nPsi), &
              IONOtoMAG_SOUTH_lat_index(IONO_nTheta,IONO_nPsi), &
              IONOtoMAG_NORTH_lon_index(IONO_nTheta,IONO_nPsi), &
              IONOtoMAG_SOUTH_lon_index(IONO_nTheta,IONO_nPsi), &
              IONOtoMAG_NORTH_lat_factor(IONO_nTheta,IONO_nPsi),&
              IONOtoMAG_SOUTH_lat_factor(IONO_nTheta,IONO_nPsi),&
              IONOtoMAG_NORTH_lon_factor(IONO_nTheta,IONO_nPsi),&
              IONOtoMAG_SOUTH_lon_factor(IONO_nTheta,IONO_nPsi),&
              stat=iError)
    IONO_NORTH_PHI=cZero
    IONO_SOUTH_PHI=cZero 
    IONO_NORTH_X=cZero
    IONO_NORTH_Y=cZero
    IONO_NORTH_Z=cZero 
    IONO_SOUTH_X=cZero
    IONO_SOUTH_Y=cZero
    IONO_SOUTH_Z=cZero 
    IONOtoMAG_NORTH_hemisphere=cZero
    IONOtoMAG_SOUTH_hemisphere=cZero
    IONOtoMAG_NORTH_lat_index=cZero
    IONOtoMAG_SOUTH_lat_index=cZero 
    IONOtoMAG_NORTH_lon_index=cZero
    IONOtoMAG_SOUTH_lon_index=cZero
    IONOtoMAG_NORTH_lat_factor=cZero
    IONOtoMAG_SOUTH_lat_factor=cZero
    IONOtoMAG_NORTH_lon_factor=cZero
    IONOtoMAG_SOUTH_lon_factor=cZero
    call check_allocate(iError,'Arrays for mapping potential')
  end subroutine init_map_potential_arrays
end module ModMapPotential

!============================================================================

module ModInnerBC

  use ModMappingParam
  implicit none
  save

  !\
  ! Ionosphere array parameters
  !/
  integer :: IONO_nTheta = -1
  integer :: IONO_nPsi   = -1

  !\
  ! Ionosphere solution array definitions
  !/
  real, dimension(:,:), allocatable ::     &
       IONO_NORTH_Theta,IONO_NORTH_Psi,             & !
       IONO_SOUTH_Theta,IONO_SOUTH_Psi


  real, dimension(:,:), allocatable ::     &
       IONO_NORTH_PHI_BC,IONO_SOUTH_PHI_BC,         & !Magnetospheric bound pot
       IONO_NORTH_ETh_BC,IONO_NORTH_EPs_BC,         & !Magnetospheric bound
       IONO_SOUTH_ETh_BC,IONO_SOUTH_EPs_BC,         & !Electric fields
       IONO_NORTH_UTh_BC,IONO_NORTH_UPs_BC,         & !Magnetosphere bound flow
       IONO_NORTH_UR_BC,                            & !Velocities
       IONO_SOUTH_UTh_BC,IONO_SOUTH_UPs_BC,         &
       IONO_SOUTH_UR_BC

end module ModInnerBC

!PHI->E->U
!============================================================================

subroutine GM_calc_iono_bcs

  ! Calculate the inner BC velocities 
  ! from the electrostatic potential PHI provided by IE on PE 0.
  ! Broadcast velocities to all PEs.

  use ModProcMH
  use ModMapPotential,ONLY:&
       IONO_nTHETA,IONO_nPSI,IONO_Radius_Mag_Boundary,&
       IONO_NORTH_Phi_BC, IONO_SOUTH_Phi_BC,         & 
       IONO_NORTH_ETh_BC, & 
       IONO_NORTH_EPs_BC, & 
       IONO_NORTH_UR_BC, &  
       IONO_NORTH_UTh_BC, & 
       IONO_NORTH_UPs_BC, & 
       IONO_SOUTH_ETh_BC, & 
       IONO_SOUTH_EPs_BC, & 
       IONO_SOUTH_UR_BC, &  
       IONO_SOUTH_UTh_BC, & 
       IONO_SOUTH_UPs_BC,&  
       IONO_NORTH_X,&       
       IONO_NORTH_Y,&       
       IONO_NORTH_Z,&       
       IONO_SOUTH_X,&       
       IONO_SOUTH_Y,&       
       IONO_SOUTH_Z,&       
       IONO_SOUTH_PSI,IONO_SOUTH_THETA,IONO_SOUTH_PHI,& !Allocated above
       IONO_NORTH_PSI,IONO_NORTH_THETA,IONO_NORTH_PHI   !Allocated above
  use ModMpi
  implicit none

  integer :: isize, iError

  logical :: IsInitialized = .false.

  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------

  call set_oktest('calc_iono_bcs',DoTest,DoTestMe)

  if(DoTestMe)write(*,*)'calc_iono_bcs starting'

  !\
  ! Send the predicted ionospheric convection velocities to
  ! all PEs for use in the application of the ionospheric boundary
  ! conditions of magnetosphere calculations.
  !/
  isize = IONO_nTheta*IONO_nPsi

  if (.not.IsInitialized)then
     if(DoTestMe)write(*,*)'calc_iono_bcs initializing'

     call MPI_Bcast(IONO_NORTH_Theta(1,1), &
          isize,MPI_REAL,0,iComm,iError)
     call MPI_Bcast(IONO_NORTH_Psi(1,1), &
          isize,MPI_REAL,0,iComm,iError) 
     call MPI_Bcast(IONO_SOUTH_Theta(1,1), &
          isize,MPI_REAL,0,iComm,iError)
     call MPI_Bcast(IONO_SOUTH_Psi(1,1), &
          isize,MPI_REAL,0,iComm,iError)

     if(DoTestMe)then
        write(*,*)'IONO_nTheta,IONO_nPsi=',IONO_nTheta,IONO_nPsi
     end if

     if(DoTestMe)write(*,*)'calc_iono_bcs call GM_map_grid'
     call GM_map_grid

     IsInitialized=.true.
  end if

  if(iProc==0)then
     if(DoTest)write(*,*)'calc_iono_bcs call ionosphere_magBCs for north'

     if(DoTestMe)write(*,*)&
          allocated(IONO_NORTH_ETh_BC),&
          allocated(IONO_NORTH_EPs_BC),&
          allocated(IONO_NORTH_UR_BC),&
          allocated(IONO_NORTH_UTh_BC),&
          allocated(IONO_NORTH_UPs_BC),&
          IONO_Radius_Mag_Boundary,&
          allocated(IONO_NORTH_PHI),&
          allocated(IONO_NORTH_X),&
          allocated(IONO_NORTH_Y),&
          allocated(IONO_NORTH_Z),&
          allocated(IONO_NORTH_Theta),maxval(abs(IONO_NORTH_Theta)),&
          allocated(IONO_NORTH_Psi),maxval(abs(IONO_NORTH_Psi)),&
          IONO_nTheta, IONO_nPsi

!Transform 

     call map_ie_potential_to_gm(IONO_NORTH_PHI,&
          IONO_Radius_Mag_Boundary,&
          IONO_NORTH_X, IONO_NORTH_Y, IONO_NORTH_Z,& 
          IONO_NORTH_Theta, IONO_NORTH_Psi, &
          IONO_nTheta, IONO_nPsi,&
          IONO_NORTH_Phi_BC) 

     call ionosphere_magBCs( IONO_NORTH_Phi_BC , IONO_NORTH_ETh_BC, &
          IONO_NORTH_EPs_BC, &
          IONO_NORTH_UR_BC, IONO_NORTH_UTh_BC, &
          IONO_NORTH_UPs_BC, &
          IONO_Radius_Mag_Boundary,&
          IONO_NORTH_Theta, IONO_NORTH_Psi, &
          IONO_nTheta, IONO_nPsi)

     if(DoTestMe)write(*,*)'calc_iono_bcs north sum(abs(UR,UTh,UPs_BC))=',&
          sum(abs(IONO_NORTH_UR_BC)),&
          sum(abs(IONO_NORTH_UTh_BC)),&
          sum(abs(IONO_NORTH_UPs_BC))

     if(DoTestMe)write(*,*)&
          allocated(IONO_SOUTH_ETh_BC),&
          allocated(IONO_SOUTH_EPs_BC),&
          allocated(IONO_SOUTH_UR_BC),&
          allocated(IONO_SOUTH_UTh_BC),&
          allocated(IONO_SOUTH_UPs_BC),&
          IONO_Radius_Mag_Boundary,&
          allocated(IONO_SOUTH_PHI),&
          allocated(IONO_SOUTH_X),&
          allocated(IONO_SOUTH_Y),&
          allocated(IONO_SOUTH_Z),&
          allocated(IONO_SOUTH_Theta),maxval(abs(IONO_SOUTH_Theta)),&
          allocated(IONO_SOUTH_Psi),maxval(abs(IONO_SOUTH_Psi)),&
          IONO_nTheta, IONO_nPsi


     call map_ie_potential_to_gm(IONO_SOUTH_PHI,&
          IONO_Radius_Mag_Boundary,&
          IONO_SOUTH_X, IONO_SOUTH_Y, IONO_SOUTH_Z,& 
          IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
          IONO_nTheta, IONO_nPsi,&
          IONO_SOUTH_Phi_BC) 

     call ionosphere_magBCs( IONO_SOUTH_Phi_BC, IONO_SOUTH_ETh_BC, &
          IONO_SOUTH_EPs_BC, &
          IONO_SOUTH_UR_BC, IONO_SOUTH_UTh_BC, &
          IONO_SOUTH_UPs_BC, &
          IONO_Radius_Mag_Boundary,&
          IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
          IONO_nTheta, IONO_nPsi)

     if(DoTestMe)write(*,*)'calc_iono_bcs south sum(abs(UR,UTh,UPs_BC))=',&
          sum(abs(IONO_SOUTH_UR_BC)),&
          sum(abs(IONO_SOUTH_UTh_BC)),&
          sum(abs(IONO_SOUTH_UPs_BC))
  endif

  if(DoTestMe)write(*,*)'calc_iono_bcs broadcasting velocities'

  call MPI_Bcast(IONO_NORTH_UR_BC(1,1), &
       isize,MPI_REAL,0,iComm,iError)
  call MPI_Bcast(IONO_NORTH_UTh_BC(1,1), &
       isize,MPI_REAL,0,iComm,iError)
  call MPI_Bcast(IONO_NORTH_UPs_BC(1,1), &
       isize,MPI_REAL,0,iComm,iError)
  call MPI_Bcast(IONO_SOUTH_UR_BC(1,1), &
       isize,MPI_REAL,0,iComm,iError)
  call MPI_Bcast(IONO_SOUTH_UTh_BC(1,1), &
       isize,MPI_REAL,0,iComm,iError)
  call MPI_Bcast(IONO_SOUTH_UPs_BC(1,1), &
       isize,MPI_REAL,0,iComm,iError)

  if(DoTestMe)write(*,*)'calc_iono_bcs finished'

end subroutine GM_calc_iono_bcs
!=============================================================================
subroutine GM_map_grid

  ! Map the ionosphere grid to the inner boundary of GM

  use ModMapPotential,ONLY:&
       IONO_NORTH_THETA,IONO_NORTH_PSI,& 
       IONO_SOUTH_THETA,IONO_SOUTH_PSI,& 
       IONOtoMAG_NORTH_hemisphere, & 
       IONOtoMAG_NORTH_lat_index , & 
       IONOtoMAG_NORTH_lon_index , & 
       IONOtoMAG_NORTH_lat_factor, & 
       IONOtoMAG_NORTH_lon_factor, & 
       IONOtoMAG_SOUTH_hemisphere, & 
       IONOtoMAG_SOUTH_lat_index, &  
       IONOtoMAG_SOUTH_lon_index, &  
       IONOtoMAG_SOUTH_lat_factor, & 
       IONOtoMAG_SOUTH_lon_factor,&  
       IONO_nTheta,IONO_nPsi, &
       IONO_Radius_Mag_Boundary,IONO_Radius,IONO_PI
  implicit none

  logical :: done

  integer :: n, i, j, k, l, i0, j0
  integer :: hem, l_lo, l_up, k_lo, k_up

  real :: Iono_Psi, Iono_Theta, r, lo, la
  real,  dimension(1:3) :: Mag_Loc, Iono_Loc

  ! Automatic arrays !!!
  real, dimension(1:IONO_nTheta,1:IONO_nPsi,2) ::                            &
                  Theta_temp, Psi_temp, phi_temp


  logical :: DoTest, DoTestMe, DoTestThis
  ! ---------------------------------------------------------------------------

  call set_oktest('map_phi',DoTest,DoTestMe)
  if(DoTestMe)then
     write(*,*)'map_phi starting with ',&
       'iono_radius_mag_boundary, iono_radius=',&
       iono_radius_mag_boundary, iono_radius
     write(*,*)'map_phi north sum(abs(theta,psi))=',&
          sum(abs(IONO_NORTH_THETA)), sum(abs(IONO_NORTH_PSI))
     write(*,*)'map_phi south sum(abs(theta,psi))=',&
          sum(abs(IONO_SOUTH_THETA)), sum(abs(IONO_SOUTH_PSI))
  end if

 
  i0 = 1
  j0 = 1

  Theta_temp(:,:,1) = IONO_NORTH_THETA
  Psi_temp(:,:,1)   = IONO_NORTH_PSI

  Theta_temp(:,:,2) = IONO_SOUTH_THETA
  Psi_temp(:,:,2)   = IONO_SOUTH_PSI

  if(DoTestMe)write(*,*)'map_phi initialized temp variables'

  do n = 1, 2

     do j = 1, IONO_nPsi
        do i = 1, IONO_nTheta

           DoTestThis = DoTestMe .and. j==1 .and. i==2

           if(DoTestThis)write(*,*)'starting n,j,i=',n,j,i
           ! Find the X,Y,Z coordinate of the point on a sphere in the
           ! magnetosphere.

           Mag_Loc(1) = IONO_Radius_Mag_Boundary * sin(Theta_temp(i,j,n))* &
                                    cos(Psi_temp(i,j,n)) / IONO_Radius
           Mag_Loc(2) = IONO_Radius_Mag_Boundary * sin(Theta_temp(i,j,n))* &
                                    sin(Psi_temp(i,j,n)) / IONO_Radius
           Mag_Loc(3) = IONO_Radius_Mag_Boundary * &
                                  cos(Theta_temp(i,j,n)) / IONO_Radius

           if(DoTestThis)&
                write(*,*)'call Get_Mapping_Point with Mag_Loc=',Mag_Loc

           ! This routine finds the X,Y,Z point in the ionosphere on the
           ! same field line.  If it returns a zero radius, that means
           ! that it took one step and was going in the wrong direction,
           ! therefore the trace has to be completed again, but with the
           ! opposite sign of the direction.

           call Get_Mapping_Point(Mag_Loc, 1.0, Iono_Loc)
           
           if(DoTestThis)&
                write(*,*)'return Get_Mapping_Point with Iono_Loc=',Iono_Loc


           ! Figure out the radius to see if we need to do the trace again.

           r = sqrt(Iono_Loc(1)**2 + Iono_Loc(2)**2 + Iono_Loc(3)**2)



           if (r > 0.0) then

              ! This means that the point is in the northern hemisphere.

              hem = 1

              ! Find the colatitude

              if (Iono_Loc(1) == 0.00 .and. Iono_Loc(2) == 0.00) then
                 Iono_Psi = 0.00
              else
                 Iono_Psi = atan2(Iono_Loc(2),Iono_Loc(1))
              end if
              if (Iono_Psi < 0.00) Iono_Psi = Iono_Psi + 2.00*IONO_PI

              ! Find the longitude

              r = sqrt(Iono_Loc(1)*Iono_Loc(1) + Iono_Loc(2)*Iono_Loc(2))
              Iono_Theta = atan2(r, Iono_Loc(3))

           else

              ! This means that the point is in the northern hemisphere.

              hem = 2
           
              ! Redo the trace.

              call Get_Mapping_Point(Mag_Loc, -1.0, Iono_Loc)

              if(DoTestThis)&
                   write(*,*)'2nd Get_Mapping_Point with Iono_Loc=',Iono_Loc

              ! Find the colatitude

              if (Iono_Loc(1) == 0.00 .and. Iono_Loc(2) == 0.00) then
                 Iono_Psi = 0.00
              else
                 Iono_Psi = atan2(Iono_Loc(2),Iono_Loc(1))
              end if
              if (Iono_Psi < 0.00) Iono_Psi = Iono_Psi + 2.00*IONO_PI

              ! Find the Longitude

              r = sqrt(Iono_Loc(1)*Iono_Loc(1) + Iono_Loc(2)*Iono_Loc(2))
              Iono_Theta = atan2(r, Iono_Loc(3))

           endif
           if(DoTestThis)write(*,*)'r, iono_theta, iono_psi=',&
                r,Iono_Theta,Iono_Psi
              

           ! Now we search the ionospheric grid for the closest point
           ! We use a binomial search in latitude and longitude to find
           ! the point.  This doesn't really assume a uniform grid in
           ! the ionosphere, but it does assume that at each longitude
           ! the latitude spacing is the same.  The same is true for
           ! the longitude spacing at each latitude.

           ! First set upper and lower limits on the search:

           l_up = IONO_nPsi
           l_lo = 1
           
           done = .false.

           if(DoTestThis)write(*,*)'start search phi'

           do while (.not.done)

              ! Get a central value

              l = (l_up + l_lo)/2

              ! Stop if central value is upper or lower value. This means
              ! that we have narrowed the search down to one value, so
              ! continuing it pointless (i.e. by process of elimination,
              ! this point must be the one we are searching for).

              if ((l == l_up).or.(l == l_lo)) done = .true.

              ! Now we figure out if we hit the value with one of our
              ! three indices that we have picked.  This is redundent,
              ! since in a binomial search, only 2 of our 3 indices
              ! change from iteration to iteration, but checking the
              ! logic to figure out which one changed would be more
              ! work than checking all of them.

              if (Iono_Psi == Psi_temp(1,l_up,hem)) then
                 done = .true.
                 l = l_up
              else
                 if (Iono_Psi == Psi_temp(1,l_lo,hem)) then
                    done = .true.
                    l = l_lo
                 else
                    if (Iono_Psi == Psi_temp(1,l,hem)) then
                       done = .true.
                    else

                       ! If we are unlucky and one of the three indices
                       ! aren't the correct one, we have to do the 
                       ! binomial split.  We just figure out if the point
                       ! is between the l_up and l or l_lo and l.  Which
                       ! ever it is, we narrow the search to that space.

                       if (Iono_Psi < Psi_temp(1,l,hem)) then
                          l_up = l
                       else
                          l_lo = l
                       endif

                    endif
                 endif
              endif

           enddo

           if(DoTestThis)write(*,*)'finish search psi'

           ! Repeat the same logic as above but with the latitudes instead
           ! of the longitudes.

           k_up = IONO_nTheta
           k_lo = 1

           done = .false.

           if(DoTestThis)write(*,*)'start search theta'
           do while (.not.done)

              k = (k_up + k_lo)/2

              if ((k == k_up).or.(k == k_lo)) done = .true.

              if (Iono_Theta == Theta_temp(k_up,1,hem)) then
                 done = .true.
                 k = k_up
              else
                 if (Iono_Theta == Theta_temp(k_lo,1,hem)) then
                    done = .true.
                    k = k_lo
                 else
                    if (Iono_Theta == Theta_temp(k,1,hem)) then
                       done = .true.
                    else

                       if (Iono_Theta < Theta_temp(k,1,hem)) then
                          k_up = k
                       else
                          k_lo = k
                       endif
                    
                    endif
                 endif
              endif

           enddo
           if(DoTestThis)write(*,*)'end search theta'

           ! Now we set the interpolation factors.  This involves knowing
           ! which index to use and the ratios to use between the point
           ! and the point next to it.

           j0 = l
           i0 = k

           if(DoTestThis)write(*,*)'j0,i0=',j0,i0

           if (j0 < IONO_nPsi) then
              lo = 1.0 - (Iono_Psi                - Psi_temp(i0,j0,hem)) / &
                         (Psi_temp(i0,j0+1,hem)   - Psi_temp(i0,j0,hem))
           else
              j0 = IONO_nPsi - 1
              lo = 0.0
           endif

           if(DoTestThis)write(*,*)'lo=',lo

           if (i0 < IONO_nTheta) then
              la = 1.0 - (Iono_Theta              - Theta_temp(i0,j0,hem)) / &
                         (Theta_temp(i0+1,j0,hem) - Theta_temp(i0,j0,hem))
           else
              i0 = IONO_nTheta - 1
              la = 0.0
           endif

           if(DoTestThis)write(*,*)'la=',la

           ! Store the data so we can use it during the next many ionosphere calls

           if(DoTestThis)write(*,*)'storage for n=',n

           if (n == 1) then

              IONOtoMAG_NORTH_hemisphere(i,j) = hem
              IONOtoMAG_NORTH_lat_index(i,j)  = i0
              IONOtoMAG_NORTH_lon_index(i,j)  = j0
              IONOtoMAG_NORTH_lat_factor(i,j) = la
              IONOtoMAG_NORTH_lon_factor(i,j) = lo

           else

              IONOtoMAG_SOUTH_hemisphere(i,j) = hem
              IONOtoMAG_SOUTH_lat_index(i,j)  = i0
              IONOtoMAG_SOUTH_lon_index(i,j)  = j0
              IONOtoMAG_SOUTH_lat_factor(i,j) = la
              IONOtoMAG_SOUTH_lon_factor(i,j) = lo

           endif

        end do
     end do

  enddo

end subroutine GM_map_grid
!==============================================================================
!-------------------------------------------------------------------------
! ionosphere_magBCs
!
!
!
!-------------------------------------------------------------------------
subroutine map_ie_potential_to_gm(PHI,&
          Radius_BC,&
          X, Y, Z,Theta_BC, Psi_BC, nTheta, nPsi,& 
          PHI_BC) 
  use ModMapPotential,ONLY:&
       IONO_nTheta,IONO_nPsi,&
       IONO_NORTH_THETA,IONO_NORTH_PSI,& !Allocated above
       IONO_SOUTH_THETA,IONO_SOUTH_PSI,& !Allocates above
       IONOtoMAG_NORTH_hemisphere, &     !Allocated above
       IONOtoMAG_NORTH_lat_index , &     !Allocated above
       IONOtoMAG_NORTH_lon_index , &     !Allocated above
       IONOtoMAG_NORTH_lat_factor, &     !Allocated above
       IONOtoMAG_NORTH_lon_factor, &     !Allocated above
       IONOtoMAG_SOUTH_hemisphere, &     !Allocated above
       IONOtoMAG_SOUTH_lat_index, &      !Allocated above
       IONOtoMAG_SOUTH_lon_index, &      !Allocated above
       IONOtoMAG_SOUTH_lat_factor, &     !Allocated above
       IONOtoMAG_SOUTH_lon_factor, &     !Allocated above
       IONO_PI,IONO_Radius, IONO_Toler,&
       IONO_NORTH_PHI, IONO_SOUTH_PHI    !Allocated above
  implicit none
  real,intent(in):: Radius_BC
  integer,intent(in):: nTheta, nPsi
  real, dimension(1:IONO_nTheta,1:IONO_nPsi), intent(in) ::    &               
       PHI, X, Y, Z, Theta_BC, Psi_BC
  real, dimension(1:IONO_nTheta,1:IONO_nPsi), intent(out)::    &
       PHI_BC

  real, dimension(1:IONO_nTheta,1:IONO_nPsi,2) ::                            &
       Theta_temp, Psi_temp, phi_temp
  integer :: i, j, i0, j0, hem
  integer :: iMin, jMin
  logical :: north,DoTest,DotestMe
  real :: la, lo
  !\
  ! For the calculated ionospheric potential solution,
  ! map the potential solution out to the magnetospheric
  ! inner boundary at R=Radius_BC.
  !/

  call set_oktest('ionosphere_magBCs',DoTest,DotestMe)

  if (Theta_BC(1,1) < IONO_PI/4.0) then
     north = .true.
  else
     north = .false.
  end if

 if(DoTestMe)write(*,*)'ionosphere_magBCs determine potential'

  if (Radius_BC/IONO_Radius < 1.05) then
     PHI_BC = PHI
  else

     Phi_temp(:,:,1)   = IONO_NORTH_PHI
     Phi_temp(:,:,2)   = IONO_SOUTH_PHI

     iMin = 1000
     jMin = 1000

     do j = 1, nPsi
        do i = 1, nTheta

           if (north) then

              hem = IONOtoMAG_NORTH_hemisphere(i,j)
              i0 = IONOtoMAG_NORTH_lat_index(i,j)
              j0 = IONOtoMAG_NORTH_lon_index(i,j)
              la = IONOtoMAG_NORTH_lat_factor(i,j)
              lo = IONOtoMAG_NORTH_lon_factor(i,j)

           else

              hem = IONOtoMAG_SOUTH_hemisphere(i,j)
              i0 = IONOtoMAG_SOUTH_lat_index(i,j)
              j0 = IONOtoMAG_SOUTH_lon_index(i,j)
              la = IONOtoMAG_SOUTH_lat_factor(i,j)
              lo = IONOtoMAG_SOUTH_lon_factor(i,j)

           endif

           if (i0 < iMin) iMin = i0
           if (j0 < jMin) jMin = j0

	   PHI_BC(i,j) = (    la)*(    lo)*PHI_temp(i0  ,j0  ,hem) + &
                         (1.0-la)*(    lo)*PHI_temp(i0+1,j0  ,hem) + &
                         (    la)*(1.0-lo)*PHI_temp(i0  ,j0+1,hem) + &
                         (1.0-la)*(1.0-lo)*PHI_temp(i0+1,j0+1,hem)

        end do
     end do

  end if

end subroutine map_ie_potential_to_gm

!=============================================================================
subroutine GM_put_from_ie(Buffer_II,iSize,jSize,NameVar)

  use ModPhysics,       ONLY: UnitSi_X, UnitSi_Electric
  use ModMapPotential,ONLY:init_map_potential_arrays,&
       IONO_nTheta,IONO_NORTH_Phi,IONO_SOUTH_Phi 
  implicit none
  character(len=*), parameter :: NameSub='GM_put_from_ie'

  integer, intent(in) :: iSize,jSize
  real, intent(in) :: Buffer_II(iSize,jSize)
  character(len=*), intent(in) :: NameVar

  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest,DoTestMe)
  if(DoTest)write(*,*)NameSub,': NameVar,iSize,jSize=',NameVar,iSize,jSize

  if(IONO_nTheta < 0)then
     if(DoTest)write(*,*)NameSub,': allocating variables'
     call init_inner_bc_arrays(iSize,jSize)
  end if
  if(.not.allocated(IONO_NORTH_Phi))call init_map_potential_arrays  

  if(DoTest)write(*,*)NameSub,': putting potential'
  select case(NameVar)
  case('PotNorth')
     IONO_NORTH_Phi = Buffer_II
  case('PotSouth')
     IONO_SOUTH_Phi = Buffer_II
  case default
     call CON_stop(NameSub//' invalid NameVar='//NameVar)
  end select

  if(DoTest)write(*,*)NameSub,': done'

end subroutine GM_put_from_ie
!==============================================================================
subroutine GM_print_variables(NameSource)

  use ModMain, ONLY: NameThisComp
  use ModNumConst
  use ModImPressure                                  !^CFG IF RCM
  use ModMapPotential,ONLY:IONO_NORTH_PHI,IONO_SOUTH_PHI
  use ModInnerBc
  use ModIoUnit, ONLY: UNITTMP_
  implicit none
  character(len=*), parameter :: NameSub='GM_print_variables'

  character(len=*),intent(in) :: NameSource
  integer            :: nFile=0
  character(len=100) :: NameFile
  character(len=100) :: NameVar
  integer            :: i,j
  !--------------------------------------------------------------------------

  if(NameThisComp/='GM') RETURN

  select case(NameSource)
  case('IM')                        !^CFG IF RCM
     NameVar='j i lon lat p'        !^CFG IF RCM
  case('IE','IE_swmf')
     NameVar='i j theta psi pot'
  case default
     write(*,*)NameSub,': incorrect NameSource=',NameSource
     RETURN
  end select

  nFile=nFile+1
  write(NameFile,'(a,i1,a)')'GM_from_'//NameSource//'_',nFile,'.dat'
  open(UNITTMP_,file=NameFile)
  write(UNITTMP_,'(a)')trim(NameVar)

  select case(NameSource)
  case('IM')                             !^CFG IF RCM BEGIN
     do i=1,iSize
        do j=1,jSize
           write(UNITTMP_,'(2i4,3G14.6)')j,i,RCM_lon(j),RCM_lat(i),RCM_p(i,j)
        end do
     end do                              !^CFG END RCM
  case('IE')
     do j=1,IONO_nPsi
        do i=1,IONO_nTheta 
           write(UNITTMP_,'(2i4,3G14.6)')i,j,&
                IONO_NORTH_Theta(i,j),IONO_NORTH_Psi(i,j),&
                IONO_NORTH_Phi(i,j)
        end do
        do i=1,IONO_nTheta
           write(UNITTMP_,'(2i4,3G14.6)')i+IONO_nTheta,j,&
                IONO_SOUTH_Theta(i,j),IONO_SOUTH_Psi(i,j),&
                IONO_SOUTH_Phi(i,j)
        end do
     end do
  case('IE_swmf')
     do j=1,IONO_nPsi
        do i=1,IONO_nTheta 
           write(UNITTMP_,'(2i4,3G14.6)')i,j,&
                IONO_NORTH_Theta(i,j),IONO_NORTH_Psi(i,j),&
                IONO_NORTH_Phi_BC(i,j)
        end do
        do i=1,IONO_nTheta
           write(UNITTMP_,'(2i4,3G14.6)')i+IONO_nTheta,j,&
                IONO_SOUTH_Theta(i,j),IONO_SOUTH_Psi(i,j),&
                IONO_SOUTH_Phi_BC(i,j)
        end do
     end do
  end select
  close(UNITTMP_)

end subroutine GM_print_variables

!=========================================================================!
!BOP
!MODULE: ionosphere_bc - converts IE  potential to GM  boundary velocity 
!DESCRIPTION:               
!             Ionosphere/magnetosphere coupling                           
!  The data for potential at the GMIE_ grid are transformed  to the       
!  ionosphere velocities by ionosphere_magBCs. Then for any point near the
!  sphere of the radius of IONO_Radius_Mag_Boundary, the subroutine       
!  calc_inner_BC_velocities gives the three cartesian components of the   
!  ionosphere velocities using bilinear interpolation                     
!EOP
!=============================================================================

!BOP
!INTERFACE:
subroutine ionosphere_magBCs(&
     PHI_BC, ETh_BC,EPs_BC,UR_BC, UTh_BC, UPs_BC, Radius_BC,              &
     Theta_BC, Psi_BC, nTheta, nPsi)
  !USES:
  use ModPhysics,ONLY:UnitSI_B,UnitSI_U
  use ModInnerBC,ONLY:&
       IONO_NORTH_THETA,IONO_NORTH_PSI,& 
       IONO_SOUTH_THETA,IONO_SOUTH_PSI,& 
       IONO_PI,IONO_Radius, IONO_Toler
!EOP
  implicit none

  integer, intent(in)  :: nTheta, nPsi
  real, intent(in)     :: Radius_BC

  real, dimension(nTheta,nPsi),intent(in)     :: &
        PHI_BC,Theta_BC, Psi_BC
  real, dimension(nTheta,nPsi), intent(out)   :: &
       ETh_BC, EPs_BC ,UR_BC, UTh_BC, UPs_BC                      
          
  real, dimension(nTheta) :: dTheta
  real, dimension(nPsi)   :: dPsi

  logical :: north
  integer :: i, j
  real :: dTheta_l, dPsi_l, dTheta_l2, dPsi_l2
  real :: dd, dd2 
  real :: cosTheta, sinTheta, cosPhi, sinPhi
  real :: xx, yy, zz
  real :: ER, Ex, Ey, Ez
  real :: Ux, Uy, Uz
  real :: BB, Bx, By,Bz
 

  real, dimension(3) :: Mag_Loc, Mag_B


  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------
  call set_oktest('ionosphere_magBCs',DoTest,DotestMe)

  if(DoTestMe)write(*,*)'ionosphere_magBCs starting with iono_pi=',IONO_PI

  if (Theta_BC(1,1) < IONO_PI/4.0) then
     north = .true.
  else
     north = .false.
  end if

  dTheta_l=(Theta_BC(nTheta,1)-Theta_BC(1,1))/real(nTheta-1)
  dPsi_l=(Psi_BC(1,nPsi)-Psi_BC(1,1))/real(nPsi-1)
  dTheta_l2=dTheta_l*dTheta_l
  dPsi_l2=dPsi_l*dPsi_l
  dd=dPsi_l
  dd2=dPsi_l2

  if(DoTestMe)write(*,*)'dTheta_l,dPsi_l=',dTheta_l,dPsi_l

  ! Determine the potential at the magnetospheric inner boundary.



  ! Compute the electric field at the magnetospheric inner boundary.

  if(DoTestMe)write(*,*)'ionosphere_magBCs compute electric field'

  do j = 1, nPsi
     if (j > 1 .and. j < nPsi ) then 
        do i = 2, nTheta-1
           sinTheta = sin(Theta_BC(i,j))
           ETh_BC(i,j) = -(PHI_BC(i+1,j)-PHI_BC(i-1,j))/ &
                         (2.00*dd*Radius_BC)
           EPs_BC(i,j) = -(PHI_BC(i,j+1)-PHI_BC(i,j-1))/ &
                         (2.00*dd*Radius_BC*sinTheta)
        end do
        ETh_BC(1,j) = -(PHI_BC(2,j)-PHI_BC(1,j))/ &
                      (dd*Radius_BC)
        EPs_BC(1,j) = EPs_BC(2,j)
        ETh_BC(nTheta,j) = -(PHI_BC(nTheta,j)-PHI_BC(nTheta-1,j))/ &
                           (dd*Radius_BC)
        EPs_BC(nTheta,j) = EPs_BC(nTheta-1,j)
     else if (j == 1) then
        do i = 2, nTheta-1
           sinTheta = sin(Theta_BC(i,j))
           ETh_BC(i,j) = -(PHI_BC(i+1,j)-PHI_BC(i-1,j))/ &
                         (2.00*dd*Radius_BC)
           EPs_BC(i,j) = -(PHI_BC(i,j+1)-PHI_BC(i,nPsi-1))/ &
                         (2.00*dd*Radius_BC*sinTheta)
        end do
        ETh_BC(1,j) = -(PHI_BC(2,j)-PHI_BC(1,j))/ &
                      (dd*Radius_BC)
        EPs_BC(1,j) = EPs_BC(2,j)
        ETh_BC(nTheta,j) = -(PHI_BC(nTheta,j)-PHI_BC(nTheta-1,j))/ &
                           (dd*Radius_BC)
        EPs_BC(nTheta,j) = EPs_BC(nTheta-1,j)
     else
        do i = 2, nTheta-1
           sinTheta = sin(Theta_BC(i,j))
           ETh_BC(i,j) = -(PHI_BC(i+1,j)-PHI_BC(i-1,j))/ &
                         (2.00*dd*Radius_BC)
           EPs_BC(i,j) = -(PHI_BC(i,2)-PHI_BC(i,j-1))/ &
                         (2.00*dd*Radius_BC*sinTheta)
        end do
        ETh_BC(1,j) = -(PHI_BC(2,j)-PHI_BC(1,j))/ &
                      (dd*Radius_BC)
        EPs_BC(1,j) = EPs_BC(2,j)
        ETh_BC(nTheta,j) = -(PHI_BC(nTheta,j)-PHI_BC(nTheta-1,j))/ &
                           (dd*Radius_BC)
        EPs_BC(nTheta,j) = EPs_BC(nTheta-1,j)
     end if
  end do

  ! Compute the convection velocities at the magnetospheric inner boundary.

  if(DoTestMe)write(*,*)'ionosphere_magBCs compute velocities'

  do j = 1, nPsi
     do i = 1, nTheta
        cosTheta = cos(Theta_BC(i,j))
        sinTheta = sin(Theta_BC(i,j))
        cosPhi = cos(Psi_BC(i,j))
        sinPhi = sin(Psi_BC(i,j))

        if (north .and. i == nTheta) then
           ER = 0.00
        else if (.not.north .and. i == 1) then
           ER = 0.00
        else
           ER = -0.50*(sinTheta/(cosTheta+IONO_Toler**2))*ETh_BC(i,j)
        end if
      
        Ex = ER*sinTheta*cosPhi + ETh_BC(i,j)*cosTheta*cosPhi - &
             EPs_BC(i,j)*sinPhi
        Ey = ER*sinTheta*sinPhi + ETh_BC(i,j)*cosTheta*sinPhi + &
             EPs_BC(i,j)*cosPhi
        Ez = ER*cosTheta - ETh_BC(i,j)*sinTheta
        
        Mag_Loc(1) = Radius_BC * sin(Theta_BC(i,j))* &
                                 cos(Psi_BC(i,j)) / IONO_Radius
        Mag_Loc(2) = Radius_BC * sin(Theta_BC(i,j))* &
                                 sin(Psi_BC(i,j)) / IONO_Radius
        Mag_Loc(3) = Radius_BC * cos(Theta_BC(i,j)) / IONO_Radius

        call get_b0(Mag_Loc(1),Mag_Loc(2),Mag_Loc(3),Mag_B)

        Mag_B = Mag_B * UnitSI_B

        BB = sum(Mag_B**2)

        bx = Mag_B(1)
        by = Mag_B(2)
        bz = Mag_B(3)
        
        Ux = (Ey*bz - Ez*by)/BB
        Uy = (Ez*bx - Ex*bz)/BB
        Uz = (Ex*by - Ey*bx)/BB
        
        xx = sinTheta*cosPhi
        yy = sinTheta*sinPhi
        zz = cosTheta

        UR_BC(i,j) = Ux*xx + Uy*yy + Uz*zz
        UTh_BC(i,j) = ((Ux*xx + Uy*yy)*zz - &
                       Uz*(xx**2+yy**2)) / &
                      sqrt(xx**2+yy**2+IONO_TOLER**2)
        UPs_BC(i,j) = (Uy*xx - Ux*yy)*sinTheta / &
                      (xx**2+yy**2+IONO_TOLER**2)

        UR_BC(i,j) = UR_BC(i,j) / UnitSI_U
        UTh_BC(i,j) = UTh_BC(i,j) / UnitSI_U
        UPs_BC(i,j) = UPs_BC(i,j) / UnitSI_U

     end do  
  end do

end subroutine ionosphere_magBCs

subroutine calc_inner_BC_velocities(iter,time_now,XFace,YFace,ZFace,&
     BxOutside,ByOutside,BzOutside,B0faceX,B0FaceY,B0FaceZ,&
     VxFace,VyFace,VzFace)
!USES:
  use ModNumConst
  use ModInnerBC
!DESCRIPTION:
  !\
  ! This subroutine should calculate -x -y and -z components of the "ionosphere" velocities
  ! as a function of the cartesian ponit coordinates x,y,z, as well as the 
  ! induction and static magnetic field
  ! Iteration number and time can be also used.
  !/
!EOP
  implicit none

  integer,intent(in) :: iter
  real, intent(in) :: time_now
  real,intent(in) :: XFace,YFace,ZFace,&
       BxOutside,ByOutside,BzOutside,B0faceX,B0FaceY,B0FaceZ
  real,intent(out) :: VxFace,VyFace,VzFace
  real:: VrFace,VThetaFace,VphiFace
  integer::i0,j0
  real:: RFace, cosTheta, sinTheta, cosPhi, sinPhi
  real :: ThetaFace, PhiFace, dTheta, dPsi, dd

  real :: xC1, yC1, xC2, yC2, xC3, yC3, xC4, yC4, &
       fC1, fC2, fC3, fC4, &
       a_bilin, b_bilin, c_bilin, d_bilin

  RFace=sqrt(XFace*XFace+YFace*YFace+ZFace*ZFace)
  cosTheta = ZFace/RFace
  sinTheta = sqrt(XFace**2+YFace**2)/RFace
  cosPhi = XFace/sqrt(XFace**2+YFace**2+cTolerance**2)
  sinPhi = YFace/sqrt(XFace**2+YFace**2+cTolerance**2)

  if (ZFace > cZero) then
     !\
     ! Northern Hemisphere Ionosphere BC.
     !/
     dTheta = (IONO_NORTH_Theta(IONO_nTheta,1)- &
          IONO_NORTH_Theta(1,1))/real(IONO_nTheta-1)
     dPsi = (IONO_NORTH_Psi(1,IONO_nPsi)- &
          IONO_NORTH_Psi(1,1))/real(IONO_nPsi-1)
     dd = dPsi

     if (XFace == cZero .and. YFace == cZero) then
        PhiFace = cZero
     else
        PhiFace = atan2(sinPhi, cosPhi)
     end if
     if (PhiFace < cZero) PhiFace = PhiFace + cTwoPi
     ThetaFace = atan2(sinTheta, cosTheta)

     i0 = ThetaFace/dd + 1
     j0 = PhiFace/dd + 1

     xC1 = IONO_NORTH_Theta(i0,j0)
     yC1 = IONO_NORTH_Psi(i0,j0)
     xC2 = IONO_NORTH_Theta(i0,j0+1)
     yC2 = IONO_NORTH_Psi(i0,j0+1)
     xC3 = IONO_NORTH_Theta(i0+1,j0+1)
     yC3 = IONO_NORTH_Psi(i0+1,j0+1)
     xC4 = IONO_NORTH_Theta(i0+1,j0)
     yC4 = IONO_NORTH_Psi(i0+1,j0)

     fC1 = IONO_NORTH_UR_BC(i0,j0)
     fC2 = IONO_NORTH_UR_BC(i0,j0+1)
     fC3 = IONO_NORTH_UR_BC(i0+1,j0+1)
     fC4 = IONO_NORTH_UR_BC(i0+1,j0)


     a_bilin = fC1
     b_bilin = (fC4-fC1)/dd
     c_bilin = (fC2-fC1)/dd
     d_bilin = (fC3+fC1-fC4-fC2)/(dd*dd)
     VrFace = a_bilin + &
          b_bilin*(ThetaFace-xC1) + &
          c_bilin*(PhiFace-yC1) + &
          d_bilin*(ThetaFace-xC1)*(PhiFace-yC1)

     fC1 = IONO_NORTH_UTh_BC(i0,j0)
     fC2 = IONO_NORTH_UTh_BC(i0,j0+1)
     fC3 = IONO_NORTH_UTh_BC(i0+1,j0+1)
     fC4 = IONO_NORTH_UTh_BC(i0+1,j0)
     a_bilin = fC1
     b_bilin = (fC4-fC1)/dd
     c_bilin = (fC2-fC1)/dd
     d_bilin = (fC3+fC1-fC4-fC2)/(dd*dd)
     VthetaFace = a_bilin + &
          b_bilin*(ThetaFace-xC1) + &
          c_bilin*(PhiFace-yC1) + &
          d_bilin*(ThetaFace-xC1)*(PhiFace-yC1)

     fC1 = IONO_NORTH_UPs_BC(i0,j0)
     fC2 = IONO_NORTH_UPs_BC(i0,j0+1)
     fC3 = IONO_NORTH_UPs_BC(i0+1,j0+1)
     fC4 = IONO_NORTH_UPs_BC(i0+1,j0)
     a_bilin = fC1
     b_bilin = (fC4-fC1)/dd
     c_bilin = (fC2-fC1)/dd
     d_bilin = (fC3+fC1-fC4-fC2)/(dd*dd)
     VphiFace = a_bilin + &
          b_bilin*(ThetaFace-xC1) + &
          c_bilin*(PhiFace-yC1) + &
          d_bilin*(ThetaFace-xC1)*(PhiFace-yC1)

  else
     !\
     ! Southern Hemisphere Ionosphere BC.
     !/
     dTheta = (IONO_SOUTH_Theta(IONO_nTheta,1)- &
          IONO_SOUTH_Theta(1,1))/real(IONO_nTheta-1)
     dPsi = (IONO_SOUTH_Psi(1,IONO_nPsi)- &
          IONO_SOUTH_Psi(1,1))/real(IONO_nPsi-1)
     dd = dPsi

     if (XFace == cZero .and. YFace == cZero) then
        PhiFace = cZero
     else
        PhiFace = atan2(sinPhi, cosPhi)
     end if
     if (PhiFace < cZero) PhiFace = PhiFace + cTwoPi
     ThetaFace = cPi - &
          atan2(sinTheta, -cosTheta)

     i0 = (ThetaFace-cHalfPi)/dd + 1
     j0 = PhiFace/dd + 1

     xC1 = IONO_SOUTH_Theta(i0,j0)
     yC1 = IONO_SOUTH_Psi(i0,j0)
     xC2 = IONO_SOUTH_Theta(i0,j0+1)
     yC2 = IONO_SOUTH_Psi(i0,j0+1)
     xC3 = IONO_SOUTH_Theta(i0+1,j0+1)
     yC3 = IONO_SOUTH_Psi(i0+1,j0+1)
     xC4 = IONO_SOUTH_Theta(i0+1,j0)
     yC4 = IONO_SOUTH_Psi(i0+1,j0)

     fC1 = IONO_SOUTH_UR_BC(i0,j0)
     fC2 = IONO_SOUTH_UR_BC(i0,j0+1)
     fC3 = IONO_SOUTH_UR_BC(i0+1,j0+1)
     fC4 = IONO_SOUTH_UR_BC(i0+1,j0)
     a_bilin = fC1
     b_bilin = (fC4-fC1)/dd
     c_bilin = (fC2-fC1)/dd
     d_bilin = (fC3+fC1-fC4-fC2)/(dd*dd)
     VrFace = a_bilin + &
          b_bilin*(ThetaFace-xC1) + &
          c_bilin*(PhiFace-yC1) + &
          d_bilin*(ThetaFace-xC1)*(PhiFace-yC1)

     fC1 = IONO_SOUTH_UTh_BC(i0,j0)
     fC2 = IONO_SOUTH_UTh_BC(i0,j0+1)
     fC3 = IONO_SOUTH_UTh_BC(i0+1,j0+1)
     fC4 = IONO_SOUTH_UTh_BC(i0+1,j0)
     a_bilin = fC1
     b_bilin = (fC4-fC1)/dd
     c_bilin = (fC2-fC1)/dd
     d_bilin = (fC3+fC1-fC4-fC2)/(dd*dd)
     VthetaFace = a_bilin + &
          b_bilin*(ThetaFace-xC1) + &
          c_bilin*(PhiFace-yC1) + &
          d_bilin*(ThetaFace-xC1)*(PhiFace-yC1)

     fC1 = IONO_SOUTH_UPs_BC(i0,j0)
     fC2 = IONO_SOUTH_UPs_BC(i0,j0+1)
     fC3 = IONO_SOUTH_UPs_BC(i0+1,j0+1)
     fC4 = IONO_SOUTH_UPs_BC(i0+1,j0)
     a_bilin = fC1
     b_bilin = (fC4-fC1)/dd
     c_bilin = (fC2-fC1)/dd
     d_bilin = (fC3+fC1-fC4-fC2)/(dd*dd)
     VphiFace = a_bilin + &
          b_bilin*(ThetaFace-xC1) + &
          c_bilin*(PhiFace-yC1) + &
          d_bilin*(ThetaFace-xC1)*(PhiFace-yC1)
  end if ! South or North

 
  VrFace = cZero
  VxFace  = VrFace*cosPhi*sinTheta + VthetaFace*cosTheta*cosPhi - VphiFace*sinPhi 
  VyFace  = VrFace*sinPhi*sinTheta + VthetaFace*cosTheta*sinPhi + VphiFace*cosPhi 
  VzFace  = VrFace*cosTheta - VthetaFace*sinTheta

end subroutine calc_inner_BC_velocities

!^CFG COPYRIGHT UM
!
!-------------------------------------------------------------------------
Module ModMappingParam
  use ModNumConst

  real,parameter::IONO_PI=cPi
  real, parameter :: IONO_TOLER = 5.0e-05
  real :: IONO_Radius_Mag_Boundary, &
       IONO_Radius, IONO_Height, rMap, DipoleSign

end Module ModMappingParam

!============================================================================

subroutine set_mapping_param

  use ModPhysics,  ONLY: unitSI_x,Rbody
  use CON_physics, ONLY: get_planet
  use ModMappingParam
  implicit none
  real::DipoleStrength

  call get_planet(&
       RadiusPlanetOut      = IONO_Radius,&
       IonosphereHeightOut  = IONO_Height,&
       DipoleStrengthOut    = DipoleStrength)

  ! final radius is fixed, store it
  rMap=(IONO_Radius+ IONO_Height)/UnitSI_x
  IONO_Radius_Mag_Boundary = Rbody*UnitSI_x
  DipoleSign=sign(cOne,DipoleStrength)

end subroutine set_mapping_param

!============================================================================

subroutine init_inner_bc_arrays(iSize,jSize)

  use CON_coupler
  use ModInnerBc,ONLY:&
       IONO_nTheta, IONO_nPsi,&            
       IONO_NORTH_Theta, IONO_NORTH_Psi, & 
       IONO_SOUTH_Theta, IONO_SOUTH_Psi, &   
       IONO_NORTH_PHI_BC,IONO_SOUTH_PHI_BC,         & 
       IONO_NORTH_ETh_BC,IONO_NORTH_EPs_BC,& 
       IONO_SOUTH_ETh_BC,IONO_SOUTH_EPs_BC,&                  
       IONO_NORTH_UR_BC, &
       IONO_NORTH_UTh_BC, &
       IONO_NORTH_UPs_BC, &
       IONO_SOUTH_UR_BC, &
       IONO_SOUTH_UTh_BC, &
       IONO_SOUTH_UPs_BC

  implicit none

  integer, intent(in) :: iSize,jSize
  character(len=*), parameter :: NameSub='GM_allocate_iono_arrays'
  integer :: nCells_D(2),iError, i
  ! Allocate and calculate coordinates
  
  nCells_D=ncells_decomposition_d(IE_)
  if(  iSize /= nCells_D(1)+1 .or. &
       jSize /= nCells_D(2)+1 ) then
     
     write(*,*)NameSub//' grid sizes do not agree iSize,jSize,nCells=',&
          iSize,jSize,nCells_D(1:2)
     call CON_stop(NameSub//' SWMF_ERROR')
  end if
  
  IONO_nTheta = iSize
  IONO_nPsi   = jSize
  
  allocate(&
       IONO_NORTH_Theta(iSize,jSize), IONO_NORTH_Psi(iSize,jSize), &
       IONO_NORTH_Phi_BC(iSize,jSize),                             &
       IONO_SOUTH_Theta(iSize,jSize), IONO_SOUTH_Psi(iSize,jSize), &
       IONO_SOUTH_Phi_BC(iSize,jSize),                             &
       IONO_NORTH_ETh_BC(IONO_nTheta,IONO_nPsi),                   &
       IONO_NORTH_EPs_BC(IONO_nTheta,IONO_nPsi),                   &
       IONO_NORTH_UR_BC(IONO_nTheta,IONO_nPsi),                    &
       IONO_NORTH_UTh_BC(IONO_nTheta,IONO_nPsi),                   &
       IONO_NORTH_UPs_BC(IONO_nTheta,IONO_nPsi),                   &
       IONO_SOUTH_ETh_BC(IONO_nTheta,IONO_nPsi),                   &
       IONO_SOUTH_EPs_BC(IONO_nTheta,IONO_nPsi),                   &
       IONO_SOUTH_UR_BC(IONO_nTheta,IONO_nPsi),                    &
       IONO_SOUTH_UTh_BC(IONO_nTheta,IONO_nPsi),                   &
       IONO_SOUTH_UPs_BC(IONO_nTheta,IONO_nPsi),                   &
       STAT = iError)
  call check_allocate(iError,NameSub//' IONO_NORTH/SOUTH_ThetaPsiPhi')
  IONO_NORTH_Phi_BC=cZero
  IONO_SOUTH_Phi_BC=cZero
  IONO_NORTH_ETh_BC=cZero
  IONO_NORTH_EPs_BC=cZero
  IONO_NORTH_UR_BC=cZero
  IONO_NORTH_UTh_BC=cZero
  IONO_NORTH_UPs_BC=cZero
  IONO_SOUTH_ETh_BC=cZero
  IONO_SOUTH_EPs_BC=cZero
  IONO_SOUTH_UR_BC=cZero
  IONO_SOUTH_UTh_BC=cZero
  IONO_SOUTH_UPs_BC=cZero
  
  do i=1,iSize
     IONO_NORTH_Theta(i,:) = Grid_C(IE_) % Coord1_I(i)
     IONO_SOUTH_Theta(i,:) = Grid_C(IE_) % Coord1_I(i+iSize-1)
     IONO_NORTH_Psi(i,:)   = Grid_C(IE_) % Coord2_I
     IONO_SOUTH_Psi(i,:)   = Grid_C(IE_) % Coord2_I
  end do

end subroutine init_inner_bc_arrays
!^CFG COPYRIGHT UM
!^CMP FILE IE
!==================================================================!
!BOP
!ROUTINE: GM_get_for_ie - gets FACurrents for coupling via SWMF
!INTERFACE:
subroutine GM_get_for_ie_swmf(&
     nPartial,iGet,Get,W,State_V,nVar)
!USES:
  use ModMain, ONLY    : TypeCoordSystem, Time_Simulation
  use CON_physics, ONLY: get_planet_field, transform_matrix
  use CON_coupler,ONLY : IndexPtrType,WeightPtrType,cOne
  use ModGeometry,ONLY : x_BLK,y_BLK,z_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModPhysics,ONLY  : UnitSI_J,UnitSI_X, RBody
  use ModAdvance,ONLY  : State_VGB, Bx_, By_, Bz_
  use ModMappingParam,ONLY:rMap,DipoleSign
  implicit none
!INPUT ARGUMENTS:  
  integer,intent(in)::nPartial,iGet,nVar
  type(IndexPtrType),intent(in)::Get
  type(WeightPtrType),intent(in)::W
!OUTPUT ARGUMENTS:
  real,dimension(nVar),intent(out)::State_V
!REVISION HISTORY:
!14AUG03 - Gabor Toth <gtoth@umich,edu> and Aaron Ridley - 
!          initial prototype/code for coupling via MPI
!21Aug03 - I.Sokolov <igorsok@umich.edu> - major revision for
!          coupling via SWMF/prolog
!EOP
  real,dimension(3)::XyzGm_D,XyzIe_D,XyzSmg_D
  real,dimension(3)::FieldGm_D,FieldIe_D
  real::dir
  integer::i,j,k,iBLK
  
  i=Get%iCB_II(1,iGet)
  j=Get%iCB_II(2,iGet)
  k=Get%iCB_II(3,iGet)
  iBLK=Get%iCB_II(4,iGet)
  
  State_V(1) = 0.50* &              
       ((State_VGB(Bz_,i,j+1,k,iBLK) - &
       State_VGB(Bz_,i,j-1,k,iBLK))/dy_BLK(iBLK) - &
       (State_VGB(By_,i,j,k+1,iBLK) - &
       State_VGB(By_,i,j,k-1,iBLK))/dz_BLK(iBLK))
  
  State_V(2) =  0.50* &
       ((State_VGB(Bx_,i,j,k+1,iBLK) - &
       State_VGB(Bx_,i,j,k-1,iBLK))/dz_BLK(iBLK) - &
       (State_VGB(Bz_,i+1,j,k,iBLK) - &
       State_VGB(Bz_,i-1,j,k,iBLK))/dx_BLK(iBLK))
  
  State_V(3) = 0.50* &
       ((State_VGB(By_,i+1,j,k,iBLK) - &
       State_VGB(By_,i-1,j,k,iBLK))/dx_BLK(iBLK) - &
       (State_VGB(Bx_,i,j+1,k,iBLK) - &
       State_VGB(Bx_,i,j-1,k,iBLK))/dy_BLK(iBLK))
  
  !Transformation to the SI system:
  State_V(1:3)=State_V(1:3)*UnitSI_J
  if(nVar==3)return
  
  XyzGm_D(1)=x_BLK(i,j,k,iBLK)
  XyzGm_D(2)=y_BLK(i,j,k,iBLK)
  XyzGm_D(3)=z_BLK(i,j,k,iBLK)
  
  ! Transform input coordinates to SMG system (used by IE)
  XyzSmg_D = matmul(transform_matrix(&
       Time_Simulation,TypeCoordSystem,'SMG'),&
       XyzGm_D)
  ! Check if the initial point is on the northern hemisphere
  If(XyzSmg_D(3)>=0)then
     call Get_Mapping_Point(XyzGm_D, 1.0, XyzIe_D)
     dir=DipoleSign
  else
     call Get_Mapping_Point(XyzGm_D, -1.0, XyzIe_D)
     dir=-DipoleSign
  end If
  call Get_MagneticField_Orient(XyzGm_D,State_V(4:6))
  call get_planet_field(Time_Simulation,XyzGm_D,&
       TypeCoordSystem//' NORM',FieldGm_D)
  call get_planet_field(Time_Simulation,XyzIe_D,&
       'SMG NORM', FieldIe_D)
  ! Get the ratio between the magnitude of the ionosphere B 
  ! and the magnetosphere B
  State_V(7)=dir*sqrt(dot_product(FieldIe_D,FieldIe_D)/&
       dot_product(FieldGm_D,FieldGm_D) )
  
  ! Get the absolute value of the cos of the angle
  !  between the magnetic field and radial direction
  State_V(8) = &
       abs(sum(FieldIe_D*XyzIe_D))/&
       (rMap*sqrt(dot_product(FieldIe_D,FieldIe_D)))
  if(nVar==11)State_V(9:11)=XyzIe_D
end subroutine GM_get_for_ie_swmf
!=====================================================================
!^CFG COPYRIGHT UM
!^CMP FILE IE
!BOP
!MODULE: GM_put_from_ie_swmf - transform IE potential to GM boundary velocities
!IROUTINE:GMIE_set_grid - set 2D spherical grid in GM_ to get IE_ potential
!INTERFACE:
subroutine GMIE_set_grid
  !USES:
  use CON_coupler
  use ModInnerBC
  use ModMain,ONLY:NameThisComp
  use ModProcMH,ONLY:iProc
!EOP
  implicit none
  integer,dimension(2)::nCell_D
  integer,parameter::GmIeGrid_=MaxComp+1
  integer::iError
  real::rBoundary=cOne
  real,allocatable,dimension(:)::Colat_I,Longit_I
  logical::DoTest,DoTestMe

  call CON_set_do_test('test_grids',DoTest,DoTestMe)

  if(.not.(done_dd_init(GM_).and.done_dd_init(IE_)))&
       call CON_stop('GMIE_grid should be intialized after GM_ and IE_')

  nCell_D=ncells_decomposition_D(IE_)
  if(NameThisComp=='GM') then
     call set_mapping_param
     call init_inner_bc_arrays(nCell_D(1)+1,nCell_D(2)+1)
  end if

  allocate(Colat_I(1:2*nCell_D(1)+1),&
       Longit_I(1:nCell_D(2)+1),stat=iError)
  call check_allocate(iError,'GMIE_set_grid')
  Colat_I=cZero;Longit_I=0



  call init_decomposition(GmIeGrid_,GM_,2)
  if(NameThisComp=='GM')then
     Colat_I(1:IONO_nTheta)=IONO_NORTH_Theta(:,1)
     Colat_I(IONO_nTheta:2*IONO_nTheta-1)=IONO_SOUTH_Theta(:,1)
     Longit_I=IONO_NORTH_Psi(1,:)
     rBoundary=IONO_Radius_Mag_Boundary/IONO_Radius
  end if
  call set_coord_system(GmIeGrid_, &
       Grid_C(GM_)%TypeCoord,&
       Colat_I,&
       Longit_I,&
       (/rBoundary/))
  if(NameThisComp=='GM'.and.iProc==0)call get_root_decomposition(&
       GridID_=GmIeGrid_,&
       iRootMapDim_D=(/2,1/),&
       XyzMin_D=(/cOne,cOne/),&
       XyzMax_D=(/real(2*IONO_nTheta-1),real(IONO_nPsi)/),&
       nCells_D=nCell_D,&
       iBlock_I=(/1,2/))
  call bcast_decomposition(GmIeGrid_)
  if(DoTest)call test_global_message_pass(GmIeGrid_)
end subroutine GMIE_set_grid

!IROUTINE: GM_put_from_ie_swmf - put the interpolated IE potential to GmIe grid
!INTERFACE:
subroutine GM_put_from_ie_swmf(nPartial,iPut,Put,W,DoAdd,State_V,nVar)
  !USES:
  use CON_coupler,ONLY:IndexPtrType, WeightPtrType
  use ModInnerBC,ONLY:IONO_NORTH_Phi_BC,IONO_SOUTH_Phi_BC
  implicit none
!INPUT ARGUMENTS:
  integer,intent(in)::nPartial,iPut,nVar
  type(IndexPtrType),intent(in)::Put
  type(WeightPtrType),intent(in)::W
  logical,intent(in)::DoAdd
  real,dimension(nVar),intent(in)::State_V
  character(LEN=*),parameter::SubName='GM_put_from_ie_swmf'
!REVISION HISTORY:
!03SEPT03     I.Sokolov <igorsok@umich.edu> - prototype/code/prolog
!EOP

  if(DoAdd)call CON_stop(SubName//' DoAdd can not be .true.')
  select case(Put%iCB_II(3,iPut))
  case(1)
     IONO_NORTH_Phi_BC(Put%iCB_II(1,iPut),Put%iCB_II(2,iPut))=State_V(1)
  case(2)
     IONO_SOUTH_Phi_BC(Put%iCB_II(1,iPut),Put%iCB_II(2,iPut))=State_V(1)
  case default
     call CON_stop(SubName//'  - wrong block number:',Put%iCB_II(3,iPut))
  end select
end subroutine GM_put_from_ie_swmf
!BOP
!IROUTINE:transform_phi_bc_to_u_bc - transform  GmIe potential  to GM boundary  velocities
!INTERFACE: 
subroutine transform_phi_bc_to_u_bc
  !USES:
  use ModProcMH
  use ModInnerBC
  use CON_coupler
!EOP
  implicit none

  integer,parameter::GmIeGrid_=MaxComp+1
  integer::nCell_D(2),nSize,iError

  if(iProc<0)return
  nCell_D=ncells_decomposition_d(GmIeGrid_)
  nSize=product(nCell_D+1)
  if(nProc>1)then
     call MPI_BCAST(IONO_NORTH_Phi_BC(1,1),nSize,MPI_REAL,0,iComm,iError)
     call MPI_BCAST(IONO_SOUTH_Phi_BC(1,1),nSize,MPI_REAL,0,iComm,iError)
  end if
  call ionosphere_magBCs(&
     IONO_NORTH_PHI_BC,&
     IONO_NORTH_ETh_BC,&
     IONO_NORTH_EPs_BC,&
     IONO_NORTH_UR_BC, &
     IONO_NORTH_UTh_BC,&
     IONO_NORTH_UPs_BC,&
     IONO_Radius_Mag_Boundary,  &
     IONO_NORTH_Theta,      &
     IONO_NORTH_Psi,        &
     IONO_nTheta,              &
     IONO_nPsi)
  call ionosphere_magBCs(&
     IONO_SOUTH_PHI_BC,&
     IONO_SOUTH_ETh_BC,&
     IONO_SOUTH_EPs_BC,&
     IONO_SOUTH_UR_BC, &
     IONO_SOUTH_UTh_BC,&
     IONO_SOUTH_UPs_BC,&
     IONO_Radius_Mag_Boundary,  &
     IONO_SOUTH_Theta,      &
     IONO_SOUTH_Psi,        &
     IONO_nTheta,              &
     IONO_nPsi)
end subroutine transform_phi_bc_to_u_bc
!^CFG COPYRIGHT UM


subroutine Get_MagneticField_Orient(qx,qb)

  implicit none

  ! Obtain normalized bb field at true location qx and put it into qb

  real, intent(in) :: qx(3)
  real, intent(out):: qb(3)

  ! Get B0

  call get_b0(qx(1),qx(2),qx(3),qb)

  ! Normalize
  qb=qb/sqrt(sum(qb**2))

end subroutine Get_MagneticField_Orient

!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------

subroutine Get_Mapping_Point(XyzIn_D, DirectionIn, XyzOut_D)

  ! Given the XyzIn_D coordinates in the magnetosphere system, 
  ! find the mapping point along the magnetic field lines for
  ! the hemisphere signaled by HemisphereIn (+1.0 for north, -1.0 for south).
  ! The coordinates of the mapping point will be in XyzOut_D in the
  ! ionosphere coordinate system. If HemisphereIn differs from the
  ! hemisphere where XyzIn_D is, then XyzOut_D is set to 0.0.

  use ModMain, ONLY: TypeCoordSystem, Time_Simulation
  use ModMappingParam, ONLY: rMap
  use CON_physics, ONLY: map_planet_field, transform_matrix

  implicit none

  real, intent(in)  :: XyzIn_D(3)
  real, intent(in)  :: DirectionIn
  real, intent(out) :: XyzOut_D(3)

  real :: XyzSmg_D(3) ! Coordinates is solar magnetic system

  integer    :: iHemisphere
  !----------------------------------------------------------------------------


  ! Transform input coordinates to SMG system (used by IE)
  XyzSmg_D = matmul(transform_matrix(Time_Simulation,TypeCoordSystem,'SMG'),&
       XyzIn_D)

  ! Check if the initial point is on the required hemisphere
  if( XyzSmg_D(3)*DirectionIn < 0.0) then
     ! We are on the wrong hemisphere, return 0. location
     XyzOut_D = 0.0
     RETURN
  end if

  ! Get the mapping point
  call map_planet_field(Time_Simulation,XyzSmg_D,'SMG NORM',rMap, &
       XyzOut_D,iHemisphere)

end subroutine Get_Mapping_Point
!==============================================================================

!-------------------------------------------------------------------------
! Get_Mapping_Information
!
subroutine GM_get_mapping_param_for_ie(rCurrentsGm,rIonosphere)
  use ModPhysics,ONLY:rCurrents
  use ModMappingParam,ONLY:rMap
  implicit none
  real,intent(out)::rCurrentsGM,rIonosphere
  rCurrentsGm=rCurrents;rIonosphere=rMap
end subroutine GM_get_mapping_param_for_ie

!=============================================================================

real function logvar_ionosphere(NameLogvar)
  use ModProcMH,  ONLY: iProc
  use ModIO,      ONLY: write_myname
  use ModInnerBc, ONLY: IONO_NORTH_PHI_BC,IONO_SOUTH_PHI_BC
  implicit none
  character (len=*), intent(in) :: NameLogvar
  integer :: nWarn = 0 ! warn multiple times to catch multiple log variables
  !---------------------------------------------------------------------------
  select case(NameLogvar)
  case('cpcpn','cpcp_n','cpcp_north','cpcpnorth')
     logvar_ionosphere = maxval(IONO_NORTH_PHI_BC)-minval(IONO_NORTH_PHI_BC)
  case('cpcps','cpcp_s','cpcp_south','cpcpsouth')
     logvar_ionosphere = maxval(IONO_SOUTH_PHI_BC)-minval(IONO_SOUTH_PHI_BC)
  case default
     if(nWarn < 2 .and. iProc==0)then
        call write_myname;
        write(*,'(a)')'WARNING in logvar_ionosphere: unknown NameLogvar='//&
             trim(NameLogvar)//', returning -777.77'
        nWarn = nWarn + 1
     end if
     logvar_ionosphere = -777.77
  end select

end function logvar_ionosphere
