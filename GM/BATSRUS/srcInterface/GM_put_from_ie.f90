!^CFG COPYRIGHT UM
!^CMP FILE IE
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

!BUFFER->PHI
subroutine GM_put_from_ie(Buffer_II,iSize,jSize,NameVar)

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
