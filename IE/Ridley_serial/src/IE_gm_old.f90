module ModIonosphere
  ...
  real, parameter    :: Max_FAC_Distance = cPi / IONO_nTheta
  integer :: IONO_NORTH_nMagBndPts=-1, IONO_SOUTH_nMagBndPts=-1

  integer, dimension(1:IONO_nTheta,1:IONO_nPsi) ::             &
       IONOtoMAG_NORTH_hemisphere, IONOtoMAG_SOUTH_hemisphere, &
       IONOtoMAG_NORTH_lat_index, IONOtoMAG_SOUTH_lat_index,   &
       IONOtoMAG_NORTH_lon_index, IONOtoMAG_SOUTH_lon_index

  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::                &
       IONOtoMAG_NORTH_lat_factor, IONOtoMAG_SOUTH_lat_factor, &
       IONOtoMAG_NORTH_lon_factor, IONOtoMAG_SOUTH_lon_factor

  real, dimension(:), allocatable ::                          &
       MAG_NORTH_X,MAG_NORTH_Y,MAG_NORTH_Z,MAG_NORTH_R,       & !Magnetospheric coordinates
       MAG_NORTH_Theta,MAG_NORTH_Psi,                         & !
       MAG_SOUTH_X,MAG_SOUTH_Y,MAG_SOUTH_Z,MAG_SOUTH_R,       & !
       MAG_SOUTH_Theta,MAG_SOUTH_Psi,                         & !
       MAG_NORTH_JR,MAG_NORTH_Jx,MAG_NORTH_Jy,MAG_NORTH_Jz,   & !Magnetospheric current
       MAG_SOUTH_JR,MAG_SOUTH_Jx,MAG_SOUTH_Jy,MAG_SOUTH_Jz

  real, dimension(:,:), allocatable ::   &
       MAG_NORTH_MagField, MAG_SOUTH_MagField

  real, dimension(:,:), allocatable ::   &
       MAG_NORTH_IONO_LOC, MAG_SOUTH_IONO_LOC

  !\
  ! Magnetosphere inner boundary current solution variable definitions.
  !/
  integer :: nMagBndPts_North, nMagBndPts_South
  integer, dimension(:), allocatable :: nMagBndPts_North_PE, nMagBndPts_South_PE
  real, dimension(:), allocatable ::             &
       Xmag_North,Ymag_North,Zmag_North,         & ! Magnetospheric coordinates
       Xmag_South,Ymag_South,Zmag_South,         & !
       JXmag_North,JYmag_North,JZmag_North,      & ! Magnetospheric current
       JXmag_South,JYmag_South,JZmag_South

  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::     &
       IONO_NORTH_PHI_BC,IONO_SOUTH_PHI_BC,         & !Magnetosphere bound pot
       IONO_NORTH_ETh_BC,IONO_NORTH_EPs_BC,         & !Magnetosphere bound E field
       IONO_SOUTH_ETh_BC,IONO_SOUTH_EPs_BC,         & !
       IONO_NORTH_UTh_BC,IONO_NORTH_UPs_BC,         & !Magnetosphere bound flow
       IONO_NORTH_UR_BC,                            & !Velocities
       IONO_SOUTH_UTh_BC,IONO_SOUTH_UPs_BC,         &
       IONO_SOUTH_UR_BC

  real, dimension(1:IONO_nTheta,1:IONO_nPsi,1:8) ::  &
       IONO_NORTH_Mapping_Distance, IONO_SOUTH_Mapping_Distance

  integer, dimension(1:IONO_nTheta,1:IONO_nPsi,1:8) ::  &
       IONO_NORTH_Mapping_Index, IONO_SOUTH_Mapping_Index

end module ModIonosphere

!=======================================================================
subroutine IE_get_for_gm_swmf(nPartial,iGetStart,Get,W,State_V,nVar)
  use CON_router
  use ModIonosphere,ONLY:IONO_NORTH_Phi,IONO_SOUTH_Phi
  implicit none
  integer,intent(in)::nPartial,iGetStart,nVar
  type(IndexPtrType),intent(in)::Get
  type(WeightPtrType),intent(in)::W
  real,dimension(nVar),intent(out)::State_V
  integer::iGet

  select case(Get%iCB_II(3,iGetStart))
  case(1)
     State_V(1)=&
          IONO_NORTH_Phi(Get%iCB_II(1,iGetStart),Get%iCB_II(2,iGetStart))*&
          W%Weight_I(iGetStart)
  case(2)
     State_V(1)=&
          IONO_SOUTH_Phi(Get%iCB_II(1,iGetStart),Get%iCB_II(2,iGetStart))*&
           W%Weight_I(iGetStart)
  case default
     call CON_stop('IE_get_for_gm_swmf: wrong block number',&      
          Get%iCB_II(3,iGetStart))
  end select
  do iGet=iGetStart+1,iGetStart+nPartial-1
     select case(Get%iCB_II(3,iGet))
     case(1)
        State_V(1)= State_V(1)+&
             IONO_NORTH_Phi(Get%iCB_II(1,iGet),Get%iCB_II(2,iGet))*&
             W%Weight_I(iGet)
     case(2)
        State_V(1)=State_V(1)+&
             IONO_SOUTH_Phi(Get%iCB_II(1,iGet),Get%iCB_II(2,iGet))*&
             W%Weight_I(iGet)
     case default
        call CON_stop('IE_get_for_gm_swmf: wrong block number',&      
             Get%iCB_II(3,iGet))
     end select
  end do
end subroutine IE_get_for_gm_swmf
!=========================================================================
subroutine IE_put_from_gm_swmf(&
     State_V,&
     nVar,&
     iBlock,&
     iPoint,&
     ColatLim)
  use CON_physics,ONLY:get_planet
  use ModIonosphere
  use IE_ModMain, ONLY: IsNewInput
  use CON_world
  implicit none
  integer,intent(in)::nVar,iBlock,iPoint
  real,dimension(nVar),intent(in)::State_V
  real,intent(inout)::ColatLim
  real,save::Radius
  real :: RadiusPlanet,IonosphereHeight
  logical,save::DoInitialize=.true.

  !  logical :: DoTest,DoTestMe
  !  character (len=*),parameter :: NameSub='IE_put_from_gm_swmf'
  !----------------------------------------------------------------------
  !  call CON_set_do_test(NameSub, DoTest, DoTestMe)
  !  if(DoTest)write(*,*)NameSub,' starting with nVar=',nVar

  if(DoInitialize)then
     call get_planet(                            &
          RadiusPlanetOut     = RadiusPlanet,    &
          IonosphereHeightOut = IonosphereHeight)
     Radius=cOne+IonosphereHeight/RadiusPlanet
     DoInitialize=.false.
  end if

  IsNewInput=.true.

  select case(iBlock)
  case(1)     !North
     MAG_NORTH_Jx(iPoint) = State_V(1)
     MAG_NORTH_Jy(iPoint) = State_V(2)
     MAG_NORTH_Jz(iPoint) = State_V(3)
     if(nVar>3) MAG_NORTH_MagField(iPoint,:) = State_V(4:8)
     if(nVar>8) then
        MAG_NORTH_IONO_LOC(iPoint,:) = State_V(9:11)
        ColatLim=max(ColatLim, acos(State_V(11)/Radius))
     end if
  case(2)
     MAG_SOUTH_Jx(iPoint) = State_V(1)
     MAG_SOUTH_Jy(iPoint) = State_V(2)
     MAG_SOUTH_Jz(iPoint) = State_V(3)
     if(nVar>3) MAG_SOUTH_MagField(iPoint,:) = State_V(4:8)
     if(nVar>8) then
        MAG_SOUTH_IONO_LOC(iPoint,:) = State_V(9:11)
        ColatLim=min(ColatLim,acos(State_V(11)/Radius))
     end if
  case default
     call CON_stop('Impossible block number in IE',iBlock)
  end select
end subroutine IE_put_from_gm_swmf

!==============================================================================

subroutine IE_check_allocation_north(nPoint)
  use ModIonosphere
  use ModUtilities, ONLY: check_allocate
  implicit none
  integer,intent(in)::nPoint
  integer::iError
  if(nPoint /= IONO_NORTH_nMagBndPts)then
     IONO_NORTH_nMagBndPts = nPoint
     if(allocated(MAG_NORTH_IONO_LOC)) deallocate( &
          MAG_NORTH_IONO_LOC, MAG_NORTH_Jx, MAG_NORTH_Jy, MAG_NORTH_Jz, &
          MAG_NORTH_MAGFIELD)
     allocate(&
          MAG_NORTH_IONO_LOC(nPoint,3), &
          MAG_NORTH_Jx(nPoint),MAG_NORTH_Jy(nPoint),MAG_NORTH_Jz(nPoint),&
          MAG_NORTH_MAGFIELD(nPoint,5),STAT=iError)
     call check_allocate(iError,&
          'IE_check_allocation_north loc3,j3,binfo5')
  end if
end subroutine IE_check_allocation_north

!==============================================================================

subroutine IE_check_allocation_south(nPoint)
  use ModIonosphere
  use ModUtilities, ONLY: check_allocate
  implicit none
  integer,intent(in)::nPoint
  integer::iError
  if(nPoint /= IONO_SOUTH_nMagBndPts)then
     IONO_SOUTH_nMagBndPts = nPoint
     if(allocated(MAG_SOUTH_IONO_LOC)) deallocate( &
          MAG_SOUTH_IONO_LOC, MAG_SOUTH_Jx, MAG_SOUTH_Jy, MAG_SOUTH_Jz, &
          MAG_SOUTH_MAGFIELD)
     allocate(&
          MAG_SOUTH_IONO_LOC(nPoint,3), &
          MAG_SOUTH_Jx(nPoint),MAG_SOUTH_Jy(nPoint),MAG_SOUTH_Jz(nPoint),&
          MAG_SOUTH_MAGFIELD(nPoint,5),STAT=iError)
     call check_allocate(iError,&
          'IE_check_allocation_south loc3,j3,binfo5')
  end if
end subroutine IE_check_allocation_south

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

!^CFG COPYRIGHT UM
subroutine IE_interpolate(iBlock,ColatLim)

  ! Establish interpolation scheme between the points provided by GM
  ! and the ionosphere polar grid.

  ! x,y,z MAG_NORTH_IONO_LOC(i,1:3) --> 
  !                            IONO_NORTH_Mapping_Index
  !                            IONO_NORTH_Mapping_Distance

  use ModIonosphere

  implicit none

  integer, intent(in) :: iBlock     ! Block index 1 for north, 2 for south
  real,    intent(in) :: ColatLim   ! Colatitude limit for the hemisphere

  integer              :: i, j, k, l, n_save, n, ierror

  real, dimension(:), allocatable :: Temp_Distance

  real :: r_max, r_min, r

  logical :: done

  logical :: DoTest, DoTestMe
  character(len=*), parameter:: NameSub='IE_interpolate'
  !---------------------------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest,DoTestMe)
  if(DoTest)write(*,*)NameSub,': starting, iBlock, ColatLim=',&
       iBlock,ColatLim


  !
  ! Limit interpolation distances to less than 10 degrees away
  !

  r_max = Max_FAC_Distance

  select case(iBlock)
  case(1) ! North hemisphere
     allocate(Temp_Distance(IONO_NORTH_nMagBndPts), stat = ierror)
     if (ierror > 0) write(*,*) "iono_coupling (north): ", &
          " allocation error for Temp_Distance (", &
          IONO_NORTH_nMagBndPts,")"
     Temp_Distance = 0.0

     do j = 1, IONO_nPsi
        do i = 1, IONO_nTheta

           IONO_NORTH_Mapping_Index(i,j,1:8) = -1
           IONO_NORTH_Mapping_Distance(i,j,1:8) = -1.0

           if (IONO_NORTH_Theta(i,j) <= ColatLim) then

              Temp_Distance(:)=                                          &
                   sqrt((MAG_NORTH_IONO_LOC(:,1)-IONO_NORTH_X(i,j))**2 +  &
                   (MAG_NORTH_IONO_LOC(:,2)-IONO_NORTH_Y(i,j))**2 +  &
                   (MAG_NORTH_IONO_LOC(:,3)-IONO_NORTH_Z(i,j))**2)

              !
              ! Find 8 closest point, which are less than r_max degrees away
              !

              do k = 1, 8

                 r_min = r_max
                 n_save = -1

                 do n = 1, IONO_NORTH_nMagBndPts

                    if ((Temp_Distance(n) < r_min) .and.  &
                         (Temp_Distance(n) > 0.0)   .and.  &
                         (Temp_Distance(n) < r_max)) then
                       r_min = Temp_Distance(n)
                       n_save = n
                    endif

                 enddo

                 if (n_save > 0) then
                    IONO_NORTH_Mapping_Index(i,j,k) = n_save
                    IONO_NORTH_Mapping_Distance(i,j,k) = Temp_Distance(n_save)
                    Temp_Distance(n_save) = -1.0
                 endif

              enddo

           endif

        end do
     end do

     deallocate(Temp_Distance)

  case(2) ! South

     allocate(Temp_Distance(IONO_SOUTH_nMagBndPts), stat = ierror)
     if (ierror > 0) write(*,*) "iono_coupling (south): ", &
          " allocation error for Temp_Distance (", &
          IONO_SOUTH_nMagBndPts,")"
     Temp_Distance = 0.0

     !
     ! South
     !

     do j = 1, IONO_nPsi
        do i = 1, IONO_nTheta

           IONO_SOUTH_Mapping_Index(i,j,1:8) = -1
           IONO_SOUTH_Mapping_Distance(i,j,1:8) = -1.0

           if (IONO_SOUTH_Theta(i,j) >= ColatLim) then

              Temp_Distance(:) =                                    &
                   sqrt((MAG_SOUTH_IONO_LOC(:,1)-IONO_SOUTH_X(i,j))**2 + &
                   (MAG_SOUTH_IONO_LOC(:,2)-IONO_SOUTH_Y(i,j))**2 + &
                   (MAG_SOUTH_IONO_LOC(:,3)-IONO_SOUTH_Z(i,j))**2)

              !
              ! Find 8 closest point, which are less than r_max degrees away
              !

              do k = 1, 8

                 r_min = r_max
                 n_save = -1

                 do n = 1, IONO_SOUTH_nMagBndPts

                    if ((Temp_Distance(n) < r_min) .and.  &
                         (Temp_Distance(n) > 0.0)   .and.  &
                         (Temp_Distance(n) < r_max)) then
                       r_min = Temp_Distance(n)
                       n_save = n
                    endif

                 enddo

                 if (n_save > 0) then
                    IONO_SOUTH_Mapping_Index(i,j,k) = n_save
                    IONO_SOUTH_Mapping_Distance(i,j,k) = Temp_Distance(n_save)
                    Temp_Distance(n_save) = -1.0
                 endif

              enddo

           endif

        end do
     end do

     if(DoTest)write(*,*)NameSub,': IONO_SOUTH_Mapping_Index(50,50,:)=',&
          IONO_SOUTH_Mapping_Index(50,50,:)


     if(DoTest)write(*,*)NameSub,&
          ': ColatLim abs(sum(IONO_SOUTH_Theta)),SOUTH_nMagBndPts',&
          ColatLim,abs(sum(IONO_SOUTH_Theta)),IONO_SOUTH_nMagBndPts

     if(DoTest)write(*,*)NameSub,': sum(abs(real(IONO_SOUTH_Mapping_Index)))',&
          sum(abs(real(IONO_SOUTH_Mapping_Index)))

     if(DoTest)write(*,*)NameSub,': sum(abs(MAG_SOUTH_IONO_LOC))',&
          sum(abs(MAG_SOUTH_IONO_LOC))

     if(DoTest)write(*,*)NameSub,': sum(abs(IONO_SOUTH_X,Y,Z))',&
          sum(abs(IONO_SOUTH_X)),sum(abs(IONO_SOUTH_Y)),sum(abs(IONO_SOUTH_Z))


     deallocate(Temp_Distance)

  end select

end subroutine IE_interpolate
!^CFG COPYRIGHT UM
subroutine ionosphere_magBCs(PHI_BC, ETh_BC, EPs_BC,                         &
                             UR_BC, UTh_BC, UPs_BC, Radius_BC,               &
                             PHI, X, Y, Z,                                   &
                             Theta_BC, Psi_BC, Radius, nTheta, nPsi,         &
                             dTheta, dPsi)

  use ModIonosphere
  implicit none

  integer :: nTheta, nPsi
  real :: Radius, Radius_BC
  real, dimension(1:IONO_nTheta,1:IONO_nPsi) ::                              &
                  PHI_BC, ETh_BC, EPs_BC,                                    &
                  UR_BC, UTh_BC, UPs_BC,                                     &
                  PHI, X, Y, Z, Theta_BC, Psi_BC

  real, dimension(1:IONO_nTheta,1:IONO_nPsi,2) ::                            &
                  Theta_temp, Psi_temp, phi_temp

  real, dimension(1:IONO_nTheta) :: dTheta
  real, dimension(1:IONO_nPsi)   :: dPsi

  ! This routine does nothing.

end subroutine ionosphere_magBCs
