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
