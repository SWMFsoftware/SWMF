module ModLevelSet

  implicit none
  
  private ! except

  public:: init_levelset        ! allocate arrays
  public:: clean_levelset       ! deallocate arrays
  public:: get_levelset         ! set levelset values at a point 
  public:: test_levelset        ! unit test

  ! Local variables

  integer, parameter:: MaxSegment = 1000

  integer:: nSegment = 0
  real,    allocatable:: LengthSegment_I(:)
  real,    allocatable:: StartSegment_DI(:,:)
  real,    allocatable:: VectorSegment_DI(:,:)
  integer, allocatable:: iLevelSegment_SI(:,:)

contains

  !===========================================================================

  subroutine init_levelset

    if(allocated(LengthSegment_I)) RETURN
    allocate(                            &
         LengthSegment_I(MaxSegment),    &
         StartSegment_DI(2,MaxSegment),  &
         VectorSegment_DI(2,MaxSegment), &
         iLevelSegment_SI(2,MaxSegment))

  end subroutine init_levelset

  !===========================================================================

  subroutine clean_levelset

    nSegment = 0

    if(.not. allocated(LengthSegment_I)) RETURN
    deallocate(            &
         LengthSegment_I,  &
         StartSegment_DI,  &
         VectorSegment_DI, &
         iLevelSegment_SI)

  end subroutine clean_levelset

  !===========================================================================

  subroutine read_levelset_param(nMaterial, NameMaterial_I, MinLevel)

    use ModReadParam, ONLY: read_var

    integer,          intent(in):: nMaterial
    character(len=*), intent(in):: NameMaterial_I(nMaterial)
    integer,          intent(in):: MinLevel

    character(len=100):: StringPoint

    integer:: i, n, iError
    character(len=100):: StringSegment
    real:: Start_D(2), End_D(2)
    character(len=2):: NameMaterial1, NameMaterial2

    character(len=*), parameter:: NameSub='read_levelset_param'
    !------------------------------------------------------------------------
    call read_var('nSegment', n)     ! Number of new segments
    do i = 1, n
       nSegment = nSegment + 1
       call read_var('Material1 Material2 x1 y1 x2 y2', StringSegment)
       read(StringSegment,*,IOSTAT=iError) &
            NameMaterial1, NameMaterial2, Start_D, End_D
       if(iError /= 0)call CON_stop(NameSub// &
            ' could not read Material1 Material2 x1 y1 x2 y2 from '// &
            trim(StringSegment))
       
       iLevelSegment_SI(1,nSegment) = i_level_material(NameMaterial1)
       iLevelSegment_SI(2,nSegment) = i_level_material(NameMaterial2)
       
       StartSegment_DI(:,nSegment)  = Start_D
       VectorSegment_DI(:,nSegment) = End_D - Start_D
       LengthSegment_I(nSegment)    = sqrt(sum( (End_D - Start_D)**2 ))
    end do

  contains
    !=========================================================================
    integer function i_level_material(NameMaterial)

      character(len=*), intent(in):: NameMaterial
      integer:: iMaterial
      !-----------------------------------------------------------------------
      do iMaterial = 1, nMaterial
         if(NameMaterial == NameMaterial_I(iMaterial)) then
            i_level_material = MinLevel - 1 + iMaterial
            RETURN
         end if
      end do

      call CON_stop(NameSub// &
           ' could not match material name='//NameMaterial)

    end function i_level_material

  end subroutine read_levelset_param

  !===========================================================================

  subroutine get_levelset(CoordIn_D, nMaterial, LevelSet_V, DoTest)

    real,    intent(in) :: CoordIn_D(2)          ! Coordinates of point
    integer, intent(in) :: nMaterial             ! Number of materials
    real,    intent(out):: LevelSet_V(nMaterial) ! Levelset functions

    logical, optional, intent(in):: DoTest

    ! Calculate levelset functions as the distance to the closest
    ! interface for all materials for a point located at CoordIn_D.

    integer:: iSegment
    integer:: iLevel        ! Level index on the side of the point
    integer:: jLevel        ! Level index on the opposite side
    real:: Coord_D(2), Segment_D(2), Point_D(2)
    real:: Length, Projection, Coeff
    real:: d                ! Distance between segment and point

    real:: dNormal         ! Normal distance square between line and point
!    real:: d2Normal_V(:)    ! Normal distance square for each material

    character(len=*), parameter:: NameSub = 'get_levelset'
    !-------------------------------------------------------------------------

!    allocate(dNormal_V(nMaterial))

    ! Initialize levelset
    LevelSet_V = -huge(1.0)

    ! Check distance from all the material interface segments
    do iSegment = 1, nSegment

       ! point coordinate relative to the starting point
       Coord_D = CoordIn_D - StartSegment_DI(:,iSegment)

       ! Segment vector
       Segment_D = VectorSegment_DI(:,iSegment)

       ! Length of the segment
       Length = LengthSegment_I(iSegment)

       ! Distance of the point projected onto the line from point 1
       Projection = sum( Segment_D*Coord_D ) / Length

       ! Normal distance
       dNormal = sqrt(max(0.0,sum(Coord_D**2) - Projection**2))

       ! Find projected point relative to Point1 limited by the segments
       Point_D = min(1.0, max(0.0, Projection/Length))*Segment_D

       ! Distance to the projected point
       d = sqrt( sum((Coord_D - Point_D)**2) )

       ! Avoid problems around segments meeting at acute angles:
       ! when two segments are at the same distance from the point
       ! because their end points coincide, the one with the LARGER
       ! normal distance should be used. We slightly modify the distance
       ! (only when d is different from dNormal) to achieve this
       d = d + 1e-6*(d-dNormal)

       ! Get the levels on the two sides of the segment.

       ! use cross product to figure out which side
       if( Coord_D(1)*Segment_D(2) - Coord_D(2)*Segment_D(1) < 0.0)then
          iLevel = iLevelSegment_SI(1,iSegment)
          jLevel = iLevelSegment_SI(2,iSegment)
       else
          iLevel = iLevelSegment_SI(2,iSegment)
          jLevel = iLevelSegment_SI(1,iSegment)
       end if

       ! Set level set function if distance is smaller than earlier values
       if(d < abs(LevelSet_V(iLevel))) LevelSet_V(iLevel) =  d
       if(d < abs(LevelSet_V(jLevel))) LevelSet_V(jLevel) = -d

       if(.not.present(DoTest)) CYCLE
       if(.not.DoTest)          CYCLE

       write(*,*)NameSub,' iSegment, CoordIn_D=', iSegment, CoordIn_D
       write(*,*)NameSub,' Length,StartSegment=', &
            Length, StartSegment_DI(:,iSegment)
       write(*,*)NameSub,' Segment_D,Coord_D  =', Segment_D, Coord_D
       write(*,*)NameSub,' Projection, Point_D=', Projection, Point_D
       write(*,*)NameSub,' d, iLevel, jLevel  =', d, iLevel, jLevel

    end do

!    deallocate(dNormal_V(nMaterial))

  end subroutine get_levelset

  !=========================================================================

  subroutine test_levelset

    use ModPlotFile, ONLY: save_plot_file
    use ModReadParam, ONLY: read_file, read_init, read_line, read_command
    use ModMpi

    integer:: iProc, iError

    integer, parameter:: nI = 100, nJ = 50, nLevel = 5
    real, allocatable:: LevelSet_VC(:,:,:)
    real :: Coord_D(2), d
    integer:: i, j, iSegment, iLevel, jLevel

    character(len=100):: NameCommand
    character(len=2):: NameMaterial_I(nLevel) = &
         (/ 'Xe', 'Be', 'Pl', 'Au', 'Ay' /)

    character(len=*), parameter:: NameSub = 'test_levelset'
    !----------------------------------------------------------------------

    call init_levelset

    allocate(LevelSet_VC(nLevel+1,nI,nJ))

    ! Read in the parameters describing the material interfaces
    call read_file('test_levelset.in')
    call read_init
    do
       if(.not.read_line()) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#MATERIALINTERFACE")
          call read_levelset_param(nLevel, NameMaterial_I, MinLevel=1)
       case default
          call CON_stop(NameSub//': unknown command='//NameCommand)
       end select
    end do

    ! Set the levelset values for all cell centers
    do j = 1, nJ; do i = 1, nI
       ! Take a simple uniform grid (does not matter really)
       Coord_D = (/ i-0.5, j-0.5 /)

       call get_levelset(Coord_D, nLevel, LevelSet_VC(1:nLevel,i,j))

       ! Set the level function based on the maximum level
       LevelSet_VC(nLevel+1,i,j) = maxloc(LevelSet_VC(1:nLevel,i,j), DIM=1)

    end do; end do

    call MPI_COMM_RANK(MPI_COMM_WORLD, iProc, iError)
    if(iProc==0) call save_plot_file('test_levelset.out', &
         NameVarIn = "x y Xe Be Pl Au Ay level nSegment", &
         ParamIn_I = (/ real(nSegment) /), &
         CoordMinIn_D = (/ 0.5, 0.5 /), &
         CoordMaxIn_D = (/ nI-0.5, nJ-0.5 /), &
         VarIn_VII    = LevelSet_VC)

    deallocate(LevelSet_VC)

  end subroutine test_levelset

end module ModLevelSet
