!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!========================================================================
module ModUser

  use ModUserEmpty,               &
       IMPLEMENTED1 => user_set_boundary_cells,         &
       IMPLEMENTED2 => user_read_inputs,                &
       IMPLEMENTED3 => user_init_session

  use ModSize
  use ModProcMH, ONLY: iProc

  include 'user_module.h' !list of public methods

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'CG Comet, G. Toth & H. Zhenguang, 2014'

  character (len=100) :: NameShapeFile

  integer:: nTriangle
  real, allocatable:: XyzTriangle_DII(:,:,:), Normal_DI(:,:)
  real :: rMinShape = 0.0, rMaxShape = 0.0

contains
  !============================================================================
  subroutine user_read_inputs

    use ModMain, ONLY: UseUserInitSession, UseExtraBoundary
    use ModReadParam

    character (len=100) :: NameCommand

    character(len=*), parameter:: NameSub = 'user_read_inputs'
    !-------------------------------------------------------------------------

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)

       case("#SHAPEFILE")
          call read_var('NameShapeFile' ,NameShapeFile)
       case('#USERINPUTEND')
          EXIT
       case default
          if(iProc==0) call stop_mpi( &
               NameSub//' invalid command='//trim(NameCommand))
       end select
    end do

    UseUserInitSession = .true.
    UseExtraBoundary   = .true.

  end subroutine user_read_inputs
  !===========================================================================
  subroutine user_init_session

    ! We need to have unit conversions before reading the shape file 
    ! which contains everything in SI units
    call read_shape_file

  end subroutine user_init_session
  !===========================================================================
  subroutine user_set_boundary_cells(iBlock)

    use ModGeometry, ONLY: ExtraBc_, IsBoundaryCell_GI, Xyz_DGB, r_BLK
    integer, intent(in):: iBlock

    integer:: i, j, k
    !------------------------------------------------------------------------

    do k = MinK, MaxK; do j = MinJ, MaxJ; do i=MinI, MaxI
       ! Check if we are close enough
       if(r_BLK(i,j,k,iBlock) > rMaxShape) then
          IsBoundaryCell_GI(i,j,k,ExtraBc_) = .false.
       elseif(r_BLK(i,j,k,iBlock) < rMinShape) then
          IsBoundaryCell_GI(i,j,k,ExtraBc_) = .true.
       else
          ! Connect cell center with origin which should be inside the shape. 
          ! If the line segment does not intersect the shape or it intersects
          ! even times then the point is inside the shape.
          IsBoundaryCell_GI(i,j,k,ExtraBc_) = .not. &
               is_segment_intersected((/0.,0.,0./), Xyz_DGB(:,i,j,k,iBlock),&
               IsOddIn=.true.)
       end if
    end do; end do; end do

  end subroutine user_set_boundary_cells

  !=========================================================================
  subroutine read_shape_file

    use ModPhysics, ONLY: Si2No_V, UnitX_
    use ModIoUnit, ONLY: UnitTmp_
    use ModCoordTransform, ONLY: cross_product

    logical:: DoReadShapeFile = .true.

    integer:: nPoint, i, j, iPoint, iTriangle, iPoint1, iPoint2, iPoint3

    real, allocatable:: Xyz_DI(:,:)

    character(len=100):: String1, String2

    character(len=*), parameter:: NameSub = 'read_shape_file'
    !-----------------------------------------------------------------------
    if(.not.DoReadShapeFile) RETURN
    DoReadShapeFile = .false.

    if(iProc==0)write(*,*) NameSub,' reading shape file ',trim(NameShapeFile)

    open(UnitTmp_, file=NameShapeFile)

    read(UnitTmp_, '(a)') String1
    read(UnitTmp_, *) String1, nPoint, String2
    if(String2 /= 'POINTS')call stop_mpi(NameSub//&
         ' no POINTS in '//trim(String2)//' of '//trim(NameShapeFile))
    read(UnitTmp_, *) String1, nTriangle, String2
    if(String2 /= 'TRIANGLES')call stop_mpi(NameSub//&
         ' no TRIANGLES in '//trim(String2)//' of '//trim(NameShapeFile))

    if(iProc==0)write(*,*) NameSub,' nPoint=', nPoint,' nTriangle=',nTriangle

    allocate(Xyz_DI(3,nPoint), &
         XyzTriangle_DII(3,3,nTriangle), Normal_DI(3,nTriangle))

    read(UnitTmp_, '(a)') String1
    do iPoint = 1, nPoint
       read(UnitTmp_,*) String1, i, j, Xyz_DI(:,iPoint)

       Xyz_DI(:,iPoint) = Xyz_DI(:,iPoint) * Si2No_V(UnitX_)
    end do
    do iTriangle = 1, nTriangle
       read(UnitTmp_,*) String1, i, j, iPoint1, iPoint2, iPoint3
       XyzTriangle_DII(:,1,iTriangle) = Xyz_DI(:,iPoint1)
       XyzTriangle_DII(:,2,iTriangle) = Xyz_DI(:,iPoint2)
       XyzTriangle_DII(:,3,iTriangle) = Xyz_DI(:,iPoint3)

       Normal_DI(:,iTriangle) = cross_product( &
            Xyz_DI(:,iPoint2) - Xyz_DI(:,iPoint1), &
            Xyz_DI(:,iPoint3) - Xyz_DI(:,iPoint1))
       Normal_DI(:,iTriangle) = Normal_DI(:,iTriangle) / &
            sqrt(sum(Normal_DI(:,iTriangle)**2))
    end do

    !write(*,*)'!!! XyzTriangle_DII(:,:,1)=',XyzTriangle_DII(:,:,1)
    !write(*,*)'!!! XyzTriangle_DII(:,:,n)=',XyzTriangle_DII(:,:,nTriangle)

    rMinShape = sqrt(minval(sum(Xyz_DI**2,DIM=1)))
    rMaxShape = sqrt(maxval(sum(Xyz_DI**2,DIM=1)))

    if(iProc==0)write(*,*) NameSub,' rMinShape, rMaxShape=', &
         rMinShape, rMaxShape

    deallocate(Xyz_DI)

    close(UnitTmp_)

  end subroutine read_shape_file

  !============================================================================
  logical function is_segment_intersected(Xyz1_D, Xyz2_D, IsOddIn)

    ! Check if a line segment connecting Xyz1_D and Xyz2_D intersects
    ! the shape. If IsEvenIn is present and true, check if the number
    ! of intersections is odd or even.

    ! Using algorithm: http://geomalgorithms.com/a06-_intersect-2.html

    real, intent(in):: Xyz1_D(3) ! segment start coordinates
    real, intent(in):: Xyz2_D(3) ! segment end   coordinates
    logical, optional:: IsOddIn  ! check for odd number of intersects

    logical:: IsOdd

    integer:: iTriangle, nIntersect

    integer, parameter:: MaxIntersect = 10
    real:: Ratio, Ratio1, Ratio2, Ratio_I(MaxIntersect)
    real, dimension(3):: v1_D, Xyz_D, u_D, v_D, w_D
    real:: nDotP2P1, nDotV1P1, u2, v2, uDotV, wDotU, wDotV, InvDenom

    character(len=*), parameter:: NameSub = 'is_segment_intersected'
    !------------------------------------------------------------------------
    ! Default is to check for first intersection only
    ! (for shadows and convex shapes)
    IsOdd = .false.
    if(present(IsOddIn)) IsOdd = IsOddIn

    nIntersect = 0
    do iTriangle = 1, nTriangle
       ! Find intersection of line segment with the plane of the triangle
       nDotP2P1 = sum(Normal_DI(:,iTriangle)*(Xyz2_D - Xyz1_D))

       ! Check if the segment is parallel to the plane
       if(abs(nDotP2P1) < 1e-20) CYCLE

       ! Vertex 1 of triangle
       v1_D = XyzTriangle_DII(:,1,iTriangle)

       nDotV1P1 = sum(Normal_DI(:,iTriangle)*(v1_D -  Xyz1_D))

       ! Intersection is at P1 + Ratio * (P2 - P1)
       Ratio = nDotV1P1 / nDotP2P1

       ! Check if point is inside the segment
       if(Ratio < 0.0) CYCLE
       if(Ratio > 1.0) CYCLE

       ! Intersection point
       Xyz_D =  Xyz1_D + Ratio*(Xyz2_D -  Xyz1_D)

       ! Calculate the barycentric coordinates Ratio1 and Ratio2
       ! The intersection point is inside if the conditions
       ! 0 < Ratio1, Ratio2 and Ratio1 + Ratio2 < 1 both hold.

       ! Vectors relative to the first vertex
       u_D = XyzTriangle_DII(:,2,iTriangle) - v1_D 
       v_D = XyzTriangle_DII(:,3,iTriangle) - v1_D 
       w_D = Xyz_D                          - v1_D

       u2 = sum(u_D**2)
       v2 = sum(v_D**2)
       uDotV = sum(u_D*v_D)
       wDotU = sum(w_D*u_D)
       wDotV = sum(w_D*v_D)

       InvDenom = 1.0/(uDotV**2 - u2*v2)

       ! First barycentric coordinate
       Ratio1 = (uDotV*wDotV - v2*wDotU)*InvDenom
       if(Ratio1 < 0.0) CYCLE
       if(Ratio1 > 1.0) CYCLE

       ! Second barycentric coordinate
       Ratio2 = (uDotV*wDotU - u2*wDotV)*InvDenom
       if(Ratio2 < 0.0) CYCLE
       if(Ratio1 + Ratio2 > 1.0) CYCLE

       if(.not. IsOdd)then
          ! The line segment intersected the triangle
          is_segment_intersected = .true.
          RETURN
       end if

       ! Check if this intersection is different from previous ones
       if(nIntersect > 0)then
          if( any(abs(Ratio - Ratio_I(1:nIntersect)) < 1e-12) ) CYCLE
       end if

       ! New intersection was found
       nIntersect =  nIntersect + 1

       if(nIntersect > MaxIntersect) call stop_mpi(NameSub// &
            ': too many intersections, increase MaxIntersect')

       ! Store the position along the segment into Ratio_I
       Ratio_I(nIntersect) = Ratio

    end do

    if(.not. IsOdd)then
       ! The line segment was not intersected by any triangle
       is_segment_intersected = .false.
       RETURN
    end if

    ! We got to the other side of the surface
    ! if there were odd number of intersections
    is_segment_intersected = modulo(nIntersect, 2) == 1

  end function is_segment_intersected

end module ModUser
