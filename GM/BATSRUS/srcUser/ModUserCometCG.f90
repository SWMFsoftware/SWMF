!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!========================================================================
module ModUser

  use ModUserEmpty,               &
       IMPLEMENTED1 => user_set_boundary_cells,         &
       IMPLEMENTED2 => user_read_inputs,                &
       IMPLEMENTED3 => user_init_session,               &
       IMPLEMENTED4 => user_set_face_boundary

  use ModSize
  use ModProcMH, ONLY: iProc
  use ModMain, ONLY: xTest, yTest, zTest
  use ModVarIndexes, ONLY: nVar
  use ModNumConst, ONLY: cPi

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

  ! Position of the sun at the start time
  real :: LatSun=43.5, LonSun=0.0, NormalSun_D(3)

  ! Rotation of the comet (changes the direction of the Sun)
  real:: RotationCometHour = 12.0
  
  ! Angular velocity
  real:: OmegaComet

  ! Minimum and maximum density
  real :: RhoCometMinDim=5.0, RhoCometMaxDim=1e6

  ! velocity amplitude and temperature of outflow in IO units
  real :: uCometDim=1.0, TempCometDim=1e3

  ! Same in normalized units
  real :: RhoCometMin, RhoCometMax, uComet, TempComet

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
       case("#SUNDIRECTION")
          call read_var('LatSun', LatSun)
          call read_var('LonSun', LonSun)
       case("#COMETSTATE")
          call read_var('RhoCometMinDim', RhoCometMinDim)
          call read_var('RhoCometMaxDim', RhoCometMaxDim)
          call read_var('uCometDim',      uCometDim)
          call read_var('TCometDim',      TempCometDim)
       case("#COMETROTATION")
          call read_var('RotationCometHour', RotationCometHour)
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

    ! Read shape file and convert units

    use ModPhysics, ONLY: Io2No_V, Si2No_V, &
         UnitRho_, UnitU_, UnitTemperature_, UnitT_
    use ModNumConst, ONLY: cTwoPi, cDegToRad
    use ModCoordTransform, ONLY: dir_to_xyz
    !------------------------------------------------------------------------
    ! We need to have unit conversions before reading the shape file 
    ! which contains everything in SI units
    call read_shape_file

    RhoCometMin  = RhoCometMinDim  * Io2No_V(UnitRho_)
    RhoCometMax  = RhoCometMaxDim  * Io2No_V(UnitRho_)
    uComet    = uCometDim    * Io2No_V(UnitU_)
    TempComet = TempCometDim * Io2No_V(UnitTemperature_)

    OmegaComet = cTwoPi / (RotationCometHour*3600 * Si2No_V(UnitT_))

    call dir_to_xyz((90-LatSun)*cDegToRad, LonSun*cDegToRad, NormalSun_D)

    if(iProc==0)then
       write(*,*)'RhoCometMinDim,  RhoCometMin =', RhoCometMinDim, RhoCometMin
       write(*,*)'RhoCometMaxDim,  RhoCometMax =', RhoCometMaxDim, RhoCometMax
       write(*,*)'uCometDim,    uComet   =', uCometDim, uComet
       write(*,*)'TempCometDim, TempComet=', TempCometDim, TempComet
       write(*,*)'RotationComet, Omega   =', RotationCometHour, OmegaComet
       write(*,*)'LatSun, LonSun, NormalSun_D=', LatSun, LonSun, NormalSun_D
    end if

  end subroutine user_init_session
  !===========================================================================
  subroutine user_set_boundary_cells(iBlock)

    use ModMain, ONLY: n_step
    use ModGeometry, ONLY: ExtraBc_, IsBoundaryCell_GI, Xyz_DGB, r_BLK
    integer, intent(in):: iBlock

    integer:: i, j, k
    real:: XyzInside_D(3)
    !------------------------------------------------------------------------
    if(n_step > 1)then
       ! This should not happen !!!
       RETURN
    end if
    ! Place a point inside rMinShape sphere with transcendent coordinates
    ! to reduce chances of hitting the edge or corner of triangles
    
    XyzInside_D = rMinShape*(/cPi/10,cPi**2/50,cPi**3/700/)

    do k = MinK, MaxK; do j = MinJ, MaxJ; do i=MinI, MaxI
       ! Check if we are close enough
       if(r_BLK(i,j,k,iBlock) > rMaxShape) then
          IsBoundaryCell_GI(i,j,k,ExtraBc_) = .false.
       elseif(r_BLK(i,j,k,iBlock) < rMinShape) then
          IsBoundaryCell_GI(i,j,k,ExtraBc_) = .true.
       else
          ! Connect cell center with a point inside.
          ! If the line segment does not intersect the shape or it intersects
          ! even times then the point is inside the shape.
          IsBoundaryCell_GI(i,j,k,ExtraBc_) = .not. is_segment_intersected( &
               XyzInside_D, Xyz_DGB(:,i,j,k,iBlock), IsOddIn=.true.)
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

  !==============================================================================

  subroutine user_set_face_boundary(VarsGhostFace_V)

    use ModMain, ONLY: n_step
    use ModVarIndexes,   ONLY: nVar, Rho_, p_, Ux_, Uz_, MassFluid_I
    use ModFaceBoundary, ONLY: iFace, jFace, kFace, FaceCoords_D, &
         iBoundary, VarsTrueFace_V, iSide, iBlockBc
    use ModGeometry,    ONLY: ExtraBc_
    use BATL_lib, ONLY: CellSize_DB
    use ModGeometry, ONLY: Xyz_DGB
    use ModConst, ONLY: cBoltzmann

    logical :: DoTestHere=.true.

    real, intent(out):: VarsGhostFace_V(nVar)

    integer:: iTrue, jTrue, kTrue, iBody, jBody, kBody

    real :: XyzIntersect_D(3), XyzStart_D(3), XyzEnd_D(3)
    real :: TestFace_D(3)
    real :: XyzTrueCell_D(3), XyzBodyCell_D(3)
    real :: Normal_D(3), CosAngle

    real, save:: VarsGhostFace_VFB(nVar,nI+1,nJ+1,nK+1,MaxBlock)

    character(len=*), parameter:: NameSub = 'user_set_face_boundary'
    !------------------------------------------------------------------------
    if(n_step > 1)then
       VarsGhostFace_V = VarsGhostFace_VFB(:,iFace,jFace,kFace,iBlockBc)
       RETURN
    end if

    if (iBoundary /= ExtraBc_) call stop_mpi(NameSub//' is implemented for extra BC only')

    ! Floating boundary condition by default
    VarsGhostFace_V = VarsTrueFace_V

    ! Default indexes for the true and body cells
    iTrue = iFace; jTrue = jFace; kTrue = kFace
    iBody = iFace; jBody = jFace; kBody = kFace

    select case(iSide)
    case(1)
       iBody = iFace - 1
    case(2)
       iTrue = iFace - 1
    case(3)
       jBody = jFace - 1
    case(4)
       jTrue = jFace - 1
    case(5)
       kBody = kFace -1
    case(6)
       kTrue = kFace -1
    end select

    XyzBodyCell_D = Xyz_DGB(:,iBody,jBody,kBody,iBlockBc)
    XyzTrueCell_D = Xyz_DGB(:,iTrue,jTrue,kTrue,iBlockBc)

    ! Find the intersection point between the true cell and the body cell
    ! that is closest to the true cell
    if (.not. is_segment_intersected(XyzTrueCell_D, XyzBodyCell_D, IsOddIn = .true., &
         XyzIntersectOut_D=XyzIntersect_D, NormalOut_D = Normal_D))then
       write(*,*) NameSub,' error for face =', iFace, jFace, kFace
       write(*,*) NameSub,' error for iside, iBlockBc=', iSide, iBlockBc
       call stop_mpi(NameSub// &
            ': No intersection points are found between true and the body cells')
    end if

    ! Fix the normal direction if it is not pointing outward
    if (sum(Normal_D*(XyzTrueCell_D - XyzBodyCell_D)) < 0.0) &
         Normal_D = -Normal_D

    ! Set outflow from shade as default that may be overwritten below
    VarsGhostFace_V(Rho_)    = RhoCometMin
    VarsGhostFace_V(Ux_:Uz_) = Normal_D*uComet
    VarsGhostFace_V(P_)      = VarsGhostFace_V(Rho_)/MassFluid_I(1)*TempComet

    ! Check if Sun light hits the shape
    CosAngle = sum(Normal_D*NormalSun_D)

    if (CosAngle > 0.0) then

       ! See whether the intersection point is in the shade by going towards the Sun
       ! and checking for intersection with the shape
       XyzStart_D = XyzIntersect_D + 1e-9*rMaxShape*NormalSun_D
       XyzEnd_D   = XyzIntersect_D +    2*rMaxShape*NormalSun_D
       if(.not.is_segment_intersected(XyzStart_D, XyzEnd_D)) then

          ! Increase density
          VarsGhostFace_V(Rho_)    = RhoCometMin + (RhoCometMax - RhoCometMin)*CosAngle

          ! Increase pressure
          VarsGhostFace_V(P_)      = VarsGhostFace_V(Rho_)/MassFluid_I(1)*TempComet

          if (DoTestHere)  then
             write(*,*) 'FaceCoords_D: ', FaceCoords_D
             write(*,*) 'TestFace_D: ', (XyzBodyCell_D + XyzTrueCell_D)/2
             write(*,*) 'XyzTrueCell_D =', XyzTrueCell_D
             write(*,*) 'XyzBodyCell_D =', XyzBodyCell_D
             write(*,*) 'XyzIntersect_D=', XyzIntersect_D
             write(*,*) 'XyzStart_D    =', XyzStart_D 
             write(*,*) 'XyzEnd_D      =', XyzEnd_D 
             write(*,*) 'Normal_D      =', Normal_D
             write(*,*) 'CosAngle      =', CosAngle
             write(*,*) 'Rho           =', VarsGhostFace_V(Rho_)
             write(*,*) 'u_D           =', VarsGhostFace_V(Ux_:Uz_)
             write(*,*) 'p             =', VarsGhostFace_V(p_)
             DoTestHere=.false.
          end if
       end if
    end if

    ! Store for future time steps
    VarsGhostFace_VFB(:,iFace,jFace,kFace,iBlockBc) = VarsGhostFace_V

  end subroutine user_set_face_boundary

  !============================================================================
  logical function is_segment_intersected(Xyz1_D, Xyz2_D, &
       IsOddIn, XyzIntersectOut_D, NormalOut_D)

    ! Check if a line segment connecting Xyz1_D and Xyz2_D intersects
    ! the shape. If IsEvenIn is present and true, check if the number
    ! of intersections is odd or even.

    ! Using algorithm: http://geomalgorithms.com/a06-_intersect-2.html

    real, intent(in):: Xyz1_D(3) ! segment start coordinates
    real, intent(in):: Xyz2_D(3) ! segment end   coordinates
    logical, optional:: IsOddIn  ! check for odd number of intersects
    real,    optional :: XyzIntersectOut_D(3)
    real,    optional :: NormalOut_D(3)

    logical:: IsOdd

    integer:: iTriangle, nIntersect, iIntersect, iMinRatio

    integer, parameter:: MaxIntersect = 10
    real:: Ratio, Ratio1, Ratio2, Ratio_I(MaxIntersect)
    integer:: iTriangle_I(MaxIntersect)
    real, dimension(3):: v1_D, Xyz_D, u_D, v_D, w_D
    real:: nDotP2P1, nDotV1P1, u2, v2, uDotV, wDotU, wDotV, InvDenom

    character(len=*), parameter:: NameSub = 'is_segment_intersected'
    !------------------------------------------------------------------------
    ! Default is to check for first intersection only
    ! (for shadows and convex shapes)
    IsOdd = .false.
    if(present(IsOddIn))   IsOdd = IsOddIn

    nIntersect = 0
    do iTriangle = 1, nTriangle
       ! Find intersection of line segment with the plane of the triangle
       nDotP2P1 = sum(Normal_DI(:,iTriangle)*(Xyz2_D - Xyz1_D))

       ! Check if the segment is parallel to the plane
       if(abs(nDotP2P1) < 1e-12) then
          if (abs(sum(Normal_DI(:,iTriangle)*(Xyz1_D - XyzTriangle_DII(:,1,iTriangle)))) &
               < 1e-12) then
             !write(*,*) 'segment lies in the same plane: iTriangle: ', Xyz2_D
             !write(*,*) 'Test:', &
             !     sum(Normal_DI(:,iTriangle)*(Xyz1_D - XyzTriangle_DII(:,1,iTriangle)))
             CYCLE
          else 
             CYCLE
          end if
       end if

       ! Vertex 1 of triangle
       v1_D = XyzTriangle_DII(:,1,iTriangle)

       nDotV1P1 = sum(Normal_DI(:,iTriangle)*(v1_D -  Xyz1_D))

       ! Intersection is at P1 + Ratio * (P2 - P1)
       Ratio = nDotV1P1 / nDotP2P1

       ! Check if point is inside the segment
       if(Ratio <    -1e-12) CYCLE
       if(Ratio > 1.0+1e-12) CYCLE

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
       if(Ratio1 <    -1e-12) CYCLE
       if(Ratio1 > 1.0+1e-12) CYCLE

       ! Second barycentric coordinate
       Ratio2 = (uDotV*wDotU - u2*wDotV)*InvDenom

!       if (abs(Xyz2_D(1) -xTest) <= 1e-4 .and. &
!            abs(Xyz2_D(2) - yTest) <= 1e-4 .and. &
!            abs(Xyz2_D(3) - zTest) <= 1e-4) then
!          if (abs(Ratio2) < 1e-10 .or. abs(Ratio2 - 1.0) < 1e-10) then
!             write(*,*) 'iTriangle: ', iTriangle, 'Ratio1: ', Ratio1, 'Ratio2: ', Ratio2
!          end if
!       end if

       if(         Ratio2 <    -1e-12) CYCLE
       if(Ratio1 + Ratio2 > 1.0+1e-12) CYCLE

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
       iTriangle_I(nIntersect) = iTriangle
    end do

    ! Only record the closeset intersection to the true cell
    if(present(NormalOut_D)) then
       iMinRatio = minloc(Ratio_I(1:nIntersect), 1)
       iTriangle = iTriangle_I(iMinRatio)
       NormalOut_D = Normal_DI(:,iTriangle)

       if (nIntersect > 1) then          
!          write(*,*) 'nIntersect: ', nIntersect, 'Ratio_I: ', Ratio_I(1:nIntersect)
!          write(*,*) 'Xyz1_D: ', Xyz1_D
!          write(*,*) 'Xyz2_D: ', Xyz2_D
!          do iIntersect = 1, nIntersect
!             iTriangle = Triangle_I(iIntersect)
!             write(*,*) XyzTriangle_DII(:,1,iTriangle)
!             write(*,*) XyzTriangle_DII(:,2,iTriangle)
!             write(*,*) XyzTriangle_DII(:,3,iTriangle)
!          end do
          write(*,*)  'Ratio_I(iMinRatio)            = ', Ratio_I(iMinRatio) 
          write(*,*)  'minval(Ratio_I(1:nIntersect)) = ', minval(Ratio_I(1:nIntersect))
          write(*,*)  'Ratio_I(1:nIntersect)', Ratio_I(1:nIntersect)
       end if
    end if
    if(present(XyzIntersectOut_D)) then
       XyzIntersectOut_D = Xyz1_D + minval(Ratio_I(1:nIntersect))*(Xyz2_D -  Xyz1_D)
    end if
    
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
