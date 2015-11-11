!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module BATL_grid

  use BATL_mpi, ONLY: iProc
  use BATL_size
  use BATL_tree
  use BATL_geometry
  use BATL_high_order

  implicit none

  SAVE

  private ! except

  public :: init_grid
  public :: clean_grid
  public :: create_grid
  public :: create_grid_block
  public :: fix_grid_res_change
  public :: average_grid_node
  public :: find_grid_block
  public :: interpolate_grid
  public :: interpolate_grid_amr
  public :: interpolate_grid_amr_gc
  public :: check_interpolate_amr_gc
  public :: test_grid
  
  ! Coordinate limits and size of domain inherited from BATL_geometry
  public:: CoordMin_D, CoordMax_D, DomainSize_D

  real, public, allocatable::   &
       CoordMin_DB(:,:),        &    ! Min gen. coordinates of a block domain
       CoordMax_DB(:,:),        &    ! Max gen. coordinates of a block domain
       CellSize_DB(:,:),        &    ! Cell size in gen. coordinates
       CellFace_DB(:,:),        &    ! Cell faces for Cartesian grids
       CellFace_DFB(:,:,:,:,:), &    ! Cell faces for general grids
       CellVolume_B(:),         &    ! Cell volume for Cartesian grids
       CellVolume_GB(:,:,:,:),  &    ! Cell volume for general grids
       Xyz_DGB(:,:,:,:,:),      &    ! Cartesian cell centers coords
       Xyz_DNB(:,:,:,:,:),      &    ! Cartesian node coordinates
       FaceNormal_DDFB(:,:,:,:,:,:),&! Normal face area vector
       CellMetrice_DDG(:,:,:,:,:), & ! Metrics at cell center. Like: dx/dXi.
       CellCoef_DDGB(:,:,:,:,:,:)    ! X,Y,Z to general coord transform coef

  ! If true, cell faces are assumed to be flat polygons formed by the nodes
  logical, public:: IsNodeBasedGrid = .true.

  ! Local variables

  logical :: DoInitializeGrid = .true.

contains
  !============================================================================
  subroutine init_grid(CoordMinIn_D, CoordMaxIn_D, UseRadiusIn, UseDegreeIn)

    use ModNumConst, ONLY: cTwoPi, cDegToRad

    ! The angular coordinate limits should be given in degrees unless
    ! UseDegreeIn is false.
    ! The radial coordinate limits should be given as true radial values even
    ! for logarithmic or strethed radial grids unless UseRadiusIn is false.

    real, intent(in):: CoordMinIn_D(nDim), CoordMaxIn_D(nDim)
    logical, optional, intent(in):: UseRadiusIn
    logical, optional, intent(in):: UseDegreeIn

    logical:: UseRadius, UseDegree
    real   :: Unit
    !-------------------------------------------------------------------------

    if(.not. DoInitializeGrid) RETURN

    DoInitializeGrid = .false.

    UseRadius = .true.
    if(present(UseRadiusIn)) UseRadius = UseRadiusIn

    UseDegree = .true.
    if(present(UseDegreeIn)) UseDegree = UseDegreeIn
    if(UseDegree)then
       Unit = 1.0
    else
       Unit = cDegToRad
    end if

    ! Make sure that the thickness is unity in the ignored dimensions
    CoordMin_D = -0.5
    CoordMax_D = +0.5
    CoordMin_D(1:nDim) = CoordMinIn_D
    CoordMax_D(1:nDim) = CoordMaxIn_D

    ! Set special boundary conditions and convert coordinates
    IsCylindricalAxis = .false.
    if(IsCylindrical .and. .not.IsLogRadius .and. .not.IsGenRadius) &
         IsCylindricalAxis = CoordMin_D(r_) == 0.0

    if(UseRadius)then
       if(IsLogRadius)then
          ! Convert rMin, rMax to log(rMin) log(rMax) for logarithmic radius
          CoordMin_D(r_) = log(CoordMin_D(r_))
          CoordMax_D(r_) = log(CoordMax_D(r_))
       elseif(IsGenRadius)then
          ! Convert rMin, rMax to generalized radial coordinates
          call radius_to_gen(CoordMin_D(r_))
          call radius_to_gen(CoordMax_D(r_))
       end if
    end if

    IsSphericalAxis = .false.
    if(IsSpherical) IsSphericalAxis = CoordMin_D(Theta_) <   0.01*Unit &
         .and.                        CoordMax_D(Theta_) > 179.99*Unit

    IsLatitudeAxis = .false.
    if(IsRLonLat) IsLatitudeAxis    = CoordMin_D(Lat_)   < -89.99*Unit &
         .and.                        CoordMax_D(Lat_)   >  89.99*Unit

    IsAnyAxis = IsCylindricalAxis .or. IsSphericalAxis .or. IsLatitudeAxis

    if(UseDegree)then
       ! Convert degrees to radians for the domain boundaries
       if(IsCylindrical .or. IsSpherical .or. IsRLonLat)then
          CoordMin_D(Phi_) = CoordMin_D(Phi_)*cDegToRad
          CoordMax_D(Phi_) = CoordMax_D(Phi_)*cDegToRad
       end if

       if(IsSpherical .or. IsRLonLat)then
          CoordMin_D(Theta_) = CoordMin_D(Theta_)*cDegToRad
          CoordMax_D(Theta_) = CoordMax_D(Theta_)*cDegToRad
       end if
    end if

    ! Set size of domain (in generalized coordinates)
    DomainSize_D = CoordMax_D - CoordMin_D

    allocate(CoordMin_DB(MaxDim,MaxBlock))
    allocate(CoordMax_DB(MaxDim,MaxBlock))
    allocate(CellSize_DB(MaxDim,MaxBlock))

    allocate(CellFace_DB(MaxDim,MaxBlock))
    if(.not.IsCartesian) &
         allocate(CellFace_DFB(MaxDim,1:nI+1,1:nJ+1,1:nK+1,MaxBlock))

    allocate(CellVolume_B(MaxBlock))
    allocate(CellVolume_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    allocate(Xyz_DGB(MaxDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    allocate(Xyz_DNB(MaxDim,nINode,nJNode,nKNode,MaxBlock))

    if(.not.IsCartesian) then
       allocate(FaceNormal_DDFB(nDim,nDim,1:nI+1,1:nJ+1,1:nK+1,MaxBlock))

    end if

    ! Periodicity in the radial direction is not possible at all
    if(r_ > 0) IsPeriodic_D(r_) = .false.

    if(Theta_ > 0 .and. (IsSphericalAxis .or. IsLatitudeAxis)) &
         IsPeriodic_D(Theta_) = .false.

    if(Phi_ > 0)then
       ! Enforce periodicity for cylindrical and spherical grids if
       ! there is a full grid in the Phi direction. 
       ! One can also have periodicity with a segment in Phi
       if( abs(DomainSize_D(Phi_) - cTwoPi) < 1e-6 )then
          IsPeriodic_D(Phi_) = .true.
          IsPeriodicCoord_D(Phi_) = .true.
       end if

       ! Set logical for the sign of the minimum Phi coordinate
       IsNegativePhiMin = CoordMin_D(Phi_) < 0.0

    end if

    if(IsRoundCube) IsPeriodic_D = .false.

  end subroutine init_grid
  !===========================================================================
  subroutine clean_grid

    if(DoInitializeGrid) RETURN

    DoInitializeGrid = .true.
    IsNodeBasedGrid  = .true.

    deallocate(CoordMin_DB, CoordMax_DB, CellSize_DB, CellFace_DB, &
         CellVolume_B, Xyz_DGB)
    if(allocated(CellVolume_GB))   deallocate(CellVolume_GB)
    if(allocated(CellFace_DFB))    deallocate(CellFace_DFB)
    if(allocated(Xyz_DNB))         deallocate(Xyz_DNB)
    if(allocated(FaceNormal_DDFB)) deallocate(FaceNormal_DDFB)
    if(allocated(CellCoef_DDGB))   deallocate(CellCoef_DDGB)

    CoordMin_D =  0.0
    CoordMax_D = -1.0
    DomainSize_D = -1.0

  end subroutine clean_grid

  !===========================================================================

  subroutine create_grid_block(iBlock, iNodeIn, DoFixFace, DoFaceOnly)

    use ModNumConst, ONLY: cHalfPi
    use ModCoordTransform, ONLY: cross_product

    ! Create geometrical information for block iBlock on the local PE

    integer, intent(in):: iBlock

    ! In case iNode_B is not set, iNodeIn can provide the node info
    integer, optional, intent(in):: iNodeIn

    ! If DoFixFace is present, fix face area (vectors) at resolution changes
    ! in 3D curvilinear grids. This only works when the neighbor information 
    ! is already known, which happens for initial refinement.
    logical, optional, intent(in):: DoFixFace

    ! If DoFaceOnly is present, then fix the face areas only. This should be
    ! used after the AMR is complete.
    logical, optional, intent(in):: DoFaceOnly

    real :: PositionMin_D(MaxDim), PositionMax_D(MaxDim), Coord_D(MaxDim), &
         FaceNormal_D(MaxDim)

    real, allocatable:: rCell_I(:), rFace_I(:), dCosTheta_I(:), &
         Xyz_DN(:,:,:,:)

    real :: Theta, Dphi, Dz, Dtheta
    integer :: iNode, i, j, k, Di, Dj, Dk

    real, parameter:: cThird = 1.0/3.0

    character(len=*), parameter:: NameSub = 'create_grid_block'
    !-------------------------------------------------------------------------
    if(present(iNodeIn))then
       iNode = iNodeIn
    else
       iNode = iNode_B(iBlock)
    end if

    if(.not.present(DoFaceOnly))then
       call get_tree_position(iNode, PositionMin_D, PositionMax_D)

       CoordMin_DB(:,iBlock) = CoordMin_D + (DomainSize_D)*PositionMin_D
       CoordMax_DB(:,iBlock) = CoordMin_D + (DomainSize_D)*PositionMax_D
       CellSize_DB(:,iBlock) = &
            (CoordMax_DB(:,iBlock) - CoordMin_DB(:,iBlock)) / nIjk_D

       ! The cell volumes and face areas in generalized coordinates.
       ! For Cartesian grid same as physical volume and area.
       CellVolume_B(iBlock)  = product(CellSize_DB(:,iBlock))
       CellFace_DB(:,iBlock) = CellVolume_B(iBlock) / CellSize_DB(:,iBlock)
    end if

    if(IsCartesianGrid .or. IsRotatedCartesian)then

       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
          Xyz_DGB(:,i,j,k,iBlock) = CoordMin_DB(:,iBlock) + &
               ( (/i, j, k/) - 0.5 ) * CellSize_DB(:,iBlock)
       end do; end do; end do

       do k = 1, nKNode; do j = 1, nJNode; do i = 1, nINode
          Xyz_DNB(:,i,j,k,iBlock) = CoordMin_DB(:,iBlock) + &
               ( (/i, j, k/) - 1 ) * CellSize_DB(:,iBlock)
       end do; end do; end do

       if(IsRzGeometry)then
          do j = MinJ, MaxJ
             ! NOTE: beyond the axis (y<0) the ghost cell volume is NEGATIVE!
             ! This allows the conservative prolongation work in BATL_amr.
             CellVolume_GB(:,j,:,iBlock) = &
                  CellVolume_B(iBlock)*Xyz_DGB(2,1,j,1,iBlock)
          end do
          CellVolume_GB(1:nI,1:nJ,1:nK,iBlock) = &
               abs(CellVolume_GB(1:nI,1:nJ,1:nK,iBlock))
          do j = 1, nJ
             CellFace_DFB(1,:,j,1:nK,iBlock) = &
                  CellFace_DB(1,iBlock)*abs(Xyz_DGB(2,1,j,1,iBlock))
          end do
          do j = 1, nJ+1
             ! Could use node coordinate here !!!
             CellFace_DFB(2,1:nI,j,1:nK,iBlock) = CellFace_DB(2,iBlock) &
                  *0.5*abs(sum(Xyz_DGB(2,1,j-1:j,1,iBlock)))
          end do
          CellFace_DFB(3,:,:,:,iBlock) = CellFace_DB(3,iBlock)

          FaceNormal_DDFB(1,1,:,:nJ,:nK,iBlock) = &
               CellFace_DFB(1,:,:nJ,:nK,iBlock)
          FaceNormal_DDFB(2,1,:,:nJ,:nK,iBlock) = 0.0

          FaceNormal_DDFB(1,2,:nI,:,:nK,iBlock) = 0.0
          FaceNormal_DDFB(2,2,:nI,:,:nK,iBlock) = &
               CellFace_DFB(2,:nI,:,:nK,iBlock)
       else
          ! Also useful for Cartesian to keep code simple
          CellVolume_GB(:,:,:,iBlock) = CellVolume_B(iBlock)
       end if

       if(IsRotatedCartesian)then
          ! Rotate coordinates
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             Xyz_DGB(:,i,j,k,iBlock) = cart_to_rot(Xyz_DGB(:,i,j,k,iBlock))
          end do; end do; end do

          do k = 1, nKNode; do j = 1, nJNode; do i = 1, nINode
             Xyz_DNB(:,i,j,k,iBlock) = cart_to_rot(Xyz_DNB(:,i,j,k,iBlock))
          end do; end do; end do

          ! Define face variables used by non-Cartesian grids
          CellFace_DFB(1,:,:,:,iBlock) = CellFace_DB(1,iBlock)
          CellFace_DFB(2,:,:,:,iBlock) = CellFace_DB(2,iBlock)
          CellFace_DFB(3,:,:,:,iBlock) = CellFace_DB(3,iBlock)

          FaceNormal_D = cart_to_rot((/CellFace_DB(1,iBlock), 0., 0./))
          do k = 1, nK+1; do j = 1, nJ+1; do i = 1, nI+1
             FaceNormal_DDFB(:,1,i,j,k,iBlock) = FaceNormal_D(1:nDim)
          end do; end do; end do

          FaceNormal_D = cart_to_rot((/0.0, CellFace_DB(2,iBlock), 0./))
          do k = 1, nK+1; do j = 1, nJ+1; do i = 1, nI+1
             FaceNormal_DDFB(:,2,i,j,k,iBlock) = FaceNormal_D(1:nDim)
          end do; end do; end do

          if(nDim==3)then
             FaceNormal_D = cart_to_rot((/0.,0.,CellFace_DB(3,iBlock)/))
             do k = 1, nK+1; do j = 1, nJ+1; do i = 1, nI+1
                FaceNormal_DDFB(:,3,i,j,k,iBlock) = FaceNormal_D
             end do; end do; end do
          end if
       end if
    else

       if(.not.present(DoFaceOnly))then
          ! Cell center positions based on generalized coordinates
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             Coord_D = CoordMin_DB(:,iBlock) + &
                  ( (/i, j, k/) - 0.5 ) * CellSize_DB(:,iBlock)
             call coord_to_xyz(Coord_D, Xyz_DGB(:,i,j,k,iBlock))
          end do; end do; end do
       end if

       if(IsRoundCube)then
          ! Allocate nodes even for ghost cells for volume calculation
          if(nDim == 2) allocate(Xyz_DN(3,MinI:MaxI+1,MinJ:MaxJ+1,1))
          if(nDim == 3) allocate(Xyz_DN(3,MinI:MaxI+1,MinJ:MaxJ+1,MinK:MaxK+1))

          ! Calculate node positions in Cartesian space
          do k = MinK, MaxK+1
             if(nDim == 2 .and. k /= 1) CYCLE
             do j = MinJ, MaxJ+1; do i = MinI, MaxI+1
                Coord_D = CoordMin_DB(:,iBlock) + &
                     ( (/i, j, k/) - 1 ) * CellSize_DB(:,iBlock)
                call coord_to_xyz(Coord_D, Xyz_DN(:,i,j,k) )
             end do; end do
          end do

       else
          ! Allocate and set the usual nodes
          allocate(Xyz_DN(MaxDim,nINode,nJNode,nKNode))

          if(present(DoFaceOnly) .and. .not.IsRoundCube)then
             ! Copy stored node values (for round cube we need extras)
             Xyz_DN(:,:,:,:) = Xyz_DNB(:,:,:,:,iBlock)
          else
             ! Calculate node positions in Cartesian space
             do k = 1, nKNode
                if(nDim == 2 .and. k /= 1) CYCLE
                do j = 1, nJNode; do i = 1, nINode
                   Coord_D = CoordMin_DB(:,iBlock) + &
                        ( (/i, j, k/) - 1 ) * CellSize_DB(:,iBlock)
                   call coord_to_xyz(Coord_D, Xyz_DN(:,i,j,k) )
                end do; end do
             end do
          end if
       end if

       ! Store node coordinates
       if(.not.present(DoFaceOnly)) &
            Xyz_DNB(:,:,:,:,iBlock) = Xyz_DN(:,1:nINode,1:nJNode,1:nKNode)

       ! Correct node positions at the resolution changes if required
       if(nDim==3 .and. present(DoFixFace) .or. present(DoFaceOnly))then
          call average_grid_node(iBlock, MaxDim, &
               Xyz_DN(:,1:nINode,1:nJNode,1:nKNode))
       end if

       if(IsNodeBasedGrid)then
          ! Calculate face area vectors assuming flat faces
          if(nDim == 2)then
             ! Calculate face vectors as 90 degree rotations of edge vectors
             do j = 1, nJ; do i = 1, nI+1
                FaceNormal_DDFB(x_,1,i,j,1,iBlock) = &
                     Xyz_DN(2,i,j+1,1) - Xyz_DN(2,i,j,1)
                FaceNormal_DDFB(y_,1,i,j,1,iBlock) = &
                     Xyz_DN(1,i,j,1) - Xyz_DN(1,i,j+1,1)

                CellFace_DFB(1,i,j,1,iBlock) = &
                     sqrt(sum(FaceNormal_DDFB(:,1,i,j,1,iBlock)**2))

             end do; end do
             do j = 1, nJ+1; do i = 1, nI
                FaceNormal_DDFB(x_,2,i,j,1,iBlock) = &
                     Xyz_DN(2,i,j,1) - Xyz_DN(2,i+1,j,1)
                FaceNormal_DDFB(y_,2,i,j,1,iBlock) = &
                     Xyz_DN(1,i+1,j,1) - Xyz_DN(1,i,j,1)

                CellFace_DFB(2,i,j,1,iBlock) = &
                     sqrt(sum(FaceNormal_DDFB(:,2,i,j,1,iBlock)**2))

             end do; end do
          else

             ! Calculate face area vectors as cross products of diagonals
             do k = 1, nK; do j = 1, nJ
                Di = 1
                if(present(DoFaceOnly).and.j>1.and.j<nJ.and.k>1.and.k<nK) &
                     Di=nI
                do i = 1, nI+1, Di
                   FaceNormal_DDFB(:,1,i,j,k,iBlock) = 0.5*cross_product( &
                        Xyz_DN(:,i,j+1,k+1) - Xyz_DN(:,i,j  ,k),           &
                        Xyz_DN(:,i,j  ,k+1) - Xyz_DN(:,i,j+1,k)          )

                   CellFace_DFB(1,i,j,k,iBlock) = &
                        sqrt(sum(FaceNormal_DDFB(:,1,i,j,k,iBlock)**2))

                end do
             end do; end do
             do k = 1, nK; do i = 1, nI 
                Dj = 1
                if(present(DoFaceOnly).and.i>1.and.i<nI.and.k>1.and.k<nK) &
                     Dj=nJ
                do j = 1, nJ+1, Dj
                   FaceNormal_DDFB(:,2,i,j,k,iBlock) = 0.5*cross_product( &
                        Xyz_DN(:,i+1,j,k+1) - Xyz_DN(:,i,j,k  ),           &
                        Xyz_DN(:,i+1,j,k  ) - Xyz_DN(:,i,j,k+1)          )

                   CellFace_DFB(2,i,j,k,iBlock) = &
                        sqrt(sum(FaceNormal_DDFB(:,2,i,j,k,iBlock)**2))

                end do; end do
             end do
             do j = 1, nJ; do i = 1, nI
                Dk = 1
                if(present(DoFaceOnly).and.i>1.and.i<nI.and.j>1.and.j<nJ) &
                     Dk=nK
                do k = 1, nK+1, Dk; 
                   FaceNormal_DDFB(:,3,i,j,k,iBlock) = 0.5*cross_product( &
                        Xyz_DN(:,i+1,j+1,k) - Xyz_DN(:,i  ,j,k),           &
                        Xyz_DN(:,i  ,j+1,k) - Xyz_DN(:,i+1,j,k)          )

                   CellFace_DFB(3,i,j,k,iBlock) = &
                        sqrt(sum(FaceNormal_DDFB(:,3,i,j,k,iBlock)**2))

                end do; end do
             end do
          end if ! if (iDim == 2)
       end if ! if (IsNodeBasedGrid)

       ! This is not actually 'node based'. 
       ! This is a FD correction. It is high-order and can preserve 
       ! free-stream solution. [Yan Jiang et al, Free-stream preserving
       ! finite difference schemes on curvilinear meshes]
       if(UseHighFDGeometry) &
            call correct_geometry_high_order

       ! Cell volumes for grids with no analytic formulas
       if(IsRoundCube)then
          if(nDim == 2)then
             ! Calculate cell volume as a sum of 2 triangle areas
             ! Also calculate cell center as the center of mass
             do j = MinJ, MaxJ; do i = MinI, MaxI
                CellVolume_GB(i,j,1,iBlock) = 0.5*abs(            &
                     (Xyz_DN(1,i+1,j+1,1) - Xyz_DN(1,i,j  ,1))*   &
                     (Xyz_DN(2,i+1,j  ,1) - Xyz_DN(2,i,j+1,1)) -  &
                     (Xyz_DN(2,i+1,j+1,1) - Xyz_DN(2,i,j  ,1))*   &
                     (Xyz_DN(1,i+1,j  ,1) - Xyz_DN(1,i,j+1,1)))
             end do; end do
          else
             ! Calculate cell volume as a sum of 6 tetrahedra
             ! The tips of the tetrahedra are at the min position of the cell
             do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
                CellVolume_GB(i,j,k,iBlock) = &
                     volume4(i,j,k, i+1,j  ,k  , i+1,j+1,k  , i+1,j  ,k+1) + &
                     volume4(i,j,k, i+1,j+1,k+1, i+1,j  ,k+1, i+1,j+1,k  ) + &
                     volume4(i,j,k, i  ,j+1,k  , i  ,j+1,k+1, i+1,j+1,k  ) + &
                     volume4(i,j,k, i+1,j+1,k+1, i+1,j+1,k  , i  ,j+1,k+1) + &
                     volume4(i,j,k, i  ,j  ,k+1, i+1,j  ,k+1, i  ,j+1,k+1) + &
                     volume4(i,j,k, i+1,j+1,k+1, i  ,j+1,k+1, i+1,j  ,k+1)
             end do; end do; end do
          end if
       elseif(.not.present(DoFaceOnly))then
          ! cylindrical, spherical, or rlonlat geometries

          allocate(rCell_I(MinI:MaxI), rFace_I(MinI:MaxI+1))
          do i = MinI, MaxI
             rCell_I(i) = CoordMin_DB(1,iBlock) + (i-0.5)*CellSize_DB(1,iBlock)
          end do
          do i = MinI, MaxI+1
             rFace_I(i) = CoordMin_DB(1,iBlock) + (i-1)*CellSize_DB(1,iBlock)
          end do
          if(IsLogRadius)then
             rCell_I = exp(rCell_I)
             rFace_I = exp(rFace_I)
          elseif(IsGenRadius)then
             do i = MinI, MaxI
                call gen_to_radius(rCell_I(i))
             end do
             do i = MinI, MaxI+1
                call gen_to_radius(rFace_I(i))
             end do
          end if

          Dphi = CellSize_DB(Phi_,iBlock)

          if(IsCylindrical)then

             Dz = CellSize_DB(3,iBlock)

             ! dV = r*dr*dphi*dz
             do i = MinI, MaxI
                ! NOTE: for ghost cells beyond the axis r=0 can be negative
                CellVolume_GB(i,:,:,iBlock) = &
                     rCell_I(i)*(rFace_I(i+1) - rFace_I(i))*Dphi*Dz
             end do

          elseif(IsSpherical)then

             Dtheta = CellSize_DB(Theta_,iBlock)

             allocate(dCosTheta_I(MinJ:MaxJ))

             do j = MinJ, MaxJ
                Theta = CoordMin_DB(Theta_,iBlock) + (j-1)*Dtheta
                ! Note the sign change
                dCosTheta_I(j) = cos(Theta) - cos(Theta + dTheta)
             end do

             ! dV = d(r^3)/3*dphi*d(cos theta)
             do j = MinJ, MaxJ; do i = MinI, MaxI
                ! NOTE: for ghost cells beyond the axis r=0 can be negative
                CellVolume_GB(i,j,:,iBlock) = &
                     cThird*(rFace_I(i+1)**3-rFace_I(i)**3)*Dphi*dCosTheta_I(j)
             end do; end do

          elseif(IsRLonLat)then
             Dtheta = CellSize_DB(Lat_,iBlock)

             allocate(dCosTheta_I(MinK:MaxK))

             do k = MinK, MaxK
                Theta = cHalfPi - (CoordMin_DB(Lat_,iBlock) + (k-1)*Dtheta)
                ! Note the sign change
                dCosTheta_I(k) = cos(Theta - dTheta) - cos(Theta)
             end do

             ! dV = d(r^3)/3*dphi*d(cos theta)
             do k = MinK, MaxK; do i = MinI, MaxI
                ! NOTE: for ghost cells beyond the axis r=0 can be negative
                CellVolume_GB(i,:,k,iBlock) = &
                     cThird*(rFace_I(i+1)**3-rFace_I(i)**3)*Dphi*dCosTheta_I(k)
             end do; end do

          else
             call CON_stop(NameSub//': '//TypeGeometry// &
                  ' geometry is not yet implemented')
          end if

          if(.not.IsNodeBasedGrid) call calc_analytic_face

       end if

       if(allocated(rCell_I))     deallocate(rCell_I, rFace_I)
       if(allocated(dCosTheta_I)) deallocate(dCosTheta_I)
       if(allocated(Xyz_DN))      deallocate(Xyz_DN)
    end if

  contains
    !=========================================================================
    subroutine calc_analytic_face

      real, allocatable:: SinThetaFace_I(:), &
           SinPhi_I(:), CosPhi_I(:), SinPhiFace_I(:), CosPhiFace_I(:)

      real:: Phi, Area, d_D(MaxDim)

      integer:: i, j, k

      character(len=*), parameter:: NameSub = 'calc_analytic_face'
      !----------------------------------------------------------------------

      ! Calculate analytic (curved) face areas
      if(IsCylindrical)then

         allocate( &
              SinPhi_I(MinJ:MaxJ), CosPhi_I(MinJ:MaxJ), &
              SinPhiFace_I(nJ+1), CosPhiFace_I(nJ+1) )

         do j = MinJ, MaxJ
            Phi = CoordMin_DB(2,iBlock) + (j-0.5)*Dphi
            SinPhi_I(j) =  sin(Phi)
            CosPhi_I(j) =  cos(Phi)
         end do

         do j = 1, nJ+1
            Phi = CoordMin_DB(2,iBlock) + (j-1)*Dphi
            SinPhiFace_I(j) =  sin(Phi)
            CosPhiFace_I(j) =  cos(Phi)
         end do

         ! dA_r = r_(i+1/2)*dphi*dz * (cos phi, sin phi, 0)
         if(nDim == 3) FaceNormal_DDFB(z_,r_,:,:,:,iBlock) = 0
         do i = 1, nI+1
            Area = rFace_I(i)*Dphi*Dz
            CellFace_DFB(r_,i,:,:,iBlock) = Area
            if(Area > 0)then
               do k = 1, nK; do j=1, nJ
                  FaceNormal_DDFB(x_,r_,i,j,k,iBlock) = Area*CosPhi_I(j)
                  FaceNormal_DDFB(y_,r_,i,j,k,iBlock) = Area*SinPhi_I(j)
               end do; end do
            else
               FaceNormal_DDFB(x_:y_,r_,i,:,:,iBlock) = 0
            end if
         end do

         ! dA_phi = dr*dz * (-sin phi, cos phi, 0) = dr*dz * (-y/r, x/r, 0)
         if(nDim == 3)FaceNormal_DDFB(z_,Phi_,:,:,:,iBlock) = 0
         do i = 1, nI
            Area = (rFace_I(i+1)-rFace_I(i))*Dz
            CellFace_DFB(Phi_,i,:,:,iBlock) = Area
            do k = 1, nK; do j=1, nJ+1
               FaceNormal_DDFB(x_,Phi_,i,j,k,iBlock) = -Area*SinPhiFace_I(j)
               FaceNormal_DDFB(y_,Phi_,i,j,k,iBlock) = +Area*CosPhiFace_I(j)
            end do; end do
         end do

         if(nDim == 3)then
            ! dA_z = r*dr*dphi * (0,0,1)
            FaceNormal_DDFB(x_:y_,z_,:,:,:,iBlock) = 0
            do i = 1, nI
               Area = rCell_I(i)*(rFace_I(i+1)-rFace_I(i))*Dphi
               CellFace_DFB(z_,i,:,:,iBlock)       = Area
               FaceNormal_DDFB(z_,z_,i,:,:,iBlock) = Area
            end do
         end if

         deallocate(SinPhi_I, CosPhi_I, SinPhiFace_I, CosPhiFace_I)

      elseif(IsSpherical)then

         allocate(SinThetaFace_I(nJ+1))

         do j = 1, nJ+1
            SinThetaFace_I(j) = sin(CoordMin_DB(Theta_,iBlock) + (j-1)*Dtheta)
         end do

         ! dA_r = r_(i+1/2)^2*dphi*d(cos theta)
         do k = 1, nK; do j = 1, nJ; do i = 1, nI+1
            ! Exact surface area
            Area = rFace_I(i)**2*Dphi*dCosTheta_I(j)
            CellFace_DFB(r_,i,j,k,iBlock) = Area

            ! Orthogonal coordinate system
            d_D = Xyz_DGB(:,i,j,k,iBlock) - Xyz_DGB(:,i-1,j,k,iBlock)
            FaceNormal_DDFB(:,r_,i,j,k,iBlock) = Area*d_D/sqrt(sum(d_D**2))

         end do; end do; end do

         ! dA_theta = r_i*sin(theta)*dr*dphi
         do k = 1, nK; do j=1, nJ+1; do i = 1, nI
            Area = rCell_I(i)*SinThetaFace_I(j)*(rFace_I(i+1)-rFace_I(i))*Dphi
            CellFace_DFB(Theta_,i,j,k,iBlock) = Area

            ! Orthogonal coordinate system
            d_D = Xyz_DGB(:,i,j,k,iBlock) - Xyz_DGB(:,i,j-1,k,iBlock)
            FaceNormal_DDFB(:,Theta_,i,j,k,iBlock)= Area*d_D/sqrt(sum(d_D**2))

         end do; end do; end do

         ! dA_phi = r*dr*dtheta
         do k = 1, nK+1; do j=1, nJ; do i = 1, nI
            Area = rCell_I(i)*(rFace_I(i+1)-rFace_I(i))*Dtheta
            CellFace_DFB(Phi_,i,j,k,iBlock) = Area

            ! Orthogonal coordinate system
            d_D = Xyz_DGB(:,i,j,k,iBlock) - Xyz_DGB(:,i,j,k-1,iBlock)
            FaceNormal_DDFB(:,Phi_,i,j,k,iBlock)= Area*d_D/sqrt(sum(d_D**2))

         end do; end do; end do

         deallocate(SinThetaFace_I)

      elseif(IsRLonLat)then

         allocate(SinThetaFace_I(nK+1))

         do k = 1, nK+1
            SinThetaFace_I(k) = cos(CoordMin_DB(Lat_,iBlock) + (k-1)*Dtheta)
         end do

         ! dA_r = r_(i+1/2)^2*dphi*d(cos theta)
         do k = 1, nK; do j = 1, nJ; do i = 1, nI+1
            ! Exact surface area
            Area = rFace_I(i)**2*Dphi*dCosTheta_I(k)
            CellFace_DFB(r_,i,j,k,iBlock) = Area

            ! Orthogonal coordinate system
            d_D = Xyz_DGB(:,i,j,k,iBlock) - Xyz_DGB(:,i-1,j,k,iBlock)
            FaceNormal_DDFB(:,r_,i,j,k,iBlock) = Area*d_D/sqrt(sum(d_D**2))

         end do; end do; end do

         ! dA_phi = r*dr*dtheta
         do k = 1, nK; do j=1, nJ+1; do i = 1, nI
            Area = rCell_I(i)*(rFace_I(i+1)-rFace_I(i))*Dtheta
            CellFace_DFB(Phi_,i,j,k,iBlock) = Area

            ! Orthogonal coordinate system
            d_D = Xyz_DGB(:,i,j,k,iBlock) - Xyz_DGB(:,i,j-1,k,iBlock)
            FaceNormal_DDFB(:,Phi_,i,j,k,iBlock)= Area*d_D/sqrt(sum(d_D**2))

         end do; end do; end do

         ! dA_lat = r_i*sin(theta)*dr*dphi
         do k = 1, nK+1; do j=1, nJ; do i = 1, nI
            Area = rCell_I(i)*SinThetaFace_I(k)*(rFace_I(i+1)-rFace_I(i))*Dphi
            CellFace_DFB(Lat_,i,j,k,iBlock) = Area

            ! Orthogonal coordinate system
            d_D = Xyz_DGB(:,i,j,k,iBlock) - Xyz_DGB(:,i,j,k-1,iBlock)
            FaceNormal_DDFB(:,Lat_,i,j,k,iBlock)= Area*d_D/sqrt(sum(d_D**2))

         end do; end do; end do

         deallocate(SinThetaFace_I)

      else
         call CON_stop(NameSub//': '//TypeGeometry// &
              ' geometry is not yet implemented')
      end if

    end subroutine calc_analytic_face
    !=========================================================================
    real function volume4(i1,j1,k1, i2,j2,k2, i3,j3,k3, i4,j4,k4)

      integer, intent(in):: i1,j1,k1, i2,j2,k2, i3,j3,k3, i4,j4,k4

      ! Return the volume of the tetrahedron enclosed by the 4 nodes.
      ! The volume can be negative for ghost cells.

      real, parameter:: cSixth = 1.0/6.0

      real, dimension(3):: a_D, b_D, c_D, d_D
      !----------------------------------------------------------------------
      a_D = Xyz_DN(:,i1,j1,k1)
      b_D = Xyz_DN(:,i2,j2,k2) - a_D
      c_D = Xyz_DN(:,i3,j3,k3) - a_D
      d_D = Xyz_DN(:,i4,j4,k4) - a_D

      ! Triple product divided by 6
      volume4 = cSixth*sum( b_D*cross_product(c_D, d_D) )

    end function volume4

    !==========================================================================
    subroutine correct_geometry_high_order
      ! Jiang, Yan, Chi-Wang Shu, and Mengping Zhang. "Free-stream preserving o
      ! finite difference schemes on curvilinear meshes." Brown University, 
      ! Scientific Computing Group, Report 10 (2013): 2013.

      call calc_metrics(iBlock)
      call coef_cart_to_noncart(iBlock)
      call calc_face_normal(iBlock)

    end subroutine correct_geometry_high_order

  end subroutine create_grid_block

  !===========================================================================

  subroutine average_grid_node(iBlock, nVar, Var_VN)

    integer, intent(in)   :: iBlock
    integer, intent(in)   :: nVar
    real,    intent(inout):: Var_VN(nVar,nINode,nJNode,nKNode)

    ! Move nodes on the fine side of resolution change to the plane
    ! defined by the coarse side. This ensures that the sum of the
    ! faces form closed surfaces, so that a uniform flow is preserved.

    ! This routine does the same as the Tecplot node fix in BATS-R-US.

    integer :: i, j, k
    integer :: i1, i2, j1, j2, k1, k2, Di, Dj, Dk
    integer :: iDir, jDir, kDir, nDir
    !----------------------------------------------------------------------
    ! Loop over neighbor directions and set index ranges
    do kDir = -1,1
       select case(kDir)
       case( 1)
          k1=1+nK; k2=1+nK; Dk=0
       case(-1)
          k1=1;    k2=1;    Dk=0
       case( 0)
          k1=2;    k2=nK;   Dk=1
       end select
       do jDir = -1,1
          select case(jDir)
          case( 1)
             j1=1+nJ; j2=1+nJ; Dj=0
          case(-1)
             j1=1;    j2=1;    Dj=0
          case( 0)
             j1=2;    j2=nJ;   Dj=1
          end select
          do iDir = -1,1

             ! Check number of non-zero directions (1:face, 2:edge, 3:corner)
             nDir = abs(iDir) + abs(jDir) + abs(kDir)

             ! ignore corners
             if(nDir == 3) CYCLE

             ! Check if there is any coarser neighbor
             if(DiLevelNei_IIIB(iDir,jDir,kDir,iBlock) /= 1) CYCLE

             select case(iDir)
             case( 1)
                i1=1+nI; i2=1+nI; Di=0
             case(-1)
                i1=1;    i2=1;    Di=0
             case( 0)
                i1=2;    i2=nI;   Di=1
             end select

             ! Correct edge nodes and some interior face nodes
             do k=k1,k2,2; do j=j1,j2,2; do i=i1,i2,2
                Var_VN(:,i,j,k) = 0.125 * ( &
                     Var_VN(:,i-Di,j-Dj,k-Dk) + &
                     Var_VN(:,i-Di,j-Dj,k+Dk) + &
                     Var_VN(:,i-Di,j+Dj,k-Dk) + &
                     Var_VN(:,i-Di,j+Dj,k+Dk) + &
                     Var_VN(:,i+Di,j-Dj,k-Dk) + &
                     Var_VN(:,i+Di,j-Dj,k+Dk) + &
                     Var_VN(:,i+Di,j+Dj,k-Dk) + &
                     Var_VN(:,i+Di,j+Dj,k+Dk) )
             end do; end do; end do

             ! Done with edge neighbors
             if(nDir == 2) CYCLE

             ! Add correction of additional interior face nodes
             if(Di==1)then
                do k=k1,k2,2; do j=j1,j2,2; do i=i1-1,i2+1,2
                   Var_VN(:,i,j,k) = 0.25 * ( &
                        Var_VN(:,i,j-Dj,k-Dk) + &
                        Var_VN(:,i,j-Dj,k+Dk) + &
                        Var_VN(:,i,j+Dj,k-Dk) + &
                        Var_VN(:,i,j+Dj,k+Dk) )
                end do; end do; end do
             end if
             if(Dj==1)then
                do k=k1,k2,2; do j=j1-1,j2+1,2; do i=i1,i2,2 
                   Var_VN(:,i,j,k) = 0.25 * ( &
                        Var_VN(:,i-Di,j,k-Dk) + &
                        Var_VN(:,i-Di,j,k+Dk) + &
                        Var_VN(:,i+Di,j,k-Dk) + &
                        Var_VN(:,i+Di,j,k+Dk) )
                end do; end do; end do
             end if
             if(Dk==1)then
                do k=k1-1,k2+1,2; do j=j1,j2,2; do i=i1,i2,2
                   Var_VN(:,i,j,k) = 0.25 * ( &
                        Var_VN(:,i-Di,j-Dj,k) + &
                        Var_VN(:,i-Di,j+Dj,k) + &
                        Var_VN(:,i+Di,j-Dj,k) + &
                        Var_VN(:,i+Di,j+Dj,k) )
                end do; end do; end do
             end if

          end do
       end do
    end do

  end subroutine average_grid_node

  !===========================================================================
  subroutine fix_grid_res_change

    integer:: iBlock, iDir, jDir, kDir
    !------------------------------------------------------------------------
    LOOPBLOCK: do iBlock = 1, nBlock
       if(Unused_B(iBlock))CYCLE

       if(iAmrChange_B(iBlock) == AmrNeiChanged_)then
          call create_grid_block(iBlock, DoFaceOnly=.true.)
       elseif(iAmrChange_B(iBlock) >= AmrMoved_)then
          do kDir = -1,1; do jDir = -1,1; do iDir = -1,1
             ! ignore corners
             if(abs(iDir) + abs(jDir) + abs(kDir) == 3) CYCLE
             if(DiLevelNei_IIIB(iDir,jDir,kDir,iBlock) == 1)then
                call create_grid_block(iBlock, DoFaceOnly=.true.)
                CYCLE LOOPBLOCK
             end if
          end do; end do; end do
       end if
    end do LOOPBLOCK

  end subroutine fix_grid_res_change
  !===========================================================================

  subroutine create_grid

    integer:: iBlock
    !------------------------------------------------------------------------
    if(nDim == 3 .and. IsNodeBasedGrid .and. .not. IsCartesianGrid)then
       do iBlock = 1, nBlock
          if(Unused_B(iBlock))CYCLE
          call create_grid_block(iBlock, DoFixFace=.true.)
       end do
    else
       do iBlock = 1, nBlock
          if(Unused_B(iBlock))CYCLE
          call create_grid_block(iBlock)
       end do
    end if

  end subroutine create_grid
  !===========================================================================

  subroutine show_grid_block(iBlock)

    integer, intent(in):: iBlock

    ! Show grid information for block iBlock

    integer:: iDim

    character(len=*), parameter:: NameSub = 'show_grid_block'
    !------------------------------------------------------------------------
    if(Unused_B(iBlock))then
       write(*,*) NameSub//' WARNING unused block ',iBlock,' on proc',iProc
       RETURN
    end if
    write(*,*)'show_grid_block for iProc, iBlock=',iProc, iBlock
    write(*,*)'CoordMin  =', CoordMin_DB(:,iBlock)
    write(*,*)'CoordMax  =', CoordMax_DB(:,iBlock)
    write(*,*)'CellSize  =', CellSize_DB(:,iBlock)
    if(IsCartesian)then
       write(*,*)'CellFace  =', CellFace_DB(:,iBlock)
       write(*,*)'CellVolume=', CellVolume_B(iBlock)
    else
       write(*,*)'CellFace(1, 1, 1)  =', CellFace_DFB(1:nDim,1,1,1,iBlock)
       write(*,*)'CellVolume(1, 1, 1)=', CellVolume_GB(1,1,1,iBlock)
       if(.not.IsRzGeometry)then
          do iDim = 1, nDim
             write(*,*)'iDim, FaceNormal_DDFB(:,iDim)=', iDim, &
                  FaceNormal_DDFB(:,iDim,1,1,1,iBlock)
          end do
       end if
    end if
    write(*,*)'Xyz( 1, 1, 1)=', Xyz_DGB(:, 1, 1, 1,iBlock)
    write(*,*)'Xyz(nI, 1, 1)=', Xyz_DGB(:,nI, 1, 1,iBlock)
    write(*,*)'Xyz( 1,nJ, 1)=', Xyz_DGB(:, 1,nJ, 1,iBlock)
    write(*,*)'Xyz( 1, 1,nK)=', Xyz_DGB(:, 1, 1,nK,iBlock)
    write(*,*)'Xyz(nI,nJ,nK)=', Xyz_DGB(:,nI,nJ,nK,iBlock)

  end subroutine show_grid_block

  !===========================================================================

  subroutine show_grid_proc

    ! Show all blocks sequentially on the calling processor

    integer:: iBlock
    !------------------------------------------------------------------------

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       call show_grid_block(iBlock)
    end do

  end subroutine show_grid_proc

  !===========================================================================

  subroutine show_grid

    use BATL_mpi, ONLY: iProc, nProc, barrier_mpi

    ! Show all blocks sequentially on all processors, ie. show_grid 
    ! must be called from all processors of the MPI communicator iComm!

    integer:: iPe
    !------------------------------------------------------------------------

    call barrier_mpi
    do iPe = 0, nProc - 1
       if(iPe == iProc) call show_grid_proc
       call barrier_mpi
    end do

  end subroutine show_grid

  !===========================================================================

  subroutine find_grid_block(XyzIn_D, &
       iProcOut, iBlockOut, iCellOut_D, DistOut_D, iNodeOut, &
       CoordMinBlockOut_D, CoordMaxBlockOut_D, CellSizeOut_D, & 
       UseGhostCell)

    ! Find the processor and block containing location XyzIn_D. 
    ! If iCellOut_D is present and UseGhostCell is not present or false,
    !   then iCell_D returns the indexes of the closest cell center:
    !   1 <= iCellOut_D <= nIjk_D
    ! If iCellOut_D is present and UseGhostCell is present and true, 
    !   then iCell_D returns the indexes of the cell to the left from XyzIn_D
    ! If present, DistOut_D returns the signed distance to the cell center
    !    given by iCell_D divided by the cell size.
    ! DistOut_D can only be present if iCellOut_D is also present.

    real,    intent(in) :: XyzIn_D(MaxDim)        ! Cartesian coords of point
    integer, intent(out):: iBlockOut, iProcOut    ! Block and proc indexes
    integer, intent(out), optional:: iCellOut_D(MaxDim) ! Closest cell indexes
    real,    intent(out), optional:: DistOut_D(MaxDim)  ! Normalized distance
    integer, intent(out), optional:: iNodeOut     ! Tree node index
    real,    intent(out), optional:: CoordMinBlockOut_D(MaxDim)! block corner
    real,    intent(out), optional:: CoordMaxBlockOut_D(MaxDim)! block corner 
    real,    intent(out), optional:: CellSizeOut_D(MaxDim) ! cell size in block
    logical, intent(in),  optional:: UseGhostCell ! Use ghost cells or not

    real:: CoordTree_D(MaxDim), Coord_D(MaxDim)
    real:: PositionMin_D(MaxDim), PositionMax_D(MaxDim)
    integer:: iNode

    character(len=*), parameter:: NameSub = 'find_grid_block'
    !------------------------------------------------------------------------
    ! Convert to generalized coordinates if necessary
    if(IsCartesianGrid)then
       Coord_D = XyzIn_D
    else
       call xyz_to_coord(XyzIn_D, Coord_D)
    end if
    ! Calculate normalized coordinates for tree search
    CoordTree_D = (Coord_D - CoordMin_D)/(DomainSize_D)

    if(any(CoordTree_D < 0.0) .or. any(CoordTree_D > 1.0))then
       iBlockOut = Unset_
       iProcOut  = Unset_
       if(present(iNodeOut)) iNodeOut = Unset_
       RETURN
    end if

    ! Find node containing the point
    if(present(iCellOut_D))then
       call find_tree_cell(CoordTree_D, iNode, iCellOut_D, DistOut_D, &
            UseGhostCell)
    else
       call find_tree_node(CoordTree_D, iNode)
    end if

    ! Check if point was found
    if(iNode > 0)then
       ! Convert to block and processor indexes
       iBlockOut = iTree_IA(Block_,iNode)
       iProcOut  = iTree_IA(Proc_, iNode)
    else
       iBlockOut = Unset_
       iProcOut  = Unset_
    end if

    if(present(iNodeOut)) iNodeOut = iNode

    if(  present(CoordMinBlockOut_D) .or. &
         present(CoordMaxBlockOut_D) .or. &
         present(CellSizeOut_D))then

       call get_tree_position(iNode, PositionMin_D, PositionMax_D)

       if(present(CoordMinBlockOut_D))&
            CoordMinBlockOut_D = CoordMin_D + PositionMin_D*DomainSize_D
       if(present(CoordMaxBlockOut_D))&
            CoordMaxBlockOut_D = CoordMin_D + PositionMax_D*DomainSize_D
       if(present(CellSizeOut_D))&
            CellSizeOut_D = (PositionMax_D-PositionMin_D)*DomainSize_D/nIjk_D

    end if

  end subroutine find_grid_block

  !===========================================================================
  subroutine interpolate_grid(Xyz_D, nCell, iCell_II, Weight_I)

    ! Find the grid cells surrounding the point Xyz_D.
    ! nCell returns the number of cells found on the processor.
    ! iCell_II returns the block+cell indexes for each cell.
    ! Weight_I returns the interpolation weights

    real,    intent(in) :: Xyz_D(MaxDim)
    integer, intent(out):: nCell  
    integer, intent(out):: iCell_II(0:nDim,2**nDim)
    real,    intent(out):: Weight_I(2**nDim)

    real:: Coord_D(MaxDim), CoordOrig_D(MaxDim), CoordCell_D(nDim)
    real:: CoordMin, CoordMax, Shift
    real:: CellSize_D(nDim), BufferLo_D(nDim), BufferHi_D(nDim)
    real:: InvSize_D(nDim), Weight_D(nDim), Weight
    integer:: iBlock, iDim, DiLevel, iCell_D(MaxDim), iCell, i_D(MaxDim)
    integer:: i, j, k, iLo, jLo, kLo, iHi, jHi, kHi

    logical, parameter:: DoTest = .false.
    character(len=*), parameter:: NameSub='BATL_grid::interpolate_grid'
    !------------------------------------------------------------------------
    ! Convert to generalized coordinates if necessary
    if(IsCartesianGrid)then
       Coord_D = Xyz_D
    else
       call xyz_to_coord(Xyz_D, Coord_D)
    end if

    ! Initialize values in case no cells are found on this processor
    nCell    = 0
    iCell_II = 0
    Weight_I = 0

    ! For periodic boundaries we may have to shift the coordinates
    ! Since this shift varies from block to block, store the original
    if(any(IsPeriodic_D)) CoordOrig_D = Coord_D


    LOOPBLOCK: do iBlock = 1, nBlock

       if(Unused_B(iBlock)) CYCLE

       CellSize_D = CellSize_DB(1:nDim,iBlock)

       ! Initialize with unshifted coordinates
       if(any(IsPeriodic_D)) Coord_D = CoordOrig_D

       ! Set buffer zone according to relative size of neighboring block
       do iDim = 1, nDim

          ! Block at the lower index side
          select case(iDim)
          case(1)
             DiLevel = DiLevelNei_IIIB(-1,0,0,iBlock)
          case(2)
             DiLevel = DiLevelNei_IIIB(0,-1,0,iBlock)
          case(3)
             DiLevel = DiLevelNei_IIIB(0,0,-1,iBlock)
          end select

          select case(DiLevel)
          case(1)
             BufferLo_D(iDim) = 0.5*CellSize_D(iDim)*iRatio_D(iDim)
          case(-1)
             BufferLo_D(iDim) = 0.5*CellSize_D(iDim)/iRatio_D(iDim)
          case default
             BufferLo_D(iDim) = 0.5*CellSize_D(iDim)
          end select

          ! Lower limit on point coordinate
          CoordMin = CoordMin_DB(iDim,iBlock) - BufferLo_D(iDim)

          ! Check if point is inside the buffer zone on the lower side
          if(Coord_D(iDim) < CoordMin .and. .not.IsPeriodic_D(iDim)) &
               CYCLE LOOPBLOCK

          ! Block at the upper index side
          select case(iDim)
          case(1)
             DiLevel = DiLevelNei_IIIB(+1,0,0,iBlock)
          case(2)
             DiLevel = DiLevelNei_IIIB(0,+1,0,iBlock)
          case(3)
             DiLevel = DiLevelNei_IIIB(0,0,+1,iBlock)
          end select

          select case(DiLevel)
          case(1)
             BufferHi_D(iDim) = 0.5*CellSize_D(iDim)*iRatio_D(iDim)
          case(-1)
             BufferHi_D(iDim) = 0.5*CellSize_D(iDim)/iRatio_D(iDim)
          case default
             BufferHi_D(iDim) = 0.5*CellSize_D(iDim)
          end select

          ! Upper limit on point coordinate
          CoordMax = CoordMax_DB(iDim,iBlock) + BufferHi_D(iDim)

          ! Check if point is inside the buffer zone on the lower side
          if(Coord_D(iDim) > CoordMax .and. .not.IsPeriodic_D(iDim)) &
               CYCLE LOOPBLOCK

          ! Done unless periodic
          if(.not.IsPeriodic_D(iDim)) CYCLE

          ! Shift coordinate and try check again
          Shift = DomainSize_D(iDim)
          if(Coord_D(iDim) < CoordMin) Coord_D(iDim) = Coord_D(iDim) + Shift
          if(Coord_D(iDim) > CoordMax) Coord_D(iDim) = Coord_D(iDim) - Shift
          if(Coord_D(iDim) < CoordMin .or. Coord_D(iDim) > CoordMax) &
               CYCLE LOOPBLOCK

       end do

       ! The coarse block SHOULD NOT contribute to the interpolation 
       ! if the point is inside the finer edge/corner neighbor cell centers
       if(nDimAmr > 1)then
          i_D = 0
          do iDim = 1, nDim
             ! No need to check non-AMR directions
             if(iRatio_D(iDim) < 2) CYCLE

             ! Check which finer edge/corner neighbor the point could belong to
             if(Coord_D(iDim) < &
                  CoordMin_DB(iDim,iBlock) - 0.25*CellSize_D(iDim)) then
                i_D(iDim) = -1
             elseif(Coord_D(iDim) > &
                  CoordMax_DB(iDim,iBlock) + 0.25*CellSize_D(iDim)) then
                i_D(iDim) = 1
             end if
          end do
          ! Check 1-2 edges
          if(iRatio_D(1) > 1 .and. iRatio_D(2) > 1) then
             if(abs(i_D(1)) + abs(i_D(2)) == 2 .and. &
                  DiLevelNei_IIIB(i_D(1),i_D(2),0,iBlock) == -1) &
                  CYCLE LOOPBLOCK
          end if
          ! Check 1-3 edges
          if(iRatio_D(1) > 1 .and. iRatio_D(3) > 1) then
             if(abs(i_D(1)) + abs(i_D(3)) == 2 .and. &
                  DiLevelNei_IIIB(i_D(1),0,i_D(3),iBlock) == -1) &
                  CYCLE LOOPBLOCK
          end if
          ! Check 2-3 edges
          if(iRatio_D(2) > 1 .and. iRatio_D(3) > 1) then
             if(abs(i_D(2)) + abs(i_D(3)) == 2 .and. &
                  DiLevelNei_IIIB(0,i_D(2),i_D(3),iBlock) == -1) &
                  CYCLE LOOPBLOCK
          end if
          ! Check corners
          if(nDimAmr == 3)then
             if(abs(i_D(1)) + abs(i_D(2)) + abs(i_D(3)) == 3 .and.&
                  DiLevelNei_IIIB(i_D(1),i_D(2),i_D(3),iBlock) == -1) &
                  CYCLE LOOPBLOCK
          end if
       end if

       ! Find closest cell center indexes towards the lower index direction
       iCell_D = 1 ! for ignored dimensions
       iCell_D(1:nDim) = floor(0.5 + &
            (Coord_D(1:nDim) - CoordMin_DB(1:nDim,iBlock))/CellSize_D )

       ! Set inverse distance between interpolation cells
       do iDim = 1, nDim
          if(iCell_D(iDim) < 1)then
             iCell_D(iDim) = 0
             InvSize_D(iDim) = 1/(0.5*CellSize_D(iDim) + BufferLo_D(iDim))
          elseif(iCell_D(iDim) >= nIjk_D(iDim))then
             iCell_D(iDim) = nIjk_D(iDim)
             InvSize_D(iDim) = 1/(0.5*CellSize_D(iDim) + BufferHi_D(iDim))
          else
             InvSize_D(iDim) = 1/CellSize_D(iDim)
          end if
       end do

       ! Index range of surrounding cell centers
       iLo = iCell_D(1); iHi = iLo + 1
       jLo = iCell_D(2); jHi = jLo + 1
       kLo = iCell_D(3); kHi = kLo + 1

       ! Index range limited to physical and outer boundary ghost cells
       if(            DiLevelNei_IIIB(-1,0,0,iBlock)/=Unset_ .or. nG==0) &
            iLo = max(1,iLo)
       if(nDim<2 .or. DiLevelNei_IIIB(0,-1,0,iBlock)/=Unset_ .or. nG==0) &
            jLo = max(1,jLo)
       if(nDim<3 .or. DiLevelNei_IIIB(0,0,-1,iBlock)/=Unset_ .or. nG==0) &
            kLo = max(1,kLo)
       if(            DiLevelNei_IIIB(+1,0,0,iBlock)/=Unset_ .or. nG==0) &
            iHi = min(nI,iHi)
       if(nDim<2 .or. DiLevelNei_IIIB(0,+1,0,iBlock)/=Unset_ .or. nG==0) &
            jHi = min(nJ,jHi)
       if(nDim<3 .or. DiLevelNei_IIIB(0,0,+1,iBlock)/=Unset_ .or. nG==0) &
            kHi = min(nK,kHi)

       if(DoTest)then
          write(*,*)'!!! iProc, iBlock, iCell_D, Coord_D=', &
               iProc, iBlock, iCell_D, Coord_D(1:nDim)
          write(*,*)'!!! iLo, iHi, jLo, jHi=',iLo, iHi, jLo, jHi
       end if

       ! Calculate the weights and store it together with index information
       do k = kLo, kHi; do j = jLo, jHi; do i = iLo, iHi
          iCell_D = (/ i, j, k /)

          CoordCell_D = CoordMin_DB(1:nDim,iBlock) &
               + ( iCell_D(1:nDim) - 0.5 )*CellSize_D

          Weight_D = 1 - InvSize_D*abs((Coord_D(1:nDim) - CoordCell_D))

          Weight = product(Weight_D)

          ! Ignore cells with 0 weight
          if(Weight <= 0.0) CYCLE

          nCell = nCell + 1

          if(nCell > 2**nDim)then
             if(nCell == 2**nDim+1)then
                write(*,*)'ERROR in ',NameSub,': too many cells!'
                write(*,*)'iProc,iBlock,nDim=',iProc,iBlock,nDim
                write(*,*)'iLo,iHi,jLo,jHi,kLo,kHi=',iLo,iHi,jLo,jHi,kLo,kHi
                write(*,*)'Coord_D=',Coord_D(1:nDim)
                i_D = 1
                do iCell = 1, 2**nDim
                   i_D(1:nDim) = iCell_II(1:nDim,iCell)
                   write(*,*)'iCell, iCell_II(:,i), weight, Xyz=', &
                        iCell, iCell_II(:,iCell), Weight_I(iCell), &
                        Xyz_DGB(1:nDim,i_D(1),i_D(2),i_D(3),iCell_II(0,iCell))
                end do
             end if
             write(*,*)'iCell, iCell_II(:,i), weight=', &
                  nCell, iBlock, iCell_D(1:nDim), Weight, &
                  Xyz_DGB(1:nDim,i,j,k,iBlock)
          else
             Weight_I(nCell)   = Weight
             iCell_II(0,nCell) = iBlock
             iCell_II(1:nDim,nCell) = iCell_D(1:nDim)
          end if
          if(DoTest)write(*,*)'!!! nCell, CoordCell, Weight=', &
               nCell, CoordCell_D(1:nDim), Weight

       end do; end do; end do
    end do LOOPBLOCK

    if(nCell > 2**nDim) call CON_stop(NameSub// &
         ': too many cells to interpolate from')

  end subroutine interpolate_grid

  !===========================================================================

  subroutine interpolate_grid_amr(XyzIn_D, nCell, iCell_II, Weight_I, &
       IsSecondOrder)

    use BATL_interpolate_amr, ONLY:interpolate_amr

    ! Find the grid cells surrounding the point Xyz_D.
    ! nCell returns the number of cells found on the processor.
    ! iCell_II returns the block+cell indexes for each cell.
    ! Weight_I returns the interpolation weights calculated 
    !                                 using AMR interpolateion procedure
    real,    intent(in) :: XyzIn_D(MaxDim)
    integer, intent(out):: nCell
    integer, intent(out):: iCell_II(0:nDim,2**nDim)
    real,    intent(out):: Weight_I(2**nDim)

    logical, optional, intent(out):: IsSecondOrder

    real   :: Coord_D(MaxDim)
    logical:: IsSecondOrderLocal
    !-----------------------------------
    ! check number of AMR dimensions:
    ! if it is 0 or 1 => call a simpler interpolation function
    if(nDimAmr <= 1)then
       call interpolate_grid(XyzIn_D, nCell, iCell_II, Weight_I)
       RETURN
    end if
    
    ! Convert to generalized coordinates if necessary
    if(IsCartesianGrid)then
       Coord_D = XyzIn_D
    else
       call xyz_to_coord(XyzIn_D, Coord_D)
    end if

    ! call the wrapper for the shared AMR interpolation procedure,
    call interpolate_amr(Coord_D, &
         nCell, iCell_II, Weight_I, IsSecondOrderLocal)

    if(present(IsSecondOrder)) IsSecondOrder = IsSecondOrderLocal

  end subroutine interpolate_grid_amr

  !===========================================================================

  subroutine interpolate_grid_amr_gc(XyzIn_D, nCell, iCell_II, Weight_I, &
       IsSecondOrder)

    use BATL_interpolate_amr, &
         ONLY:interpolate_amr_gc, find_block_to_interpolate_gc

    ! Find the grid cells surrounding the point Xyz_D.
    ! nCell returns the number of cells found on the processor.
    ! iCell_II returns the block+cell indexes for each cell.
    ! Weight_I returns the interpolation weights calculated 
    !                                 using AMR interpolateion procedure
    ! Interpolation is performed using cells (including ghost) of single block
    real,    intent(in) :: XyzIn_D(MaxDim)
    integer, intent(out):: nCell
    integer, intent(out):: iCell_II(0:nDim,2**nDim)
    real,    intent(out):: Weight_I(2**nDim)

    logical, optional, intent(out):: IsSecondOrder

    real   :: Coord_D(MaxDim), DCoord_D(MaxDim), CoordMin_D(MaxDim)
    logical:: IsSecondOrderLocal
    integer:: iBlockOut, iProcOut
    !-----------------------------------
    ! check number of AMR dimensions:
    ! if it is 0 or 1 => call a simpler interpolation function
    if(nDimAmr <= 1)then
       call interpolate_grid(XyzIn_D, nCell, iCell_II, Weight_I)
       RETURN
    end if
    
    ! Convert to generalized coordinates if necessary
    if(IsCartesianGrid)then
       Coord_D = XyzIn_D
    else
       call xyz_to_coord(XyzIn_D, Coord_D)
    end if

    ! find a block suitable for interpolation
    call find_block_to_interpolate_gc(Coord_D, iProcOut, iBlockOut)
    
    ! check if it is on the current processor
    if (iProcOut /= iProc)then
       !call CON_stop("Can't perform interpolation on this processor")
       nCell = 0
       RETURN
    end if

    ! get corner coordinates and cel size of the block
    CoordMin_D =  CoordMin_DB(:,iBlockOut)
    DCoord_D   = (CoordMax_DB(:,iBlockOut) - CoordMin_DB(:,iBlockOut)) / nIJK_D

    ! call the wrapper for the shared AMR interpolation procedure
    call interpolate_amr_gc(Coord_D, CoordMin_D, DCoord_D, &
         DiLevelNei_IIIB(:, :, :, iBlockOut), &
         nCell, iCell_II, Weight_I, IsSecondOrderLocal)

    ! return block number as well
    iCell_II(0,:) = iBlockOut

    if(present(IsSecondOrder)) IsSecondOrder = IsSecondOrderLocal

  end subroutine interpolate_grid_amr_gc

  !==================================

  subroutine check_interpolate_amr_gc(Xyz_D, iBlock, IsPossible, iPeOut,& 
       iBlockOut, IsBoundary)
    !\
    ! Checks if a point with the Cartesian coordinates, Xyz_D,
    ! can be interpolated using the data in the block iBlock at the
    ! given processor, with the ghost cell values included, if needed
    !/
    use BATL_interpolate_amr, ONLY:find_block_to_interpolate_gc
    !\
    ! INPUTS:
    !/
    real,    intent(in) :: Xyz_D(MaxDim) !Point coordinates
    integer, intent(in) :: iBlock        !# of block
    !\
    ! OUTPUTS:
    !/
    logical, intent(out):: IsPossible    !.true. if interpolation is possible  
    integer, intent(out):: iPeOut, iBlockOut !else, Pe and Block to be used
    logical, intent(out):: IsBoundary    !.true. for the point out of domain
    !\
    !Generalized Coordinates
    !/
    real    :: Coord_D(MaxDim) 
    !\
    !Direction along which the point goes out of the block inner part
    !/
    integer :: iShift_D(MaxDim) 
    !-----------------------------------   
    ! Convert to generalized coordinates if necessary
    if(IsCartesianGrid)then
       Coord_D = Xyz_D
    else
       call xyz_to_coord(Xyz_D, Coord_D)
    end if
   
    !\
    ! For periodic boundary conditions fix the input coordinate if
    ! beyond the tree bounadaries
    !/
    where(IsPeriodic_D(1:nDim)) Coord_D(1:nDim) = CoordMin_D(1:nDim) + &
         modulo(Coord_D(1:nDim) - CoordMin_D(1:nDim), DomainSize_D(1:nDim))

    !\
    ! Figure out if the point goes out of the computational domain
    !/
    IsBoundary = any(Coord_D(1:nDim) < CoordMin_D(1:nDim)&
         .or. Coord_D(1:nDim) >= CoordMax_D(nDim))
    if(IsBoundary)then
       IsPossible = .false.
       iPeOut = Unset_; iBlockOut = Unset_
       RETURN
    end if
    if(Unused_B(iBlock))then
       !\
       !This may happen if the block is moved to another PE in the 
       !course of load balancing (while particles are not moved) 
       IsPossible = .false.
    else
       !\
       ! Check if the block is suitable to interpolate with ghost cells
       !/
       iShift_D(1:nDim) = floor((Coord_D(1:nDim) - CoordMin_DB(1:nDim,iBlock)&
            - 0.50*CellSize_DB(1:nDim,iBlock))/(CellSize_DB(1:nDim,iBlock) &
            *(nIJK_D(1:nDim) - 1)))
       if(all(iShift_D(1:nDim)==0))then
          !The point falls into an inner part of block
          IsPossible = .true.
       elseif(any(DiLevelNei_IIIB(min(0,iShift_D(1)):max(0,iShift_D(1)),&
                                  min(0,iShift_D(2)):max(0,iShift_D(2)),&
                                  min(0,iShift_D(3)):max(0,iShift_D(3)),&
                                  iBlock)==1))then
          !\
          ! The required information in the ghost cells is missing.
          !/ 
          IsPossible = .false.
       elseif(any(DiLevelNei_IIIB(min(0,iShift_D(1)):max(0,iShift_D(1)),&
                                  min(0,iShift_D(2)):max(0,iShift_D(2)),&
                                  min(0,iShift_D(3)):max(0,iShift_D(3)),&
                                  iBlock)==-1))then
          !\
          ! We need to interpolate using ghostcells even if the point is 
          ! inside the first layer of ghostcells
          !/
          IsPossible = all(&
               Coord_D(1:nDim) >= CoordMin_DB(1:nDim,iBlock)    &
               -CellSize_DB(1:nDim,iBlock) .and.&
               Coord_D(1:nDim) <  CoordMax_DB(1:nDim,iBlock)    &
               +CellSize_DB(1:nDim,iBlock))
       else
          !\
          ! We interpolate using ghostcells only if the point is 
          ! inside the physical block
          !/
          IsPossible = all(&
               Coord_D(1:nDim) >= CoordMin_DB(1:nDim,iBlock).and.&
               Coord_D(1:nDim) <  CoordMax_DB(1:nDim,iBlock))
       end if
    end if
    if(IsPossible)then
       iPeOut = Unset_; iBlockOut = Unset_
       RETURN
    end if
    !\
    ! call the wrapper for the shared interpolation procedure
    !/
    call find_block_to_interpolate_gc(Coord_D, iPeOut, iBlockOut)
  end subroutine check_interpolate_amr_gc
  !===========================================================================
  subroutine calc_face_normal(iBlock)
    ! Interpolate dx3/dx1 to the face, where x3=hat(Xi,Eta,Zeta), x1=x,y,z.
    integer, intent(in):: iBlock
    integer:: iFace, jFace, kFace
    integer:: iDimCart
    integer:: nIFace, nJFace, nKFace
    real   :: Area
    !----------------------------------------------------------------------

    nIFace = nI + 1; nJFace = nJ + 1; nKFace = nK + 1
    do kFace = 1, nK; do jFace = 1, nJ; do iFace = 1, nIFace
       Area = CellFace_DFB(1,iFace,jFace,kFace,iBlock)
       if(Area < 1e-15) then
          ! For singular point, the face area should be zero. 
          ! Now, the area calculated with a flat plane assumation is used
          ! to determine weather it is a singular point, and FaceNormal_DDFB
          ! and CellFace_DFB are calculated twice. The high order face area 
          ! in this subroutine can be used for determination, then only one
          ! calculation is needed, but the tolerance (now it is 1e-15) may 
          ! need to change. 
          FaceNormal_DDFB(:,1,iFace,jFace,kFace,iBlock) = 0
          CYCLE
       endif

       do iDimCart = 1, nDim
          FaceNormal_DDFB(iDimCart,1,iFace,jFace,kFace,iBlock) = &
               calc_face_value(CellCoef_DDGB(Xi_,iDimCart, &
               iFace-3:iFace+2,jFace,kFace,iBlock))               
       enddo

       CellFace_DFB(1,iFace,jFace,kFace,iBlock) = &
            sqrt(sum(FaceNormal_DDFB(:,1,iFace,jFace,kFace,iBlock)**2))
    enddo; enddo; enddo

    do kFace = 1, nK; do jFace = 1, nJFace; do iFace = 1, nI 
       Area = CellFace_DFB(2,iFace,jFace,kFace,iBlock)
       if(Area < 1e-15) then
          FaceNormal_DDFB(:,2,iFace,jFace,kFace,iBlock) = 0
          CYCLE
       endif

       do iDimCart = 1, nDim
          FaceNormal_DDFB(iDimCart,2,iFace,jFace,kFace,iBlock) = &
               calc_face_value(CellCoef_DDGB(Eta_,iDimCart, &
               iFace,jFace-3:jFace+2,kFace,iBlock))          
       enddo

       CellFace_DFB(2,iFace,jFace,kFace,iBlock) = &
            sqrt(sum(FaceNormal_DDFB(:,2,iFace,jFace,kFace,iBlock)**2))
    enddo; enddo; enddo

    if(nK > 1) then
       do kFace = 1, nKFace; do jFace = 1, nJ; do iFace = 1, nI
          Area = CellFace_DFB(3,iFace,jFace,kFace,iBlock)
          if(Area < 1e-15) then
             FaceNormal_DDFB(:,3,iFace,jFace,kFace,iBlock) = 0
             CYCLE
          endif

          do iDimCart = 1, nDim
             FaceNormal_DDFB(iDimCart,3,iFace,jFace,kFace,iBlock) = &
                  calc_face_value(CellCoef_DDGB(Zeta_, iDimCart, &
                  iFace,jFace,kFace-3:kFace+2,iBlock))
          enddo

          CellFace_DFB(3,iFace,jFace,kFace,iBlock) = &
               sqrt(sum(FaceNormal_DDFB(:,3,iFace,jFace,kFace,iBlock)**2))
       enddo; enddo; enddo
    endif

  end subroutine calc_face_normal
  !===========================================================================
  subroutine coef_cart_to_noncart(iBlock)
    ! Eq (26).
    ! Calc dx3/dx1 at cell center, where x3=hat(Xi,Eta,Zeta), x1=x,y,z.

    integer, intent(in):: iBlock
    integer:: iDimCart, iDimNonCart, iSub1, iSub2, iCart1, iCart2, i, j, k
    real:: CartValue1_I(1:7), CartValue2_I(1:7)
    !----------------------------------------------------------------------

    if(.not. allocated(CellCoef_DDGB)) then
       allocate(CellCoef_DDGB(&
            nDim,nDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
       CellCoef_DDGB = 0
    endif

    
    if(nK > 1) then ! nDim == 3

       !  ~       ~        ~
       ! dXi/dx, dXi/dy, dXi/dz
       iDimNonCart = 1
       iSub1 = mod(iDimNonCart, nDim) + 1
       iSub2 = mod(iDimNonCart+1, nDim) + 1
       do iDimCart = 1, nDim
          do k = 1, nK; do j = 1, nJ; do i = MinI, MaxI
             iCart1 = mod(iDimCart, nDim) + 1
             iCart2 = mod(iDimCart+1, nDim) + 1

             CartValue1_I = CellMetrice_DDG(iCart1,iSub1,i,j,k-3:k+3)* &
                  Xyz_DGB(iCart2,i,j,k-3:k+3,iBlock)
             CartValue2_I = CellMetrice_DDG(iCart1,iSub2,i,j-3:j+3,k)* &
                  Xyz_DGB(iCart2,i,j-3:j+3,k,iBlock)

             CellCoef_DDGB(iDimNonCart,iDimCart,i,j,k,iBlock) = &
                  calc_center_first_derivate(CartValue1_I, &
                  DxIn = CellSize_DB(iSub2,iBlock)) - &
                  calc_center_first_derivate(CartValue2_I, &
                  DxIn = CellSize_DB(iSub1,iBlock))
          enddo; enddo; enddo
       enddo
       CellCoef_DDGB(iDimNonCart,:,MinI:MaxI,1:nJ,1:nK,iBlock) = &
            CellCoef_DDGB(iDimNonCart,:,MinI:MaxI,1:nJ,1:nK,iBlock)&
            *CellSize_DB(iSub1,iBlock)&
            *CellSize_DB(iSub2,iBlock)


       !  ~         ~        ~
       ! dEta/dx, dEta/dy, dEta/dz
       iDimNonCart = 2
       iSub1 = mod(iDimNonCart, nDim) + 1
       iSub2 = mod(iDimNonCart+1, nDim) + 1
       do iDimCart = 1, nDim
          do k = 1, nK; do j = MinJ, MaxJ; do i = 1, nI
             iCart1 = mod(iDimCart, nDim) + 1
             iCart2 = mod(iDimCart+1, nDim) + 1

             CartValue1_I = CellMetrice_DDG(iCart1,iSub1,i-3:i+3,j,k)* &
                  Xyz_DGB(iCart2,i-3:i+3,j,k,iBlock)
             CartValue2_I = CellMetrice_DDG(iCart1,iSub2,i,j,k-3:k+3)* &
                  Xyz_DGB(iCart2,i,j,k-3:k+3,iBlock)

             CellCoef_DDGB(iDimNonCart,iDimCart,i,j,k,iBlock) = &
                  calc_center_first_derivate(CartValue1_I, &
                  DxIn = CellSize_DB(iSub2,iBlock)) - &
                  calc_center_first_derivate(CartValue2_I, &
                  DxIn = CellSize_DB(iSub1,iBlock))
          enddo; enddo; enddo
       enddo
       CellCoef_DDGB(iDimNonCart,:,1:nI,MinJ:MaxJ,1:nK,iBlock) = &
            CellCoef_DDGB(iDimNonCart,:,1:nI,MinJ:MaxJ,1:nK,iBlock)&
            *CellSize_DB(iSub1,iBlock)&
            *CellSize_DB(iSub2,iBlock)


       !   ~         ~         ~
       ! dZeta/dx, dZeta/dy, dZeta/dz
       iDimNonCart = 3
       iSub1 = mod(iDimNonCart, nDim) + 1
       iSub2 = mod(iDimNonCart+1, nDim) + 1
       do iDimCart = 1, nDim
          do k = MinK, MaxK; do j = 1, nJ; do i = 1, nI
             iCart1 = mod(iDimCart, nDim) + 1
             iCart2 = mod(iDimCart+1, nDim) + 1

             CartValue1_I = CellMetrice_DDG(iCart1,iSub1,i,j-3:j+3,k)* &
                  Xyz_DGB(iCart2,i,j-3:j+3,k,iBlock)
             CartValue2_I = CellMetrice_DDG(iCart1,iSub2,i-3:i+3,j,k)* &
                  Xyz_DGB(iCart2,i-3:i+3,j,k,iBlock)

             CellCoef_DDGB(iDimNonCart,iDimCart,i,j,k,iBlock) = &
                  calc_center_first_derivate(CartValue1_I, &
                  DxIn = CellSize_DB(iSub2,iBlock)) - &
                  calc_center_first_derivate(CartValue2_I, &
                  DxIn = CellSize_DB(iSub1,iBlock))
          enddo; enddo; enddo
       enddo
       CellCoef_DDGB(iDimNonCart,:,1:nI,1:nJ,MinK:MaxK,iBlock) = &
            CellCoef_DDGB(iDimNonCart,:,1:nI,1:nJ,MinK:MaxK,iBlock)&
            *CellSize_DB(iSub1,iBlock)&
            *CellSize_DB(iSub2,iBlock)

    elseif(nJ > 1) then ! nDim == 2       
       CellCoef_DDGB(Xi_,x_,:,:,:,iBlock)  =   &
            CellMetrice_DDG(y_,Eta_,:,:,:)*CellSize_DB(Eta_,iBlock)
       CellCoef_DDGB(Xi_,y_,:,:,:,iBlock)  = &
            - CellMetrice_DDG(x_,Eta_,:,:,:)*CellSize_DB(Eta_,iBlock)
       CellCoef_DDGB(Eta_,x_,:,:,:,iBlock) = &
            - CellMetrice_DDG(y_,Xi_,:,:,:)*CellSize_DB(Xi_,iBlock)
       CellCoef_DDGB(Eta_,y_,:,:,:,iBlock) = &
            CellMetrice_DDG(x_,Xi_,:,:,:)*CellSize_DB(Xi_,iBlock)
    else
       write(*,*) &
            'Warning: high-order does not available for 1D non-cartesian case!'
    endif

  end subroutine coef_cart_to_noncart
  !==========================================================================

  subroutine calc_metrics(iBlock) 
    ! Eq (10).

    ! Should be called from create_grid_block.
    ! Calc dx1/dx2 at cell center, where x1=x,y,z and x2=Xi,Eta,Zeta.

    integer, intent(in):: iBlock
    integer:: i, j, k, iDimCart, iDimNonCart
    real:: CellValue_I(7)
    !----------------------------------------------------------------------
    if(.not.allocated(CellMetrice_DDG)) then
       allocate(CellMetrice_DDG(nDim,nDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK))
       CellMetrice_DDG = 0.0
    endif

    ! dx/dXi, dy/dXi, dz/dXi
    iDimNonCart = 1
    do iDimCart = 1, nDim
       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = 1, nI
          CellValue_I = Xyz_DGB(iDimCart,i-3:i+3,j,k,iBlock)

          CellMetrice_DDG(iDimCart,iDimNonCart,i,j,k) = &
               calc_center_first_derivate(CellValue_I,&
               DxIn=CellSize_DB(iDimNonCart,iBlock))
       enddo; enddo; enddo ! i, j, k
    enddo ! iDimCart

    ! dx/dEta, dy/dEta, dz/dEta
    iDimNonCart = 2
    do iDimCart = 1, nDim
       do k = MinK, MaxK; do j = 1, nJ; do i = MinI, MaxI
          CellValue_I = Xyz_DGB(iDimCart,i,j-3:j+3,k,iBlock)

          CellMetrice_DDG(iDimCart,iDimNonCart,i,j,k) = &
               calc_center_first_derivate(CellValue_I,&
               DxIn=CellSize_DB(iDimNonCart,iBlock))
       enddo; enddo; enddo ! i, j, k
    enddo ! iDimCart

    if(nK >1) then
       ! dx/dZeta, dy/dZeta, dz/dZeta,
       iDimNonCart = 3
       do iDimCart = 1, nDim
          do k = 1, nK; do j = MinJ, MaxJ; do i = MinI, MaxI
             CellValue_I = Xyz_DGB(iDimCart,i,j,k-3:k+3,iBlock)

             CellMetrice_DDG(iDimCart,iDimNonCart,i,j,k) = &
                  calc_center_first_derivate(CellValue_I,&
                  DxIn=CellSize_DB(iDimNonCart,iBlock))
          enddo; enddo; enddo ! i, j, k
       enddo ! iDimCart
    endif
  end subroutine calc_metrics
  !==========================================================================

  subroutine test_grid

    use BATL_mpi, ONLY: iProc, nProc, iComm
    use BATL_geometry, ONLY: init_geometry
    use ModNumConst, ONLY: i_DD
    use ModMpi, ONLY: MPI_reduce, MPI_real, MPI_sum

    integer :: iBlock

    integer, parameter:: MaxBlockTest            = 50
    integer, parameter:: nRootTest_D(MaxDim)     = (/3,3,3/)
    logical, parameter:: IsPeriodicTest_D(MaxDim)= (/.true., .true., .false./)
    real:: DomainMin_D(MaxDim) = (/ 3.0, 2.0, 1.0 /)
    real:: DomainMax_D(MaxDim) = (/ 9.0, 6.0, 4.0 /)

    ! number of points in each dimension to test interpolation
    integer, parameter:: &
         nPointI = 91, &
         nPointJ = 1 + 90*min(1,nDim-1), &
         nPointK = 1 + 90*max(0,nDim-2)

    integer, parameter:: nPoint_D(MaxDim) = (/ nPointI, nPointJ, nPointK /)
    real, allocatable:: Point_VIII(:,:,:,:), PointAll_VIII(:,:,:,:)
    integer, parameter:: nVarPoint = nDim
    real:: XyzPoint_D(MaxDim), Point_V(nVarPoint), Weight
    integer:: iPoint, jPoint, kPoint, iPoint_D(MaxDim), iCell, nCell, iError
    integer:: iCell_II(0:nDim,2**nDim)
    logical:: IsSecondOrder
    integer:: iDiscr_D(MaxDim)
    real   :: Weight_I(2**nDim)

    real:: Tolerance

    integer:: i, j, k, Di, Dj, Dk, iDim, iBlockOut, iProcOut, iCell_D(MaxDim)
    integer:: iNodeCenter
    real:: Radius, Phi, Xyz_D(MaxDim), Coord_D(MaxDim), Distance_D(MaxDim)
    real:: Good, Good_D(MaxDim)
    real, allocatable:: CellVolumeCart_B(:), CellFaceCart_DB(:,:)

    logical:: DoTestMe
    character(len=*), parameter :: NameSub = 'test_grid'
    !-----------------------------------------------------------------------
    DoTestMe = iProc == 0

    if(DoTestMe)then
       write(*,*)'Starting ',NameSub
       write(*,*)'Testing init_grid'
       write(*,*)'nDimAmr, nIJK_D=', nDimAmr, nIJK_D
    end if
    ! Set Cartesian grid geometry before initializing tree and grid
    call init_geometry( IsPeriodicIn_D = IsPeriodicTest_D(1:nDim) )
    call init_tree(MaxBlockTest)
    call init_grid( DomainMin_D(1:nDim), DomainMax_D(1:nDim) )
    call set_tree_root( nRootTest_D(1:nDim))

    call find_tree_node((/0.5,0.5,0.5/),iNodeCenter)
    call refine_tree_node(iNodeCenter)
    call distribute_tree(.true.)
    if(DoTestMe) call show_tree('After distribute_tree')

    if(DoTestMe) write(*,*)'Testing create_grid'
    call create_grid

    if(iProc==0)call show_grid_proc

    if(DoTestMe) write(*,*)'Testing find_grid_block'
    Xyz_D = 0.0
    Xyz_D(1:nDim) = DomainMin_D(1:nDim)
    call find_grid_block(Xyz_D, iProcOut, iBlockOut, &
         iCell_D, Distance_D, UseGhostCell=.true.)
    if(iProc == iProcOut) then
       Xyz_D = Xyz_DGB(:,iCell_D(1),iCell_D(2),iCell_D(3),iBlockOut) &
            + 0.5*CellSize_DB(:,iBlockOut)
       if(any(abs(DomainMin_D(1:nDim) - Xyz_D(1:nDim)) > 1e-6)) then
          write(*,*) 'Error: DomainMin_D, Xyz_D=', &
               DomainMin_D, Xyz_D
          write(*,*) 'iProcOut, iBlockOut, iCell_D = ',&
               iProcOut, iBlockOut, iCell_D
       end if
    end if

    if(any(iCell_D(1:nDim) /= 0)) then
       write(*,*) 'Error: iCell_D=', iCell_D(1:nDim),' should be 0'
       write(*,*) 'iProcOut, iBlockOut, Distance_D = ',&
            iProcOut, iBlockOut, Distance_D
    end if

    if(any(abs(Distance_D(1:nDim) - 0.5) > 1e-6)) then
       write(*,*) 'Error: Distance_D=', Distance_D(1:nDim),' should be -0.5'
       write(*,*) 'iProcOut, iBlockOut, iCell_D = ',&
            iProcOut, iBlockOut, iCell_D
    end if

    Xyz_D = 0.0
    Xyz_D(1:nDim) = DomainMax_D(1:nDim)
    call find_grid_block(Xyz_D, iProcOut, iBlockOut, &
         iCell_D, Distance_D, UseGhostCell=.true.)
    if(iProc == iProcOut) then
       Xyz_D = Xyz_DGB(:,iCell_D(1),iCell_D(2),iCell_D(3),iBlockOut) &
            + 0.5*CellSize_DB(:,iBlockOut)
       if(any(abs(DomainMax_D(1:nDim) - Xyz_D(1:nDim)) > 1e-6)) then
          write(*,*) 'Error: DomainMax_D, Xyz_D=', &
               DomainMax_D, Xyz_D
          write(*,*) 'iProcOut, iBlockOut, iCell_D = ',&
               iProcOut, iBlockOut, iCell_D
       end if
    end if

    if(any(iCell_D(1:nDim) /= nIJK_D(1:nDim))) then
       write(*,*) 'Error: iCell_D=', iCell_D(1:nDim), &
            ' should be ', nIJK_D(1:nDim)
       write(*,*) 'iProcOut, iBlockOut, Distance_D = ',&
            iProcOut, iBlockOut, Distance_D
    end if

    if(any(abs(Distance_D(1:nDim) - 0.5) > 1e-6)) then
       write(*,*) 'Error: Distance_D=', Distance_D(1:nDim),' should be +0.5'
       write(*,*) 'iProcOut, iBlockOut, iCell_D = ',&
            iProcOut, iBlockOut, iCell_D
    end if

    if(DoTestMe) write(*,*)'Testing interpolate_grid'

    Xyz_D = 0.0
    if(.not.allocated(Point_VIII)) &
         allocate(Point_VIII(0:nVarPoint,nPointI,nPointJ,nPointK))
    Point_VIII = 0.0
    do kPoint = 1, nPointK; do jPoint = 1, nPointJ; do iPoint = 1, nPointI
       iPoint_D = (/ iPoint, jPoint, kPoint /)
       XyzPoint_D(1:nDim) = CoordMin_D(1:nDim) + (iPoint_D(1:nDim)-0.5) &
            *DomainSize_D(1:nDim)/nPoint_D(1:nDim)

       call interpolate_grid(XyzPoint_D, nCell, iCell_II, Weight_I)

       do iCell = 1, nCell
          Point_VIII(0,iPoint,jPoint,kPoint) = &
               Point_VIII(0,iPoint,jPoint,kPoint) + Weight_I(iCell)
          iBlock = iCell_II(0,iCell)
          iCell_D = 1
          iCell_D(1:nDim) = iCell_II(1:nDim,iCell)

          ! Interpolate the coordinates to check order of accuracy
          ! Note: Using array syntax in combination with the indirect
          ! iCell_D index fails with optimization for NAG v5.1
          do iDim = 1, nDim
             Xyz_D(iDim) = &
                  Xyz_DGB(iDim,iCell_D(1),iCell_D(2),iCell_D(3),iBlock)
          end do

          ! Take care of periodic dimensions: shift coordinates as necessary
          do iDim = 1, nDim
             if(.not.IsPeriodicTest_D(iDim)) CYCLE
             if(XyzPoint_D(iDim) < Xyz_D(iDim) - 2*CellSize_DB(iDim,iBlock)) &
                  Xyz_D(iDim) = Xyz_D(iDim) - DomainSize_D(iDim)

             if(XyzPoint_D(iDim) > Xyz_D(iDim) + 2*CellSize_DB(iDim,iBlock)) &
                  Xyz_D(iDim) = Xyz_D(iDim) + DomainSize_D(iDim)
          end do

          Point_VIII(1:nDim,iPoint,jPoint,kPoint) = &
               Point_VIII(1:nDim,iPoint,jPoint,kPoint) &
               + Weight_I(iCell)*Xyz_D(1:nDim)

          if(.false..and.iPoint==7.and.jPoint==1) write(*,*) &
               '!!!iProc,ijkPoint,iCell,iCell_D,iBlock,Weight,Xyz,Point=',&
               iProc, iPoint_D(1:nDim), iCell, iCell_D(1:nDim), &
               iBlock, Weight_I(iCell), &
               Xyz_DGB(1:nDim,iCell_D(1),iCell_D(2),iCell_D(3),iBlock), &
               Point_VIII(1:nDim,iPoint,jPoint,kPoint)

       end do
    end do; end do; end do

    ! Collect contributions from all processors to proc 0
    if(nProc > 1)then
       allocate(PointAll_VIII(0:nVarPoint,nPointI,nPointJ,nPointK))
       call MPI_reduce(Point_VIII, PointAll_VIII, size(PointAll_VIII),&
            MPI_REAL, MPI_SUM, 0, iComm, iError)
       Point_VIII = PointAll_VIII
       deallocate(PointAll_VIII)
    end if

    if(iProc == 0)then
       ! Check interpolated coordinate values against point coordinates
       do kPoint = 1, nPointK; do jPoint = 1, nPointJ; do iPoint = 1, nPointI
          iPoint_D = (/ iPoint, jPoint, kPoint /)
          Xyz_D(1:nDim) = CoordMin_D(1:nDim) + (iPoint_D(1:nDim)-0.5) &
               *DomainSize_D(1:nDim)/nPoint_D(1:nDim)

          Weight  = Point_VIII(0,iPoint,jPoint,kPoint)
          Point_V = Point_VIII(1:nDim,iPoint,jPoint,kPoint)/Weight

          if(abs(Weight - 1.0) < 1e-6)then
             Tolerance = 1e-6
          else
             Tolerance = 3e-2
          end if

          if(any(abs(Xyz_D(1:nDim) - Point_V) > Tolerance))then
             write(*,*) 'ERROR: Point_V=',Point_V(1:nDim),&
                  ' should be ',Xyz_D(1:nDim)
             write(*,*) 'Total weight=',Weight
             write(*,*) 'i,j,kPoint=', iPoint_D(1:nDim)
             write(*,*) 'CoordMin,Max=',CoordMin_D(1:nDim),CoordMax_D(1:nDim)
          end if
       end do; end do; end do
    end if

    if(DoTestMe) write(*,*)'Testing interpolate_grid_amr'
    Xyz_D = 0.0
    if(.not.allocated(Point_VIII)) &
         allocate(Point_VIII(0:nVarPoint,nPointI,nPointJ,nPointK))
    Point_VIII = 0.0
    do kPoint = 1, nPointK; do jPoint = 1, nPointJ; do iPoint = 1, nPointI
       iPoint_D = (/ iPoint, jPoint, kPoint /)
       XyzPoint_D(1:nDim) = CoordMin_D(1:nDim) + (iPoint_D(1:nDim)-0.5) &
            *DomainSize_D(1:nDim)/nPoint_D(1:nDim)

       call interpolate_grid_amr(XyzPoint_D, nCell, iCell_II, Weight_I,&
            IsSecondOrder)

       do iCell = 1, nCell
          Point_VIII(0,iPoint,jPoint,kPoint) = &
               Point_VIII(0,iPoint,jPoint,kPoint) + Weight_I(iCell)
          iBlock = iCell_II(0,iCell)
          iCell_D = 1
          iCell_D(1:nDim) = iCell_II(1:nDim,iCell)

          ! Interpolate the coordinates to check order of accuracy
          ! Note: Using array syntax in combination with the indirect
          ! iCell_D index fails with optimization for NAG v5.1
          do iDim = 1, nDim
             Xyz_D(iDim) = &
                  Xyz_DGB(iDim,iCell_D(1),iCell_D(2),iCell_D(3),iBlock)
          end do

          ! Take care of periodic dimensions: shift coordinates as necessary
          do iDim = 1, nDim
             if(.not.IsPeriodicTest_D(iDim)) CYCLE
             if(XyzPoint_D(iDim) < Xyz_D(iDim) - 2*CellSize_DB(iDim,iBlock)) &
                  Xyz_D(iDim) = Xyz_D(iDim) - DomainSize_D(iDim)

             if(XyzPoint_D(iDim) > Xyz_D(iDim) + 2*CellSize_DB(iDim,iBlock)) &
                  Xyz_D(iDim) = Xyz_D(iDim) +  DomainSize_D(iDim)
          end do

          ! if point is close to boundary => interpolation isn't of 2nd order
          where( .not.( IsSecondOrder.or.IsPeriodicTest_D(1:nDim) ) )&
               Xyz_D(1:nDim) = XyzPoint_D(1:nDim)

          Point_VIII(1:nDim,iPoint,jPoint,kPoint) = &
               Point_VIII(1:nDim,iPoint,jPoint,kPoint) &
               + Weight_I(iCell)*Xyz_D(1:nDim)
                    
       end do

    end do; end do; end do

    ! Collect contributions from all processors to proc 0
    if(nProc > 1)then
       allocate(PointAll_VIII(0:nVarPoint,nPointI,nPointJ,nPointK))
       call MPI_reduce(Point_VIII, PointAll_VIII, size(PointAll_VIII),&
            MPI_REAL, MPI_SUM, 0, iComm, iError)
       Point_VIII = PointAll_VIII
       deallocate(PointAll_VIII)
    end if

    if(iProc == 0)then
       ! Check interpolated coordinate values against point coordinates
       do kPoint = 1, nPointK; do jPoint = 1, nPointJ; do iPoint = 1, nPointI
          iPoint_D = (/ iPoint, jPoint, kPoint /)
          Xyz_D(1:nDim) = CoordMin_D(1:nDim) + (iPoint_D(1:nDim)-0.5) &
               *DomainSize_D(1:nDim)/nPoint_D(1:nDim)

          Weight  = Point_VIII(0,iPoint,jPoint,kPoint)
          Point_V = Point_VIII(1:nDim,iPoint,jPoint,kPoint)/Weight

          if(abs(Weight - 1.0) < 1e-6)then
             Tolerance = 1e-6
          else
             Tolerance = 3e-2
          end if

          if(any(abs(Xyz_D(1:nDim) - Point_V) > Tolerance))then
             write(*,*) 'ERROR: Point_V=',Point_V(1:nDim),&
                  ' should be ',Xyz_D(1:nDim)
             write(*,*) 'Total weight=',Weight
             write(*,*) 'i,j,kPoint=', iPoint_D(1:nDim)
             write(*,*) 'CoordMin,Max=',CoordMin_D(1:nDim),CoordMax_D(1:nDim)
          end if
       end do; end do; end do
    end if

    if(nDim==nDimAmr)then
       if(DoTestMe) write(*,*)'Testing interpolate_grid_amr_gc'
       Xyz_D = 0.0
       if(.not.allocated(Point_VIII)) &
            allocate(Point_VIII(0:nVarPoint,nPointI,nPointJ,nPointK))
       Point_VIII = 0.0
       do kPoint = 1, nPointK; do jPoint = 1, nPointJ; do iPoint = 1, nPointI
          iPoint_D = (/ iPoint, jPoint, kPoint /)
          XyzPoint_D(1:nDim) = CoordMin_D(1:nDim) + (iPoint_D(1:nDim)-0.5) &
               *DomainSize_D(1:nDim)/nPoint_D(1:nDim)
          call interpolate_grid_amr_gc(XyzPoint_D, nCell, iCell_II, Weight_I,&
               IsSecondOrder)

          do iCell = 1, nCell
             Point_VIII(0,iPoint,jPoint,kPoint) = &
                  Point_VIII(0,iPoint,jPoint,kPoint) + Weight_I(iCell)
             iBlock = iCell_II(0,iCell)
             iCell_D = 1
             iCell_D(1:nDim) = iCell_II(1:nDim,iCell)

             ! Interpolate the coordinates to check order of accuracy
             ! Note: Using array syntax in combination with the indirect
             ! iCell_D index fails with optimization for NAG v5.1
             do iDim = 1, nDim
                Xyz_D(iDim) = &
                     Xyz_DGB(iDim,iCell_D(1),iCell_D(2),iCell_D(3),iBlock)
             end do

             ! Take care of periodic dimensions: shift coordinates as necessary
             do iDim = 1, nDim
                if(.not.IsPeriodicTest_D(iDim)) CYCLE
                if(XyzPoint_D(iDim) < Xyz_D(iDim) - 2*CellSize_DB(iDim,iBlock)) &
                     Xyz_D(iDim) = Xyz_D(iDim) - DomainSize_D(iDim)

                if(XyzPoint_D(iDim) > Xyz_D(iDim) + 2*CellSize_DB(iDim,iBlock)) &
                     Xyz_D(iDim) = Xyz_D(iDim) +  DomainSize_D(iDim)
             end do

             ! Take care of ghost cells: shift coordinates for coarser neighbor
             iDiscr_D = 0
             where(iCell_D(1:nDim) < 1)
                iDiscr_D(1:nDim) =-1
             elsewhere(iCell_D(1:nDim) > nIJK_D(1:nDim))
                iDiscr_D(1:nDim) = 1
             end where
             ! check that neighbor is coarser
             if(DiLevelNei_IIIB(iDiscr_D(1),iDiscr_D(2),iDiscr_D(3),iBlock)==1)&
                  Xyz_D = Xyz_D + &
                  0.5 * (2*modulo(iCell_D, 2) - 1) * CellSize_DB(:,iBlock)
             
          ! if point is close to boundary => interpolation isn't of 2nd order
          where( .not.( IsSecondOrder.or.IsPeriodicTest_D(1:nDim) ) )&
               Xyz_D(1:nDim) = XyzPoint_D(1:nDim)

             Point_VIII(1:nDim,iPoint,jPoint,kPoint) = &
                  Point_VIII(1:nDim,iPoint,jPoint,kPoint) &
                  + Weight_I(iCell)*Xyz_D(1:nDim)

          end do
       end do; end do; end do

       ! Collect contributions from all processors to proc 0
       if(nProc > 1)then
          allocate(PointAll_VIII(0:nVarPoint,nPointI,nPointJ,nPointK))
          call MPI_reduce(Point_VIII, PointAll_VIII, size(PointAll_VIII),&
               MPI_REAL, MPI_SUM, 0, iComm, iError)
          Point_VIII = PointAll_VIII
          deallocate(PointAll_VIII)
       end if

       if(iProc == 0)then
          ! Check interpolated coordinate values against point coordinates
          do kPoint = 1, nPointK; do jPoint = 1, nPointJ; do iPoint = 1, nPointI
             iPoint_D = (/ iPoint, jPoint, kPoint /)
             Xyz_D(1:nDim) = CoordMin_D(1:nDim) + (iPoint_D(1:nDim)-0.5) &
                  *DomainSize_D(1:nDim)/nPoint_D(1:nDim)

             Weight  = Point_VIII(0,iPoint,jPoint,kPoint)
             Point_V = Point_VIII(1:nDim,iPoint,jPoint,kPoint)/Weight

             if(abs(Weight - 1.0) < 1e-6)then
                Tolerance = 1e-6
             else
                Tolerance = 3e-2
             end if

             if(any(abs(Xyz_D(1:nDim) - Point_V) > Tolerance))then
                write(*,*) 'ERROR: Point_V=',Point_V(1:nDim),&
                     ' should be ',Xyz_D(1:nDim)
                write(*,*) 'Total weight=',Weight
                write(*,*) 'i,j,kPoint=', iPoint_D(1:nDim)
                write(*,*) 'CoordMin,Max=',CoordMin_D(1:nDim),CoordMax_D(1:nDim)
             end if
          end do; end do; end do
       end if
    end if

       if(nDim == 2)then
          if(DoTestMe) write(*,*)'Testing create_grid in RZ geometry'

        ! Store Cartesian values for checking
        allocate(CellVolumeCart_B(MaxBlock), CellFaceCart_DB(MaxDim,MaxBlock))
        CellFaceCart_DB = CellFace_DB
        CellVolumeCart_B= CellVolume_B

          ! Clean Cartesian grid
          call clean_grid

          ! Initialize RZ grid
          call init_geometry(TypeGeometryIn = 'rz')
          call init_grid( DomainMin_D(1:nDim), DomainMax_D(1:nDim) )
          call create_grid
          if(iProc==0)call show_grid_proc

          ! Check relative to Cartesian
          Tolerance = 1e-6
          do iBlock = 1, nBlock
             do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
                if(abs( CellVolume_GB(i,j,k,iBlock) &
                     - abs(Xyz_DGB(2,i,j,k,iBlock))*CellVolumeCart_B(iBlock)) &
                     < Tolerance) CYCLE
                write(*,*)NameSub,' ERROR: incorrect cell volume=', &
                     CellVolume_GB(i,j,k,iBlock),' should be', &
                     abs(Xyz_DGB(2,i,j,k,iBlock))*CellVolumeCart_B(iBlock), &
                     ' at i,j,k,iBlock,iProc=', i, j, k, iBlock, iProc
             end do; end do; end do
             do iDim = 1, nDim
                Di = i_DD(1,iDim); Dj = i_DD(2,iDim)
                do k = 1, nK; do j = 1, nJ+Dj; do i = 1, nI+Di
                   Radius = 0.5*sum(abs(Xyz_DGB(2,i-Di:i,j-Dj:j,k,iBlock)))
                   if(abs(CellFace_DFB(iDim,i,j,k,iBlock) - &
                        Radius*CellFaceCart_DB(iDim,iBlock)) &
                        < Tolerance) CYCLE
                   write(*,*)NameSub,' ERROR: incorrect face area=', &
                        CellFace_DFB(iDim,i,j,k,iBlock),' should be', &
                        Radius*CellFaceCart_DB(iDim,iBlock), &
                        ' at iDim,i,j,k,iBlock,iProc=', &
                        iDim, i, j, k, iBlock, iProc

                end do; end do; end do
             end do
          end do
       end if

       if(nDim >= 2)then
          if(DoTestMe) write(*,*)'Testing create_grid in cylindrical geometry'

          ! Clean  grid
          call clean_grid

          ! Initialize cylindrical grid
          call init_geometry(TypeGeometryIn = 'cylindrical')

          DomainMin_D = (/1., 0., -0.5/)
          DomainMax_D = (/3., 90., 0.5/)

          ! This is temporary solution to keep the test working
          IsNodeBasedGrid = .false.
          call init_grid( DomainMin_D(1:nDim), DomainMax_D(1:nDim) )
          call create_grid
          if(iProc==0)call show_grid_proc

          ! Check relative to generalized coordinate volumes and areas
          Tolerance = 1e-6
          do iBlock = 1, nBlock
             do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
                Good = sqrt(sum(Xyz_DGB(1:2,i,j,k,iBlock)**2)) &
                     *CellVolume_B(iBlock)
                if(abs( CellVolume_GB(i,j,k,iBlock) - Good) < Tolerance) CYCLE
                write(*,*)NameSub,' ERROR: incorrect cell volume=', &
                     CellVolume_GB(i,j,k,iBlock),' should be', Good, &
                     ' at i,j,k,iBlock,iProc=', i, j, k, iBlock, iProc
             end do; end do; end do
             do iDim = 1, nDim
                Di = i_DD(1,iDim); Dj = i_DD(2,iDim); Dk = i_DD(3,iDim)
                do k = 1, nK+Dk; do j = 1, nJ+Dj; do i = 1, nI+Di
                   ! Calculate face center in generalized coordinates
                   Coord_D = CoordMin_DB(:,iBlock) + CellSize_DB(:,iBlock) &
                        *(/i-0.5*(1+Di),j-0.5*(1+Dj),k-0.5*(1+Dk)/)

                   Good = CellFace_DB(iDim,iBlock)
                   if(iDim /= 2) Good = Good*Coord_D(1)
                   if(abs(CellFace_DFB(iDim,i,j,k,iBlock) - Good) > Tolerance) &
                        write(*,*)NameSub,' ERROR: incorrect face area=', &
                        CellFace_DFB(iDim,i,j,k,iBlock),' should be', Good, &
                        ' at iDim,i,j,k,iBlock,iProc=', &
                        iDim, i, j, k, iBlock, iProc

                   Phi = Coord_D(2)
                   if(iDim == 1)then
                      Good_D = (/ cos(Phi), sin(Phi), 0.0 /)
                   elseif(iDim == 2)then
                      Good_D = (/ -sin(Phi), cos(Phi), 0.0 /)
                   else
                      Good_D = (/ 0.0, 0.0, 1.0 /)
                   end if
                   ! Multiply by area (for now)
                   Good_D = Good_D*CellFace_DFB(iDim,i,j,k,iBlock)

                   if(any( Tolerance < abs(FaceNormal_DDFB(:,iDim,i,j,k,iBlock) &
                        - Good_D(1:nDim)))) &
                        write(*,*)NameSub,' ERROR: incorrect face area=', &
                        FaceNormal_DDFB(:,iDim,i,j,k,iBlock),' should be',Good_D,&
                        ' at iDim,i,j,k,iBlock,iProc=', &
                        iDim, i, j, k, iBlock, iProc

                end do; end do; end do
             end do
          end do
       end if

    if(nDim == 3)then
       if(DoTestMe) write(*,*)'Testing create_grid in spherical geometry'

       ! Clean  grid
       call clean_grid

       ! Initialize cylindrical grid
       call init_geometry(TypeGeometryIn = 'spherical')

       DomainMin_D = (/1.,  45.0,  0.0 /)
       DomainMax_D = (/3., 135.0, 90.0/)

       ! This is temporary solution to keep the test working
       IsNodeBasedGrid = .false.
       call init_grid( DomainMin_D, DomainMax_D )
       call create_grid
       if(iProc==0)call show_grid_proc

       if(DoTestMe) write(*,*)'Testing create_grid in rlonlat geometry'

       ! Clean  grid
       call clean_grid

       ! Initialize cylindrical grid
       call init_geometry(TypeGeometryIn = 'rlonlat')

       DomainMin_D = (/1.,  0.0, -45.0 /)
       DomainMax_D = (/3., 90.0,  45.0 /)

       ! This is temporary solution to keep the test working
       IsNodeBasedGrid = .false.
       call init_grid( DomainMin_D, DomainMax_D )
       call create_grid
       if(iProc==0)call show_grid_proc

    end if

    if(DoTestMe) write(*,*)'Testing clean_grid'
    call clean_grid
    call clean_tree

  end subroutine test_grid

end module BATL_grid
