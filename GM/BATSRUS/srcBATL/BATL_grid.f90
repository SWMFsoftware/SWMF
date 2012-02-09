module BATL_grid

  use BATL_mpi, ONLY: iProc
  use BATL_size
  use BATL_tree
  use BATL_geometry

  implicit none

  SAVE

  private ! except

  public :: init_grid
  public :: clean_grid
  public :: create_grid
  public :: create_grid_block
  public :: find_grid_block
  public :: interpolate_grid
  public :: test_grid

  real, public:: &
       CoordMin_D(MaxDim),      & ! Min gen. coords of domain
       CoordMax_D(MaxDim)         ! Max gen. coordinates of domain

  real, public, allocatable::   &
       CoordMin_DB(:,:),        &    ! Min gen. coordinates of a block domain
       CoordMax_DB(:,:),        &    ! Max gen. coordinates of a block domain
       CellSize_DB(:,:),        &    ! Cell size in gen. coordinates
       CellFace_DB(:,:),        &    ! Cell faces for Cartesian grids
       CellFace_DFB(:,:,:,:,:), &    ! Cell faces for general grids
       CellVolume_B(:),         &    ! Cell volume for Cartesian grids
       CellVolume_GB(:,:,:,:),  &    ! Cell volume for general grids
       Xyz_DGB(:,:,:,:,:),      &    ! Cartesian cell centers coords
       FaceNormal_DDFB(:,:,:,:,:,:)  ! Normal face area vector

  logical, public:: IsNodeBasedGrid = .true.

  ! Local variables

  logical :: DoInitializeGrid = .true.

contains
  !============================================================================
  subroutine init_grid(CoordMinIn_D, CoordMaxIn_D, UseRadiusIn, UseDegreeIn)

    use ModNumConst, ONLY: cPi, cTwoPi, cDegToRad

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

    allocate(CoordMin_DB(MaxDim,MaxBlock))
    allocate(CoordMax_DB(MaxDim,MaxBlock))
    allocate(CellSize_DB(MaxDim,MaxBlock))

    allocate(CellFace_DB(MaxDim,MaxBlock))
    if(.not.IsCartesian) &
         allocate(CellFace_DFB(MaxDim,1:nI+1,1:nJ+1,1:nK+1,MaxBlock))

    allocate(CellVolume_B(MaxBlock))
    allocate(CellVolume_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    allocate(Xyz_DGB(MaxDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))

    if(.not.IsCartesian .and. .not. IsRzGeometry) &
         allocate(FaceNormal_DDFB(nDim,nDim,1:nI+1,1:nJ+1,1:nK+1,MaxBlock))

    ! Periodicity in the radial direction is not possible at all
    if(r_ > 0) IsPeriodic_D(r_) = .false.

    if(Phi_ > 0)then
       ! Enforce periodicity for cylindrical and spherical grids if
       ! there is a full grid in the Phi direction. 
       ! One can also have periodicity with a segment in Phi
       if( abs(CoordMax_D(Phi_) - CoordMin_D(Phi_) - cTwoPi) < 1e-6 ) &
            IsPeriodic_D(Phi_) = .true.
    end if

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
    if(allocated(FaceNormal_DDFB)) deallocate(FaceNormal_DDFB)

    CoordMin_D =  0.0
    CoordMax_D = -1.0

  end subroutine clean_grid

  !===========================================================================

  subroutine create_grid_block(iBlock, iNodeIn)

    use ModNumConst, ONLY: cHalfPi
    use ModCoordTransform, ONLY: cross_product

    ! Create geometrical information for block iBlock on the local PE

    integer, intent(in):: iBlock

    ! In case iNode_B is not set, iNodeIn can provide the node info
    integer, optional, intent(in):: iNodeIn

    real :: PositionMin_D(MaxDim), PositionMax_D(MaxDim), Coord_D(MaxDim)

    real, allocatable:: &
         rCell_I(:), rFace_I(:),  &
         SinPhi_I(:), CosPhi_I(:), SinPhiFace_I(:), CosPhiFace_I(:), &
         SinThetaFace_I(:), dCosTheta_I(:), &
         Xyz_DN(:,:,:,:)

    real :: Area, Phi, Theta, Dphi, Dz, Dtheta, Xyz_D(MaxDim), d_D(MaxDim)
    integer :: iNode, i, j, k

    real, parameter:: cThird = 1.0/3.0
    
    character(len=*), parameter:: NameSub = 'create_grid_block'
    !-------------------------------------------------------------------------
    if(present(iNodeIn))then
       iNode = iNodeIn
    else
       iNode = iNode_B(iBlock)
    end if
    call get_tree_position(iNode, PositionMin_D, PositionMax_D)

    CoordMin_DB(:,iBlock)= CoordMin_D + (CoordMax_D - CoordMin_D)*PositionMin_D
    CoordMax_DB(:,iBlock)= CoordMin_D + (CoordMax_D - CoordMin_D)*PositionMax_D

    CellSize_DB(:,iBlock) = (CoordMax_DB(:,iBlock) - CoordMin_DB(:,iBlock)) &
         / nIjk_D

    ! The cell volumes and face areas in generalized coordinates.
    ! For Cartesian grid same as physical volume and area.
    CellVolume_B(iBlock)  = product(CellSize_DB(:,iBlock))
    CellFace_DB(:,iBlock) = CellVolume_B(iBlock) / CellSize_DB(:,iBlock)

    if(IsCartesian .or. IsRzGeometry)then

       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
          Xyz_DGB(:,i,j,k,iBlock) = CoordMin_DB(:,iBlock) + &
               ( (/i, j, k/) - 0.5 ) * CellSize_DB(:,iBlock)
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
       else
          ! Also useful for Cartesian to keep code simple
          CellVolume_GB(:,:,:,iBlock) = CellVolume_B(iBlock)
       end if

    elseif(IsNodeBasedGrid)then

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

       ! Calculate face area vectors
       if(nDim == 2)then
          ! Calculate face area vectors as 90 degree rotations of edge vectors
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
          do k = 1, nK; do j = 1, nJ; do i = 1, nI+1
             FaceNormal_DDFB(:,x_,i,j,k,iBlock) = 0.5*cross_product( &
                  Xyz_DN(:,i,j+1,k+1) - Xyz_DN(:,i,j  ,k),           &
                  Xyz_DN(:,i,j  ,k+1) - Xyz_DN(:,i,j+1,k)          )

             CellFace_DFB(1,i,j,k,iBlock) = &
                  sqrt(sum(FaceNormal_DDFB(:,1,i,j,k,iBlock)**2))

          end do; end do; end do
          do k = 1, nK; do j = 1, nJ+1; do i = 1, nI
             FaceNormal_DDFB(:,y_,i,j,k,iBlock) = 0.5*cross_product( &
                  Xyz_DN(:,i+1,j,k+1) - Xyz_DN(:,i,j,k  ),           &
                  Xyz_DN(:,i+1,j,k  ) - Xyz_DN(:,i,j,k+1)          )

             CellFace_DFB(2,i,j,k,iBlock) = &
                  sqrt(sum(FaceNormal_DDFB(:,2,i,j,k,iBlock)**2))

          end do; end do; end do
          do k = 1, nK+1; do j = 1, nJ; do i = 1, nI
             FaceNormal_DDFB(:,z_,i,j,k,iBlock) = 0.5*cross_product( &
                  Xyz_DN(:,i+1,j+1,k) - Xyz_DN(:,i  ,j,k),           &
                  Xyz_DN(:,i  ,j+1,k) - Xyz_DN(:,i+1,j,k)          )

             CellFace_DFB(3,i,j,k,iBlock) = &
                  sqrt(sum(FaceNormal_DDFB(:,3,i,j,k,iBlock)**2))

          end do; end do; end do
       end if

       ! Calculate cell volumes
       if(nDim == 2)then
          ! Calculate cell volume as a sum of 2 triangle areas
          ! Also calculate cell center as the center of mass
          do j = MinJ, MaxJ; do i = MinI, MaxI
             Xyz_D = 0.0
             CellVolume_GB(i,j,1,iBlock) = &
                  volume3(i,j, i+1,j, i+1,j+1) + &
                  volume3(i,j, i+1,j+1, i,j+1)
             Xyz_DGB(:,i,j,1,iBlock) = Xyz_D/CellVolume_GB(i,j,1,iBlock)
          end do; end do
       else
          ! Calculate cell volume as a sum of 6 tetrahedra
          ! The tips of the tetrahedra are at the min position of the cell
          ! Also calculate cell center as the center of mass
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             Xyz_D = 0.0
             CellVolume_GB(i,j,k,iBlock) = &
                  volume4(i,j,k, i+1,j  ,k  , i+1,j+1,k  , i+1,j  ,k+1) + &
                  volume4(i,j,k, i+1,j+1,k+1, i+1,j  ,k+1, i+1,j+1,k  ) + &
                  volume4(i,j,k, i  ,j+1,k  , i  ,j+1,k+1, i+1,j+1,k  ) + &
                  volume4(i,j,k, i+1,j+1,k+1, i+1,j+1,k  , i  ,j+1,k+1) + &
                  volume4(i,j,k, i  ,j  ,k+1, i+1,j  ,k+1, i  ,j+1,k+1) + &
                  volume4(i,j,k, i+1,j+1,k+1, i  ,j+1,k+1, i+1,j  ,k+1)
             Xyz_DGB(:,i,j,k,iBlock) = Xyz_D/CellVolume_GB(i,j,k,iBlock)
          end do; end do; end do
       end if

       ! Cell center positions based on generalized coordinates
       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
          Coord_D = CoordMin_DB(:,iBlock) + &
               ( (/i, j, k/) - 0.5 ) * CellSize_DB(:,iBlock)
          call coord_to_xyz(Coord_D, Xyz_DGB(:,i,j,k,iBlock))
       end do; end do; end do

    else

       allocate(rCell_I(MinI:MaxI), rFace_I(MinI:MaxI+1))

       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
          Coord_D = CoordMin_DB(:,iBlock) + &
               ( (/i, j, k/) - 0.5 ) * CellSize_DB(:,iBlock)
          call coord_to_xyz(Coord_D, Xyz_DGB(:,i,j,k,iBlock))
       end do; end do; end do

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

          Dz = CellSize_DB(3,iBlock)

          ! dV = r*dr*dphi*dz
          do i = MinI, MaxI
             ! NOTE: for ghost cells beyond the axis r=0 can be negative
             CellVolume_GB(i,:,:,iBlock) = &
                  rCell_I(i)*(rFace_I(i+1)-rFace_I(i))*Dphi*Dz
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

       elseif(IsSpherical)then

          Dtheta = CellSize_DB(Theta_,iBlock)

          allocate( SinThetaFace_I(nJ+1), dCosTheta_I(MinJ:MaxJ))

          do j = 1, nJ+1
             SinThetaFace_I(j) = sin(CoordMin_DB(Theta_,iBlock) + (j-1)*Dtheta)
          end do

          do j = MinJ, MaxJ
             Theta = CoordMin_DB(Theta_,iBlock) + (j-1)*Dtheta
             ! Note the sign change
             dCosTheta_I(j) = cos(Theta) - cos(Theta+dTheta)
          end do

          ! dV = d(r^3)/3*dphi*d(cos theta)
          do j = MinJ, MaxJ; do i = MinI, MaxI
             ! NOTE: for ghost cells beyond the axis r=0 can be negative
             CellVolume_GB(i,j,:,iBlock) = &
                  cThird*(rFace_I(i+1)**3-rFace_I(i)**3)*Dphi*dCosTheta_I(j)
          end do; end do

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

       elseif(IsRLonLat)then
          Dtheta = CellSize_DB(Lat_,iBlock)

          allocate(SinThetaFace_I(nK+1), dCosTheta_I(MinK:MaxK))

          do k = 1, nK+1
             SinThetaFace_I(k) = cos(CoordMin_DB(Lat_,iBlock) + (k-1)*Dtheta)
          end do

          do k = MinK, MaxK
             Theta = cHalfPi - (CoordMin_DB(Lat_,iBlock) + (k-1)*Dtheta)
             ! Note the sign change
             dCosTheta_I(k) = cos(Theta-dTheta) - cos(Theta)
          end do

          ! dV = d(r^3)/3*dphi*d(cos theta)
          do k = MinK, MaxK; do i = MinI, MaxI
             ! NOTE: for ghost cells beyond the axis r=0 can be negative
             CellVolume_GB(i,:,k,iBlock) = &
                  cThird*(rFace_I(i+1)**3-rFace_I(i)**3)*Dphi*dCosTheta_I(k)
          end do; end do

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

       else
          call CON_stop(NameSub//': '//TypeGeometry// &
               ' geometry is not yet implemented')
       end if

       if(allocated(rCell_I)) &
            deallocate(rCell_I, rFace_I)
       if(allocated(SinPhi_I)) &
            deallocate(SinPhi_I, CosPhi_I, SinPhiFace_I, CosPhiFace_I)
       if(allocated(SinThetaFace_I)) &
            deallocate(SinThetaFace_I, dCosTheta_I)
       if(allocated(Xyz_DN)) &
            deallocate(Xyz_DN)

    end if

  contains
    !=========================================================================
    real function volume3(i1,j1, i2,j2, i3,j3)

      integer, intent(in):: i1,j1, i2,j2, i3,j3

      real, parameter:: cThird = 1.0/3.0

      real, dimension(2):: a_D, b_D, c_D
      !----------------------------------------------------------------------
      a_D = Xyz_DN(1:2,i1,j1,1)
      b_D = Xyz_DN(1:2,i2,j2,1) - a_D
      c_D = Xyz_DN(1:2,i3,j3,1) - a_D

      volume3 = 0.5*( b_D(1)*c_D(2) - b_D(2)*c_D(1))

      Xyz_D(1:2) = Xyz_D(1:2) + volume3*(a_D + cThird*(b_D + c_D))

    end function volume3
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

      Xyz_D = Xyz_D + volume4*(a_D + 0.25*(b_D + c_D + d_D))

    end function volume4

  end subroutine create_grid_block

  !===========================================================================
  subroutine create_grid

    integer:: iBlock
    !------------------------------------------------------------------------
    do iBlock = 1, nBlock
       if(Unused_B(iBlock))CYCLE
       call create_grid_block(iBlock)
    end do

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

  subroutine show_grid

    use BATL_mpi, ONLY: iProc, nProc, barrier_mpi

    ! Show all blocks sequentially on all processors, ie. show_grid 
    ! must be called from all processors of the MPI communicator iComm!

    integer:: iBlock, iPe
    !------------------------------------------------------------------------

    call barrier_mpi
    do iPe = 0, nProc - 1
       if(iPe == iProc) then
          do iBlock = 1, nBlock
             if(Unused_B(iBlock)) CYCLE
             call show_grid_block(iBlock)
          end do
       end if
       call barrier_mpi
    end do

  end subroutine show_grid

  !===========================================================================

  subroutine find_grid_block(XyzIn_D, &
       iProcOut, iBlockOut, iCellOut_D, DistOut_D)

    ! Find the processor and block containing location XyzIn_D. 
    ! If present, iCellOut_D returns the indexes of the closest cell center:
    ! 1 <= iCellOut_D <= nIjk_D
    ! If present, DistOut_D returns the signed distance to the closest 
    ! cell center divided by the cell size.
    ! DistOut_D can only be present if iCellOut_D is also present.

    real,    intent(in) :: XyzIn_D(MaxDim)        ! Cartesian coords of point
    integer, intent(out):: iBlockOut, iProcOut    ! Block and proc indexes

    integer, intent(out), optional:: iCellOut_D(MaxDim) ! Closest cell indexes
    real,    intent(out), optional:: DistOut_D(MaxDim)  ! Normalized distance

    real:: CoordTree_D(MaxDim), Coord_D(MaxDim)
    integer:: iNode

    character(len=*), parameter:: NameSub = 'find_grid_block'
    !------------------------------------------------------------------------
    ! Convert to generalized coordinates if necessary
    if(IsCartesian .or. IsRzGeometry)then
       Coord_D = XyzIn_D
    else
       call xyz_to_coord(XyzIn_D, Coord_D)
    end if
    ! Calculate normalized coordinates for tree search
    CoordTree_D = (Coord_D - CoordMin_D)/(CoordMax_D - CoordMin_D)

    ! Find node containing the point
    if(present(iCellOut_D))then
       call find_tree_cell(CoordTree_D, iNode, iCellOut_D, DistOut_D)
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
    integer:: iBlock, iDim, DiLevel, iCell_D(MaxDim)
    integer:: i, j, k, iLo, jLo, kLo, iHi, jHi, kHi

    logical, parameter:: DoTest = .false.
    !------------------------------------------------------------------------
    ! Convert to generalized coordinates if necessary
    if(IsCartesian .or. IsRzGeometry)then
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
          Shift = CoordMax_D(iDim) - CoordMin_D(iDim)
          if(Coord_D(iDim) < CoordMin) Coord_D(iDim) = Coord_D(iDim) + Shift
          if(Coord_D(iDim) > CoordMax) Coord_D(iDim) = Coord_D(iDim) - Shift
          if(Coord_D(iDim) < CoordMin .or. Coord_D(iDim) > CoordMax) &
               CYCLE LOOPBLOCK

       end do

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
       if(            DiLevelNei_IIIB(-1,0,0,iBlock)/=Unset_) iLo = max(1,iLo)
       if(nDim<2 .or. DiLevelNei_IIIB(0,-1,0,iBlock)/=Unset_) jLo = max(1,jLo)
       if(nDim<3 .or. DiLevelNei_IIIB(0,0,-1,iBlock)/=Unset_) kLo = max(1,kLo)

       if(            DiLevelNei_IIIB(+1,0,0,iBlock)/=Unset_) iHi = min(nI,iHi)
       if(nDim<2 .or. DiLevelNei_IIIB(0,+1,0,iBlock)/=Unset_) jHi = min(nJ,jHi)
       if(nDim<3 .or. DiLevelNei_IIIB(0,0,+1,iBlock)/=Unset_) kHi = min(nK,kHi)

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
             write(*,*)'ERROR in interpolate_grid: too many cells!'
             write(*,*)'iProc,iBlock,nDim=',iProc,iBlock,nDim
             write(*,*)'iLo,iHi,jLo,jHi,kLo,kHi=',iLo,iHi,jLo,jHi,kLo,kHi
             write(*,*)'Coord_D=',Coord_D(1:nDim)
          end if

          Weight_I(nCell)   = Weight
          iCell_II(0,nCell) = iBlock
          iCell_II(1:nDim,nCell) = iCell_D(1:nDim)

          if(DoTest)write(*,*)'!!! nCell, CoordCell, Weight=', &
               nCell, CoordCell_D(1:nDim), Weight

       end do; end do; end do
    end do LOOPBLOCK

  end subroutine interpolate_grid
  !===========================================================================

  subroutine test_grid

    use BATL_mpi, ONLY: iProc, nProc, iComm
    use BATL_geometry, ONLY: init_geometry
    use ModNumConst, ONLY: i_DD
    use ModMpi, ONLY: MPI_reduce, MPI_real, MPI_sum

    integer :: iBlock

    integer, parameter:: MaxBlockTest            = 50
    integer, parameter:: nRootTest_D(MaxDim)     = (/3,2,1/)
    logical, parameter:: IsPeriodicTest_D(MaxDim)= (/.true., .true., .false./)
    real:: DomainMin_D(MaxDim) = (/ 3.0, 2.0, 1.0 /)
    real:: DomainMax_D(MaxDim) = (/ 9.0, 6.0, 4.0 /)

    ! number of points in each dimension to test interpolation
    integer, parameter:: &
         nPointI = 51, &
         nPointJ = 1 + 50*min(1,nDim-1), &
         nPointK = 1 + 50*max(0,nDim-2)

    integer, parameter:: nPoint_D(MaxDim) = (/ nPointI, nPointJ, nPointK /)
    real, allocatable:: Point_VIII(:,:,:,:), PointAll_VIII(:,:,:,:)
    integer, parameter:: nVarPoint = nDim
    real:: XyzPoint_D(MaxDim), Point_V(nVarPoint), Weight
    integer:: iPoint, jPoint, kPoint, iPoint_D(MaxDim), iCell, nCell, iError
    integer:: iCell_II(0:nDim,2**nDim)
    real   :: Weight_I(2**nDim)

    real:: Tolerance

    integer:: i, j, k, Di, Dj, Dk, iDim, iBlockOut, iProcOut, iCell_D(MaxDim)
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

    call refine_tree_node(3)
    call distribute_tree(.true.)
    if(DoTestMe) call show_tree('After distribute_tree')

    if(DoTestMe) write(*,*)'Testing create_grid'
    call create_grid

    call show_grid

    if(DoTestMe) write(*,*)'Testing find_grid_block'
    call find_grid_block(DomainMin_D, iProcOut, iBlockOut, &
         iCell_D, Distance_D)
    if(iProc == iProcOut) then
       Xyz_D = Xyz_DGB(:,iCell_D(1),iCell_D(2),iCell_D(3),iBlockOut) &
            - 0.5*CellSize_DB(:,iBlockOut)
       if(any(abs(DomainMin_D(1:nDim) - Xyz_D(1:nDim)) > 1e-6)) then
          write(*,*) 'Error: DomainMin_D, Xyz_D=', &
               DomainMin_D, Xyz_D
          write(*,*) 'iProcOut, iBlockOut, iCell_D = ',&
               iProcOut, iBlockOut, iCell_D
       end if
    end if

    if(any(iCell_D(1:nDim) /= 1)) then
       write(*,*) 'Error: iCell_D=', iCell_D(1:nDim),' should be 1'
       write(*,*) 'iProcOut, iBlockOut, Distance_D = ',&
            iProcOut, iBlockOut, Distance_D
    end if

    if(any(abs(Distance_D(1:nDim) + 0.5) > 1e-6)) then
       write(*,*) 'Error: Distance_D=', Distance_D(1:nDim),' should be -0.5'
       write(*,*) 'iProcOut, iBlockOut, iCell_D = ',&
            iProcOut, iBlockOut, iCell_D
    end if

    call find_grid_block(DomainMax_D, iProcOut, iBlockOut, &
         iCell_D, Distance_D)
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
            *(CoordMax_D(1:nDim) - CoordMin_D(1:nDim))/nPoint_D(1:nDim)

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
                  Xyz_D(iDim) = Xyz_D(iDim) - &
                  (CoordMax_D(iDim) - CoordMin_D(iDim))

             if(XyzPoint_D(iDim) > Xyz_D(iDim) + 2*CellSize_DB(iDim,iBlock)) &
                  Xyz_D(iDim) = Xyz_D(iDim) + &
                  (CoordMax_D(iDim) - CoordMin_D(iDim))
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
               *(CoordMax_D(1:nDim) - CoordMin_D(1:nDim))/nPoint_D(1:nDim)

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
       call show_grid

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

    ! We should probably have a single option.

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
       call show_grid

       ! Check relative to generalized coordinate volumes and areas
       Tolerance = 1e-6
       do iBlock = 1, nBlock
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             Good = sqrt(sum(Xyz_DGB(1:2,i,j,k,iBlock)**2))*CellVolume_B(iBlock)
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
       call show_grid

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
       call show_grid

    end if

    if(DoTestMe) write(*,*)'Testing clean_grid'
    call clean_grid
    call clean_tree
    
  end subroutine test_grid

end module BATL_grid
