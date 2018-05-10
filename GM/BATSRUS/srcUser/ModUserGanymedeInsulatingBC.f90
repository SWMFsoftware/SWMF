!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModUser

  use BATL_lib, ONLY: &
       test_start, test_stop, iTest, jTest, kTest, iBlockTest
  use ModProcMH, ONLY: iProc,iComm
  ! User module for Ganymede, insulating boundary model from Duling et.al [2014]
  ! B1 is used instead of B1+B0 as input/output to the function.
  ! Must compile with MHDHypPe equation set.
  use ModUserEmpty,                           &
       IMPLEMENTED1 => user_set_ics,          &
       IMPLEMENTED2 => user_set_cell_boundary,&
       IMPLEMENTED3 => user_action
       
  include 'user_module.h' ! list of public methods
  
  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = 'Ganymede insulating, Hongyang Zhou'

  real, parameter :: PlanetRadius = 1.

  integer :: iProcInner  ! Processor rank for inner boundary
  integer :: nProcInner  ! Number of processors containing the inner boundary
  integer :: iCommInner  ! MPI communicator for inner boundary group
  integer :: iError
  
contains
  !============================================================================

  subroutine user_set_ics(iBlock)

    use ModAdvance,    ONLY: State_VGB
    use ModGeometry,   ONLY: Rmin_BLK, R_BLK
    use ModSize,       ONLY: nI, nJ, nK, nG
    use ModVarIndexes, ONLY: Bx_, Bz_, Rho_, P_, Pe_
    use ModMultiFluid, ONLY: select_fluid, iFluid, nFluid, iP, &
         iRho, iRhoUx, iRhoUz, iRhoUx_I, iRhoUz_I
    use ModPhysics,    ONLY: CellState_VI
    
    integer, intent(in) :: iBlock

    integer :: i,j,k
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_set_ics'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)     
    
    ! Modify near planet condition only
    if(Rmin_BLK(iBlock) > 5.0*PlanetRadius) RETURN
    
    do k=1,nK; do j=1,nJ; do i=1,nI
       do iFluid = 1,nFluid
          iRhoUx = iRhoUx_I(iFluid); iRhoUz = iRhoUz_I(iFluid)
          if(R_BLK(i,j,k,iBlock) < 2.5*PlanetRadius)then
             State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = 0.0
          elseif(R_BLK(i,j,k,iBlock) < 5.0*PlanetRadius)then
             State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = &
                  CellState_VI(iRhoUx:iRhoUz,2) * (R_BLK(i,j,k,iBlock)/2.5 - 1)
          end if
          
       end do
    end do; end do; end do
    
    call test_stop(NameSub, DoTest, iBlock)
  end subroutine user_set_ics
  !============================================================================

  subroutine user_set_cell_boundary(iBlock, iSide, TypeBc, IsFound)

    use ModAdvance,    ONLY: State_VGB
    use ModGeometry,   ONLY: r_BLK
    use ModSize,       ONLY: nI, nJ, nK, nG
    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUz_, Pe_, Bx_, Bz_, &
                             iRho_I,iRhoUx_I, iRhoUy_I, iRhoUz_I, iP_I
    use ModEnergy,     ONLY: calc_energy_cell
    use BATL_lib,      ONLY: Xyz_DGB, CellSize_DB, nRoot_D
    use ModMain,       ONLY: nBlock, MaxBlock 
    use ModPhysics,    ONLY: No2Si_V, UnitB_, CellState_VI
    use ModB0,         ONLY: B0_DGB
    use ModMultiFluid, ONLY: iFluid, nFluid, iRhoUx, iRhoUz
    use ModMpi
    use ModCellBoundary, ONLY: iMin, iMax, jMin, jMax, kMin, kMax
    use ModCoordTransform, ONLY: rot_xyz_sph, xyz_to_sph
    use ModNumConst,   ONLY: cRadToDeg, cPi
    use ModParallel,   ONLY: NOBLK, NeiLev
    
    integer,          intent(in)  :: iBlock, iSide
    character(len=*), intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    real :: r_D(3), b_D(3), u_D(3), RhoUr, XyzSph_DD(3,3)

    ! Insulating BC setup
    real :: BoundaryCellCoord_DG(3,nJ,nK)
    real :: GhostCellCoord_D(3,-1:0)

    real, allocatable, save :: Br_S(:,:)
    real, allocatable, save :: BrAllProc_S(:,:)
    real, allocatable, save :: BrGhost_S(:,:), BtGhost_S(:,:), BpGhost_S(:,:)
    real, allocatable, save :: t(:), p(:)  ! theta and phi array
    real, allocatable, save :: dt(:), dp(:)
    ! Lookup table for grid mapping 
    integer, allocatable, save :: index_lookuptable(:,:,:)
    
    ! Number of cells in theta-direction of the simulation grid (N_theta)
    integer, save :: Nt
    ! Number of cells in phi-direction of the simulation grid (N_phi)
    integer, save :: Np    
    ! Maximum transformation degree (N_l)
    integer, parameter :: Nl = 20
    ! Radius of boundary surface (R_0)
    real :: R0 = 1.0 ! I may need to change this to cell center!
       
    logical, save :: IsFirstCall=.true.
    integer, save :: MaxInnerBCBlock ! Max block index on each proc for inner BC
    integer :: i,j,k, jIndex, kIndex, blockIndex, iVar

    ! Switch to save if constants are calculated
    integer, save :: constants = 0
    
    ! Constants for transformation in spatial-spectral space
    real, allocatable, save :: lmc(:,:,:,:) ! cos part
    real, allocatable, save :: lms(:,:,:,:) ! sin part
    
    ! Constants for back-transformation
    real, allocatable, save :: lmb1c(:,:,:,:) ! B_r     cos part
    real, allocatable, save :: lmb1s(:,:,:,:) ! B_r     sin part
    real, allocatable, save :: lmb2c(:,:,:,:) ! B_theta cos part
    real, allocatable, save :: lmb2s(:,:,:,:) ! B_theta sin part
    real, allocatable, save :: lmb3c(:,:,:,:) ! B_phi   cos part
    real, allocatable, save :: lmb3s(:,:,:,:) ! B_phi   sin part
    
    ! Gauss coefficients of internal field
    real, dimension(1:Nl,0:Nl), save :: Glm, Hlm
    
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_set_cell_boundary'
    !---------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    
    ! Outer boundary condition
    if(iSide == 2)then
       do k=kMin,kMax; do j=jMin,jMax; do i=iMin,iMax
          if(Xyz_DGB(1,i,j,k,iBlock) < 0)then
             ! Upstream fixed
             State_VGB(:,i,j,k,iBlock) = CellState_VI(:,iSide)
          else
             ! Downstream float
             State_VGB(:,i,j,k,iBlock) = State_VGB(:,nI,j,k,iBlock)
          end if
       end do; end do; end do
       
       IsFound = .true.
       RETURN
    end if
    
    ! Inner boundary condition
    ! The only requirement is that boundary grid must be uniform!
    
    ! Initial setup and calculation when first being called
    if(IsFirstCall) then
       ! Get the maximum block index of the inner boundary block
       MaxInnerBCBlock = nBlock
       do while( NeiLev(1,MaxInnerBCBlock) /= NOBLK)
          MaxInnerBCBlock = MaxInnerBCBlock - 1
       end do

       ! Get the number of grid points on the surface
       Nt = nRoot_D(2) * nJ
       Np = nRoot_D(3) * nK

       ! Initialization on first call
       if(iProcInner == 0)then
          if(.not.allocated(lmc)) then
             allocate(lmc(Nt,Np,1:Nl,0:Nl))
             allocate(lms(Nt,Np,1:Nl,0:Nl))
             allocate(lmb1c(1:Nl,0:Nl,Nt,Np))
             allocate(lmb1s(1:Nl,0:Nl,Nt,Np))
             allocate(lmb2c(1:Nl,0:Nl,Nt,Np))
             allocate(lmb2s(1:Nl,0:Nl,Nt,Np))
             allocate(lmb3c(1:Nl,0:Nl,Nt,Np))
             allocate(lmb3s(1:Nl,0:Nl,Nt,Np))
             allocate(BrAllProc_S(Nt,Np))
             lmc   = 0.0; lms   = 0.0
             lmb1c = 0.0; lmb1s = 0.0
             lmb2c = 0.0; lmb2s = 0.0
             lmb3c = 0.0; lmb3s = 0.0
             BrAllProc_S = 0.0
          end if
       end if
       
       if(.not.(allocated(Br_S))) then
          allocate(Br_S(Nt,Np))
          allocate(BrGhost_S(Nt,Np))
          allocate(BtGhost_S(Nt,Np))
          allocate(BpGhost_S(Nt,Np))
          Br_S = 0.0
          BrGhost_S = 0.0
          BtGhost_S = 0.0
          BpGhost_S = 0.0
          
          allocate(index_lookuptable(3,Nt,Np))
          index_lookuptable = 0
       end if
        
       if(.not.allocated(t)) then
          allocate(t(nRoot_D(2)*nJ))
          allocate(p(nRoot_D(3)*nK))
          allocate(dt(nRoot_D(2)*nJ))
          allocate(dp(nRoot_D(3)*nK))

          ! Theta grid and interval
          do j = 1,Nt
             t(j)  = cPi/Nt*(j-0.5)
             dt(j) = cPi/Nt
          end do

          ! Phi grid and interval
          do k = 1,Np
             p(k)  = 2*cPi/Np*(k-0.5)
             dp(k) = 2*cPi/Np
          end do
          
       end if
       
       IsFirstCall = .false.
    end if
    
    ! Get the closest layer of physical cells' positions and values
    ! and store into BoundaryCellCoord_DG
    do k = 1, nK; do j = 1, nJ
       call xyz_to_sph(Xyz_DGB(:,1,j,k,iBlock),BoundaryCellCoord_DG(:,j,k))
    end do; end do

    ! Set the radius to the cell center radius (turn it on later)
    !R0 = BoundaryCellCoord_DG(1,1,1)
    
    ! Get the ghost cell positions (there are duplicated calculations)
    call xyz_to_sph(Xyz_DGB(:,iMin,1,1,iBlock),GhostCellCoord_D(:,iMin))
    call xyz_to_sph(Xyz_DGB(:,iMax,1,1,iBlock),GhostCellCoord_D(:,iMax))
        
    ! From coords. find index and calculate Br
    do k=1,nK; do j=1,nJ
       ! Get the global boundary indexes from coords.
       jIndex = BoundaryCellCoord_DG(2,j,k) / dt(1) + 1
       kIndex = BoundaryCellCoord_DG(3,j,k) / dp(1) + 1

       ! Mapping the grid
       Index_LookupTable(:,jIndex,kIndex) = [j,k,iBlock]
       
       ! Get Br from (Bx,By,Bz), and use B1 as input and output
       r_D = Xyz_DGB(x_:z_,1,j,k,iBlock) / r_BLK(1,j,k,iBlock)
       Br_S(jIndex,kIndex) = dot_product( State_VGB(Bx_:Bz_,1,j,k,iBlock),r_D )
    end do; end do    

    ! Message pass when all the boundary blocks are found per proc
    if(iBlock == MaxInnerBCBlock)then
       
       call MPI_Reduce(Br_S,BrAllProc_S,Nt*Np,MPI_REAL,MPI_SUM,&
            0,iCommInner,iError)

       do i=iMin,iMax
          ! Calculate boundary values at ghost cell radius
          if(iProc == 0) then          
             call insulating_bvalues(BrAllProc_S,t,p,dt,dp,&
                  GhostCellCoord_D(1,i), BrGhost_S,BtGhost_S,BpGhost_S)
          end if
       
          call MPI_Bcast(BrGhost_S,Nt*Np,MPI_REAL,0,iCommInner,iError)
          call MPI_Bcast(BtGhost_S,Nt*Np,MPI_REAL,0,iCommInner,iError)
          call MPI_Bcast(BpGhost_S,Nt*Np,MPI_REAL,0,iCommInner,iError)
       
          ! (Br,Bt,Bp) --> (Bx,By,Bz) and set layer i ghost cell values
          do kIndex=1,Np; do jIndex=1,Nt
             j = Index_LookupTable(1,jIndex,kIndex)
             k = Index_LookupTable(2,jIndex,kIndex)
             BlockIndex = Index_LookupTable(3,jIndex,kIndex)
             if(BlockIndex == 0) cycle
             
             ! Matrix for converting between Cartesian and spherical Coords.
             XyzSph_DD = rot_xyz_sph(Xyz_DGB(:,i,j,k,BlockIndex))
             
             ! Spherical --> Cartesian
             State_VGB(Bx_:Bz_,i,j,k,BlockIndex) = matmul(XyzSph_DD,&
                  [BrGhost_S(jIndex,kIndex),BtGhost_S(jIndex,kIndex),&
                  BpGhost_S(jIndex,kIndex)])
          end do; end do
       end do     
    end if

    ! Set density, pressure and velocity boundary conditions
    do k=kMin,kMax; do j=jMin,jMax; do i=iMin,iMax      
       ! Velocity perpendicular to B is float, parallel to B is 0
!       b_D = B0_DGB(:,1,j,k,iBlock) + State_VGB(Bx_:Bz_,1,j,k,iBlock)
!       b_D = b_D / sqrt(sum(b_D**2))
!       do iFluid = 1, nFluid
!          iRhoUx = iRhoUx_I(iFluid); iRhoUz = iRhoUz_I(iFluid)
!          State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = &
!               State_VGB(iRhoUx:iRhoUz,1,j,k,iBlock) - &
!               sum(b_D*State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock))*b_D
!       end do

       ! radial vector
       r_D = Xyz_DGB(x_:z_,1,j,k,iBlock)/ r_BLK(1,j,k,iBlock)

       ! float density
       !State_VGB(iRho_I,i,j,k,iBlock) = State_VGB(iRho_I,1,j,k,iBlock)
       
       ! No outflow, inflow float
       do iFluid = 1, nFluid
          iRhoUx = iRhoUx_I(iFluid); iRhoUz = iRhoUz_I(iFluid)

          ! Get radial momentum          
          RhoUr = dot_product(State_VGB(iRhoUx:iRhoUz,1,j,k,iBlock),r_D)
          if(RhoUr>0)then
             ! Push it harder
             State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = &
                  State_VGB(iRhoUx:iRhoUz,1,j,k,iBlock) - 2*RhoUr*r_D
             
             !State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = &
             !     State_VGB(iRhoUx:iRhoUz,1,j,k,iBlock) - RhoUr*r_D
             
             u_D = State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) / &
                  State_VGB(Rho_,i,j,k,iBlock)
             
             ! Float density
             State_VGB(iRho_I,i,j,k,iBlock) = State_VGB(iRho_I,1,j,k,iBlock)
             
             State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = &
                  u_D*State_VGB(Rho_,i,j,k,iBlock)
             
          else
             ! inflow float
             State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = &
                  State_VGB(iRhoUx:iRhoUz,1,j,k,iBlock)
             ! Float density
             State_VGB(iRho_I,i,j,k,iBlock) = State_VGB(iRho_I,1,j,k,iBlock)
          end if
       end do
       
       ! Float pressure
       State_VGB(iP_I,i,j,k,iBlock) = State_VGB(iP_I,1,j,k,iBlock)
       State_VGB(Pe_,i,j,k,iBlock) = State_VGB(Pe_,1,j,k,iBlock)
       
    end do; end do; end do

    IsFound = .true.
    
    call calc_energy_cell(iBlock)
    
    call test_stop(NameSub, DoTest, iBlock)

  contains
    !--------------------------------------------------------------
    subroutine insulating_bvalues(br,t,p,dt,dp,r_b,br_b,bt_b,bp_b)
      
      implicit none
      
      real, intent(in) :: br(Nt,Np)
      real, intent(in) :: t(Nt)
      real, intent(in) :: p(Np)
      real, intent(in) :: dt(Nt)
      real, intent(in) :: dp(Np)
      real, intent(in) :: r_b
      
      real, intent(out) :: br_b(Nt,Np)
      real, intent(out) :: bt_b(Nt,Np)
      real, intent(out) :: bp_b(Nt,Np)
      
      ! Temporal sum variables
      real :: summc
      real :: summs
      real :: summb1
      real :: summb2
      real :: summb3
      
      ! Cos and Sin (real and imagenary) part of poloidal coefficients at R_0
      real :: pc(1:Nl,0:Nl)
      real :: ps(1:Nl,0:Nl)
      
      ! Cos and Sin (real and imagenary) part of poloidal coefficients at r_b
      real :: pc_b(1:Nl,0:Nl)
      real :: ps_b(1:Nl,0:Nl)
      ! Gradient of Cos and Sin part of poloidal coefficients at r_b
      real :: dpc_b(1:Nl,0:Nl)
      real :: dps_b(1:Nl,0:Nl)
      
      ! Loop variables
      integer :: l, m, j, k

      
      ! Calculate time indepedent constants only once during the first call
      ! of insulating_bvalues to save computational time.
      if (constants .eq. 0) then
        call insulating_constants(t,p,dt,dp)
        constants = 1
      endif
      
! ------------------------- Step 1: Transformation --------------------------
! Here the magnetic field at the boundary surface is transformed into 
! spatial-spectral space und the poloidal coefficients are calculated.

      do l=1,Nl
        do m=0,l
      summc=0.0
      summs=0.0
        
      ! Integration over grid cells
      do k=1,Np
        do j=1,Nt
      summc = summc + br(j,k)*lmc(j,k,l,m)
      summs = summs + br(j,k)*lms(j,k,l,m)
        enddo
      enddo

      ! poloidal coefficients
      pc(l,m) = summc
      ps(l,m) = summs

        enddo
      enddo
      
! ----------------------- Step 2: Boundary coefficients ---------------------
! Here the poloidal coefficients p_lm and their gradients are calculated 
! where the boundary values are needed.

      do l=1,Nl
        do m=0,l
            ! Poloidal coefficients at r_b
            pc_b(l,m) = pc(l,m)*(r_b/R0)**(l+1)
            ps_b(l,m) = ps(l,m)*(r_b/R0)**(l+1)

            ! Add intrinsic fields
            pc_b(l,m) = pc_b(l,m)+Glm(l,m)/l* &
                 ( (R0/r_b)**l - (r_b/R0)**(l+1) )
            ps_b(l,m) = ps_b(l,m)+Hlm(l,m)/l* &
                 ( (R0/r_b)**l - (r_b/R0)**(l+1) )

            ! Gradients of polodial coefficients at r_b
            dpc_b(l,m) = (l+1)/R0    *pc (l,m)*(r_b/R0)**l
            dps_b(l,m) = (l+1)/R0    *ps (l,m)*(r_b/R0)**l

            ! Add intrinsic fields
            dpc_b(l,m) = dpc_b(l,m)+Glm(l,m)/(l*R0) &
                 *(-l*(R0/r_b)**(l+1) - (l+1)*(r_b/R0)**l )
            dps_b(l,m) = dps_b(l,m)+Hlm(l,m)/(l*R0) &         
                 *(-l*(R0/r_b)**(l+1) - (l+1)*(r_b/R0)**l )
        enddo
      enddo
      
! ------------------------ Step 3: Back-Transformation ----------------------
! Here all three components of the boundary values are calculated

      do k=1,Np
        do j=1,Nt
      summb1 = 0.0
      summb2 = 0.0
      summb3 = 0.0
      ! Summ over spherical harmonics
      do l=1,Nl
        do m=0,l
      summb1 = summb1 + lmb1c(l,m,j,k)*pc_b(l,m) & ! B_b,r
                      + lmb1s(l,m,j,k)*ps_b(l,m)
      summb2 = summb2 + lmb2c(l,m,j,k)*dpc_b(l,m)& ! B_b,theta
                      + lmb2s(l,m,j,k)*dps_b(l,m)
      summb3 = summb3 + lmb3c(l,m,j,k)*dpc_b(l,m)& ! B_b,phi
                      + lmb3s(l,m,j,k)*dps_b(l,m)
        enddo
      enddo
      br_b(j,k) = summb1*(R0     /r_b)**2.0
      bt_b(j,k) = summb2*(R0**2.0/r_b)
      bp_b(j,k) = summb3*(R0**2.0/r_b)

        enddo
      enddo
      
      end subroutine insulating_bvalues
! ===========================================================================

      subroutine insulating_constants(t,p,dt,dp)

        use ModHyperGeometric, ONLY: factorial!, lgndrp
        use ModPhysics, ONLY: Bdp,DipoleStrengthSi
        use CON_axes,   ONLY: MagAxis_D
        
        real, intent(in) :: t(Nt)
        real, intent(in) :: p(Np)
        real, intent(in) :: dt(Nt)
        real, intent(in) :: dp(Np)
        
        integer :: l, m, j, k
        real :: klm, costh, temp
       
        do l=1,Nl
           do m=0,l
              Glm(l,m)=0.0
              Hlm(l,m)=0.0
           enddo
        enddo

        ! Due to the separation of B0 and B1 in BATSRUS, we only calculate
        ! the boundary condition for B1 here, so there's no need for including
        ! the inner dipole field.
        ! If used, the discretization may introduce weird solution.
        !Glm(1,0)=Bdp*MagAxis_D(3)
        !Glm(1,1)=Bdp*MagAxis_D(1)
        !Hlm(1,1)=Bdp*MagAxis_D(2)
        
        ! If higher degrees of internal fields are used add them here.
        
        ! ---------------------------------- for Step 1 -----------------------
        do l=1,Nl
           do m=0,l
              ! constant K_lm
              klm=(2*l+1)*factorial(l-m)/(4*cPi*R0*l*(l+1)*factorial(l+m))
              if (m .gt. 0) then
                 klm=klm*2.0
              endif
              
              klm = klm * R0
              
              do j=1,Nt
                 ! theta part
                 costh = cos(t(j))
                 temp = klm*lgndrp(l,m,costh)*sin(t(j))*dt(j)
                 do k=1,Np
                    ! phi part
                    lmc(j,k,l,m)=temp*cos(m*p(k))*dp(k)
                    lms(j,k,l,m)=temp*sin(m*p(k))*dp(k)
                 enddo
              enddo
              
           enddo
        enddo
        
        ! ---------------------------------- for Step 3 ----------------------
        ! ------- for B_r ------
        do l=1,Nl
           do m=0,l
              klm = l*(l+1)
              
              do j=1,Nt
                 costh = cos(t(j))
                 temp = klm*lgndrp(l,m,costh)
                 do k=1,Np
                    lmb1c(l,m,j,k)=temp*cos(m*p(k))
                    lmb1s(l,m,j,k)=temp*sin(m*p(k))
                 enddo
              enddo
              
           enddo
        enddo
        
      ! ------- for B_theta -------
        do l=1,Nl
           do m=0,l
              
              do j=1,Nt
                 costh = cos(t(j))
                 ! Derivation of associated legendre polynomials
                 if (m .eq. 0) then
                    temp = lgndrp(l,1,costh)
                 elseif (m .eq. l) then
                    temp = -0.5*(l+l)*lgndrp(l,m-1,costh)
                 else
                    temp = -0.5*((l+m)*(l-m+1)*lgndrp(l,m-1,costh) &
                         -lgndrp(l,m+1,costh))
                 endif
                 do k=1,Np
                    lmb2c(l,m,j,k)=temp*cos(m*p(k))
                    lmb2s(l,m,j,k)=temp*sin(m*p(k))
                 enddo
              enddo
              
           enddo
        enddo
        
        ! ------- for B_phi ------- 
        do l=1,Nl
           do m=0,l
              klm = m
              
              do j=1,Nt
                 costh = cos(t(j))
                 temp = klm*lgndrp(l,m,costh)/sin(t(j))
                 do k=1,Np
                    lmb3c(l,m,j,k)=-temp*sin(m*p(k))
                    lmb3s(l,m,j,k)=temp*cos(m*p(k))
                 enddo
              enddo
              
           enddo
        enddo
        
      end subroutine insulating_constants

! ===========================================================================
! ---------------------------------------------------------------------------
!                                   LGNDRP
! ---------------------------------------------------------------------------
! Calculates associated legendre polynomials 
! (recursion algorithm of 'Numerical Recipes')
! ---------------------------------------------------------------------------
      FUNCTION lgndrp(l,m,x)
      
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: l,m
      REAL, INTENT(IN) :: x
      REAL :: lgndrp

      REAL :: p0,p1,p2,m_eins

      INTEGER :: i

      IF(m<0 .OR. l<0 .OR. m>l .OR. abs(x)>1.0D0) THEN
            WRITE(*,*) l,m,x
            call stop_mpi("Lgndrp: wrong arguments!")
      ENDIF

      !P_mm, (l=m) 
      p0=1.0

      m_eins=-1.0*SQRT(1.0-x*x)
      DO i=1,m ! only m>0
         p0=p0*m_eins*(2*i-1.0)
      ENDDO

      IF(l==m) THEN
         lgndrp=p0
         RETURN
      ENDIF

      ! P_ml, l=m+1
      p1=x*(2*m+1)*p0

      IF(l==m+1) THEN
         lgndrp=p1
         RETURN
      ENDIF

      !P_ml, l>m+1
      DO i=m+2,l
         p2=(x*(2*i-1)*p1-(i+m-1)*p0)/(i-m)
         p0=p1
         p1=p2
      ENDDO

      lgndrp=p2
      RETURN
      END FUNCTION lgndrp
      
    end subroutine user_set_cell_boundary
!==============================================================================

    subroutine user_action(NameAction)

      use ModMain,       ONLY: nBlock
      use ModGeometry,   ONLY: Rmin_BLK
      
      character(len=*), intent(in):: NameAction

      integer :: iTypeProc ! Proc type: 1 if inner BC contained; 0 otherwise
      logical :: DoTest
      character(len=*), parameter:: NameSub = 'user_action'
      !------------------------------------------------------------------------
      call test_start(NameSub, DoTest)

      select case(NameAction)
      case('initial condition done','loadbalancedone')
         iTypeProc = 0
         if(minval(Rmin_BLK(1:nBlock)) < 1.1*PlanetRadius) &
              iTypeProc = 1
         
         ! Split the communicator based on inner boundary processor
         call MPI_Comm_split(iComm,iTypeProc,iProc,iCommInner,iError)
         call MPI_Comm_rank(iCommInner,iProcInner,iError)
         call MPI_Comm_size(iCommInner,nProcInner,iError)
      end select
      
      call test_stop(NameSub, DoTest)
    end subroutine user_action
    !==========================================================================
    
end module ModUser
!==============================================================================
