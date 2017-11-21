!  Copyright (C) 2002 Regents of the University of Michigan
!  Portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModPIC

  use BATL_lib, ONLY: &
       test_start, test_stop

  ! Variables and methods for coupling BATSRUS with a PIC code

  implicit none

  SAVE

  private ! except

  public:: pic_read_param
  public:: pic_init_region
  public:: pic_find_node
  public:: pic_find_region
  public:: mhd_to_pic_vec

  logical, public:: UsePic = .false.

  logical, public:: DoBalancePicBlock=.true.

  ! Local variables

  ! Coupling parameters
  integer, public :: nGhostPic   = 2  ! Number of ghost cells around PIC region

  ! Conversion to PIC units
  real, public  :: xUnitPicSi  = 1.0
  real, public  :: uUnitPicSi  = 1.0
  real, public  :: mUnitPicSi  = 1.0

  ! File sent by the PIC code

  ! PIC regions
  integer, public :: nRegionPic = 0

  ! 1) LenPic_D, XyzMinPic_DI and XyzMaxPic_DI include the ghost cell layers.
  !    An 2D example: XyzMinPic_DI is the bottom left point, XyzMaxPic_DI
  !      is the bottom right point. Note the PIC box can be rotated.
  ! 2) R_DDI: mhd coordinates to pic coordinates.
  real, public, allocatable:: XyzMinPic_DI(:,:), XyzMaxPic_DI(:,:), &
       LenPic_DI(:,:), DxyzPic_DI(:,:), R_DDI(:,:,:)

  ! The PIC regon can rotate around x, y and z axis. The rotation rule is the
  ! same as the rotation for #REGION command.
  ! In the rotated coordinates, XyzMinPic_DI instead of
  ! XyzPic0_DI is the origin point.
  real, public, allocatable:: XyzPic0_DI(:,:)
  logical, public :: DoRotatePIC = .false.

  ! Is the node overlaped with PIC domain?
  logical, public,allocatable:: IsPicNode_A(:)

  logical :: IsPicRegionInitialized = .false.
contains
  !============================================================================
  subroutine pic_read_param(NameCommand)

    use ModProcMH,    ONLY: iProc
    use ModReadParam, ONLY: read_var
    use BATL_lib,     ONLY: x_, y_, z_, nDim
    use ModNumConst,  ONLY: cDegToRad
    use ModCoordTransform, ONLY: rot_matrix_x, rot_matrix_y, rot_matrix_z, &
         rot_matrix

    character(len=*), intent(in):: NameCommand

    integer:: iRegion

    real :: xRotate, yRotate, zRotate

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'pic_read_param'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    select case(NameCommand)
    case ("#PICGHOST")
      call read_var('nGhostPic',nGhostPic)
    case("#PICUNIT")
       call read_var('xUnitPicSi', xUnitPicSi)
       call read_var('uUnitPicSi', uUnitPicSi)
    case("#PICBALANCE")
       call read_var('DoBalancePicBlock', DoBalancePicBlock)
    case("#PICREGION")
       call read_var('nPicRegion', nRegionPic)
       UsePic = nRegionPic > 0
       if(allocated(XyzMinPic_DI)) deallocate( &
            XyzMinPic_DI, XyzMaxPic_DI, DxyzPic_DI,XyzPic0_DI)
       allocate( &
            XyzMinPic_DI(nDim,nRegionPic), &
            XyzMaxPic_DI(nDim,nRegionPic), &
            XyzPic0_DI(nDim,nRegionPic), &
            LenPic_DI(nDim,nRegionPic), &
            R_DDI(3,3,nRegionPic),       &
            DxyzPic_DI(nDim,nRegionPic))
       do iRegion = 1, nRegionPic
          call              read_var('xMinPic', XyzMinPic_DI(x_,iRegion))
          call              read_var('xMaxPic', XyzMaxPic_DI(x_,iRegion))
          if(nDim > 1) call read_var('yMinPic', XyzMinPic_DI(y_,iRegion))
          if(nDim > 1) call read_var('yMaxPic', XyzMaxPic_DI(y_,iRegion))
          if(nDim > 2) call read_var('zMinPic', XyzMinPic_DI(z_,iRegion))
          if(nDim > 2) call read_var('zMaxPic', XyzMaxPic_DI(z_,iRegion))
          call              read_var('DxPic',   DxyzPic_DI(x_,iRegion))
          if(nDim > 1) call read_var('DyPic',   DxyzPic_DI(y_,iRegion))
          if(nDim > 2) call read_var('DzPic',   DxyzPic_DI(z_,iRegion))

          ! Origo point for the IPIC3D grid for this region
          XyzPic0_DI(1:nDim,iRegion) = XyzMinPic_DI(1:nDim,iRegion)

          LenPic_DI(1:nDim, iRegion) = XyzMaxPic_DI(1:nDim,iRegion) - &
               XyzMinPic_DI(1:nDim,iRegion)
          R_DDI(:,:,iRegion) = 0
          R_DDI(1,1,iRegion) = 1
          R_DDI(2,2,iRegion) = 1
          R_DDI(3,3,iRegion) = 1

       end do
    case("#PICREGIONROTATE")       
       DoRotatePIC = .true.
       call read_var('nPicRegion', nRegionPic)
       UsePic = nRegionPic > 0
       if(allocated(XyzMinPic_DI)) deallocate( &
            XyzMinPic_DI, XyzMaxPic_DI, LenPic_DI, DxyzPic_DI,XyzPic0_DI)
       allocate( &
            XyzMinPic_DI(nDim,nRegionPic), &
            XyzMaxPic_DI(nDim,nRegionPic), &
            XyzPic0_DI(nDim,nRegionPic), &
            DxyzPic_DI(nDim,nRegionPic), &
            LenPic_DI(nDim,nRegionPic),  &
            R_DDI(3,3,nRegionPic))
       XyzMinPic_DI = 0
       XyzMaxPic_DI = 0
       XyzPic0_DI = 0
       DxyzPic_DI = 0
       LenPic_DI = 0
       R_DDI=0
       do iRegion = 1, nRegionPic
          call              read_var('xMinPic', XyzMinPic_DI(x_,iRegion))
          call              read_var('xLenPic', LenPic_DI(x_,iRegion))
          if(nDim > 1) call read_var('yMinPic', XyzMinPic_DI(y_,iRegion))
          if(nDim > 1) call read_var('yLenPic', LenPic_DI(y_,iRegion))
          if(nDim > 2) call read_var('zMinPic', XyzMinPic_DI(z_,iRegion))
          if(nDim > 2) call read_var('zLenPic', LenPic_DI(z_,iRegion))
          call              read_var('DxPic',   DxyzPic_DI(x_,iRegion))
          if(nDim > 1) call read_var('DyPic',   DxyzPic_DI(y_,iRegion))
          if(nDim > 2) call read_var('DzPic',   DxyzPic_DI(z_,iRegion))
          if(nDim==2 ) then
             xRotate = 0; yRotate = 0
             call read_var('zRotate', zRotate)
             R_DDI(1:2,1:2,iRegion) = rot_matrix(-zRotate*cDegToRad)
          elseif(nDim==3) then
             call              read_var('xRotate', xRotate)
             call              read_var('yRotate', yRotate)
             call              read_var('zRotate', zRotate)
          endif

          XyzPic0_DI(1:nDim,iRegion) = XyzMinPic_DI(1:nDim,iRegion)

          ! Rotation matrix rotates around X, Y and Z axes in this order
          R_DDI(:,:,iRegion) = matmul( matmul( &
               rot_matrix_z(-zRotate*cDegToRad),&
               rot_matrix_y(-yRotate*cDegToRad)),&
               rot_matrix_x(-xRotate*cDegToRad))
       end do
    case default
       if(iProc==0) call stop_mpi(NameSub//': unknown command='//NameCommand)
    end select

    call test_stop(NameSub, DoTest)
  end subroutine pic_read_param
  !============================================================================

  subroutine pic_init_region

    use ModProcMH,    ONLY: iProc
    use BATL_lib,     ONLY: nDim
    use ModPhysics,   ONLY: No2Si_V, UnitMass_, UnitCharge_
    use ModHallResist, ONLY: HallFactorMax
    use ModPhysics,   ONLY: IonMassPerCharge
    use ModMultiFluid, ONLY: nIonFLuid, MassIon_I
    use ModVarIndexes, ONLY: IsMhd

    ! PIC grid indexes
    integer:: iRegion

    ! MHD grid indexes
    integer:: i

    ! mass per charge SI
    real:: IonMassPerChargeSi

    ! Region check
    integer :: nCell

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'pic_init_region'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    if(nIonFluid > 1 .and. IsMhd) then
       if(iProc == 0) then
          write(*,*) ' '
          write(*,*) "Error!!!!!"
          write(*,*) 'Multi-fluid MHD with total density and total ',&
               'momentum variables is used for  the current simulation ', &
               '(see ModEquation.f90). ', &
               'MHD-EPIC does not support this case now. Please eliminate ', &
               'the total density and total momentum variables and try again.'
          write(*,*) ' '
       endif
       call stop_mpi(NameSub)
    endif

    ! Normalizing the system so q/(mc) == 1 in IPIC3D.
    !
    ! In CGS units the Hall speed is uH_CGS = j/(nq) = c/4pi curlB m/(q rho)
    !
    ! In SI  units the Hall speed is uH_SI  = j/nq = curlB/mu0 m/(q rho)
    !
    ! The CGS dimension of curlB/(uH rho) is
    !   1/[L] * [U]sqrt[RHO]/([U][RHO]) = 1/([L] sqrt[RHO]) = sqrt[L]/sqrt[M]
    ! which can be used to scale the PIC units to true CGS units.
    !
    ! CGS->SI: 1kg=1000g, 1m=100cm, 1T=10000G and mu0 = 4pi*1e-7
    !
    ! 1  = q_PIC/(m_PIC c_pic)
    !    = curlB_PIC /  (4pi uH_PIC rho_PIC)
    !    = curlB_CGS /  (4pi uH_CGS rho_CGS) * sqrt( [M]_CGS / [L]_CGS )
    !    = curlB_SI*100/(4pi uH_SI*100 rho_SI*0.001)
    !                                        * sqrt{ [M]_SI*1000 / ([L]_SI*100)
    !    = 10^3.5/4pi * curlB_SI/(uH_SI rho_SI) * sqrt( [M]_SI / [L]_SI )
    !    = 10^3.5/4pi * mu0_SI*q_SI/m_SI        * sqrt( [M]_SI / [L]_SI )
    !    = 10^(-3.5)  *        q_SI/m_SI        * sqrt( [M]_SI / [L]_SI )
    !
    ! Then we can solve for the mass unit
    !
    !   [M]_SI = 10^7 * [L]_SI * (m_SI/q_SI)^2

    ! iPIC3D does not suppor multi-sessions now. If the pic region related
    ! parameters has been initialized, then skip this subroutine in
    ! the following sessions.
    if(IsPicRegionInitialized) RETURN
    IsPicRegionInitialized = .true.
    IonMassPerChargeSi = IonMassPerCharge* &
         No2Si_V(UnitMass_)/No2Si_V(UnitCharge_)

    if(nIonFluid == 1) IonMassPerChargeSi = IonMassPerChargeSi/MassIon_I(1)

    mUnitPicSi = 1e7*xUnitPicSi * (IonMassPerChargeSi*HallFactorMax)**2

    if(iProc==0)then
       write(*,*) NameSub,': IonMassPerChargeSi=', IonMassPerChargeSi
       write(*,*) NameSub,': xUnitPicSi = ',xUnitPicSi
       write(*,*) NameSub,': mUnitPicSi = ',mUnitPicSi
       write(*,*) NameSub,': uUnitPicSi = ',uUnitPicSi
    end if

    do iRegion = 1, nRegionPic
       ! extending the region sizes with 1 ghost cell
       do i=1, nDim
          ! Number of nodes (not cells! )
          nCell = nint(LenPic_DI(i,iRegion)/DxyzPic_DI(i,iRegion))
          ! 1 ghost cell layer at each side
          LenPic_DI(i,iRegion) = (nCell+2)*DxyzPic_DI(i,iRegion)
       end do

       XyzMinPic_DI(:,iRegion) = XyzPic0_DI(:,iRegion) - DxyzPic_DI(:,iRegion)
    end do

    if(iProc == 0) then
       write(*,*) "Corrected IPIC3D  regions"
       do iRegion = 1, nRegionPic
          write(*,*) "  Region : ", iRegion
          write(*,*) "     Min Cordinate : ", &
               XyzMinPic_DI(1:nDim,iRegion) + DxyzPic_DI(1:nDim,iRegion)
          write(*,*) "     Max Cordinate : ", &
               XyzMaxPic_DI(1:nDim,iRegion) - DxyzPic_DI(1:nDim,iRegion)
       end do
    end if
    call test_stop(NameSub, DoTest)
  end subroutine pic_init_region
  !============================================================================

  subroutine pic_find_node()
    ! Find out the blocks that overlaped with PIC region(s).
    use BATL_lib,  ONLY: nDim, MaxDim, find_grid_block, &
         x_, y_, z_, MaxNode
    integer:: nIJK_D(1:MaxDim), ijk_D(1:MaxDim)
    real:: PIC_D(1:MaxDim), MHD_D(1:MaxDim)
    integer:: iRegion, iBlock, i, j, k, iProcFound, iNode

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'pic_find_node'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(DoTest)write(*,*) NameSub,' is called'

    if(.not.allocated(IsPicNode_A)) allocate(IsPicNode_A(MaxNode))
    IsPicNode_A = .false.

    nIJK_D = 1; PIC_D = 0; MHD_D=0

    do iRegion = 1, nRegionPic
       ! nIJK_D includes two ghost cells.
       nIJK_D(1:nDim) = int(&
            LenPic_DI(1:nDim,iRegion)/DxyzPic_DI(1:nDim,iRegion) + 0.5)

       if(DoTest) write(*,*) NameSub,' iRegion = ',iRegion, &
            ' nIJK_D = ',nIJK_D(1:nDim)

       do k=1, nIJK_D(z_); do j=1, nIJK_D(y_); do i=1, nIJK_D(x_)
          ! Loop through all the PIC node points.
          ijk_D(x_) = i - 1; ijk_D(y_) = j - 1; ijk_D(z_) = k - 1
          PIC_D(1:nDim) = ijk_D(1:nDim)*DxyzPic_DI(1:nDim,iRegion)
          call pic_to_mhd_vec(iRegion,PIC_D,MHD_D)
          call find_grid_block(MHD_D, iProcFound, iBlock, iNodeOut=iNode)
          IsPicNode_A(iNode) = .true.
       enddo; enddo; enddo

    enddo

    if(DoTest) write(*,*)'IsPicNode= ', IsPicNode_A(:)
    call test_stop(NameSub, DoTest)
  end subroutine pic_find_node
  !============================================================================

  integer function pic_find_region(iBlock,i,j,k)
    ! If a cell is inside the PIC region, return 1;
    ! otherwise, return 0;
    use BATL_lib, ONLY: nDim, Xyz_DGB

    integer, intent(in) :: iBlock,i,j,k

    integer:: iStatus
    integer:: iRegion
    real:: xyz_D(nDim), PIC_D(nDim)

    character(len=*), parameter:: NameSub = 'pic_find_region'
    !--------------------------------------------------------------------------
    xyz_D = Xyz_DGB(1:nDim,i,j,k,iBlock)

    iStatus=0
    do iRegion = 1, nRegionPic
       call mhd_to_pic_vec(iRegion, xyz_D, PIC_D)

       if(  all(PIC_D > 0 ).and.&
            all(PIC_D < LenPic_DI(:,iRegion))) & ! Not accurate here. --Yuxi
            iStatus = 1
    enddo

    pic_find_region=iStatus
  end function pic_find_region
  !============================================================================

  subroutine pic_to_mhd_vec(iRegion, CoordIn_D, CoordOut_D, OriginIn_D)
    ! Transfer Pic coordinates to Mhd coordinates. Origin_D
    ! is the origin of the PIC coordinates.

    use BATL_lib, ONLY: nDim

    integer, intent(in) :: iRegion
    real, intent(in)    :: CoordIn_D(nDim)
    real, intent(out)   :: CoordOut_D(nDim)
    real, intent(in), optional :: OriginIn_D(nDim)
    real :: Origin_D(nDim), coord_D(nDim), R_DD(3, 3)

    integer:: iDim, jDim

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'pic_to_mhd_vec'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(DoTest)write(*,*) NameSub,' is called'

    Origin_D = XyzMinPic_DI(:,iRegion)
    if(present(OriginIn_D)) Origin_D = OriginIn_D

    Coord_D = 0

    ! R_pic2mhd = transpose(R_mhd2pic)
    coord_D = CoordIn_D
    do iDim = 1, nDim; do jDim = 1, nDim
       R_DD(iDim,jDim) = R_DDI(jDim,iDim,iRegion)
    enddo; enddo

    CoordOut_D = 0
    do iDim = 1, nDim
       CoordOut_D(iDim) = sum(R_DD(iDim,1:nDim)*coord_D)
    enddo

    CoordOut_D = CoordOut_D + Origin_D

    if(DoTest) then
       write(*,*) 'Origin_D   = ', Origin_D
       write(*,*) 'CoordIn_D  = ', CoordIn_D
       write(*,*) 'CoordOut_D = ', CoordOut_D
       do iDim = 1, nDim
          write(*,*) 'iDim = ', iDim, 'R = ', R_DD(iDim,:)
       enddo
    endif

    call test_stop(NameSub, DoTest)
  end subroutine pic_to_mhd_vec
  !============================================================================

  subroutine mhd_to_pic_vec(iRegion, CoordIn_D, CoordOut_D, OriginIn_D)
    ! DoMhd2Pic == true: transfer Mhd coordinates to a coordinates
    !   that is parallel to the PIC coordinates but the origin point is
    !   defined by Origin_D.

    use BATL_lib, ONLY: nDim

    integer, intent(in) :: iRegion
    real, intent(in)    :: CoordIn_D(nDim)
    real, intent(out)   :: CoordOut_D(nDim)
    real, intent(in), optional :: OriginIn_D(nDim)
    real :: Origin_D(nDim), coord_D(nDim), R_DD(3, 3)

    integer :: iDim

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'mhd_to_pic_vec'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(DoTest)write(*,*) NameSub,' is called'

    Origin_D = XyzMinPic_DI(:,iRegion)
    if(present(OriginIn_D)) Origin_D = OriginIn_D

    Coord_D = CoordIn_D - Origin_D

    CoordOut_D = 0
    do iDim = 1, nDim
       CoordOut_D(iDim) = sum(R_DDI(iDim,1:nDim,iRegion)*coord_D)
    enddo

    if(DoTest) then
       write(*,*) 'Origin_D   = ', Origin_D
       write(*,*) 'CoordIn_D  = ', CoordIn_D
       write(*,*) 'CoordOut_D = ', CoordOut_D
       do iDim = 1, nDim
          write(*,*) 'iDim = ', iDim, 'R = ', R_DD(iDim,:)
       enddo
    endif

    call test_stop(NameSub, DoTest)
  end subroutine mhd_to_pic_vec
  !============================================================================

end module ModPIC
!==============================================================================
