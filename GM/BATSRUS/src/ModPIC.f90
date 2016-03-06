!  Copyright (C) 2002 Regents of the University of Michigan
!  Portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModPIC

  use ModMain,      ONLY: UseB0
  use ModB0,        ONLY: B0_DGB
  ! Variables and methods for coupling BATSRUS with a PIC code

  implicit none

  SAVE

  private ! except

  public:: pic_read_param
  public:: pic_init_region
  public:: pic_find_node
  
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
  character(len=100):: NameFilePic = 'GM/IO2/ipic3d.dat'

  ! PIC regions
  integer, public :: nRegionPic = 0
  real, public, allocatable:: XyzMinPic_DI(:,:), XyzMaxPic_DI(:,:), DxyzPic_DI(:,:)
  real, public, allocatable:: XyzPic0_DI(:,:)

  ! Is the node overlaped with PIC domain?
  logical, public,allocatable:: IsPicNode_A(:)

  logical :: IsPicRegionInitialized = .false.
contains
  !===========================================================================
  subroutine pic_read_param(NameCommand)

    use ModProcMH,    ONLY: iProc
    use ModReadParam, ONLY: read_var
    use BATL_lib,     ONLY: x_, y_, z_, nDim

    character(len=*), intent(in):: NameCommand

    integer:: iRegion

    character(len=*), parameter:: NameSub = 'read_pic_param'
    !------------------------------------------------------------------------
    select case(NameCommand)
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
       end do
    case default
       if(iProc==0) call stop_mpi(NameSub//': unknown command='//NameCommand)
    end select

  end subroutine pic_read_param

  !===========================================================================
  subroutine pic_init_region

    use ModProcMH,    ONLY: iProc
    use ModMain,      ONLY: Dt
    use BATL_lib,     ONLY: nDim, MaxDim, Xyz_DGB, CellSize_DB, find_grid_block
    use ModPhysics,   ONLY: No2Si_V, UnitT_, UnitMass_, UnitCharge_
    use ModHallResist,ONLY: HallFactorMax
    use ModPhysics,   ONLY: IonMassPerCharge
    use ModMultiFluid,ONLY: nIonFLuid, MassIon_I

    ! Assuming ideal/aniso MHD for now !!! Add Pe, PePar multi-ion???
    integer, parameter:: RhoPic_=1, UxPic_=2, UzPic_=4, BxPic_=5, BzPic_=7, &
         PparPic_=8, pPic_=9, JxPic_=10, JzPic_=12, nVarPic = 12

    ! Coordinate, variable and parameter names
    character(len=*), parameter:: NameVarPic = &
         'x y rho ux uy uz bx by bz ppar p jx jy jz dt xUnitPic uUnitPic mUnitPic'

    ! PIC grid indexes
    integer:: iPic, jPic, kPic, nPic_D(MaxDim), iRegion

    ! MHD grid indexes
    integer:: i, j, k, iBlock, iProcFound

    ! Cell indexes in an array
    integer:: iCell_D(MaxDim)

    ! Location of PIC node
    real:: XyzPic_D(MaxDim)

    ! The PIC variable array
    real, allocatable:: StatePic_VC(:,:,:,:)

    ! mass per charge SI
    real:: IonMassPerChargeSi 

    ! Region check
    integer :: nCell

    character(len=*), parameter:: NameSub = 'pic_init_region'
    !-------------------------------------------------------------------------

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

    XyzPic_D = 0.0
    nPic_D = 1

    do iRegion = 1, nRegionPic
       ! extending the region sizes with 1 ghost cell
       do i=1, nDim
          ! Number of nodes (not cells!)
          nCell = nint( &
               (XyzMaxPic_DI(i,iRegion) - XyzMinPic_DI(i,iRegion)) &
               /DxyzPic_DI(i,iRegion) )

          ! Adding 1 ghost cell layer at max and 1 at min
          XyzMaxPic_DI(i,iRegion) = XyzPic0_DI(i,iRegion) &
               + (nCell + 1)*DxyzPic_DI(i,iRegion)
          XyzMinPic_DI(i,iRegion) = XyzPic0_DI(i,iRegion) &
               - DxyzPic_DI(i,iRegion)
       end do
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
  end subroutine pic_init_region
  !===========================================================================

  subroutine pic_find_node()
    ! Find out the blocks that overlaped with PIC region(s). 
    use ModProcMH, ONLY: iProc
    use BATL_lib,  ONLY: nDim, MaxDim, nBlock, MaxBlock, find_grid_block, &
         x_, y_, z_, iTree_IA, MaxNode, Block_
    integer:: nIJK_D(1:MaxDim), ijk_D(1:MaxDim)
    real:: Xyz_D(1:MaxDim)
    integer:: iRegion, iBlock, i, j, k, iProcFound, iNode

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'pic_find_node'
    !----------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*) NameSub,' is called'

    if(.not.allocated(IsPicNode_A)) allocate(IsPicNode_A(MaxNode))
    IsPicNode_A = .false.
    nIJK_D = 1; Xyz_D = 0
    
    do iRegion = 1, nRegionPic
       ! nIJK_D includes two ghost cells. 
       nIJK_D(1:nDim) = int(&
            (XyzMaxPic_DI(1:nDim,iRegion)- XyzMinPic_DI(1:nDim,iRegion)) &
            /DxyzPic_DI(1:nDim,iRegion) + 0.5)

       if(DoTestMe) write(*,*) NameSub,' iRegion = ',iRegion, &
            ' nIJK_D = ',nIJK_D(1:nDim)

       do k=1, nIJK_D(z_); do j=1, nIJK_D(y_); do i=1, nIJK_D(x_)
          ! Loop through all the PIC node points. 
          ijk_D(x_) = i - 1; ijk_D(y_) = j - 1; ijk_D(z_) = k - 1
          Xyz_D(1:nDim) = XyzMinPic_DI(1:nDim,iRegion) &
               + ijk_D(1:nDim)*DxyzPic_DI(1:nDim,iRegion)
          call find_grid_block(Xyz_D, iProcFound, iBlock, iNodeOut=iNode)
          IsPicNode_A(iNode) = .true.
       enddo; enddo; enddo
       
    enddo
    
    if(DoTestMe) write(*,*)'IsPicNode= ', IsPicNode_A(:)
  end subroutine pic_find_node
  !===========================================================================
end module ModPIC
