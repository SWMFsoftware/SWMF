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

  logical, public:: UsePic = .false.

  ! Local variables

  ! Coupling parameters
  integer, public :: nGhostPic   = 1  ! Number of ghost cells around PIC region
  integer:: nOverlapPic = 0  ! Overlap region with linear interpolation
  integer, public :: DnCouplePic = 1  ! Coupling frequency

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
    case("#PIC")
       call read_var('UsePic',      UsePic)
       call read_var('nGhostPic',   nGhostPic)
       call read_var('nOverlapPic', nOverlapPic)

    case("#PICCOUPLE")
       call read_var('DnCouplePic', DnCouplePic)

    case("#PICUNIT")
       call read_var('xUnitPicSi', xUnitPicSi)
       call read_var('uUnitPicSi', uUnitPicSi)

    case("#PICREGION")
       call read_var('nPicRegion', nRegionPic)
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

    ! Time step in SI units
    real:: DtSi

    ! mass per charge SI
    real:: IonMassPerChargeSi 

    ! Region check
    integer :: nNode

    character(len=*), parameter:: NameSub = 'pic_init_region'
    !-------------------------------------------------------------------------

    ! Save first step and then every DnCouplePic steps

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

    IonMassPerChargeSi = IonMassPerCharge* &
         No2Si_V(UnitMass_)/No2Si_V(UnitCharge_)

    ! if(nIonFluid == 1) IonMassPerChargeSi = IonMassPerChargeSi/MassIon_I(1)

    mUnitPicSi = 1e7*xUnitPicSi * (IonMassPerChargeSi*HallFactorMax)**2

    if(iProc==0)then
       write(*,*) NameSub,': IonMassPerChargeSi=', IonMassPerChargeSi
       write(*,*) NameSub,': xUnitPicSi = ',xUnitPicSi
       write(*,*) NameSub,': mUnitPicSi = ',mUnitPicSi
       write(*,*) NameSub,': uUnitPicSi = ',uUnitPicSi
    end if

    DtSi = Dt*No2Si_V(UnitT_)*DnCouplePic
    XyzPic_D = 0.0
    nPic_D = 1

    do iRegion = 1, nRegionPic
       ! checking and corecting that the  domain size and Dx 
       ! match what IPIC3D needs ( odd number of cells)
       do i=1, nDim
          nNode = 1 + nint( &
               (XyzMaxPic_DI(i,iRegion) - XyzMinPic_DI(i,iRegion)) &
               /DxyzPic_DI(i,iRegion) )

          ! Make sure there is an odd number of nodes
          if(mod(nNode, 2) == 0) nNode = nNode + 1

          ! Adding ghost cell layars (1 at Min, 1 or 2 at max)
          XyzMaxPic_DI(i,iRegion) = XyzPic0_DI(i,iRegion) &
               + nNode*DxyzPic_DI(i,iRegion)
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
end module ModPIC
