!This code is a copyright protected software (c) 2002- University of Michigan
module ModPIC

  ! Variables and methods for coupling BATSRUS with a PIC code

  implicit none

  SAVE

  private ! except

  public:: pic_read_param
  public:: pic_save_region
  public:: pic_update_states
  public:: pic_param

  logical, public:: UsePic = .false.

  ! Local variables
  integer:: nGhostPic   = 3
  integer:: nOverlapPic = 0
  real   :: TimeUnitPic = 1.0
  real   :: xUnitPicSi  = 1.0
  real   :: uUnitPicSi  = 1.0
  real   :: mUnitPicSi  = 1.0

  character(len=100):: NameFilePic = 'GM/IO2/ipic3d.dat'

  integer:: nRegionPic = 0
  real, allocatable:: XyzMinPic_DI(:,:), XyzMaxPic_DI(:,:), DxyzPic_DI(:,:)

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
       call read_var('TimeUnitPic', TimeUnitPic)

    case("#PICUNIT")
       call read_var('xUnitPicSi', xUnitPicSi)
       call read_var('uUnitPicSi', uUnitPicSi)
       call read_var('mUnitPicSi', mUnitPicSi)

    case("#PICREGION")
       call read_var('nPicRegion', nRegionPic)
       if(allocated(XyzMinPic_DI)) deallocate( &
            XyzMinPic_DI, XyzMaxPic_DI, DxyzPic_DI)
       allocate( &
            XyzMinPic_DI(nDim,nRegionPic), &
            XyzMaxPic_DI(nDim,nRegionPic), &
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
       end do
    case default
       if(iProc==0) call stop_mpi(NameSub//': unknown command='//NameCommand)
    end select

  end subroutine pic_read_param
  !===========================================================================
  subroutine pic_save_region

    use ModProcMH,  ONLY: iProc, nProc, iComm
    use ModAdvance, ONLY: Rho_, RhoUx_, RhoUz_, Bx_, Bz_, p_, State_VGB
    use ModMain,    ONLY: UseB0, Dt, time_simulation, n_step
    use ModB0,      ONLY: B0_DGB
    use ModCurrent, ONLY: get_current
    use ModMpi,     ONLY: MPI_reduce, MPI_SUM, MPI_REAL
    use BATL_lib,   ONLY: nDim, MaxDim, Xyz_DGB, CellSize_DB, find_grid_block
    use ModIO,      ONLY: NamePlotDir
    use ModPhysics, ONLY: No2Si_V, UnitT_, UnitX_, UnitRho_, UnitU_, UnitB_, &
         UnitP_, UnitJ_
    use ModPlotFile,ONLY: save_plot_file

    ! Assuming ideal MHD for now !!! Add Pe, Ppar, PePar multi-ion???
    integer, parameter:: RhoPic_=1, UxPic_=2, UzPic_=4, BxPic_=5, BzPic_=7, &
         pPic_=8, JxPic_=9, JzPic_=11, nVarPic = 11

    ! Coordinate, variable and parameter names
    character(len=*), parameter:: NameVarPic = &
         'x y rho ux uy uz bx by bz p jx jy jz dt xUnitPic uUnitPic mUnitPic'

    ! PIC grid indexes
    integer:: iPic, jPic, kPic, nPic_D(MaxDim), iRegion

    ! MHD grid indexes
    integer:: i, j, k, iBlock, iProcFound

    ! Cell indexes in an array
    integer:: iCell_D(MaxDim)

    ! Location of PIC node
    real:: XyzPic_D(MaxDim)

    ! Current in normalized units
    real:: Current_D(MaxDim)

    ! The PIC variable array
    real, allocatable:: StatePic_VC(:,:,:,:), StateAllPic_VC(:,:,:,:)

    ! Time step in SI units
    real:: DtSi

    ! MPI error
    integer:: iError

    character(len=100):: NameFile

    character(len=*), parameter:: NameSub = 'pic_save_region'
    !-------------------------------------------------------------------------
    DtSi = dt*No2Si_V(UnitT_)
    XyzPic_D = 0.0
    nPic_D = 1
    do iRegion = 1, nRegionPic

       nPic_D(1:nDim) = nint( &
            (XyzMaxPic_DI(:,iRegion) - XyzMinPic_DI(:,iRegion)) &
            / DxyzPic_DI(:,iRegion) )

       allocate(StatePic_VC(nVarPic, nPic_D(1), nPic_D(2), nPic_D(3)))
       StatePic_VC = 0.0

       do kPic = 1, nPic_D(3); do jPic = 1, nPic_D(2); do iPic = 1, nPic_D(1)

          ! Set location of PIC node
          iCell_D = (/iPic, jPic, kPic/)
          XyzPic_D(1:nDim) = XyzMinPic_DI(:,iRegion) &
               + (iCell_D(1:nDim) - 0.5)*DxyzPic_DI(:,iRegion)

          call find_grid_block(XyzPic_D, iProcFound, iBlock, iCell_D)
          if(iProcFound /= iProc) CYCLE

          ! Find corresponding MHD cell center
          i = iCell_D(1); j = iCell_D(2); k = iCell_D(3)

          ! For now we only support coinciding grids...
          if(any(abs(Xyz_DGB(:,i,j,k,iBlock) - XyzPic_D) > &
               1e-5*CellSize_DB(:,iBlock))) then
             write(*,*)'ERROR in ',NameSub
             write(*,*)'XyzPic_D = ', XyzPic_D
             write(*,*)'XyzMhd_D = ', Xyz_DGB(:,i,j,k,iBlock)
             write(*,*)'iProc, iBlock, i, j, k=', iProc, iBlock, i, j, k
             call stop_mpi('PIC cell center does not match MHD grid!')
          end if

          StatePic_VC(RhoPic_,iPic,jPic,kPic) = No2Si_V(UnitRho_)* &
               State_VGB(Rho_,i,j,k,iBlock)
          StatePic_VC(UxPic_:UzPic_,iPic,jPic,kPic) = No2Si_V(UnitU_)* &
               State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) &
               / StatePic_VC(RhoPic_,iPic,jPic,kPic)
          if(UseB0)then
             StatePic_VC(BxPic_:BzPic_,iPic,jPic,kPic) = No2Si_V(UnitB_)* &
                  (State_VGB(Bx_:Bz_,i,j,k,iBlock) + B0_DGB(:,i,j,k,iBlock))
          else
             StatePic_VC(BxPic_:BzPic_,iPic,jPic,kPic) = No2Si_V(UnitB_)* &
                  State_VGB(Bx_:Bz_,i,j,k,iBlock)
          end if
          StatePic_VC(pPic_,iPic,jPic,kPic) = No2Si_V(UnitP_)* &
               State_VGB(p_,i,j,k,iBlock)

          ! Put current into the last three elemets
          call get_current(i, j, k, iBlock, Current_D)
          StatePic_VC(JxPic_:JzPic_,iPic,jPic,kPic) = No2Si_V(UnitJ_)*Current_D

       end do; end do; end do

       if(nProc > 1)then
          ! Collect information from all processors
          allocate(StateAllPic_VC(nVarPic,nPic_D(1),nPic_D(2),nPic_D(3)))
          call MPI_reduce(StatePic_VC, StateAllPic_VC, size(StatePic_VC), &
               MPI_REAL, MPI_SUM, 0, iComm, iError)
          if(iProc == 0) StatePic_VC = StateAllPic_VC
          deallocate(StateAllPic_VC)
       endif

       if(iProc == 0)then
          write(NameFile,'(a,i2.2,a)') &
               trim(NamePlotDir)//'mhd_to_pic_',iRegion,'.dat'
          call save_plot_file(NameFile, &
               StringHeaderIn='PIC region in SI units', &
               nStepIn = n_step, TimeIn = time_simulation, &
               ParamIn_I = (/ DtSi, xUnitPicSi, uUnitPicSi, mUnitPicSi /), &
               NameVarIn = NameVarPic, &
               nDimIn = nDim, &
               CoordMinIn_D = No2Si_V(UnitX_)* &
               (XyzMinPic_DI(:,iRegion) + 0.5*DxyzPic_DI(:,iRegion)), &
               CoordMaxIn_D = No2Si_V(UnitX_)* &
               (XyzMaxPic_DI(:,iRegion) - 0.5*DxyzPic_DI(:,iRegion)), &
               VarIn_VIII = StatePic_VC)
       end if

       deallocate(StatePic_VC)

    end do ! iRegion

  end subroutine pic_save_region
  !===========================================================================
  real function pic_param(NameParam)

    use ModProcMH, ONLY: iProc

    character(len=*), intent(in):: NameParam

    character(len=*), parameter:: NameSub = 'pic_param'
    !------------------------------------------------------------------------
    select case(NameParam)
    case('tunitpic')
       pic_param = TimeUnitPic
    case('noverlap')
       pic_param = nOverlapPic
    case('nghostpic')
       pic_param = nGhostPic
    case default
       if(iProc==0)call stop_mpi(NameSub//': unknown NameParam='//NameParam)
    end select

  end function pic_param

  !===========================================================================

  subroutine pic_update_states(iBlock)

    ! Overwrite the PIC region with the PIC solution

    use ModAdvance,    ONLY: State_VGB
    use ModVarIndexes

    use ModSize, ONLY: nDim, MaxDim, nI, nJ, nK
    use ModProcMH, ONLY: iProc
    use ModMain, ONLY: n_step, time_simulation
    use BATL_lib, ONLY: Xyz_DGB, CellSize_DB

    use ModEnergy,      ONLY: calc_energy_cell
    use ModPlotFile,    ONLY: read_plot_file
    use ModInterpolate, ONLY: bilinear
    use ModUtilities,   ONLY: sleep
    use ModIoUnit,      ONLY: UnitTmp_

    integer, intent(in):: iBlock

    integer:: i, j, k, iError
    real:: XyzNorm_D(MaxDim)

    ! PIC coupling related local variables
    integer:: nStepLast = -1000

    ! PIC grid size
    integer, save:: nXPic, nYPic, nCellPic_D(nDim)
    real,    save:: CoordMinPic_D(nDim), CoordMaxPic_D(nDim), DxyzPic_D(nDim)

    ! PIC variables
    integer, save:: nVarPic
    real,    save, allocatable:: StatePic_VC(:,:,:), StatePic_V(:)

    integer:: Dn
    real:: WeightMhd, WeightPic

    character(len=*), parameter :: NameSub = 'pic_update_states'
    !--------------------------------------------------------------------------

    ! There is no PIC solution at the very beginning
    if(time_simulation == 0.0) RETURN

    ! Check if we should read in a new PIC file
    if(n_step > nStepLast)then
       nStepLast  = n_step

       if(iProc == 0)write(*,*) NameSub,' trying to read ', NameFilePic
       ! Wait until file exists
       do
          open(UnitTmp_, FILE=NameFilePic, STATUS='OLD', IOSTAT=iError)
          if(iError /= 0)then
             ! If not successful, wait a bit and try open again
             call sleep(0.1)
             CYCLE
          end if
          close(UnitTmp_)
          EXIT
       end do

       ! Check if this is the first time 
       if(.not.allocated(StatePic_VC))then
          ! Get size of PIC grid and allocate array
          do
             call read_plot_file(NameFilePic, &
                  nVarOut = nVarPic, nOut_D = nCellPic_D, iErrorOut=iError)
             if(iError /= 0)then
                write(*,*) NameSub,': could not read header from ', &
                     trim(NameFilePic), ' on processor', iProc
                call sleep(0.5)
                CYCLE
             end if
             EXIT
          end do
          nXPic = nCellPic_D(1)
          nYPic = nCellPic_D(min(nDim,2))
          allocate(StatePic_VC(nVarPic,nXPic,nYPic), StatePic_V(nVarPic))
          if(iProc == 0)write(*,*) NameSub, &
               ' allocated StatePic_VC with nXPic, nYPic=', nXPic, nYPic

          ! Read first PIC data and coordinate limits
          do
             call read_plot_file(NameFilePic, VarOut_VII=StatePic_VC, &
                  CoordMinOut_D = CoordMinPic_D, &
                  CoordMaxOut_D = CoordMaxPic_D, &
                  iErrorOut=iError)
             if(iError /= 0)then
                write(*,*) NameSub,': could not read first data from ', &
                     trim(NameFilePic), ' on processor', iProc
                call sleep(0.5)
                CYCLE
             end if
             EXIT
          end do

          DxyzPic_D = (CoordMaxPic_D - CoordMinPic_D)/(nCellPic_D - 1)

          if(iProc == 0)then
             write(*,*) NameSub, ' CoordMinPic_D=', CoordMinPic_D
             write(*,*) NameSub, ' CoordMaxPic_D=', CoordMaxPic_D
             write(*,*) NameSub, ' DxyzPic_D    =', DxyzPic_D
          end if
       else
          ! Read in PIC data
          do
             call read_plot_file(NameFilePic, VarOut_VII=StatePic_VC, &
                  iErrorOut=iError)
             if(iError /= 0)then
                write(*,*) NameSub,': could not read data from ', &
                     trim(NameFilePic), ' on processor', iProc
                call sleep(0.5)
                CYCLE
             end if
             EXIT
          end do

       end if

       ! We reuse the same filename, so delete file after it was read in
       call barrier_mpi
       if(iProc ==0 ) then
          open(UnitTmp_, FILE=NameFilePic, STATUS='OLD', IOSTAT=iError)
          close(UnitTmp_,STATUS='DELETE')
       end if

    end if

    ! Overwrite cells inside the PIC domain
    ! Check if block overlaps with PIC domain (Cartesian ONLY for now!)
!!!    if(any(Xyz_DGB(1:nDim, 1, 1, 1,iBlock) > CoordMaxPic_D + 0.1*CellSize_DB(1:nDim,iBlock)) &
!!!         .or. &
!!!       any(Xyz_DGB(1:nDim,nI,nJ,nK,iBlock) < CoordMinPic_D  0.1*CellSize_DB(1:nDim,iBlock)))&
!!!       RETURN

    if(  all(Xyz_DGB(1:nDim, 1, 1, 1,iBlock) <= CoordMaxPic_D + 0.1*CellSize_DB(1:nDim,iBlock)) .and. &
         all(Xyz_DGB(1:nDim,nI,nJ,nK,iBlock) >= CoordMinPic_D - 0.1*CellSize_DB(1:nDim,iBlock))) then


    do k = 1, nK; do j = 1, nJ; do i = 1, nI

       ! Normalized PIC grid coordinates (1...nCellPic_D)
       XyzNorm_D(1:nDim) = 1 + (Xyz_DGB(1:nDim,i,j,k,iBlock) - CoordMinPic_D) &
            /DxyzPic_D

       ! Distance from edge
       Dn = minval( min(nint(XyzNorm_D(1:nDim) - 1),  &
            nint(nCellPic_D - XyzNorm_D(1:nDim))) )

       ! Nothing to do within PIC ghost region
       if(Dn < nGhostPic) CYCLE

       ! Distance from ghost layers
       Dn = Dn - nGhostPic + 1

       if(Dn <= nOverlapPic)then
          ! For nOverlapPic=1, Dn = 1, so use 0.5 as weight
          ! For nOverlapPic=2, Dn = 1, 2, so use 1/3 and 2/3 weights
          ! ...
          WeightPic = Dn/(nOverlapPic + 1.0)
       else
          WeightPic = 1.0
       end if

       WeightMhd = 1.0 - WeightPic

       StatePic_V = &
            bilinear(StatePic_VC, nVarPic, 1, nXPic, 1, nYPic, XyzNorm_D)

       ! Convert velocity to momentum
       StatePic_V(RhoUx_:RhoUz_) = StatePic_V(Rho_)*StatePic_V(Ux_:Uz_)

       ! Interpolate MHD and PIC states. 
       ! Skip hyperbolic scalar if present (Hyp=Bz_+1)
       State_VGB(Rho_:Bz_,i,j,k,iBlock) = &
            WeightMhd*State_VGB(Rho_:Bz_,i,j,k,iBlock) + &
            WeightPic*StatePic_V(Rho_:Bz_)

       State_VGB(p_,i,j,k,iBlock) = &
            WeightMhd*State_VGB(p_,i,j,k,iBlock) &
            + WeightPic*StatePic_V(nVarPic)

       ! Set hyperbolic scalar to zero if present
       if(Hyp_>1) State_VGB(Hyp_,i,j,k,iBlock) = 0.0

    end do; end do; end do

    call calc_energy_cell(iBlock)

    end if !!!

  end subroutine pic_update_states

end module ModPIC
