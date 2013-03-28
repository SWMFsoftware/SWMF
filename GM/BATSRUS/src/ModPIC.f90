!This code is a copyright protected software (c) 2002- University of Michigan
module ModPIC

  ! Variables and methods for coupling BATSRUS with a PIC code

  implicit none

  SAVE

  private ! except

  public:: pic_read_param
  public:: pic_update_states
  public:: pic_param

  logical, public:: UsePic = .false.

  ! Local variables
  integer:: nGhostPic   = 3
  integer:: nOverlapPic = 0
  real   :: TimeUnitPic = 1.0
  character(len=100):: NameFilePic = 'GM/IO2/ipic3d.dat'

contains
  !===========================================================================
  subroutine pic_read_param(NameCommand)

    use ModProcMH,    ONLY: iProc
    use ModReadParam, ONLY: read_var

    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'read_pic_param'
    !------------------------------------------------------------------------
    select case(NameCommand)
    case("#PIC")
       call read_var('UsePic',      UsePic)
       call read_var('nGhostPic',   nGhostPic)
       call read_var('nOverlapPic', nOverlapPic)
       call read_var('TimeUnitPic', TimeUnitPic)
    case default
       if(iProc==0) call stop_mpi(NameSub//': unknown command='//NameCommand)
    end select

  end subroutine pic_read_param

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
