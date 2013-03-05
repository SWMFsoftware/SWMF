!#NOTPUBLIC  email:gtoth@umich.edu  expires:12/31/2099
!This code is a copyright protected software (c) 2002- University of Michigan
!=============================================================================
module ModUser

  use ModUserEmpty,               &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_set_ics,                    &
       IMPLEMENTED3 => user_get_log_var,                &
       IMPLEMENTED4 => user_update_states

  use ModNumConst, ONLY: cTwoPi

  include 'user_module.h' !list of public methods

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'GEM and PIC coupling, G. Toth and L. Daldorff'

  ! GEM challenge parameters
  real:: Tp=0.01           ! plasma temperature
  real:: B0=0.0014         ! Background field
  real:: Lambda0=0.5       ! Width of current sheet
  real:: Apert = 0.2       ! amplitude of perturbation
  real:: GaussXInv = 2.0   ! X size of Gaussian perturbation in Az
  real:: GaussYInv = 2.0   ! Y size of Gaussian perturbation in Az
  real:: Kx = cTwoPi/25.6  ! X wave number of perturbation
  real:: Ky = cTwoPi/12.8  ! Y wave number of perturbation

  ! PIC coupling related variables
  integer:: nGhostPic   = 3
  integer:: nOverlapPic = 0
  character(len=100):: NameFilePic

contains

  subroutine user_read_inputs

    use ModProcMH,    ONLY: iProc
    use ModReadParam
    use ModMain,      ONLY: UseUserUpdateStates, UseUserIcs, UseUserLogFiles

    real:: WaveLengthX, WaveLengthY, GaussX, GaussY

    character(len=100) :: NameCommand
    !-------------------------------------------------------------------------
    UseUserLogFiles = .true.
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case('#GEM')
          UseUserIcs = .true.
          call read_var('Amplitude', Apert)

       case('#GEMPARAM')
          call read_var('B0', B0)
          call read_var('Tp', Tp)
          call read_var('CurrentSheetWidth', Lambda0)

       case('#GEMPERTURB')
          call read_var('GaussWidthX', GaussX)
          call read_var('GaussWidthY', GaussY)
          call read_var('WaveLengthX', WaveLengthX)
          call read_var('WaveLengthY', WaveLengthY)

          if(GaussX <= 0)then
             GaussXInv = 0.0
          else
             GaussXInv = 1.0/GaussX
          end if

          if(GaussY <= 0)then
             GaussYInv = 0.0
          else
             GaussYInv = 1.0/GaussY
          end if

          if(WaveLengthX <= 0.0)then
             Kx = 0.0
          else
             Kx = cTwoPi/WaveLengthX
          end if
          if(WaveLengthY <= 0.0)then
             Ky = 0.0
          else
             Ky = cTwoPi/WaveLengthY
          end if
       case('#PIC')
          UseUserUpdateStates = .true.
          call read_var('nGhostPic'  , nGhostPic)
          call read_var('nOverlapPic', nOverlapPic)
          call read_var('NameFilePic', NameFilePic)
       case('#USERINPUTEND')
          if(iProc==0) write(*,*)'USERINPUTEND'
          EXIT

       case default
          if(iProc==0) call stop_mpi( &
               'read_inputs: unrecognized command: '//NameCommand)
       end select
    end do
  end subroutine user_read_inputs
  !============================================================================
  subroutine user_set_ics(iBlock)

    use ModGeometry, ONLY: Xyz_DGB

    use ModPhysics,  ONLY: ShockLeftState_V

    use ModAdvance,  ONLY: State_VGB, Bx_, By_, rho_, Ppar_, p_, Pe_, &
         UseElectronPressure, UseAnisoPressure
    use ModSize,     ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK

    integer, intent(in) :: iBlock

    real                :: x, y, a
    integer             :: i, j, k

    character(len=*), parameter :: NameSub = 'user_set_ics'
    !--------------------------------------------------------------------------

    !write(*,*)'GEM problem set up'
    State_VGB(Bx_,:,:,:,iBlock) = B0*tanh(Xyz_DGB(y_,:,:,:,iBlock)/Lambda0)

    ! Modify pressure(s) to balance magnetic pressure
    if(UseElectronPressure) then
       ! Distribute the correction proportionally between electrons and ions
       State_VGB(Pe_,:,:,:,iBlock) = ShockLeftState_V(Pe_)*(1.0 &
            + 0.5*(B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2)      &
            /(ShockLeftState_V(Pe_) + ShockLeftState_V(p_)))

       State_VGB(p_,:,:,:,iBlock) = ShockLeftState_V(p_)*(1.0   &
            + 0.5*(B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2)      &
            /(ShockLeftState_V(Pe_) + ShockLeftState_V(p_)))
    else
       State_VGB(p_,:,:,:,iBlock)  = ShockLeftState_V(p_) &
            + 0.5*(B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2)
    end if

    if(UseAnisoPressure) &
         ! parallel pressure
         State_VGB(Ppar_,:,:,:,iBlock) = ShockLeftState_V(Ppar_)*(1.0 &
         + 0.5*(B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2)               &
         /ShockLeftState_V(p_))

    ! Get density from the uniform temperature assumption
    State_VGB(rho_,:,:,:,iBlock) = State_VGB(p_,:,:,:,iBlock)/Tp

    ! set intial perturbation Az = exp(-x^2/GaussX^2-y^2/Gauss^2)*cos(Kx*x)*cos(Ky*y)
    do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
       x = Xyz_DGB(x_,i,j,k,iBlock)
       y = Xyz_DGB(y_,i,j,k,iBlock)
       a = Apert*B0*exp(-x**2*GaussXInv**2 - y**2*GaussYInv**2)
       ! Bx = dAz/dy
       State_VGB(Bx_,i,j,k,iBlock) = State_VGB(Bx_,i,j,k,iBlock) + &
            a*(-2*y*GaussYInv**2*cos(Kx*x)*cos(Ky*y) - Ky*cos(Kx*x)*sin(Ky*y))

       ! By = -dAz/dx
       State_VGB(By_,i,j,k,iBlock) = State_VGB(By_,i,j,k,iBlock) + &
            a*(+2*x*GaussXInv**2*cos(Kx*x)*cos(Ky*y) + Kx*sin(Kx*x)*cos(Ky*y))
    end do; end do; end do


  end subroutine user_set_ics

  !=====================================================================
  subroutine user_get_log_var(VarValue, TypeVar, Radius)

    use ModMain,     ONLY: nI, nJ, nK, nBlock, Unused_B
    use ModAdvance,  ONLY: By_, State_VGB
    use ModGeometry, ONLY: z2, z1
    use BATL_lib,    ONLY: CellFace_DB, CellSize_DB, Xyz_DGB

    real, intent(out)            :: VarValue
    character (len=*), intent(in):: TypeVar
    real, intent(in), optional :: Radius

    character (len=*), parameter :: Name='user_get_log_var'

    integer :: k1, k2, iBlock
    real:: y1, y2, dy1, dy2, HalfInvWidth, Flux
    !-------------------------------------------------------------------
    HalfInvWidth = 0.5/(z2-z1)
    VarValue=0.0
    select case(TypeVar)
    case('byflux')
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          y1 = Xyz_DGB(y_,1,0,1,iBlock)
          y2 = Xyz_DGB(y_,1,nJ+1,1,iBlock)

          if(y1*y2 > 0) CYCLE
          k1 = -y1/CellSize_DB(y_,iBlock)
          k2 = k1 + 1
          dy1 = abs(Xyz_DGB(y_,1,k1,1,iBlock))/CellSize_DB(y_,iBlock)
          dy2 = 1.0 - dy1
          Flux = CellFace_DB(2,iBlock)*HalfInvWidth* &
               ( dy2*sum(abs(State_VGB(By_,1:nI,k1,1:nK,iBlock))) &
               + dy1*sum(abs(State_VGB(By_,1:nI,k2,1:nK,iBlock))))
          if(k1==0 .or. k2==nJ+1) Flux = 0.5*Flux
          VarValue = VarValue + Flux
       end do
    case default
       call stop_mpi('Unknown user logvar='//TypeVar)
    end select

  end subroutine user_get_log_var

  !============================================================================

  subroutine user_update_states(iStage, iBlock)

    use ModAdvance,    ONLY: State_VGB
    use ModVarIndexes

    use ModSize, ONLY: nDim, nI, nJ, nK
    use ModProcMH, ONLY: iProc
    use ModMain, ONLY: n_step
    use BATL_lib, ONLY: Xyz_DGB

    use ModEnergy, ONLY: calc_energy_cell
    use ModPlotFile, ONLY: read_plot_file
    use ModInterpolate, ONLY: bilinear
    use ModUtilities,   ONLY: sleep
    use ModIoUnit,      ONLY: UnitTmp_

    integer,intent(in):: iStage, iBlock

    ! PIC coupling related local variables
    integer:: nStepLast = -1000
    integer:: iCouplePic = 0

    character(len=100):: NameFile
    integer:: i, j, k, iError
    real:: XyzNorm_D(nDim)

    ! PIC grid size
    integer, save:: nXPic, nYPic, nCellPic_D(nDim)
    real,    save:: CoordMinPic_D(nDim), CoordMaxPic_D(nDim), DxyzPic_D(nDim)

    ! PIC variables
    integer, save:: nVarPic
    real,    save, allocatable:: StatePic_VC(:,:,:), StatePic_V(:)

    integer:: Dn
    real:: WeightMhd, WeightPic

    character(len=*), parameter :: NameSub = 'user_update_states'
    !--------------------------------------------------------------------------

    call update_states_mhd(iStage, iBlock)

    if(n_step >= 1)then
       ! Overwrite the region with the PIC solution

       ! Check if we should read in a new PIC file
       if(n_step > nStepLast)then
          nStepLast  = n_step
          iCouplePic = iCouplePic + 1

          ! Construct file name
          !write(NameFile,'(a,i7.7,a)') trim(NameFilePic), iCouplePic, '.out'
          !write(NameFile,'(a,i7.7,a)') trim(NameFilePic), n_step, '.out'
          write(NameFile,'(a,a)') trim(NameFilePic), '.out'

          if(iProc == 0)write(*,*) NameSub,' trying to read ',NameFile
          ! Wait until file exists
          do
             open(UnitTmp_, FILE=NameFile, STATUS='OLD', IOSTAT=iError)
             ! If successful, wait a bit so that file is fully written
             ! If not successful, wait a bit and try open again
             if(iError /= 0)then
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
                call read_plot_file(NameFile, &
                     nVarOut = nVarPic, nOut_D = nCellPic_D, iErrorOut=iError)
                if(iError /= 0)then
                   write(*,*) NameSub,': could not read header from ', &
                        trim(NameFile), ' on processor', iProc
                   call sleep(0.5)
                   CYCLE
                end if
                EXIT
             end do
             nXPic = nCellPic_D(1)
             nYPic = nCellPic_D(2)
             allocate(StatePic_VC(nVarPic,nXPic,nYPic), StatePic_V(nVarPic))
             if(iProc == 0)write(*,*) NameSub, &
                  ' allocated StatePic_VC with nXPic, nYPic=', nXPic, nYPic

             ! Read first PIC data and coordinate limits
             do
                call read_plot_file(NameFile, VarOut_VII=StatePic_VC, &
                     CoordMinOut_D = CoordMinPic_D, &
                     CoordMaxOut_D = CoordMaxPic_D, &
                     iErrorOut=iError)
                if(iError /= 0)then
                   write(*,*) NameSub,': could not read first data from ', &
                        trim(NameFile), ' on processor', iProc
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
             do
                call read_plot_file(NameFile, VarOut_VII=StatePic_VC, &
                     iErrorOut=iError)
                if(iError /= 0)then
                   write(*,*) NameSub,': could not read data from ', &
                        trim(NameFile), ' on processor', iProc
                   call sleep(0.5)
                   CYCLE
                end if
                EXIT
             end do

          end if
       end if


       ! As we resuse the same filename, we will need to delete data we have read  
       call barrier_mpi
       if(iProc ==0 ) then
          open(UnitTmp_, FILE=NameFile, STATUS='OLD', IOSTAT=iError)
          close(UnitTmp_,STATUS='DELETE')
       end if
       
       ! Overwrite cells inside the PIC domain
       if(  all(Xyz_DGB(1:nDim, 1, 1, 1,iBlock) <= CoordMaxPic_D) .and. &
            all(Xyz_DGB(1:nDim,nI,nJ,nK,iBlock) >= CoordMinPic_D)) then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI

             ! Normalized PIC grid coordinates (1...nCellPic_D)
             XyzNorm_D = 1 + (Xyz_DGB(1:nDim,i,j,k,iBlock) - CoordMinPic_D) &
                  /DxyzPic_D

             ! Distance from edge
             Dn = minval( min(nint(XyzNorm_D - 1), nint(nCellPic_D - XyzNorm_D)) )

             ! Nothing to do within PIC ghost region
             if(Dn < nGhostPic) CYCLE

             ! Distance from ghost layers
             Dn = Dn - nGhostPic + 1

             if(Dn <= nOverlapPic)then
                ! For nOverlapPic=1, Dn = 1, so use 0.5 as weight
                ! For nOverlapPic=2, Dn = 1, 2, so use 1/3 and 2/3 weights for PIC
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

             ! Interpolate MHD and PIC states. Skip hyperbolic scalar if present (Hyp=Bz_+1)
             State_VGB(Rho_:Bz_,i,j,k,iBlock) = WeightMhd*State_VGB(Rho_:Bz_,i,j,k,iBlock) &
                  + WeightPic*StatePic_V(Rho_:Bz_)

             State_VGB(p_,i,j,k,iBlock) = WeightMhd*State_VGB(p_,i,j,k,iBlock) &
                  + WeightPic*StatePic_V(nVarPic)

             ! Set hyperbolic scalar to zero if present
             if(Hyp_>1) State_VGB(Hyp_,i,j,k,iBlock) = 0.0

          end do; end do; end do

          call calc_energy_cell(iBlock)
       end if
    end if

  end subroutine user_update_states

  !============================================================================


end module ModUser
