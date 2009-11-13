!^CFG COPYRIGHT UM
!==============================================================================
module ModUser
  use ModMain, ONLY: nBLK
  use ModSize, ONLY: nI,nJ,nK
  use ModReadParam, ONLY: lStringLine
  use ModVarIndexes
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED3 => user_set_ics,                    &
       IMPLEMENTED4 => user_initial_perturbation,       &
       IMPLEMENTED6 => user_get_log_var,                &
       IMPLEMENTED8 => user_update_states,              &
       IMPLEMENTED10=> user_set_boundary_cells       
  
  include 'user_module.h' !list of public methods
  
  real, parameter               :: VersionUserModule = 1.0
  character (len=*), parameter  :: NameUserModule = 'Multigroup Frequency Advection'
  character(len=lStringLine)    :: NameModel

  logical                       :: IsInitWave = .false.
  real                          :: LogFreqCutOff, SpectralIndex
  integer                       :: LowestFreqNum,SpectrumWidth
  integer                       :: nWaveHalf ! number of waves in each direction
  real,dimension(nWave)         :: LogFreq_I ! frequency grid

contains
  !============================================================================
  subroutine user_read_inputs
    use ModMain
    use ModProcMH,      ONLY: iProc
    use ModReadParam,   ONLY: read_line, read_command, read_var
    use ModIO,          ONLY: write_prefix, write_myname, iUnitOut,NamePlotDir
    use ModWaves,       ONLY: read_alfven_speed,read_wave_pressure,read_frequency
    implicit none

    character (len=100) :: NameCommand
    !--------------------------------------------------------------------------
    UseUserInitSession = .false.

    if(iProc == 0 .and. lVerbose > 0)then
       call write_prefix;
       write(iUnitOut,*)'User read_input for log advection starts'
    endif

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)

       case("#ALFVENSPEED")
          call read_alfven_speed

       case("#WAVEPRESSURE")
          call read_wave_pressure

       case("#FREQUENCY")
          call read_frequency

       case("#SPECTRUM")
          call read_var('LowestFreqNum',LowestFreqNum)
          call read_var('SpectrumWidth',SpectrumWidth)
          call read_var('SpectralIndex',SpectralIndex)

       case('#USERINPUTEND')
          if(iProc == 0 .and. lVerbose > 0)then
             call write_prefix;
             write(iUnitOut,*)'User read_input TURBULENCE CORONA ends'
          endif
          EXIT
       case default
          if(iProc == 0) then
             call write_myname; write(*,*) &
                  'ERROR: Invalid user defined #COMMAND in user_read_inputs. '
             write(*,*) '--Check user_read_inputs for errors'
             write(*,*) '--Check to make sure a #USERINPUTEND command was used'
             write(*,*) '  *Unrecognized command was: '//NameCommand
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if
       end select
    end do

  end subroutine user_read_inputs
  !============================================================================
  subroutine user_initial_perturbation
    use ModMain, ONLY: nBLK,unusedBLK,x_,y_,z_,n_step
    use ModGeometry
    use ModProcMH

    implicit none

    logical :: oktest,oktest_me

    !--------------------------------------------------------------------------
 
        call set_oktest('user_initial_perturbation',oktest,oktest_me)

        call init_wave_spectrum

        if (iProc==0) then
           
           write(*,*) 'SC: Finished initializing wave spectrum'
        end if
     
   end subroutine user_initial_perturbation
  !============================================================================
  subroutine user_set_ics
    use ModMain,      ONLY: globalBLK,nI,nJ,nK
    use ModVarIndexes
    use ModAdvance,   ONLY: State_VGB 
    use ModPhysics,   ONLY: inv_gm1,BodyTDim_I
    use ModGeometry
    implicit none

    integer :: i,j,k,iBLK
    logical :: oktest,oktest_me
    !--------------------------------------------------------------------------
    call set_oktest('user_set_ics',oktest,oktest_me)
    iBLK = globalBLK

    do k=1,nK; do j=1,nJ; do i=1,nI
       State_VGB(Bx_:Bz_,i,j,k,iBLK) = 0.0
       State_VGB(Uy_:Uz_,i,j,k,iBLK) = 0.0
       State_VGB(Ux_,i,j,k,iBLK) = real(i)
    end do; end do; end do

  end subroutine user_set_ics
  !============================================================================
   subroutine user_update_states(iStage,iBlock)
     use ModVarIndexes
     use ModSize
     use ModAdvance
     use ModWaves,   ONLY: UseWavePressure
    implicit none

    integer,intent(in)           :: iStage,iBlock
    real                         :: DensCell,PresCell,GammaCell,Beta,WavePres
    character(len=*),parameter   :: NameSub='user_update_states'
    !--------------------------------------------
    
    !Disable advection of waves with medium
    Source_VC(:,:,:,:) = 0.0
    Flux_VX(:,:,:,:) = 0.0
    Flux_VY(:,:,:,:) = 0.0
    Flux_VZ(:,:,:,:) = 0.0

    ! Advect solution in frequency dimension only
    if(any(State_VGB(WaveFirst_:WaveLast_,:,:,:,iBlock)<0.0)) then
       write(*,*) NameSub,' : negative wave energy before MHD'
    end if
    call update_states_MHD(iStage,iBlock)
    if(any(State_VGB(WaveFirst_:WaveLast_,:,:,:,iBlock)<0.0)) then
       write(*,*) NameSub, ': negative wave energy after MHD'
    end if

  end subroutine user_update_states
  !=======================================================================
   subroutine write_spectrogram
   
    use ModProcMH
    use ModMain,   ONLY: iteration_number, nBLK,unusedBLK,nBlockALL
    use ModSize,   ONLY: nI,nJ,nK
    use ModGeometry, ONLY: x_BLK,y_BLK, z_BLK, dz_BLK,dx_BLK
    use ModIoUnit, ONLY: io_unit_new
    use ModVarIndexes
    use ModAdvance, ONLY: State_VGB
    use ModWaves

    implicit none
    
    real, allocatable,dimension(:,:) :: Cut_II ! Array to store log variables
    integer                          :: nCell,nRow,iRow 
    real                             :: dx, dz, x,y,z, IwPlusSi,IwMinusSi
    integer                          :: iFreq,i,j,k,iBLK
    integer                          :: iUnit,iError,aError
    character(len=40)                :: FileNameTec 
    character(len=11)                :: NameStage
    character(len=7)                 :: NameProc
    character(len=*),parameter       :: NameSub='write_spectrogram'
    !-------------------------------------------------------------------
    !\
    ! count cells in cut x=0, z=0
    !/
    nCell=0
    do iBLK=1,nBLK
       if(unusedBLK(iBLK)) CYCLE
       do k=1,nK ; do j=1,nJ ; do i=1,nI
          x=x_BLK(i,j,k,iBLK)
          dx=dx_BLK(iBLK)
          z=z_BLK(i,j,k,iBLK)
          dz=dz_BLK(iBLK)
          if((z< dz) .and. (z >=0.0) .and. (x<dx) .and. (x>=0.0)) then
             nCell=nCell+1
          end if
       end do; end do ; end do
    end do
    nRow=nCell*nWaveHalf
    !\
    ! Allocate plot arrays
    !/
    ALLOCATE(Cut_II(nRow,4),STAT=aError)
    !\
    ! Fill plot array
    !/
    if (aError .ne. 0) then
       call stop_mpi('Allocation failed for spectrogram array')
    else
       iRow=1
       do iBLK=1,nBLK
          !if(unusedBLK(iBLK)) CYCLE
          do k=1,nK ; do j=1,nJ ; do i=1,nI
             x=x_BLK(i,j,k,iBLK)
             y=y_BLK(i,j,k,iBLK)
             z=z_BLK(i,j,k,iBLK)
             dx=dx_BLK(iBLK)
             dz=dz_BLK(iBLK)
             if((z< dz) .and. (z >=0.0) .and. (x<dx) .and. (x>=0.0)) then
                do iFreq=1,nWaveHalf
                   IwPlusSi  = State_VGB(AlfvenSpeedPlusFirst_+iFreq-1,i,j,k,iBLK)
                   IwMinusSi = State_VGB(AlfvenSpeedMinusFirst_+iFreq-1,i,j,k,iBLK)
                   Cut_II(iRow,1) = y
                   Cut_II(iRow,2) = LogFreq_I(iFreq)
                   Cut_II(iRow,3) = IwPlusSi
                   Cut_II(iRow,4) = IwMinusSi
                   iRow=iRow+1
                end do
             end if
          end do; end do ; end do
       end do
    end if
    !\
    ! write data file
    !/
    write(NameStage,'(i5.5)') iteration_number
    write(NameProc,'(a,i4.4)') "_pe",iProc
    FileNameTec='SC/IO2/Spectrum_n_'//trim(NameStage)//trim(NameProc)//'.tec'
    if (iProc==0) then
       write(*,*) 'SC:  writing file ', FileNameTec
    end if
    iUnit=io_unit_new()
    open(unit=iUnit, file=FileNameTec,form='formatted',access='sequential',&
         status='replace',iostat=iError)
    write(iUnit, '(a)') 'Title: BATSRUS SC Spectrogram'
    write(iUnit, '(a)') 'Variables = "Y[R]", "log(w)","I+[Jm-3]","I-[Jm-3]" '
    write(iUnit,'(a,i3.3,a,i5.5,a)') 'Zone I= ',nWaveHalf,' J= ',nCell,' F=point'
    do iRow=1,nRow
       write(iUnit, fmt="(30(e14.6))") &
            Cut_II(iRow,:)
    end do
    close(iUnit)

    if(allocated(Cut_II)) deallocate(Cut_II,STAT=aError)
    if(aError .ne. 0) then
       write(*,*) NameSub, 'Deallocation of spectrogram array failed'
       call stop_mpi(NameSub)
    end if
  end subroutine write_spectrogram
  !=====================================================================
  subroutine calc_cutoff_freq(i,j,k,iBLK , LogFreqCutOff)

    ! Arbitrary value, set to be LogFreq_I(nWaveHalf-1) f
       
    use ModVarIndexes

    implicit none
     
    real, intent(out)           :: LogFreqCutOff
    integer, intent(in)         :: i,j,k,iBLK

    character(len=*),parameter  :: NameSub = 'calc_cutoff_freq'
    ! -----------------------------------------------------------------

    LogFreqCutOff = LogFreq_I(nWaveHalf-1)

  end subroutine calc_cutoff_freq
  !=======================================================================
  subroutine set_freq_grid

    use ModVarIndexes
    use ModNumConst, ONLY: cPi
    use ModWaves,    ONLY: FreqMinSI, DeltaLogFrequency
    implicit none
    
    integer                    :: iFreq
    real                       :: LogFreqMin
    character(len=*),parameter :: NameSub='set_freq_grid'
    !-----------------------------------------------------------------
    nWaveHalf = max(nWave/2,1)
    ! Minimal frequency in frequency grid
    LogFreqMin = log(2*cPi*FreqMinSI) 

    ! calculate frequencies
    ! Plus waves (+Va)
    do iFreq = 1,nWaveHalf
       LogFreq_I(iFreq)=LogFreqMin+(iFreq-1)*DeltaLogFrequency
    end do
    ! Minus waves (-Va)
    do iFreq = 1,nWaveHalf
       LogFreq_I(nWaveHalf+iFreq)=LogFreqMin+(iFreq-1)*DeltaLogFrequency
    end do
   end subroutine set_freq_grid
  !===================================================================
  subroutine init_wave_spectrum
  
    use ModMain,        ONLY: unusedBLK,nBLK
    use ModVarIndexes
    use ModAdvance,     ONLY: State_VGB
    use ModSize,        ONLY: nI,nJ,nK
    use ModWaves,       ONLY: DeltaLogFrequency, &
         AlfvenSpeedPlusFirst_,&
         AlfvenSpeedMinusFirst_
    
    implicit none
    integer                    :: i,j,k,iBLK, iWave
    character(len=*),parameter :: NameSub= 'init_wave_spectrum'
    ! ------------------------------------------------------------------
    IsInitWave=.true.
    call set_freq_grid
    State_VGB(WaveFirst_:WaveLast_,:,:,:,:) = 0.0
    SpectralIndex = SpectralIndex + 1 ! State_VGB(iWave) represents I*w
    do iBLK=1,nBLK
       do k=1,nK ; do j = 1,nJ ; do i=1,nI
          
          do iWave=1,nWaveHalf
             if( (iWave .le. LowestFreqNum) .or. &
                  (iWave .gt. LowestFreqNum+SpectrumWidth)) then
                State_VGB(AlfvenSpeedPlusFirst_+iWave-1,i,j,k,iBLK) = 0.0
                State_VGB(AlfvenSpeedMinusFirst_+iWave-1,i,j,k,iBLK)= 0.0
             else
                State_VGB(AlfvenSpeedPlusFirst_+iWave-1,i,j,k,iBLK) = &
                     exp(LogFreq_I(iWave)*SpectralIndex)
                State_VGB(AlfvenSpeedMinusFirst_+iWave-1,i,j,k,iBLK)= &
                     exp(LogFreq_I(iWave)*SpectralIndex)
             end if
          end do
       end do; end do ; end do
    end do

  end subroutine init_wave_spectrum
  !========================================================================
  subroutine user_get_log_var(VarValue,TypeVar,Radius)
    
    use ModIO,         ONLY: write_myname
   
    character (LEN=10), intent(in):: TypeVar
    real,intent(out)              :: VarValue
    real,intent(in),optional      :: Radius
    !--------------------------------------------------------------------------
    !\
    ! Define log variable to be saved::
    !/
    select case(TypeVar)
    case('spec')
       call write_spectrogram
    case default
       VarValue = -7777.
       call write_myname;
       write(*,*) 'Warning in set_user_logvar: unknown logvarname = ',TypeVar
    end select
  end subroutine user_get_log_var
 !===========================================================================
  subroutine user_set_boundary_cells(iBLK)
    use ModGeometry,      ONLY: ExtraBc_, IsBoundaryCell_GI, r_Blk
    use ModBoundaryCells, ONLY: SaveBoundaryCells
    use ModPhysics,       ONLY: rBody
    implicit none

    integer, intent(in):: iBLK

    character (len=*), parameter :: Name='user_set_boundary_cells'
    !--------------------------------------------------------------------------
    IsBoundaryCell_GI(:,:,:,ExtraBc_) = r_Blk(:,:,:,iBLK) < rBody

    if(SaveBoundaryCells) return
    call stop_mpi('Set SaveBoundaryCells=.true. in PARAM.in file')

  end subroutine user_set_boundary_cells

end module ModUser

