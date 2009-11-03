!^CFG COPYRIGHT UM
!==============================================================================
module ModUser
  use ModMain, ONLY: nBLK
  use ModSize, ONLY: nI,nJ,nK
  use ModReadParam, ONLY: lStringLine
  use ModVarIndexes
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_init_session,               &
       IMPLEMENTED3 => user_set_ics,                    &
       IMPLEMENTED4 => user_initial_perturbation,       &
       IMPLEMENTED5 => user_face_bcs,                   &
       IMPLEMENTED6 => user_get_log_var,                &
       IMPLEMENTED8 => user_update_states,              &
       IMPLEMENTED10=> user_set_boundary_cells,         &
  
  include 'user_module.h' !list of public methods
  
  real, parameter               :: VersionUserModule = 1.0
  character (len=*), parameter  :: &
       NameUserModule = 'Multigroup Frequency Advection'
  character(len=lStringLine)    :: NameModel

  !Global variables - frequency grid
  logical                       :: IsInitWave = .false.

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
    UseUserInitSession = .true.

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
  subroutine user_init_session
    use ModIO,          ONLY: write_prefix, iUnitOut,NamePlotDir
    use ModMain,        ONLY: UseUserB0
    use ModPhysics,     ONLY: BodyNDim_I,BodyTDim_I,g
    use ModProcMH,      ONLY: iProc
    use ModReadParam,   ONLY: i_line_command
    implicit none
    !--------------------------------------------------------------------------
    if(iProc == 0)then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'user_init_session:'
       call write_prefix; write(iUnitOut,*) ''
    end if

    

    if(iProc == 0)then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'user_init_session finished'
       call write_prefix; write(iUnitOut,*) ''
    end if

  end subroutine user_init_session
  !============================================================================
  subroutine user_face_bcs(VarsGhostFace_V)
    use ModSize,       ONLY: East_,West_,South_,North_,Bot_,Top_,nDim
    use ModMain,       ONLY: time_accurate,x_,y_,z_, n_step, Iteration_Number
    use ModVarIndexes 
    use ModAdvance,    ONLY: State_VGB
    use ModPhysics,    ONLY: No2Si_V,Si2No_V,UnitU_,UnitRho_,UnitP_,UnitX_
    use ModNumConst,   ONLY: cTolerance,cTiny
    use ModFaceBc,     ONLY: FaceCoords_D, VarsTrueFace_V, TimeBc, &
                             iFace, jFace, kFace, iSide, iBlockBc
    use ModWaves,      ONLY: AlfvenSpeedPlusFirst_,  &
                             AlfvenSpeedPlusLast_,   &
                             AlfvenSpeedMinusFirst_, &
                             AlfvenSpeedMinusLast_ , &
                             DeltaLogFrequency
    implicit none

    real, intent(out):: VarsGhostFace_V(nVar)

    integer:: iCell,jCell,kCell, iWave

    real:: DensCell,PresCell,GammaCell,TBase,B1dotR  
    real, dimension(3):: RFace_D,B1_D,U_D,B1t_D,B1n_D
    real :: BRCell, vAlfvenSi, wEnergyDensSi, wEnergyDens
  
    !--------------------------------------------------------------------------

    RFace_D  = FaceCoords_D/sqrt(sum(FaceCoords_D**2))

    U_D (x_:z_)  = VarsTrueFace_V(Ux_:Uz_)
    B1_D(x_:z_)  = VarsTrueFace_V(Bx_:Bz_)
    B1dotR       = dot_product(RFace_D,B1_D)
    B1n_D(x_:z_) = B1dotR*RFace_D(x_:z_)
    B1t_D        = B1_D-B1n_D

    !\
    ! Update BCs for velocity and induction field::
    !/
    VarsGhostFace_V(Ux_:Uz_) = -U_D(x_:z_)
    VarsGhostFace_V(Bx_:Bz_) = B1t_D(x_:z_)!-B1n_D(x_:z_)

    !\
    ! Update BCs for the mass density, EnergyRL, 
    ! and pressure::
    !/
    iCell = iFace; jCell = jFace; kCell = kFace
    select case(iSide)
    case(East_)
       iCell  = iFace
    case(West_)
       iCell  = iFace-1
    case(South_)
       jCell  = jFace
    case(North_)
       jCell  = jFace-1
    case(Bot_)
       kCell  = kFace
    case(Top_)
       kCell  = kFace-1
    case default
       write(*,*)'ERROR: iSide = ',iSide
       call stop_mpi('incorrect iSide value in user_face_bcs')
    end select

    !\
    ! Update BCs for wave spectrum
    !/
    if(.not.IsInitWave) then
       BRCell = sum(RFace_D * &
            (State_VGB(Bx_:Bz_,iCell,jCell,kCell,iBlockBc) + &
            B0_DGB(:,iCell,jCell,kCell,iBlockBc)) )
       vAlfvenSi = (BRCell/sqrt(VarsTrueFace_V(Rho_))) * No2Si_V(UnitU_)
      
       call get_total_wave_energy_dens(&
            FaceCoords_D(x_),&
            FaceCoords_D(y_),&
            FaceCoords_D(z_),&
            vAlfvenSi, wEnergyDensSi)
      
       wEnergyDens = wEnergyDensSI * Si2No_V(UnitP_)
       do iWave = AlfvenSpeedPlusFirst_,AlfvenSpeedPlusLast_
          VarsGhostFace_V(iWave) = (2.0/3.0) * DeltaLogFrequency * wEnergyDens* &
               exp((LogFreq_I(iWave-WaveFirst_+1)-LogFreq_I(1))*(-2.0/3.0))
       end do
    end if
  end subroutine user_face_bcs
  !===========================================================================
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
     use ModAdvance, ONLY: State_VGB, Source_VC
     use ModWaves,   ONLY: UseWavePressure
    implicit none

    integer,intent(in)           :: iStage,iBlock
    real                         :: DensCell,PresCell,GammaCell,Beta,WavePres
    character(len=*),parameter   :: NameSub='user_update_states'
    !--------------------------------------------
    
    !Disable advection of waves with medium
    Source_VC(:,:,;,:) = 0.0
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
    use ModPhysics, ONLY: No2Si_V, UnitX_, UnitP_
    use ModWaves

    implicit none
    
    real, allocatable,dimension(:,:) :: Cut_II ! Array to store log variables
    integer                          :: nCell,nRow,iRow 
    real                             :: dx, dz, x,y,z, IwPlusSi,IwMinusSi
    integer                          :: iFreq,i,j,k,iBLK
    integer                          :: iUnit,iError,aError
    character(len=40)                :: FileNameTec,FileNameHead,FileNameGrid 
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
    nRow=nCell*nWave/2
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
          if(unusedBLK(iBLK)) CYCLE
          do k=1,nK ; do j=1,nJ ; do i=1,nI
             x=x_BLK(i,j,k,iBLK)
             y=y_BLK(i,j,k,iBLK)
             z=z_BLK(i,j,k,iBLK)
             dx=dx_BLK(iBLK)
             dz=dz_BLK(iBLK)
             if((z< dz) .and. (z >=0.0) .and. (x<dx) .and. (x>=0.0)) then
                do iFreq=1,nWave/2
                   IwPlusSi  = No2Si_V(UnitP_)*State_VGB(AlfvenSpeedPlusFirst_+iFreq-1,i,j,k,iBLK)
                   IwMinusSi = No2Si_V(UnitP_)*State_VGB(AlfvenSpeedMinusFirst_+iFreq-1,i,j,k,iBLK)
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
    ! write header file (iProc==0)
    !/
    write(NameStage,'(i5.5)') iteration_number
    write(NameProc,'(a,i4.4)') "_pe",iProc
   
    if (iProc==0) then
       FileNameHead='SC/IO2/Spectrum_n_'//trim(NameStage)//'.H'
       write(*,*) 'SC:  writing file ', FileNameHead
       iUnit=io_unit_new()
       open(unit=iUnit, file=FileNameHead, form='formatted', access='sequential',&
         status='replace',iostat=iError)
       write(iUnit, '(a)') 'Title: BATSRUS SC Spectrogram'
       write(iUnit, '(i4.4)') nProc
       write(iUnit, '(a)') 'Variables = "Y[R]", "w","I+[Jm-3]","I-[Jm-3]" '
       close(iUnit)
    end if
    !\
    ! Write grid size file
    !/
    FileNameGrid='SC/IO2/Spectrum_n_'//trim(NameStage)//trim(NameProc)//'.G'
    iUnit=io_unit_new()
    open(unit=iUnit, file=FileNameGrid,form='formatted',access='sequential',&
         status='replace',iostat=iError)
    write(iUnit,'(i6.6)') nWave/2
    write(iUnit,'(i6.6)') nCell

    !\
    ! Write data file
    !/
    FileNameTec='SC/IO2/Spectrum_n_'//trim(NameStage)//trim(NameProc)//'.tec'
    iUnit=io_unit_new()
    open(unit=iUnit, file=FileNameTec, form='formatted',access='sequential',&
         status='replace',iostat=iError)
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

    ! This subroutine calculates the cut-off frequency for Alfven waves, 
    ! which is the ion cyclotron frequency (in radians) 
    ! \omega_{c.o.}=zeB/m. For ions z=1. 
     
    use ModVarIndexes
    use ModAdvance,           ONLY: State_VGB
    use ModConst,             ONLY: cElectronCharge, cProtonMass,cTiny
    use ModSize,              ONLY: nI,nJ,nK
    use ModPhysics,           ONLY: No2Si_V,UnitB_

    implicit none
     
    real, intent(out)           :: LogFreqCutOff
    integer, intent(in)         :: i,j,k,iBLK
    real                        :: BxNo, ByNo, BzNo,BtotSi 
    character(len=*),parameter  :: NameSub = 'calc_cutoff_freq'
    ! -----------------------------------------------------------------

    BxNo = State_VGB(Bx_,i,j,k,iBLK)
    ByNo = State_VGB(By_,i,j,k,iBLK)
    BzNo = State_VGB(Bz_,i,j,k,iBLK)

    BtotSi = No2Si_V(UnitB_)*sqrt(BxNo**2 + ByNo**2 + BzNo**2)
    BtotSi = max(BtotSi,cTiny**2)
    LogFreqCutOff = log((cElectronCharge*BtotSi)/cProtonMass)

  end subroutine calc_cutoff_freq
  !=======================================================================
  subroutine set_freq_grid

    use ModVarIndexes
    use ModNumConst, ONLY: cPi
    use ModWaves,    ONLY: FreqMinSI,FreqMaxSI,&
         DeltaLogFrequency
    implicit none
    
    integer                    :: iFreq,nWaveHalf
    real                       :: LogFreqMin, LogFreqMax
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
    use ModGeometry,    ONLY: R_BLK
    use ModSize,        ONLY: nI,nJ,nK
    use ModNumConst,    ONLY: cTiny,cPi
    use ModPhysics,     ONLY: No2Si_V,Si2No_V,UnitB_,UnitRho_,UnitX_,UnitP_,UnitU_
    use ModWaves,       ONLY: DeltaLogFrequency
    
    implicit none

    ! Spatial grid variables
    integer                      :: i,j,k,iBLK, iWave
    character(len=*),parameter    :: NameSub= 'init_wave_spectrum'
    ! ------------------------------------------------------------------
    IsInitWave=.true.

    call set_freq_grid
 
    State_VGB(WaveFirst_:WaveLast_,:,:,:,:) = 0.0

    do iBLK=1,nBLK
       if (unusedBLK(iBLK)) CYCLE
 
       do k=1,nK ; do j = i,nJ ; do i=1,nI

          call calc_cutoff_freq(i,j,k,iBLK, LogFreqCutOff)
          do iWave=1,nWave
             if (LogFreq_I(iWave) .ge. LogFreqCutOff) then
                WaveEnergy_I(iWave) = 0.0
             else
                WaveEnergy_I(iWave)=exp(LogFreq_I(iWave) * FreqPower)
             end if
         
             State_VGB(WaveFirst_+iWave-1,i,j,k,iBLK) =&
                  WaveEnergy_I(iWave)*DeltaLogFrequency*Si2No_V(UnitP_)
             if(State_VGB(WaveFirst_+iWave-1,i,j,k,iBLK) < 0.0) then
                write(*,*) 'Negative energy at init spectrum'
                write(*,*) 'dLogFreq= ', DeltaLogFrequency
                write(*,*) 'WaveEnergy_I', WaveEnergy_I(iWave)
             end if
          end do
       end do; end do ; end do
    end do

  end subroutine init_wave_spectrum
  !========================================================================
  subroutine user_get_log_var(VarValue,TypeVar,Radius)
    
    use ModIO,         ONLY: write_myname
   
    character (LEN=10), intent(in):: TypeVar 
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

