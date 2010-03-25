!^CFG COPYRIGHT UM
!==============================================================================
module ModUser
  use ModMain, ONLY: nBLK
  use ModSize, ONLY: nI,nJ,nK
  use ModVarIndexes
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_init_session,               &
       IMPLEMENTED3 => user_set_ics,                    &
       IMPLEMENTED4 => user_face_bcs,                   &
       IMPLEMENTED5 => user_get_log_var,                &
       IMPLEMENTED6 => user_get_b0,                     &
       IMPLEMENTED7 => user_update_states,              &
       IMPLEMENTED8 => user_specify_refinement,         &
       IMPLEMENTED9=> user_set_boundary_cells,         &
       IMPLEMENTED10=> user_set_plot_var
  
  include 'user_module.h' !list of public methods
  
   real, parameter               :: VersionUserModule = 1.0
  character (len=*), parameter  :: &
       NameUserModule = 'Alfven Waves Driven Solar Wind'

  logical                       :: IsInitWave = .false.
  logical                       :: DoDampCutOff = .false., DoDampSurface = .false.
  real,dimension(nWave)         :: LogFreq_I ! frequency grid
  real                          :: LogFreqCutOff,WaveInnerBcFactor
  real                          :: SpectralIndex, SpectralCoeff, dLogFreq
  integer                       :: WaveFirstBc_, WaveLastBc_, nWaveHalf
  real                          :: xTrace = 0.0, zTrace = 0.0
  real                          :: xTestSpec,yTestSpec,zTestSpec
  character(len=10)             :: TypeWaveInnerBc
contains
  !============================================================================
  subroutine user_read_inputs
    use ModMain
    use ModProcMH,      ONLY: iProc
    use ModReadParam,   ONLY: read_line, read_command, read_var
    use ModIO,          ONLY: write_prefix, write_myname, iUnitOut,NamePlotDir
 
    character (len=100) :: NameCommand
    !--------------------------------------------------------------------------
    UseUserInitSession = .true.

    if(iProc == 0 .and. lVerbose > 0)then
       call write_prefix;
       write(iUnitOut,*)'User read_input TURBULENCE CORONA starts'
    endif

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)
          
       case("#INITSPECTRUM")
          call read_var('IsInitWave',IsInitWave)
       
       case("#WAVEINNERBC")
          call read_var('TypeWaveInnerBc', TypeWaveInnerBc)
          call read_var('WaveInnerBcFactor',WaveInnerBcFactor)
          call read_var('WaveFirstBc_',WaveFirstBc_)
          call read_var('WaveLastBc_',WaveLastBc_)
          call read_var('SpectralIndex',SpectralIndex)
         
       case("#DAMPWAVES")
          call read_var('DoDampCutoff',DoDampCutoff)
          call read_var('DoDampSurface',DoDampSurface)

       case("#SPECTROGRAM")
          call read_var('xTrace',xTrace)
          call read_var('zTrace',zTrace)
       case("#TESTSPEC")
          call read_var('xTestSpec',xTestSpec)
          call read_var('yTestSpec',yTestSpec)
          call read_var('zTestSpec',zTestSpec)
      
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

    use ModIO,          ONLY: write_prefix, iUnitOut
    use ModProcMH,      ONLY: iProc
    use ModReadParam,   ONLY: i_line_command
 
    !--------------------------------------------------------------------------
    if(iProc == 0)then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'user_init_session:'
       call write_prefix; write(iUnitOut,*) ''
    end if

    if(i_line_command("#INITSPECTRUM") > 0 .and. IsInitWave) then
       call init_wave_spectrum
    end if

      if(iProc == 0)then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'user_init_session finished'
       call write_prefix; write(iUnitOut,*) ''
    end if

  end subroutine user_init_session
  !====================================================================
  subroutine set_freq_grid

    use ModProcMH,      ONLY: iProc
    use ModNumConst,    ONLY: cPi
    use ModWaves,       ONLY: FreqMinSI, FreqMaxSI
    
    integer                    :: iFreq
    real                       :: LogFreqMinRadian
    character(len=*),parameter :: NameSub='set_freq_grid'
    !-----------------------------------------------------------------
    nWaveHalf = max(nWave/2, 1)
    LogFreqMinRadian = log(2*cPi*FreqMinSI) 
    dLogFreq = log(FreqMaxSI/FreqMinSI)/(nWaveHalf-1)
    if(iProc == 0) write(*,*) &
         'Setting Frequency grid, dLogfreq = ',dLogFreq

    ! calculate frequencies
    ! Plus waves - parallel to B
    do iFreq = 1,nWaveHalf
       LogFreq_I(iFreq)=LogFreqMinRadian+(iFreq-1)*dLogFreq
    end do
    ! Minus waves- antiparallel to B
    do iFreq = 1,nWaveHalf
       LogFreq_I(nWaveHalf+iFreq)=LogFreqMinRadian +(iFreq-1)*dLogFreq
    end do
   end subroutine set_freq_grid
  !===================================================================
  subroutine init_wave_spectrum
   
    use ModProcMH,      ONLY: iProc

    character(len=*),parameter    :: NameSub= 'init_wave_spectrum'
    ! ------------------------------------------------------------------
    call set_freq_grid
    SpectralIndex = SpectralIndex + 1 ! State_VGB(iWave) represents I*w
    SpectralCoeff = -(1.0/SpectralIndex)
    if (iProc == 0) &
         write(*,*) 'Spectral coefficient is ',SpectralCoeff

  end subroutine init_wave_spectrum
  !============================================================================
  subroutine user_face_bcs(VarsGhostFace_V)
   
    use ModSize,             ONLY: East_,West_,South_,North_,Bot_,Top_,nDim
    use ModMain,             ONLY: x_,y_,z_, UseRotatingFrame
    !use ModExpansionFactors, ONLY: ExpansionFactorInv
    use ModAdvance,          ONLY: State_VGB, B0_DGB
    use ModPhysics,          ONLY: OmegaBody,No2Si_V,Si2No_V,UnitU_,UnitP_
    use ModNumConst,         ONLY: cTolerance
    use ModFaceBc,           ONLY: FaceCoords_D, VarsTrueFace_V, &
                                   iFace, jFace, kFace, iSide, iBlockBc, B0Face_D
    use ModWaves,            ONLY: AlfvenWavePlusFirst_,AlfvenWavePlusLast_,&
                                   AlfvenWaveMinusFirst_,AlfvenWaveMinusLast_ , &
                                   UseWavePressureLtd      

    real, intent(out)           :: VarsGhostFace_V(nVar)
    real                        :: WavesGhostFace_I(nWaveHalf)
    integer                     :: iCell,jCell,kCell
    real                        :: DensCell,PresCell,TBase,FullBr, FullB  
    real, dimension(3)          :: RFace_D,B1r_D,B1t_D,FullB_D
    real                        :: vAlfvenSi, wEnergyDensSi, wEnergyDensBc
    character(len=*),parameter  :: NameSub = "user_face_bc"
    !--------------------------------------------------------------------------
    ! vector normal to the face
    RFace_D  = FaceCoords_D/sqrt(sum(FaceCoords_D**2))

    !\
    ! Magnetic field reltaed quantities
    !/
    ! B1 normal to the face
    B1r_D = sum(RFace_D*VarsTrueFace_V(Bx_:Bz_))*RFace_D(x_:z_)
   
    ! B1 tangential
    B1t_D = VarsTrueFace_V(Bx_:Bz_) - B1r_D

    ! Total mangnetic field - add B0
    FullB_D = B0Face_D + B1t_D
    FullBr = sum(RFace_D*FullB_D)
    FullB = sum(FullB_D**2)
   
    !\
    ! Update BCs for velocity and induction field: reflective BC
    !/
    VarsGhostFace_V(Ux_:Uz_) = -VarsTrueFace_V(Ux_:Uz_)
    VarsGhostFace_V(Bx_:Bz_) = B1t_D
  
    !\
    ! Update BCs for the mass density and pressure::
    !/
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

    call get_plasma_parameters_cell(iCell,jCell,kCell,iBlockBc,&
         DensCell,PresCell)
    VarsGhostFace_V(Rho_) = &
         max(-VarsTrueFace_V(Rho_) + 2.0*(DensCell), &
         VarsTrueFace_V(Rho_))
    TBase = PresCell/DensCell
    VarsGhostFace_V(P_) = &
         max(VarsGhostFace_V(Rho_)*TBase,VarsTrueFace_V(P_))
    
    !\
    ! Update BCs for wave spectrum
    !/
    if(IsInitWave) then
       !if(ExpansionFactorInv < cTolerance) then
       !   ! no wave energy in closed field lines
       !   wEnergyDens = 1.0e-30
       !else
          select case(TypeWaveInnerBc)
          case('WSA')
             vAlfvenSi = (FullB/sqrt(VarsGhostFace_V(Rho_))) * No2Si_V(UnitU_)
             call get_total_wave_energy_dens(&
                  FaceCoords_D(x_),&
                  FaceCoords_D(y_),&
                  FaceCoords_D(z_),&
                  vAlfvenSi, wEnergyDensSi)
             wEnergyDensBc = wEnergyDensSi * Si2No_V(UnitP_)*WaveInnerBcFactor
  
          case('Turb')
             wEnergyDensBc = FullB*WaveInnerBcFactor
          end select
          if(wEnergyDensBc < 0.0)then
             write(*,*) 'Negative total wave energy density at inner BC'
             call stop_MPI('Error in face_bcs')
          end if
       !end if
       !\
       ! Deconstruct total energy into frequency groups
       !/
       call get_wave_group_bc(wEnergyDensBc,WavesGhostFace_I)
        
       if (FullBr .le. 0.0) then
          
          ! full absorption of "Plus"  waves - float bc
          VarsGhostFace_V(AlfvenWavePlusFirst_:AlfvenWavePlusLast_) =&
               VarsTrueFace_V(AlfvenWavePlusFirst_:AlfvenWavePlusLast_)
          ! Set spectrum of "Minus" waves
          VarsGhostFace_V(AlfvenWaveMinusFirst_:AlfvenWaveMinusLast_) = &
               WavesGhostFace_I(1:nWaveHalf)
      
       else
   
          ! full absorption of "Minus" waves - float bc
          VarsGhostFace_V(AlfvenWaveMinusFirst_:AlfvenWaveMinusLast_) =&
               VarsTrueFace_V(AlfvenWaveMinusFirst_:AlfvenWaveMinusLast_)
          ! Set spectrum of "Plus" waves
          VarsGhostFace_V(AlfvenWavePlusFirst_:AlfvenWavePlusLast_) = &
               WavesGhostFace_I
   
       end if
       if(any(VarsGhostFace_V(WaveFirst_:WaveLast_) <= 0.0)) then
          write(*,*) 'Negative wave energy at inner boundary'
          call stop_MPI('ERROR in plus waves BC')
       end if
     
       !\
       ! BC for total wave pressure
       !/
       if(UseWavePressureLtd) then
          VarsGhostFace_V(Ew_) = 0.5*sum(VarsGhostFace_V(WaveFirst_:WaveLast_))
       else
          VarsGhostFace_V(Ew_) = 1.0e-30
       end if
    end if

    !\
    ! Apply corotation
    !/
    if (.not.UseRotatingFrame) then
       VarsGhostFace_V(Ux_) = VarsGhostFace_V(Ux_) &
            - 2*OmegaBody*FaceCoords_D(y_)
       VarsGhostFace_V(Uy_) = VarsGhostFace_V(Uy_) &
            + 2*OmegaBody*FaceCoords_D(x_)
    end if

  end subroutine user_face_bcs
  !===========================================================================
  subroutine get_wave_group_bc(wEnergyDensBc,WavesGhostFace_I)
   
  ! This subroutine decomposes the total wave energy into individual frequency
  ! groups according to the spectral index
  !
  ! Note that this subroutine is designed for the inner boundary. The actual
  ! emmitted frequency range can be narrower than frequency grid. 
    
    real, intent(out)           :: WavesGhostFace_I(1:nWaveHalf)
    real,intent(in)             :: wEnergyDensBc
    integer                     :: iWave
    character(len=*),parameter  :: NameSub = "get_wave_group_bc"
    !--------------------------------------------------------------------------
    do iWave = 1, nWaveHalf
       if(iWave .lt. WaveFirstBc_ .or. iWave .gt. WaveLastBc_) then
          ! no wave energy outside of emitted freq. range
          WavesGhostFace_I(iWave)=1e-30
       else
          ! spectral decomposition
          WavesGhostFace_I(iWave) = SpectralCoeff * dLogFreq * wEnergyDensBc* &
               exp((LogFreq_I(iWave-WaveFirst_+1))*SpectralIndex)
       end if
    end do
  end subroutine get_wave_group_bc
  !===========================================================================
  subroutine get_plasma_parameters_cell(iCell,jCell,kCell,iBlock,&
       DensCell,PresCell)
     
    ! This subroutine computes the cell values for density and pressure 
    ! assuming an isothermal atmosphere
    
    use ModGeometry,   ONLY: x_BLK,y_BLK,z_BLK,R_BLK
    use ModNumConst
    use ModPhysics,    ONLY: GBody,BodyRho_I,Si2No_V,UnitTemperature_
    use ModExpansionFactors,  ONLY: UMin,CoronalT0Dim
    implicit none

    integer, intent(in)  :: iCell,jCell,kCell,iBlock
    real, intent(out)    :: DensCell,PresCell
    real :: UFinal       !The solar wind speed at the far end of the Parker spiral,
                         !which originates from the given cell
    real :: URatio       !The coronal based values for temperature density 
                         !are scaled as functions of UFinal/UMin ratio
    real :: Temperature
    !--------------------------------------------------------------------------

    !call get_gamma_emp(x_BLK(iCell,jCell,kCell,iBlock),&
    !     y_BLK(iCell,jCell,kCell,iBlock),&
    !     z_BLK(iCell,jCell,kCell,iBlock),&
    !     GammaCell)
    call get_bernoulli_integral(x_BLK(iCell,jCell,kCell,iBlock)/&
         R_BLK(iCell,jCell,kCell,iBlock),&
         y_BLK(iCell,jCell,kCell,iBlock)/R_BLK(iCell,jCell,kCell,iBlock),&
         z_BLK(iCell,jCell,kCell,iBlock)/R_BLK(iCell,jCell,kCell,iBlock),UFinal)
    URatio=UFinal/UMin

    !This is the temperature variation
    Temperature = CoronalT0Dim*Si2No_V(UnitTemperature_)/(min(URatio,2.0))

    DensCell  = (1.0/URatio) &          !This is the density variation
         *BodyRho_I(1)*exp(-GBody/Temperature &
         *(1.0/max(R_BLK(iCell,jCell,kCell,iBlock),0.90)-1.0))

    PresCell = DensCell*Temperature

  end subroutine get_plasma_parameters_cell
  !============================================================================
  subroutine user_set_ics

    use ModMain,      ONLY: globalBLK
    use ModAdvance,   ONLY: State_VGB 
    use ModPhysics,   ONLY: inv_gm1,BodyTDim_I
    use ModGeometry
    implicit none

    integer :: i,j,k,iBLK
    logical :: oktest,oktest_me
    real    :: Dens_BLK,Pres_BLK
    real    :: x,y,z,R,ROne,Rmax,U0
    !--------------------------------------------------------------------------
    call set_oktest('user_set_ics',oktest,oktest_me)

    iBLK = globalBLK

    select case(TypeGeometry)
    case('cartesian')
       Rmax = max(2.1E+01,sqrt(x2**2+y2**2+z2**2))
    case('spherical_lnr')
       Rmax = max(2.1E+01,exp(XyzMax_D(1)))
    end select

    ! The sqrt is for backward compatibility with older versions of the Sc
    U0 = 4.0*sqrt(2.0E+6/BodyTDim_I(1))
    
    State_VGB(:,1:nI,1:nJ,1:nK,iBLK) = 1.0e-30 !Initialize the wave spectrum
    do k=1,nK; do j=1,nJ; do i=1,nI
       x = x_BLK(i,j,k,iBLK)
       y = y_BLK(i,j,k,iBLK)
       z = z_BLK(i,j,k,iBLK)
       R = max(R_BLK(i,j,k,iBLK),cTolerance)
       ROne = max(1.0,R)
       State_VGB(Bx_:Bz_,i,j,k,iBLK) = 0.0
       call get_plasma_parameters_cell(i,j,k,iBLK,&
            Dens_BLK,Pres_BLK)
       State_VGB(rho_,i,j,k,iBLK) = Dens_BLK
       State_VGB(P_,i,j,k,iBLK)   = Pres_BLK
       State_VGB(RhoUx_,i,j,k,iBLK) = Dens_BLK &
            *U0*((ROne-1.0)/(Rmax-1.0))*x/R
       State_VGB(RhoUy_,i,j,k,iBLK) = Dens_BLK &
            *U0*((ROne-1.0)/(Rmax-1.0))*y/R
       State_VGB(RhoUz_,i,j,k,iBLK) = Dens_BLK &
            *U0*((ROne-1.0)/(Rmax-1.0))*z/R
       State_VGB(Ew_,i,j,k,iBLK) = nWave*1.0e-30

    end do; end do; end do

  end subroutine user_set_ics
  !============================================================================
  subroutine user_get_b0(xInput,yInput,zInput,B0_D)
 
    use ModPhysics,     ONLY: Si2No_V,UnitB_
    use ModMagnetogram, ONLY: get_magnetogram_field

    real, intent(in):: xInput,yInput,zInput
    real, intent(out), dimension(3):: B0_D
    !---------------------------------------------

    call get_magnetogram_field(xInput,yInput,zInput,B0_D)
    B0_D = B0_D*Si2No_V(UnitB_)

   end subroutine user_get_b0
  !===========================================================================
  subroutine user_update_states(iStage,iBlock)

    use ModAdvance, ONLY: State_VGB, B0_DGB

    use ModGeometry,ONLY: R_BLK
    use ModEnergy,  ONLY: calc_energy_cell
    use ModWaves,   ONLY: UseWavePressure, UseWavePressureLtd, &
                          AlfvenWavePlusFirst_,  AlfvenWavePlusLast_, &
                          AlfvenWaveMinusFirst_, AlfvenWaveMinusLast_ 

    integer,intent(in)           :: iStage,iBlock
    integer                      :: i,j,k
    real                         :: DensCell,PresCell
    character(len=*),parameter   :: NameSub='user_update_states'
    !--------------------------------------------
    !\
    ! Advect solution
    !/
    if(any(State_VGB(WaveFirst_:WaveLast_,1:nI,1:nJ,1:nK,iBlock)<0.0)) then
       write(*,*) NameSub,' : negative wave energy before MHD'
    end if
    
    call update_states_MHD(iStage,iBlock)
    
    if(any(State_VGB(WaveFirst_:WaveLast_,1:nI,1:nJ,1:nK,iBlock)<0.0)) then
       write(*,*) NameSub, ': negative wave energy after MHD'
    end if

    if (UseWavePressure .and. (.not. UseWavePressureLtd)) &
         State_VGB(Ew_,i,j,k,iBlock) = nWave*1.0e-30

    !\
    ! Dissipate wave energy after advection
    !/
    if (IsInitWave .and. DoDampCutOff) call dissipate_waves(iBlock)
   
    !call calc_energy_cell(iBlock)
 
  end subroutine user_update_states
  !=======================================================================
   subroutine dissipate_waves(iBlock)

    use ModAdvance, ONLY: State_VGB
    use ModPhysics, ONLY: Gamma0

    integer,intent(in)      :: iBlock
    integer                 :: i, j, k, iWave
    real                    :: dWavePres = 0.0
    real                    :: LogFreqCutOff
    character(len=*),parameter :: NameSub="dissipate_waves"
    ! -----------------------------------------------------------------
    if (DoDampCutOff) then
       do k=1,nK ; do j=1, nJ ; do i=1,nI
          call calc_cutoff_freq(i,j,k,iBlock,LogFreqCutOff)
          do iWave=WaveFirst_, WaveLast_
             if (LogFreq_I(iWave) .ge. LogFreqCutOff .and. &
                  State_VGB(iWave,i,j,k,iBlock) > 0.0) then
                dWavePres = -State_VGB(iWave,i,j,k,iBlock)
             end if
          end do
       end do; end do; end do
    end if

    ! Add more dissipation mechanisms here
  

    ! Remove dissipated energy density from spectrum
    State_VGB(iWave,i,j,k,iBlock) = &
         State_VGB(iWave,i,j,k,iBlock) + dWavePres

    ! pass dissipated wave energy density to the MHD pressure
    State_VGB(p_,i,j,k,iBlock) = State_VGB(p_,i,j,k,iBlock) - &
        0.5*(Gamma0-1)*dWavePres

  end subroutine dissipate_waves
  !=====================================================================
   subroutine write_cell_spectrum
    
    use ModProcMH
    use ModMain,     ONLY: iteration_number,unusedBLK,nBlockALL
    use ModGeometry, ONLY: x_BLK,y_BLK,z_BLK
    use ModIoUnit,   ONLY: io_unit_new
    use ModAdvance,  ONLY: State_VGB
    use ModPhysics,  ONLY: No2Si_V, UnitX_, UnitP_
    use ModWaves,    ONLY: AlfvenWaveMinusFirst_,AlfvenWavePlusFirst_

    implicit none
    
    real,dimension(nWaveHalf,3)      :: Spectrum_II 
    real                             :: x,y,z, IwPlusSi,IwMinusSi
    integer                          :: iFreq,i,j,k,iBLK
    integer                          :: iUnit,iError,aError
    logical                          :: DoSaveCellSpec = .false.
    character(len=40)                :: FileNameTec 
    character(len=11)                :: NameStage
    character(len=7)                 :: NameProc
    character(len=*),parameter       :: NameSub='write_cell_spectrum'
    !-------------------------------------------------------------------
   
    do iBLK=1,nBLK
       if(unusedBLK(iBLK)) CYCLE
       do k=1,nK ; do j=1,nJ ; do i=1,nI
          x=x_BLK(i,j,k,iBLK)
          y=y_BLK(i,j,k,iBLK)
          z=z_BLK(i,j,k,iBLK)
          if((x == xTestSpec) .and. (y == yTestSpec) .and. (z == zTestSpec)) then
             DoSaveCellSpec = .true.
             do iFreq=1,nWaveHalf
                IwPlusSi  = No2Si_V(UnitP_)* &
                     State_VGB(AlfvenWavePlusFirst_+iFreq-1,i,j,k,iBLK)
                IwMinusSi = No2Si_V(UnitP_)* &
                     State_VGB(AlfvenWaveMinusFirst_+iFreq-1,i,j,k,iBLK)
                Spectrum_II(iFreq,1) = LogFreq_I(iFreq)
                Spectrum_II(iFreq,2) = IwPlusSi
                Spectrum_II(iFreq,3) = IwMinusSi
                
             end do
          end if
       end do; end do ; end do
    end do
    
    if (DoSaveCellSpec) then
       !\
       ! write data file
       !/
       write(NameStage,'(i6.6)') iteration_number
       write(NameProc,'(a,i4.4)') "_pe",iProc
       FileNameTec = 'SC/IO2/Cell_Spectrum_n_'//trim(NameStage)//trim(NameProc)//'.spec'
       if(iProc == 0) write(*,*) 'SC: writing file ',FileNameTec
      
       iUnit=io_unit_new()
       open(unit=iUnit, file=FileNameTec, form='formatted', access='sequential',&
            status='replace',iostat=iError)
       write(iUnit, '(a)') 'Title: BATSRUS SC Spectrogram'
       write(iUnit, '(a)') 'Variables = "log(w)","I+[Jm-3]","I-[Jm-3]" '
       write(iUnit,'(a,i3.3,a,i5.5,a)') 'Zone I= ',nWaveHalf,' F=point'
       do iFreq=1,nWaveHalf
          write(iUnit, fmt="(30(e14.6))") &
               Spectrum_II(iFreq,:)
       end do
       close(iUnit)
    end if
     
  end subroutine write_cell_spectrum
   !=====================================================================
   subroutine write_spectrogram
    
    use ModProcMH
    use ModMain,      ONLY: iteration_number,unusedBLK,nBlockALL
    use ModGeometry,  ONLY: x_BLK,y_BLK, z_BLK, dz_BLK,dx_BLK
    use ModIoUnit,    ONLY: io_unit_new
    use ModAdvance,   ONLY: State_VGB
    use ModPhysics,   ONLY: No2Si_V, UnitX_, UnitP_
    use ModWaves,     ONLY: AlfvenWaveMinusFirst_,AlfvenWavePlusFirst_

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
    write(NameStage,'(i6.6)') iteration_number
    write(NameProc,'(a,i4.4)') "_pe",iProc
    FileNameTec = 'SC/IO2/Spectrum_n_'//trim(NameStage)//trim(NameProc)//'.tec'
    if(iProc == 0) write(*,*) 'SC: writing file ',FileNameTec
    !\
    ! count cells along x=xTrace, z=zTrace line parallel to y axis
    !/
    nCell=0
    do iBLK=1,nBLK
       if(unusedBLK(iBLK)) CYCLE
       do k=1,nK ; do j=1,nJ ; do i=1,nI
          x=x_BLK(i,j,k,iBLK)
          dx=dx_BLK(iBLK)
          z=z_BLK(i,j,k,iBLK)
          dz=dz_BLK(iBLK)
          if((z < zTrace+dz) .and. (z >= zTrace) .and. &
             (x < xTrace+dx) .and. (x >= xTrace)) then
             nCell=nCell+1
          end if
       end do; end do ; end do
    end do
    if (nCell > 0 )then !write spectrogram output file
     
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
             if(unusedBLK(iBLK)) CYCLE
             do k=1,nK ; do j=1,nJ ; do i=1,nI
                x=x_BLK(i,j,k,iBLK)
                y=y_BLK(i,j,k,iBLK)
                z=z_BLK(i,j,k,iBLK)
                dx=dx_BLK(iBLK)
                dz=dz_BLK(iBLK)
                if((z < zTrace+dz) .and. (z >= zTrace) .and. &
                   (x < xTrace+dx) .and. (x >= xTrace)) then
                   do iFreq=1,nWaveHalf
                      IwPlusSi  = No2Si_V(UnitP_)* &
                           State_VGB(AlfvenWavePlusFirst_+iFreq-1,i,j,k,iBLK)
                      IwMinusSi = No2Si_V(UnitP_)* &
                           State_VGB(AlfvenWaveMinusFirst_+iFreq-1,i,j,k,iBLK)
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
       iUnit=io_unit_new()
       open(unit=iUnit, file=FileNameTec, form='formatted', access='sequential',&
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
    end if
  end subroutine write_spectrogram
  !=====================================================================
  subroutine calc_cutoff_freq(i,j,k,iBLK , LogFreqCutOff)

    ! This subroutine calculates the cut-off frequency for Alfven waves, 
    ! which is the ion cyclotron frequency (in radians) 
    ! \omega_{c.o.}=zeB/m. For ions z=1. 
     
    use ModMain,                ONLY: x_,y_,z_
    use ModAdvance,             ONLY: State_VGB, B0_DGB
    use ModConst,               ONLY: cElectronCharge, cProtonMass,cTiny
    use ModPhysics,             ONLY: No2Si_V,UnitB_

    real, intent(out)           :: LogFreqCutOff
    integer, intent(in)         :: i,j,k,iBLK
    real                        :: BtotSi
    real,dimension(3)           :: FullB_D
    character(len=*),parameter  :: NameSub = 'calc_cutoff_freq'
    ! -----------------------------------------------------------------
    FullB_D = State_VGB(Bx_:Bz_,i,j,k,iBLK) + B0_DGB(x_:z_,i,j,k,iBLK)
   
    BtotSi = No2Si_V(UnitB_)*sqrt(sum(FullB_D**2))
    if(BtotSi <=0.0) then
       write(*,*) 'Btot negative, ',BtotSi
    end if
    LogFreqCutOff = log((cElectronCharge*BtotSi)/cProtonMass)

  end subroutine calc_cutoff_freq
  !=======================================================================
  subroutine calc_poynt_flux(i,j,k,iBLK, UseUr, PoyntFluxSi)

    Use ModMain,       ONLY: x_,y_,z_
    use ModAdvance,    ONLY: State_VGB, B0_DGB
    use ModGeometry,   ONLY: x_BLK,y_BLK,z_BLK,R_BLK
    use ModPhysics,    ONLY: No2Si_V, UnitB_, UnitRho_,UnitU_,UnitP_
    use ModConst,      ONLY: cMu
    implicit none
    
    integer,intent(in)         :: i,j,k,iBLK
    logical,intent(in)         :: UseUr
    real,intent(out)           :: PoyntFluxSi
    real,dimension(3)          :: r_D,BSi_D,USi_D
    real                       :: x,y,z,BrSi,UrSi, RhoSi  
    real                       :: vAlfvenRadialSi ! in radial direction
    character(len=*),parameter :: NameSub='calc_poynt_flux'
    ! -----------------------------------------------------------------
    !  Poynting flux is calculated in SI UNITS
    !  This subroutine is used for finding flux leaving a spherical surface
    ! (radius depends on calling routine), thus only the radial component is calculated.
    
    ! get radial unit vector
    x = x_BLK(i,j,k,iBLK)
    y = y_BLK(i,j,k,iBLK)
    z = z_BLK(i,j,k,iBLK)
    r_D = (/x,y,z/)
    r_D = r_D/sqrt(sum(r_D**2)) ! unit vector in radial direction

    ! get B, Br vectors in SI units
    BSi_D = (State_VGB(Bx_:Bz_,i,j,k,iBLK) + B0_DGB(x_:z_,i,j,k,iBLK)) &
         * No2Si_V(UnitB_)
    BrSi = sum(r_D*BSi_D)
    
    ! get density in SI units
    RhoSi = No2Si_V(UnitRho_)*State_VGB(rho_,i,j,k,iBLK)
    ! calculate radial component of Alfven speed
    vAlfvenRadialSi = abs(BrSi/sqrt(RhoSi*cMu)) ! for outgoing waves
    
    ! add Ur to vAlfven if bulk velocity is not neglected
    if(UseUr) then
       USi_D = No2Si_V(UnitU_)*State_VGB(Ux_:Uz_,i,j,k,iBLK)
       UrSi = sum(r_D*USi_D)
        
       vAlfvenRadialSi = UrSi+vAlfvenRadialSi
    end if
     
    PoyntFluxSi=vAlfvenRadialSi*No2Si_V(UnitP_)* & 
         sum(State_VGB(WaveFirst_:WaveLast_,i,j,k,iBLK))

  end subroutine calc_poynt_flux
  !========================================================================
  subroutine user_get_log_var(VarValue,TypeVar,Radius)
    
    use ModIO,         ONLY: write_myname
    use ModMain,       ONLY: unusedBLK,x_,y_,z_
    use ModVarIndexes !ONLY: Ew_,Bx_,By_,Bz_,rho_,rhoUx_,rhoUy_,rhoUz_,P_ 
    use ModGeometry,   ONLY: R_BLK
    use ModAdvance,    ONLY: State_VGB,tmp1_BLK,B0_DGB
    use ModPhysics,    ONLY: inv_gm1,No2Si_V,UnitEnergydens_, &
                             UnitX_,UnitU_,UnitRho_

    real, intent(out)              :: VarValue
    character (LEN=10), intent(in) :: TypeVar 
    real, intent(in), optional     :: Radius
    integer                        :: iBLK
    real                           :: unit_energy,unit_mass
    real, external                 :: integrate_BLK
    !--------------------------------------------------------------------------
    unit_energy = 1.0e7*No2Si_V(UnitEnergydens_)*No2Si_V(UnitX_)**3
    unit_mass   = 1.0e3*No2Si_V(UnitRho_)*No2Si_V(UnitX_)**3
    !\
    ! Define log variable to be saved::
    !/
    select case(TypeVar)
    case('spec')
       call write_spectrogram
    case('cellspec')
       call write_cell_spectrum
    case('em_t','Em_t','em_r','Em_r')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(:,:,:,iBLK) = & 
               (B0_DGB(x_,:,:,:,iBLK)+State_VGB(Bx_,:,:,:,iBLK))**2+&
               (B0_DGB(y_,:,:,:,iBLK)+State_VGB(By_,:,:,:,iBLK))**2+&
               (B0_DGB(z_,:,:,:,iBLK)+State_VGB(Bz_,:,:,:,iBLK))**2
       end do
       VarValue = unit_energy*0.5*integrate_BLK(1,tmp1_BLK)
    case('ek_t','Ek_t','ek_r','Ek_r')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(:,:,:,iBLK) = &
               (State_VGB(rhoUx_,:,:,:,iBLK)**2 +&
               State_VGB(rhoUy_,:,:,:,iBLK)**2 +&
               State_VGB(rhoUz_,:,:,:,iBLK)**2)/&
               State_VGB(rho_  ,:,:,:,iBLK)             
       end do
       VarValue = unit_energy*0.5*integrate_BLK(1,tmp1_BLK)
    case('et_t','Et_t','et_r','Et_r')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(:,:,:,iBLK) = State_VGB(P_,:,:,:,iBLK)
       end do
       VarValue = unit_energy*inv_gm1*integrate_BLK(1,tmp1_BLK)
    case('ew_t','Ew_t','ew_r','Ew_r')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(:,:,:,iBLK) = State_VGB(Ew_,:,:,:,iBLK)
       end do
       VarValue = unit_energy*integrate_BLK(1,tmp1_BLK)
    case('ms_t','Ms_t')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(:,:,:,iBLK) = &
               State_VGB(rho_,:,:,:,iBLK)/R_BLK(:,:,:,iBLK)
       end do
       VarValue = unit_mass*integrate_BLK(1,tmp1_BLK)
    case('vol','Vol')
       tmp1_BLK(:,:,:,iBLK) = 1.0
       VarValue = integrate_BLK(1,tmp1_BLK)
    case default
       VarValue = -7777.
       call write_myname;
       write(*,*) 'Warning in set_user_logvar: unknown logvarname = ',TypeVar
    end select
  end subroutine user_get_log_var
 !===========================================================================
  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModAdvance, ONLY: State_VGB
    use ModPhysics, ONLY: NameTecUnit_V, NameIdlUnit_V, &
         No2Si_V,UnitPoynting_,UnitP_,UnitX_,UnitEnergyDens_

    integer,          intent(in)   :: iBlock
    character(len=*), intent(in)   :: NameVar
    logical,          intent(in)   :: IsDimensional
    real,             intent(out)  :: PlotVar_G(-1:nI+2, -1:nJ+2, -1:nK+2)
    real,             intent(out)  :: PlotVarBody
    logical,          intent(out)  :: UsePlotVarBody
    character(len=*), intent(inout):: NameTecVar
    character(len=*), intent(inout):: NameTecUnit
    character(len=*), intent(inout):: NameIdlUnit
    logical,          intent(out)  :: IsFound

    character (len=*), parameter :: NameSub = 'user_set_plot_var'
    real                         :: PoyntFlux, UnitEnergy
    integer                      :: i,j,k, I_Index
    logical                      :: IsError
    !-------------------------------------------------------------------    
    !UsePlotVarBody = .true. 
    !PlotVarBody = 0.0 
    IsFound=.true.

    UnitEnergy=1.0e7*No2Si_V(UnitEnergydens_)*No2Si_V(UnitX_)**3
    !\                                                                              
    ! Define plot variable to be saved::
    !/ 
    !
    select case(NameVar)
       !Always use lower case !!
    case('wpres')
       do k=-1,nK+2 ; do j=-1,nJ+2 ; do i=-1,nI+2
          PlotVar_G(i,j,k) = 0.5*sum(State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock))
       end do ; end do ; end do
       PlotVar_G=No2Si_V(UnitP_)*PlotVar_G
       NameTecVar = 'wPres'
       NameTecUnit = NameTecUnit_V(UnitP_)
       NameIdlUnit = NameIdlUnit_V(UnitP_)
    !case('wdiss')
    !   do k=-1,nK+2 ; do j=-1,nJ+2 ; do  i=-1,nI+2
    !      PlotVar_G(i,j,k)=WaveDissip_CB(i,j,k,iBlock)
    !   end do ; end do ; end do
    !   PlotVar_G=UnitEnergy*PlotVar_G
    !   NameTecVar = 'wDissip'
    case('poynt')
       ! neglect radial bulk velociy
       do i=-1,nI+2 ; do j=-1,nJ+2; do k=-1,nK+2
          call calc_poynt_flux(i,j,k,iBlock,.false.,PoyntFlux)
          PlotVar_G(i,j,k)=PoyntFlux
       end do; end do ; end do
       NameTecVar = 'S_r(u=0)'
       NameTecUnit = NameTecUnit_V(UnitPoynting_)
       NameIdlUnit = NameIdlUnit_V(UnitPoynting_)
    case('poyntur')
       ! include radial bulk velocity
       do i=-1,nI+2; do j=-1,nJ+2; do k=-1,nK+2
          call calc_poynt_flux(i,j,k,iBlock,.true.,PoyntFlux)
          PlotVar_G(i,j,k)=PoyntFlux
       end do; end do ; end do
       NameTecVar = 'S_r(total)'
       NameTecUnit = NameTecUnit_V(UnitPoynting_)
       NameIdlUnit = NameIdlUnit_V(UnitPoynting_)
    case default
       IsFound= .false.
    end select
  end subroutine user_set_plot_var
  !=========================================================================== 
  subroutine user_specify_refinement(iBlock, iArea, DoRefine)

    use ModAdvance,  ONLY: State_VGB, Bx_, By_, Bz_, B0_DGB
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, far_field_BCs_BLK
    use ModNumConst, ONLY: cTiny
    use ModMain,     ONLY: x_, y_, z_

    integer, intent(in) :: iBlock, iArea
    logical,intent(out) :: DoRefine

    real :: rDotB_G(1:nI,1:nJ,0:nK+1)
    integer :: i,j,k
    character (len=*), parameter :: NameSub = 'user_specify_refinement'
    !-------------------------------------------------------------------

    if(far_field_BCs_BLK(iBlock))then
       DoRefine = .false.
       RETURN
    end if

    ! Calculate r.B in all physical cells and ghost cells 
    ! in the Z/Theta direction to find current sheet 
    ! passing between blocks
    do k=0, nK+1; do j=1, nJ; do i=1, nI
       rDotB_G(i,j,k) = x_BLK(i,j,k,iBlock)   &
            * (B0_DGB(x_,i,j,k,iBlock) + State_VGB(Bx_,i,j,k,iBlock)) &
            +              y_BLK(i,j,k,iBlock)   &
            * (B0_DGB(y_,i,j,k,iBlock) + State_VGB(By_,i,j,k,iBlock)) &
            +              z_BLK(i,j,k,iBlock)   &
            * (B0_DGB(z_,i,j,k,iBlock) + State_VGB(Bz_,i,j,k,iBlock))
    end do; end do; end do;

    DoRefine = maxval(rDotB_G) > cTiny .and. minval(rDotB_G) < -cTiny

  end subroutine user_specify_refinement

  !============================================================================
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

