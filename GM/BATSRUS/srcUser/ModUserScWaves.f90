
!^CFG COPYRIGHT UM
!==============================================================================
Module ModUser

! User Module for Solar Corona with a spectrum of Alfven waves
! Description (Oran, Jun 20, 2010)

! The following user routines are implemented:
! --------------------------------------------
! 1. user_read_inputs : read input commands for wave related quantities
!
!  #INITSPECTRUM : allows the user to set the logical flag IsInitWave. 
!                  If true, the wave spectrum will be initialized. Should be used only in the
!                  first session of the simulation.
!  #WAVEINNERBC  : allows the user to choose type of inner boundary condition for Alfven waves.
!                  TypeWaveInnerBc (string) -  Two options are implemented:
!                  'WSA' - uses the WSA model to get the solar wind velocity at 1AU, then calculates
!                          the total wave energy at coronal base from the Bernoulli integral.
!                  'turb' - calculates the total wave energy from the local magnetic field.
!                  WaveInnerbcFactor (real) -  allows scaling of the resulting total 
!                  wave energy by the same factor for all boundary faces.
!  #DAMPWAVES    : allows the user to choose frequency dependent wave damping mechanism.
!                  DoDampCutoff (logical) - if true, Alfven waves will be totaly damped at
!                                             and above the local ion cyclotron frequency.
!                  DoDampSurface (logical) - if true, Alfven waves will be damped according
!                                            to the local surface waves dissipation length.
!  #SPECTROGRAM  : allows the user to choose a line parallel to one of the axes along which
!                  the full wave spectrum will be extracted and written to a file for plotting.
!                  >>>
!                  xTrace, yTrace, zTrace are the coordinates of the line
!  #CELLSPECTRUM : allows the user to output the spectrum in a single cell for testing purposes.
! -----------------------------------------------------------------------------
! 2. user_set_ic : intilize MHD parameters according to an isothermal atmosphere solution.
! -----------------------------------------------------------------------------
! 3. user_face_bc : set inner boundary conditions for MHD variables and Alfven waves.
!                   Magnetic field B1 is set to equal its tangential component.
!                   Reflective boundary conditions for velocity.
!                   Density and pressure set according to isothermal atmosphere.
!                   Alfven Waves: wave energy is present in the inner boundary only
!                                 in the open field lines region (to avoid blowing out
!                                 of Helmet streamers). The total wave energy at the face
!                                 is calculated according to the type of boundary condition
!                                 (see desctription of #WAVEINNERBC in section 1 above).   
!                                 Once the total wave energy is obtained, it is deconstructed
!                                 into the frequency bins according to the assumed spectral shape
!                                 by calling set_wave_state (in ModWaves.f90).
!                   Rotating frame : if the non rotating frame is used, the x and y velocity components
!                                    are adjusted.
! -----------------------------------------------------------------------------                            
! 4. user_get_log_var : allows the user to output spectral data to a dedicated log file.
!                       In order to use this option, the user must include the command #SAVELOGFILE
!                       in the PARAM.in file, with the following format:
!                       #SAVELOGFILE
!                       T
!                       var
!                       dn
!                       dt
!                       'string'
!
!                       where 'string' can be set to 'spectrogram' or 'cellspecrum' (refer to the the SWMF
!                       user manual for further details on this command).
!                       'spectrogram' : the full wave spectrum along a line is extracted
!                       according to the parameters set in #SPECTROGRAM (see above). The data is
!                       written to a file called Spectrum_n_xxx_peyyy.tec, where xxx stands
!                       for simulation time/iteration number and yyy to the processor number. 
!                       Files from different processors should be combined later.
!                       
!                       'cellspectrum' : the full spectrum in a single cell is written to a file
!                       The cell is chosen by its x,y,z coordinates set in #CELLSPECTRUM (see above).
!                       The output is written to a file Cell_Spectrum_n_xxx_peyyy.tec. Since only
!                       one processor is involved, no post-processing is necessary.
! -----------------------------------------------------------------------------------------------
! 5. user_update_states : calls update states_MHD, followed by a call to dissipate_waves, which
!                         removes wave energy from the spectrum and adds it to the plasma pressure.
! -----------------------------------------------------------------------------------------------
! 6. user_specify_refinement : improve refinement in the current sheet region.
! -----------------------------------------------------------------------------------------------
! 7. user_set_boundary_cells :required when "extra" boundary conditions are used.
! -----------------------------------------------------------------------------------------------
! 8. user_set_plot_var : implement plot varaibles related to Alfven waves.
!                        'wpres' - total wave pressure (summed over spectrum).
!                        'poynt' - Poynting flux of Alfven waves neglecting plasma speed.
!                        'poyntur' - Poynting flux of Alfven waves.
! -----------------------------------------------------------------------------------------------
  use ModMain, ONLY: nBLK
  use ModSize, ONLY: nI,nJ,nK
  use ModVarIndexes
  use ModWaves,  ONLY: DeltaLogFrequency
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_set_ics,                    &
       IMPLEMENTED3 => user_face_bcs,                   &
       IMPLEMENTED4 => user_get_log_var,                &
       IMPLEMENTED5 => user_update_states,              &
       IMPLEMENTED6 => user_specify_refinement,         &
       IMPLEMENTED7 => user_set_boundary_cells,         &
       IMPLEMENTED8 => user_set_plot_var
  
  include 'user_module.h' !list of public methods

  real, parameter               :: VersionUserModule = 1.0
  character (len=*), parameter  :: NameUserModule = 'Alfven Waves Driven Solar Wind'
  
  ! varaibles read by user_read_inputs and used by other subroutines
  logical                       :: IsInitWave = .false.
  logical                       :: DoDampCutOff = .false., DoDampSurface = .false.
  real                          :: WaveInnerBcFactor
  real                          :: xTrace = 0.0, zTrace = 0.0
  real                          :: xTestSpectrum, yTestSpectrum, zTestSpectrum
  character(len=10)             :: TypeWaveInnerBc

contains
  !============================================================================
  subroutine user_read_inputs

    !  user_read_inputs : read input commands for wave related quantities
    !
    !  #INITSPECTRUM : allows the user to set the logical flag IsInitWave. 
    !                  If true, the wave spectrum will be initialized. Should be used only in the
    !                  first session of the simulation.
    !  #WAVEINNERBC  : allows the user to choose type of inner boundary condition for Alfven waves.
    !                  TypeWaveInnerBc (string) -  Two options are implemented:
    !                  'WSA' - uses the WSA model to get the solar wind velocity at 1AU, then calculates
    !                          the total wave energy at coronal base from the Bernoulli integral.
    !                  'turb' - calculates the total wave energy from the local magnetic field.
    !                  WaveInnerbcFactor (real) -  allows scaling of the resulting total 
    !                  wave energy by the same factor for all boundary faces.
    !  #DAMPWAVES    : allows the user to choose frequency dependent wave damping mechanism.
    !                  DoDampCutoff (logical) - if true, Alfven waves will be totaly damped at
    !                                             and above the local ion cyclotron frequency.
    !                  DoDampSurface (logical) - if true, Alfven waves will be damped according
    !                                            to the local surface waves dissipation length.
    !  #SPECTROGRAM  : allows the user to choose a line parallel to one of the axes along which
    !                  the full wave spectrum will be extracted and written to a file for plotting.
    !                  >>>
    !                  xTrace, yTrace, zTrace are the coordinates of the line
    !  #CELLSPECTRUM : allows the user to output the spectrum in a single cell for testing purposes.
    ! -----------------------------------------------------------------------------
    use ModMain,        ONLY: UseUserInitSession, lVerbose
    use ModProcMH,      ONLY: iProc
    use ModReadParam,   ONLY: read_line, read_command, read_var
    use ModIO,          ONLY: write_prefix, write_myname, iUnitOut, NamePlotDir
 
    character (len=100) :: NameCommand
    !--------------------------------------------------------------------------
    !UseUserInitSession = .true.

    if(iProc == 0 .and. lVerbose > 0)then
       call write_prefix;
       write(iUnitOut,*)'User read_input TURBULENCE CORONA starts'
    endif

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)
          
       case("#INITSPECTRUM")
          call read_var('IsInitWave', IsInitWave)
       
       case("#WAVEINNERBC")
          call read_var('TypeWaveInnerBc', TypeWaveInnerBc)
          call read_var('WaveInnerBcFactor', WaveInnerBcFactor)
         
       case("#DAMPWAVES")
          call read_var('DoDampCutoff', DoDampCutoff)
          call read_var('DoDampSurface', DoDampSurface)

       case("#SPECTROGRAM")
          call read_var('xTrace', xTrace)
          call read_var('zTrace', zTrace)
       case("#CELLSPECTRUM")
          call read_var('xTestSpectrum', xTestSpectrum)
          call read_var('yTestSpectrum', yTestSpectrum)
          call read_var('zTestSpectrum', zTestSpectrum)
      
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
  subroutine user_face_bcs(VarsGhostFace_V)
   
    ! user_face_bc : set inner boundary conditions for MHD variables and Alfven waves.
    !                   Magnetic field B1 is set to equal to its tangential component.
    !                   Reflective boundary conditions for velocity.
    !                   Density and pressure set according to isothermal atmosphere.
    !                   Alfven Waves: wave energy is present in the inner boundary only
    !                                 in the open field lines region (to avoid blowing out
    !                                 of Helmet streamers). The total wave energy at the face
    !                                 is calculated according to the type of boundary condition
    !                                 (see desctription of #WAVEINNERBC in section 1 above).   
    !                                 Once the total wave energy is obtained, it is deconstructed
    !                                 into the frequency bins according to the assumed spectral shape
    !                                 by calling set_wave_state (in ModWaves.f90).
    !                   Rotating frame : if the non rotating frame is used, the x and y velocity components
    !                                    are adjusted.
    ! -----------------------------------------------------------------------------                            
    use ModSize,             ONLY: East_, West_, South_, North_, Bot_, Top_, nDim
    use ModMain,             ONLY: x_, y_, z_, UseRotatingFrame
    use ModExpansionFactors, ONLY: ExpansionFactorInv_N, get_interpolated
    use ModAdvance,          ONLY: State_VGB
    use ModPhysics,          ONLY: OmegaBody, No2Si_V, Si2No_V, UnitU_, UnitP_
    use ModNumConst,         ONLY: cTolerance
    use ModFaceBc,           ONLY: FaceCoords_D, VarsTrueFace_V, B0Face_D, &
                                   iFace, jFace, kFace, iSide, iBlockBc
    use ModWaves,            ONLY: set_wave_state
    
    real, intent(out)           :: VarsGhostFace_V(nVar)



    integer                     :: iCell, jCell, kCell
    real                        :: DensCell, PresCell, TBase, TotalB  
    real, dimension(3)          :: RFace_D, B1r_D, B1t_D, TotalB_D
    real                        :: vAlfvenSi, wEnergyDensSi, wEnergyDensBc
    real                        :: ExpansionFactorInv
    logical                     :: IsClosedWSA = .false. 
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
    
    !\
    ! Update BCs for velocity and induction field: reflective BC
    !/
    VarsGhostFace_V(Ux_:Uz_) = -VarsTrueFace_V(Ux_:Uz_)
    VarsGhostFace_V(Bx_:Bz_) = B1t_D
  
    !\
    ! Update BCs for the mass density and pressure::
    !/
    iCell = iFace ; jCell = jFace ; kCell = kFace
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
       write(*,*)'ERROR: iSide = ', iSide
       call stop_mpi('incorrect iSide value in user_face_bcs')
    end select

    call get_plasma_parameters_cell(iCell, jCell, kCell, iBlockBc,&
         DensCell, PresCell)
    VarsGhostFace_V(Rho_) = &
         max(-VarsTrueFace_V(Rho_) + 2.0*(DensCell), &
         VarsTrueFace_V(Rho_))
    TBase = PresCell/DensCell
    VarsGhostFace_V(P_) = &
         max(VarsGhostFace_V(Rho_)*TBase, VarsTrueFace_V(P_))
    
    !\
    ! Update BCs for wave spectrum
    !/
    if(IsInitWave) then
       ! Check if this is closed or open field region
       call get_interpolated(ExpansionFactorInv_N, FaceCoords_D(x_),&
            FaceCoords_D(y_), FaceCoords_D(z_), ExpansionFactorInv)
       if(ExpansionFactorInv < cTolerance) IsClosedWSA = .true.
       
       ! set wave energy at inner boundary in open field region only
       if (IsClosedWSA) then  
          VarsGhostFace_V(WaveFirst_:WaveLast_) = 0.0
       else
          ! total wave energy depends on magnetic field magnitude
          TotalB_D = B0Face_D + VarsTrueFace_V(Bx_:Bz_)  
          TotalB = sum(TotalB_D**2)

          select case(TypeWaveInnerBc)
          case('WSA')
             
             vAlfvenSi = (TotalB/sqrt(VarsGhostFace_V(Rho_))) * No2Si_V(UnitU_)
             call get_total_wave_energy_dens(FaceCoords_D(x_), FaceCoords_D(y_),&
                  FaceCoords_D(z_), vAlfvenSi, wEnergyDensSi)
             wEnergyDensBc = wEnergyDensSi * Si2No_V(UnitP_)*WaveInnerBcFactor
         
          case('Turb')

             wEnergyDensBc = TotalB*WaveInnerBcFactor

          end select
       
          if(wEnergyDensBc < 0.0)then
             write(*,*) 'Negative TOTAL wave energy at inner BC'
             call stop_MPI('Error in user_face_bcs')
          end if
       
          ! Set BC for each frequency group
          call set_wave_state(wEnergyDensBc, VarsGhostFace_V, RFace_D, B0Face_D) 
       
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
  subroutine get_plasma_parameters_cell(iCell, jCell, kCell, iBlock,&
       DensCell, PresCell)
     
    ! This subroutine computes the cell values for density and pressure 
    ! assuming an isothermal atmosphere
    
    use ModGeometry,   ONLY: x_BLK, y_BLK, z_BLK, R_BLK
    use ModNumConst
    use ModPhysics,    ONLY: GBody, BodyRho_I, Si2No_V, UnitTemperature_
    use ModExpansionFactors,  ONLY: UMin, CoronalT0Dim
    implicit none

    integer, intent(in)  :: iCell, jCell, kCell, iBlock
    real, intent(out)    :: DensCell, PresCell
    real :: UFinal       !The solar wind speed at the far end of the Parker spiral,
                         !which originates from the given cell
    real :: URatio       !The coronal based values for temperature density 
                         !are scaled as functions of UFinal/UMin ratio
    real :: Temperature
    !--------------------------------------------------------------------------
    call get_bernoulli_integral(x_BLK(iCell, jCell, kCell, iBlock)/&
         R_BLK(iCell, jCell, kCell, iBlock),&
         y_BLK(iCell, jCell, kCell, iBlock)/R_BLK(iCell, jCell, kCell, iBlock),&
         z_BLK(iCell, jCell, kCell, iBlock)/R_BLK(iCell, jCell, kCell, iBlock), UFinal)
    URatio=UFinal/UMin

    !This is the temperature variation
    Temperature = CoronalT0Dim*Si2No_V(UnitTemperature_)/(min(URatio, 1.5))

    DensCell  = (1.0/URatio) &          !This is the density variation
         *BodyRho_I(1)*exp(-GBody/Temperature &
         *(1.0/max(R_BLK(iCell, jCell, kCell, iBlock), 0.90)-1.0))

    PresCell = DensCell*Temperature

  end subroutine get_plasma_parameters_cell
  !============================================================================
  subroutine user_set_ics

    ! user_set_ic : intilize MHD parameters according to an isothermal atmosphere solution.

    use ModMain,      ONLY: globalBLK
    use ModAdvance,   ONLY: State_VGB 
    use ModPhysics,   ONLY: inv_gm1, BodyTDim_I
    use ModGeometry
    implicit none

    integer :: i, j, k, iBLK
    logical :: oktest,oktest_me
    real    :: Dens, Pres
    real    :: x, y, z, R, ROne, Rmax, U0
    !--------------------------------------------------------------------------
    call set_oktest('user_set_ics', oktest, oktest_me)

    iBLK = globalBLK

    ! The sqrt is for backward compatibility with older versions of the Sc
    U0 = 4.0*sqrt(2.0E+6/BodyTDim_I(1))
    
    State_VGB(:,1:nI,1:nJ,1:nK,iBLK) = 1.0e-30 !Initialize the wave spectrum
    do k=1,nK; do j=1,nJ; do i=1,nI
       r = r_BLK(i,j,k,iBLK)
       State_VGB(Bx_:Bz_,i,j,k,iBLK) = 0.0
       call get_plasma_parameters_cell(i,j,k,iBLK,&
            Dens,Pres)
       State_VGB(rho_,i,j,k,iBLK) = Dens
       State_VGB(P_,i,j,k,iBLK)   = Pres
       State_VGB(RhoUx_,i,j,k,iBLK) = Dens *U0*x/R
       State_VGB(RhoUy_,i,j,k,iBLK) = Dens *U0*y/R
       State_VGB(RhoUz_,i,j,k,iBLK) = Dens *U0*z/R
       State_VGB(Ew_,i,j,k,iBLK) = nWave*1.0e-30

    end do; end do; end do

  end subroutine user_set_ics
  !============================================================================
  subroutine user_update_states(iStage,iBlock)

    ! user_update_states : calls update states_MHD, followed by a call to dissipate_waves, which
    !                      removes wave energy from the spectrum and adds it to the plasma pressure.
    ! ---------------------------------------------------------------------------------------------

    use ModAdvance, ONLY: State_VGB, B0_DGB

    use ModGeometry,ONLY: R_BLK
    use ModEnergy,  ONLY: calc_energy_cell
    use ModWaves,   ONLY: UseWavePressure, UseWavePressureLtd, &
                          AlfvenWavePlusFirst_,  AlfvenWavePlusLast_, &
                          AlfvenWaveMinusFirst_, AlfvenWaveMinusLast_ 

    integer,intent(in)           :: iStage, iBlock
    integer                      :: i, j, k
    character(len=*),parameter   :: NameSub='user_update_states'
    !--------------------------------------------
    !\
    ! Advect solution
    !/
    if(any(State_VGB(WaveFirst_:WaveLast_,1:nI,1:nJ,1:nK,iBlock)<0.0)) then
       write(*,*) NameSub,' : negative wave energy before MHD'
    end if
    
    call update_states_MHD(iStage, iBlock)
    
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
  contains
    subroutine dissipate_waves(iBlock)

      ! Implements frequency dependent wave damping mechanism.  This subroutine removes wave
      ! energy from the spectrum and adds it to the plasma pressure.
      ! Usage: The user can control the damping by setting logical flags via the input
      !        command #DAMPWAVES. Currently two mechanisms are available:
      !       DoDampCutoff (logical) - if true, Alfven waves will be totaly damped at
      !                                and above the local ion cyclotron frequency.
      !       DoDampSurface (logical) - if true, Alfven waves will be damped according
      !                                 to the local surface waves dissipation length.
      !                                 (implemented by R. Evans)

      use ModAdvance, ONLY: State_VGB
      use ModPhysics, ONLY: Gamma0
      use ModWaves,   ONLY: FrequencySi_W
      use ModNumConst,ONLY: cPi

      integer,intent(in)      :: iBlock
      integer                 :: i, j, k, iWave
      real                    :: dWavePres = 0.0, LogFreqRadian
      real                    :: LogFreqCutOff 
      character(len=*),parameter :: NameSub="dissipate_waves"
      ! -----------------------------------------------------------------
      if (DoDampCutOff) then
         do k=1,nK ; do j=1, nJ ; do i=1,nI
            call calc_cutoff_freq(i,j,k,iBlock,LogFreqCutOff)
            do iWave=WaveFirst_, WaveLast_
               LogFreqRadian = log(2*cPi*FrequencySi_W(iWave))
               if (LogFreqRadian .ge. LogFreqCutOff .and. &
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
    subroutine calc_cutoff_freq(i,j,k,iBLK , LogFreqCutOff)

      ! This subroutine calculates the cut-off frequency for Alfven waves, 
      ! which is the ion cyclotron frequency (in radians) 
      ! \omega_{c.o.}=zeB/m. For ions z=1. 
      
      use ModMain,                ONLY: x_, y_, z_
      use ModAdvance,             ONLY: State_VGB, B0_DGB
      use ModConst,               ONLY: cElectronCharge, cProtonMass,cTiny
      use ModPhysics,             ONLY: No2Si_V, UnitB_
      
      real, intent(out)           :: LogFreqCutOff
      integer, intent(in)         :: i, j, k, iBLK
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
    
  end subroutine user_update_states
  !========================================================================
  subroutine user_get_log_var(VarValue, TypeVar, Radius)
    
    ! user_get_log_var : allows the user to output spectral data to a dedicated log file.
    !                       In order to use this option, the user must modify the PARAM.in file:
    !                       - set the user flag UseUserLogVar to T                       
    !                       - include the command #SAVELOGFILE, with the following format:
    !                       #SAVELOGFILE
    !                       T
    !                       var
    !                       dn
    !                       dt
    !                       'string'
    !
    !                       where 'string' can be set to 'spectrogram' or 'cellspectrum' (refer to the the SWMF
    !                       user manual for further details on this command).
    !                       'spectrogram' : the full wave spectrum along a line is extracted
    !                       according to the parameters set in #SPECTROGRAM (see above). The data is
    !                       written to a file called Spectrum_n_xxx_peyyy.tec, where xxx stands
    !                       for simulation time/iteration number and yyy to the processor number. 
    !                       Files from different processors should be combined later.
    !                       
    !                       'cellspectrum' : the full spectrum in a single cell is written to a file
    !                       The cell is chosen by its x,y,z coordinates set in #CELLSPECTRUM (see above).
    !                       The output is written to a file Cell_Spectrum_n_xxx_peyyy.tec. Since only
    !                       one processor is involved, no post-processing is necessary.
    ! -----------------------------------------------------------------------------------------------

    use ModIO,                 ONLY: write_myname
    use ModVarIndexes,         ONLY: nWave

    real, intent(out)              :: VarValue 
    character (LEN=10), intent(in) :: TypeVar 
    real, intent(in), optional     :: Radius
    integer                        :: nWaveHalf
 
    !--------------------------------------------------------------------------
    ! This subroutine is used to output the wave spectrum to a file.
    ! No extravariables are defined for the log file here.
    VarValue = 0.0
    nWaveHalf = max(nWave/2,1)
   
    select case(TypeVar)
    case('spectrogram')
       ! write spectrum to file, extracted along a line (parallel to an axis).
       ! See description of #SPECTROGRAM command in the beginning of this module.
       call write_spectrogram
    case('cellspectrum')
       ! write spectrum to file, extracted from a single cell. Cell is chosen by #CELLSPECTRUM
       ! command (see description in the biginning of this module).
       call write_cell_spectrum
    case default
       VarValue = -7777.
       call write_myname;
       write(*,*) 'Warning in set_user_logvar: unknown logvarname = ',TypeVar
    end select

  contains
      subroutine write_cell_spectrum
    
        use ModProcMH
        use ModMain,     ONLY: iteration_number, unusedBLK, nBlockALL
        use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, dx_BLK, dy_BLK, dz_BLK
        use ModIoUnit,   ONLY: io_unit_new
        use ModAdvance,  ONLY: State_VGB
        use ModPhysics,  ONLY: No2Si_V, UnitX_, UnitP_
        use ModWaves,    ONLY: AlfvenWaveMinusFirst_,AlfvenWavePlusFirst_, &
             FrequencySi_W
    
        real,dimension(nWaveHalf,3)      :: Spectrum_II 
        real                             :: x, y, z, dxHalf, dyHalf,dzHalf
        real                             :: IwPlusSi,IwMinusSi
        integer                          :: iFreq, i, j, k, iBLK
        integer                          :: iUnit, iError, aError
        logical                          :: DoSaveCellSpectrum = .false.
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
              dxHalf = dx_BLK(iBLK)/2.0
              dyHalf = dy_BLK(iBLK)/2.0
              dzHalf = dz_BLK(iBLK)/2.0
              
              if((xTestSpectrum >= x - dxHalf ) .and. (xTestSpectrum <= x + dxHalf) .and. &
                 (yTestSpectrum >= y - dyHalf ) .and. (yTestSpectrum <= y + dyHalf) .and. &
                 (zTestSpectrum >= z - dzHalf ) .and. (zTestSpectrum <= z + dzHalf) ) then

                 DoSaveCellSpectrum = .true.
                 do iFreq=1,nWaveHalf
                    IwPlusSi  = No2Si_V(UnitP_)* &
                         State_VGB(AlfvenWavePlusFirst_+iFreq-1,i,j,k,iBLK)
                    IwMinusSi = No2Si_V(UnitP_)* &
                         State_VGB(AlfvenWaveMinusFirst_+iFreq-1,i,j,k,iBLK)
                    Spectrum_II(iFreq,1) = log(FrequencySi_W(iFreq))
                    Spectrum_II(iFreq,2) = IwPlusSi
                    Spectrum_II(iFreq,3) = IwMinusSi
                    
                 end do
              end if
           end do; end do ; end do
        end do
    
        if (DoSaveCellSpectrum) then
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
 
      subroutine write_spectrogram
    
        use ModProcMH
        use ModMain,      ONLY: iteration_number, unusedBLK, nBlockALL
        use ModGeometry,  ONLY: x_BLK, y_BLK, z_BLK, dz_BLK, dx_BLK
        use ModIoUnit,    ONLY: io_unit_new
        use ModAdvance,   ONLY: State_VGB
        use ModPhysics,   ONLY: No2Si_V, UnitX_, UnitP_
        use ModWaves,     ONLY: AlfvenWaveMinusFirst_, AlfvenWavePlusFirst_, FrequencySi_W
        
        implicit none
        
        real, allocatable,dimension(:,:) :: Cut_II ! Array to store log variables
        integer                          :: nCell, nRow, iRow 
        real                             :: dx, dz, x, y, z, IwPlusSi, IwMinusSi
        integer                          :: iFreq,i,j,k,iBLK
        integer                          :: iUnit, iError, aError
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
                          Cut_II(iRow,2) = log(FrequencySi_W(iFreq))
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

  end subroutine user_get_log_var
!===========================================================================
  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    ! user_set_plot_var : implement plot varaibles related to Alfven waves.
    !                     'wpres' - total wave pressure (summed over spectrum).
    !                     'poynt' - Poynting flux of Alfven waves neglecting plasma speed.
    !                     'poyntur' - Poynting flux of Alfven waves.
    ! ----------------------------------------------------------------------------------

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
    real                         :: PoyntFlux
    integer                      :: i, j, k
    logical                      :: IsError
    !-------------------------------------------------------------------    
    !UsePlotVarBody = .true. 
    !PlotVarBody = 0.0 
    IsFound=.true.

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

  contains
    subroutine calc_poynt_flux(i,j,k,iBLK, UseUr, PoyntFluxSi)

      Use ModMain,       ONLY: x_, y_, z_
      use ModAdvance,    ONLY: State_VGB, B0_DGB
      use ModGeometry,   ONLY: x_BLK, y_BLK, z_BLK, R_BLK
      use ModPhysics,    ONLY: No2Si_V, UnitB_, UnitRho_, UnitU_, UnitP_
      use ModConst,      ONLY: cMu
      implicit none
      
      integer,intent(in)         :: i, j, k, iBLK
      logical,intent(in)         :: UseUr
      real,intent(out)           :: PoyntFluxSi
      real,dimension(3)          :: r_D, BSi_D, USi_D
      real                       :: x, y, z, BrSi, UrSi, RhoSi  
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
 
  end subroutine user_set_plot_var
  !=========================================================================== 
  subroutine user_specify_refinement(iBlock, iArea, DoRefine)

    ! user_specify_refinement : improve refinement in the current sheet region.
    ! ------------------------------------------------------------------------

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

    ! user_set_boundary_cells :required when "extra" boundary conditions are used.
    ! --------------------------------------------------------------------------

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

