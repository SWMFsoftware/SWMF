!instrumental broadening
!check boxes

program spectrum

  use ModConst, ONLY           : cLightSpeed, cBoltzmann, cProtonMass, cAU, &
       cPi, rSun

  implicit none

  ! Logical variables of running modes
  logical                     :: IsVerbose =  .false., IsDebug = .false.
  logical                     :: IsDataFile = .false. ! Use input file or not

  integer                     :: iError 

  !One line option
  logical                     :: IsOneLine = .false.
  real                        :: OneLineWavelength

  ! Variables for output file
  character(len=200)          :: NameSpectrumFile = 'spectrum.out'
  character(len=200)          :: NameLabelFile = 'label.out'

  ! Variables for input files
  character(len=200)          :: StringLine 
  character(len=200)          :: NameDataFile, NameTableFile
  character(len=200)          :: TypeDataFile

  ! Variables for instrument (if any)
  logical                     :: IsInstrument = .false.
  character(len=200)          :: NameInstrument
  integer                     :: nPixel ! Number of pixels along slit
  integer                     :: nWavelengthBin, iWavelengthBin
  real,allocatable            :: DLambdaInstr_I(:) ! Instrumental broadening
  real                        :: SizeWavelengthBin ! Resolution in wavelength

  ! Variables for solar wind input data file (if any)
  logical                     :: IsDataBlock = .false. ! Overwrite data size
  logical                     :: IsOnePixel = .false. ! Select one pixel only
  logical                     :: IsPe = .false.
  logical                     :: IsPpar = .false., IsPperp = .false.
  logical                     :: IsNoAlfven = .false. ! Ignore Alfven waves
  logical                     :: IsDoppler = .true. ! Calculate Doppler shift
  integer                     :: n1Block, n2Block, n3Block ! Define data size
  integer                     :: iOnePixel, jOnePixel, kOnePixel
  ! Rotation angles for box variable rotation
  integer                     :: iDimLOS=1, iDimVertical=2, iDimHorizontal=3

  ! Variables to read solar wind data file in
  integer                     :: nVar   ! Number of variables   
  integer                     :: nDim   ! Number of dimensions  
  integer                     :: nParam ! Number of parameters  
  integer                     :: n1, n2, n3 ! Data box size
  real                        :: CoordMin_D(3), CoordMax_D(3)        
  real,allocatable            :: Var_VIII(:,:,:,:)
  ! For H:He 10:1 fully ionized plasma the proton:electron ratio is 1/(1+2*0.1)
  real                        :: ProtonElectronRatio = 0.83
  real                        :: dx = 1., dy = 1., dz = 1., dVperd2 = 1.
  
  ! Variables for uniform data
  logical                     :: IsUniData = .false. ! Overwrite data w/ const
  real                        :: RhoUni, UxUni, UyUni, UzUni, BxUni, ByUni
  real                        :: BzUni
  real                        :: TparUni, TperpUni, TeUni, I01Uni, I02Uni

  ! Indexes for the solar wind variables
  integer, parameter          :: &
       rho_  =  1, &              
       ux_   =  2, &                  
       uy_   =  3, &               
       uz_   =  4, & 
       bx_   =  5, &
       by_   =  6, &
       bz_   =  7, &
       tpar_ =  8, &
       tperp_ =  9, &
       te_   = 10, &
       t_    = 11, &
       I01_  = 12, &
       I02_  = 13

  integer                     :: nLocalVar = 13 ! Number of used variables

  ! Variables for the wavelengths of interest
  logical                     :: IsAllLines = .false. ! Ignore unobserved lines
  integer                     :: iWavelengthInterval, nWavelengthInterval
  integer                     :: nMaxLine, nLineFound
  real, allocatable           :: WavelengthInterval_II(:,:) 
  ! For Doppler-shift wavelength intervals are shifted by 10% of speed of light
  real, allocatable           :: WavelengthIntervalShifted_II(:,:) 

  ! Temperature and density grid for contrbution function (G)
  real                        :: DLogN, DLogT
  real                        :: MinLogN, MinLogT, MaxLogN, MaxLogT

  ! Maximum index range in density and temperature
  integer                     :: MaxI, MinI, MaxJ, MinJ

  ! Variables for lines
  integer                     :: nBin ! Wavelength bin number  
  integer                     :: i, iLine, jPixel, kPixel 
  integer                     :: nLineAll ! All lines of interest

  ! Derived type to read tabulated G values
  type LineTableType
     character(len=6)         :: NameIon
     real                     :: Aion
     integer                  :: nLevelFrom, nLevelTo ! Levels of transition
     ! Indexes on density-temperature grid 
     ! i = density index, j = temperature index 
     integer                  :: iMin, jMin, iMax, jMax  
     real                     :: LineWavelength 
     real, allocatable        :: g_II(:,:)
     real                     :: StartLogT
  end type LineTableType

  type(LineTableType), allocatable :: LineTable_I(:)

  ! Derived type of output 
  type SpectrumTableType
     integer                  :: nBin     
     real,allocatable         :: Spectrum_III(:,:,:), SpectrumGrid_I(:)
  end type SpectrumTableType

  type(SpectrumTableType), allocatable :: SpectrumTable_I(:)

  ! Derived type for the Labels
  type LabelTableType
     character(len=6)         :: NameIon
     real                     :: Lambda
     real                     :: FluxMono
  end type LabelTableType

  type(LabelTableType), allocatable :: LabelTable_III(:,:,:)

  character(len=*), parameter :: NameSub = 'spectrum.f90'

  !---------------------------------------------------------------------------
  write(*,*)'Spectrum.exe starting'

!!!  call MPI_init(iError)

  call read_param

  if(IsDataFile)then
     call read_data
  else
     call set_data_block
  endif

  call read_table

  nLineAll = min(nMaxLine,nLineFound)

  ! Allocate table for labels
  allocate(LabelTable_III(nLineAll,n2,n3))
  if(IsVerbose)write(*,*)'allocated LabelTable_III'

  ! Loop over Chianti lines
  if(IsVerbose)write(*,*)'IsOneLine = ', IsOneLine
  do iLine = 1,nLineAll
     if (.not. IsOneLine)call calc_flux
     if (IsOneLine .and.  &
          OneLineWavelength == LineTable_I(iLine)%LineWavelength)call calc_flux
  end do


  call save_all
  if(IsVerbose)write(*,*)'done with save_all'
  call save_label
  if(IsVerbose)write(*,*)'done with save_label'
  write(*,*)'Spectrum.exe ending'

!!!  call MPI_finalize(iError)

contains

  !==========================================================================
  subroutine set_data_block
    ! When no data file input is used set up uniform data values in defined box
    integer, parameter             :: iUnitOut = 18
    real                           :: MinWavelength, Maxwavelength
    character(len=*), parameter    :: NameSub = 'set_data_block'
    !------------------------------------------------------------------------

    n1 = n1Block
    n2 = n2Block
    n3 = n3Block
    allocate(Var_VIII(nLocalVar,n1,n2,n3))

    if (IsUniData) then
       ! cm^-3 --> kg/m^3
       Var_VIII(rho_,1:n1,1:n2,1:n3)  = RhoUni * 1e6 * cProtonMass
       ! km/s --> m/s
       Var_VIII(ux_,1:n1,1:n2,1:n3)   = UxUni * 1e3
       Var_VIII(uy_,1:n1,1:n2,1:n3)   = UyUni * 1e3
       Var_VIII(uz_,1:n1,1:n2,1:n3)   = UzUni * 1e3
       ! G --> T
       Var_VIII(bx_,1:n1,1:n2,1:n3)   = BxUni * 1e-4
       Var_VIII(by_,1:n1,1:n2,1:n3)   = ByUni * 1e-4
       Var_VIII(bz_,1:n1,1:n2,1:n3)   = BzUni * 1e-4
       ! K
       Var_VIII(tpar_,1:n1,1:n2,1:n3) = TparUni
       Var_VIII(tperp_,1:n1,1:n2,1:n3) = TperpUni
       Var_VIII(te_,1:n1,1:n2,1:n3)   = TeUni
       ! erg cm^-3 --> J m^-3
       Var_VIII(I01_,1:n1,1:n2,1:n3)  = I01Uni * 1e-1
       Var_VIII(I02_,1:n1,1:n2,1:n3)  = I02Uni * 1e-1

       dVperd2 = 1.0
    else 
       call CON_stop( &
            NameSub//' need data input!!! ')
    end if

    ! Set up bins for Spectrum
    allocate(SpectrumTable_I(nWavelengthinterval))
    nBin = 0
    do iWavelengthInterval=1, nWavelengthInterval
       MinWavelength = WavelengthInterval_II(1,iWavelengthInterval)
       MaxWavelength = WavelengthInterval_II(2,iWavelengthInterval)
       nWavelengthBin = nint((MaxWavelength-MinWavelength)/SizeWavelengthBin)
       nBin = nWavelengthBin
       allocate(SpectrumTable_I(iWavelengthinterval)%Spectrum_III(n2,n3, &
            nWavelengthBin))
       allocate(SpectrumTable_I(iWavelengthinterval)%SpectrumGrid_I(&
            nWavelengthBin+1))
       SpectrumTable_I(iWavelengthinterval)%Spectrum_III(:,:,:)=0.0
       SpectrumTable_I(iWavelengthinterval)%nBin=nBin
       do iWavelengthBin=1,nWavelengthBin+1
          ! Values of each wavelength bin correspond to the center of the bin
          SpectrumTable_I(iWavelengthinterval)%SpectrumGrid_I(iWavelengthBin)&
               = MinWavelength + (iWavelengthBin-.5)*SizeWavelengthBin
       end do
    end do
  end subroutine set_data_block
  !==========================================================================
  subroutine save_label
    integer, parameter             :: iUnitOut = 18

    character(len=*), parameter    :: NameSub = 'save_label'
    !------------------------------------------------------------------------

    open(unit=iUnitOut,file=NameLabelFile,action='write')
    write(iUnitOut,*)"SPECTRUM output file - labels of lines"
    write(iUnitOut,*)"found spectral lines: nLineAll = ",nLineAll
    write(iUnitOut,*)"pixels in grid: n2 X n3 = ",n2 ," X ", n3
    write(iUnitOut,*)"jPixel kPixel NameIon Lambda FluxMono"
    write(iUnitOut,*)"Number of data lines in this file = ",nLineAll*n2*n3

    do jPixel = 1,n2
       do kPixel = 1,n3
          do iLine = 1,nLineAll
             write(iUnitOut,*)jPixel,kPixel, &
                  LabelTable_III(iLine,jPixel,kPixel)%NameIon, &
                  LabelTable_III(iLine,jPixel,kPixel)%Lambda, &
                  LabelTable_III(iLine,jPixel,kPixel)%FluxMono
          end do
       end do
    end do

    close(unit=iUnitOut)

  end subroutine save_label
  !==========================================================================
  subroutine save_all

    use ModPlotFile, ONLY          : save_plot_file

    character                      :: Namefile
    !------------------------------------------------------------------------

    ! Variables for output file
    character(len=5)     :: TypeFileSpectrum = 'ascii'
    character(len=200)   :: StringHeaderSpectrum = &
         '[A] [erg sr^-1 cm^-2 A^-1]'
    real,allocatable            :: &
         Intensity_VIII(:,:,:,:), CoordWave_I(:), CoordPixel_DII(:,:,:)

    integer:: nWave, iWaveInterval, nWaveInterval, iWaveBin, iWave, & 
         nWaveBin, nWaveAll
    real                        :: DWaveBin, WavelengthMin, WavelengthMax
    character(len=*), parameter :: NameSub = 'save_all'
    !------------------------------------------------------------------------

    nWaveAll = 0
    nWaveInterval = nWavelengthInterval
    do iWavelengthInterval = 1,nWavelengthInterval
       nWaveAll = nWaveAll + SpectrumTable_I(iWavelengthInterval)%nBin
    end do

    allocate( &
         Intensity_VIII(1,nWaveAll,n2,n3), &
         CoordWave_I(nWaveAll), CoordPixel_DII(2,n2,n3))

    ! Number of wave length processed so far
    nWave = 0
    ! Loop over intervals
    if(IsVerbose) write(*,*)"  nWaveInterval = ",nWaveInterval
    do iWaveInterval = 1, nWaveInterval

       ! Number of wave length in this interval
       nWaveBin = SpectrumTable_I(iWaveInterval)%nBin
       WavelengthMin = WavelengthInterval_II(1,iWaveInterval)
       WavelengthMax = WavelengthInterval_II(2,iWaveInterval)
       DWaveBin = (WavelengthMax - WavelengthMin)/nWaveBin

       ! Loop over wave lengths in this interval
       do iWaveBin = 1, nWaveBin
          ! Global wave length index
          iWave = nWave + iWaveBin

          ! Save intensity and coordinate into global array
          Intensity_VIII(1,iWave,:,:) = &
               SpectrumTable_I(iWaveInterval)%Spectrum_III(:,:,iWaveBin)
          CoordWave_I(iWave) = WavelengthMin + (iWaveBin - 0.5)*DWaveBin
       end do
       ! Finished processing this interval
       nWave = nWave + nWaveBin
    end do

    ! nDim = 1 if n2 = n3 = 1, nDim = 2 if n3 = 1 /= n2 otherwise nDim= 3
    if(n3==1)then
       if(n2==1)then
          call save_plot_file(NameFile = NameSpectrumFile, &
               TypeFileIn     = TypeFileSpectrum,      &
               StringHeaderIn = StringHeaderSpectrum,  &
               NameVarIn      = "wavelength flux",     &
               nDimIn         = 1,                     &
               Coord1In_I     = CoordWave_I,           &
               VarIn_VIII      = Intensity_VIII)
       else
          call save_plot_file(NameFile = NameSpectrumFile, &
               TypeFileIn     = TypeFileSpectrum,      &
               StringHeaderIn = StringHeaderSpectrum,  &
               NameVarIn      = "wavelength y flux",   &
               nDimIn         = 2,                     &
               CoordMinIn_D   = CoordMin_D,            &
               CoordMaxIn_D   = CoordMax_D,            &
               Coord1In_I     = CoordWave_I,           &
               VarIn_VIII     = Intensity_VIII)
       endif
    else
       call save_plot_file(NameFile = NameSpectrumFile, &
            TypeFileIn     = TypeFileSpectrum,      &
            StringHeaderIn = StringHeaderSpectrum,  &
            NameVarIn      = "wavelength x y flux", &
            nDimIn         = 3,                     &
            CoordMinIn_D   = CoordMin_D,            &
            CoordMaxIn_D   = CoordMax_D,            &
            Coord1In_I     = CoordWave_I,           &
            VarIn_VIII      = Intensity_VIII)
    endif

    deallocate(Intensity_VIII, CoordWave_I, CoordPixel_DII)

  end subroutine save_all

  !==========================================================================
  subroutine calc_flux

    use ModInterpolate, ONLY: bilinear

    integer                        :: i, iInterval, iBin, iCenter
    integer                        :: iNMin, jTMin, iNMax, jTMax

    real                           :: FluxMono
    real                           :: Lambda, LambdaSI, Lambda0SI, DLambdaSI
    real                           :: DLambda, DLambdaSI2, Lambda_shifted
    real                           :: DLambdaInstr2 = 0.0 !!! do it later
    real                           :: Zplus2, Zminus2, CosAlpha, SinAlpha
    real                           :: B_D(3), Bnorm_D(3)
    real                           :: Unth2, Uth2
    real                           :: Gint, LogNe, LogTe, Rho
    real, allocatable              :: Glambda_II(:,:)
    real                           :: Tlos
    real                           :: Aion
    real                           :: TShift
    character(len=*), parameter    :: NameSub='calc_flux'
    !------------------------------------------------------------------------

    allocate(Glambda_II(MinI:MaxI,MinJ:MaxJ))

    Lambda   = LineTable_I(iLine)%LineWavelength
    Aion     = LineTable_I(iLine)%Aion
    ! Fill in LabelTabel%NameIon
    LabelTable_III(iLine,:,:)%NameIon = LineTable_I(iLine)%NameIon

    if(IsOnePixel)write(*,*)'n1 n2 n3 i j k =',n1,n2,n3,i,jPixel,kPixel

    do kPixel=1,n3
       do jPixel=1,n2
          do i=1,n1

             if(IsOnePixel .and. (kPixel/=kOnePixel .or. &
                  jPixel/=jOnePixel .or. i/=iOnePixel))CYCLE

             ! Cells inside body or behind the solar disk
             if(all(Var_VIII(1:7,i,jPixel,kPixel)==0))CYCLE

             ! Doppler shift while x axis is oriented towards observer
             if(IsDoppler)then
                Lambda_shifted = (-Var_VIII(ux_,i,jPixel,kPixel)/cLightSpeed &
                     + 1 ) * Lambda
                Lambda = Lambda_shifted
             endif

             ! Find which interval wavelength belongs to
             iInterval = -1
             do iWavelengthInterval = 1, nWavelengthInterval
                if((Lambda <= WavelengthInterval_II(2,iWavelengthInterval)) &
                     .and.&
                     (Lambda >= WavelengthInterval_II(1,iWavelengthInterval)))&
                     then
                   iInterval = iWavelengthInterval
                   EXIT
                end if
             end do
             if(iInterval == -1)then
                write(*,*)"Lambda = ",Lambda," will be left out!"
                write(*,*)"Intervals begin and end = ", &
                     WavelengthInterval_II(:,:)
                EXIT
             endif

             iBin =0
             do
                iBin=iBin+1
                if ((Lambda > &
                     SpectrumTable_I(iInterval)%SpectrumGrid_I(iBin)) &
                     .and.&
                     (Lambda < & 
                     SpectrumTable_I(iInterval)%SpectrumGrid_I(iBin+1)))&
                     then
                   iCenter=iBin
                   EXIT
                end if
                if (SpectrumTable_I(iInterval)%SpectrumGrid_I(iBin+1) > &
                     WavelengthInterval_II(2,iInterval))EXIT
                if (iBin==nBin)then
                   write(*,*)"Lambda = ",Lambda," will be left out!"
                   EXIT
                endif
             end do

             ! Calculate Elzasser variables
             Rho = Var_VIII(rho_,i,jPixel,kPixel)
             Zplus2   = Var_VIII(I01_,i,jPixel,kPixel) * 4.0 / Rho
             Zminus2  = Var_VIII(I02_,i,jPixel,kPixel) * 4.0 / Rho

             ! Calculate angle between LOS and B directions
             B_D      = Var_VIII(bx_:bz_,i,jPixel,kPixel)
             CosAlpha = abs(B_D(1)/sqrt(max(sum(B_D**2), 1e-30)))

             ! Calculate temperature relative to the LOS direction
             SinAlpha = sqrt(1 - CosAlpha**2)
             Tlos = SinAlpha**2 * Var_VIII(tperp_,i,jPixel,kPixel) &
                  + CosAlpha**2 * Var_VIII(tpar_,i,jPixel,kPixel)

             ! Calculate thermal and non-thermal broadening
             Unth2    = 1.0/16.0 * (Zplus2 + Zminus2) * SinAlpha**2
             Uth2     = cBoltzmann * Tlos/(cProtonMass * Aion)  

             ! Convert from kg m^-3 to kg cm^-3 (*1e-6)
             ! and divide by cProtonMass in kg so Ne is in cm^-3
             ! 1 : 0.83  electron to proton ratio is assumed
             LogNe = log10(Rho*1e-6/cProtonMass/ProtonElectronRatio)
             LogTe = log10(Var_VIII(te_,i,jPixel,kPixel))

             ! Convert to SI
             LambdaSI = LineTable_I(iLine)%LineWavelength * 1e-10 

             ! Add thermal and non-thermal broadening
             DLambdaSI2 = LambdaSI**2 * (Uth2 + Unth2)/cLightSpeed**2

             ! Add instrumental broadening (if any)
             if(IsInstrument)DLambdaSI2 = DLambdaSI2 + DLambdaInstr2

             ! Convert [m] --> [A]
             DLambdaSI = sqrt(DLambdaSI2)
             DLambda   = DLambdaSI * 1e10

             ! Get the contribution function
             iNMin  = LineTable_I(iLine)%iMin
             jTMin  = LineTable_I(iLine)%jMin
             iNMax  = LineTable_I(iLine)%iMax
             jTMax  = LineTable_I(iLine)%jMax
             Glambda_II = LineTable_I(iLine)%g_II(:,:)

             Gint = bilinear(Glambda_II, iNMin, iNMax, jTMin, jTMax, &
                  (/ LogNe/DLogN , LogTe/DLogT /),DoExtrapolate=.true.)

             ! When Gint becomes negative due to extrapolation -> move to next
             if(Gint<=0)CYCLE 

             ! Calculate flux and spread it on the Spectrum_II grids
             ! Intensity calculation according to Aschwanden p.58 Eq(2.8.4)
             FluxMono = Gint * (10.0**LogNe)**2 / (4*cPi) * dVperd2

             if(IsDebug)then
                write(*,*)'                                                   '
                write(*,*)'      ',LineTable_I(iLine)%NameIon,'         '
                write(*,*)'LineWavelength = ',LineTable_I(iLine)%LineWavelength
                write(*,*)'Gint = ',log10(Gint),' FluxMono = ',FluxMono
                write(*,*)'Aion = ',Aion
                write(*,*)'LambdaSI =',LambdaSI
                write(*,*)'DLambdaSI, DLambda =', DLambdaSI, DLambda
                write(*,*)'Uth, Unth = ',sqrt(Uth2), sqrt(Unth2)
                write(*,*)'LogTlos = ',log10(Tlos)
                write(*,*)'iLine = ',iLine
                write(*,*)'iCenter , iInterval = ',iCenter, iInterval
                write(*,*)'jTMin, jTMax = ', jTMin, jTMax
                write(*,*)'iNMin, iNMax = ', iNMin, iNMax
                write(*,*)'LogNe,DLogN , LogTe,DLogT = ', &
                     LogNe,DLogN , LogTe,DLogT
                write(*,*)'LogNe/DLogN , LogTe/DLogT = ', &
                     LogNe/DLogN , LogTe/DLogT
                write(*,*)'Gint = ',Gint
                write(*,*)'                                                   '
             endif

             ! Fill in LabelTabel%Lambda,FluxMono
             LabelTable_III(iLine,:,:)%Lambda = Lambda
             LabelTable_III(iLine,jPixel,kPixel)%FluxMono = &
                  FluxMono /(sqrt(2*cPi) * DLambda)

             call disperse_line(iInterval, iCenter, Lambda, DLambda, FluxMono)
          end do
       end do
    end do

  end subroutine calc_flux

  !==========================================================================
  subroutine disperse_line(iInterval,iCenter,Lambda,DLambda,FluxMono)

    real, intent(in)            :: Lambda, DLambda, FluxMono
    integer, intent(in)         :: iInterval, iCenter

    integer                     :: iStep, iWave, nWaveBin    
    integer                     :: iBin, iBegin, iEnd
    real                        :: Flux, Phi, InvNorm, InvSigma2
    real                        :: LambdaBin, LambdaBegin, LambdaEnd
    real                        :: LambdaDist

    character(len=*), parameter :: NameSub='disperse_line'
    !------------------------------------------------------------------------

    nWaveBin = SpectrumTable_I(iInterval)%nBin

    ! Beginning and end of gaussian truncated to +/-5 sigma in [A]
    LambdaBegin = (Lambda - 5*DLambda)
    LambdaEnd   = (Lambda + 5*DLambda)

    ! Find the corresponding wavelength bin for starting wavelength
    do iBegin = 1, nWaveBin - 1
       if (LambdaBegin < SpectrumTable_I(iInterval)%SpectrumGrid_I(iBegin+1))&
            EXIT
    end do

    ! Start at iBegin for efficiency and
    ! stop at nWaveBin-1 so iEnd = mWaveBin if no EXIT was performed
    do iEnd = iBegin, nWaveBin-1
       if (LambdaEnd < SpectrumTable_I(iInterval)%SpectrumGrid_I(iEnd)) EXIT
    end do

    InvNorm   = 1/(sqrt(2*cPi) * DLambda)
    InvSigma2 = 1/(2*DLambda**2) 

    ! Update bins between begin and end indices by adding the Gaussian 
    ! distribution
    do iBin = iBegin , iEnd
       ! Get wavelength from the center of the bin
       LambdaBin = SpectrumTable_I(iInterval)%SpectrumGrid_I(iBin)

       ! Get distance from peak wavelength in SI
       LambdaDist = Lambda - LambdaBin 

       ! Calculate Gaussian of the line
       Phi = InvNorm * exp(-LambdaDist**2 * InvSigma2)

       ! Calculate total monochromatic flux 
       Flux = FluxMono*Phi

       ! Update bin with flux
       SpectrumTable_I(iInterval)%Spectrum_III(jPixel,kPixel,iBin) = &
            SpectrumTable_I(iInterval)%Spectrum_III(jPixel,kPixel,iBin) + Flux
    end do

  end subroutine disperse_line

  !==========================================================================
  subroutine read_param

    use ModMPI, ONLY: MPI_COMM_SELF
    use ModReadParam

    character(len=200)      :: NameCommand
    logical                     :: IsNoInstrument = .false., DoEcho = .false.
    integer                     :: iPixel
    character(len=*), parameter :: NameSub='read_param'
    !------------------------------------------------------------------------
    call read_file('SPECTRUM.in', iCommIn=MPI_COMM_SELF)

    call read_init('  ') 

    ! Read SPECTRUM.in

    READPARAM: do
       if(.not.read_line(StringLine) ) EXIT READPARAM

       if(.not.read_command(NameCommand)) CYCLE READPARAM

       select case(NameCommand)
       case("#ECHO")
          call read_var('DoEcho',DoEcho)
          call read_echo_set(DoEcho)

       case("#VERBOSE")
          call read_var('IsVerbose', IsVerbose)

       case("#DEBUG")
          call read_var('IsDebug', IsDebug)

       case("#OUTFILE")
          call read_var('NameSpectrumFile',NameSpectrumFile)
          call read_var('NameLabelFile',NameLabelFile)

       case("#DATAFILE")
          IsDataFile = .true.
          call read_var('NameDataFile',NameDataFile)
          call read_var('TypeDataFile',TypeDataFile)

       case("#TABLEFILE")
          call read_var('NameTableFile',NameTableFile)
          call read_var('nMaxLine',nMaxLine)

       case("#WAVELENGTHINTERVAL")
          IsNoInstrument = .true.

          call read_var('nWavelengthInterval',nWavelengthInterval)
          if(IsInstrument)then
             deallocate(WavelengthInterval_II)
             write(*,*)'INSTRUMENT intervals changed to WAVELENGTHINTERVALS'
          endif
          allocate(WavelengthInterval_II(2,nWavelengthInterval))
          do iWavelengthInterval=1, nWavelengthInterval
             call read_var('IntervalMin', &
                  WavelengthInterval_II(1,iWavelengthInterval))
             call read_var('IntervalMax', &
                  WavelengthInterval_II(2,iWavelengthInterval))
          end do
          call read_var('SizeWavelengthBin',SizeWavelengthBin)

       case("#DATABLOCK")
          IsDataBlock = .true.
          call read_var('n1',n1Block)
          call read_var('n2',n2Block)
          call read_var('n3',n3Block)

       case("#INSTRUMENT")
          IsInstrument = .true.
          iPixel = 0
          nPixel = 0
          call read_var('NameInstrument',NameInstrument)

          select case(NameInstrument)
          case("EIS")
             nWavelengthInterval = 2
             nPixel = 512
             allocate(DLambdaInstr_I(nPixel))

             do iPixel=1,nPixel
                DLambdaInstr_I(iPixel) = 0.0
             end do

             SizeWavelengthBin = 0.0223
             if(IsNoInstrument)then
                write(*,*)'INSTRUMENT interval changed to WAVELENGTHINTERVALS'
             else
                allocate(WavelengthInterval_II(2,nWavelengthInterval))
                WavelengthInterval_II(:,1) = (/ 170 ,210 /)
                WavelengthInterval_II(:,2) = (/ 250 ,290 /)
             endif
          case default
             write(*,*) NameSub // ' WARNING: unknown #INSTRUMENT '
          end select

       case("#UNIFORMDATA")
          IsUniData = .true.
          call read_var('IsPpar',IsPpar)
          call read_var('IsPe',IsPe)
          call read_var('RhoUni',RhoUni)
          call read_var('UxUni',UxUni)
          call read_var('UyUni',UyUni)
          call read_var('UzUni',UzUni)
          call read_var('BxUni',BxUni)
          call read_var('ByUni',ByUni)
          call read_var('BzUni',BzUni)
          ! One temperature
          if(.not. IsPe .and. .not. IsPpar)then
             call read_var('TparUni',TparUni)
             TperpUni = TparUni
             TeUni = TparUni
             ! Proton +  electron temperatures
          elseif(IsPe .and. .not. IsPpar)then
             call read_var('TparUni',TparUni)
             TperpUni = TparUni
             call read_var('TeUni',TeUni)
             ! Electron + anisotropic proton temperatures
          elseif(IsPe .and. IsPpar)then
             IsPperp = .true.
             call read_var('TparUni',TparUni)
             call read_var('TperpUni',TperpUni)
             call read_var('TeUni',TeUni)
          else 
             write(*,*) NameSub // ' WARNING: check temperature setting' // &
                  'Select either 1-, 2-(proton+electron), or 3 temperature' // & 
                  '(electron+anisotropic proton) model!' 
          endif
          call read_var('I01Uni',I01Uni)
          call read_var('I02Uni',I02Uni)

       case("#NOALFVEN")
          call read_var('IsNoAlfven',IsNoAlfven)

       case("#ISONEPIXEL")
          call read_var('IsOnePixel',IsOnePixel)
          call read_var('iOnePixel',iOnePixel)
          call read_var('jOnePixel',jOnePixel)
          call read_var('kOnePixel',kOnePixel)

       case("#ISALLLINES")
          call read_var('IsAllLines',IsAllLines) ! Unobserved lines included

       case("#ISDOPPLER")
          call read_var('IsDoppler',IsDoppler)

       case("#ONELINE")
          IsOneLine = .true.
          call read_var('OneLineWavelength',OneLineWavelength)

       case default
          write(*,*) NameSub // ' WARNING: unknown #COMMAND '

       end select
    end do READPARAM

    ! In case of Doppler shift applied some lines might appear that would not
    ! otherwise. For this reason we introduce the extended WaveLengthIntervals
    ! for the sake of searching of lines of interest, assuming plasma moves no
    ! faster than 10% of the speed of light.
    if(IsDoppler)then
       allocate(WavelengthIntervalShifted_II(2,nWavelengthInterval))
       WavelengthIntervalShifted_II(1,:) =  WavelengthInterval_II(1,:) * 0.9
       WavelengthIntervalShifted_II(2,:) =  WavelengthInterval_II(2,:) * 1.1
    endif

  end subroutine read_param

  !==========================================================================
  subroutine read_data

    use ModPlotFile,         ONLY: read_plot_file
    use ModUtilities,        ONLY: split_string, lower_case
    use ModCoordTransform,   ONLY: rot_matrix_x, rot_matrix_y, rot_matrix_z

    character(len=200)          :: StringHeader
    character(len=200)          :: NameVar

    integer, parameter          :: MaxNameVar = 100
    character(len=20)           :: NameVar_V(MaxNameVar)

    integer                     :: nStep, i, j, k
    integer                     :: k1, k2
    integer                     :: nVarName, iVar

    real                        :: MinWavelength, MaxWavelength
    real                        :: Xangle, Yangle, Zangle, Rot_DD(3,3)
    real                        :: Param_I(3)
    real,allocatable            :: VarIn_VIII(:,:,:,:)
    real                        :: Coord, Dz1, Dz2

    character(len=*), parameter :: NameSub = 'read_data'
    !------------------------------------------------------------------------
    call read_plot_file(NameFile=NameDataFile, &
         TypeFileIn=TypeDataFile,              &
         StringHeaderOut=StringHeader,         &
         nStepOut=nStep,                       &
         nDimOut=nDim,                         &
         nParamOut=nParam,                     &
         ParamOut_I = Param_I,                 &
         nVarOut=nVar,                         &
         n1Out=n1,                             &
         n2Out=n2,                             &
         n3Out=n3,                             &
         NameVarOut=NameVar,                   &
         iErrorOut = iError)

    if(iError /= 0) call CON_stop( &
         NameSub//' could not header from '//trim(NameDataFile))

    allocate(VarIn_VIII(nVar,n1,n2,n3))

    call read_plot_file(NameFile=NameDataFile, &
         TypeFileIn = TypeDataFile,            &
         VarOut_VIII = VarIn_VIII,             &
         CoordMinOut_D = CoordMin_D,           &
         CoordMaxOut_D = CoordMax_D,           &   
         iErrorOut = iError)

    if(iError /= 0) call CON_stop( &
         NameSub//' could not read data from '//trim(NameDataFile))

    ! Assign angles of rotated box
    Xangle = Param_I(1)
    Yangle = Param_I(2)
    Zangle = Param_I(3)

    ! Assign var names to indexes, drop unused data, convert to SI
    call split_string(NameVar, MaxNameVar, NameVar_V, nVarName)

    if(IsDataBlock)then
       n1 = n1Block
       n2 = n2Block
       n3 = n3Block
       allocate(Var_VIII(nLocalVar,n1,n2,n3))
    else
       allocate(Var_VIII(nLocalVar,n1,n2,n3))
    endif

    ! Set up bins for Spectrum
    allocate(SpectrumTable_I(nWavelengthinterval))
    nBin = 0
    do iWavelengthInterval=1, nWavelengthInterval
       MinWavelength = WavelengthInterval_II(1,iWavelengthInterval)
       MaxWavelength = WavelengthInterval_II(2,iWavelengthInterval)
       nWavelengthBin = nint((MaxWavelength-MinWavelength)/SizeWavelengthBin)
       nBin = nWavelengthBin

       allocate(SpectrumTable_I(iWavelengthinterval)%Spectrum_III(n2,n3, &
            nWavelengthBin))
       allocate(SpectrumTable_I(iWavelengthinterval)%SpectrumGrid_I(&
            nWavelengthBin+1))

       SpectrumTable_I(iWavelengthinterval)%Spectrum_III(:,:,:)=0.0
       SpectrumTable_I(iWavelengthinterval)%nBin=nBin
       do iWavelengthBin=1,nWavelengthBin+1
          ! Values of each wavelength bin correspond to the center of the bin
          SpectrumTable_I(iWavelengthinterval)%SpectrumGrid_I(iWavelengthBin)&
               = MinWavelength + (iWavelengthBin-.5)*SizeWavelengthBin
       end do
    end do

    if (IsUniData) then
       ! cm^-3 --> kg/m^3
       Var_VIII(rho_,1:n1,1:n2,1:n3)   = rhoUni * 1e6 * cProtonMass
       ! km/s --> m/s
       Var_VIII(ux_,1:n1,1:n2,1:n3)    = uxUni * 1e3
       Var_VIII(uy_,1:n1,1:n2,1:n3)    = uyUni * 1e3
       Var_VIII(uz_,1:n1,1:n2,1:n3)    = uzUni * 1e3
       ! G --> T
       Var_VIII(bx_,1:n1,1:n2,1:n3)    = bxUni * 1e-4
       Var_VIII(by_,1:n1,1:n2,1:n3)    = byUni * 1e-4
       Var_VIII(bz_,1:n1,1:n2,1:n3)    = bzUni * 1e-4
       ! K
       Var_VIII(tpar_,1:n1,1:n2,1:n3)  = tparUni
       Var_VIII(tperp_,1:n1,1:n2,1:n3) = tperpUni
       Var_VIII(te_,1:n1,1:n2,1:n3)    = teUni
       ! erg cm^-3 --> J m^-3
       Var_VIII(I01_,1:n1,1:n2,1:n3)   = I01Uni * 1e-1
       Var_VIII(I02_,1:n1,1:n2,1:n3)   = I02Uni * 1e-1

    else
       do iVar=1, nVar
          if(IsVerbose)write(*,*)'NameVar_V(iVar+nDim) = ',&
               NameVar_V(iVar+nDim)

          call lower_case(NameVar_V(iVar+nDim))

          select case(NameVar_V(iVar+nDim))
          case('rho')
             ! g/cm3 --> kg/m3
             Var_VIII(rho_,1:n1,1:n2,1:n3)  = VarIn_VIII(iVar,1:n1,1:n2,1:n3) &
                  * 1e3
          case('ux')
             ! km/s --> m/s
             Var_VIII(ux_,1:n1,1:n2,1:n3)   = VarIn_VIII(iVar,1:n1,1:n2,1:n3) &
                  * 1e3
          case('uy')
             Var_VIII(uy_,1:n1,1:n2,1:n3)   = VarIn_VIII(iVar,1:n1,1:n2,1:n3) &
                  * 1e3
          case('uz')
             Var_VIII(uz_,1:n1,1:n2,1:n3)   = VarIn_VIII(iVar,1:n1,1:n2,1:n3) &
                  * 1e3
          case('bx')
             ! G --> T
             Var_VIII(bx_,1:n1,1:n2,1:n3)   = VarIn_VIII(iVar,1:n1,1:n2,1:n3) &
                  *1e-4
          case('by')
             Var_VIII(by_,1:n1,1:n2,1:n3)   = VarIn_VIII(iVar,1:n1,1:n2,1:n3) &
                  *1e-4
          case('bz')
             Var_VIII(bz_,1:n1,1:n2,1:n3)   = VarIn_VIII(iVar,1:n1,1:n2,1:n3) &
                  *1e-4
          case('ppar')
             ! dyne/cm2 --> N/m2
             ! T = p/rho * M / kB
             Var_VIII(tpar_,1:n1,1:n2,1:n3) = VarIn_VIII(iVar,1:n1,1:n2,1:n3) &
                  * 1e-1 / max(1e-30,Var_VIII(rho_,1:n1,1:n2,1:n3)) * &
                  cProtonMass &
                  / cBoltzmann
             IsPpar = .true.
          case('pperp')
             Var_VIII(tperp_,1:n1,1:n2,1:n3) = VarIn_VIII(iVar,1:n1,1:n2,1:n3)&
                  * 1e-1 / max(1e-30,Var_VIII(rho_,1:n1,1:n2,1:n3)) * &
                  cProtonMass &
                  / cBoltzmann
             IsPperp = .true.
          case('pe')
             Var_VIII(te_,1:n1,1:n2,1:n3)   = VarIn_VIII(iVar,1:n1,1:n2,1:n3) &
                  * 1e-1 / max(1e-30,Var_VIII(rho_,1:n1,1:n2,1:n3)) * &
                  cProtonMass &
                  / cBoltzmann
             IsPe = .true.
          case('p')
             Var_VIII(t_,1:n1,1:n2,1:n3) = VarIn_VIII(iVar,1:n1,1:n2,1:n3) &
                  * 1e-1 / max(1e-30,Var_VIII(rho_,1:n1,1:n2,1:n3)) * &
                  cProtonMass &
                  / cBoltzmann
          case('i01')
             ! erg/cm^3 --> J/m^3
             Var_VIII(I01_,1:n1,1:n2,1:n3)  = VarIn_VIII(iVar,1:n1,1:n2,1:n3) &
                  * 1e-1
          case('i02')
             Var_VIII(I02_,1:n1,1:n2,1:n3)  = VarIn_VIII(iVar,1:n1,1:n2,1:n3) &
                  * 1e-1
          case default
             write(*,*) NameSub // ' unused NameVar = '&
                  // NameVar_V(iVar+nDim)
          end select
       end do

       if(Xangle /= 0.0 .or. Yangle /= 0.0 .or. Zangle /= 0.0)then
          ! Apply rotation on vector variables
          Rot_DD = matmul(rot_matrix_z(-Zangle), &
               matmul(rot_matrix_y(-Yangle), rot_matrix_x(-Xangle)))

          if(IsDebug)then
             write(*,*)'angles = ',Param_I
             write(*,*)'ux_,uz_,bx_,bz_=', ux_,uz_,bx_,bz_
             write(*,*)'Rot_DD=', Rot_DD
             write(*,*)'Before: u_D=', Var_VIII(ux_:uz_,1,1,1)
          endif

          do k = 1, n3; do j = 1, n2; do i = 1, n1
             Var_VIII(ux_:uz_,i,j,k) = matmul(Rot_DD, Var_VIII(ux_:uz_,i,j,k))
             Var_VIII(bx_:bz_,i,j,k) = matmul(Rot_DD, Var_VIII(bx_:bz_,i,j,k))
          end do; end do; end do

          if(IsDebug)write(*,*)'After: u_D=', Var_VIII(ux_:uz_,1,1,1)

       end if
    endif

    deallocate(VarIn_VIII)

    ! Copy temperatures if not all components are given
    ! tpar = tper
    if(.not. IsPpar .and. .not. IsPperp)then
       Var_VIII(tpar_,1:n1,1:n2,1:n3) = Var_VIII(t_,1:n1,1:n2,1:n3)
       Var_VIII(tperp_,1:n1,1:n2,1:n3) = Var_VIII(t_,1:n1,1:n2,1:n3)
       ! tperp = (3*t - tpar)/2
    elseif(IsPpar .and. .not. IsPperp)then
       Var_VIII(tperp_,1:n1,1:n2,1:n3) = (3*Var_VIII(t_,1:n1,1:n2,1:n3) &
            - Var_VIII(tpar_,1:n1,1:n2,1:n3))/2.0
       ! tpar = 3*t - 2*tperp
    elseif(.not. IsPpar .and.IsPperp)then
       Var_VIII(tpar_,1:n1,1:n2,1:n3) = 3*Var_VIII(t_,1:n1,1:n2,1:n3) &
            - 2*Var_VIII(tperp_,1:n1,1:n2,1:n3)
    endif
    ! te = t
    if(.not. IsPe)Var_VIII(te_,1:n1,1:n2,1:n3) = Var_VIII(t_,1:n1,1:n2,1:n3)

    if(IsInstrument .and. nPixel /= n3)write(*,*) &
         '!!! nPixel= ',nPixel,' /= n3  = ',n3,' -> we use n3 '

    if(IsNoAlfven)then
       Var_VIII(I01_,1:n1,1:n2,1:n3) = 0
       Var_VIII(I02_,1:n1,1:n2,1:n3) = 0
       write(*,*)"IsNoAlfven ON !!!"
    endif
    
    ! Cells behind the solar disk
    do kPixel = 1, n3
       do jPixel = 1, n2
          do i = 1, n1
             if(CoordMin_D(iDimLOS) + &
                  (i-0.5) * (CoordMax_D(iDimLOS)-CoordMin_D(iDimLOS))/n1 < 0 .and. &
                  (CoordMin_D(iDimVertical) + (jPixel-0.5) * &
                  (CoordMax_D(iDimVertical)-CoordMin_D(iDimVertical))/n2)**2+&
                  (CoordMin_D(iDimHorizontal) + (kPixel-0.5) * &
                  (CoordMax_D(iDimHorizontal)-&
                  CoordMin_D(iDimHorizontal))/n3)**2 < 1 ) &
                  Var_VIII(:,i,jPixel,kPixel) = 0.
          end do
       end do
    end do

    dx = (CoordMax_D(iDimLOS)-CoordMin_D(iDimLOS))/n1 * rSun
    dy = (CoordMax_D(iDimVertical)-CoordMin_D(iDimVertical))/n2 * rSun
    dz = (CoordMax_D(iDimHorizontal)-CoordMin_D(iDimHorizontal))/n3 * rSun

    ! Constant multiplier converted from SI [m] to CGS units [cm]
    ! dV_cell / (1AU)^2 * 1e2
    dVperd2 = dx*dy*dz /(1.496e11)**2 * 1e2

  end subroutine read_data

  !==========================================================================
  subroutine read_table

    use ModIoUnit, ONLY: UnitTmp_
    use ModUtilities, ONLY: open_file, close_file

    ! Data read from the file
    character(len=6)            :: NameIon
    integer                     :: nLevelFrom, nLevelTo
    integer                     :: nFirstLevelFrom, nFirstLevelTo
    real                        :: LineWavelength, FirstLineWavelength
    real                        :: LogN, LogT, LogG, Aion 

    ! End of file indicated by iError /= 0
    integer                     :: iError

    ! Intensity table sized to the maximum
    real, allocatable           :: g_II(:,:)

    ! Size of table for a given wave
    integer                     :: iMin, jMin, iMax, jMax

    ! Density, temperature and wave indexes
    integer                     :: iN, iT, iLine

    ! Switches for header and information for a wavelength of interest
    logical                     :: IsHeader
    logical                     :: DoStore  

    character(len=*), parameter :: NameSub='read_table'
    !------------------------------------------------------------------------
    if(IsVerbose) write(*,*)'reading table file=', trim(NameTableFile)

    ! Read only wavelength of interest into a nice table
    allocate(LineTable_I(nMaxLine))
    nLineFound          = 0
    iLine               = 0
    nFirstLevelFrom     = -1
    nFirstLevelTo       = -1
    FirstLineWavelength = -1.0
    DoStore             = .false.
    IsHeader            = .true.

    ! Start to read data file
    call open_file(FILE=NameTableFile, STATUS='old')

    ! Read grid size information from header
    READGRID: do
       read(UnitTmp_,'(a)',iostat=iError) StringLine
       if(iError  /= 0)then
          write(*,*)'iError = ',iError
          call CON_stop('failed reading header of chianti table')
       end if
       if(StringLine == "#GRID")then
          read(UnitTmp_,*)MinLogN, MaxLogN, DLogN
          read(UnitTmp_,*)MinLogT, MaxLogT, DLogT
          if(IsVerbose)then
             write(*,*)'MinLogN, MaxLogN, DLogN = ',MinLogN, MaxLogN, DLogN 
             write(*,*)'MinLogT, MaxLogT, DLogT = ',MinLogT, MaxLogT, DLogT
          end if
          EXIT READGRID
       end if
    end do READGRID

    ! Set up maximum size grid to store tabulated values by wavelength
    MinI = nint(MinLogN/DLogN)
    MaxI = nint(MaxLogN/DLogN)
    MinJ = nint(MinLogT/DLogT)
    MaxJ = nint(MaxLogT/DLogT)
    allocate(g_II(MinI:MaxI,MinJ:MaxJ))

    ! Read remaining portion of table file
    READLOOP: do
       ! Read remaining header lines of table file
       if(IsHeader)then
          read(UnitTmp_,'(a)',iostat=iError) StringLine
          if(iError  /= 0)then
             write(*,*)'iError = ',iError
             call CON_stop('failed reading remaining header of chianti table')
          end if
          if(IsVerbose) write(*,'(a)') StringLine
          if(StringLine == "#START") IsHeader = .false.
          CYCLE READLOOP
       end if

       ! Read data line
       read(UnitTmp_,*,iostat=iError) &
            NameIon, Aion, nLevelFrom, nLevelTo, LineWavelength, &
            LogN, LogT, LogG

       if(iError  /= 0 .and. iError /= -1)then
          write(*,*)'iError = ',iError
          write(*,*)'last line = ',NameIon, Aion, nLevelFrom, nLevelTo, &
               LineWavelength, LogN, LogT, LogG
          call CON_stop('failed reading chianti table')
       end if

       ! Unobserved lines are stored with negative wavelength
       ! If interested in unobserved lines, use absolute value of wavelength
       if(IsAllLines)then
          if(IsVerbose .and. LineWavelength < 0)&
               write(*,*)'Unobserved line read : ',LineWavelength
          LineWavelength = abs(LineWavelength)
       endif

       ! Check if current line belongs to the same wavelength as previous one
       if(LineWavelength == FirstLineWavelength .and. &
            nLevelFrom == nFirstLevelFrom .and. &
            nLevelTo == nFirstLevelTo .and. &
            iError == 0) then
          ! Calculate index and store extra elements of LogG in oversized g_II
          iN = nint(LogN/DLogN)
          iT = nint(LogT/DLogT)
          g_II(iN,iT) = 10.0**LogG
          CYCLE READLOOP
       end if

       if(DoStore)then
          ! Store last indexes as the maximum indexes for N and T
          iMax = iN
          jMax = iT
          LineTable_I(iLine)%iMax = iMax
          LineTable_I(iLine)%jMax = jMax

          ! Extract min indexes into scalars
          iMin = LineTable_I(iLine)%iMin
          jMin = LineTable_I(iLine)%jMin

          ! Create intensity table
          allocate(LineTable_I(iLine)%g_II(iMin:iMax,jMin:jMax))
          LineTable_I(iLine)%g_II = g_II(iMin:iMax,jMin:jMax)

          ! Storage is done
          DoStore = .false.
       end if

       ! If wavelength is different than previous line
       ! pass only if wavelength is inside of wavelengthintervals of interest
       do iWavelengthInterval = 1, nWavelengthInterval
          if(IsDoppler)then
             ! If Doppler shift is calculated, use shifted intervals
             if(LineWavelength < &
                  WavelengthIntervalShifted_II(1,iWavelengthInterval))&
                  CYCLE
             if(LineWavelength > &
                  WavelengthIntervalShifted_II(2,iWavelengthInterval))&
                  CYCLE
          else
             ! If no Doppler shift is calculated, use original intervals
             if(LineWavelength < WavelengthInterval_II(1,iWavelengthInterval))&
                  CYCLE
             if(LineWavelength > WavelengthInterval_II(2,iWavelengthInterval))&
                  CYCLE
          endif

          ! New line of interest found, decide to store it
          iLine      = iLine + 1
          if(IsVerbose)write(*,*)iLine,NameIon,LineWavelength
          ! Check if there are too many lines already
          if(iLine > nMaxLine)&
               call CON_stop('Too many lines are found, increase MaxWave')

          ! Change reference line 
          FirstLineWavelength = LineWavelength
          nFirstLevelFrom = nLevelFrom
          nFirstLEvelTo = nLevelTo
          DoStore    = .true.
          nLineFound = nLineFound + 1

          ! Store ion name and wavelength 
          LineTable_I(iLine)%NameIon        = NameIon
          LineTable_I(iLine)%Aion           = Aion
          LineTable_I(iLine)%nLevelFrom     = nLevelFrom
          LineTable_I(iLine)%nLevelTo       = nLevelTo
          LineTable_I(iLine)%LineWavelength = LineWavelength

          if(IsVerbose)write(*,*)'NameIon, LineWavelength = ',NameIon,LineWavelength

          ! Calculate indexes and store as the minimum indexes
          iN                      = nint(LogN/DLogN)
          iT                      = nint(LogT/DLogT)
          LineTable_I(iLine)%StartLogT = LogT
          LineTable_I(iLine)%iMin = iN
          LineTable_I(iLine)%jMin = iT

          ! Initially zero out the input array
          g_II = 0.0

          ! Store first element
          g_II(iN,iT) = 10.0**LogG

          ! Go on reading the lines corresponding to the wavelength above
          EXIT

       end do

       ! When reached end of file, exit loop
       if(iError /= 0) EXIT READLOOP

    end do READLOOP
    call close_file

    deallocate(g_II)

    if(IsVerbose)write(*,*)'nLineFound = ',nLineFound
    if(IsVerbose)write(*,*)'nMaxLine = ',nMaxLine

  end subroutine read_table
  !==========================================================================
end program spectrum

!============================================================================

subroutine CON_stop(String)
  use ModMpi, ONLY: MPI_abort, MPI_COMM_WORLD
  implicit none
  integer                       :: iError, nError
  character (len=*), intent(in) :: String
  !--------------------------------------------------------------------------
  write(*,*)'CON_stop called on with String='
  write(*,*) String

!!!  call MPI_abort(MPI_COMM_WORLD, nError, iError)
  stop
  
end subroutine CON_stop
