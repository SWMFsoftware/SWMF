!instrumental broadening
!comments!!!
!verbose levels
! isothermal.pro

program spectrum

  use ModConst, ONLY           : cLightSpeed, cBoltzmann, cProtonMass, cAU, &
       cPi, rSun

  implicit none
  logical                     :: IsVerbose =  .false., IsNoAlfven = .false.
  logical                     :: IsDataFile = .false., IsAllLines = .false.
  logical                     :: IsPermuteAxis = .false.

  ! Variables for output file
  character(len=200)   :: NameSpectrumFile = 'spectrum.out'
  character(len=200)   :: NameLabelFile = 'label.out'

  ! Variables for input files
  integer                     :: lString = 200 
  character(len=200)          :: StringLine
  character(len=200)          :: NameDataFile,NameTableFile
  character(len=200)          :: TypeDataFile

  ! Variables for instrument (if any)
  logical                     :: IsInstrument = .false.
  character(len=200)          :: NameInstrument
  integer                     :: nPixel
  real                        :: Ainstrument
  real,allocatable            :: dLambdaInstr_I(:)
  real                        :: SizeWavelengthBin
  integer                     :: nWavelengthBin, iWavelengthBin

  ! Variables for solar wind input data file
  logical                     :: IsDataBlock = .false. ! read part of data
  logical                     :: IsUniData = .false. ! overwrite data w/ const
  integer                     :: n1Block,n2Block,n3Block ! defines the size
  integer                     :: iError
  integer                     :: nVar   ! number of variables   
  integer                     :: nDim   ! number of dimensions  
  integer                     :: nParam ! number of parameters  
  integer                     :: n1, n2, n3 ! grid size
  real                        :: CoordMin_D(3), CoordMax_D(3)        
  real,allocatable            :: Var_VIII(:,:,:,:)
  real                        :: rhoUni, uxUni, uyUni, uzUni, bxUni, byUni
  real                        :: bzUni
  real                        :: tparUni, tperUni, teUni, I01Uni, I02Uni
  integer                     :: nTemperature = 2

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
       tper_ =  9, &
       te_   = 10, &
       I01_  = 11, &
       I02_  = 12
  
  integer                     ::  iDimLOS,iDimVertical,iDimHorizontal

  ! Variables for the wavelengths of interest
  integer                     :: iWavelengthInterval, nWavelengthInterval
  integer                     :: nMaxLine, nLineFound
  real, allocatable           :: WavelengthInterval_II(:,:), & 
       WavelengthIntervalShifted_II(:,:)

  ! Temperature and density grid for G
  real                        :: dLogN, dLogT, MaxLogN, MaxLogT, &
       MinLogN, MinLogT

  ! Maximum index range in density and temperature
  integer                     :: MaxI, MinI, MaxJ, MinJ

  ! Derived type to read tabulated G values
  type LineTableType
     character(len=6)         :: NameIon
     integer                  :: nLevelFrom, nLevelTo
     real                     :: LineWavelength
     integer                  :: iMin, jMin, iMax, jMax
     real, allocatable        :: g_II(:,:)
  end type LineTableType

  type(LineTableType), allocatable :: LineTable_I(:)

  real                        :: Dist = 0.99*cAU * 1e2 ! Sun-L1 distance in cm
  real                        :: LOS_D(3), LOSnorm_D(3)
  real                        :: dx, dA ! dx is plasma thickness along the LOS
  real                        :: dy, dz ! in case of boxdata used
  ! Derived type of output 
  type SpectrumTableType
     real,allocatable         :: Spectrum_III(:,:,:), SpectrumGrid_I(:)
     real                     :: nBin
  end type SpectrumTableType

  type(SpectrumTableType), allocatable :: SpectrumTable_I(:)
  integer                     :: i, nBin
  real,allocatable            :: Spectrum_II(:,:), SpectrumGrid_I(:)

  integer                     :: iLine, jPixel, kPixel

  logical                     :: IsDoppler = .true.

  ! Derived type for the Labels
  type LabelTableType
     character(len=6)         :: NameIon
     real                     :: Lambda
     real                     :: FluxMono
  end type LabelTableType

  type(LabelTableType), allocatable :: LabelTable_III(:,:,:)

  integer                     :: nLineAll

  character(len=*), parameter :: NameSub = 'spectrum.f90'

  !---------------------------------------------------------------------------
  write(*,*)'Spectrum.exe starting'

  call MPI_init(iError)

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
  do iLine = 1,nLineAll
     call calc_flux 
     if(IsVerbose)write(*,*)'calc_flux: line ',iLine,' out of nLineAll = ',nLineAll
  end do

  call save_all
  if(IsVerbose)write(*,*)'done with save_all'
  call save_label
  if(IsVerbose)write(*,*)'done with save_label'
  write(*,*)'Spectrum.exe ending'

  call MPI_finalize(iError)

contains


  !==========================================================================
  subroutine set_data_block
    integer, parameter             :: iUnitOut = 18
    integer                        :: MinWavelength, Maxwavelength
    character(len=*), parameter    :: NameSub = 'set_data_block'
    !------------------------------------------------------------------------

    if(IsVerbose)write(*,*)'dx, dA = ',dx,dA,' inside set_data_block'

    n1 = n1Block
    n2 = n2Block
    n3 = n3Block
    allocate(Var_VIII(12,n1,n2,n3))

    if (IsUniData) then
       ! cm^-3 --> kg/m^3
       Var_VIII(rho_,1:n1,1:n2,1:n3)  = rhoUni * 1e6 * cProtonMass
       ! km/s --> m/s
       Var_VIII(ux_,1:n1,1:n2,1:n3)   = uxUni * 1e3
       Var_VIII(uy_,1:n1,1:n2,1:n3)   = uyUni * 1e3
       Var_VIII(uz_,1:n1,1:n2,1:n3)   = uzUni * 1e3
       ! G --> T
       Var_VIII(bx_,1:n1,1:n2,1:n3)   = bxUni * 1e-4
       Var_VIII(by_,1:n1,1:n2,1:n3)   = byUni * 1e-4
       Var_VIII(bz_,1:n1,1:n2,1:n3)   = bzUni * 1e-4
       ! K
       Var_VIII(tpar_,1:n1,1:n2,1:n3) = tparUni
       Var_VIII(tper_,1:n1,1:n2,1:n3) = tperUni
       Var_VIII(te_,1:n1,1:n2,1:n3)   = teUni
       ! erg cm^-3 --> J m^-3
       Var_VIII(I01_,1:n1,1:n2,1:n3)  = I01Uni * 1e-1
       Var_VIII(I02_,1:n1,1:n2,1:n3)  = I02Uni * 1e-1
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
       nWavelengthBin = int((MaxWavelength-MinWavelength)/SizeWavelengthBin)
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
    real                        :: dWaveBin, WavelengthMin,WavelengthMax
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
       dWaveBin = (WavelengthMax - WavelengthMin)/nWaveBin

       ! Loop over wave lengths in this interval
       do iWaveBin = 1, nWaveBin
          ! Global wave length index
          iWave = nWave + iWaveBin

          ! Save intensity and coordinate into global array
          Intensity_VIII(1,iWave,:,:) = &
               SpectrumTable_I(iWaveInterval)%Spectrum_III(:,:,iWaveBin)
          CoordWave_I(iWave) = WavelengthMin + (iWaveBin - 0.5)*dWaveBin
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

    real                           :: FluxMono
    integer                        :: i, iInterval, iBin, iCenter
    integer                        :: iNMin, jTMin, iNMax, jTMax
    real                           :: LambdaSI, Lambda0SI, dLambdaSI, &
         dLambda, Lambda, dLambdaSI2, Lambda_shifted
    real                           :: dLambdaInstr2 = 0.0 !!! do it later
    real                           :: zPlus2, zMinus2, cosAlpha, sinAlpha
    real                           :: B_D(3), Bnorm_D(3)
    real                           :: uNth2, uTh2
    real                           :: Gint, LogNe, LogTe, Rho
    real, allocatable              :: gLambda_II(:,:)
    real                           :: FluxConst
    real                           :: Tlos

    character(len=*), parameter    :: NameSub='calc_flux'
    !------------------------------------------------------------------------

    allocate(gLambda_II(MinI:MaxI,MinJ:MaxJ))

    Lambda   = LineTable_I(iLine)%LineWavelength

    ! Fill in LabelTabel%NameIon
    LabelTable_III(iLine,:,:)%NameIon = LineTable_I(iLine)%NameIon

    ! Calculate flux from 1AU distant plasma parcel
    ! Convert arsec^-2 --> sr^-1 to match Chianti output
    ! 1 str = 4.25e10 arcsec^2
    ! 1 str^-1 = 2.3529412e-11 arcsec^-2 
    FluxConst = dA * dx / (4 * cPi * Dist**2)/2.3529412e-11

    do kPixel=1,n3
       do jPixel=1,n2
          do i=1,n1

             if(IsVerbose)write(*,*)'inside calc_flux loop begins i,j,k =',i,jPixel,kPixel

             ! Cells inside body
             if(all(Var_VIII(:,i,jPixel,kPixel)==0))EXIT 

             if(IsVerbose)write(*,*)'inside calc_flux loop after body-check'

             ! Doppler shift while x axis is oriented towards observer
             if(IsDoppler)then
                Lambda_shifted = (-Var_VIII(ux_,i,jPixel,kPixel)/cLightSpeed &
                     + 1 ) * Lambda
                Lambda = Lambda_shifted
             endif

             if(IsVerbose)write(*,*)'inside calc_flux loop after doppler-check'

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

             if(IsVerbose)write(*,*)'inside calc_flux loop after interval-check'

             iBin =0

             do
                iBin=iBin+1
                if(IsVerbose)write(*,*)'inside calc_flux loop iBin =',iBin
                if(IsVerbose)write(*,*)'inside calc_flux loop iInterval =',iInterval
                if(IsVerbose)write(*,*)'inside calc_flux loop nBin',nBin
                if ((Lambda>SpectrumTable_I(iInterval)%SpectrumGrid_I(iBin)) &
                     .and.&
                     (Lambda<SpectrumTable_I(iInterval)%SpectrumGrid_I(iBin+1)))&
                     then
                   iCenter=iBin
                   if(IsVerbose)write(*,*)'inside calc_flux loop iCenter =',iCenter
                   EXIT
                end if
                if (SpectrumTable_I(iInterval)%SpectrumGrid_I(iBin+1) > &
                     WavelengthInterval_II(2,iInterval))EXIT
                if (iBin==nBin)then
                   write(*,*)"Lambda = ",Lambda," will be left out!"
                   EXIT
                endif
             end do

             if(IsVerbose)write(*,*)'in calc_flux, before calculations'

             ! Calculate Elzasser variables
             Rho = Var_VIII(rho_,i,jPixel,kPixel)
             zPlus2   = Var_VIII(I01_,i,jPixel,kPixel) * 4.0 / Rho
             zMinus2  = Var_VIII(I02_,i,jPixel,kPixel) * 4.0 / Rho

             ! Calculate angle between LOS and B directions
             B_D      = Var_VIII(bx_:bz_,i,jPixel,kPixel)
             Bnorm_D  = B_D/sqrt(max(sum(B_D**2), 1e-30))
             cosAlpha = sum(LOSnorm_D*Bnorm_D)

             ! Calculate temperature relative to the LOS direction
             sinAlpha = sqrt(1 - cosAlpha**2)
             Tlos = sinAlpha * Var_VIII(tper_,i,jPixel,kPixel) &
                  + cosAlpha * Var_VIII(tpar_,i,jPixel,kPixel)

             ! Calculate thermal and non-thermal broadening
             uNth2    = 1.0/16.0 * (zPlus2 + zMinus2) * abs(cosAlpha)
             uTh2     = cBoltzmann * Tlos/cProtonMass

             ! Convert from kg m^-3 to kg cm^-3 (*1e-6)
             ! and divide by cProtonMass in kg so Ne is in cm^-3  
             LogNe = log10(Rho*1e-6/cProtonMass)
             LogTe = log10(Var_VIII(te_,i,jPixel,kPixel))

             ! Convert to SI
             LambdaSI = LineTable_I(iLine)%LineWavelength * 1e-10 

             ! Add instrumental broadening if there is any
             dLambdaSI2 = (LambdaSI/cLightSpeed)**2 * (uTh2 + uNth2)
             if(IsInstrument)dLambdaSI2 = dLambdaSI2 + dLambdaInstr2
             dLambdaSI = sqrt(dLambdaSI2)

             ! Convert [m] --> [A]
             dLambda   = dLambdaSI * 1e10

             ! Get the contribution function
             iNMin  = LineTable_I(iLine)%iMin
             jTMin  = LineTable_I(iLine)%jMin
             iNMax  = LineTable_I(iLine)%iMax
             jTMax  = LineTable_I(iLine)%jMax
             gLambda_II = LineTable_I(iLine)%g_II(:,:)
             Gint = bilinear(gLambda_II, iNMin, iNMax, jTMin, jTMax, &
                  (/ LogNe/dLogN , LogTe/dLogT /),DoExtrapolate=.true.)

             ! When Gint becomes negative due to extrapolation -> move to next
             if(Gint<0.0)CYCLE 

             ! Calculate flux and spread it on the Spectrum_II grids
             ! FluxMono = 1 / (4 * cPi * Dist**2.0) * &
             !      Gint * (10.0**LogNe)**2.0 * dV
             FluxMono = FluxConst * Gint * (10.0**LogNe)**2

             ! Fill in LabelTabel%Lambda,FluxMono
             LabelTable_III(iLine,:,:)%Lambda = Lambda
             LabelTable_III(iLine,jPixel,kPixel)%FluxMono = FluxMono

             call disperse_line(iInterval, iCenter, Lambda, dLambda, FluxMono)

          end do
       end do
    end do

    if(IsVerbose)then
       write(*,*)'iLine = ',iLine
       write(*,*)'NameIon = ',LineTable_I(iLine)%NameIon
       write(*,*)'LineWaveLength = ',LineTable_I(iLine)%LineWavelength
       write(*,*)'iMin = ',LineTable_I(iLine)%iMin, &
            'jMin = ',LineTable_I(iLine)%jMin
       write(*,*)'iMax = ',LineTable_I(iLine)%iMax, &
            'jMax = ',LineTable_I(iLine)%jMax
       write(*,*)'int(LogNe/dLogN)',int(LogNe/dLogN),'int(LogTe/dLogT)', & 
            int(LogTe/dLogT)
    endif

  end subroutine calc_flux

  !==========================================================================
  subroutine disperse_line(iInterval,iCenter,Lambda,dLambda,FluxMono)

    real, intent(in)            :: Lambda, dLambda, FluxMono
    integer, intent(in)         :: iInterval, iCenter

    ! 
    integer                     :: iBin, iBegin, iEnd
    real                        :: Flux, Phi, InvNorm, InvSigma2
    real                        :: LambdaBin, LambdaBegin, LambdaEnd, &
         LambdaDist
    integer                     :: iStep, iWave, nWaveBin

    character(len=*), parameter :: NameSub='disperse_line'
    !------------------------------------------------------------------------

    nWaveBin = SpectrumTable_I(iInterval)%nBin

    ! Beginning and end of gaussian truncated to +/-5 sigma in [A]
    LambdaBegin = (Lambda - 5*dLambda)
    LambdaEnd   = (Lambda + 5*dLambda)

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

    InvNorm   = 1/(sqrt(2*cPi) * dLambda)
    InvSigma2 = 1/(2*dLambda**2) 

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

    use ModReadParam

    character(len=lString)      :: NameCommand
    logical                     :: IsNoInstrument = .false., DoEcho = .false.
    integer                     :: iPixel
    character(len=*), parameter :: NameSub='read_param'
    !------------------------------------------------------------------------
    call read_file('SPECTRUM.in')

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
          call read_var('dx',dx)
          call read_var('dy',dy)
          call read_var('dz',dz)
          dx = dx*rSun*1e2 ! Convert to CGS
          dA = dy*rSun*1e2 *dz*rSun*1e2

       case("#INSTRUMENT")
          IsInstrument = .true.
          iPixel = 0
          nPixel = 0
          call read_var('NameInstrument',NameInstrument)

          select case(NameInstrument)
          case("EIS")
             if(.not.IsDataBlock)then
                nPixel = 1024
             endif
             allocate(dLambdaInstr_I(nPixel))

             do iPixel=1,nPixel
                dLambdaInstr_I(iPixel) = 0.0
             end do
             ! 1 arcsec x 1 arcsec on solar surface in cm^-2
             Ainstrument = (14e10/1800.)**2.

             SizeWavelengthBin = 0.0223
             if(IsNoInstrument)then
                write(*,*)'INSTRUMENT interval changed to WAVELENGTHINTERVALS'
             else
                allocate(WavelengthInterval_II(2,2))
                WavelengthInterval_II(:,1) = (/ 170 ,210 /)
                WavelengthInterval_II(:,2) = (/ 250 ,290 /)
             endif
          case default
             write(*,*) NameSub // ' WARNING: unknown #INSTRUMENT '
          end select

       case("#UNIFORMDATA")
          IsUniData = .true.
          call read_var('rhoUni',rhoUni)
          call read_var('uxUni',uxUni)
          call read_var('uyUni',uyUni)
          call read_var('uzUni',uzUni)
          call read_var('bxUni',bxUni)
          call read_var('byUni',byUni)
          call read_var('bzUni',bzUni)
          call read_var('tparUni',tparUni)
          call read_var('tperUni',tperUni)
          call read_var('teUni',teUni)
          call read_var('I01Uni',I01Uni)
          call read_var('I02Uni',I02Uni)

       case("#NOALFVEN")
          call read_var('IsNoAlfven',IsNoAlfven)

       case("#ISPERMUTEAXIS")
          call read_var('IsPermuteAxis',IsPermuteAxis)
          call read_var('iDimLOS',iDimLOS)
          call read_var('iDimVertical',iDimVertical)
          call read_var('iDimHorizontal',iDimHorizontal)

       case("#ISALLLINES")
          call read_var('IsAllLines',IsAllLines) ! Unobserved lines included

       case("#ISDOPPLER")
          call read_var('IsDoppler',IsDoppler)

       case("NTEMPERATURE")
          call read_var('nTemperature',nTemperature)

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

    character(len=lString)      :: StringHeader
    character(len=lString)      :: NameVar

    integer                     :: nStep, i, j, k, Size_D(3)

    integer                     :: nVarName, iVar
    integer, parameter          :: MaxNameVar = 100
    character(len=20)           :: NameVar_V(MaxNameVar)

    real                        :: MinWavelength, MaxWavelength

    real                        :: xAngle, yAngle, zAngle, Rot_DD(3,3)
    real                        :: Param_I(3)
    real,allocatable            :: VarIn_VIII(:,:,:,:)

    real                        :: Coord, Dz1, Dz2
    integer                     :: k1, k2

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
    xAngle = Param_I(1)
    yAngle = Param_I(2)
    zAngle = Param_I(3)

    ! Assign var names to indexes, drop unused data, convert to SI
    call split_string(NameVar, MaxNameVar, NameVar_V, nVarName)

    if(IsPermuteAxis)then
       Size_D(1:3) = (/n1,n2,n3/)
       CoordMax_D = (/ CoordMax_D(iDimLOS), CoordMax_D(iDimVertical), &
            CoordMax_D(iDimHorizontal)/)
       CoordMin_D = (/ CoordMin_D(iDimLOS), CoordMin_D(iDimVertical), &
            CoordMin_D(iDimHorizontal)/)
       n1 = Size_D(iDimLOS)
       n2 = Size_D(iDimVertical)
       n3 = Size_D(iDimHorizontal)
    endif

    if(IsDataBlock)then
       n1 = n1Block
       n2 = n2Block
       n3 = n3Block
       allocate(Var_VIII(12,n1,n2,n3))
    else
       allocate(Var_VIII(12,n1,n2,n3))
    endif

    ! Set up bins for Spectrum
    allocate(SpectrumTable_I(nWavelengthinterval))
    nBin = 0
    do iWavelengthInterval=1, nWavelengthInterval
       MinWavelength = WavelengthInterval_II(1,iWavelengthInterval)
       MaxWavelength = WavelengthInterval_II(2,iWavelengthInterval)
       nWavelengthBin = int((MaxWavelength-MinWavelength)/SizeWavelengthBin)
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
       Var_VIII(rho_,1:n1,1:n2,1:n3) = rhoUni * 1e6 * cProtonMass
       ! km/s --> m/s
       Var_VIII(ux_,1:n1,1:n2,1:n3)   = uxUni * 1e3
       Var_VIII(uy_,1:n1,1:n2,1:n3)   = uyUni * 1e3
       Var_VIII(uz_,1:n1,1:n2,1:n3)   = uzUni * 1e3
       ! G --> T
       Var_VIII(bx_,1:n1,1:n2,1:n3)   = bxUni * 1e-4
       Var_VIII(by_,1:n1,1:n2,1:n3)   = byUni * 1e-4
       Var_VIII(bz_,1:n1,1:n2,1:n3)   = bzUni * 1e-4
       ! K
       Var_VIII(tpar_,1:n1,1:n2,1:n3) = tparUni
       Var_VIII(tper_,1:n1,1:n2,1:n3) = tperUni
       Var_VIII(te_,1:n1,1:n2,1:n3)   = teUni
       ! erg cm^-3 --> J m^-3
       Var_VIII(I01_,1:n1,1:n2,1:n3)  = I01Uni * 1e-1
       Var_VIII(I02_,1:n1,1:n2,1:n3)  = I02Uni * 1e-1

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
                  * 1e-1 / Var_VIII(rho_,1:n1,1:n2,1:n3) * cProtonMass / cBoltzmann
          case('p')
             Var_VIII(tper_,1:n1,1:n2,1:n3) = VarIn_VIII(iVar,1:n1,1:n2,1:n3) &
                  * 1e-1 / Var_VIII(rho_,1:n1,1:n2,1:n3) * cProtonMass / cBoltzmann
          case('pe')
             Var_VIII(te_,1:n1,1:n2,1:n3)   = VarIn_VIII(iVar,1:n1,1:n2,1:n3) &
                  * 1e-1 / Var_VIII(rho_,1:n1,1:n2,1:n3) * cProtonMass / cBoltzmann
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

       if(xAngle /= 0.0 .or. yAngle /= 0.0 .or. zAngle /= 0.0)then
          ! Apply rotation on vector variables
          Rot_DD = matmul(rot_matrix_z(-zAngle), &
               matmul(rot_matrix_y(-yAngle), rot_matrix_x(-xAngle)))

          if(IsVerbose)then
             write(*,*)'angles = ',Param_I
             write(*,*)'ux_,uz_,bx_,bz_=', ux_,uz_,bx_,bz_
             write(*,*)'Rot_DD=', Rot_DD
             write(*,*)'Before: u_D=', Var_VIII(ux_:uz_,1,1,1)
          endif

          do k = 1, n3; do j = 1, n2; do i = 1, n1
             Var_VIII(ux_:uz_,i,j,k) = matmul(Rot_DD, Var_VIII(ux_:uz_,i,j,k))
             Var_VIII(bx_:bz_,i,j,k) = matmul(Rot_DD, Var_VIII(bx_:bz_,i,j,k))
          end do; end do; end do

          if(IsVerbose)write(*,*)'After: u_D=', Var_VIII(ux_:uz_,1,1,1)

       end if
    endif

    deallocate(VarIn_VIII)

    ! Copy temperatures if no three-temperature model is in use
    ! tpar = tper
    if(nTemperature == 2)Var_VIII(tpar_,1:n1,1:n2,1:n3) = Var_VIII(tper_,1:n1,1:n2,1:n3)
    ! tpar = te = tper
    if(nTemperature == 1)then
       Var_VIII(tpar_,1:n1,1:n2,1:n3) = Var_VIII(tper_,1:n1,1:n2,1:n3)
       Var_VIII(te_,1:n1,1:n2,1:n3)   = Var_VIII(tper_,1:n1,1:n2,1:n3)
    endif

    if(IsInstrument .and. nPixel /= n3) then
       write(*,*)'interpolate from n3 to nPixel! nPixel=',nPixel,' /= n3=',n3
       allocate(VarIn_VIII(nVar,n1,n2,n3))

       VarIn_VIII =  Var_VIII
       deallocate(Var_VIII)
       allocate(Var_VIII(nVar,n1,n2,nPixel))

       ! Interpolate to new grid
       do k = 1, nPixel
          Coord = 1 + real(n3-1)/real(nPixel-1)*real(k-1)
          k1 = floor(Coord)
          k2 = k1+1
          Dz1 = Coord - real(k1)
          Dz2 = 1.0 - Dz1

          Var_VIII(:,:,:,k) = Dz2*VarIn_VIII(:,:,:,k1)+Dz1*VarIn_VIII(:,:,:,k2)
       end do

       n3 = nPixel
       deallocate(VarIn_VIII)
    endif

    if(IsNoAlfven)then
       Var_VIII(I01_,1:n1,1:n2,1:n3) = 0
       Var_VIII(I02_,1:n1,1:n2,1:n3) = 0
       write(*,*)"IsNoAlfven ON !!!"
    endif

    ! LOS is set align with the x axis
    LOSnorm_D = (/1,0,0/) 

    ! Convert to CGS
    if(.not. IsDataBlock)then
       dx = (CoordMax_D(1) - CoordMin_D(1))/n1 *rSun*1e2
       dy = (CoordMax_D(2) - CoordMin_D(2))/n2 *rSun*1e2
       dz = (CoordMax_D(3) - CoordMin_D(3))/n3 *rSun*1e2
       if(dy*dz /= 0) then 
          dA = dy*dz 
       else 
          dA = max(dy**2,dz**2)
       end if
    endif

    if(IsVerbose)write(*,*)'dx, dA = ',dx,dA

  end subroutine read_data

  !==========================================================================
  subroutine read_table

    use ModIoUnit, ONLY: UnitTmp_

    ! Data read from the file
    character(len=6)            :: NameIon
    real                        :: LineWavelength, FirstLineWavelength
    integer                     :: nLevelFrom, nLevelTo
    real                        :: LogN, LogT, LogG 

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
    nLineFound      = 0
    iLine           = 0
    FirstLineWavelength = -1.0
    DoStore         = .false.
    IsHeader        = .true.

    ! Start to read data file
    open(UnitTmp_, FILE=NameTableFile)

    READGRID: do
       read(UnitTmp_,'(a)',iostat=iError) StringLine
       if(StringLine == "#GRID")then
          read(UnitTmp_,*)MinLogN, MaxLogN, dLogN
          read(UnitTmp_,*)MinLogT, MaxLogT, dLogT
          EXIT READGRID
       end if
    end do READGRID

    MinI = nint(MinLogN/dLogN)
    MaxI = nint(MaxLogN/dLogN)
    MinJ = nint(MinLogT/dLogT)
    MaxJ = nint(MaxLogT/dLogT)
    allocate(g_II(MinI:MaxI,MinJ:MaxJ))
    if(IsVerbose)&
         write(*,*)'MinLogN, MaxLogN, dLogN, MinLogT, MaxLogT, dLogT = ', &
         MinLogN, MaxLogN, dLogN,MinLogT, MaxLogT, dLogT

    READLOOP: do
       if(IsHeader)then
          read(UnitTmp_,'(a)',iostat=iError) StringLine
          if(IsVerbose) write(*,'(a)') StringLine
          if(StringLine == "#START") IsHeader = .false.
          CYCLE READLOOP
       end if

       read(UnitTmp_,*,iostat=iError) &
            NameIon, nLevelFrom, nLevelTo, LineWavelength, LogN, LogT, LogG
       ! If unobserved line, it is stored as a negative wavelength
       ! If interested in unobserved lines, use absolute value of wavelength
       if(IsAllLines)then
          if(IsVerbose .and. LineWavelength < 0)write(*,*)'Unobserved line read : ',LineWavelength
          LineWavelength = abs(LineWavelength)
       endif

       ! Check if this belongs to the same line
       if(LineWavelength == FirstLineWavelength .and. iError == 0) then
          ! Calculate indexes and store extra elements of LogG
          iN = nint(LogN/dLogN)
          iT = nint(LogT/dLogT)
          g_II(iN,iT) = 10.0**LogG
          CYCLE READLOOP
       end if

       if(DoStore)then
          ! Store last indexes as the maximum indexes for N and T
          iMax = iN; jMax = iT
          LineTable_I(iLine)%iMax = iMax; LineTable_I(iLine)%jMax = jMax

          ! Extract min indexes into scalars
          iMin = LineTable_I(iLine)%iMin; jMin = LineTable_I(iLine)%jMin

          ! Create intensity table
          allocate(LineTable_I(iLine)%g_II(iMin:iMax,jMin:jMax))
          LineTable_I(iLine)%g_II = g_II(iMin:iMax,jMin:jMax)

          ! Storage is done
          DoStore = .false.
       end if

       if(iError /= 0) EXIT READLOOP

       ! Check if wavelength is inside any of the intervals
       do iWavelengthInterval = 1, nWavelengthInterval
          if(IsDoppler)then
             if(LineWavelength < &
                  WavelengthIntervalShifted_II(1,iWavelengthInterval))&
                  CYCLE
             if(LineWavelength > &
                  WavelengthIntervalShifted_II(2,iWavelengthInterval))&
                  CYCLE
          else
             if(LineWavelength < WavelengthInterval_II(1,iWavelengthInterval))&
                  CYCLE
             if(LineWavelength > WavelengthInterval_II(2,iWavelengthInterval))&
                  CYCLE
          endif
          ! New line of interest found
          iLine = iLine + 1
          if(iLine > nMaxLine)&
               call CON_stop('Too many waves, increase MaxWave')

          FirstLineWavelength = LineWavelength
          DoStore = .true.
          nLineFound = nLineFound + 1

          ! Store ion name and wavelength 
          LineTable_I(iLine)%NameIon    = NameIon
          LineTable_I(iLine)%nLevelFrom = nLevelFrom
          LineTable_I(iLine)%nLevelTo = nLevelTo
          LineTable_I(iLine)%LineWavelength = LineWavelength

          ! Calculate indexes and store as the minimum indexes
          iN = nint(LogN/dLogN); iT = nint(LogT/dLogT)
          LineTable_I(iLine)%iMin = iN; LineTable_I(iLine)%jMin = iT

          ! To be safe zero out the input array
          g_II = 0.0

          ! Store first element
          g_II(iN,iT) = 10.0**LogG

          ! Wavelength was found to be interesting.
          EXIT

       end do

    end do READLOOP
    close(UnitTmp_)

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

  call MPI_abort(MPI_COMM_WORLD, nError, iError)
  stop
  
end subroutine CON_stop
