!comments!!!
program spectrum

  use ModConst, ONLY           : cLightSpeed, cBoltzmann, cProtonMass, cAU, cPi

  implicit none

  logical                     :: IsVerbose =  .false., IsNoAlfven = .false.

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
  logical                     :: IsDataBlock = .false. ! read only part of data
  logical                     :: IsUniData = .false. ! overwrite data with const
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
  real                        :: tiUni, teUni, I01Uni, I02Uni

  ! Indexes for the solar wind variables
  integer, parameter          :: &
       rho_ = 1, &              
       ux_  = 2, &                  
       uy_  = 3, &               
       uz_  = 4, & 
       bx_  = 5, &
       by_  = 6, &
       bz_  = 7, &
       ti_  = 8, &
       te_  = 9, &
       I01_ =10, &
       I02_ =11

  ! Variables for the wavelengths of interest
  integer                     :: iWavelengthInterval, nWavelengthInterval
  integer                     :: nMaxLine, nLineFound
  real, allocatable           :: WavelengthInterval_II(:,:)

  ! Temperature and density grid for G
  real                        :: dLogN, dLogT, MaxLogN, MaxLogT, &
       MinLogN, MinLogT

  ! Maximum index range in density and temperature
  integer                     :: MaxI, MinI, MaxJ, MinJ

  ! Derived type to read tabulated G values
  type LineTableType
     character(len=6)         :: NameIon
     real                     :: LineWavelength
     integer                  :: iMin, jMin, iMax, jMax
     real, allocatable        :: g_II(:,:)
  end type LineTableType

  type(LineTableType), allocatable :: LineTable_I(:)

  real                        :: Dist = 0.99*cAU ! Sun-L1 distance
  real                        :: LOS_D(3), LOSnorm_D(3), dx, A

  ! Derived type of output 
  type SpectrumTableType
     real,allocatable         :: Spectrum_III(:,:,:), SpectrumGrid_I(:)
     real                     :: nBin
  end type SpectrumTableType

  type(SpectrumTableType), allocatable :: SpectrumTable_I(:)
  integer                     :: i, nBin
  real,allocatable            :: Spectrum_II(:,:), SpectrumGrid_I(:)

  integer                     :: iLine, jPixel, kPixel

  character(len=*), parameter :: NameSub = 'spectrum.f90'

  !---------------------------------------------------------------------------
  write(*,*)'Spectrum.exe starting'

  call MPI_init(iError)

  call read_param

  call read_data

  call read_table

  ! Loope over Chianti lines and pixels in the orthogonal direction
  do iLine = 1,min(nMaxLine,nLineFound)
     do kPixel = 1, n3
        do jPixel = 1, n2
           call calc_flux 
        end do
     end do
  end do

  call save_all

  write(*,*)'Spectrum.exe ending'

  call MPI_finalize(iError)

contains
  !==========================================================================
  subroutine save_all

    use ModPlotFile, ONLY          : save_plot_file

    character                      :: Namefile

    !------------------------------------------------------------------------

    ! Variables for output file
    character(len=200)   :: NameSpectrumFile = 'spectrum.out'
    character(len=5)     :: TypeFileSpectrum = 'ascii'
    character(len=200)   :: StringHeaderSpectrum = '[A] [erg sr^-1 cm^-2 A^-1]'
    real,allocatable            :: &
         Intensity_VIII(:,:,:,:), CoordWave_I(:), CoordPixel_DII(:,:,:)

    integer:: nWave, iWaveInterval, nWaveInterval, iWaveBin, iWave, nWaveBin, &
         nWaveAll
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
          ! Convert intensity to [erg cm^-2 sr^-1 s^-1 A^-1]
          Intensity_VIII(1,iWave,:,:) = &
               SpectrumTable_I(iWaveInterval)%Spectrum_III(:,:,iWaveBin) * &
               10.**7
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
    real                           :: LambdaSI, Lambda0SI, DeltaLambda, &
         DeltaLambda2
    real                           :: dLambdaInstr2 = 0.0 !!! do it later
    real                           :: zPlus2, zMinus2, cosAlpha
    real                           :: B_D(3), Bnorm_D(3)
    real                           :: uNth2, uTh2, Lambda
    real                           :: Gint, LogNe, LogTe, Rho
    real, allocatable              :: gLambda_II(:,:)

    character(len=*), parameter    :: NameSub='calc_flux'
    !------------------------------------------------------------------------

    allocate(gLambda_II(MinI:MaxI,MinJ:MaxJ))
       
    ! Convert to SI
    Lambda = LineTable_I(iLine)%LineWavelength
    
    LambdaSI = LineTable_I(iLine)%LineWavelength * 1e-10 

    ! Find which interval wavelength belongs to
    iInterval = -1
    do iWavelengthInterval = 1, nWavelengthInterval
       if((Lambda <= WavelengthInterval_II(2,iWavelengthInterval)) &
            .and.(Lambda >= WavelengthInterval_II(1,iWavelengthInterval))) then
          iInterval = iWavelengthInterval
          EXIT
       end if
    end do
    if(iInterval == -1)then
       write(*,*)"Lambda = ",Lambda
       write(*,*)"Intervals begin and end = ", WavelengthInterval_II(:,:)
       call CON_stop(NameSub // " no interval founnd")
    endif
       
    iBin =0

    do
       iBin=iBin+1
       if ((Lambda>SpectrumTable_I(iInterval)%SpectrumGrid_I(iBin)) &
            .and.(Lambda<SpectrumTable_I(iInterval)%SpectrumGrid_I(iBin+1)))then
          iCenter=iBin
          EXIT
       end if
       if (SpectrumTable_I(iInterval)%SpectrumGrid_I(iBin+1) > &
            WavelengthInterval_II(2,iInterval))EXIT
    end do

    do kPixel=1,n3
       do jPixel=1,n2
          do i=1,n1
             ! Calculate thermal and non-thermal broadening
             if(all(Var_VIII(:,i,jPixel,kPixel)==0))CYCLE ! cells inside body
             Rho = Var_VIII(rho_,i,jPixel,kPixel)
             zPlus2   = Var_VIII(I01_,i,jPixel,kPixel) * 4.0 / Rho
             zMinus2  = Var_VIII(I02_,i,jPixel,kPixel) * 4.0 / Rho
             B_D      = Var_VIII(bx_:bz_,i,jPixel,kPixel)
             Bnorm_D  = B_D/sqrt(max(sum(B_D**2), 1e-30))
             CosAlpha = sum(LOSnorm_D*Bnorm_D)
             uNth2    = 1.0/16.0 * (zPlus2 + zMinus2) * abs(cosAlpha)
             uTh2     = cBoltzmann * Var_VIII(ti_,i,jPixel,kPixel)/cProtonMass
             
             ! Convert from kg m^-3 to kg cm^-3 (*1e-6)
             ! and divide by cProtonMass in kg so Ne is in cm^-3  
             LogNe = log10(Rho*1e-6/cProtonMass)
             LogTe = log10(max(Var_VIII(te_,i,jPixel,kPixel),1e-30))

             ! Add instrumental broadening if there is any
             DeltaLambda2 = (LambdaSI/cLightSpeed)**2.0 * (uTh2 + uNth2)
             if(IsInstrument)DeltaLambda2 = DeltaLambda2 + dLambdaInstr2
             DeltaLambda = sqrt(DeltaLambda2)
             
             ! Get the contribution function
             iNMin  = LineTable_I(iLine)%iMin
             jTMin  = LineTable_I(iLine)%jMin
             iNMax  = LineTable_I(iLine)%iMax
             jTMax  = LineTable_I(iLine)%jMax
             gLambda_II = LineTable_I(iLine)%g_II(:,:)
             Gint = bilinear(gLambda_II, iNMin, iNMax, jTMin, jTMax, &
                  (/ LogNe/dLogN , LogTe/dLogT /),DoExtrapolate=.true.)

             ! Calculate flux and spread it on the Spectrum_II grids
             FluxMono = Ainstrument/(4*cPi*Dist**2.0)*Gint*(10.0**LogNe)**2.0*dx
             
             call disperse_line(iInterval, iCenter, LambdaSI, DeltaLambda, FluxMono)

          end do
       end do
    end do

  end subroutine calc_flux

  !==========================================================================
  subroutine disperse_line(iInterval,iCenter,LambdaSI,dLambdaSI,FluxMono)

    real, intent(in)            :: LambdaSI, dLambdaSI, FluxMono
    integer, intent(in)         :: iInterval, iCenter

    ! 
    


    integer                     :: iBin, iBegin, iEnd
    real                        :: Flux, Phi, InvNorm, InvSigma2
    real                        :: LambdaBin, LambdaBegin, LambdaEnd, LambdaDistSI
    integer                     :: iStep, iWave, nWaveBin

    character(len=*), parameter :: NameSub='disperse_line'
    !------------------------------------------------------------------------

    nWaveBin = SpectrumTable_I(iInterval)%nBin

    ! Beginning and end of gaussian truncated to +/-5 sigma in [A]
    LambdaBegin = (LambdaSI - 5*dLambdaSI)*1e10
    LambdaEnd   = (LambdaSI + 5*dLambdaSI)*1e10

    ! Find the corresponding wavelength bin for starting wavelength
    do iBegin = 1, nWaveBin - 1
       if (LambdaBegin < SpectrumTable_I(iInterval)%SpectrumGrid_I(iBegin+1)) EXIT
    end do

    ! Start at iBegin for efficiency and
    ! stop at nWaveBin-1 so iEnd = mWaveBin if no EXIT was performed
    do iEnd = iBegin, nWaveBin-1
       if (LambdaEnd < SpectrumTable_I(iInterval)%SpectrumGrid_I(iEnd)) EXIT
    end do

    InvNorm   = 1/(sqrt(2*cPi) * dLambdaSI)
    InvSigma2 = 1/(2*dLambdaSI**2) 
    
    ! Update bins between begin and end indices by adding the Gaussian distribution
    do iBin = iBegin , iEnd
       ! Get wavelength from the center of the bin
       LambdaBin = SpectrumTable_I(iInterval)%SpectrumGrid_I(iBin)
    
       ! Get distance from peak wavelength in SI
       LambdaDistSI = LambdaSI - LambdaBin*1e-10

       ! Calculate Gaussian of the line
       Phi = InvNorm * exp(-LambdaDistSI**2 * InvSigma2)
            
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
    logical                     :: IsNoInstrument = .false.
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
       case("#VERBOSE")
          call read_var('IsVerbose', IsVerbose)

       case("#DATAFILE")
          call read_var('NameDataFile',NameDataFile)
          call read_var('TypeDataFile',TypeDataFile)

       case("#TABLEFILE")
          call read_var('NameTableFile',NameTableFile)
          call read_var('nMaxLine',nMaxLine)

       case("#WAVELENGTHINTERVAL")
          IsNoInstrument = .true.
          Ainstrument = 7.7778e5**2
          call read_var('nWavelengthInterval',nWavelengthInterval)
          if(IsInstrument)then
             deallocate(WavelengthInterval_II)
             write(*,*)'INSTRUMENT intervals are changed to given WAVELENGTHINTERVALS'
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
             if(.not.IsDataBlock)then
                nPixel = 1024
             endif
             allocate(dLambdaInstr_I(nPixel))

             do iPixel=1,nPixel
                dLambdaInstr_I(iPixel) = 0.0
             end do
             Ainstrument = 7.7778e5**2

             SizeWavelengthBin = 0.0223
             if(IsNoInstrument)then
                write(*,*)'INSTRUMENT intervals are changed to given WAVELENGTHINTERVALS'
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
      call read_var('tiUni',tiUni)
      call read_var('teUni',teUni)
      call read_var('I01Uni',I01Uni)
      call read_var('I02Uni',I02Uni)

      case("#NOALFVEN")
         call read_var('IsNoAlfven',IsNoAlfven)

       case default
          write(*,*) NameSub // ' WARNING: unknown #COMMAND '

       end select
    end do READPARAM

  end subroutine read_param

  !==========================================================================
  subroutine read_data

    use ModPlotFile, ONLY: read_plot_file
    use ModUtilities, ONLY: split_string, lower_case

    character(len=lString)      :: StringHeader
    character(len=lString)      :: NameVar

    integer                     :: nStep

    real                        :: Time

    integer                     :: nVarName, iVar
    integer, parameter          :: MaxNameVar = 100
    character(len=20)           :: NameVar_V(MaxNameVar)
    
    real                        :: MinWavelength, MaxWavelength

    real,allocatable            :: VarIn_VIII(:,:,:,:)

    character(len=*), parameter :: NameSub = 'read_data'
    !------------------------------------------------------------------------
    call read_plot_file(NameFile=NameDataFile, &
         TypeFileIn=TypeDataFile,              &
         StringHeaderOut=StringHeader,         &
         nStepOut=nStep,                       &
         TimeOut=Time,                         &
         nDimOut=nDim,                         &
         nParamOut=nParam,                     &
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
         CoordMinOut_D=CoordMin_D,             &
         CoordMaxOut_D=CoordMax_D,             &   
         iErrorOut = iError)

    if(iError /= 0) call CON_stop( &
         NameSub//' could not read data from '//trim(NameDataFile))

    ! Assign var names to indexes, drop unused data, convert to SI
    call split_string(NameVar, MaxNameVar, NameVar_V, nVarName)

    if(IsDataBlock)then
       n1 = n1Block
       n2 = n2Block
       n3 = n3Block
       allocate(Var_VIII(11,n1,n2,n3))
    else
       allocate(Var_VIII(11,n1,n2,n3))
    endif

!    if(IsInstrument .and. nPixel /= n3)call CON_stop( &
!         NameSub//' incorrect number of vertical pixels for instrument in data file'//trim(NameInstrument))
    if(IsInstrument .and. nPixel /= n3)write(*,*)&
         '!!! nPixel /= n3 !!! nPixel = ',nPixel,' and n3 = ',n3
    

    ! Set up bins for Spectrum
    allocate(SpectrumTable_I(nWavelengthinterval))
    nBin = 0
    do iWavelengthInterval=1, nWavelengthInterval
       MinWavelength = WavelengthInterval_II(1,iWavelengthInterval)
       MaxWavelength = WavelengthInterval_II(2,iWavelengthInterval)
       nWavelengthBin = int((MaxWavelength-MinWavelength)/SizeWavelengthBin)
       nBin = nWavelengthBin
       
       allocate(SpectrumTable_I(iWavelengthinterval)%Spectrum_III(n2,n3,nWavelengthBin))
       allocate(SpectrumTable_I(iWavelengthinterval)%SpectrumGrid_I(nWavelengthBin+1))

       SpectrumTable_I(iWavelengthinterval)%Spectrum_III(:,:,:)=0.0
       SpectrumTable_I(iWavelengthinterval)%nBin=nBin
       do iWavelengthBin=1,nWavelengthBin+1
          ! Values of each wavelength bin correspond to the center of the bin
          SpectrumTable_I(iWavelengthinterval)%SpectrumGrid_I(iWavelengthBin) = &
               MinWavelength + (iWavelengthBin-.5)*SizeWavelengthBin
       end do
    end do

    if (IsUniData) then
       ! cm^-3 --> kgm^-3
       Var_VIII(rho_,1:n1,1:n2,1:n3) = rhoUni * 1e6 * cProtonMass
       ! km/s --> m/s
       Var_VIII(ux_,1:n1,1:n2,1:n3)  = uxUni * 1e3
       Var_VIII(uy_,1:n1,1:n2,1:n3)  = uyUni * 1e3
       Var_VIII(uz_,1:n1,1:n2,1:n3)  = uzUni * 1e3
       ! G --> T
       Var_VIII(bx_,1:n1,1:n2,1:n3)  = bxUni * 1e-4
       Var_VIII(by_,1:n1,1:n2,1:n3)  = byUni * 1e-4
       Var_VIII(bz_,1:n1,1:n2,1:n3)  = bzUni * 1e-4
       ! K
       Var_VIII(ti_,1:n1,1:n2,1:n3)  = tiUni
       Var_VIII(te_,1:n1,1:n2,1:n3)  = teUni
       ! erg/cm^3 --> J/m^3
       Var_VIII(I01_,1:n1,1:n2,1:n3) = I01Uni * 1e-1
       Var_VIII(I02_,1:n1,1:n2,1:n3) = I02Uni * 1e-1

    else
       do iVar=1, nVar
          if(IsVerbose)write(*,*)'NameVar_V(iVar+nDim) = ',NameVar_V(iVar+nDim)

          call lower_case(NameVar_V(iVar+nDim))

          select case(NameVar_V(iVar+nDim))
          case('rho')
             ! g/cm3 --> kg/m3
             Var_VIII(rho_,1:n1,1:n2,1:n3) = VarIn_VIII(iVar,1:n1,1:n2,1:n3)*1e3
          case('ux')
             ! km/s --> m/s
             Var_VIII(ux_,1:n1,1:n2,1:n3)  = VarIn_VIII(iVar,1:n1,1:n2,1:n3)*1e3
          case('uy')
             Var_VIII(uy_,1:n1,1:n2,1:n3)  = VarIn_VIII(iVar,1:n1,1:n2,1:n3)*1e3
          case('uz')
             Var_VIII(uz_,1:n1,1:n2,1:n3)  = VarIn_VIII(iVar,1:n1,1:n2,1:n3)*1e3
          case('bx')
             ! G --> T
             Var_VIII(bx_,1:n1,1:n2,1:n3)  = VarIn_VIII(iVar,1:n1,1:n2,1:n3) &
                  *1e-4
          case('by')
             Var_VIII(by_,1:n1,1:n2,1:n3)  = VarIn_VIII(iVar,1:n1,1:n2,1:n3) &
                  *1e-4
          case('bz')
             Var_VIII(bz_,1:n1,1:n2,1:n3)  = VarIn_VIII(iVar,1:n1,1:n2,1:n3) &
                  *1e-4
          case('ti')
             ! K
             Var_VIII(ti_,1:n1,1:n2,1:n3)  = VarIn_VIII(iVar,1:n1,1:n2,1:n3)
          case('te')
             Var_VIII(te_,1:n1,1:n2,1:n3)  = VarIn_VIII(iVar,1:n1,1:n2,1:n3)
          case('i01')
             ! erg/cm^3 --> J/m^3
             Var_VIII(I01_,1:n1,1:n2,1:n3) = VarIn_VIII(iVar,1:n1,1:n2,1:n3)*1e-1
          case('i02')
             Var_VIII(I02_,1:n1,1:n2,1:n3) = VarIn_VIII(iVar,1:n1,1:n2,1:n3)*1e-1
          case default
             write(*,*) NameSub // ' unused NameVar = ' // NameVar_V(iVar+nDim)
          end select
       end do
    endif

    if(IsNoAlfven)then
       Var_VIII(I01_,1:n1,1:n2,1:n3) = 0
       Var_VIII(I02_,1:n1,1:n2,1:n3) = 0
       write(*,*)"IsNoAlfven ON !!!"
    endif

    deallocate(VarIn_VIII)

    LOS_D = CoordMax_D-CoordMin_D
    LOSnorm_D  = LOS_D/sqrt(max(sum(LOS_D**2), 1e-30))

    dx = (CoordMax_D(1)-CoordMin_D(1))/n1

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

    MinI = nint(MinLogN/dLogN) + 1
    MaxI = nint(MaxLogN/dLogN) + 1
    MinJ = nint(MinLogT/dLogT) + 1
    MaxJ = nint(MaxLogT/dLogT) + 1
    allocate(g_II(MinI:MaxI,MinJ:MaxJ))
    if(IsVerbose)&
         write(*,*)'MinLogN, MaxLogN, dLogN, MinLogT, MaxLogT, dLogT = ',&
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

       ! Check if this belongs to the same line
       if(LineWavelength == FirstLineWavelength .and. iError == 0) then
          ! Calculate indexes and store extra elements of LogG
          iN = nint(LogN/dLogN) + 1
          iT = nint(LogT/dLogT) + 1
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

          if(IsVerbose)then
             write(*,*)'LineWavelength = ', LineTable_I(iLine)%LineWavelength
             write(*,*)'log(g_II(iMin,jMin:jMax)) = ', &
                  log10(g_II(iMin,jMin:jMax))
             write(*,*)'log(g_II(iMax,jMin:jMax)) = ', &
                  log10(g_II(iMax,jMin:jMax))
          endif
       end if

       if(iError /= 0) EXIT READLOOP

       ! Check if wavelength is inside any of the intervals
       do iWavelengthInterval = 1, nWavelengthInterval
          if(LineWavelength < WavelengthInterval_II(1,iWavelengthInterval)) CYCLE
          if(LineWavelength > WavelengthInterval_II(2,iWavelengthInterval)) CYCLE

          ! New line of interest found
          iLine = iLine + 1
          if(iLine > nMaxLine) call CON_stop('Too many waves, increase MaxWave')
          
          FirstLineWavelength = LineWavelength
          DoStore = .true.
          nLineFound = nLineFound + 1
          ! Store ion name and wavelength 
          LineTable_I(iLine)%NameIon    = NameIon
          LineTable_I(iLine)%LineWavelength = LineWavelength
          
          ! Calculate indexes and store as the minimum indexes
          iN = nint(LogN/dLogN) + 1; iT = nint(LogT/dLogT) + 1
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
