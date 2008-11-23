!^GFG COPYRIGHT UM
!============================================================================
subroutine write_plot_radiowave(iFile)

  !
  ! Purpose:  creates radio telescope images of the radiowaves at several
  !     frequencies by inegrating the plasma emissivity along the refracting
  !     rays.
  !     The plasma emissivity is considered here a function or the plasma 
  !     density only.
  ! Written by Leonid Benkevitch.
  !

  use ModProcMH, ONLY: iProc
  use ModMain, ONLY: Time_Accurate, n_Step, Time_Simulation
  use ModPhysics, ONLY: UnitX_, NO2IO_V
  use ModNumConst
  use ModConst
  use ModIO

  implicit none

  !
  ! Arguments
  !
  integer, intent(in) :: iFile

  !
  ! Local variables
  !
  real :: XyzObserv_D(3)                     
  ! Observ position
  real :: ImageRange_I(4)    
  ! Image plane: XLower, YLower, XUpper, YUpper
  real :: HalfImageRangeX, HalfImageRangeY
  real :: rIntegration           
  ! Radius of "integration sphere"
  integer ::  nXPixel, nYPixel
  ! Dimensions of the raster in pixels
  real, dimension(nPlotRfrFreqMax) :: RadioFrequency_I
  character (LEN=20) :: strFreqNames_I(nPlotRfrFreqMax)
  real, allocatable, dimension(:,:,:) :: Intensity_III
  ! The result of the emissivity integration
  integer :: nFreq
  ! Number of frequencies read from StringRadioFrequency_I(iFile)
  integer :: i, j, iFreq, iPixel, jPixel
  real :: XPixel, XPixelSize, YPixel, YPixelSize, nRay
  real :: XLower, YLower, XUpper, YUpper
  real :: ImagePlaneDiagRadius, RadiusCr, DensityCr
  real, parameter :: ProtonChargeSGSe = 4.8e-10 !SGSe
  real, parameter :: DensityAtSolarSurface = 3.3452e-16 ! g/cm^3
  character (LEN=120) :: allnames, StringHeadLine, strFreq
  character (LEN=500) :: unitstr_TEC, unitstr_IDL
  character (LEN=4) :: file_extension
  character (LEN=40) :: file_format
  logical :: oktest, oktest_me, DoTiming, DoTimingMe
  integer :: neqpar=1, nplotvar=1
  !--------------------------------------------------------

  !
  ! Initialize
  !
  call set_oktest('write_plot_radiowave', oktest,oktest_me)
  call set_oktest('rfr_timing', DoTiming, DoTimingMe)
  call timing_start('write_plot_radiowave')

  !
  ! Set file specific parameters 
  !
  XyzObserv_D = ObsPos_DI(:,iFile)
  nXPixel = n_Pix_X(iFile)
  nYPixel = n_Pix_Y(iFile)
  HalfImageRangeX = cHalf*X_Size_Image(iFile)
  HalfImageRangeY = cHalf*Y_Size_Image(iFile)
  ImageRange_I = (/-HalfImageRangeX, -HalfImageRangeY, &
       HalfImageRangeX, HalfImageRangeY/)
  nRay = nXPixel*nYPixel

  !
  ! Determine the image plane inner coordinates of pixel centers
  !
  XLower = ImageRange_I(1)
  YLower = ImageRange_I(2)
  XUpper = ImageRange_I(3)
  YUpper = ImageRange_I(4)
  XPixelSize = (XUpper - XLower)/nXPixel
  YPixelSize = (YUpper - YLower)/nYPixel

  call parse_freq_string(StringRadioFrequency_I(iFile), RadioFrequency_I, &
       strFreqNames_I, nFreq) 
  strFreq = ''
  do iFreq = 1, nFreq
     strFreq = trim(strFreq)//' "'//trim(adjustl(strFreqNames_I(iFreq)))//'"'
     if (iFreq .lt. nFreq) strFreq = trim(strFreq)//','
     if (iProc .eq. 0) then
        write(*,*) 'iFreq = ', iFreq, &
             ', trim(adjustl(strFreqNames_I(iFreq))) = ', &
             trim(adjustl(strFreqNames_I(iFreq))), ', strFreq = "', &
             trim(strFreq), '"'
     end if
  end do

  if (iProc .eq. 0) write(*,*) 'Preved 4'


  if (oktest .and. (iProc .eq. 0)) then
     write(*,*) 'XyzObserv_D     =', XyzObserv_D
     write(*,*) 'ImageRange_I   =', ImageRange_I, &
          '(XLower, YLower, XUpper, YUpper)'
     write(*,*) 'nXPixel        =', nXPixel
     write(*,*) 'nYPixel        =', nYPixel
     write(*,*) 'StringRadioFrequency_I(iFile): '
     write(*,*) StringRadioFrequency_I(iFile)
     write(*,*) ''
     do i = 1, nFreq
        write(*,*) RadioFrequency_I(i)
     end do
  end if

  unitstr_TEC = ''
  unitstr_IDL = ''

  plot_type1=plot_type(ifile)
  plot_vars1 = plot_vars(ifile)
  plot_pars1 = plot_pars(ifile)

  if (iProc .eq. 0) then
     write(*,*) 'iFile = ', iFile
     write(*,*) 'nFreq = ', nFreq
     write(*,*) 'StringRadioFrequency_I(iFile) = ', &
          StringRadioFrequency_I(iFile)
     write(*,*) 'strFreqNames_I:'
     do iFreq = 1, nFreq
        write(*,*) strFreqNames_I(iFreq)
     end do
     write(*,*) 'strFreq = ', strFreq
     write(*,*) 'RadioFrequency_I = '
     do iFreq = 1, nFreq
        write(*,*) RadioFrequency_I(iFreq)
     end do
  end if

  if (oktest_me) then
     write(*,*) 'iFile = ', iFile, ', plot_type = ', 'rfr '//strFreq, &
          ', form = ', plot_form(iFile)
     write(*,*) 'nFreq = ', nFreq
  end if

  ! Get the headers that contain variable names and units
  select case(plot_form(ifile))
  case('tec')
     !call get_TEC_rfr_variables(ifile,nplotvar,plotvarnames,unitstr_TEC)
     unitstr_TEC = 'VARIABLES = "X", "Y",'//strFreq
     if(oktest .and. iProc==0) write(*,*) unitstr_TEC
  case('idl')
     !call get_IDL_rfr_units(ifile,nplotvar,plotvarnames,unitstr_IDL)
     if(oktest .and. iProc==0) write(*,*) unitstr_IDL
  end select

  allocate(Intensity_III(nYPixel,nXPixel,nFreq))
  Intensity_III = cZero

  if (DoTiming) call timing_start('rfr_raytrace_loop')

  do iFreq = 1, nFreq
     ! Calculate approximate radius of the  critical surface around the sun
     ! from the frequency
     DensityCr = cPi*cProtonMass*cElectronMass*1e6* &
          (RadioFrequency_I(iFreq)/ProtonChargeSGSe)**2
     if (oktest_me) write(*,*) 'DensityCr = ', DensityCr
     RadiusCr = sqrt(DensityAtSolarSurface/DensityCr)
     if (oktest_me) write(*,*) 'RadiusCritical = ', RadiusCr
     ImagePlaneDiagRadius = sqrt(HalfImageRangeX**2 + HalfImageRangeY**2)
     rIntegration = ceiling(max(ImagePlaneDiagRadius+1.0, RadiusCr+1.0, 5.0))

     if (oktest_me) write(*,*) 'rIntegration = ', rIntegration

     if (iProc .eq. 0) write(*,*) 'RAYTRACE START: RadioFrequency = ', &
          RadioFrequency_I(iFreq)
     if (iProc .eq. 0) write(*,*) 'RAYTRACE START: DensityCr = ', DensityCr
     if (iProc .eq. 0) write(*,*) 'RAYTRACE START: RadiusCr = ', RadiusCr
     if (iProc .eq. 0) write(*,*) 'RAYTRACE START: ImagePlaneDiagRadius = ', &
          ImagePlaneDiagRadius 
     if (iProc .eq. 0) write(*,*) 'RAYTRACE START: rIntegration = ', &
          rIntegration

     call get_ray_bunch_intensity(XyzObserv_D, RadioFrequency_I(iFreq), &
          ImageRange_I, rIntegration, &
          nXPixel, nYPixel, Intensity_III(:,:,iFreq))

     if (iProc .eq. 0) write(*,*) 'RAYTRACE END'
     !if (iProc .eq. 0) then
     !   do iPixel = 1, nXPixel
     !      XPixel = XLower + (real(iPixel) - cHalf)*XPixelSize
     !      do jPixel = 1, nYPixel
     !         YPixel = YLower + (real(jPixel) - cHalf)*YPixelSize
     !         write(*,fmt="(30(E14.6))") XPixel, YPixel, &
     !              Intensity_III(jPixel,iPixel,1:nFreq)
     !      end do
     !   end do
     !end if
  end do

  if (DoTiming) call timing_stop('rfr_raytrace_loop')

  !
  ! Save results on file
  !
  if (DoTiming) call timing_start('rfr_save_plot')
  if (iProc==0) then

     select case(plot_form(ifile))
     case('tec')
        file_extension='.dat'
     case('idl')
        file_extension='.out'
     end select

     if (ifile-plot_ > 9) then
        file_format='("' // trim(NamePlotDir) // '",a,i2,a,i7.7,a)'
     else
        file_format='("' // trim(NamePlotDir) // '",a,i1,a,i7.7,a)'
     end if

     if(time_accurate)then
        call get_time_string
        write(filename,file_format) &
             trim(plot_type1)//"_",&
             ifile-plot_,"_t"//trim(StringDateOrTime)//"_n",n_step,&
             file_extension
     else
        write(filename,file_format) &
             trim(plot_type1)//"_",&
             ifile-plot_,"_n",n_step,file_extension
     end if

     write(*,*) 'filename = ', filename

     open(unit_tmp, file=filename, status="replace", err=999) !OPEN OPEN

     !
     ! Write the file header
     !
     select case(plot_form(ifile))
     case('tec')
        write(unit_tmp,*) 'TITLE="BATSRUS: Radiotelescope Image"'
        write(unit_tmp,'(a)') trim(unitstr_TEC)
        write(unit_tmp,*) 'ZONE T="RFR Image"', &
             ', I=',nXPixel,', J=',nYPixel,', F=POINT'
        ! Write point values
        do iPixel = 1, nXPixel
           XPixel = XLower + (real(iPixel) - cHalf)*XPixelSize
           do jPixel = 1, nYPixel
              YPixel = YLower + (real(jPixel) - cHalf)*YPixelSize
              write(unit_tmp,fmt="(30(E14.6))") XPixel, YPixel, &
                   Intensity_III(jPixel,iPixel,1:nFreq)
           end do
        end do

     case('idl')
        ! description of file contains units, physics and dimension
        StringHeadLine = 'RFR Radiorelescope Image'
        write(unit_tmp,"(a)") StringHeadLine

        ! 2 in the next line means 2 dimensional plot, 1 in the next line
        !  is for a 2D cut, in other words one dimension is left out)
        write(unit_tmp,"(i7,1pe13.5,3i3)") &
             n_step, time_simulation, 2, neqpar, nplotvar

        ! Grid size
        write(unit_tmp,"(2i4)") nXPixel, nYPixel

        ! Equation parameters
        ! write(unit_tmp,"(100(1pe13.5))") eqpar(1:neqpar)

        ! Coordinate, variable and equation parameter names
        write(unit_tmp,"(a)") allnames

        ! Data
        do iPixel = 1, nXPixel
           YPixel = YLower + (real(jPixel) - cHalf)*YPixelSize
           do jPixel = 1, nYPixel
              XPixel = XLower + (real(iPixel) - cHalf)*XPixelSize
              write(unit_tmp,fmt="(30(1pe13.5))") &
                   XPixel, YPixel, Intensity_III(jPixel,iPixel,1:nFreq)
           end do
        end do
     end select

     close(unit_tmp)  !CLOSE CLOSE CLOSE CLOSE CLOSE CLOSE CLOSE CLOSE CLOSE

  end if  !iProc ==0

  if (DoTiming) call timing_stop('rfr_save_plot')

  deallocate(Intensity_III)

  if (oktest_me) write(*,*) 'write_plot_radiowave finished'

  call timing_stop('write_plot_radiowave')

  return

999 continue

  call stop_mpi("Error in opening or writing file in write_plot_radiowave")

end subroutine write_plot_radiowave

!==========================================================================

subroutine parse_freq_string(StringRadioFrequency, Value_I, strValue_I, nValue)
  use ModIO, ONLY: nPlotRfrFreqMax
  implicit none
  ! String read from PARAM.in, like '1500kHz, 11MHz, 42.7MHz, 1.08GHz'
  character(len=*), intent(in) :: StringRadioFrequency
  real, intent(out) :: Value_I(nPlotRfrFreqMax)
  integer, intent(out) :: nValue
  character(len=*), intent(out) :: strValue_I(nPlotRfrFreqMax)
  character(len=50) :: cTmp, StrScale
  integer :: i, j, lenStr, iChar, iTmp


  lenStr = len(trim(StringRadioFrequency))
  nValue = 1
  iChar = 1

  ! Skip spaces, commas, or semicolons
  if (is_delim(StringRadioFrequency(1:1))) then
     do while(is_delim(StringRadioFrequency(iChar:iChar)) &
          .and. (iChar .le. lenStr))
        iChar = iChar + 1
     end do
  end if

  do while (iChar .le. lenStr)

     iTmp = 0
     do while(is_num(StringRadioFrequency(iChar:iChar)) &
          .and. (iChar .le. lenStr))
        iTmp = iTmp + 1
        cTmp(iTmp:iTmp) = StringRadioFrequency(iChar:iChar)
        iChar = iChar + 1
     end do

     read(cTmp(1:iTmp),*) Value_I(nValue)

     do while(is_delim(StringRadioFrequency(iChar:iChar)) &
          .and. (iChar .le. lenStr))
        iChar = iChar + 1
     end do

     iTmp = 0
     do while((.not. is_delim(StringRadioFrequency(iChar:iChar))) &
          .and. (iChar .le. lenStr))
        iTmp = iTmp + 1
        cTmp(iTmp:iTmp) = StringRadioFrequency(iChar:iChar)
        iChar = iChar + 1
     end do

     read(cTmp(1:iTmp),*) StrScale

     select case(trim(StrScale))
     case('Hz', 'HZ', 'hz')
        ! Do not scale
     case('kHz', 'kHZ', 'khz', 'KHz', 'Khz')
        Value_I(nValue) = 1e3*Value_I(nValue)
     case('MHz', 'MHZ', 'Mhz')
        Value_I(nValue) = 1e6*Value_I(nValue)
     case('GHz', 'GHZ', 'Ghz')
        Value_I(nValue) = 1e9*Value_I(nValue)
     case default
        write(*,*) '+++ Unrecognized frequency unit "'//trim(StrScale) &
             //'". Use only Hz, kHz, MHz, or GHz'
        stop
     end select

     do while(is_delim(StringRadioFrequency(iChar:iChar)) &
          .and. (iChar .le. lenStr))
        iChar = iChar + 1
     end do

     nValue = nValue + 1

  end do
  nValue = nValue - 1
  Value_I = abs(Value_I)    ! Just in case: make all the frequencies positive

  !
  ! Create standard frequency value array
  !
  do i = 1, nValue
     if ((Value_I(i) .gt. 0.0) .and. (Value_I(i) .lt. 1e3)) then
        write(strValue_I(i),'(a,f6.2,a)') Value_I(i), '_Hz' 
     else if ((Value_I(i) .ge. 1e3) .and. (Value_I(i) .lt. 1e6)) then
        write(strValue_I(i),'(a,f6.2,a)') Value_I(i)/1e3, '_kHz' 
     else if ((Value_I(i) .ge. 1e6) .and. (Value_I(i) .lt. 1e9)) then
        write(strValue_I(i),'(f6.2,a)') Value_I(i)/1e6, '_MHz' 
     else if ((Value_I(i) .ge. 1e9) .and. (Value_I(i) .lt. 1e12)) then
        write(strValue_I(i),'(f6.2,a)') Value_I(i)/1e9, '_GHz' 
     else if ((Value_I(i) .ge. 1e12) .and. (Value_I(i) .lt. 1e15)) then
        write(strValue_I(i),'(f6.2,a)') Value_I(i)/1e12, '_THz' 
     end if
  end do

  do i = 1, nValue
     strValue_I(i) = 'f_'//trim(adjustl(strValue_I(i)))
  end do

contains !===============================

  function is_num(c) result(yesno)
    character(len=1) :: c
    logical :: yesno
    yesno = (lge(c, '0') .and. lle(c, '9')) .or. (c .eq. '.') &
         .or. (c .eq. 'e') .or. (c .eq. 'E') &
         .or. (c .eq. 'd') .or. (c .eq. 'D') &
         .or. (c .eq. '+') .or. (c .eq. '-')
  end function is_num

  function is_delim(c) result(yesno)
    character(len=1) :: c
    logical :: yesno
    yesno = (c .eq. ' ') .or. (c .eq. ',') .or. (c .eq. ';')
  end function is_delim

end subroutine parse_freq_string
