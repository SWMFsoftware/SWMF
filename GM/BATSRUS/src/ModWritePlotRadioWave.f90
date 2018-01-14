!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModWritePlotRadiowave
  use ModRadioWaveImage, ONLY: ray_bunch_intensity, nRay,      &
       Intensity_I, check_allocate,  rIntegration2, StateIn_VI,&
       SlopeX_, SlopeZ_
  use BATL_lib, ONLY: test_start, test_stop, x_, y_, z_
  implicit none
  SAVE

  private ! Except
  public:: write_plot_radiowave
contains
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
    use ModCoordTransform, ONLY: cross_product
    use ModProcMH, ONLY: iProc
    use ModMain, ONLY: Time_Accurate, n_Step, Time_Simulation
    use ModIO, ONLY: StringRadioFrequency_I, plot_type1, &
         plot_type, plot_form, plot_, ObsPos_DI, &
         n_Pix_X, n_Pix_Y, X_Size_Image, Y_Size_Image, &
         NamePlotDir, StringDateOrTime, nPlotRfrFreqMax
    use ModUtilities, ONLY: open_file, close_file
    use ModIoUnit, ONLY: UnitTmp_
    use ModPlotFile, ONLY: save_plot_file
    use ModCellGradient, ONLY: calc_gradient_ghost
    use ModAdvance, ONLY: State_VGB
    use ModVarIndexes, ONLY: Rho_
    use BATL_lib, ONLY:  nI, nJ, nK, MaxDim
    use ModUtilities, ONLY: join_string
    use ModWaves,     ONLY: WaveFirst_, WaveLast_, FrequencySi_W
    !
    ! Arguments
    !
    integer, intent(in) :: iFile

    !\
    ! Local variables
    !/
    !\
    ! Observer's location; actually it is the location 
    !     of the radiotelescope.
    !/
    real :: XyzObserver_D(MaxDim)
    !\
    !Misc functions of inputs 
    !/     
    real :: ObserverDistance !=sqrt(sum(XyzObs_D**2))
    !\
    ! Radius of "integration sphere"
    !/
    real :: rIntegration
    !\
    ! Dimensions of the raster in pixels
    !/
    integer ::  nXPixel, nYPixel
    !\
    !Pixel grid
    !/ 
    real :: XLower, YLower, XUpper, YUpper
    real :: DxPixel, DyPixel
    ! Pixel coordinates INSIDE the image plane
    real :: XPixel, YPixel  
    ! Unity vector normal to the image plane 
    real :: Normal_D(MaxDim)                                          
    real :: ZAxisOrt_D(MaxDim) = (/0, 0, 1/)
    ! Image plane inner Cartesian orts
    real :: Tau_D(MaxDim), Xi_D(MaxDim) 
    !\
    ! Pixel 3D Cartesian coordinates
    real, dimension(MaxDim) :: XyzPixel_D
    !\
    ! Initial directions of rays
    !/
    real :: SlopeUnscaled_D(MaxDim)
    !\
    ! Distance from the radiotelescope to the integration sphere  
    real :: ObservToIntSphereDist 
    !\
    ! Intersection point of the line of sight with the integration sphere
    !/
    real :: XyzAtIntSphere_D(MaxDim)
    !\
    !  Number of frequencies read from StringRadioFrequency_I(iFile)
    !/
    integer :: nFreq
    !\
    ! Frequencies in Hertz:
    !/
    real               :: RadioFrequency_I(nPlotRfrFreqMax)
 
    character (LEN=20) :: NameVar_I(nPlotRfrFreqMax)
    !\
    ! The result of the emissivity integration
    !/
    real, allocatable, dimension(:,:,:) :: Intensity_IIV
    !Loop variables
    integer :: iFreq, i, j, iRay
    ! Radius of "integration sphere"
    real    :: ImagePlaneDiagRadius

    character (LEN=120) :: NameVarAll, NameFile
    character (LEN=4)   :: NameFileExtension
    character (LEN=40)  :: NameFileFormat
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'write_plot_radiowave'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    !
    ! Initialize
    !
    call timing_start('write_plot_radiowave')

    !
    ! Set file specific parameters
    !
    XyzObserver_D = ObsPos_DI(:,iFile)
    nXPixel = n_Pix_X(iFile)
    nYPixel = n_Pix_Y(iFile)
    !
    ! Total number of rays and pixels in the raster
    !
    nRay = nXPixel*nYPixel
    call check_allocate
    !
    ! Determine the image plane inner coordinates of pixel centers
    !
    XLower = -0.50*X_Size_Image(iFile)
    YLower = -0.50*Y_Size_Image(iFile)
    XUpper =  0.50*X_Size_Image(iFile)
    YUpper =  0.50*Y_Size_Image(iFile)
    DxPixel = (XUpper - XLower)/nXPixel
    DyPixel = (YUpper - YLower)/nYPixel

    !
    ! Determune the orts, Tau and Xi, of the inner coordinate system
    ! of the image plane
    !
    ObserverDistance = sqrt(sum(XyzObserver_D**2))
    Normal_D = XyzObserver_D/ObserverDistance
    Tau_D = cross_product(ZAxisOrt_D, Normal_D)
    Tau_D = Tau_D/sqrt(sum(Tau_D**2))
    Xi_D = cross_product(Normal_D, Tau_D)
    !\
    ! Find the radius of the integration domain
    ImagePlaneDiagRadius = sqrt(xUpper**2 + yUpper**2)
    rIntegration = ceiling(max(ImagePlaneDiagRadius+1.0, 5.0))
    iRay = 0
    !
    ! Calculate coordinates of all the pixels 
    !
    do j = 1, nYPixel
       YPixel = YLower + (real(j) - 0.5)*DyPixel
       do i = 1, nXPixel
          XPixel = XLower + (real(i) - 0.5)*DxPixel
          XyzPixel_D = XPixel*Tau_D + YPixel*Xi_D
          SlopeUnscaled_D = XyzPixel_D - XyzObserver_D
          iRay = iRay + 1
          !Store direction vectors
          StateIn_VI(SlopeX_:SlopeZ_,iRay) = &
               SlopeUnscaled_D/sqrt(sum(SlopeUnscaled_D**2)) 
          !\
          ! Find the points on the integration sphere where it intersects
          ! with the straight "rays" 
          !/
          !\
          ! Solve a quadratic equation,
          !   XyzObs_D + ObservToIntSphereDist*Slope_DI || = rIntegration
          !or  ObservToIntSphereDist**2 + 
          !  + 2*ObservToIntSphereDist*XyzObs_D
          !  + ObserverDistance**2 - rIntegration**2 = 0
          !/  
          ObservToIntSphereDist = -sum(XyzObserver_D*&
               StateIn_VI(SlopeX_:SlopeZ_,iRay))&
               - sqrt(rIntegration**2 - ObserverDistance**2 &
               + sum(StateIn_VI(SlopeX_:SlopeZ_,iRay)&
               *XyzObserver_D)**2)
          !So the points at the boundary of integration sphere
          StateIn_VI(x_:z_,iRay) = XyzObserver_D &
               + StateIn_VI(SlopeX_:SlopeZ_,iRay)&
               *ObservToIntSphereDist
       end do
    end do
 
    call parse_freq_string(StringRadioFrequency_I(iFile), RadioFrequency_I, &
         NameVar_I, nFreq)

    if (iProc == 0) then
       if(DoTest)then
          write(*,*) 'XyzObserv_D     =', XyzObserver_D
          write(*,*)  '(XLower, YLower, XUpper, YUpper)=',&
                XLower, YLower, XUpper, YUpper    
          write(*,*) 'nXPixel        =', nXPixel
          write(*,*) 'nYPixel        =', nYPixel
       end if
       select case(plot_form(ifile))
       case('tec')
          NameFileExtension='.dat'
       case('idl')
          NameFileExtension='.out'
       end select
       call join_string(NameVar_I(1:nFreq), NameVarAll)

       plot_type1=plot_type(ifile)

       write(*,*) 'iFile = ', iFile
       write(*,*) 'nFreq = ', nFreq
       write(*,*) 'NameVarAll = ', NameVarAll
       write(*,*) 'RadioFrequency_I = ',RadioFrequency_I(1:nFreq)
    end if
    allocate(Intensity_IIV(nXPixel,nYPixel,nFreq))
    Intensity_IIV = 0.0

    ! Get density gradient
    call calc_gradient_ghost(State_VGB(Rho_,1:nI,1:nJ,1:nK,:))
    if (DoTest) write(*,*) 'rIntegration = ', rIntegration 
    rIntegration2 = rIntegration**2 + 0.01
    do iFreq = 1, nFreq
       ! Calculate approximate radius of the  critical surface around the sun
       ! from the frequency
       if (iProc == 0)then
          write(*,*) 'RAYTRACE START: RadioFrequency = ', &
               RadioFrequency_I(iFreq)
          write(*,*) 'RAYTRACE START: ImagePlaneDiagRadius = ', &
               ImagePlaneDiagRadius
          write(*,*) 'RAYTRACE START: rIntegration = ', &
               rIntegration
       end if
       FrequencySi_W(WaveFirst_) = RadioFrequency_I(iFreq)
       FrequencySi_W(WaveLast_ ) = 1
       call ray_bunch_intensity(RadioFrequency_I(iFreq))
       Intensity_IIV(:,:,iFreq) = &
            reshape(Intensity_I(1:nRay), (/nXPixel,nYPixel/))

       if (iProc == 0) write(*,*) 'RAYTRACE END'
    end do

    if (iProc==0) then
       if (ifile-plot_ > 9) then
          NameFileFormat='(a,i2,a,i7.7,a)'
       else
          NameFileFormat='(a,i1,a,i7.7,a)'
       end if

       if(time_accurate)then
          call get_time_string
          write(NameFile,NameFileFormat) &
               trim(NamePlotDir)//trim(plot_type1)//"_",&
               ifile-plot_,"_t"//trim(StringDateOrTime)//"_n",n_step,&
               NameFileExtension
       else
          write(NameFile,NameFileFormat) &
               trim(NamePlotDir)//trim(plot_type1)//"_",&
               ifile-plot_,"_n",n_step,NameFileExtension
       end if

       write(*,*) 'filename = ', NameFile
       !
       ! Write the file header
       !
       select case(plot_form(ifile))
       case('tec')
         ! description of file contains units, physics and dimension
          call save_plot_file(NameFile=NameFile,          &
               TypeFileIn = 'tec',                        &
               StringHeaderIn='BATSRUS: Radiotelescope Image', &
               nStepIn=n_step, TimeIn=Time_Simulation,    &
               NameVarIn='  X  Y '//trim(NameVarAll),     &
               CoordMinIn_D=&
               (/XLower + 0.5*DxPixel, YLower + 0.5*DyPixel/),&
               CoordMaxIn_D=&
               (/XUpper - 0.5*DxPixel, YUpper - 0.5*DyPixel/),&
               VarIn_IIV=Intensity_IIV, StringFormatIn = '(30(E14.6))')
       case('idl')
          ! description of file contains units, physics and dimension
          call save_plot_file(NameFile=NameFile,          &
               StringHeaderIn='RFR Radiorelescope Image', &
               nStepIn=n_step, TimeIn=Time_Simulation,    &
               NameVarIn='  X  Y '//trim(NameVarAll),     &
               CoordMinIn_D=&
               (/XLower + 0.5*DxPixel, YLower + 0.5*DyPixel/),&
               CoordMaxIn_D=&
               (/XUpper - 0.5*DxPixel, YUpper - 0.5*DyPixel/),&
               VarIn_IIV=Intensity_IIV, StringFormatIn = '(30es13.5)')
       end select

    end if  ! iProc ==0

    deallocate(Intensity_IIV)

    if (DoTest) write(*,*) 'write_plot_radiowave finished'

    call timing_stop('write_plot_radiowave')

    call test_stop(NameSub, DoTest)
  end subroutine write_plot_radiowave
  !============================================================================

  subroutine parse_freq_string(NameVarAll, Frequency_I, NameVar_I, nFreq)

    use ModIO, ONLY: nPlotRfrFreqMax

    !\
    ! INPUT
    !/
    ! String read from PARAM.in, like '1500kHz, 11MHz, 42.7MHz, 1.08GHz'

    character(len=*), intent(in) :: NameVarAll
    real,    intent(out) :: Frequency_I(nPlotRfrFreqMax)
    integer, intent(out) :: nFreq
    character(len=*), intent(out) :: NameVar_I(nPlotRfrFreqMax)
    character(len=50) :: cTmp, NameFreqUnit
    integer :: iFreq, lNameVarAll, iChar, iTmp
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'parse_freq_string'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    lNameVarAll = len(trim(NameVarAll))
    nFreq = 1
    iChar = 1

    ! Skip spaces, commas, or semicolons
    if (is_delim(NameVarAll(1:1))) then
       do while(is_delim(NameVarAll(iChar:iChar)) &
            .and. (iChar <= lNameVarAll))
          iChar = iChar + 1
       end do
    end if

    do while (iChar <= lNameVarAll)

       iTmp = 0
       do while(is_num(NameVarAll(iChar:iChar)) &
            .and. (iChar <= lNameVarAll))
          iTmp = iTmp + 1
          cTmp(iTmp:iTmp) = NameVarAll(iChar:iChar)
          iChar = iChar + 1
       end do

       read(cTmp(1:iTmp),*) Frequency_I(nFreq)

       do while(is_delim(NameVarAll(iChar:iChar)) &
            .and. (iChar <= lNameVarAll))
          iChar = iChar + 1
       end do

       iTmp = 0
       do while((.not. is_delim(NameVarAll(iChar:iChar))) &
            .and. (iChar <= lNameVarAll))
          iTmp = iTmp + 1
          cTmp(iTmp:iTmp) = NameVarAll(iChar:iChar)
          iChar = iChar + 1
       end do

       read(cTmp(1:iTmp),*) NameFreqUnit

       select case(trim(NameFreqUnit))
       case('Hz', 'HZ', 'hz')
          ! Do not scale
       case('kHz', 'kHZ', 'khz', 'KHz', 'Khz')
          Frequency_I(nFreq) = 1e3*Frequency_I(nFreq)
       case('MHz', 'MHZ', 'Mhz')
          Frequency_I(nFreq) = 1e6*Frequency_I(nFreq)
       case('GHz', 'GHZ', 'Ghz')
          Frequency_I(nFreq) = 1e9*Frequency_I(nFreq)
       case default
          write(*,*) '+++ Unrecognized frequency unit "'//trim(NameFreqUnit) &
               //'". Use only Hz, kHz, MHz, or GHz'
          stop
       end select

       do while(is_delim(NameVarAll(iChar:iChar)) &
            .and. (iChar <= lNameVarAll))
          iChar = iChar + 1
       end do

       nFreq = nFreq + 1

    end do
    nFreq = nFreq - 1
    !\
    ! Just in case: make all the frequencies positive
    !/
    Frequency_I = abs(Frequency_I)

    !
    ! Create standard frequency value array
    !
    do iFreq = 1, nFreq
       if ((Frequency_I(iFreq) > 0.0) .and. (Frequency_I(iFreq)<1e3)) then
          write(NameVar_I(iFreq),'(f6.2,a)') Frequency_I(iFreq), '_Hz'
       else if((Frequency_I(iFreq)>=1e3) .and. (Frequency_I(iFreq)<1e6)) then
          write(NameVar_I(iFreq),'(f6.2,a)') Frequency_I(iFreq)/1e3, '_kHz'
       else if((Frequency_I(iFreq) >= 1e6) .and. (Frequency_I(iFreq)<1e9)) then
          write(NameVar_I(iFreq),'(f6.2,a)') Frequency_I(iFreq)/1e6, '_MHz'
       else if((Frequency_I(iFreq) >= 1e9).and.(Frequency_I(iFreq)<1e12)) then
          write(NameVar_I(iFreq),'(f6.2,a)') Frequency_I(iFreq)/1e9, '_GHz'
       else if((Frequency_I(iFreq) >= 1e12).and.(Frequency_I(iFreq)<1e15))then
          write(NameVar_I(iFreq),'(f6.2,a)') Frequency_I(iFreq)/1e12, '_THz'
       end if
    end do

    do iFreq = 1, nFreq
       NameVar_I(iFreq) = 'f='//trim(adjustl(NameVar_I(iFreq)))
    end do

    call test_stop(NameSub, DoTest)
  contains
    !==========================================================================
    logical function is_num(c)
      character, intent(in) :: c

      !------------------------------------------------------------------------
      is_num = (lge(c, '0') .and. lle(c, '9')) .or. (c == '.') &
           .or. (c == 'e') .or. (c == 'E') &
           .or. (c == 'd') .or. (c == 'D') &
           .or. (c == '+') .or. (c == '-')

    end function is_num
    !==========================================================================
    logical function is_delim(c)
      character, intent(in) :: c

      !------------------------------------------------------------------------
      is_delim = (c == ' ') .or. (c == ',') .or. (c == ';')

    end function is_delim
    !==========================================================================
  end subroutine parse_freq_string
  !============================================================================

end module ModWritePlotRadiowave
!==============================================================================
