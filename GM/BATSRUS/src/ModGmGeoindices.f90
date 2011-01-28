!==============================================================================
module ModGmGeoindices

  ! ModGmGeoindices is a module for calculating and writing geomagnetic 
  ! indices such as Kp, Dst, AsymH, etc. through the use of simulated
  ! ground-based measurments via ModGroundMagPerturb.  While some indices,
  ! such as Dst, can be calculated straight-forwardly (e.g. a Biort-Savart
  ! integration), this module attempts to create the index in a manner
  ! similar to how they are calculated in the real world.  

  implicit none
  save

  ! Should a geomagnetic indices file be written?
  logical :: DoWriteIndices = .false.

  ! What indices should be calculated and written?
  ! As more are added in the future (e.g. Dst, Asym-H, etc.), more options
  ! should be added to the #COMMAND so that desired indices are included
  ! but unnecessary ones may be shut off.
  logical :: DoCalcKp = .false., DoCalcDst = .false.
  
  ! K-index is evaluated over a rolling time window, typically three hours.
  ! It may be desirable to reduce this window, changing the Kp index so
  ! that it is more indicative of dB/dt.  This requires careful rescaling
  ! of the K-index conversion tables, however.  Units are left in minutes
  ! because the size of the window should be large, making seconds cumbersome.
  integer :: nKpMins = 180
  
  integer :: iUnitOut

  ! Output frequency.  The absense of a dn-type variable stems from how
  ! kp is calculated -- using a floating time window that is hard to define
  ! when kp is calculated on a constant iteration cadence instead of a 
  ! constant time cadence.
  real :: dtWriteIndices

  ! Has GeoIndices been initialized?
  logical :: IsInitialized = .false.
  
  ! VARIABLES FOR KP CALCULATION
  integer, parameter :: nKpMag = 24 ! Currently, only the faKe_p stations.
  integer            :: iSizeKpWindow = 0      ! Size of MagPerturb_II
  real, parameter    :: faKepLat = 50.0        ! Fake Kp geomag. latitude.
  real, allocatable  :: MagPerturb_II(:,:)     ! Magnetometer time history.
  real               :: MagPerbIE_DI(3,nKpMag)=0.0! IE contribution.
  real               :: XyzKp_DI(3, nKpMag)    ! Locations of kp mags, SMG.
  real               :: faKeP=0.0, LocalK(nKpMag)
  logical            :: IsFirstCalc=.true.
  ! K CONVERSION TABLES
  real, dimension(9), parameter :: &
       table50_i=(/5,10,20,40,70,120,200,330,500/) ! faKe_p, lat=50 deg.

  character(len=*), parameter :: NameKpVars = &
       'Kp K_12 K_13 K_14 K_15 K_16 K_17 K_18 K_19 K_20 K_21 K_22 K_23 '//&
       'K_00 K_01 K_02 K_03 K_04 K_05 K_06 K_07 K_08 K_09 K_10 K_11 '

contains

  !===========================================================================
  subroutine init_geoindices
    ! Initialize variables, arrays, and output file.
    use ModNumConst,  ONLY: cDegToRad, cTwoPi
    use ModUtilities, ONLY: flush_unit
    use ModMain,      ONLY: n_step
    use ModProcMH,    ONLY: iProc
    use ModIoUnit,    ONLY: io_unit_new
    use ModIO,        ONLY: NamePlotDir

    integer            :: i
    real               :: radXY, phi
    character(len=100) :: NameFile

    character(len=*), parameter :: NameSub='init_geoindices'
    logical :: DoTest, DoTestMe
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    IsInitialized=.true.

    ! Initialize grid and arrays.  FaKe_p uses stations fixed in SMG coords.
    XyzKp_DI(3,:) = sin(fakepLat * cDegToRad) ! SMG Z for all stations.
    radXY         = cos(fakepLat * cDegToRad) ! radial dist. from z-axis.
    do i=1, nKpMag
       phi = cTwoPi * (i-1)/24.0
       XyzKp_DI(1,i) = radXY * cos(phi)
       XyzKp_DI(2,i) = radXY * sin(phi)
       !if(iProc==0) &
       !write(*,'(a, 3(1x, e13.3))') 'Coords = ', XyzKp_DI(:,i)
    end do

    ! Allocate array to follow time history of magnetometer readings.
    iSizeKpWindow = int(nKpMins / (dtWriteIndices / 60.0))

    ! Allocate MagPerturb_II, open index file, write header.
    if (iProc==0) then
       allocate(MagPerturb_II(nKpMag, iSizeKpWindow))
       MagPerturb_II = 0.0

       write(NameFile, '(a, a, i8.8, a)') trim(NamePlotDir), &
            'GeoIndices_n', n_step, '.txt'
       iUnitOut = io_unit_new()
       open(iUnitOut, file=NameFile, status='replace')

       write(iUnitOut, '(2a,f8.2,a,i4.4)') 'Synthetic Geomagnetic Indices', &
            ' DtOutput=', DtWriteIndices, '   SizeKpWindow(Mins)=', &
            iSizeKpWindow
       write(iUnitOut, '(a)', advance='NO') &
            'it year mo dy hr mn sc msc '
       if (DoCalcKp) write(iUnitOut, '(a)', advance='NO') NameKpVars
       write(iUnitOut, '(a)') '' ! Close out header line.
       call flush_unit(iUnitOut)
    end if

  end subroutine init_geoindices

  !===========================================================================
  subroutine write_geoindices

    use ModMain,  ONLY: n_step, time_simulation
    use ModProcMH,ONLY: iProc, iComm
    use ModUtilities, ONLY: flush_unit
    use ModMpi

    integer :: i, iError, iTime_I(7)

    character(len=*), parameter :: NameSub='write_geoindices'
    logical :: DoTest, DoTestMe
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Calculate indices on all nodes.
    if(DoCalcKp) call calc_kp
    !if(DoCalcDst) call calc_dst ...etc...

    if(iProc > 0)RETURN ! Write only on head node.

    ! Write date and time.
    call get_date_time(iTime_I)
    write(iUnitOut, '(i7.7, 1x, i4.4, 5(1x, i2.2), 1x, i3.3)', ADVANCE='NO') &
         n_step, iTime_I

    if(DoCalcKp) then
       write(iUnitOut, '(1x, f4.2)', ADVANCE='NO') faKeP
       do i=1, nKpMag
          write(iUnitOut, '(1x, f4.2)', ADVANCE='NO') LocalK(i)
       end do
    end if

    write(iUnitOut, '(a)') ' ' ! Add carriage return.
    call flush_unit(iUnitOut)

  end subroutine write_geoindices

  !===========================================================================
  subroutine calc_kp
    use ModProcMH,ONLY: iProc, nProc, iComm
    use CON_axes, ONLY: transform_matrix
    use ModPhysics,        ONLY: No2Io_V, UnitB_
    use ModMain,           ONLY: time_simulation,TypeCoordSystem
    use ModMpi
    use ModGroundMagPerturb, ONLY: ground_mag_perturb, ground_mag_perturb_fac

    integer :: i, iError, nTmpMag
    real, dimension(3,3)       :: SmgToGsm_DD, GsmToSmg_DD, XyzToSph_DD
    real, dimension(3, nKpMag) :: Bmag_DI, Bfac_DI, Bsum_DI=0.0, XyzGsm_DI
    real :: deltaH
    character(len=*), parameter :: NameSub='calc_kp'
    logical :: DoTest, DoTestMe
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Obtain locations in correct (GSM) coordinates.
    SmgToGsm_DD = transform_matrix(Time_simulation, 'SMG', TypeCoordSystem)
    GsmToSmg_DD = transform_matrix(Time_simulation, TypeCoordSystem, 'SMG')
    do i=1, nKpMag
       XyzGsm_DI(:,i) = matmul(SmgToGsm_DD, XyzKp_DI(:,i))
    end do
  
    ! Obtain geomagnetic pertubation. Output is in NED Coordinates.
    call ground_mag_perturb(    nKpMag, XyzGsm_DI, Bmag_DI)
    call ground_mag_perturb_fac(nKpMag, XyzKp_DI,  Bfac_DI)

    ! Sum contributions.
    do i=1, nKpMag
       Bmag_DI(:,i) = (Bmag_DI(:,i) + Bfac_DI(:,i)) * No2Io_V(UnitB_)
    end do

    ! MPI Reduce to head node.
    if(nProc>1) call MPI_reduce(Bmag_DI, Bsum_DI, 3*nKpMag, &
         MPI_REAL, MPI_SUM, 0, iComm, iError)

    ! Head node calculates K-values and shares them with all other nodes.
    if(iProc==0)then
       ! Shift MagPerturb to make room for new measurements.
       MagPerturb_II(:,1:iSizeKpWindow-1) = MagPerturb_II(:,2:iSizeKpWindow)

       ! Add IE component of pertubation.
       do i=1, nKpMag
          ! Store H-component; add IE component.
          if (IsFirstCalc) then
             MagPerturb_II(i,:)=Bsum_DI(1,i)+MagPerbIE_DI(1,i)
          else
             MagPerturb_II(i,iSizeKpWindow)=Bsum_DI(1,i)+MagPerbIE_DI(1,i)
          end if

          ! Calculate deltaH, convert to K.
          deltaH = maxval(MagPerturb_II(i,:)) - minval(MagPerturb_II(i,:))
          LocalK(i) = convert_to_k(deltaH, table50_I)
       end do
    end if
    
    IsFirstCalc=.false.

    ! Kp is average of Ks.
    faKeP = sum(LocalK)/nKpMag
    ! Quantize to -/+ levels.
    faKeP = (nint(faKeP * 3.0))/3.0

  end subroutine calc_kp

  !===========================================================================
  subroutine finalize_geoindices
    ! De-allocate arrays, close files.
    use ModProcMH, ONLY: iProc
    character(len=*), parameter :: NameSub='finalize_geoindices'
    logical :: DoTest, DoTestMe
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    
    if(iProc==0)then
       close(iUnitOut)

       ! Clean up allocatables.
       deallocate(MagPerturb_II)
    end if

  end subroutine finalize_geoindices
  !===========================================================================
  real function convert_to_k(deltaB, table_I)!, Kout)
    ! Convert a deltaB value (max-min over window) to a K-value using given
    ! conversion table 'table'.  Table should be a 9-element vector that
    ! lists the upper limit for each K window.  For example, a K of 0 is
    ! given for a deltaB .le. 5, table(1) = 5.
    
    !real :: convert_to_k = 9
    real, intent(in) :: deltaB, table_I(9)
    !real, intent(out):: Kout=9 
    
    integer :: i
    !------------------------------------------------------------------------
    do i=1, 9
       if( (deltaB-table_I(i)) < 0) then
          convert_to_k=i-1
          exit
       end if
    end do
    
  end function convert_to_k
  
  !===========================================================================
 
end module ModGmGeoindices

!==============================================================================
