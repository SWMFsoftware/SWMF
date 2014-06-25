!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!This code is a copyright protected software (c) 2002- University of Michigan

module GM_couple_rb

  ! Coupling with Radiation Belt component

  use ModMpi
  use CON_coupler, ONLY: Grid_C, ncells_decomposition_d
  use ModProcMH, ONLY: iProc
  use ModMain, ONLY: n_step
  use ModPhysics, ONLY: No2Si_V, Si2No_V, UnitP_, UnitRho_, UnitTemperature_

  implicit none
  save

  private ! except

  public:: GM_get_for_rb_trace
  public:: GM_get_for_rb
  public:: GM_satinit_for_rb
  public:: GM_get_sat_for_rb

  ! Local variables

  character (len=*), parameter :: NameMod='GM_couple_rb'

  ! RB Grid size
  integer :: nCells_D(2), iSize,jSize

  ! Information about the RB grid ! 2D non-uniform regular grid only !!!
  real, allocatable, dimension(:) :: RB_lat, RB_lon

  integer :: i,j

  real, allocatable :: MHD_lat_boundary(:)
  real, dimension(:,:), allocatable :: &
       MHD_SUM_vol, &
       MHD_SUM_rho, &
       MHD_SUM_p, &
       MHD_Beq, &
       MHD_Xeq, &
       MHD_Yeq

  real, parameter :: noValue=-99999.

  integer, parameter :: maxMessages=10
  integer :: iError


  integer, parameter :: vol_=1, z0x_=2, z0y_=3, bmin_=4, rho_=5, p_=6

contains
  !===========================================================================
  subroutine allocate_gm_rb(iSizeIn,jSizeIn)
    use CON_comp_param, ONLY: RB_

    integer, intent(in) :: iSizeIn, jSizeIn
    character(len=*), parameter:: NameSub=NameMod//'::allocate_gm_rb'
    !--------------------------------------------------------------------------
    iSize = iSizeIn
    jSize = jSizeIn

    if(.not.allocated(RB_lat))then
       nCells_D=ncells_decomposition_d(RB_)
       if(  iSize /= nCells_D(1) .or. &
            jSize /= nCells_D(2) ) then
          write(*,*)NameSub//' grid sizes do not agree iSize,jSize,nCells=',&
               iSize,jSize, nCells_D(1:2)
          call CON_stop(NameSub//' ERROR')
       end if
       allocate(RB_lat(iSize), RB_lon(jSize))
       ! Convert colat, lon to lat-lon in degrees
       RB_lat = Grid_C(RB_) % Coord1_I
       RB_lon = Grid_C(RB_) % Coord2_I
    end if

    if(.not.allocated(MHD_SUM_vol))then
       allocate( MHD_SUM_vol(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_SUM_vol")
       MHD_SUM_vol = 0.
    end if

    if(.not.allocated(MHD_SUM_rho))then
       allocate( MHD_SUM_rho(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_SUM_rho")
       MHD_SUM_rho = 0.
    end if

    if(.not.allocated(MHD_SUM_p))then
       allocate( MHD_SUM_p(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_SUM_p")
       MHD_SUM_p = 0.
    end if

    if(.not.allocated(MHD_Beq))then
       allocate( MHD_Beq(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_Beq")
       MHD_Beq = 0.
    end if

    if(.not.allocated(MHD_Xeq))then
       allocate( MHD_Xeq(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_Xeq")
       MHD_Xeq = 0.
    end if

    if(.not.allocated(MHD_Yeq))then
       allocate( MHD_Yeq(isize,jsize), stat=iError )
       call alloc_check(iError,"MHD_Yeq")
       MHD_Yeq = 0.
    end if

    if(.not.allocated(MHD_lat_boundary))then
       allocate( MHD_lat_boundary(jsize), stat=iError )
       call alloc_check(iError,"MHD_lat_boundary")
       MHD_lat_boundary = 0
    end if

  end subroutine allocate_gm_rb

  !============================================================================
  subroutine write_integrated_data_tec

    use ModIoUnit, ONLY: UNITTMP_
    character(LEN=80) :: filename
    integer :: j2, nCall=0
    real :: tmpT, tmpV1,tmpV2, lonShift
    !-------------------------------------------------------------------------

    nCall=nCall+1

    !write values to plot file
    write(filename,'(a,i6.6,a,i4.4,a)')"RB/RbValues_n=",n_step,"_",nCall,".dat"

    OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown')
    write(UNITTMP_,'(a)') 'TITLE="Raytrace Values"'
    write(UNITTMP_,'(a)') &
         'VARIABLES="J", "I", "Lon", "Lat", "Lat Boundary (I)"', &
         ', "Xeq", "Yeq"', &
         ', "Volume", "Volume**(-2/3)"', &
         ', "MHD `r", "MHD p", "MHD T", "Beq"'
    write(UNITTMP_,'(a,i3.3,a,i4,a,i4,a)') &
         'ZONE T="PE=',iProc,'", I=',jsize+1,', J=',isize,', K=1, F=POINT'
    do i=1,isize
       do j2=1,jsize+1
          j=j2; if(j2==jsize+1) j=1
          lonShift=0.; if(j2==jsize+1) lonShift=360.
          tmpT=-1.; if(MHD_SUM_rho(i,j)>0.) &
               tmpT = ((MHD_SUM_p(i,j)*Si2No_V(UnitP_)) / &
               (MHD_SUM_rho(i,j)*Si2No_V(UnitRho_))) &
               * No2Si_V(UnitTemperature_)
          tmpV1=0.; if(MHD_SUM_vol(i,j)>0.) &
               tmpV1 = (MHD_SUM_vol(i,j)/1.e9)
          tmpV2=0.; if(MHD_SUM_vol(i,j)>0.) &
               tmpV2 = (MHD_SUM_vol(i,j)/1.e9)**(-2./3.)
          write(UNITTMP_,'(2i4,12G14.6)') j2,i,RB_lon(j)+lonShift,RB_lat(i), &
               MHD_lat_boundary(j), &
               MHD_Xeq(i,j),MHD_Yeq(i,j), &
               tmpV1,tmpV2, &
               MHD_SUM_rho(i,j),MHD_SUM_p(i,j),tmpT,MHD_Beq(i,j)
       end do
    end do
    CLOSE(UNITTMP_)

  end subroutine write_integrated_data_tec

  !============================================================================
  subroutine write_integrated_data_idl
    use ModPhysics, ONLY: No2Si_V, UnitB_
    use ModIoUnit, ONLY: UNITTMP_
    use ModMain,   ONLY: time_simulation
    CHARACTER (LEN=100) :: filename
    integer :: nCall = 0
    !-------------------------------------------------------------------------

    !write values to plot file
    nCall = nCall+1
    write(filename,'(a,i6.6,a,i4.4,a)')"RB/RbValues_n=",n_step,"_",nCall,".out"

    OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown', &
         iostat =iError)
    if (iError /= 0) call CON_stop("Can not open raytrace File "//filename)
    write(UNITTMP_,'(a79)')            'Raytrace Values_var22'
    write(UNITTMP_,'(i7,1pe13.5,3i3)') n_step,time_simulation,2,1,6
    write(UNITTMP_,'(3i4)')            jSize,iSize
    write(UNITTMP_,'(100(1pe13.5))')   0.0
    write(UNITTMP_,'(a79)')            'Lon Lat Xeq Yeq vol rho p Beq nothing'

    do i=isize,1,-1
       do j=1,jsize
          write(UNITTMP_,'(100(1pe18.10))') &
               modulo(RB_lon(j)-180.0,360.0),RB_lat(i), &
               MHD_Xeq(i,j),MHD_Yeq(i,j),&
               MHD_SUM_vol(i,j), &
               MHD_SUM_rho(i,j),MHD_SUM_p(i,j),MHD_Beq(i,j) * No2Si_V(UnitB_)
       end do
    end do
    CLOSE(UNITTMP_)
    
  end subroutine write_integrated_data_idl

  !==========================================================================
  subroutine GM_get_for_rb_trace(iSizeIn, jSizeIn, NameVar, nVarLine, &
       nPointLine)

    ! Do ray tracing for RB. 
    ! Provide total number of points along rays 
    ! and the number of variables to pass to RB
    use ModRayTrace, ONLY: DoExtractUnitSi

    use CON_line_extract, ONLY: line_get

    integer, intent(in)           :: iSizeIn, jSizeIn
    character (len=*), intent(in) :: NameVar
    integer, intent(out)          :: nVarLine, nPointLine
    real :: Radius

    character (len=*), parameter :: NameSub='GM_get_for_rb_trace'
    !---------------------------------------------------------------------

    if(NameVar /= 'x:y:bmin:I_I:S_I:R_I:B_I:rho:p') &
         call CON_stop(NameSub//' invalid NameVar='//NameVar)

    ! Allocate arrays
    call allocate_gm_rb(iSizeIn, jSizeIn)

    ! The RB ionosphere radius in normalized units
    Radius = (6378.+100.)/6378.  !!! could be derived from Grid_C ?
    DoExtractUnitSi = .true.
    call integrate_ray_accurate(iSizeIn, jSizeIn, RB_lat, RB_lon, Radius, &
         NameVar)

    call line_get(nVarLine, nPointLine)

    nVarLine = 4 ! We only pass line index, length, B and radial distance to RB

  end subroutine GM_get_for_rb_trace

  !==========================================================================

  subroutine GM_get_for_rb(Buffer_IIV, iSizeIn, jSizeIn, nVarIn, &
       BufferLine_VI, nVarLine, nPointLine, NameVar)

    !call stop_mpi('RAYTRACE is OFF')


    use ModGeometry,ONLY: x2
    use ModProcMH,  ONLY: iProc

    use ModMain, ONLY: Time_Simulation,TypeCoordSystem


    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUy_, RhoUz_, Bx_, By_, Bz_, p_,&
         MassFluid_I, IonFirst_, nVar

    use ModPhysics, ONLY: No2Si_V, UnitN_, UnitU_, UnitB_, UnitP_, rBody
    use ModSolarwind, ONLY: get_solar_wind_point
    use ModConst, ONLY: cProtonMass

    use CON_line_extract, ONLY: line_get, line_clean
    use CON_axes,         ONLY: transform_matrix
    use CON_planet,       ONLY: RadiusPlanet

    character (len=*), parameter :: NameSub='GM_get_for_rb'

    integer, intent(in)                                :: iSizeIn, jSizeIn, nVarIn
    real, intent(out), dimension(iSizeIn,jSizeIn,nVarIn) :: Buffer_IIV

    integer, intent(in) :: nPointLine, nVarLine
    real, intent(out)   :: BufferLine_VI(nVarLine, nPointLine)
    character (len=*), intent(in):: NameVar

    integer :: nVarExtract, nPoint, iPoint, iStartPoint
    real, allocatable :: Buffer_VI(:,:)

    logical :: DoTestTec, DoTestIdl
    logical :: DoTest, DoTestMe

    integer :: iLat,iLon,iLine, iLocBmin
    real    :: SolarWind_V(nVar), SmGm_DD(3,3), XyzBminSm_D(3)
    !--------------------------------------------------------------------------

    if(NameVar /= 'x:y:bmin:I_I:S_I:R_I:B_I:rho:p') &
         call CON_stop(NameSub//' invalid NameVar='//NameVar)

    if(iProc /= 0)then
       ! Clean and return
       call line_clean
       RETURN
    end if

    call CON_set_do_test(NameSub//'_tec', DoTestTec, DoTestMe)
    call CON_set_do_test(NameSub//'_idl', DoTestIdl, DoTestMe)
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Initialize buffer_iiv
    Buffer_IIV = 0.0

    ! Put the extracted data into BufferLine_VI
    call line_get(nVarExtract, nPoint)
    if(nPoint /= nPointLine)call stop_mpi(NameSub//': nPointLine error')
    if(nVarExtract < nVarLine)call stop_mpi(NameSub//': nVarLine error')
    allocate(Buffer_VI(0:nVarExtract, nPoint))
    call line_get(nVarExtract, nPoint, Buffer_VI, DoSort=.true.)

    ! Transformation matrix between CRCM(SM) and GM coordinates
    SmGm_DD = transform_matrix(Time_Simulation,TypeCoordSystem,'SMG')

    ! The first field line starts from iPoint = 1
    iStartPoint = 1
    do iPoint = 1, nPoint

       iLine =  Buffer_VI(0,iPoint)     ! line index
       iLat = mod(iLine-1, iSizeIn) + 1
       iLon = (iLine-1)/iSizeIn + 1

       BufferLine_VI(1,iPoint) = iLine
       BufferLine_VI(2,iPoint) = Buffer_VI(1,iPoint)                 ! Length
       BufferLine_VI(3,iPoint) = sqrt(sum(Buffer_VI(2:4,iPoint)**2)) ! |r|
       BufferLine_VI(4,iPoint) = &
            sqrt(sum(Buffer_VI(4+Bx_:4+Bz_,iPoint)**2))       ! |B|

       ! Find the location of minimum B, Bmin, and other variables at Bmin 
       ! for each field line
       if(Buffer_VI(0,min(nPoint,iPoint+1)) /= Buffer_VI(0,iPoint) &
            .or. iPoint == nPoint)then
          ! Exclude open field lines by checking the radial 
          ! distance of the last point on a field line
          if(BufferLine_VI(3,iPoint) > 1.0001*rBody*RadiusPlanet)then
             ! set line index to -1. for non-closed field line
             ! It would be much better not to send the line at all !!!
             BufferLine_VI(1,iStartPoint:iPoint) = -1
             Buffer_IIV(iLat, iLon, 1:2) = NoValue*1e6    ! x, y
             Buffer_IIV(iLat, iLon, 3) = NoValue       ! Bmin
             Buffer_IIV(iLat, iLon, 4:5) = 0.0         ! rho, p

          else
             ! For closed field lines 
             ! Location of Bmin for this field line
             iLocBmin = minloc(BufferLine_VI(4,iStartPoint:iPoint), dim=1) &
                  + iStartPoint - 1
             Buffer_IIV(iLat, iLon, 3) = BufferLine_VI(4, iLocBmin)    ! Bmin

             ! Convert location from GM to SMG coordinates
             XyzBminSm_D = matmul(SmGm_DD, Buffer_VI(2:4,iLocBmin))
             Buffer_IIV(iLat, iLon, 1) = XyzBminSm_D(1) ! x
             Buffer_IIV(iLat, iLon, 2) = XyzBminSm_D(2) ! y
             Buffer_IIV(iLat, iLon, 4) = Buffer_VI(4+Rho_,iLocBmin) &
                  /cProtonMass                          ! rho in [#/m^3]
             Buffer_IIV(iLat, iLon, 5) = Buffer_VI(4+p_,  iLocBmin)    ! p
          end if

          ! Set the start point for the next field line
          iStartPoint = iPoint +1
       end if
    end do

    deallocate(Buffer_VI)
    call line_clean

    ! Output before processing
    if(DoTest .or. DoTestTec)call write_integrated_data_tec
    if(DoTest .or. DoTestIdl)call write_integrated_data_idl


    ! Send solar wind values in the array of the extra integral
    ! This is a temporary solution. RB should use MHD_SUM_rho and MHD_SUM_p

    call get_solar_wind_point(Time_Simulation, (/x2, 0.0, 0.0/), SolarWind_V)

    Buffer_IIV(1,:,6) = SolarWind_V(Rho_)/MassFluid_I(IonFirst_) &
         *No2Si_V(UnitN_)
    Buffer_IIV(2,:,6) = SolarWind_V(RhoUx_) * No2Si_V(UnitU_)
    Buffer_IIV(3,:,6) = SolarWind_V(RhoUy_) * No2Si_V(UnitU_)
    Buffer_IIV(4,:,6) = SolarWind_V(RhoUz_) * No2Si_V(UnitU_)
    Buffer_IIV(5,:,6) = SolarWind_V(Bx_) * No2Si_V(UnitB_)
    Buffer_IIV(6,:,6) = SolarWind_V(By_) * No2Si_V(UnitB_)
    Buffer_IIV(7,:,6) = SolarWind_V(Bz_) * No2Si_V(UnitB_)
    Buffer_IIV(8,:,6) = SolarWind_V(p_)  * No2Si_V(UnitP_)


  end subroutine GM_get_for_rb

  !==========================================================================

  subroutine GM_satinit_for_rb(nSats)

    !This subroutine collects the number of satellite files for use in 
    !SWMF GM and RB coupling.

    !Module variables to use:
    use ModMain,   ONLY: DoRbSatTrace
    use ModSatelliteFile, ONLY: nSatellite

    !Subroutine Arguments:
    integer,           intent(out) :: nSats
    !--------------------------------------------------------------------------

    !If RB sat tracing is on, collect the number of satellites to trace.
    !If RB sat tracing is off, set nSats to zero.
    if (DoRbSatTrace) then
       nSats = nSatellite
    else 
       nSats = 0
    endif

  end subroutine GM_satinit_for_rb

  !==========================================================================
  subroutine GM_get_sat_for_rb(Buffer_III, Buffer_I, nSats)

    ! Subroutine to update and collect satellite locations for RB tracing

    !Modules
    use ModProcMH,        ONLY: iProc
    use ModSatelliteFile, ONLY: NameSat_I, XyzSat_DI, &
         gm_trace_sat
    use ModMain,          ONLY: UseB0, nBlock
    use ModPhysics,       ONLY: No2Si_V, UnitB_
    use ModVarIndexes,    ONLY: nVar, Bx_, By_, Bz_
    use ModB0,            ONLY: get_b0
    use ModMPI

    !Arguments
    integer, intent(in)               :: nSats
    real, intent(out)                 :: Buffer_III(4,2,nSats)
    character (len=100), intent(out)  :: Buffer_I(nSats)

    !Internal variables
    character (len=*), parameter :: NameSub='GM_get_sat_for_rb'

    real ::SatRay_D(3)

    real :: StateSat_V(0:nVar+3), B0Sat_D(3)  
    real :: Bx,By,Bz,B2
    integer :: iSat
    !--------------------------------------------------------------------------
    ! Store satellite names in Buffer_I
    Buffer_I = NameSat_I(1:nSats)

    do iSat=1, nSats
       ! Update satellite position.
       !call set_satellite_flags(iSat)
       !call get_satellite_ray(iSat, sat_RayVars)
       !
       !! Reduce values from all 
       !call MPI_reduce(sat_RayVars, sat_RayVarsSum, 5, MPI_REAL, MPI_SUM, &
       !     0, iComm, iError)
       !
       !write(*,*) 'sat_RayVars',sat_RayVars
       call GM_trace_sat(XyzSat_DI(1:3,iSat),SatRay_D)
       ! Determine magnetic field magnitude at satellite B=B0+B1
       if(UseB0)then
          call get_b0(XyzSat_DI(:,iSat), B0Sat_D)
       else
          B0Sat_D=0.00
       end if
       call get_point_data(0.0,XyzSat_DI(:,iSat),1,nBlock,1,nVar+3,StateSat_V)
       call collect_satellite_data(XyzSat_DI(:,iSat),StateSat_V)

       Bx = StateSat_V(Bx_)+B0Sat_D(1)
       By = StateSat_V(By_)+B0Sat_D(2)
       Bz = StateSat_V(Bz_)+B0Sat_D(3)

       B2 = (Bx**2.0 + By**2.0 + Bz**2.0) * (No2Si_V(UnitB_))**2.0 

       ! Store results in Buffer_III
       if (iProc == 0) then 
          Buffer_III(1:3,1,iSat)  = XyzSat_DI(1:3,iSat)
          !Buffer_III(1:3,2,iSat) = sat_RayVarsSum(1:3)
          Buffer_III(1:3,2,iSat)  = SatRay_D
          Buffer_III(4,2,iSat)    = B2
       end if
    end do

  end subroutine GM_get_sat_for_rb

end module GM_couple_rb
