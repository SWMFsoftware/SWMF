!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE IM

module GM_couple_im

  use ModMpi
  use ModNumConst, ONLY: cRadToDeg, cDegToRad
  use CON_coupler, ONLY: Grid_C, ncell_id

  use ModProcMH
  use ModMain, ONLY: n_step, &
       DoMultiFluidIMCoupling, DoAnisoPressureIMCoupling
  use ModPhysics, ONLY: No2Si_V, Si2No_V, &
       UnitP_, UnitRho_, UnitTemperature_, UnitB_, &
       Bdp, DipoleStrengthSi, rCurrents, rBody

  implicit none

  private ! except

  public:: GM_get_for_im_trace_crcm ! for IM/CRCM
  public:: GM_get_for_im_crcm       ! for IM/CRCM
  public:: GM_get_sat_for_im_crcm   ! for IM/CRCM
  public:: GM_get_for_im_trace      ! for IM/RAM
  public:: GM_get_for_im_line       ! for IM/RAM
  public:: GM_get_for_im            ! for IM/RCM
  public:: GM_satinit_for_im        ! initialize satellite
  public:: GM_get_sat_for_im        ! get satellite info
  public:: GM_put_from_im           ! from IM

  character(len=*), parameter :: NameMod='GM_couple_im'

  ! IM Grid size
  integer :: nCells_D(2), iSize,jSize

  ! Information about the IM grid ! 2D non-uniform regular grid only !!!
  real, allocatable, dimension(:) :: IM_lat, IM_lon

  integer :: i,j, i0

  real, save, dimension(:), allocatable :: &
       MHD_lat_boundary
  real, save, dimension(:,:), allocatable :: &
       MHD_SUM_vol, MHD_tmp, &
       MHD_SUM_rho, MHD_HpRho, MHD_OpRho, &
       MHD_SUM_p, MHD_HpP, MHD_OpP, &
       MHD_Beq, &
       MHD_Xeq, &
       MHD_Yeq, &
       MHD_Fluxerror
  real, parameter :: noValue=-99999.
  real :: eqB
  real :: colat,Ci,Cs,FCiCs,factor,Vol,Ri,s2,s6,s8,factor1,factor2

  integer :: iError


  logical :: DoTestTec, DoTestIdl

  ! This is for the GM-IM/RAM coupling                                      
  real, allocatable :: StateLine_VI(:,:)

contains
  !==========================================================================
  ! FOR CRCM COUPLING
  !==========================================================================
  subroutine GM_get_for_im_trace_crcm(iSizeIn, jSizeIn, NameVar, nVarLine, &
       nPointLine)

    ! Do field tracing for IM. 
    ! Provide total number of points along rays 
    ! and the number of variables to pass to IM
    use ModFieldTrace, ONLY: DoExtractUnitSi, integrate_field_from_sphere
    use CON_line_extract, ONLY: line_get

    integer, intent(in)           :: iSizeIn, jSizeIn
    character (len=*), intent(in) :: NameVar
    integer, intent(out)          :: nVarLine, nPointLine
    real :: Radius

    character (len=*), parameter :: NameSub='GM_get_for_im_trace'
    !---------------------------------------------------------------------

    if(DoMultiFluidIMCoupling)then
       if(NameVar /= 'x:y:bmin:I_I:S_I:R_I:B_I:rho:p:Hprho:Oprho:Hpp:Opp')&
            call CON_stop(NameSub//' invalid NameVar='//NameVar)
    else if(DoAnisoPressureIMCoupling)then
       if(NameVar /= 'x:y:bmin:I_I:S_I:R_I:B_I:rho:p:ppar') &
            call CON_stop(NameSub//' invalid NameVar='//NameVar)
    else
       if(NameVar /= 'x:y:bmin:I_I:S_I:R_I:B_I:rho:p') &
            call CON_stop(NameSub//' invalid NameVar='//NameVar)
    end if

    ! Allocate arrays
    call allocate_gm_im(iSizeIn, jSizeIn)

    ! The CRCM ionosphere radius in normalized units
    Radius = (6378.+100.)/6378.  !!! could be derived from Grid_C ?
    DoExtractUnitSi = .true.

    call integrate_field_from_sphere(&
         iSizeIn, jSizeIn, IM_lat, IM_lon, Radius, NameVar)

    call line_get(nVarLine, nPointLine)

    nVarLine = 4 ! We only pass line index, length, B and radial distance to RB

  end subroutine GM_get_for_im_trace_crcm

  !==========================================================================

  subroutine GM_get_for_im_crcm(Buffer_IIV, KpOut,iSizeIn, jSizeIn, nVarIn, &
       BufferLine_VI, nVarLine, nPointLine, NameVar)

    use ModGeometry,ONLY: x2
    use ModProcMH,  ONLY: iProc
    use ModIoUnit, ONLY: UNITTMP_
    use ModMain, ONLY: Time_Simulation, TypeCoordSystem
    use ModVarIndexes, ONLY: &
         Rho_, RhoUx_, RhoUy_, RhoUz_, Bx_, By_, Bz_, p_, Ppar_, &
         iRho_I, iP_I, MassFluid_I, IonFirst_, IonLast_, nVar
    use ModPhysics, ONLY: No2Si_V, &
         UnitN_, UnitU_, UnitB_, UnitP_, rBody
    use ModSolarwind, ONLY: get_solar_wind_point
    use ModConst, ONLY: cProtonMass
    use ModGroundMagPerturb, ONLY: DoCalcKp, Kp
    
    use CON_line_extract, ONLY: line_get, line_clean
    use CON_axes,         ONLY: transform_matrix
    use CON_planet,       ONLY: RadiusPlanet

    integer, intent(in) :: iSizeIn, jSizeIn, nVarIn
    real,    intent(out):: Buffer_IIV(iSizeIn,jSizeIn,nVarIn), KpOut

    integer, intent(in) :: nPointLine, nVarLine
    real, intent(out)   :: BufferLine_VI(nVarLine, nPointLine)
    character (len=*), intent(in):: NameVar

    integer :: nVarExtract, nPoint, iPoint, iStartPoint
    real, allocatable :: Buffer_VI(:,:)

    integer :: iLat,iLon,iLine,iLocBmin
    integer :: iIonSecond
    real    :: SmGm_DD(3,3), XyzBminSm_D(3)
    real    :: SolarWind_V(nVar)
    character(len=100) :: NameOut

    logical :: DoTest, DoTestMe
    character (len=*), parameter :: NameSub='GM_get_for_im_crcm'
    !--------------------------------------------------------------------------
    if(iProc /= 0) RETURN
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! If Kp is being calculated, share it.  Otherwise, share -1.
    if(DoCalcKp) then
       KpOut = Kp
    else
       KpOut = -1
    endif

    
    if(DoMultiFluidIMCoupling)then
       if(NameVar /= 'x:y:bmin:I_I:S_I:R_I:B_I:rho:p:Hprho:Oprho:Hpp:Opp')&
            call CON_stop(NameSub//' invalid NameVar='//NameVar)
    else if(DoAnisoPressureIMCoupling)then
       if(NameVar /= 'x:y:bmin:I_I:S_I:R_I:B_I:rho:p:ppar') &
            call CON_stop(NameSub//' invalid NameVar='//NameVar)
    else
       if(NameVar /= 'x:y:bmin:I_I:S_I:R_I:B_I:rho:p') &
            call CON_stop(NameSub//' invalid NameVar='//NameVar)
    end if

    ! Initialize buffer_iiv
    Buffer_IIV = 0.0

    ! Put the extracted data into BufferLine_I
    call line_get(nVarExtract, nPoint)
    if(nPoint /= nPointLine)call stop_mpi(NameSub//': nPointLine error')
    if(nVarExtract < nVarLine)call stop_mpi(NameSub//': nVarLine error')
    allocate(Buffer_VI(0:nVarExtract, nPoint))
    call line_get(nVarExtract, nPoint, Buffer_VI, DoSort=.true.)

    ! Transformation matrix between CRCM(SM) and GM coordinates
    SmGm_DD = transform_matrix(Time_Simulation,TypeCoordSystem,'SMG')

    ! For multifluid coupling
    if(DoMultiFluidIMCoupling) &
         iIonSecond = min(IonFirst_+1, IonLast_)

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
             if(DoMultiFluidIMCoupling) &
                  Buffer_IIV(iLat, iLon, 7:10) = 0.0
             if(DoAnisoPressureIMCoupling) &
                  Buffer_IIV(iLat, iLon, 7) = 0.0

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
             if(DoMultiFluidIMCoupling)then
                Buffer_IIV(iLat, iLon, 7:8) &                      ! HpRho, OpRho
                     = Buffer_VI(4+iRho_I(IonFirst_:iIonSecond), iLocBmin) 
                Buffer_IIV(iLat, iLon, 9:10) &                         ! HpP,OpP
                     = Buffer_VI(4+iP_I(IonFirst_:iIonSecond),   iLocBmin) 
             end if
             ! Ppar and HpRho share the same index in buffer_iiv for now
             ! does not work for multifluid MHD with anisotropic pressure
             if(DoAnisoPressureIMCoupling) &                           ! Ppar
                  Buffer_IIV(iLat, iLon, 7) = Buffer_VI(4+Ppar_,  iLocBmin) 
          end if

          ! Set the start point for the next field line
          iStartPoint = iPoint +1
       end if
    end do

    deallocate(Buffer_VI)
    call line_clean

    ! Send solar wind values in the array of the extra integral
    ! This is a temporary solution. IM should use MHD_SUM_rho and MHD_SUM_p

    call get_solar_wind_point(Time_Simulation, (/x2, 0.0, 0.0/), SolarWind_V)

    Buffer_IIV(1,:,6) = SolarWind_V(Rho_)/MassFluid_I(IonFirst_)*No2Si_V(UnitN_)
    Buffer_IIV(2,:,6) = SolarWind_V(RhoUx_) * No2Si_V(UnitU_)
    Buffer_IIV(3,:,6) = SolarWind_V(RhoUy_) * No2Si_V(UnitU_)
    Buffer_IIV(4,:,6) = SolarWind_V(RhoUz_) * No2Si_V(UnitU_)
    Buffer_IIV(5,:,6) = SolarWind_V(Bx_) * No2Si_V(UnitB_)
    Buffer_IIV(6,:,6) = SolarWind_V(By_) * No2Si_V(UnitB_)
    Buffer_IIV(7,:,6) = SolarWind_V(Bz_) * No2Si_V(UnitB_)
    Buffer_IIV(8,:,6) = SolarWind_V(p_)  * No2Si_V(UnitP_)



    ! test only for anisop
    if(iProc == 0 .and. DoTest .and. DoAnisoPressureIMCoupling)then
       write(NameOut,"(a,f6.1)") 'GM_get_for_IM_t_', Time_Simulation
       open(UnitTmp_,FILE=NameOut)
       write(UnitTmp_,"(a)") &
            'GM_get_for_im_crcm, Buffer_IIV, last index 1:5 and 7 '
       write(UnitTmp_,"(a)") 'Xeq Yeq Beq rho p ppar'
       do iLat =iSizeIn,1,-1
          do iLon =1,jSizeIn
             write(UnitTmp_,"(100es18.10)") Buffer_IIV(iLat,iLon,1:5), &
                  Buffer_IIV(iLat,iLon,7)
          enddo
       enddo
       close(UnitTmp_)
    end if

  end subroutine GM_get_for_im_crcm

  !==========================================================================
  subroutine GM_get_sat_for_im_crcm(Buffer_III, Name_I, nSats)

    ! Subroutine to update and collect satellite locations for IM tracing

    !Modules
    use ModProcMH,        ONLY: iProc
    use ModSatelliteFile, ONLY: NameSat_I, XyzSat_DI, gm_trace_sat
    use ModWriteLogSatFile, ONLY: collect_satellite_data
    use ModMain,          ONLY: UseB0, nBlock
    use ModPhysics,       ONLY: No2Si_V, UnitB_
    use ModVarIndexes,    ONLY: nVar, Bx_, Bz_
    use ModB0,            ONLY: get_b0
    use ModCurrent,       ONLY: get_point_data
    use ModMPI

    !Arguments
    integer, intent(in)               :: nSats
    real, intent(out)                 :: Buffer_III(4,2,nSats)
    character (len=100), intent(out)  :: Name_I(nSats)

    !Internal variables
    real ::SatRay_D(3)
    real :: StateSat_V(0:nVar+3), B0Sat_D(3)
    integer :: iSat

    character (len=*), parameter :: NameSub='GM_get_sat_for_im'
    !--------------------------------------------------------------------------
    ! Store satellite names in Buffer_I (known on all processors)
    Name_I = NameSat_I(1:nSats)

    do iSat = 1, nSats
       ! Update satellite position.
       call GM_trace_sat(XyzSat_DI(1:3,iSat), SatRay_D)

       call get_point_data( &
            0.0, XyzSat_DI(:,iSat), 1, nBlock, Bx_, Bz_, StateSat_V)
       call collect_satellite_data(XyzSat_DI(:,iSat), StateSat_V)

       ! Store results in Buffer_III from processor zero
       if (iProc == 0) then 
          Buffer_III(1:3,1,iSat) = XyzSat_DI(1:3,iSat)
          Buffer_III(1:3,2,iSat) = SatRay_D

          ! Determine total magnetic field squared at satellite: B2 = (B0+B1)^2
          B0Sat_D = 0.0
          if(UseB0) call get_b0(XyzSat_DI(:,iSat), B0Sat_D)          
          Buffer_III(4,2,iSat)   = &
               sum( (StateSat_V(Bx_:Bz_) + B0Sat_D)**2 ) * No2Si_V(UnitB_)**2
       end if
    end do

  end subroutine GM_get_sat_for_im_crcm

  !==========================================================================
  subroutine GM_get_for_im_trace(nRadius, nLon, nVarLine, nPointLine, NameVar)

    ! Do ray tracing for IM/RAM_SCB. 
    ! Provide total number of points along rays 
    ! and the number of variables to pass to IM

    use ModProcMH,     ONLY: iProc
    use ModMain,       ONLY: Time_Simulation, TypeCoordSystem
    use ModVarIndexes, ONLY: Bx_, Bz_
    use ModIO,         ONLY: NamePrimitiveVarOrig
    use CON_comp_param,   ONLY: lNameVersion
    use CON_world,        ONLY: get_comp_info
    use CON_line_extract, ONLY: line_get, line_clean
    use CON_coupler,      ONLY: Grid_C, IM_
    use CON_axes,         ONLY: transform_matrix
    use ModMultiFluid,    ONLY: iFluid, nFluid, iUx_I, iUz_I
    use ModFieldTrace,    ONLY: DoExtractBGradB1, trace_field_equator

    integer, intent(in)           :: nRadius, nLon
    integer, intent(out)          :: nVarLine, nPointLine
    character (len=*), intent(out):: NameVar

    real, allocatable, save :: RadiusIm_I(:), LongitudeIm_I(:)

    integer :: nVarExtract, iPoint, iUx5, iUz5
    real    :: SmGm_DD(3,3)

    character(len=lNameVersion) :: NameVersionIm
    character(len=*), parameter :: NameSub='GM_get_for_im_trace'
    !---------------------------------------------------------------------

    if(.not.allocated(RadiusIm_I))then
       if(   nRadius /= Grid_C(IM_) % nCoord_D(1) .or. &
            nLon    /= Grid_C(IM_) % nCoord_D(2) )then
          write(*,*)NameSub//' grid sizes do not agree nRadius,nLon,nCells=',&
               nRadius, nLon, Grid_C(IM_) % nCoord_D(1:2)
          call CON_stop(NameSub//' ERROR')
       end if
       allocate(RadiusIm_I(nRadius), LongitudeIm_I(nLon))
       RadiusIm_I    = Grid_C(IM_) % Coord1_I
       LongitudeIm_I = Grid_C(IM_) % Coord2_I
    end if

    ! Trace field lines starting from IM equatorial grid (do message pass)
    ! Include BGradB1 only if needed:
    call get_comp_info(IM_, NameVersion=NameVersionIm)
    if(NameVersionIm(1:7) == 'RAM-SCB')then
       DoExtractBGradB1 = .false. ! RAM-SCB does not need BGradB1
    else
       DoExtractBGradB1 = .true.  ! Other IM models do need BGradB1
    end if

    ! The variables to be passed: line index, length along line, 
    ! coordinatess and primitive variables. Total is 5 + nVar.
    NameVar = 'iLine Length x y z '//NamePrimitiveVarOrig
    if(DoExtractBGradB1) NameVar = trim(NameVar)//' bgradb1x bgradb1y bgradb1z'

    call trace_field_equator(nRadius, nLon, RadiusIm_I, LongitudeIm_I, .true.)

    if(iProc /= 0) RETURN

    ! Put the extracted data into BufferLine_VI
    call line_get(nVarExtract, nPointLine)
    ! The +1 is for the line index
    nVarLine = nVarExtract + 1
    allocate(StateLine_VI(nVarLine, nPointLine))
    ! Get all the line data and clean up
    call line_get(nVarExtract, nPointLine, StateLine_VI, DoSort=.true.)
    call line_clean

    ! Transformation matrix between the SM and GM coordinates
    SmGm_DD = transform_matrix(time_simulation,TypeCoordSystem,'SMG')
    do iPoint = 1, nPointLine
       StateLine_VI(3:5,iPoint)  = &
            matmul(SmGm_DD, StateLine_VI(3:5,iPoint)) ! X,Y,Z
       do iFluid = 1, nFluid
          iUx5 = iUx_I(iFluid)+5; iUz5 = iUz_I(iFluid)+5
          StateLine_VI(iUx5:iUz5,iPoint)  = &
               matmul(SmGm_DD, StateLine_VI(iUx5:iUz5,iPoint)) ! velocity
       end do
       StateLine_VI(Bx_+5:Bz_+5,iPoint)= &
            matmul(SmGm_DD, StateLine_VI(Bx_+5:Bz_+5,iPoint)) ! mag field

       ! b.grad(B1)
       if(DoExtractBGradB1) &
            StateLine_VI(nVarLine-2:nVarLine,iPoint) = &
            matmul(SmGm_DD, StateLine_VI(nVarLine-2:nVarLine,iPoint)) 

    end do

  end subroutine GM_get_for_im_trace

  !==========================================================================

  subroutine GM_get_for_im_line(nRadius, nLon, MapOut_DSII, &
       nVarLine, nPointLine, BufferLine_VI)

    !call stop_mpi('RAYTRACE is OFF')


    use ModProcMH,  ONLY: iProc
    use ModFieldTrace, ONLY: RayMap_DSII

    character (len=*), parameter :: NameSub='GM_get_for_im_line'

    integer, intent(in) :: nRadius, nLon
    real,    intent(out):: MapOut_DSII(3,2,nRadius,nLon)
    integer, intent(in) :: nPointLine, nVarLine
    real, intent(out)   :: BufferLine_VI(nVarLine, nPointLine)

    logical :: DoTest, DoTestMe
    !--------------------------------------------------------------------------

    ! This routine should be called from processor 0 only
    if(iProc /= 0) RETURN

    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Transfer extracted rays
    BufferLine_VI = StateLine_VI

    ! Transfer mapping info
    MapOut_DSII = RayMap_DSII

    deallocate(RayMap_DSII, StateLine_VI)


  end subroutine GM_get_for_im_line

  !==========================================================================

  subroutine GM_get_for_im(Buffer_IIV,KpOut,iSizeIn,jSizeIn,nVar,NameVar)

    use ModProcMH, ONLY: iProc
    use ModFieldTrace, ONLY: RayResult_VII, RayIntegral_VII, &
         InvB_, Z0x_, Z0y_, Z0b_, RhoInvB_, pInvB_,  &
         HpRhoInvB_, OpRhoInvB_, HpPInvB_, OpPInvB_, xEnd_, CLOSEDRAY, &
         integrate_field_from_sphere
    use ModGroundMagPerturb, ONLY: DoCalcKp, Kp

    integer,          intent(in) :: iSizeIn, jSizeIn, nVar
    real,             intent(out):: Buffer_IIV(iSizeIn,jSizeIn,nVar), KpOut
    character(len=*), intent(in) :: NameVar

    real :: Radius

    logical :: DoTest, DoTestMe
    character (len=*), parameter :: NameSub='GM_get_for_im'
    !--------------------------------------------------------------------------

    if(DoMultiFluidIMCoupling)then
       if(NameVar /= 'vol:z0x:z0y:bmin:Hprho:Oprho:Hpp:Opp') &
            call CON_stop(NameSub//' invalid NameVar='//NameVar)
    else
       if(NameVar /= 'vol:z0x:z0y:bmin:rho:p') &
            call CON_stop(NameSub//' invalid NameVar='//NameVar)
    end if
    call CON_set_do_test(NameSub//'_tec', DoTestTec, DoTestMe)
    call CON_set_do_test(NameSub//'_idl', DoTestIdl, DoTestMe)
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Allocate arrays
    call allocate_gm_im(iSizeIn, jSizeIn)

    ! The RCM ionosphere radius in normalized units
    Radius = (6378.+100.)/6378.  !!! could be derived from Grid_C ?
    if(.not. DoMultiFluidIMCoupling)then
       call integrate_field_from_sphere(&
            iSizeIn, jSizeIn, IM_lat, IM_lon, Radius, &
            'InvB,RhoInvB,pInvB,Z0x,Z0y,Z0b')
    else
       call integrate_field_from_sphere(&
            iSizeIn, jSizeIn, IM_lat, IM_lon, Radius, &
            'InvB,RhoInvB,pInvB,Z0x,Z0y,Z0b,HpRhoInvB,OpRhoInvB,HppInvB,OppInvB')
       ! but not pass Rhoinvb, Pinvb to IM
    end if

    if(iProc==0)then
       ! Copy RayResult into small arrays
       MHD_SUM_vol = RayResult_VII(InvB_   ,:,:)
       MHD_Xeq     = RayResult_VII(Z0x_    ,:,:)
       MHD_Yeq     = RayResult_VII(Z0y_    ,:,:)
       MHD_Beq     = RayResult_VII(Z0b_    ,:,:)
       if(.not.DoMultiFluidIMCoupling)then
          MHD_SUM_rho = RayResult_VII(RhoInvB_,:,:)
          MHD_SUM_p   = RayResult_VII(pInvB_  ,:,:)
       else
          MHD_HpRho= RayResult_VII(HpRhoInvB_,:,:)
          MHD_OpRho= RayResult_VII(OpRhoInvB_,:,:)
          MHD_HpP= RayResult_VII(HpPInvB_,:,:)
          MHD_OpP= RayResult_VII(OpPInvB_,:,:)
       end if
       ! Put impossible values if the ray is not closed
       if(.not.DoMultiFluidIMCoupling)then
          where(RayResult_VII(xEnd_,:,:) <= CLOSEDRAY)
             MHD_Xeq     = NoValue
             MHD_Yeq     = NoValue
             MHD_SUM_vol = 0.0
             MHD_SUM_rho = 0.0
             MHD_SUM_p   = 0.0
             MHD_Beq     = NoValue
          end where
       else
          where(RayResult_VII(xEnd_,:,:) <= CLOSEDRAY)
             MHD_Xeq     = NoValue
             MHD_Yeq     = NoValue
             MHD_SUM_vol = 0.0
             MHD_Hprho = 0.0
             MHD_Oprho = 0.0
             MHD_HpP   = 0.0
             MHD_OpP   = 0.0
             MHD_Beq     = NoValue
          end where
       end if

    end if

    deallocate(RayIntegral_VII, RayResult_VII)

    ! If Kp is being calculated, share it.  Otherwise, share -1.
    if(DoCalcKp) then
       KpOut = Kp
    else
       KpOut = -1
    endif
   
    
    if (iProc == 0) then
       ! Output before processing
       if(DoTestTec)call write_integrated_data_tec
       if(DoTestIdl)call write_integrated_data_idl

       ! Convert the units of the MHD_variables to SI units!!                 
       call process_integrated_data

       ! Output after processing
       if(DoTestTec)call write_integrated_data_tec
       if(DoTestIdl)call write_integrated_data_idl

       ! Put results into output buffer
       Buffer_IIV(:,:,InvB_)    = MHD_SUM_vol
       Buffer_IIV(:,:,Z0x_)     = MHD_Xeq
       Buffer_IIV(:,:,Z0y_)     = MHD_Yeq
       Buffer_IIV(:,:,Z0b_)     = MHD_Beq

       if(.not.DoMultiFluidIMCoupling)then
          Buffer_IIV(:,:,RhoInvB_) = MHD_SUM_rho
          Buffer_IIV(:,:,pInvB_)   = MHD_SUM_p
       else
          ! the index is not continuous as in ModFieldTrace                
          Buffer_IIV(:,:,HpRhoInvB_) = MHD_HpRho
          Buffer_IIV(:,:,OpRhoInvB_) = MHD_Oprho
          Buffer_IIV(:,:,HpPInvB_)   = MHD_HpP
          Buffer_IIV(:,:,OpPInvB_)   = MHD_OpP
       end if

    end if

  end subroutine GM_get_for_im

  !==========================================================================
  subroutine GM_satinit_for_im(nSats)

    !This subroutine collects the number of satellite files for use in 
    !SWMF GM and IM coupling.   !!!DTW 2007

    !Module variables to use:
    use ModMain,   ONLY: DoImSatTrace
    use ModSatelliteFile, ONLY: nSatellite

    !Subroutine Arguments:
    integer,           intent(out) :: nSats
    !--------------------------------------------------------------------------

    !If IM sat tracing is on, collect the number of satellites to trace.
    !If IM sat tracing is off, set nSats to zero.
    if (DoImSatTrace) then
       nSats = nSatellite
    else 
       nSats = 0
    endif

  end subroutine GM_satinit_for_im

  !==========================================================================
  subroutine GM_get_sat_for_im(Buffer_III, Name_I, nSats)

    ! Subroutine to update and collect satellite locations for IM tracing
    ! !!!DTW 2007

    !Modules
    use ModProcMH, ONLY: iProc, iComm
    use ModSatelliteFile, ONLY: NameSat_I, XyzSat_DI, &
         get_satellite_ray, set_satellite_flags
    use ModMPI

    !Arguments
    integer, intent(in)               :: nSats
    real, intent(out)                 :: Buffer_III(3,2,nSats)
    character (len=100), intent(out)  :: Name_I(nSats)

    !Internal variables
    character (len=*), parameter :: NameSub='GM_get_sat_for_im'

    real :: sat_RayVars(5), sat_RayVarsSum(5)

    integer :: iSat, iError
    !--------------------------------------------------------------------------
    ! Store satellite names in Buffer_I
    Name_I = NameSat_I(1:nSats)

    do iSat = 1, nSats
       ! Update satellite position.
       call set_satellite_flags(iSat)
       call get_satellite_ray(iSat, sat_RayVars)

       ! Reduce values from all 
       call MPI_reduce(sat_RayVars, sat_RayVarsSum, 5, MPI_REAL, MPI_SUM, &
            0, iComm, iError)

       ! Store results in Buffer_III
       if (iProc == 0) then 
          Buffer_III(:,1,iSat) = XyzSat_DI(:,iSat)
          Buffer_III(:,2,iSat) = sat_RayVarsSum(1:3)
       end if
    end do

  end subroutine GM_get_sat_for_im

  !==========================================================================

  subroutine GM_put_from_im(Buffer_IIV,iSizeIn,jSizeIn,nVar,NameVar)

    !call stop_mpi('RCM is OFF')

    use CON_coupler
    use CON_world,      ONLY: get_comp_info
    use CON_comp_param, ONLY: lNameVersion
    use ModImCoupling                              ! Storage for IM pressure
    use ModMain, ONLY : n_step,time_simulation
    use ModIoUnit, ONLY: UNITTMP_
    use ModProcMH, ONLY: iProc
    use ModFieldTrace, ONLY: UseAccurateTrace, DoMapEquatorRay

    character(len=80):: filename

    integer, intent(in) :: iSizeIn,jSizeIn,nVar
    real, intent(in) :: Buffer_IIV(iSizeIn,jSizeIn,nVar)
    character(len=*), intent(in) :: NameVar
    character(len=lNameVersion) :: NameVersionIm
    integer :: nCells_D(2), iError, i,j
    integer, parameter :: pres_=1, dens_=2, parpres_=3, bmin_=4, &
         Hpres_=3,Opres_=4,Hdens_=5,Odens_=6

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='GM_put_from_im'
    !--------------------------------------------------------------------------

    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoMultiFluidIMCoupling)then
       if(NameVar /= 'p:rho:Hpp:Opp:Hprho:Oprho') &
            call CON_stop(NameSub//' invalid NameVar='//NameVar)
    else if(DoAnisoPressureIMCoupling)then
       if(NameVar /= 'p:rho:ppar:bmin') &
            call CON_stop(NameSub//' invalid NameVar='//NameVar)
    else
       if(NameVar /= 'p:rho') &
            call CON_stop(NameSub//' invalid NameVar='//NameVar)
    end if

    nCells_D=ncell_id(IM_)
    if( iSizeIn /= nCells_D(1) .or. jSizeIn /= nCells_D(2) ) then

       write(*,*)NameSub//' grid sizes do not agree iSize,jSize,nCells=',&
            iSizeIn,jSizeIn,nCells_D(1:2)
       call CON_stop(NameSub//' SWMF_ERROR')
    end if

    if(.not.allocated(IM_lat))then
       ! Allocate IM_lat, IM_lon, IM_p, IM_dens
       call im_pressure_init(iSizeIn, jSizeIn)
       ! Set up IM ionospheric grid and store.
       ! Latitude specification is module specific, we must set it up
       ! according to the IM module we have selected.
       ! Determine version of IM:
       call get_comp_info(IM_, NameVersion=NameVersionIm)
       if(NameVersionIm(1:3) == 'RAM')then
          ! HEIDI and RAM-SCB have similar equatorial grids.
          IM_lat = Grid_C(IM_) % Coord1_I

          ! Coupling with RAM/HEIDI grid requires accurate raytrace
          ! that stops at the equator
          UseAccurateTrace= .true.
          DoMapEquatorRay = .true.
       else
          ! RCM uses a CoLat based grid, information is stored in 
          ! module grid information.
          IM_lat = (cHalfPi - Grid_C(IM_) % Coord1_I) * cRadToDeg
       end if
       IM_lon = Grid_C(IM_)% Coord2_I * cRadToDeg
    end if

    ! Store IM variable for internal use
    IM_p    = Buffer_IIV(:,:,pres_)
    IM_dens = Buffer_IIV(:,:,dens_)
    iNewPIm  = iNewPIm + 1

    ! for multifluid                                               
    if(DoMultiFluidIMCoupling)then
       IM_Hpp = Buffer_IIV(:,:,Hpres_)
       IM_Opp = Buffer_IIV(:,:,Opres_)
       IM_Hpdens = Buffer_IIV(:,:,Hdens_)
       IM_Opdens = Buffer_IIV(:,:,Odens_)
    endif

    ! for anisotropic pressure                                                                                       
    if(DoAnisoPressureIMCoupling)then
       IM_ppar = Buffer_IIV(:,:,parpres_)
       IM_bmin = Buffer_IIV(:,:,bmin_)
    end if

    if(DoTest)call write_IMvars_tec  ! TecPlot output
    if(DoTest)call write_IMvars_idl  ! IDL     output

  contains

    !============================================================================
    subroutine write_IMvars_tec
      integer :: j2
      real :: lonShift
      !-------------------------------------------------------------------------
      if(iProc /= 0)RETURN

      !write values to plot file
      write(filename,'(a,i6.6,a)')"IMp_n=",n_step,".dat"
      OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown')
      write(UNITTMP_,'(a)') 'TITLE="Raytrace Values"'
      if(DoMultiFluidIMCoupling)then
         write(UNITTMP_,'(a)') 'VARIABLES="J", "I", "Lon", "Lat",&                
              &"IM pressure", "IM density", &                                   
              &"IM Hp pressure", "IM Hp density", &                             
              &"IM Op pressure", "IM Op density"'
      else
         write(UNITTMP_,'(a)') 'VARIABLES="J", "I", "Lon", "Lat",&                
              &"IM pressure", "IM density"'
      end if
      write(UNITTMP_,'(a,i4,a,i4,a)') &
           'ZONE T="IM Pressure", I=',jSizeIn+1,', J=',iSizeIn,', K=1, F=POINT'
      do i=1,iSizeIn
         do j2=1,jSizeIn+1
            j=j2; if(j2==jSizeIn+1) j=1
            lonShift=0.; if(j2==jSizeIn+1) lonShift=360.
            if(DoMultiFluidIMCoupling)then
               write(UNITTMP_,'(2i4,8G14.6)') j2,i,IM_lon(j)+lonShift,IM_lat(i), &
                    IM_p(i,j),IM_dens(i,j), &
                    IM_Hpp(i,j),IM_Hpdens(i,j),IM_Opp(i,j), IM_Opdens(i,j)
            else
               write(UNITTMP_,'(2i4,4G14.6)') j2,i,IM_lon(j)+lonShift,IM_lat(i), &
                    IM_p(i,j),IM_dens(i,j)
            endif
         end do
      end do
      CLOSE(UNITTMP_)

    end subroutine write_IMvars_tec

    !==========================================================================
    subroutine write_IMvars_idl
      if(iProc /= 0)RETURN

      !write values to plot file
      write(filename,'(a,i6.6,a)')"IMp_n=",n_step,".out"
      OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown', &
           iostat =iError)
      if (iError /= 0) call CON_stop("Can not open file "//filename)
      write(UNITTMP_,'(a79)')            'IM pressure'
      write(UNITTMP_,'(i7,1pe13.5,3i3)') n_step,time_simulation,2,0,2
      write(UNITTMP_,'(3i4)')            jSizeIn,iSizeIn
      if(DoMultiFluidIMCoupling)then
         write(UNITTMP_,'(a79)')'Lon Lat p rho Hpp Hprho Opp Oprho'
      else
         write(UNITTMP_,'(a79)')'Lon Lat p rho'
      endif
      do i=iSizeIn,1,-1
         do j=1,jSizeIn
            if(DoMultiFluidIMCoupling)then
               write(UNITTMP_,'(100(1pe18.10))') &
                    IM_lon(j),IM_lat(i),IM_p(i,j),IM_dens(i,j), &
                    IM_Hpp(i,j), IM_Hpdens(i,j), IM_Opp(i,j), IM_Opdens(i,j)
            else
               write(UNITTMP_,'(100(1pe18.10))') &
                    IM_lon(j),IM_lat(i),IM_p(i,j),IM_dens(i,j)
            endif
         end do
      end do
      CLOSE(UNITTMP_)

    end subroutine write_IMvars_idl

  end subroutine GM_put_from_im

  !============================================================================

  subroutine allocate_gm_im(iSizeIn,jSizeIn)
    use CON_comp_param, ONLY: IM_

    integer, intent(in) :: iSizeIn, jSizeIn
    character(len=*), parameter:: NameSub=NameMod//'::allocate_gm_im'

    if(allocated(MHD_lat_boundary)) deallocate(MHD_lat_boundary)
    if(allocated(MHD_SUM_vol))      deallocate(MHD_SUM_vol)
    if(allocated(MHD_tmp))          deallocate(MHD_tmp)
    if(allocated(MHD_Beq))          deallocate(MHD_Beq)
    if(allocated(MHD_Xeq))          deallocate(MHD_Xeq)
    if(allocated(MHD_Yeq))          deallocate(MHD_Yeq)
    if(allocated(MHD_Fluxerror))    deallocate(MHD_Fluxerror)

    if(DoMultiFluidIMCoupling)then
       if(allocated(MHD_HpRho))        deallocate(MHD_HpRho)
       if(allocated(MHD_HpP))          deallocate(MHD_HpP)
       if(allocated(MHD_OpRho))        deallocate(MHD_OpRho)
       if(allocated(MHD_OpP))          deallocate(MHD_OpP)
    else
       if(allocated(MHD_SUM_rho))      deallocate(MHD_SUM_rho)
       if(allocated(MHD_SUM_p))        deallocate(MHD_SUM_p)
    end if

    iSize = iSizeIn
    jSize = jSizeIn

    if(.not.allocated(IM_lat))then
       nCells_D=ncell_id(IM_)
       if(  iSize /= nCells_D(1) .or. &
            jSize /= nCells_D(2) ) then
          write(*,*)NameSub//' grid sizes do not agree iSize,jSize,nCells=',&
               iSize,jSize, nCells_D(1:2)
          call CON_stop(NameSub//' ERROR')
       end if
       allocate(IM_lat(iSize), IM_lon(jSize))
       ! Convert colat, lon to lat-lon in degrees
       IM_lat = 90.0 - Grid_C(IM_) % Coord1_I * cRadToDeg
       IM_lon =        Grid_C(IM_) % Coord2_I * cRadToDeg
    end if

    ! Arrays needed for the field line integrals
    allocate( MHD_SUM_vol(isize,jsize))
    MHD_SUM_vol = 0.

    if(DoMultiFluidIMCoupling)then
       allocate( MHD_Hprho(isize,jsize), MHD_Oprho(isize,jsize), &
            MHD_Hpp(isize,jsize), MHD_Opp(isize,jsize))
       MHD_Hprho = 0.
       MHD_Oprho = 0.
       MHD_Hpp = 0.
       MHD_Opp = 0.
    else
       allocate(MHD_SUM_rho(isize,jsize), MHD_SUM_p(isize,jsize))
       MHD_SUM_rho = 0.
       MHD_SUM_p = 0.
    end if

    allocate( MHD_Beq(isize,jsize))
    MHD_Beq = 0.

    allocate( MHD_Xeq(isize,jsize))
    MHD_Xeq = 0.

    allocate( MHD_Yeq(isize,jsize))
    MHD_Yeq = 0.

    allocate( MHD_tmp(isize,jsize))
    MHD_tmp = 0.

    allocate( MHD_Fluxerror(isize,jsize))
    MHD_Fluxerror = 0.

    allocate( MHD_lat_boundary(jsize))
    MHD_lat_boundary = 0

  end subroutine allocate_gm_im

  !============================================================================
  subroutine write_integrated_data_tec
    use ModIoUnit, ONLY: UNITTMP_
    CHARACTER (LEN=80) :: filename
    integer :: j2, nCall=0
    real :: tmpT, tmpV1,tmpV2, lonShift,tmpHpT,tmpOpT
    !-------------------------------------------------------------------------

    nCall=nCall+1

    !write values to plot file
    write(filename,'(a,i6.6,a,i4.4,a)')"rayValues_n=",n_step,"_",nCall,".dat"

    OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown')
    write(UNITTMP_,'(a)') 'TITLE="Raytrace Values"'
    if(DoMultiFluidIMCoupling)then
       write(UNITTMP_,'(a)') 'VARIABLES="J", "I", "Lon", "Lat", "Lat Boundary (I)"', &
            ', "Xeq", "Yeq"', &
            ', "Volume", "Volume**(-2/3)"', &
            ', "MHD `r", "MHD p", "MHD T"', &
            ', "MHD Hp`r", "MHD Hp p", "MHD Hp T"', &
            ', "MHD Op`r", "MHD Op p", "MHD Op T", "Beq"', &
            ', "FluxError"'
    else
       write(UNITTMP_,'(a)') 'VARIABLES="J", "I", "Lon", "Lat", "Lat Boundary (I)"', &
            ', "Xeq", "Yeq"', &
            ', "Volume", "Volume**(-2/3)"', &
            ', "MHD `r", "MHD p", "MHD T", "Beq"', &
            ', "FluxError"'
    end if
    write(UNITTMP_,'(a,i3.3,a,i4,a,i4,a)') &
         'ZONE T="PE=',iProc,'", I=',jsize+1,', J=',isize,', K=1, F=POINT'
    do i=1,isize
       do j2=1,jsize+1
          j=j2; if(j2==jsize+1) j=1
          lonShift=0.; if(j2==jsize+1) lonShift=360.
          if(DoMultiFluidIMCoupling)then
             tmpHpT=-1.; if(MHD_Hprho(i,j)>0.) &
                  tmpHpT = ((MHD_Hpp(i,j)*Si2No_V(UnitP_))/(MHD_Hprho(i,j)*Si2No_V(UnitRho_))) &
                  * No2Si_V(UnitTemperature_)
             tmpOpT=-1.; if(MHD_Oprho(i,j)>0.) &
                  tmpOpT = ((MHD_Opp(i,j)*Si2No_V(UnitP_))/(MHD_Oprho(i,j)*Si2No_V(UnitRho_))) &
                  * No2Si_V(UnitTemperature_)
          else
             tmpT=-1.; if(MHD_SUM_rho(i,j)>0.) &
                  tmpT = ((MHD_SUM_p(i,j)*Si2No_V(UnitP_))/(MHD_SUM_rho(i,j)*Si2No_V(UnitRho_))) &
                  * No2Si_V(UnitTemperature_)
          end if
          tmpV1=0.; if(MHD_SUM_vol(i,j)>0.) &
               tmpV1 = (MHD_SUM_vol(i,j)/1.e9)
          tmpV2=0.; if(MHD_SUM_vol(i,j)>0.) &
               tmpV2 = (MHD_SUM_vol(i,j)/1.e9)**(-2./3.)
          if(DoMultiFluidIMCoupling)then
             write(UNITTMP_,'(2i4,18G14.6)') j2,i,IM_lon(j)+lonShift,IM_lat(i),&
                  MHD_lat_boundary(j), &
                  MHD_Xeq(i,j),MHD_Yeq(i,j), &
                  tmpV1,tmpV2, &
                  MHD_Hprho(i,j), MHD_Hpp(i,j), tmpHpT, &
                  MHD_Oprho(i,j), MHD_Opp(i,j), tmpOpT,MHD_Beq(i,j), &
                  MHD_Fluxerror(i,j)
          else
             write(UNITTMP_,'(2i4,12G14.6)') j2,i,IM_lon(j)+lonShift,IM_lat(i),&
                  MHD_lat_boundary(j), &
                  MHD_Xeq(i,j),MHD_Yeq(i,j), &
                  tmpV1,tmpV2, &
                  MHD_SUM_rho(i,j),MHD_SUM_p(i,j),tmpT,MHD_Beq(i,j), &
                  MHD_Fluxerror(i,j)
          end if
       end do
    end do
    CLOSE(UNITTMP_)

  end subroutine write_integrated_data_tec

  !============================================================================
  subroutine write_integrated_data_idl

    use ModIoUnit, ONLY: UNITTMP_
    use ModMain,   ONLY: time_simulation
    CHARACTER (LEN=100) :: filename
    integer :: nCall = 0
    !-------------------------------------------------------------------------

    !write values to plot file
    nCall = nCall+1
    write(filename,'(a,i6.6,a,i4.4,a)')"rayValues_n=",n_step,"_",nCall,".out"

    OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown', &
         iostat =iError)
    if (iError /= 0) call CON_stop("Can not open raytrace File "//filename)
    write(UNITTMP_,'(a79)')            'Raytrace Values_var22'
    write(UNITTMP_,'(i7,1pe13.5,3i3)') n_step,time_simulation,2,1,7
    write(UNITTMP_,'(3i4)')            jSize+1,iSize
    write(UNITTMP_,'(100(1pe13.5))')   0.0
    write(UNITTMP_,'(a79)') 'Lon Lat Xeq Yeq vol rho p Beq FluxError nothing'
    do i=isize,1,-1
       do j=1,jsize
          write(UNITTMP_,'(100(1pe18.10))') &
               IM_lon(j),       &
               IM_lat(i),       &
               MHD_Xeq(i,j),     &
               MHD_Yeq(i,j),     &
               MHD_SUM_vol(i,j), &
               MHD_SUM_rho(i,j), &
               MHD_SUM_p(i,j),   &
               MHD_Beq(i,j),     &
               MHD_FluxError(i,j)
       end do
       write(UNITTMP_,'(100(1pe18.10))') &
            IM_lon(1)+360.0, &
            IM_lat(i),       &
            MHD_Xeq(i,1),     &
            MHD_Yeq(i,1),     &
            MHD_SUM_vol(i,1), &
            MHD_SUM_rho(i,1), &
            MHD_SUM_p(i,1),   &
            MHD_Beq(i,1),     &
            MHD_FluxError(i,1)               
    end do
    CLOSE(UNITTMP_)

  end subroutine write_integrated_data_idl

  !===========================================================================
  subroutine process_integrated_data

    integer :: iLoc_I(1),iEquator,iNorthPole

    if(DoMultiFluidIMCoupling)then
       where(MHD_SUM_vol>0.)
          MHD_Hprho = MHD_Hprho/MHD_SUM_vol
          MHD_Oprho = MHD_Oprho/MHD_SUM_vol
          MHD_Hpp   = MHD_Hpp/MHD_SUM_vol
          MHD_Opp   = MHD_Opp/MHD_SUM_vol
       end where
    else
       where(MHD_SUM_vol>0.)
          MHD_SUM_p   = MHD_SUM_p/MHD_SUM_vol
          MHD_SUM_rho = MHD_SUM_rho/MHD_SUM_vol
       end where
    end if

    !Set volume floor
    MHD_SUM_vol = max(1.E-8,MHD_SUM_vol)

    ! If the field-line tracer returned a good value, we may need to
    ! add contribution to V from the part of the field line that
    ! goes inside the body. If the field-line tracer did not
    ! return a good value, we will compute total V assuming dipole.
    Ri=(6378.+100.)/6378.
    Factor = 2. * (Ri**4) / abs(Bdp)
    do i=1,isize
       Colat = (90.0 - IM_lat(i))*cDegToRad
       s2    = sin(colat)**2
       s6    = s2**3
       s8    = s6*s2

       Ci=abs(cos(colat))

       if( s2 < Ri/Rbody )then
          !Fieldline goes beyond Rbody, add piece of fieldline volume
          Cs=sqrt(1.-(Rbody/Ri)*s2)
          FCiCs = (Ci-Cs) - (Ci**3-Cs**3) + (3./5.)*(Ci**5-Cs**5) &
               - (1./7.)*(Ci**7-Cs**7)
!!!CHANGE!!!
          if (s8 /= 0.0) then
             Vol = Factor*FCiCs/s8
          else
             Vol = 0.0
          endif
          where(MHD_SUM_vol(i,:)>1.1E-8)
             MHD_SUM_vol(i,:)=MHD_SUM_vol(i,:) + Vol
          end where
       end if

       if( s2 > Ri/Rcurrents )then
          !Fieldline stays inside of Rcurrents, recompute some values

          !Compute the full analytic volume
          FCiCs = Ci - Ci**3 + (3./5.)*Ci**5 - (1./7.)*Ci**7
!!!CHANGE!!!
          if (s8 /= 0.0) then
             Vol = factor*FCiCs/s8
          else
             Vol = 0.0
          endif

          !Compute equatorial B value for dipole at this latitude
          eqB = abs(Bdp)*s6/Ri**3

          if( s2 > Ri/Rbody )then
             ! Fieldline stays inside of Rbody

             ! Recompute exact volume
             MHD_SUM_vol(i,:)=Vol

             ! Fix the grid inside Rbody
             MHD_Xeq(i,:) = (Ri/s2)*cos(IM_lon(:)*cDegToRad)
             MHD_Yeq(i,:) = (Ri/s2)*sin(IM_lon(:)*cDegToRad)

             ! Fix the equatorial B value
             MHD_Beq(i,:) = eqB

          else
             ! Fieldline stays inside of Rcurrents but outside Rbody
             ! Weight analytic formula proportional to the distance 
             ! of the field line from rBody within the rBody-rCurrents range

             Factor1= (Ri/Rbody - s2)/ (Ri/Rbody - Ri/Rcurrents)
             Factor2= 1.0 - Factor1

             ! Check if numerical volume exists
             where(MHD_SUM_vol(i,:)>1.1E-8)
                ! Blend numerical volume with exact volume
                MHD_SUM_vol(i,:) = Factor1*MHD_SUM_vol(i,:) + Factor2*Vol
             end where

          end if
       end if
    end do

    ! find index for latitude closest to equator
    iLoc_I   = minloc(abs(IM_lat))
    iEquator = iLoc_I(1)

    ! find index for latitude closest to north pole
    iLoc_I = maxloc(IM_lat)
    iNorthPole = iLoc_I(1)

    iEquator = iEquator + sign(1, iNorthPole - iEquator) 
    !set open fieldline values
    do j=1,jsize

       ! Initialize the index of the last closed field line to be at
       ! the highest latitude in the RCM grid
       i0 = iNorthPole

       do i = iEquator, iNorthPole, sign(1, iNorthPole - iEquator)

          if(MHD_SUM_vol(i,j) < 1.1E-8 .or. &
               abs(MHD_Xeq(i,j)) > 200.0 .or. &
               abs(MHD_Yeq(i,j)) > 200.0 ) then
             i0=i+1
             exit
          end if
       end do

       ! Save the index of the "last" closed field line into MHD_lat_boundary
       MHD_lat_boundary(j)=i0
       ! write(*,*) "finding closed field-line ",j,i0,MHD_lat_boundary(j)

    end do

    ! Set impossible values for open fieldlines
    ! except for Xeq and Yeq where the last closed values are used
    ! which is useful when the equatorial grid is plotted
    do j=1,jsize
       i = int(MHD_lat_boundary(j))
       MHD_Beq    (1:i-1,j) = -1.
       MHD_Xeq    (1:i-1,j) = MHD_Xeq(i,j)
       MHD_Yeq    (1:i-1,j) = MHD_Yeq(i,j)
       MHD_SUM_vol(1:i-1,j) = -1.       
       if(DoMultiFluidIMCoupling)then
          MHD_HpRho(1:i-1,j) = -1.
          MHD_OpRho(1:i-1,j) = -1.
          MHD_HpP  (1:i-1,j) = -1.
          MHD_OpP  (1:i-1,j) = -1.
       else
          MHD_SUM_rho(1:i-1,j) = -1.
          MHD_SUM_p  (1:i-1,j) = -1.
       endif
    end do

    ! Dimensionalize values
    ! Note: dimensions of "MHD_sum_vol" is Distance/Magnetic field.
    ! The distance unit is planetary radius in GM, and the same is
    ! used in RCM, so only the magnetic unit is converted to SI units.
    ! Similarly Xeq and Yeq (equatorial crossing coords) remain in 
    ! normalized units.
    if(DoMultiFluidIMCoupling)then
       where(MHD_SUM_vol > 0.)
          MHD_SUM_vol = MHD_SUM_vol / No2Si_V(UnitB_)
          MHD_HpRho = MHD_HpRho * No2Si_V(UnitRho_)
          MHD_OpRho = MHD_OpRho * No2Si_V(UnitRho_)
          MHD_HpP   = MHD_HpP * No2Si_V(UnitP_)
          MHD_OpP   = MHD_OpP * No2Si_V(UnitP_)
       elsewhere
          MHD_SUM_vol = -1.
          MHD_HpRho = -1.
          MHD_OpRho = -1.
          MHD_HpP   = -1.
          MHD_OpP   = -1.
       end where
    else
       where(MHD_SUM_vol > 0.)
          MHD_SUM_vol = MHD_SUM_vol / No2Si_V(UnitB_)
          MHD_SUM_rho = MHD_SUM_rho * No2Si_V(UnitRho_)
          MHD_SUM_p   = MHD_SUM_p   * No2Si_V(UnitP_)
       elsewhere
          MHD_SUM_vol = -1.
          MHD_SUM_rho = -1.
          MHD_SUM_p   = -1.
       end where
    end if

    where(MHD_Beq > 0.)
       MHD_Beq = MHD_Beq * No2Si_V(UnitB_)
    elsewhere
       MHD_Beq = -1.
    end where

    if(DoTestTec .or. DoTestIdl)then
       do j=1,jsize-1; do i=1,isize-1
          if ( MHD_SUM_vol(i  ,j) > 0.0 .AND. &
               MHD_SUM_vol(i+1,j) > 0.0 .AND. &
               MHD_SUM_vol(i,j+1) > 0.0 .AND. &
               MHD_SUM_vol(i+1,j+1)> 0.0) THEN
             MHD_Fluxerror(i,j) = &
                  0.25E+9*(MHD_Beq(i,j  ) + MHD_Beq(i+1,j  ) + &
                  MHD_Beq(i,j+1) + MHD_Beq(i+1,j+1)) * &
                  0.5*(ABS((MHD_Xeq(i,j+1)-MHD_Xeq(i,j))* &
                  (MHD_Yeq(i+1,j)-MHD_Yeq(i,j))  &
                  - (MHD_Yeq(i,j+1)-MHD_Yeq(i,j)) &
                  * (MHD_Xeq(i+1,j)-MHD_Xeq(i,j))) &
                  + ABS((MHD_Xeq(i+1,j)-MHD_Xeq(i+1,j+1)) &
                  *(MHD_Yeq(i,j+1)-MHD_Yeq(i+1,j+1)) &
                  -(MHD_Yeq(i+1,j)-MHD_Yeq(i+1,j+1)) &
                  *(MHD_Xeq(i,j+1)-MHD_Xeq(i+1,j+1))))/&
                  (ABS(DipoleStrengthSi)*(SIN(Im_lat(i)*cDegToRad)**2 &
                  -SIN(Im_lat(i+1)*cDegToRad)**2)* &
                  (IM_lon(j+1)-IM_lon(j))*cDegToRad )- 1.0
          ELSE
             MHD_Fluxerror (i,j) = 0.0
          END IF
       end do; end do
       MHD_fluxerror(:,jsize) = MHD_Fluxerror (:,1)
       MHD_fluxerror(isize,:) = MHD_FLuxerror(isize-1,:)
    end if

  end subroutine process_integrated_data

end module GM_couple_im
