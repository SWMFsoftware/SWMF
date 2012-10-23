!This code is a copyright protected software (c) 2002- University of Michigan
!^CMP FILE IM
!==========================================================================
! FOR CRCM COUPLING
!==========================================================================
subroutine GM_get_for_im_trace_crcm(iSizeIn, jSizeIn, NameVar, nVarLine, &
     nPointLine)

  ! Do ray tracing for IM. 
  ! Provide total number of points along rays 
  ! and the number of variables to pass to IM
  use ModGmImCoupling, ONLY: allocate_gm_im, IM_lat, IM_lon
  use ModRayTrace, ONLY: DoExtractUnitSi
  use ModMain, ONLY: DoMultiFluidIMCoupling, DoAnisoPressureIMCoupling
  use CON_line_extract, ONLY: line_get
  implicit none
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

  ! The RB ionosphere radius in normalized units
  Radius = (6378.+100.)/6378.  !!! could be derived from Grid_C ?
  DoExtractUnitSi = .true.

  call integrate_ray_accurate(iSizeIn, jSizeIn, IM_lat, IM_lon, Radius, &
       NameVar)

  call line_get(nVarLine, nPointLine)

  nVarLine = 4 ! We only pass line index, length, B and radial distance to RB

end subroutine GM_get_for_im_trace_crcm

!==========================================================================

subroutine GM_get_for_im_crcm(Buffer_IIV, iSizeIn, jSizeIn, nVarIn, &
     BufferLine_VI, nVarLine, nPointLine, NameVar)

  !call stop_mpi('RAYTRACE is OFF')


  use ModGeometry,ONLY: x2
  use ModProcMH,  ONLY: iProc
  use ModIoUnit, ONLY: UNITTMP_
  use ModMain, ONLY: Time_Simulation, TypeCoordSystem, &
       DoMultiFluidIMCoupling, DoAnisoPressureIMCoupling
  use ModGmImCoupling, ONLY: NoValue
  use ModVarIndexes, ONLY: &
       Rho_, RhoUx_, RhoUy_, RhoUz_, Bx_, By_, Bz_, p_, Ppar_, &
       iRho_I, iP_I, MassFluid_I, IonFirst_, IonLast_, nVar
  use ModPhysics, ONLY: No2Si_V, &
       UnitN_, UnitU_, UnitB_, UnitP_, rBody
  use ModSolarwind, ONLY: get_solar_wind_point
  use ModConst, ONLY: cProtonMass

  use CON_line_extract, ONLY: line_get, line_clean
  use CON_axes,         ONLY: transform_matrix
  use CON_planet,       ONLY: RadiusPlanet

  implicit none

  character (len=*), parameter :: NameSub='GM_get_for_im_crcm'

  integer, intent(in) :: iSizeIn, jSizeIn, nVarIn
  real,    intent(out):: Buffer_IIV(iSizeIn,jSizeIn,nVarIn)

  integer, intent(in) :: nPointLine, nVarLine
  real, intent(out)   :: BufferLine_VI(nVarLine, nPointLine)
  character (len=*), intent(in):: NameVar

  integer :: nVarExtract, nPoint, iPoint, iStartPoint
  real, allocatable :: Buffer_VI(:,:)

  logical :: DoTest, DoTestMe

  integer :: iLat,iLon,iLine,iLocBmin
  integer :: iIonSecond
  real    :: SmGm_DD(3,3), XyzBminSm_D(3)
  real    :: SolarWind_V(nVar)
  character(len=100) :: NameOut
  !--------------------------------------------------------------------------

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

  if(iProc /= 0)then
     ! Clean and return
     call line_clean
     RETURN
  end if

  call CON_set_do_test(NameSub, DoTest, DoTestMe)

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
subroutine GM_get_sat_for_im_crcm(Buffer_III, Buffer_I, nSats)

  ! Subroutine to update and collect satellite locations for IM tracing

  !Modules
  use ModProcMH,        ONLY: iProc
  use ModSatelliteFile, ONLY: NameSat_I, XyzSat_DI, gm_trace_sat
  use ModMain,          ONLY: UseB0, nBlock
  use ModPhysics,       ONLY: No2Si_V, UnitB_
  use ModVarIndexes,    ONLY: nVar, Bx_, By_, Bz_
  use ModB0,            ONLY: get_b0
  use ModMPI

  implicit none

  !Arguments
  integer, intent(in)               :: nSats
  real, intent(out)                 :: Buffer_III(4,2,nSats)
  character (len=100), intent(out)  :: Buffer_I(nSats)

  !Internal variables
  character (len=*), parameter :: NameSub='GM_get_sat_for_im'

  real ::SatRay_D(3)

  real :: StateSat_V(0:nVar+3), B0Sat_D(3)  
  real :: Bx,By,Bz,B2
  integer :: iSat
  !--------------------------------------------------------------------------
  ! Store satellite names in Buffer_I
  if (iProc == 0) then
     Buffer_I = NameSat_I(1:nSats)
  end if

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
        Buffer_III(1:3,1,iSat)   = XyzSat_DI(1:3,iSat)
        !Buffer_III(1:3,2,iSat) = sat_RayVarsSum(1:3)
        Buffer_III(1:3,2,iSat) = SatRay_D
        Buffer_III(4,2,iSat)   = B2
     end if
  end do

end subroutine GM_get_sat_for_im_crcm

!==========================================================================
subroutine GM_get_for_im_trace(nRadius, nLon, nVarLine, nPointLine, NameVar)

  ! Do ray tracing for IM/RAM_SCB. 
  ! Provide total number of points along rays 
  ! and the number of variables to pass to IM

  use ModProcMH, ONLY: iProc
  use ModMain, ONLY: Time_Simulation, TypeCoordSystem
  use ModVarIndexes, ONLY: NamePrimitiveVar, Bx_, Bz_
  use CON_line_extract, ONLY: line_get, line_clean
  use CON_coupler, ONLY: Grid_C, IM_
  use CON_axes,         ONLY: transform_matrix
  use ModGmImCoupling, ONLY: StateLine_VI
  use ModMultiFluid, ONLY: iFluid, nFluid, iUx_I, iUz_I
  use ModRaytrace, ONLY: DoExtractBGradB1

  implicit none

  integer, intent(in)           :: nRadius, nLon
  integer, intent(out)          :: nVarLine, nPointLine
  character (len=*), intent(out):: NameVar

  real, allocatable, save :: RadiusIm_I(:), LongitudeIm_I(:)

  integer :: nVarExtract, iPoint, iUx5, iUz5
  real    :: SmGm_DD(3,3)

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
  DoExtractBGradB1 = .true. !!! Maybe check NameVersionIm here ???

  ! The variables to be passed: line index, length along line, 
  ! coordinatess and primitive variables. Total is 5 + nVar.
  NameVar = 'iLine Length x y z '//NamePrimitiveVar
  if(DoExtractBGradB1) NameVar = trim(NameVar) //' bgradb1x bgradb1y bgradb1z'

  call trace_ray_equator(nRadius, nLon, RadiusIm_I, LongitudeIm_I, .true.)

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
          matmul(SmGm_DD,         StateLine_VI(3:5        ,iPoint)) ! X,Y,Z
     do iFluid = 1, nFluid
        iUx5 = iUx_I(iFluid)+5; iUz5 = iUz_I(iFluid)+5
        StateLine_VI(iUx5:iUz5,iPoint)  = &
             matmul(SmGm_DD,       StateLine_VI(iUx5:iUz5,iPoint)) ! velocity
     end do
     StateLine_VI(Bx_+5:Bz_+5,iPoint)= &
          matmul(SmGm_DD,         StateLine_VI(Bx_+5:Bz_+5,iPoint)) ! mag field

     if(DoExtractBGradB1) &
          StateLine_VI(nVarLine-2:nVarLine,iPoint)= &
          matmul(SmGm_DD, StateLine_VI(nVarLine-2:nVarLine,iPoint))! b.grad(B1)

  end do

end subroutine GM_get_for_im_trace

!==========================================================================

subroutine GM_get_for_im_line(nRadius, nLon, MapOut_DSII, &
     nVarLine, nPointLine, BufferLine_VI)

  !call stop_mpi('RAYTRACE is OFF')


  use ModProcMH,  ONLY: iProc
  use ModGmImCoupling, ONLY: StateLine_VI
  use ModRayTrace, ONLY: RayMap_DSII

  implicit none

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

subroutine GM_get_for_im(Buffer_IIV,iSizeIn,jSizeIn,nVar,NameVar)

  !call stop_mpi('RCM is OFF')

  use ModProcMH, ONLY: iProc
  use ModMain, ONLY: DoMultiFluidIMCoupling
  use ModGmImCoupling, ONLY: &
       allocate_gm_im, &
       process_integrated_data, DoTestTec, DoTestIdl, &
       write_integrated_data_tec, write_integrated_data_idl, &
       IM_lat, IM_lon, &
       MHD_SUM_vol, MHD_Xeq, MHD_Yeq, MHD_Beq, MHD_SUM_rho, MHD_SUM_p, &
       NoValue, MHD_HpRho, MHD_OpRho, MHD_Hpp, MHD_Opp

  use ModRaytrace, ONLY: RayResult_VII, RayIntegral_VII, &
       InvB_, Z0x_, Z0y_, Z0b_, RhoInvB_, pInvB_,  &
       HpRhoInvB_, OpRhoInvB_, HpPInvB_, OpPInvB_, xEnd_, CLOSEDRAY

  implicit none

  character (len=*), parameter :: NameSub='GM_get_for_im'

  integer, intent(in)                                :: iSizeIn, jSizeIn, nVar
  real, intent(out), dimension(iSizeIn,jSizeIn,nVar) :: Buffer_IIV
  character (len=*), intent(in)                      :: NameVar

  real :: Radius

  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------

  if(DoMultiFluidIMCoupling)then
     if(NameVar /= 'vol:z0x:z0y:bmin:rho:p:Hprho:Oprho:Hpp:Opp') &
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
     call integrate_ray_accurate(iSizeIn, jSizeIn, IM_lat, IM_lon, Radius, &
          'InvB,RhoInvB,pInvB,Z0x,Z0y,Z0b')
  else
     call integrate_ray_accurate(iSizeIn, jSizeIn, IM_lat, IM_lon, Radius, &
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
        MHD_SUM_rho = RayResult_VII(RhoInvB_,:,:)
        MHD_SUM_p   = RayResult_VII(pInvB_  ,:,:)
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
           MHD_SUM_rho = 0.0
           MHD_SUM_p   = 0.0
           MHD_Hprho = 0.0
           MHD_Oprho = 0.0
           MHD_HpP   = 0.0
           MHD_OpP   = 0.0
           MHD_Beq     = NoValue
        end where
     end if

  end if

  deallocate(RayIntegral_VII, RayResult_VII)

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
        ! the index is not continuous as in ModRayTrace                
        Buffer_IIV(:,:,RhoInvB_) = MHD_SUM_rho
        Buffer_IIV(:,:,pInvB_)   = MHD_SUM_p
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

  implicit none

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
subroutine GM_get_sat_for_im(Buffer_III, Buffer_I, nSats)

  ! Subroutine to update and collect satellite locations for IM tracing
  ! !!!DTW 2007

  !Modules
  use ModProcMH, ONLY: iProc, iComm
  use ModSatelliteFile, ONLY: NameSat_I, XyzSat_DI, &
       get_satellite_ray, set_satellite_flags
  use ModMPI

  implicit none

  !Arguments
  integer, intent(in)               :: nSats
  real, intent(out)                 :: Buffer_III(3,2,nSats)
  character (len=100), intent(out)  :: Buffer_I(nSats)

  !Internal variables
  character (len=*), parameter :: NameSub='GM_get_sat_for_im'

  real :: sat_RayVars(5), sat_RayVarsSum(5)

  integer :: iSat, iError
  !--------------------------------------------------------------------------
  ! Store satellite names in Buffer_I
  if (iProc == 0) then
     Buffer_I = NameSat_I(1:nSats)
  end if

  do iSat=1, nSats
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

subroutine GM_get_multi_for_im(DoMultiFluidIM)

  !This subroutine return the logical value of DoMultiFluidIMCoupling        

  use ModMain,   ONLY: DoMultiFluidIMCoupling
  implicit none

  !Subroutine Arguments:                                                      
  logical,           intent(out) :: DoMultiFluidIM
  !--------------------------------------------------------------------------  

  if (DoMultiFluidIMCoupling) then
     DoMultiFluidIM = .true.
  else
     DoMultiFluidIM = .false.
  endif

end subroutine GM_get_multi_for_im

!==========================================================================    

