!^CFG COPYRIGHT UM
!^CMP FILE RB

!==========================================================================
subroutine GM_get_for_rb_trace(iSizeIn, jSizeIn, NameVar, nVarLine, nPointLine)

  ! Do ray tracing for RB. 
  ! Provide total number of points along rays 
  ! and the number of variables to pass to RB
  use ModGmRbCoupling, ONLY: allocate_gm_rb, RB_lat, RB_lon
  use ModRayTrace, ONLY: DoExtractUnitSi

  use CON_line_extract, ONLY: line_get
  implicit none
  integer, intent(in)           :: iSizeIn, jSizeIn
  character (len=*), intent(in) :: NameVar
  integer, intent(out)          :: nVarLine, nPointLine
  real :: Radius

  character (len=*), parameter :: NameSub='GM_get_for_rb_trace'
  !---------------------------------------------------------------------

  if(NameVar /= 'Z0x:Z0y:Z0b:I_I:S_I:R_I:B_I:IMF') &
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

  !call stop_mpi('RAYTRACE is OFF') !^CFG UNCOMMENT IF NOT RAYTRACE
  !^CFG IF RAYTRACE BEGIN

  use ModGeometry,ONLY: x2
  use ModProcMH,  ONLY: iProc

  use ModMain, ONLY: Time_Simulation

  use ModGmRbCoupling, ONLY: &
       RB_lat, RB_lon, &
       write_integrated_data_tec, write_integrated_data_idl, &
       MHD_SUM_vol, MHD_Xeq, MHD_Yeq, MHD_Beq, MHD_SUM_rho, MHD_SUM_p, NoValue

  use ModRaytrace, ONLY: RayResult_VII, RayIntegral_VII, &
       InvB_, Z0x_, Z0y_, Z0b_, RhoInvB_, pInvB_, xEnd_, CLOSEDRAY

  use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUy_, RhoUz_, Bx_, By_, Bz_, p_,&
                           MassFluid_I, IonFirst_, nVar

  use ModPhysics, ONLY: No2Si_V, UnitN_, UnitU_, UnitB_, UnitP_,UnitRho_
  use ModSolarwind, ONLY: get_solar_wind_point

  use CON_line_extract, ONLY: line_get, line_clean

  implicit none

  character (len=*), parameter :: NameSub='GM_get_for_rb'

  integer, intent(in)                                :: iSizeIn, jSizeIn, nVarIn
  real, intent(out), dimension(iSizeIn,jSizeIn,nVarIn) :: Buffer_IIV

  integer, intent(in) :: nPointLine, nVarLine
  real, intent(out)   :: BufferLine_VI(nVarLine, nPointLine)
  character (len=*), intent(in):: NameVar

  integer :: nVarExtract, nPoint, iPoint
  real, allocatable :: Buffer_VI(:,:)
  real :: Rho, Ux, Uy, Uz, Bx, By, Bz, p

  logical :: DoTestTec, DoTestIdl
  logical :: DoTest, DoTestMe

  integer :: iLat,iLon,iLine
  real    :: SolarWind_V(nVar)
  !--------------------------------------------------------------------------

  if(NameVar /= 'Z0x:Z0y:Z0b:I_I:S_I:R_I:B_I:IMF') &
       call CON_stop(NameSub//' invalid NameVar='//NameVar)

  if(iProc /= 0)then
     ! Clean and return
     deallocate(RayIntegral_VII, RayResult_VII)
     call line_clean
     RETURN
  end if

  call CON_set_do_test(NameSub//'_tec', DoTestTec, DoTestMe)
  call CON_set_do_test(NameSub//'_idl', DoTestIdl, DoTestMe)
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  ! Copy RayResult into small arrays used in old algorithm
  MHD_SUM_vol = RayResult_VII(InvB_   ,:,:)
  MHD_Xeq     = RayResult_VII(Z0x_    ,:,:)
  MHD_Yeq     = RayResult_VII(Z0y_    ,:,:)
  MHD_Beq     = RayResult_VII(Z0b_    ,:,:)
  MHD_SUM_rho = RayResult_VII(RhoInvB_,:,:)
  MHD_SUM_p   = RayResult_VII(pInvB_  ,:,:)

  ! Put impossible values if ray was not found for a lat-lon grid cell
  where(RayResult_VII(xEnd_,:,:) <= CLOSEDRAY)
     MHD_Xeq     = NoValue
     MHD_Yeq     = NoValue
     MHD_SUM_vol = 0.0
     MHD_SUM_rho = 0.0
     MHD_SUM_p   = 0.0
     MHD_Beq     = NoValue
  end where

  ! Put the extracted data into BufferLine_VI

  call line_get(nVarExtract, nPoint)
  if(nPoint /= nPointLine)call stop_mpi(NameSub//': nPointLine error')
  if(nVarExtract < nVarLine)call stop_mpi(NameSub//': nVarLine error')
  allocate(Buffer_VI(0:nVarExtract, nPoint))
  call line_get(nVarExtract, nPoint, Buffer_VI, DoSort=.true.)
  
  do iPoint = 1, nPoint

     iLine =  Buffer_VI(0,iPoint)     ! line index
     iLat = mod(iLine-1, iSizeIn) + 1
     iLon = (iLine-1)/iSizeIn + 1

     ! exclude open field lines by setting impossible line index
     if(MHD_Xeq(iLat, iLon) == NoValue)then
        iLine = -1
     endif

     BufferLine_VI(1,iPoint) = iLine
     BufferLine_VI(2,iPoint) = Buffer_VI(1,iPoint)                 ! Length
     BufferLine_VI(3,iPoint) = sqrt(sum(Buffer_VI(2:4,iPoint)**2)) ! |r|
     BufferLine_VI(4,iPoint) = &
          sqrt(sum(Buffer_VI(4+Bx_:4+Bz_,iPoint)**2))       ! |B|
  end do
  
  deallocate(RayIntegral_VII, RayResult_VII, Buffer_VI)
  call line_clean

  ! Output before processing
  if(DoTest .or. DoTestTec)call write_integrated_data_tec
  if(DoTest .or. DoTestIdl)call write_integrated_data_idl

  ! Put results into output buffer
  Buffer_IIV(:,:,1)  = MHD_Xeq
  Buffer_IIV(:,:,2)  = MHD_Yeq
  Buffer_IIV(:,:,3)  = MHD_Beq * No2Si_V(UnitB_)
  
  ! Send solar wind values in the array of the extra integral
  ! This is a temporary solution. RB should use MHD_SUM_rho and MHD_SUM_p

  call get_solar_wind_point(Time_Simulation, (/x2, 0.0, 0.0/), SolarWind_V)

  Buffer_IIV(1,:,4) = SolarWind_V(Rho_)/MassFluid_I(IonFirst_)*No2Si_V(UnitN_)
  Buffer_IIV(2,:,4) = SolarWind_V(RhoUx_) * No2Si_V(UnitU_)
  Buffer_IIV(3,:,4) = SolarWind_V(RhoUy_) * No2Si_V(UnitU_)
  Buffer_IIV(4,:,4) = SolarWind_V(RhoUz_) * No2Si_V(UnitU_)
  Buffer_IIV(5,:,4) = SolarWind_V(Bx_) * No2Si_V(UnitB_)
  Buffer_IIV(6,:,4) = SolarWind_V(By_) * No2Si_V(UnitB_)
  Buffer_IIV(7,:,4) = SolarWind_V(Bz_) * No2Si_V(UnitB_)
  Buffer_IIV(8,:,4) = SolarWind_V(p_)  * No2Si_V(UnitP_)

  !^CFG END RAYTRACE
end subroutine GM_get_for_rb


!==========================================================================
subroutine GM_satinit_for_rb(nSats)

  !This subroutine collects the number of satellite files for use in 
  !SWMF GM and RB coupling.

  !Module variables to use:
  use ModMain,   ONLY: DoRbSatTrace
  use ModSatelliteFile, ONLY: nSatellite
  use ModProcMH, ONLY: iProc

  implicit none

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
  use ModProcMH, ONLY: iProc, iComm
  use ModSatelliteFile, ONLY: Satellite_name, Xsatellite, &
       get_satellite_ray, set_satellite_flags, gm_trace_sat
  use ModMain,   ONLY: UseB0, nBlock
  use ModMPI
  use ModPhysics, ONLY: No2Si_V, UnitB_
  use ModVarIndexes, ONLY : nVar,Bx_,By_,Bz_
  implicit none
  
  !Arguments
  integer, intent(in)               :: nSats
  real, intent(out)                 :: Buffer_III(4,2,nSats)
  character (len=100), intent(out)  :: Buffer_I(nSats)
  
  !Internal variables
  character (len=*), parameter :: NameSub='GM_get_sat_for_rb'

  real :: sat_RayVars(5), sat_RayVarsSum(5),SatRay_D(3)

  real :: StateSat_V(0:nVar+3), B0Sat_D(3)  
  real :: Bx,By,Bz,B2
  integer :: iSat, iError
  !--------------------------------------------------------------------------
  ! Store satellite names in Buffer_I
  if (iProc == 0) then
     Buffer_I = Satellite_name(1:nSats)
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
     call GM_trace_sat(xSatellite(iSat,1:3),SatRay_D)
     ! Determine magnetic field magnitude at satellite B=B0+B1
     if(UseB0)then
        call get_b0(xSatellite(iSat,1),xSatellite(iSat,2),xSatellite(iSat,3),&
             B0Sat_D)
     else
        B0Sat_D=0.00
     end if
     call get_point_data(0.0,xSatellite(iSat,:),1,nBlock,1,nVar+3,StateSat_V)
     call collect_satellite_data(xSatellite(iSat,:),StateSat_V)

     Bx = StateSat_V(Bx_)+B0Sat_D(1)
     By = StateSat_V(By_)+B0Sat_D(2)
     Bz = StateSat_V(Bz_)+B0Sat_D(3)
     
     B2 = (Bx**2.0 + By**2.0 + Bz**2.0) * (No2Si_V(UnitB_))**2.0 
     
     ! Store results in Buffer_III
     if (iProc == 0) then 
        Buffer_III(1:3,1,iSat)   = Xsatellite(iSat,1:3)
        !Buffer_III(1:3,2,iSat) = sat_RayVarsSum(1:3)
        Buffer_III(1:3,2,iSat) = SatRay_D
        Buffer_III(4,2,iSat)   = B2
     end if
  end do

end subroutine GM_get_sat_for_rb

