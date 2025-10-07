module RIM_grid_comp

  ! This is the RIM Gridded Component, which acts as an interface
  ! to the SWMF/IE/RIM model.

  ! ESMF Framework module
  use ESMF
  use NUOPC

  use NUOPC_Model, only: NUOPC_ModelGet
  use NUOPC_Model, only: modelSS => SetServices
  use NUOPC_Model, only: model_label_Advance => label_Advance
  use NUOPC_Model, only: model_label_Finalize => label_Finalize

  use ESMFSWMF_variables, ONLY: &
       NameFieldEsmf_V, nVarEsmf, &
       add_fields, write_log, write_error, &
       DoTest, FieldTest_V, CoordCoefTest, dHallPerDtTest

  ! Size of ionosphere grid in SWMF/IE model
  use IE_ModSize, ONLY: nLat => IONO_nTheta, nLon => IONO_nPsi

  use CON_axes, ONLY: transform_matrix
  use ModNumConst, ONLY: cRadToDeg, cDegToRad

  ! Conversion to radians

  implicit none

  private

  public:: set_services

  ! Local variables

  ! Grid related to the SWMF Ridley Ionosphere Model
  ! This is a 2D spherical grid representing the height integrated ionosphere.
  ! RIM is in SM coordinates (aligned with Sun-Earth direction):
  ! +Z points to north magnetic dipole and the Sun is in the +X-Z halfplane.

  ! Coordinate arrays. They store the RIM mesh node coordinates 
  ! in MAG (used by IPE) coordinates. They are used for coupling.
  real(ESMF_KIND_R8), pointer, save:: Lon_I(:), Lat_I(:)  

  integer:: MinLon, MaxLon, MinLat, MaxLat

contains
  !============================================================================
  subroutine set_services(gComp, iError)
    type(ESMF_GridComp) :: gComp
    integer, intent(out):: iError

    !--------------------------------------------------------------------------
    call NUOPC_CompDerive(gComp, modelSS, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_CompDerive')
    call ESMF_GridCompSetEntryPoint(gComp, ESMF_METHOD_INITIALIZE, &
         userRoutine=my_init_p0, phase=0, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_GridCompSetEntryPoint')
    call NUOPC_CompSetEntryPoint(gComp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=["IPDv01p1"], userRoutine=my_init_advertise, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_CompSetEntryPoint')
    call NUOPC_CompSetEntryPoint(gComp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=["IPDv01p3"], userRoutine=my_init_realize, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_CompSetEntryPoint')
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
          specRoutine=my_run, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_CompSetEntryPoint')
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=my_final, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_CompSetEntryPoint')

  end subroutine set_services
  !============================================================================
  subroutine my_init_p0(gComp, ImportState, ExportState, ExternalClock, iError)

    type(ESMF_GridComp) :: gComp
    type(ESMF_State) :: ImportState
    type(ESMF_State) :: ExportState
    type(ESMF_Clock) :: ExternalClock
    integer, intent(out):: iError

    !--------------------------------------------------------------------------
    call NUOPC_CompFilterPhaseMap(gComp, ESMF_METHOD_INITIALIZE, &
         acceptStringList=["IPDv01p"], rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_CompFilterPhaseMap')

  end subroutine my_init_p0
  !============================================================================
  subroutine my_init_advertise(gComp, ImportState, ExportState, ExternalClock, iError)

    type(ESMF_GridComp) :: gComp
    type(ESMF_State) :: ImportState
    type(ESMF_State) :: ExportState
    type(ESMF_Clock) :: ExternalClock
    integer, intent(out):: iError

    integer :: n
    !--------------------------------------------------------------------------
    call write_log("RIM_grid_comp:init_advertise routine called")
    iError = ESMF_FAILURE

    do n = 1, nVarEsmf
      call NUOPC_Advertise(ImportState, standardName=trim(NameFieldEsmf_V(n)), &
           TransferOfferGeomObject='will provide', rc=iError)
      if(iError /= ESMF_SUCCESS) call my_error('NUOPC_Advertise')
    end do

    iError = ESMF_SUCCESS
    call write_log("RIM_grid_comp:init_advertise routine returned")

  end subroutine my_init_advertise
  !============================================================================
  subroutine my_init_realize(gComp, ImportState, ExportState, ExternalClock, iError)

    type(ESMF_GridComp) :: gComp
    type(ESMF_State) :: ImportState
    type(ESMF_State) :: ExportState
    type(ESMF_Clock) :: ExternalClock
    integer, intent(out):: iError

    type(ESMF_VM)    :: Vm
    type(ESMF_Grid)  :: Grid
    integer          :: PetCount
    type(ESMF_TimeInterval) :: SimTime, RunDuration
    
    integer(ESMF_KIND_I4)   :: iSec, iMilliSec
    real(ESMF_KIND_R8) :: tCurrent

    real(ESMF_KIND_R8), allocatable :: LonSM_I(:)

    real(ESMF_KIND_R8) :: SmToMag_DD(3,3)

    ! Phi is the rotation angle from SM to MAG coordinates 
    ! in the counter-clockwise direction.
    real(ESMF_KIND_R8) :: CosPhi, SinPhi, Phi

    integer:: i, j
    !--------------------------------------------------------------------------
    call write_log("RIM_grid_comp:init_realize routine called")
    iError = ESMF_FAILURE

    ! Obtain the VM for the IE gridded component
    call ESMF_GridCompGet(gComp, vm=Vm, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_GridCompGet')

    ! Obtain the PET count for the VM
    call ESMF_VMGet(Vm, petCount=PetCount, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_VMGet')

    ! RIM grid is node based. Internally it is Colat-Lon grid, but we pretend
    ! here that it is a Lat-Lon grid, so ESMF can use it.
    ! Lon from 0 to 360-dPhi (periodic), Lat from -90 to +90
    Grid = ESMF_GridCreateNoPeriDim(maxIndex=[nLon-1, nLat-1], &
         regDecomp=[1, petCount], coordDep1=[1], coordDep2=[2], &
         coordSys=ESMF_COORDSYS_CART, indexflag=ESMF_INDEX_GLOBAL, &
         name="RIM grid", rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCreateNoPeriDim')

    call ESMF_GridAddCoord(Grid, staggerloc=ESMF_STAGGERLOC_CORNER, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridAddCoord')

    nullify(Lon_I)
    call ESMF_GridGetCoord(Grid, CoordDim=1, &
         staggerLoc=ESMF_STAGGERLOC_CORNER, farrayPtr=Lon_I, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridGetCoord 1')
    write(*,*)'ESMF_GridComp size(Lon_I)=', size(Lon_I)

    nullify(Lat_I)
    call ESMF_GridGetCoord(Grid, CoordDim=2, &
         staggerLoc=ESMF_STAGGERLOC_CORNER, farrayPtr=Lat_I, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridGetCoord 2')
    write(*,*)'ESMF_GridComp size(Lat_I)=', size(Lat_I)

    ! Uniform longitude grid from -180 to 180 (to match IPE)
    MinLon = lbound(Lon_I, dim=1)
    MaxLon = ubound(Lon_I, dim=1)
    write(*,'(a,2i4)')'RIM grid: Lon_I Min, Max=', MinLon, MaxLon

    allocate(LonSM_I(MinLon:MaxLon))
    
    do i = MinLon, MaxLon
       LonSM_I(i) = (i-1)*(360.0/(nLon-1)) - 180
    end do

    write(*,*)'RIM grid: LonSM_I(Min,Min+1,Max)=', LonSM_I([MinLon,MinLon+1,MaxLon])

    ! Get the current time from the clock
    call ESMF_ClockGet(ExternalClock, CurrSimTime=SimTime, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_ClockGet')
    call ESMF_TimeIntervalGet(SimTime, s=iSec, ms=iMilliSec, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_TimeIntervalGet current')
    tCurrent = iSec + 0.001*iMilliSec
    write(*,*)'RIM_grid_comp init_realize: current time, isec, msec=', &
         tCurrent, iSec, iMilliSec

    SmToMag_DD = transform_matrix(tCurrent, 'SMG', 'MAG') 
    write(*,*)'RIM_grid_comp: SM to MAG matrix='
    do i = 1, 3
       write(*,'(3f12.6)') SmToMag_DD(i,:)
    end do
    CosPhi = SmToMag_DD(1, 1)
    SinPhi = SmToMag_DD(1, 2)
    Phi = atan2(SinPhi, CosPhi)*cRadToDeg
    write(*,*)'RIM_grid_comp: rotation angle from SM to MAG=', Phi

    ! Make sure the range is [-180, 180] after the shift by Theta.
    Lon_I = modulo(LonSM_I + 180 - Phi, 360.0) - 180
    write(*,*)'RIM grid: Lon_I(Min,Min+1,Max)=', Lon_I([MinLon,MinLon+1,MaxLon])

    !do i = MinLon, MaxLon
     !write(*,*)'lonsm lonmag lon sm=', i, LonSM_I(i), Lon_I(i)
    !end do


    ! Uniform latitude grid
    MinLat = lbound(Lat_I, dim=1)
    MaxLat = ubound(Lat_I, dim=1)
    write(*,'(a,2i4)')'RIM grid: Lat_I Min, Max=', MinLat, MaxLat
    do j = MinLat, MaxLat
       Lat_I(j) = (j-1)*(180./(nLat-1)) - 90
    end do
    write(*,*)'RIM grid: Lat_I(Min,Min+1,Max)=', Lat_I([MinLat,MinLat+1,MaxLat])

    ! Add fields to the RIM import state
    call add_fields(Grid, ImportState, IsFromEsmf=.true., iError=iError)
    if(iError /= ESMF_SUCCESS) call my_error('add_fields')

    if(allocated(LonSM_I)) deallocate(LonSM_I)

    iError = ESMF_SUCCESS
    call write_log("RIM_grid_comp:init_realize routine returned")

  end subroutine my_init_realize
  !============================================================================
  subroutine my_run(gComp, iError)

    type(ESMF_GridComp):: gComp
    integer, intent(out):: iError

    ! Access to the data
    type(ESMF_State):: ImportState
    type(ESMF_Clock):: Clock
    real(ESMF_KIND_R8), pointer     :: Ptr_II(:,:)
    real(ESMF_KIND_R8), allocatable :: Data_VII(:,:,:)
    integer                         :: iVar

    ! Access to time
    type(ESMF_TimeInterval) :: SimTime, TimeStep
    integer(ESMF_KIND_I4)   :: iSec, iMilliSec

    ! Current time (needed for test only)
    real(ESMF_KIND_R8) :: tCurrent

    ! Misc variables
    type(ESMF_Field):: Field
    character(len=4):: NameField
    integer:: i, j
    real(ESMF_KIND_R8):: Exact_V(2)
    !--------------------------------------------------------------------------
    call write_log("RIM_grid_comp:run routine called")
    iError = ESMF_FAILURE

    ! Query component
    call NUOPC_ModelGet(gComp, modelClock=Clock, &
         importState=ImportState, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_ModelGet')

    ! Get the current time from the clock
    call ESMF_ClockGet(Clock, CurrSimTime=SimTime, TimeStep=TimeStep, &
         rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_ClockGet')
    call ESMF_TimeIntervalGet(SimTime, s=iSec, ms=iMilliSec, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_TimeIntervalGet current')
    tCurrent = iSec + 0.001*iMilliSec

    ! Obtain pointer to the data obtained from the ESMF component
    allocate(Data_VII(nVarEsmf,MinLon:MaxLon,MinLat:MaxLat), stat=iError)
    if(iError /= 0) call my_error('allocate(Data_VII)')

    ! Copy fields into an array
    do iVar = 1, nVarEsmf
       nullify(Ptr_II)
       NameField = NameFieldEsmf_V(iVar)
       call ESMF_StateGet(ImportState, itemName=NameField, &
            field=Field, rc=iError)
       if(iError /= ESMF_SUCCESS) call my_error("ESMF_StateGet")
       call ESMF_FieldGet(Field, farrayPtr=Ptr_II, rc=iError)
       if(iError /= ESMF_SUCCESS) call my_error("ESMF_FieldGet")

       Data_VII(iVar,:,:) = Ptr_II
    end do
    if(DoTest)then
       write(*,*)'SWMF_GridComp shape of Ptr =', shape(Ptr_II)
       ! Do not check the poles
       do j = MinLat, MaxLat; do i = MinLon, MaxLon
          ! Calculate exact solution
          Exact_V = FieldTest_V
          ! add time dependence for Hall field
          Exact_V(1) = Exact_V(1) + tCurrent*dHallPerDtTest
          ! add spatial dependence
          Exact_V = Exact_V + CoordCoefTest &
               *abs(Lon_I(i))*(90-abs(Lat_I(j)))
          ! *sin(Lon_I(i)*cDegToRad)*cos(Lat_I(j)*cDegToRad)
          if(abs(Data_VII(1,i,j) - Exact_V(1)) > 1e-10) &
               write(*,*) 'ERROR in SWMF_GridComp ', &
               'at i, j, Lon, Lat, Hall, Exact, Error=', &
               i, j, Lon_I(i), Lat_I(j), Data_VII(1,i,j), &
               Exact_V(1), Data_VII(1,i,j) - Exact_V(1)
          if(abs(Data_VII(2,i,j) - Exact_V(2)) > 1e-10) &
               write(*,*) 'ERROR in SWMF_GridComp ', &
               'at i, j, Lon, Lat, Pede, Exact, Error=', &
               i, j, Lon_I(i), Lat_I(j), Data_VII(2,i,j), &
               Exact_V(2), Data_VII(2,i,j) - Exact_V(2)
       end do; end do
       write(*,*)'SWMF_GridComp value of Data(MidLon,MidLat)=', &
            Data_VII(:,(MinLon+MaxLon)/2,(MinLat+MaxLat)/2)
    end if
    deallocate(Data_VII)

    call write_log("RIM_grid_comp:run routine returned")

    iError = ESMF_SUCCESS

  end subroutine my_run
  !============================================================================
  subroutine my_final(gComp, iError)

    type(ESMF_GridComp) :: gComp
    integer, intent(out) :: iError
    !--------------------------------------------------------------------------
    call write_log("RIM_finalize routine called")

    call write_log("RIM_finalize routine returned")

    iError = ESMF_SUCCESS

  end subroutine my_final
  !============================================================================
  subroutine my_error(String)

    ! Write out error message and stop

    character(len=*), intent(in) :: String
    !--------------------------------------------------------------------------

    call write_error('RIM_grid_comp '//String)

  end subroutine my_error
  !============================================================================
end module RIM_grid_comp
!==============================================================================
