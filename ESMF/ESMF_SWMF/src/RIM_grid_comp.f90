module RIM_grid_comp

  ! This is the RIM Gridded Component, which acts as an interface
  ! to the SWMF/IE/RIM model.

  ! ESMF Framework module
  use ESMF

  use ESMFSWMF_variables, ONLY: &
       nProcSwmfComp, NameFieldEsmf_V, nVarEsmf, &
       add_fields, write_log, write_error, &
       DoTest, FieldTest_V, CoordCoefTest, dHallPerDtTest

  ! Size of ionosphere grid in SWMF/IE model
  use IE_ModSize, ONLY: nLat => IONO_nTheta, nLon => IONO_nPsi

  ! Conversion to radians

  implicit none

  private

  public:: set_services

  ! Local variables

  ! Grid related to the SWMF Ridley Ionosphere Model
  ! This is a 2D spherical grid representing the height integrated ionosphere.
  ! RIM is in SM coordinates (aligned with Sun-Earth direction):
  ! +Z points to north magnetic dipole and the Sun is in the +X-Z halfplane.

  ! Coordinate arrays
  real(ESMF_KIND_R8), pointer, save:: Lon_I(:), Lat_I(:)

  integer:: MinLat = 0, MaxLat = 0
  integer, allocatable:: iPetMap_III(:,:,:)
  integer:: iPet = 0

contains
  !============================================================================
  subroutine set_services(gComp, iError)
    type(ESMF_GridComp) :: gComp
    integer, intent(out):: iError

    !--------------------------------------------------------------------------
    call ESMF_GridCompSetEntryPoint(gComp, ESMF_METHOD_INITIALIZE, &
         userRoutine=my_init, rc=iError)
    call ESMF_GridCompSetEntryPoint(gComp, ESMF_METHOD_RUN, &
         userRoutine=my_run, rc=iError)
    call ESMF_GridCompSetEntryPoint(gComp, ESMF_METHOD_FINALIZE, &
         userRoutine=my_final, rc=iError)

  end subroutine set_services
  !============================================================================
  subroutine my_init(gComp, ImportState, ExportState, ExternalClock, iError)

    type(ESMF_GridComp) :: gComp
    type(ESMF_State) :: ImportState
    type(ESMF_State) :: ExportState
    type(ESMF_Clock) :: ExternalClock
    integer, intent(out):: iError

    type(ESMF_VM)    :: Vm
    type(ESMF_Grid)  :: Grid
    integer          :: iComm
    type(ESMF_TimeInterval) :: SimTime, RunDuration

    integer:: i
    !--------------------------------------------------------------------------
    call write_log("RIM_grid_comp:init routine called")
    iError = ESMF_FAILURE

    ! Obtain the VM for the IE gridded component
    call ESMF_GridCompGet(gComp, vm=Vm, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_GridCompGet')

    ! Obtain the MPI communicator for the VM
    call ESMF_VMGet(Vm, mpiCommunicator=iComm, localPet=iPet, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_VMGet')

    ! RIM grid is node based. Internally it is Colat-Lon grid, but we pretend
    ! here that it is a Lat-Lon grid, so ESMF can use it.
    ! Lon from 0 to 360-dPhi (periodic), Lat from -90 to +90
    if(nProcSwmfComp == 1)then
       allocate(iPetMap_III(1,1,1))
       iPetMap_III(1,1,1) = 0
    else
       allocate(iPetMap_III(1,2,1))
       iPetMap_III(1,:,1) = [ 0, 1 ]
    end if

    Grid = ESMF_GridCreateNoPeriDim(maxIndex=[nLon-1, nLat-1], &
         RegDecomp = [1, nProcSwmfComp], &
         coordDep1=[1], coordDep2=[2], coordSys=ESMF_COORDSYS_CART, &
         name="RIM grid", rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCreateNoPeriDim')
    deallocate(iPetMap_III)

    call ESMF_GridAddCoord(Grid, staggerloc=ESMF_STAGGERLOC_CORNER, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridAddCoord')

    nullify(Lon_I)
    call ESMF_GridGetCoord(Grid, CoordDim=1, &
         staggerLoc=ESMF_STAGGERLOC_CORNER, farrayPtr=Lon_I, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridGetCoord 1')
    write(*,*)'ESMF_GridComp size(Lon_I)=', size(Lon_I)

    call ESMF_GridGetCoord(Grid, CoordDim=2, &
         staggerLoc=ESMF_STAGGERLOC_CORNER, farrayPtr=Lat_I, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridGetCoord 2')

    MinLat = lbound(Lat_I, DIM=1); MaxLat = ubound(Lat_I, DIM=1)
    write(*,'(a,2i4)')'RIM grid: MinLat, MaxLat=', MinLat, MaxLat

    ! Uniform longitude grid from -180 to 180
    do i = 1, nLon
       Lon_I(i) = (i-1)*(360.0/(nLon-1)) - 180
    end do
    write(*,*)'RIM grid: Lon_I(1,2,last)=', Lon_I([1,1,nLon])
    ! Uniform latitude grid
    do i = MinLat, MaxLat
       Lat_I(i) = (i-1)*(180./(nLat-1)) - 90
    end do
    write(*,*)'RIM grid: Lat_I(Min,Min+1,Max)=', &
         Lat_I([MinLat,MinLat+1,MaxLat])

    ! Add fields to the RIM import state
    call add_fields(Grid, ImportState, IsFromEsmf=.true., iError=iError)
    if(iError /= ESMF_SUCCESS) call my_error('add_fields')

    iError = ESMF_SUCCESS
    call write_log("RIM_grid_comp:init routine returned")

  end subroutine my_init
  !============================================================================
  subroutine my_run(gComp, ImportState, ExportState, Clock, iError)

    type(ESMF_GridComp):: gComp
    type(ESMF_State):: ImportState
    type(ESMF_State):: ExportState
    type(ESMF_Clock):: Clock
    integer, intent(out):: iError

    ! Access to the data
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

    ! Get the current time from the clock
    call ESMF_ClockGet(Clock, CurrSimTime=SimTime, TimeStep=TimeStep, &
         rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_ClockGet')
    call ESMF_TimeIntervalGet(SimTime, s=iSec, ms=iMilliSec, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_TimeIntervalGet current')
    tCurrent = iSec + 0.001*iMilliSec

    ! Obtain pointer to the data obtained from the ESMF component
    allocate(Data_VII(nVarEsmf,1:nLon,MinLat:MaxLat), stat=iError)
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
       do j = MinLat, MaxLat; do i = 1, nLon
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
            Data_VII(:,nLon/2,(MinLat+MaxLat)/2)
    end if
    deallocate(Data_VII)

    call write_log("RIM_grid_comp:run routine returned")

    iError = ESMF_SUCCESS

  end subroutine my_run
  !============================================================================
  subroutine my_final(gComp, ImportState, ExportState, Externalclock, iError)

    type(ESMF_GridComp) :: gComp
    type(ESMF_State) :: ImportState
    type(ESMF_State) :: ExportState
    type(ESMF_Clock) :: Externalclock
    integer, intent(out) :: iError
    !--------------------------------------------------------------------------
    call write_log("RIM_finalize routine called")

    call write_log("RIM_finalize routine returned")

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
