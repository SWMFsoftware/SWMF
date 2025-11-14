module IPE_grid_comp

  ! ESMF Framework module
  use ESMF
  use NUOPC

  use NUOPC_Model, only: NUOPC_ModelGet
  use NUOPC_Model, only: modelSS => SetServices
  use NUOPC_Model, only: model_label_DataInitialize => label_DataInitialize
  use NUOPC_Model, only: model_label_Advance => label_Advance
  use NUOPC_Model, only: model_label_Finalize => label_Finalize

  use ESMFSWMF_variables, ONLY: &
       NameFieldImport_V, nVarImport, NameFieldExport_V, nVarExport, &
       add_fields, DoTest, FieldTest_V, CoordCoefTest, &
       dHallPerDtTest, write_log, write_error

  implicit none
  private

  public:: set_services

  ! ESMF "dynamo" grid and domain. The real grid is not uniform in
  ! latitude, which will be implemented in the near future !!!
  ! This is a 2D spherical grid in MAG coordinates (rotates with Earth):
  ! +Z points to north magnetic dipole, +Y is towards rotational Omega x Z

  integer, parameter:: nLon=81, nLat=97 ! Default ESMF grid size

  ! Coordinate arrays
  integer:: MinLon, MaxLon, MinLat, MaxLat
  real(ESMF_KIND_R8), pointer:: Lon_I(:), Lat_I(:)

  ! IPE dynamo grid latitudes
  real, parameter:: LatIpe_I(nLat) = [ &
       -90.0000, -88.1238, -86.2386, -84.3344, -82.4013, -80.4296, -78.4095, &
       -76.3318, -74.1877, -71.9690, -69.6682, -67.2793, -64.7977, -62.2208, &
       -59.5484, -56.7835, -53.9323, -51.0045, -48.0138, -44.9776, -41.9167, &
       -38.8546, -35.8165, -32.8285, -29.9165, -27.1046, -24.4146, -21.8655, &
       -19.4724, -17.2473, -15.1984, -13.3307, -11.6462, -10.1443,  -8.8219, &
       -7.6733,   -6.6900,  -5.8603,  -5.1688,  -4.5959,  -4.1191,  -3.7133, &
       -3.3532,   -3.0142,  -2.6728,  -2.3049,  -1.8786,  -1.3276,   0.0000, &
       1.3276,   1.8786,  2.3049,  2.6728,  3.0142,  3.3532,  3.7133,  &
       4.1191,   4.5959,  5.1688,  5.8603,  6.6900,  7.6733,  8.8219,  &
       10.1443, 11.6462, 13.3307, 15.1984, 17.2473, 19.4724, 21.8655,  &
       24.4146, 27.1046, 29.9165, 32.8285, 35.8165, 38.8546, 41.9167,  &
       44.9776, 48.0138, 51.0045, 53.9323, 56.7835, 59.5484, 62.2208,  &
       64.7977, 67.2793, 69.6682, 71.9690, 74.1877, 76.3318, 78.4095,  &
       80.4296, 82.4013, 84.3344, 86.2386, 88.1238, 90.0000 ]

  ! Coupling interval
  integer :: iCoupleFreq

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
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_DataInitialize, &
         specRoutine=my_data_init, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_CompSpecialize')
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=my_run, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_CompSpecialize')
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=my_final, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_CompSpecialize')

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
    call write_log("IPE_grid_comp:init_advertise routine called")
    iError = ESMF_FAILURE

    ! Note that NameFieldImport_V and NameFieldExport_V are defined for RIM
    ! RIM export is IPE import and RIM import is IPE export
    do n = 1, nVarImport
       ! IPE -> RIM coupling
       call NUOPC_Advertise(ExportState, standardName=trim(NameFieldImport_V(n)), &
            TransferOfferGeomObject='will provide', rc=iError)
       if(iError /= ESMF_SUCCESS) call my_error('NUOPC_Advertise')
    end do

    do n = 1, nVarExport
       ! RIM -> IPE coupling
       call NUOPC_Advertise(ImportState, standardName=trim(NameFieldExport_V(n)), &
            TransferOfferGeomObject='will provide', rc=iError)
       if(iError /= ESMF_SUCCESS) call my_error('NUOPC_Advertise')
    end do

    iError = ESMF_SUCCESS
    call write_log("IPE_grid_comp:init_advertise routine returned")

  end subroutine my_init_advertise
  !============================================================================
  subroutine my_init_realize(gComp, ImportState, ExportState, Clock, iError)

    type(ESMF_GridComp):: gComp
    type(ESMF_State)   :: ImportState
    type(ESMF_State)   :: ExportState
    type(ESMF_Clock)   :: Clock
    integer, intent(out):: iError

    type(ESMF_Grid):: Grid
    type(ESMF_Field):: Field
    type(ESMF_VM):: Vm
    real(ESMF_KIND_R8), pointer :: Ptr_II(:,:)
    integer:: iVar, i, j, PetCount
    character(len=4):: NameField
    !--------------------------------------------------------------------------
    call write_log("IPE_grid_comp init_realize called")
    iError = ESMF_FAILURE

    ! Query component
    call ESMF_GridCompGet(gComp, vm=Vm, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCompGet')

    call ESMF_VMGet(Vm, petCount=PetCount, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_VMGet')

    ! Create Lon-Lat grid where -180<=Lon<=180-dLon, -90<=Lat<=90
    Grid = ESMF_GridCreateNoPeriDim(maxIndex=[nLon-1, nLat-1], &
         regDecomp=[1, petCount], coordDep1=[1], coordDep2=[2], &
         coordSys=ESMF_COORDSYS_CART, indexflag=ESMF_INDEX_GLOBAL, &
         name="dynamo grid", rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCreateNoPeriDim')

    call ESMF_GridAddCoord(Grid, staggerloc=ESMF_STAGGERLOC_CORNER, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridAddCoord')

    nullify(Lon_I)
    call ESMF_GridGetCoord(Grid, CoordDim=1, &
         staggerLoc=ESMF_STAGGERLOC_CORNER, farrayPtr=Lon_I, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridGetCoord 1')
    write(*,*)'IPE size(Lon_I)=', size(Lon_I)

    nullify(Lat_I)
    call ESMF_GridGetCoord(Grid, CoordDim=2, &
         staggerLoc=ESMF_STAGGERLOC_CORNER, farrayPtr=Lat_I, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridGetCoord 2')
    write(*,*)'IPE size(Lat_I)=', size(Lat_I)

    ! Uniform longitude grid from -180 to 180
    MinLon = lbound(Lon_I, dim=1)
    MaxLon = ubound(Lon_I, dim=1)
    do i = MinLon, MaxLon
       Lon_I(i) = (i-1)*(360.0/(nLon-1)) - 180
    end do
    write(*,*)'IPE grid: Lon_I(MinLon,MinLon+1,MaxLon)=', Lon_I([MinLon, MinLon+1, MaxLon])
    ! Nonuniform latitude grid
    MinLat = lbound(Lat_I, dim=1)
    MaxLat = ubound(Lat_I, dim=1)
    do j = MinLat, MaxLat
       Lat_I(j) = LatIpe_I(j)
    end do
    write(*,*)'IPE grid: Lat_I(MinLat,MinLat+1,MaxLat)=', Lat_I([MinLat, MinLat+1, MaxLat])

    ! Note that NameFieldImport_V and NameFieldExport_V are defined for RIM
    ! RIM export is IPE import and RIM import is IPE export

    ! Add fields to the export state
    call add_fields(Grid, ExportState, nVarImport, NameFieldImport_V, iError=iError)
    if(iError /= ESMF_SUCCESS) call my_error("add_fields")

    ! Add fields to the import state
    call add_fields(Grid, ImportState, nVarExport, NameFieldExport_V, iError=iError)
    if(iError /= ESMF_SUCCESS) call my_error("add_fields")

    iError = ESMF_SUCCESS
    call write_log("IPE_grid_comp init_realize returned")
    call ESMF_LogFlush()

  end subroutine my_init_realize
  !============================================================================
  subroutine my_data_init(gComp, iError)

    type(ESMF_GridComp):: gComp
    integer, intent(out):: iError

    type(ESMF_Clock) :: Clock
    type(ESMF_TimeInterval) :: TimeStep
    type(ESMF_State):: ExportState
    type(ESMF_Field):: Field
    real(ESMF_KIND_R8), pointer :: Ptr_II(:,:)
    integer:: i, j, iVar
    character(len=4):: NameField
    !--------------------------------------------------------------------------
    call write_log("IPE_grid_comp data_init called")
    iError = ESMF_FAILURE

    call NUOPC_ModelGet(gComp, modelClock=Clock, exportState=ExportState, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error("NUOPC_ModelGet")

    call ESMF_ClockGet(Clock, timeStep=TimeStep, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_ClockGet')

    call ESMF_TimeIntervalGet(TimeStep, s=iCoupleFreq, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_TimeIntervalGet')

    do iVar = 1, nVarImport
       ! Get pointers to the variables in the export state
       nullify(Ptr_II)
       NameField = NameFieldImport_V(iVar)
       call ESMF_StateGet(ExportState, itemName=NameField, field=Field, &
            rc=iError)
       if(iError /= ESMF_SUCCESS) call my_error("ESMF_StateGet "//NameField)

       call ESMF_FieldGet(Field, farrayPtr=Ptr_II, rc=iError)
       if(iError /= ESMF_SUCCESS) call my_error("ESMF_FieldGet "//NameField)

       ! Get dimension extents
       MinLon = lbound(Ptr_II, dim=1)
       MaxLon = ubound(Ptr_II, dim=1)
       MinLat = lbound(Ptr_II, dim=2)
       MaxLat = ubound(Ptr_II, dim=2)

       if(iError /= ESMF_SUCCESS) RETURN
       select case(NameField)
       case('Hall')
          Ptr_II(MinLon:MaxLon,MinLat:MaxLat) = FieldTest_V(1)
       case('Ped')
          Ptr_II(MinLon:MaxLon,MinLat:MaxLat) = FieldTest_V(2)
       case default
          write(*,*)'ERROR in ESMF_GridComp:init: unknown NameField=',&
               NameField,' for iVar=',iVar
          iError = ESMF_FAILURE; RETURN
       end select

       ! Add coordinate dependence
       do j = MinLat, MaxLat; do i = MinLon, MaxLon
          Ptr_II(i,j) = Ptr_II(i,j) + CoordCoefTest &
               *abs(Lon_I(i))*(90-abs(Lat_I(j)))
       end do; end do

    end do ! iVar

    iError = ESMF_SUCCESS
    call write_log("IPE_grid_comp data_init returned")
    call ESMF_LogFlush()

  end subroutine my_data_init
  !============================================================================
  subroutine my_run(gComp, iError)

    type(ESMF_GridComp):: gComp
    integer, intent(out):: iError

    type(ESMF_Field):: Field
    type(ESMF_State):: ExportState

    type(ESMF_State):: ImportState
    character(len=4):: NameField
    real(ESMF_KIND_R8), pointer     :: Ptr_II(:,:)
    real(ESMF_KIND_R8), allocatable :: Data_VII(:,:,:)
    integer                         :: iVar

    integer :: i, j
    real(ESMF_KIND_R8) :: Coef, ExactValue
    !--------------------------------------------------------------------------
    call write_log("IPE_grid_comp run called")
    iError = ESMF_FAILURE

    ! Query component
    call NUOPC_ModelGet(gcomp, exportState=ExportState, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_ModelGet')

    ! We should execute the ESMF code here and put the result into
    ! the fields of the ExportState
    if(DoTest)then
       ! Get pointers to the MHD variables in the export state
       !!! This could be done in the initialization ?!
       nullify(Ptr_II)
       call ESMF_StateGet(ExportState, itemName='Hall', field=Field, rc=iError)
       if(iError /= ESMF_SUCCESS) call my_error("ESMF_StateGet for Hall")

       call ESMF_FieldGet(Field, farrayPtr=Ptr_II, rc=iError)
       if(iError /= ESMF_SUCCESS) call my_error("ESMF_FieldGet for Hall")

       ! Update state by changing Hall conductivity
       ! Get dimension extents
       MinLon = lbound(Ptr_II, dim=1)
       MaxLon = ubound(Ptr_II, dim=1)
       MinLat = lbound(Ptr_II, dim=2)
       MaxLat = ubound(Ptr_II, dim=2)
       write(*,*)'IPE_grid_comp:run old Hall=', &
            Ptr_II((MaxLon+MinLon)/2,(MaxLat+MinLat)/2)
       Ptr_II = Ptr_II + iCoupleFreq*dHallPerDtTest
       write(*,*)'IPE_grid_comp:run new Hall=', &
            Ptr_II((MaxLon+MinLon)/2,(MaxLat+MinLat)/2)
    end if

    ! Query component
    call NUOPC_ModelGet(gComp, importState=ImportState, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_ModelGet')

    allocate(Data_VII(nVarExport,MinLon:MaxLon,MinLat:MaxLat), stat=iError)
    if(iError /= 0) call my_error('allocate(Data_VII)')

    do iVar = 1, nVarExport
       nullify(Ptr_II)
       NameField = NameFieldExport_V(iVar)
       call ESMF_StateGet(ImportState, itemName=NameField, &
            field=Field, rc=iError)
       if(iError /= ESMF_SUCCESS) call my_error("ESMF_StateGet")
       call ESMF_FieldGet(Field, farrayPtr=Ptr_II, rc=iError)
       if(iError /= ESMF_SUCCESS) call my_error("ESMF_FieldGet")

       Data_VII(iVar,:,:) = Ptr_II

       Coef = 10**iVar
       ! Skip boundary longitudes (+180 and -180), where the interpolated 
       ! values and the ExactValue are different.  
       do i = MinLon+1, MaxLon-1; do j = MinLat, MaxLat
          ExactValue = Lon_I(i)*Lat_I(j)*Coef
          if( abs(Data_VII(iVar,i,j)-ExactValue) > 1.0e-6 ) then
             write(*,*)'Error in RIM->IPE coupling for ', NameField, &
                  ' at lon=', Lon_I(i), ' lat=', Lat_I(j), &
                  ' value=', Data_VII(iVar,i,j), ' ExactValue=', ExactValue
          end if
       end do; end do
    end do

    ! Clean memory
    deallocate(Data_VII, stat=iError)
    if(iError /= 0) call my_error('deallocate(Data_VII)')

    iError = ESMF_SUCCESS
    call write_log("IPE_grid_comp run returned")

  end subroutine my_run
  !============================================================================
  subroutine my_final(gComp, iError)

    type(ESMF_GridComp) :: gComp
    integer, intent(out):: iError

    !--------------------------------------------------------------------------
    call write_log("IPE_grid_comp finalize called")
    call write_log("IPE_grid_comp finalize returned")

    iError = ESMF_SUCCESS

  end subroutine my_final
  !============================================================================
  subroutine my_error(String)

    ! Write out error message and stop

    character(len=*), intent(in) :: String
    !--------------------------------------------------------------------------
    call write_error('IPE_grid_comp '//String)

  end subroutine my_error
  !============================================================================
end module IPE_grid_comp
!==============================================================================
