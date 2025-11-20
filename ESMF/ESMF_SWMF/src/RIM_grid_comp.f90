module RIM_grid_comp

  ! This is the RIM Gridded Component, which acts as an interface
  ! to the SWMF/IE/RIM model.

  ! ESMF Framework module
  use ESMF
  use NUOPC

  use NUOPC_Model, only: NUOPC_ModelGet
  use NUOPC_Model, only: modelSS => SetServices
  use NUOPC_Model, only: model_label_DataInitialize => label_DataInitialize
  use NUOPC_Model, only: model_label_Advance => label_Advance
  use NUOPC_Model, only: model_label_Finalize => label_Finalize

  use ESMFSWMF_variables, ONLY: &
       NameFieldIpe2Rim_V, nVarIpe2Rim, NameFieldRim2Ipe_V, nVarRim2Ipe, &
       add_fields, write_log, write_error, &
       DoTest, FieldTest_V, CoordCoefTest, dHallPerDtTest, &
       get_coords, update_coordinates, DoShiftDataCoupling, &
       get_sm_to_mag_angle

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
  real(ESMF_KIND_R8), allocatable :: LonSM_I(:), LatSM_I(:)
  integer:: MinLon, MaxLon, MinLat, MaxLat

  ! dPhiSm2Mag is the rotation angle from SM to MAG coordinates 
  ! in the counter-clockwise direction.
  real(ESMF_KIND_R8) :: dPhiSm2Mag

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

    do n = 1, nVarIpe2Rim
       ! IPE -> RIM coupling 
       call NUOPC_Advertise(ImportState, standardName=trim(NameFieldIpe2Rim_V(n)), &
            TransferOfferGeomObject='will provide', rc=iError)
       if(iError /= ESMF_SUCCESS) call my_error('NUOPC_Advertise - import')
    end do

    do n = 1, nVarRim2Ipe
       ! RIM -> IPE coupling
       call NUOPC_Advertise(ExportState, standardName=trim(NameFieldRim2Ipe_V(n)), &
            TransferOfferGeomObject='will provide', rc=iError)
       if(iError /= ESMF_SUCCESS) call my_error('NUOPC_Advertise - export')
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
    type(ESMF_Grid)  :: ImportGrid, ExportGrid
    integer          :: PetCount, i, j

    real(ESMF_KIND_R8), pointer :: Lon_I(:), Lat_I(:)  
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
    ImportGrid = ESMF_GridCreateNoPeriDim(maxIndex=[nLon-1, nLat-1], &
         regDecomp=[1, petCount], coordDep1=[1], coordDep2=[2], &
         coordSys=ESMF_COORDSYS_CART, indexflag=ESMF_INDEX_GLOBAL, &
         name="RIM grid", rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCreateNoPeriDim')

    ! Add corner coordinates to the grid
    call ESMF_GridAddCoord(ImportGrid, staggerloc=ESMF_STAGGERLOC_CORNER, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridAddCoord')

    !----- Get Lon/Lat ranges ------------------------------------------------
    call get_coords(ImportGrid, Lon_I, Lat_I, iError)
    if(iError /= ESMF_SUCCESS) call my_error('get_coords')

    MinLon = lbound(Lon_I, dim=1)
    MaxLon = ubound(Lon_I, dim=1)
    MinLat = lbound(Lat_I, dim=1)
    MaxLat = ubound(Lat_I, dim=1)
    write(*,'(a,2i4)')'RIM grid: Lon_I Min, Max=', MinLon, MaxLon
    write(*,'(a,2i4)')'RIM grid: Lat_I Min, Max=', MinLat, MaxLat
    !--------------------------------------------------------------------------

    allocate(LonSM_I(MinLon:MaxLon))
    allocate(LatSM_I(MinLat:MaxLat))

    do i = MinLon, MaxLon
       LonSM_I(i) = (i-1)*(360.0/(nLon-1)) - 180
    end do
    write(*,*)'RIM grid: LonSM_I(Min,Min+1,Max)=', LonSM_I([MinLon,MinLon+1,MaxLon])

    do j = MinLat, MaxLat
       LatSM_I(j) = (j-1)*(180./(nLat-1)) - 90
    end do
    write(*,*)'RIM grid: LatSM_I(Min,Min+1,Max)=', LatSM_I([MinLat,MinLat+1,MaxLat])

    if(DoShiftDataCoupling) then 
       Lon_I = LonSM_I
       Lat_I = LatSM_I
       call get_sm_to_mag_angle(ExternalClock, dPhiSm2Mag, iError)
    else 
       ! Sets the corner coordinates
       call update_coordinates(ImportGrid, ExternalClock, LonSM_I, LatSM_I, &
            .true., dPhiSm2Mag = dPhiSm2Mag, iError=iError)
       if(iError /= ESMF_SUCCESS) call my_error('update_coordinates')
    end if

    ! Add fields to the RIM import state
    call add_fields(ImportGrid, ImportState, nVarIpe2Rim, NameFieldIpe2Rim_V, iError=iError)
    if(iError /= ESMF_SUCCESS) call my_error('add_fields - import')

    !---- Create Export grid (same as Import grid) -----------------------------
    ExportGrid = ESMF_GridCreateNoPeriDim(maxIndex=[nLon-1, nLat-1], &
         regDecomp=[1, petCount], coordDep1=[1], coordDep2=[2], &
         coordSys=ESMF_COORDSYS_CART, indexflag=ESMF_INDEX_GLOBAL, &
         name="RIM grid", rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCreateNoPeriDim')

    ! Add corner coordinates to the grid
    call ESMF_GridAddCoord(ExportGrid, staggerloc=ESMF_STAGGERLOC_CORNER, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridAddCoord')

    call get_coords(ExportGrid, Lon_I, Lat_I, iError)
    if(iError /= ESMF_SUCCESS) call my_error('get_coords')

    Lon_I = LonSM_I
    Lat_I = LatSM_I
    !--------------------------------------------------------------------------

    ! Add fields to the RIM export state
    call add_fields(ExportGrid, ExportState, nVarRim2Ipe, NameFieldRim2Ipe_V, iError=iError)
    if(iError /= ESMF_SUCCESS) call my_error('add_fields - export')


    iError = ESMF_SUCCESS
    call write_log("RIM_grid_comp:init_realize routine returned")

  end subroutine my_init_realize
  !============================================================================
  subroutine my_data_init(gComp, iError)
    type(ESMF_GridComp):: gComp
    integer, intent(out):: iError
    !--------------------------------------------------------------------------

    call update_export_state(gComp, iError)

  end subroutine my_data_init
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
    type(ESMF_Time) :: CurrTime
    type(ESMF_TimeInterval) :: SimTime, TimeStep
    integer(ESMF_KIND_I4)   :: iSec, iMilliSec

    ! Current time (needed for test only)
    real(ESMF_KIND_R8) :: tCurrent

    ! Grid 
    type(ESMF_Grid) :: Grid

    ! Misc variables
    type(ESMF_Field):: Field
    character(len=4):: NameField
    character(len=ESMF_MAXSTR) :: timeStr
    integer:: i, j, iLeft, iRight
    real(ESMF_KIND_R8):: Exact_V(2), LonMag, dLon, CoefL, CoefR    
    !--------------------------------------------------------------------------
    call write_log("RIM_grid_comp:run routine called")
    iError = ESMF_FAILURE

    ! Query component
    call NUOPC_ModelGet(gComp, modelClock=Clock, &
         importState=ImportState, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('NUOPC_ModelGet')

    ! Get the current time from the clock
    call ESMF_ClockGet(Clock, CurrSimTime=SimTime, TimeStep=TimeStep, &
         currTime=CurrTime, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_ClockGet')
    call ESMF_TimeIntervalGet(SimTime, s=iSec, ms=iMilliSec, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_TimeIntervalGet current')
    tCurrent = iSec + 0.001*iMilliSec

    ! Get time in ISO format
    call ESMF_TimeGet(CurrTime, timeStringISOFrac=timeStr , rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error('ESMF_TimeGet ISO')

    ! Obtain pointer to the data obtained from the ESMF component
    allocate(Data_VII(nVarIpe2Rim,MinLon:MaxLon,MinLat:MaxLat), stat=iError)
    if(iError /= 0) call my_error('allocate(Data_VII)')

    ! Copy fields into an array
    do iVar = 1, nVarIpe2Rim
       nullify(Ptr_II)
       NameField = NameFieldIpe2Rim_V(iVar)
       call ESMF_StateGet(ImportState, itemName=NameField, &
            field=Field, rc=iError)
       if(iError /= ESMF_SUCCESS) call my_error("ESMF_StateGet")
       !   call ESMF_FieldWrite(Field, "rim_import_"//trim(timeStr)//".nc", &
       !        overwrite=.true., rc=iError)
       !   if(iError /= ESMF_SUCCESS) call my_error("ESMF_FieldWrite")
       call ESMF_FieldGet(Field, farrayPtr=Ptr_II, rc=iError)
       if(iError /= ESMF_SUCCESS) call my_error("ESMF_FieldGet")

       if(DoShiftDataCoupling) then        
          ! Ptr_II is the data in IPE(MAG) coordinates, need to shift and interpolate 
          ! to RIM(SM) coordinates. 
          do i = MinLon, MaxLon
             dLon = LonSM_I(2) - LonSM_I(1)
             LonMag = modulo(LonSM_I(i) + 180 - dPhiSm2Mag, 360.0) - 180

             ! The starting position of Ptr_II lon (in MAG) is 
             ! the same as LonSM_I(1)
             iLeft = floor((LonMag - LonSM_I(1))/dLon) + 1           
             iRight = iLeft + 1
             CoefL = (LonSM_I(iRight) - LonMag)/dLon
             CoefR = 1 - CoefL
             Data_VII(iVar,i,:) = Ptr_II(iLeft,:)*CoefL + Ptr_II(iRight,:)*CoefR
          end do
       else 
          Data_VII(iVar,:,:) = Ptr_II(:,:)
       end if
    end do
    write(*,*)'RIM received = ', Data_VII(:, MinLon, MinLat)

    if(DoTest)then
       write(*,*)'SWMF_GridComp shape of Ptr =', shape(Ptr_II)
       ! Do not check the poles
       do j = MinLat, MaxLat; do i = MinLon, MaxLon
          ! Calculate exact solution
          Exact_V = FieldTest_V
          ! add time dependence for Hall field
          Exact_V(1) = Exact_V(1) + tCurrent*dHallPerDtTest

          LonMag = modulo(LonSM_I(i) + 180 - dPhiSm2Mag, 360.0) - 180
          ! add spatial dependence
          Exact_V = Exact_V + CoordCoefTest &
               *abs(LonMag)*(90-abs(LatSM_I(j)))

          ! Exact_V = LonMag
          ! *sin(Lon_I(i)*cDegToRad)*cos(Lat_I(j)*cDegToRad)
          if(abs(Data_VII(1,i,j) - Exact_V(1)) > 1e-10) &
               write(*,*) 'ERROR in SWMF_GridComp ', &
               'at i, j, Lon, Lat, Hall, Exact, Error=', &
               i, j, LonSM_I(i), LatSM_I(j), Data_VII(1,i,j), &
               Exact_V(1), Data_VII(1,i,j) - Exact_V(1)
          if(abs(Data_VII(2,i,j) - Exact_V(2)) > 1e-10) &
               write(*,*) 'ERROR in SWMF_GridComp ', &
               'at i, j, Lon, Lat, Pede, Exact, Error=', &
               i, j, LonSM_I(i), LatSM_I(j), Data_VII(2,i,j), &
               Exact_V(2), Data_VII(2,i,j) - Exact_V(2)
       end do; end do
       write(*,*)'SWMF_GridComp value of Data(MidLon,MidLat)=', &
            Data_VII(:,(MinLon+MaxLon)/2,(MinLat+MaxLat)/2)
    end if

    ! Clean memory
    deallocate(Data_VII, stat=iError)
    if(iError /= 0) call my_error('deallocate(Data_VII)')

    if(DoShiftDataCoupling) then 
       call get_sm_to_mag_angle(Clock, dPhiSm2Mag, iError)

       call update_export_state(gComp, iError)

    else 

       call ESMF_FieldGet(Field, grid=Grid, rc=iError)
       if(iError /= ESMF_SUCCESS) call my_error('ESMF_FieldGetGrid')

       call update_coordinates(Grid, Clock, LonSM_I, LatSM_I, .true., &
            dPhiSm2Mag=dPhiSm2Mag, iError=iError)
       if(iError /= ESMF_SUCCESS) call my_error('update_coordinates')
    end if

    call write_log("RIM_grid_comp:run routine returned")

    iError = ESMF_SUCCESS

  end subroutine my_run
  !============================================================================
  subroutine update_export_state(gComp, iError)
    type(ESMF_GridComp):: gComp
    integer, intent(out):: iError

    type(ESMF_State):: ExportState
    type(ESMF_Field):: Field
    type(ESMF_Grid):: Grid

    real(ESMF_KIND_R8), pointer :: Ptr_II(:,:), Lon_I(:), Lat_I(:)  
    real(ESMF_KIND_R8), allocatable :: Data_VII(:,:,:)
    integer:: i, j, iVar, itemCount
    character(len=4):: NameField

    integer:: iLeft, iRight
    real(ESMF_KIND_R8):: Coef, dLon, LonSM, CoefL, CoefR
    !--------------------------------------------------------------------------
    call write_log("RIM_grid_comp:update_export_state routine called")

    call NUOPC_ModelGet(gComp, exportState=ExportState, rc=iError)
    if(iError /= ESMF_SUCCESS) call my_error("NUOPC_ModelGet")

    allocate(Data_VII(nVarRim2Ipe,MinLon:MaxLon,MinLat:MaxLat))

    do iVar = 1, nVarRim2Ipe
       ! Get pointers to the variables in the export state
       nullify(Ptr_II)
       NameField = NameFieldRim2Ipe_V(iVar)

       ! Check field is available or not
       ! Export fields may not be defined in case of uni-directional coupling
       call ESMF_StateGet(ExportState, itemSearch=trim(NameField), &
            itemCount=itemCount, rc=iError)

       if (itemCount /= 0) then
          call ESMF_StateGet(ExportState, itemName=NameField, field=Field, &
               rc=iError)
          if(iError /= ESMF_SUCCESS) call my_error("ESMF_StateGet "//NameField)

          call ESMF_FieldGet(Field, farrayPtr=Ptr_II, rc=iError)
          if(iError /= ESMF_SUCCESS) call my_error("ESMF_FieldGet "//NameField)

          call ESMF_FieldGet(Field, grid=Grid, rc=iError)
          if(iError /= ESMF_SUCCESS) call my_error('ESMF_FieldGet Grid')

          call get_coords(Grid, Lon_I, Lat_I, iError)
          if(iError /= ESMF_SUCCESS) call my_error('get_coords')

          Coef = 10**iVar

          ! Add coordinate dependence
          ! With abs(Lon)*abs(Lat) dependence, the coupling does not pass
          ! correct values. To be investigated.
          do j = MinLat, MaxLat; do i = MinLon, MaxLon
             Data_VII(iVar,i,j) = abs(Lon_I(i))*abs(Lat_I(j))*Coef
          end do; end do

          if(DoShiftDataCoupling) then 
             do i = MinLon, MaxLon
                dLon = LonSM_I(2) - LonSM_I(1)

                ! For Ptr_II, its coordinates are in MAG system. Get the corresponding
                ! Lon in SM system.
                LonSM = modulo(LonSM_I(i) + 180 + dPhiSm2Mag, 360.0) - 180

                iLeft = floor((LonSM - LonSM_I(1))/dLon) + 1
                iRight = iLeft + 1
                CoefL = (LonSM_I(iRight) - LonSM)/dLon
                CoefR = 1.0 - CoefL
                Ptr_II(i,:) = Data_VII(iVar,iLeft,:)*CoefL + Data_VII(iVar,iRight,:)*CoefR
             end do
          else 
             Ptr_II(:,:) = Data_VII(iVar,:,:)
          end if


       end if
    end do ! iVar



    deallocate(Data_VII)

    iError = ESMF_SUCCESS
    call write_log("RIM_grid_comp:update_export_state routine returned")
    call ESMF_LogFlush()

  end subroutine update_export_state
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
end module RIM_grid_comp
!==============================================================================
