module ESMF_grid_comp

  ! Code for the ESMFSWMF Gridded Component which creates 4 child Components:
  ! IPE, SWMF, IE and IPE_ie_coupler.

  ! ESMF Framework module
  use ESMF
  use NUOPC

  use NUOPC_Driver, &
      Driver_routine_SS             => SetServices, &
      Driver_label_SetModelServices => label_SetModelServices, &
      Driver_label_SetRunSequence   => label_SetRunSequence, &
      Driver_label_SetRunClock      => label_SetRunClock
  use NUOPC_Driver, only: NUOPC_DriverAddComp

  ! Various variables
  use ESMFSWMF_variables, ONLY: write_log, write_error, NameParamFile

  ! User Component registration routines
  use SWMF_grid_comp, ONLY: swmf_set_services    => set_services
  use IPE_grid_comp,  ONLY: ipe_set_services     => set_services
  use RIM_grid_comp,  ONLY: rim_set_services     => set_services

  implicit none
  private

  public:: ESMF_set_services

  type(ESMF_GridComp), save :: IpeComp, SwmfComp, RimComp

contains
  !============================================================================
  subroutine ESMF_set_services(gComp, iError)

    type(ESMF_GridComp) :: gComp
    integer, intent(out):: iError

    type(ESMF_Config)   :: config
    !--------------------------------------------------------------------------
    call NUOPC_CompDerive(gComp, Driver_routine_SS, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('NUOPC_CompDerive')
    call NUOPC_CompSpecialize(gComp, specLabel=Driver_label_SetModelServices, &
         specRoutine=set_model_services, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('NUOPC_CompSpecialize - init')
    call NUOPC_CompSpecialize(gComp, specLabel=Driver_label_SetRunSequence, &
         specRoutine=set_run_sequence, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('NUOPC_CompSpecialize - set run seq.')
    call ESMF_MethodRemove(gComp, Driver_label_SetRunClock, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_MethodRemove')
    call NUOPC_CompSpecialize(gComp, specLabel=Driver_label_SetRunClock, &
         specRoutine=NUOPC_NoOp, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('NUOPC_CompSpecialize - set run clock')
    config = ESMF_ConfigCreate(rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_ConfigCreate')
    call ESMF_ConfigLoadFile(config, trim(NameParamFile), rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_ConfigLoadFile')
    call ESMF_GridCompSet(gComp, config=config, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCompSet')
    call NUOPC_FieldDictionarySetup('fd_swmf.yaml', rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('NUOPC_FieldDictionarySetup')

  end subroutine ESMF_set_services
  !============================================================================
  subroutine set_run_sequence(gComp, iError)

    type(ESMF_GridComp) :: gComp
    integer, intent(out):: iError

    type(ESMF_Config) :: config
    type(NUOPC_FreeFormat) :: runSeqFF
    !--------------------------------------------------------------------------

    ! Read free format run sequence from config
    call ESMF_GridCompGet(gComp, config=config, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCompGet')
    runSeqFF = NUOPC_FreeFormatCreate(config, label="runSeq::", rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('NUOPC_FreeFormatCreate')

    ! Ingest FreeFormat run sequence
    call NUOPC_DriverIngestRunSequence(gComp, runSeqFF, &
         autoAddConnectors=.true., rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('NUOPC_DriverIngestRunSequence')
    call NUOPC_DriverPrint(gComp, orderflag=.true., rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('NUOPC_DriverPrint')

  end subroutine set_run_sequence
  !============================================================================
  subroutine set_model_services(gComp, iError)

    type(ESMF_GridComp) :: gComp
    integer, intent(out):: iError

    type(ESMF_Config) :: config
    type(NUOPC_FreeFormat) :: attrFF
    integer :: i, j, petCount, compCount
    integer :: petListBounds(2)
    integer, allocatable :: petList(:)
    character(len=32) :: model, prefix
    character(len=32), allocatable :: compLabels(:)
    !--------------------------------------------------------------------------
    call write_log("ESMF_gric_comp set_model_services called")
    iError = ESMF_FAILURE

    ! Query gridded component
    call ESMF_GridCompGet(gComp, petCount=petCount, config=config, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCompGet')

    ! Read and ingest free format driver attributes
    attrFF = NUOPC_FreeFormatCreate(config, label='DRV_attributes::', &
             relaxedflag=.true., rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('NUOPC_FreeFormatCreate')
    call NUOPC_CompAttributeIngest(gComp, attrFF, addFlag=.true., rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('NUOPC_CompAttributeIngest')
    call NUOPC_FreeFormatDestroy(attrFF, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('NUOPC_FreeFormatDestroy')

    ! Determine the generic component labels
    compCount = ESMF_ConfigGetLen(config, label='DRV_component_list:', rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_ConfigGetLen')

    allocate(compLabels(compCount))
    call ESMF_ConfigGetAttribute(config, valueList=compLabels, &
         label="DRV_component_list:", count=compCount, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_ConfigGetAttribute - drv')

    ! Determine information for each component and add to the driver
    do i = 1, compCount
       ! Construct component prefix
       prefix = trim(compLabels(i))
       ! Read in petList bounds
       call ESMF_ConfigGetAttribute(config, petListBounds, &
            label=trim(prefix)//"_petlist_bounds:", default=-1, rc=iError)
       if(iError /= ESMF_SUCCESS)call my_error('ESMF_ConfigGetAttribute - '//trim(prefix))
       ! Handle negative values
       if (petListBounds(1) < 0) petListBounds(1) = petCount + petListBounds(1)
       if (petListBounds(2) < 0) petListBounds(2) = petCount + petListBounds(2)
       petListBounds = min(petCount-1, petListBounds)
       ! Read in model instance name
       call ESMF_ConfigGetAttribute(config, model, &
            label=trim(prefix)//"_model:", default="none", rc=iError)
       if(iError /= ESMF_SUCCESS)call my_error('ESMF_ConfigGetAttribute - '//trim(model))
       ! Set petList for this component
       allocate(petList(petListBounds(2)-petListBounds(1)+1))
       do j = petListBounds(1), petListBounds(2)
          petList(j-petListBounds(1)+1) = j
       end do
       ! Add component/s
       ! SWMF
       if (trim(model) == 'swmf') then
          call NUOPC_DriverAddComp(gComp, trim(prefix), &
               swmf_set_services, petlist=petList, comp=SwmfComp, rc=iError)
          if(iError /= ESMF_SUCCESS)call my_error('NUOPC_DriverAddComp - SWMF')
       end if
       ! IPE
       if (trim(model) == 'ipe') then
          call NUOPC_DriverAddComp(gComp, trim(prefix), &
               ipe_set_services, petlist=petList, comp=IpeComp, rc=iError)
          if(iError /= ESMF_SUCCESS)call my_error('NUOPC_DriverAddComp - IPE')
       end if
       ! RIM
       if (trim(model) == 'rim') then
          call NUOPC_DriverAddComp(gComp, trim(prefix), &
               rim_set_services, petlist=petList, comp=RimComp, rc=iError)
          if(iError /= ESMF_SUCCESS)call my_error('NUOPC_DriverAddComp - RIM')
       end if
       ! Clear memory
       deallocate(petList)
    end do

    iError = ESMF_SUCCESS
    call write_log("ESMF_grid_comp set_model_services finished")

  end subroutine set_model_services
  !============================================================================
  subroutine my_error(String)

    character(len=*), intent(in) :: String
    !--------------------------------------------------------------------------
    call write_error("ESMF_grid_comp "//String)

  end subroutine my_error
  !============================================================================
end module ESMF_grid_comp
!==============================================================================
