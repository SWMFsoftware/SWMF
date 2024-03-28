module EsmfGridCompMod

  ! ESMF Framework module
  use ESMF
  use ESMF_SWMF_Mod, ONLY: add_fields, nVarEsmf, NameFieldEsmf_V, &
       DoTest, FieldTest_V, CoordCoefTest, dHallPerDtTest, iCoupleFreq
  ! Conversion to radians
  use ModNumConst, ONLY: cDegToRad

  implicit none
  private

  public:: SetServices

  ! ESMF "dynamo" grid and domain. The real grid is not uniform in
  ! latitude, which will be implemented in the near future !!!
  ! This is a 2D spherical grid in MAG coordinates (rotates with Earth):
  ! +Z points to north magnetic dipole, +Y is towards rotational Omega x Z
  
  integer, public:: nLonEsmf=81, nLatEsmf=97 ! Default ESMF grid size
  real(ESMF_KIND_R8), parameter:: LonMinEsmf = -180.! Min longitude
  real(ESMF_KIND_R8), parameter:: LonMaxEsmf = +180.0 ! Max longitude
  real(ESMF_KIND_R8), parameter:: LatMinEsmf = -90.0 ! Min latitude
  real(ESMF_KIND_R8), parameter:: LatMaxEsmf = +90.0 ! Max latitude

  ! Coordinate arrays
  real(ESMF_KIND_R8), pointer, save:: Lon_I(:), Lat_I(:)
  
  integer:: MinLon = 0, MaxLon = 0, MinLat = 0, MaxLat = 0

  ! Define which ESMF procesor contains the full grid
  integer:: PetMap(1,1,1) = 0

  ! Only true for the PET having the grid
  logical:: IsGridPet = .false.

contains
  !============================================================================
  subroutine SetServices(gcomp, rc)

    type(ESMF_GridComp) :: gcomp
    integer, intent(out):: rc

    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=my_init, rc=rc)
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
         userRoutine=my_run, rc=rc)
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_FINALIZE, &
         userRoutine=my_final, rc=rc)

  end subroutine SetServices
  !============================================================================
  subroutine my_init(gComp, importState, exportState, Clock, rc)

    type(ESMF_GridComp):: gComp
    type(ESMF_State)   :: importState
    type(ESMF_State)   :: exportState
    type(ESMF_Clock)   :: Clock
    integer, intent(out):: rc

    type(ESMF_Grid):: Grid
    type(ESMF_Field):: Field
    type(ESMF_VM):: Vm
    real(ESMF_KIND_R8), pointer :: Ptr_II(:,:)
    integer                     :: iVar, i, j, iPet
    character(len=4):: NameField
    !-------------------------------------------------------------------------
    call ESMF_LogWrite("ESMFGridComp init called", ESMF_LOGMSG_INFO)
    rc = ESMF_FAILURE

    ! Create Lon-Lat grid where -180<=Lon<=180-dLon, -90<=Lat<=90
    Grid = ESMF_GridCreate1PeriDimUfrm(maxIndex=[nLonEsmf-1, nLatEsmf-1], &
         minCornerCoord=[LonMinEsmf, LatMinEsmf], &
         maxCornerCoord=[LonMaxEsmf, LatMaxEsmf], &
         staggerLocList=[ESMF_STAGGERLOC_CORNER, ESMF_STAGGERLOC_CORNER], &
         petMap=PetMap, name="ESMF grid", rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_GridCreate1PeriDimUfrm')

    call ESMF_VMGetCurrent(Vm, rc=rc)
    if (rc /= ESMF_SUCCESS) call my_error('ESMF_VMGetCurrent')
    call ESMF_VMGet(Vm, localPet=iPet, rc=rc)
    if (rc /= ESMF_SUCCESS) call my_error('ESMF_VMGet')

    IsGridPet = iPet == 0
    
    if(IsGridPet)then
       nullify(Lon_I)
       call ESMF_GridGetCoord(Grid, coordDim=1, &
            staggerloc=ESMF_STAGGERLOC_CORNER, farrayPtr=Lon_I, rc=rc)
       if(rc /= ESMF_SUCCESS)call	my_error('ESMF_GridGetCoord 1')
       MinLon = lbound(Lon_I, DIM=1); MaxLon = ubound(Lon_I, DIM=1)
       write(*,'(a,2i4)')'ESMF grid: MinLon, MaxLon=', MinLon, MaxLon
       if(MaxLon > MinLon) write(*,'(a,3f8.2)') &
            'ESMF grid: Lon_I(Min,Min+1,Max)=', Lon_I([MinLon,MinLon+1,MaxLon])

       nullify(Lat_I)
       call ESMF_GridGetCoord(Grid, coordDim=2, &
            staggerloc=ESMF_STAGGERLOC_CORNER, farrayPtr=Lat_I, rc=rc)
       if(rc /= ESMF_SUCCESS)call	my_error('ESMF_GridGetCoord 2')
       MinLat = lbound(Lat_I, DIM=1); MaxLat = ubound(Lat_I, DIM=1)
       write(*,'(a,2i4)')'ESMF grid: MinLat, MaxLat=', MinLat, MaxLat
       if(MaxLat > MinLat) write(*,'(a,3f8.2)') &
            'ESMF grid: Lat_I(1,2,last)=', Lat_I([MinLat,MinLat+1,MaxLat])

       ! Add fields to the export state
       call add_fields(Grid, ExportState, IsFromEsmf=.true., rc=rc)
       if(rc /= ESMF_SUCCESS) call my_error("add_fields")

       ! Initialize the data
       do iVar = 1, nVarEsmf
          ! Get pointers to the variables in the export state
          nullify(Ptr_II)
          NameField = NameFieldEsmf_V(iVar)
          call ESMF_StateGet(ExportState, itemName=NameField, field=Field, &
               rc=rc)
          if(rc /= ESMF_SUCCESS) call my_error("ESMF_StateGet for "//NameField)
            
          call ESMF_FieldGet(Field, farrayPtr=Ptr_II, rc=rc) 
          if(rc /= ESMF_SUCCESS) call my_error("ESMF_FieldGet for "//NameField)

          if(rc /= ESMF_SUCCESS) RETURN
          select case(NameField)
          case('Hall')
             Ptr_II = FieldTest_V(1)
          case('Ped')
             Ptr_II = FieldTest_V(2)
          case default
             write(*,*)'ERROR in ESMF_GridComp:init: unknown NameField=',&
                  NameField,' for iVar=',iVar
             rc = ESMF_FAILURE; return
          end select

          ! Add coordinate dependence
          do j = MinLat, MaxLat; do i = MinLon, MaxLon
             Ptr_II(i,j) = Ptr_II(i,j) + CoordCoefTest &
                  *sin(Lon_I(i)*cDegToRad)*cos(Lat_I(j)*cDegToRad)
          end do; end do

       end do ! iVar
    end if ! IsGridPet
    rc = ESMF_SUCCESS
    call ESMF_LogWrite("ESMFGridComp init returned", ESMF_LOGMSG_INFO)

  end subroutine my_init
  !============================================================================
  subroutine my_run(gComp, ImportState, ExportState, Clock, rc)

    type(ESMF_GridComp):: gComp
    type(ESMF_State)   :: ImportState
    type(ESMF_State)   :: ExportState
    type(ESMF_Clock)   :: Clock
    integer, intent(out):: rc

    type(ESMF_Field):: Field
    
    ! Access to the MHD data
    real(ESMF_KIND_R8), pointer :: Ptr_II(:,:)
    !--------------------------------------------------------------------------
    call ESMF_LogWrite("ESMFGridComp run called", ESMF_LOGMSG_INFO)
    rc = ESMF_FAILURE

    ! We should execute the ESMF code here and put the result into
    ! the fields of the ExportState
    if(DoTest .and. IsGridPet)then
       ! Get pointers to the MHD variables in the export state
       !!! This could be done in the initialization ?!
       nullify(Ptr_II)
       call ESMF_StateGet(ExportState, itemName='Hall', field=Field, rc=rc)
       if(rc /= ESMF_SUCCESS) call my_error("ESMF_StateGet for Hall")
       call ESMF_FieldGet(Field, farrayPtr=Ptr_II, rc=rc) 
       if(rc /= ESMF_SUCCESS) call my_error("ESMF_FieldGet for Hall")
       ! Update state by changing Hall conductivity
       write(*,*)'ESMFGridComp:run old Hall=', &
            Ptr_II((MinLon+MaxLon)/2,(MinLat+MaxLat)/2)
       Ptr_II = Ptr_II + iCoupleFreq*dHallPerdtTest
       write(*,*)'ESMFGridComp:run new Hall=', &
            Ptr_II((MinLon+MaxLon)/2,(MinLat+MaxLat)/2)
    end if

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("ESMFGridComp run returned", ESMF_LOGMSG_INFO)

  end subroutine my_run
  !============================================================================
  subroutine my_final(gComp, ImportState, ExportState, Clock, rc)

    type(ESMF_GridComp) :: gComp
    type(ESMF_State) :: ImportState
    type(ESMF_State) :: ExportState
    type(ESMF_Clock) :: Clock
    integer, intent(out):: rc

    call ESMF_LogWrite("ESMFGridComp finalize called", ESMF_LOGMSG_INFO)
    call ESMF_LogWrite("ESMFGridComp finalize returned", ESMF_LOGMSG_INFO)

  end subroutine my_final
  !============================================================================
  subroutine my_error(String)

    ! Write out error message and stop
    
    character(len=*), intent(in) :: String
    !--------------------------------------------------------------------------
    write(*,*)'ERROR in EsmfGridCompMod: ',String
    
    call ESMF_finalize

  end subroutine my_error
  !============================================================================
end module EsmfGridCompMod
