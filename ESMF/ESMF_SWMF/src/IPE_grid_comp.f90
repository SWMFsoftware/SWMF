module IPE_grid_comp

  ! ESMF Framework module
  use ESMF
  use ESMFSWMF_variables, ONLY: add_fields, nVarEsmf, NameFieldEsmf_V, &
       DoTest, FieldTest_V, CoordCoefTest, dHallPerDtTest, iCoupleFreq, &
       write_log, write_error
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

  ! Coordinate arrays
  real(ESMF_KIND_R8), pointer:: Lon_I(:), Lat_I(:)

contains
  !============================================================================
  subroutine SetServices(gComp, rc)

    type(ESMF_GridComp) :: gComp
    integer, intent(out):: rc

    call ESMF_GridCompSetEntryPoint(gComp, ESMF_METHOD_INITIALIZE, &
         userRoutine=my_init, rc=rc)
    call ESMF_GridCompSetEntryPoint(gComp, ESMF_METHOD_RUN, &
         userRoutine=my_run, rc=rc)
    call ESMF_GridCompSetEntryPoint(gComp, ESMF_METHOD_FINALIZE, &
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
    integer                     :: iVar, i, j
    character(len=4):: NameField
    !-------------------------------------------------------------------------
    call write_log("IPE_grid_comp init called")
    rc = ESMF_FAILURE

    ! Create Lon-Lat grid where -180<=Lon<=180-dLon, -90<=Lat<=90
    Grid = ESMF_GridCreateNoPeriDim(maxIndex=[nLonEsmf-1, nLatEsmf-1], &
         coordDep1=[1], coordDep2=[2], coordSys=ESMF_COORDSYS_CART, &
         name="dynamo grid", rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_GridCreateNoPeriDim')

    call ESMF_GridAddCoord(Grid, staggerloc=ESMF_STAGGERLOC_CORNER, rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_GridAddCoord')

    nullify(Lon_I)
    call ESMF_GridGetCoord(Grid, CoordDim=1, &
         staggerLoc=ESMF_STAGGERLOC_CORNER, farrayPtr=Lon_I, rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_GridGetCoord 1')
    write(*,*)'ESMF_GridComp size(Lon_I)=', size(Lon_I)

    call ESMF_GridGetCoord(Grid, CoordDim=2, &
         staggerLoc=ESMF_STAGGERLOC_CORNER, farrayPtr=Lat_I, rc=rc)
    if(rc /= ESMF_SUCCESS)call my_error('ESMF_GridGetCoord 2')

    write(*,*)'ESMF_GridComp size(Lat_I)=', size(Lat_I)
    ! Uniform longitude grid
    do i = 1, nLonEsmf
       Lon_I(i) = (i-1)*(360.0/(nLonEsmf-1)) - 180
    end do
    write(*,*)'IPE grid: Lon_I(1,2,last)=', Lon_I([1,2,nLonEsmf])
    ! Uniform latitude grid (for now!!!)
    do i = 1, nLatEsmf
       Lat_I(i) = (i-1)*(180./(nLatEsmf-1)) - 90
    end do
    write(*,*)'IPE grid: Lat_I(1,2,last)=', Lat_I([1,2,nLatEsmf])

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
       do j = 1, nLatEsmf; do i = 1, nLonEsmf
          Ptr_II(i,j) = Ptr_II(i,j) + CoordCoefTest &
               *abs(Lon_I(i))*(90-abs(Lat_I(j)))
       end do; end do

    end do ! iVar

    rc = ESMF_SUCCESS
    call write_log("IPE_grid_comp init returned")
    call ESMF_LogFlush()

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
    call write_log("IPE_grid_comp run called")
    rc = ESMF_FAILURE

    ! We should execute the ESMF code here and put the result into
    ! the fields of the ExportState
    if(DoTest)then
       ! Get pointers to the MHD variables in the export state
       !!! This could be done in the initialization ?!
       nullify(Ptr_II)
       call ESMF_StateGet(ExportState, itemName='Hall', field=Field, rc=rc)
       if(rc /= ESMF_SUCCESS) call my_error("ESMF_StateGet for Hall")

       call ESMF_FieldGet(Field, farrayPtr=Ptr_II, rc=rc) 
       if(rc /= ESMF_SUCCESS) call my_error("ESMF_FieldGet for Hall")

       ! Update state by changing Hall conductivity
       write(*,*)'IPE_grid_comp:run old Hall=', Ptr_II(nLonEsmf/2,nLatEsmf/2)
       Ptr_II = Ptr_II + iCoupleFreq*dHallPerdtTest
       write(*,*)'IPE_grid_comp:run new Hall=', Ptr_II(nLonEsmf/2,nLatEsmf/2)
    end if

    rc = ESMF_SUCCESS
    call write_log("IPE_grid_comp run returned")

  end subroutine my_run
  !============================================================================
  subroutine my_final(gComp, ImportState, ExportState, Clock, rc)

    type(ESMF_GridComp) :: gComp
    type(ESMF_State) :: ImportState
    type(ESMF_State) :: ExportState
    type(ESMF_Clock) :: Clock
    integer, intent(out):: rc

    call write_log("IPE_grid_comp finalize called")
    call write_log("IPE_grid_comp finalize returned")

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
