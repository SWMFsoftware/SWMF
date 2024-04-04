module IPE_grid_comp

  ! ESMF Framework module
  use ESMF
  use ESMFSWMF_variables, ONLY: add_fields, nVarEsmf, NameFieldEsmf_V, &
       DoTest, FieldTest_V, CoordCoefTest, dHallPerDtTest, iCoupleFreq, &
       write_log, write_error
  ! Conversion to radians

  implicit none
  private

  public:: set_services

  ! ESMF "dynamo" grid and domain. The real grid is not uniform in
  ! latitude, which will be implemented in the near future !!!
  ! This is a 2D spherical grid in MAG coordinates (rotates with Earth):
  ! +Z points to north magnetic dipole, +Y is towards rotational Omega x Z

  integer, parameter:: nLon=81, nLat=97 ! Default ESMF grid size

  ! Coordinate arrays
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
  subroutine my_init(gComp, ImportState, ExportState, Clock, iError)

    type(ESMF_GridComp):: gComp
    type(ESMF_State)   :: ImportState
    type(ESMF_State)   :: ExportState
    type(ESMF_Clock)   :: Clock
    integer, intent(out):: iError

    type(ESMF_Grid):: Grid
    type(ESMF_Field):: Field
    type(ESMF_VM):: Vm
    real(ESMF_KIND_R8), pointer :: Ptr_II(:,:)
    integer:: iVar, i, j
    character(len=4):: NameField
    !--------------------------------------------------------------------------
    call write_log("IPE_grid_comp init called")
    iError = ESMF_FAILURE

    ! Create Lon-Lat grid where -180<=Lon<=180-dLon, -90<=Lat<=90
    Grid = ESMF_GridCreateNoPeriDim(maxIndex=[nLon-1, nLat-1], &
         coordDep1=[1], coordDep2=[2], coordSys=ESMF_COORDSYS_CART, &
         name="dynamo grid", rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridCreateNoPeriDim')

    call ESMF_GridAddCoord(Grid, staggerloc=ESMF_STAGGERLOC_CORNER, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridAddCoord')

    nullify(Lon_I)
    call ESMF_GridGetCoord(Grid, CoordDim=1, &
         staggerLoc=ESMF_STAGGERLOC_CORNER, farrayPtr=Lon_I, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridGetCoord 1')
    write(*,*)'IPE size(Lon_I)=', size(Lon_I)

    call ESMF_GridGetCoord(Grid, CoordDim=2, &
         staggerLoc=ESMF_STAGGERLOC_CORNER, farrayPtr=Lat_I, rc=iError)
    if(iError /= ESMF_SUCCESS)call my_error('ESMF_GridGetCoord 2')

    write(*,*)'IPE size(Lat_I)=', size(Lat_I)
    ! Uniform longitude grid from -180 to 180
    do i = 1, nLon
       Lon_I(i) = (i-1)*(360.0/(nLon-1)) - 180
    end do
    write(*,*)'IPE grid: Lon_I(1,2,last)=', Lon_I([1, 2, nLon])
    ! Nonuniform latitude grid
    Lat_I = LatIpe_I
    write(*,*)'IPE grid: Lat_I(1,2,last)=', Lat_I([1, 2, nLat])

    ! Add fields to the export state
    call add_fields(Grid, ExportState, IsFromEsmf=.true., iError=iError)
    if(iError /= ESMF_SUCCESS) call my_error("add_fields")

    ! Initialize the data
    do iVar = 1, nVarEsmf
       ! Get pointers to the variables in the export state
       nullify(Ptr_II)
       NameField = NameFieldEsmf_V(iVar)
       call ESMF_StateGet(ExportState, itemName=NameField, field=Field, &
            rc=iError)
       if(iError /= ESMF_SUCCESS) call my_error("ESMF_StateGet "//NameField)

       call ESMF_FieldGet(Field, farrayPtr=Ptr_II, rc=iError)
       if(iError /= ESMF_SUCCESS) call my_error("ESMF_FieldGet "//NameField)

       if(iError /= ESMF_SUCCESS) RETURN
       select case(NameField)
       case('Hall')
          Ptr_II = FieldTest_V(1)
       case('Ped')
          Ptr_II = FieldTest_V(2)
       case default
          write(*,*)'ERROR in ESMF_GridComp:init: unknown NameField=',&
               NameField,' for iVar=',iVar
          iError = ESMF_FAILURE; RETURN
       end select

       ! Add coordinate dependence
       do j = 1, nLat; do i = 1, nLon
          Ptr_II(i,j) = Ptr_II(i,j) + CoordCoefTest &
               *abs(Lon_I(i))*(90-abs(Lat_I(j)))
       end do; end do

    end do ! iVar

    iError = ESMF_SUCCESS
    call write_log("IPE_grid_comp init returned")
    call ESMF_LogFlush()

  end subroutine my_init
  !============================================================================
  subroutine my_run(gComp, ImportState, ExportState, Clock, iError)

    type(ESMF_GridComp):: gComp
    type(ESMF_State)   :: ImportState
    type(ESMF_State)   :: ExportState
    type(ESMF_Clock)   :: Clock
    integer, intent(out):: iError

    type(ESMF_Field):: Field

    ! Access to the MHD data
    real(ESMF_KIND_R8), pointer :: Ptr_II(:,:)
    !--------------------------------------------------------------------------
    call write_log("IPE_grid_comp run called")
    iError = ESMF_FAILURE

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
       write(*,*)'IPE_grid_comp:run old Hall=', Ptr_II(nLon/2,nLat/2)
       Ptr_II = Ptr_II + iCoupleFreq*dHallPerdtTest
       write(*,*)'IPE_grid_comp:run new Hall=', Ptr_II(nLon/2,nLat/2)
    end if

    iError = ESMF_SUCCESS
    call write_log("IPE_grid_comp run returned")

  end subroutine my_run
  !============================================================================
  subroutine my_final(gComp, ImportState, ExportState, Clock, iError)

    type(ESMF_GridComp) :: gComp
    type(ESMF_State) :: ImportState
    type(ESMF_State) :: ExportState
    type(ESMF_Clock) :: Clock
    integer, intent(out):: iError

    !--------------------------------------------------------------------------
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
!==============================================================================
