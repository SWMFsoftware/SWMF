!^CFG COPYRIGHT UM
!^CMP FILE IM
!==========================================================================
subroutine GM_put_from_im(Buffer_II,iSizeIn,jSizeIn,NameVar)

  !call stop_mpi('RCM is OFF') !^CFG UNCOMMENT IF NOT RCM
  !^CFG IF RCM BEGIN
  use CON_coupler
  use ModImPressure                              ! Storage for IM pressure
  use ModNumConst
  use ModMain, ONLY : n_step,time_simulation
  use ModIoUnit, ONLY: UNITTMP_
  implicit none
  CHARACTER (LEN=80) :: filename
  character(len=*), parameter :: NameSub='GM_put_from_im'

  integer, intent(in) :: iSizeIn,jSizeIn
  real, intent(in) :: Buffer_II(iSizeIn,jSizeIn)
  character(len=*), intent(in) :: NameVar
  integer :: nCells_D(2), iError, i,j
  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------

  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  if(NameVar /= 'p') &
       call CON_stop(NameSub//' invalid NameVar='//NameVar)

  nCells_D=ncells_decomposition_d(IM_)
  if( iSizeIn /= nCells_D(1) .or. jSizeIn /= nCells_D(2) ) then

     write(*,*)NameSub//' grid sizes do not agree iSize,jSize,nCells=',&
          iSizeIn,jSizeIn,nCells_D(1:2)
     call CON_stop(NameSub//' SWMF_ERROR')
  end if

  if(.not.allocated(RCM_lat))then
     ! Allocate RCM_lat, RCM_lon, RCM_p
     call im_pressure_init(iSizeIn, jSizeIn)
     ! Convert colat, lon to lat-lon in degrees and store
     RCM_lat = (cHalfPi - Grid_C(IM_) % Coord1_I) * cRadToDeg
     RCM_lon = Grid_C(IM_)% Coord2_I              * cRadToDeg
  end if

  ! Store IM variable for internal use
  RCM_p = Buffer_II

  ! Make sure that ray tracing has been done, so the pressure can be applied
  call ray_trace

  if(DoTest)call write_pressure_tec  ! TecPlot output
  if(DoTest)call write_pressure_idl  ! IDL     output

contains

  !============================================================================
  subroutine write_pressure_tec
    integer :: j2
    real :: lonShift
    !-------------------------------------------------------------------------

    !write values to plot file
    write(filename,'(a,i6.6,a)')"IMp_n=",n_step,".dat"
    OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown')
    write(UNITTMP_,'(a)') 'TITLE="Raytrace Values"'
    write(UNITTMP_,'(a)') 'VARIABLES="J", "I", "Lon", "Lat", "RCM pressure"'
    write(UNITTMP_,'(a,i4,a,i4,a)') &
         'ZONE T="IM Pressure", I=',jSizeIn+1,', J=',iSizeIn,', K=1, F=POINT'
    do i=1,iSizeIn
       do j2=1,jSizeIn+1
          j=j2; if(j2==jSizeIn+1) j=1
          lonShift=0.; if(j2==jSizeIn+1) lonShift=360.
          write(UNITTMP_,'(2i4,3G14.6)') j2,i,RCM_lon(j)+lonShift,RCM_lat(i), &
               RCM_p(i,j)
       end do
    end do
    CLOSE(UNITTMP_)

  end subroutine write_pressure_tec

  !============================================================================
  subroutine write_pressure_idl

    !write values to plot file

    write(filename,'(a,i6.6,a)')"IMp_n=",n_step,".out"

    OPEN (UNIT=UNITTMP_, FILE=filename, STATUS='unknown', &
         iostat =iError)
    if (iError /= 0) call CON_stop("Can not open raytrace File "//filename)
    write(UNITTMP_,'(a79)')            'Raytrace Values_var22'
    write(UNITTMP_,'(i7,1pe13.5,3i3)') n_step,time_simulation,2,1,1
    write(UNITTMP_,'(3i4)')            jSizeIn,iSizeIn
    write(UNITTMP_,'(100(1pe13.5))')   0.0
    write(UNITTMP_,'(a79)')            'Lon Lat p nothing'
    do i=iSizeIn,1,-1
       do j=1,jSizeIn
          write(UNITTMP_,'(100(1pe18.10))') &
               RCM_lon(j),RCM_lat(i),RCM_p(i,j)
       end do
    end do
    CLOSE(UNITTMP_)

  end subroutine write_pressure_idl

  !^CFG END RCM

end subroutine GM_put_from_im
