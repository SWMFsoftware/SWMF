!^CFG COPYRIGHT UM
!==============================================================================
module ModSatelliteFile
  use ModIO,     ONLY: MaxSatelliteFile
  implicit none
  save
  private !Except

  logical :: IsOpen_I(MaxSatelliteFile) = .false.
  character(LEN=80),dimension(MaxSatelliteFile) :: NameFile_I

  public:: set_satellite_file_status !Open, open to append or close the file
contains
  !============================================================================
  subroutine set_satellite_file_status(iSat,TypeStatus)
    use ModMain,   ONLY: n_step
    use ModIoUnit, ONLY: io_unit_new
    use ModIO,     ONLY: Satellite_name, &
         NamePlotDir, iUnitSat_I

    integer, intent(in) :: iSat
    character(LEN=*),intent(in) :: TypeStatus

    integer :: l1, l2
    logical :: oktest, oktest_me
    !-----------------------------------------
    select case(TypeStatus)
    case('open')
       call set_oktest('open_satellite_output_files', oktest, oktest_me)

       l1 = index(Satellite_name(iSat), '/', back=.true.) + 1
       l2 = index(Satellite_name(iSat), '.') - 1
       if (l1-1<=0) l1=1
       if (l2+1<=0) l2=len_trim(Satellite_name(iSat))

       write(NameFile_I(iSat),'(a,i6.6,a)')trim(NamePlotDir)//&
         'sat_'//Satellite_Name(iSat)(l1:l2)//'_n',n_step,'.sat'

       if(oktest) then
          write(*,*) 'open_satellite_output_files: satellitename:', &
               Satellite_name(iSat), 'status ='//TypeStatus
          write(*,*) 'iSat,l1,l2: ',iSat,l1,l2
          write(*,*) 'open_satellite_output_files: NameFile_I(iSat):', NameFile_I(iSat)
       end if

       iUnitSat_I(iSat)=io_unit_new()
       open(iUnitSat_I(iSat), file=trim(NameFile_I(iSat)), status='replace')

       IsOpen_I(iSat) = .true.
    case('append')
       if(.not.IsOpen_I(iSat))then
          iUnitSat_I(iSat)=io_unit_new()
       open(iUnitSat_I(iSat), file=trim(NameFile_I(iSat)), status='old',POSITION='append')
          IsOpen_I(iSat) = .true.
       end if
    case('close')
       if(IsOpen_I(iSat))close(iUnitSat_I(iSat))
       IsOpen_I(iSat) = .false.
    case default
       call stop_mpi(&
            'Unknown TypeStatus='//TypeStatus//' in set_satellite_file_status')
    end select
  end subroutine set_satellite_file_status
  !=============================================================================
  subroutine open_satellite_output_files(iSat,TypeStatus)

   
    implicit none
    integer, intent(in) :: iSat
    character(LEN=*),intent(in) :: TypeStatus
   

    !--------------------------------------------------------------------------
   

   

  end subroutine open_satellite_output_files

  !============================================================================

  subroutine close_satellite_output_files(iSat)

    use ModIO, ONLY : nSatellite,iUnitSat_I
    implicit none

    integer, intent(in) :: iSat

    

  end subroutine close_satellite_output_files
  !============================================================================
end module ModSatelliteFile
