subroutine write_output_RIM(iFile)
  use ModProcIE
  use ModRIM
  use ModIoRIM

  implicit none

  integer, intent(in) :: iFile

  character (len=100) :: Filename

  integer, parameter :: tec_type = 1
  integer, parameter :: idl_type = 2
  integer, parameter :: all_vars = 1
  integer, parameter :: min_vars = 2
  integer, parameter :: uam_vars = 3                     !^CFG  IF TIEGCM
  integer, parameter :: aur_vars = 4

  integer :: output_type, variables, year

  integer :: iLon, iLat
  character(len=*), parameter :: NameSub='IE_write_output'
  character(len=4) :: IO_ext
  character (len=23) :: textNandT
  !------------------------------------------------------------------------

  variables = min_vars
  output_type = idl_type

  select case(plot_vars(ifile))
  case('minimum')
     variables = min_vars
  case('maximum')
     variables = all_vars
  case('uam')                                   !^CFG  IF TIEGCM
     variables = uam_vars                       !^CFG  IF TIEGCM
  case('aur')
     variables = aur_vars
  end select

  select case(plot_form(ifile))
  case('idl')
     output_type = idl_type
  case('tec')
     output_type = tec_type
  case default
     call CON_stop(NameSub//' IE_ERROR invalid plot_form='//plot_form(ifile))
  end select

  if (output_type == idl_type) then
     IO_ext = ".idl"
  elseif (output_type == tec_type) then
     IO_ext = ".tec"
  endif

  if (.not.IsTimeAccurate) then

     write(Filename,'(a,i6.6,a,i2.2,a)')trim(NameOutputDir)//"in",nSolve,&
          "_b",iProc,IO_ext

     write(textNandT,'(a,i7.7)') "N=",nSolve

  else
     write(Filename,'(a,3i2.2,"_",3i2.2,"_",i3.3,"_b",i1,a)') &
          trim(NameOutputDir)//"it",&
          mod(TimeArray(1),100),TimeArray(2:7),&
          iProc,IO_ext

     write(textNandT,'(a,i7.7,a,i4.4,a,i2.2,a,i2.2)') &
          "N=",nSolve," T=", &
          TimeArray(4),":",TimeArray(5),":",TimeArray(6)

  end if

  open(unit=iUnit,file=Filename,status="unknown")

  do iLon=0,nLons+1
     do iLat=1,nLats

        write(iUnit,*) Longitude(iLon,iLat),Latitude(iLon,iLat),Potential(iLon,iLat)

     end do
  end do

  close(iUnit)
end subroutine write_output_RIM
