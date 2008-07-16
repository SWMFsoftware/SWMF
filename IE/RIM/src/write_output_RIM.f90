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

  character (len=5)  :: cBlock
  character (len=14) :: cTime

  !------------------------------------------------------------------------

  variables = min_vars
  output_type = idl_type

  select case(plot_vars(ifile))
  case('min')
     variables = min_vars
  case('max')
     variables = all_vars
  case('uam')
     variables = uam_vars
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

  write(cBlock,'(a1,i4.4)') "b",iProc+1

  if (IsTimeAccurate) then
     write(cTime,'("t",3i2.2,"_",3i2.2)') mod(TimeArray(1),100), TimeArray(2:6)
  else
     write(cTime,'("i",i7.7)') nSolve
  endif

  ! Write out header

  if (iProc == 0) then

     open(unit=iUnit,status="unknown", &
          file = trim(NameOutputDir)//trim(cTime)//"_"//trim(plot_vars(ifile))//".header")

     write(iUnit,*) "BLOCKS"
     write(iUnit,"(I7,A)") 1, " nBlocksAlt"
     write(iUnit,"(I7,A)") 1, " nBlocksLat"
     write(iUnit,"(I7,A)") nProc, " nBlocksLon"
     write(iUnit,*) ""

     write(iUnit,*) "ITERATION"
     write(iUnit,"(I7,A)") nSolve, " iteration"
     write(iUnit,*) ""

     write(iUnit,*) "TIME"
     write(iUnit,"(I7,A)") TimeArray(1), " Year"
     write(iUnit,"(I7,A)") TimeArray(2), " Month"
     write(iUnit,"(I7,A)") TimeArray(3), " Day"
     write(iUnit,"(I7,A)") TimeArray(4), " Hour"
     write(iUnit,"(I7,A)") TimeArray(5), " Minute"
     write(iUnit,"(I7,A)") TimeArray(6), " Second"
     write(iUnit,"(I7,A)") TimeArray(7), " Millisecond"
     write(iUnit,*) ""
  
     write(iUnit,*) "VERSION"
     write(iUnit,*) Version
     write(iUnit,*) ""

     write(iUnit,*) "NUMERICAL VALUES"
     select case(plot_vars(ifile))
     case('min')
        write(iUnit,"(I7,A)") 12, " nVars"
     case('max')
        write(iUnit,"(I7,A)") 12, " nVars"
     case('uam')
        write(iUnit,"(I7,A)") 12, " nVars"
     case('aur')
        write(iUnit,"(I7,A)") 12, " nVars"
     end select
     write(iUnit,"(I7,A)")       1, " nAltitudes"
     write(iUnit,"(I7,A)") nLats  , " nLatitude"
     write(iUnit,"(I7,A)") nLons+2, " nLongitudes"
     write(iUnit,*) ""

     write(iUnit,*) "NGHOSTCELLS"
     write(iUnit,"(I7,A)") 0, " nAltitudes"
     write(iUnit,"(I7,A)") 0, " nLatitude"
     write(iUnit,"(I7,A)") 1, " nLongitudes"
     write(iUnit,*) ""

     write(iUnit,*) "VARIABLE LIST"
     write(iUnit,"(I7,A1,a)")  1, " ", "MagneticLocalTime"
     write(iUnit,"(I7,A1,a)")  2, " ", "Latitude"
     write(iUnit,"(I7,A1,a)")  3, " ", "Potential (kV)"
     write(iUnit,"(I7,A1,a)")  4, " ", "Jr (microA/m2)"
     write(iUnit,"(I7,A1,a)")  5, " ", "AveE (keV)"
     write(iUnit,"(I7,A1,a)")  6, " ", "EFlux (ergs/cm2/s)"
     write(iUnit,"(I7,A1,a)")  7, " ", "SigmaH (mhos)"
     write(iUnit,"(I7,A1,a)")  8, " ", "SigmaP (mhos)"
     write(iUnit,"(I7,A1,a)")  9, " ", "MHD-Rho (unknown)"
     write(iUnit,"(I7,A1,a)") 10, " ", "MHD-P (Pa)"
     write(iUnit,"(I7,A1,a)") 11, " ", "MHD-InvB (unknown)"
     write(iUnit,"(I7,A1,a)") 12, " ", "MHD-T (unknown)"

     write(iUnit,*) ""

     close(iUnit)

  endif

  open(unit=iUnit,status="unknown", form="unformatted", &
       file = trim(NameOutputDir)//trim(cTime)//"_"//trim(plot_vars(ifile))//"."//cBlock)

  do iLat=1,nLats
     do iLon=0,nLons+1
        write(iUnit) Longitude(iLon,iLat),Latitude(iLon,iLat),&
             Potential(iLon,iLat), Jr(iLon,iLat), &
             AveE(iLon,iLat), EFlux(iLon,iLat), &
             SigmaH(iLon,iLat), SigmaP(iLon,iLat), &
             OuterMagRho(iLon,iLat), OuterMagP(iLon,iLat), &
             OuterMagInvB(iLon,iLat), OuterMagT(iLon,iLat)
     end do
  end do

  close(iUnit)

end subroutine write_output_RIM
