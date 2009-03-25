program earth_traj

  ! earth_traj
  !
  ! Generation of the earth trajectoy file for input by SWMF
  ! as a thajectory of a satellite in the inner heliosphere.
  !
  ! The tragectory comprises the earth locations in HAE, HGI,
  ! or HGR coordinate system.
  !
  ! The input parameters are read from EARTH_TRAJ.in, and the result is
  ! written into STDOUT.
  ! 

  use CON_geopack_internal
  use ModUT
  use ModNumConst
  use ModTimeConvert,   ONLY: time_int_to_real, time_real_to_int
  use ModConst

  implicit none

  real :: R, msec
  real :: StartTimeRel=0.0, EndTimeRel=0.0, TimeStep=0.0
  real(Real8_) :: StartTime, EndTime, CurrentTime
  character(len=128), allocatable, dimension(:)  :: arg
  integer, allocatable, dimension(:) :: len_arg
  character(len=32) :: NameCommand
  character(len=32) :: NameFile
  character(len=3)  :: CoordSys
  character(len=32) :: cDate
  character(len=32) :: cTimeLoop
  integer, dimension(7) :: iStartTime, iCT
  integer :: i, j, k, iFile
  integer ::  ios, ic1, ic2
  integer iarg, errno, argc, oUnit
  R=cAU/Rsun

  !
  ! Command line argument processing
  !


  !
  ! Parse the command line arguments
  !
  NameFile='EARTH_TRAJ.in'
  iFile = 10
  open(iFile,file=NameFile,status='old')
  READPARAM: do 
     read(iFile,'(a)',end=9000,err=9000) NameCommand
     if(index(NameCommand,'#').ne.1)NameCommand=''
     select case(NameCommand)
     case("#COOR")
        read(iFile,*,end=9000,err=9000) CoordSys
     case("#INCLUDE")
        read(iFile,*,end=9000,err=9000) NameFile
        iFile = iFile + 1
        open(iFile,file=trim(NameFile),status='old')
     case("#DATE","#STARTTIME")
        do i = 1, 6
           read(iFile,*,end=9000,err=9000) iStartTime(i)
        end do
        read(iFile,*,end=9000,err=9000)R   !Auxiliary
        iStartTime(7) = 0
     case("#TIMELOOP")
        read(iFile,*,end=9000,err=9000) StartTimeRel, EndTimeRel, TimeStep
     case("#START","#END")
        if(iFile>10) then
           close(iFile)
           iFile = iFile-1
        else
           exit READPARAM
        end if
     case default
        if (len_trim(NameCommand) .ne. 0) then  ! i.e., just skip an empty line
           write(*,*) &
                '+++ Wrong command "'//trim(NameCommand)//'" encountered +++'
           stop
        end if
     end select
  end do READPARAM
  close(10)

  !
  ! Check correctness of arguments
  !
  if ((CoordSys .ne. 'HAE') .and. (CoordSys .ne. 'HGI') .and. (CoordSys .ne. 'HGR')) &
       call CON_stop('+++ Unknown coordinate system "'//trim(CoordSys)//'" specified +++')
  if (StartTimeRel > EndTimeRel) &
       call CON_stop('+++ StartTime MUST be less than EndTime +++')
  if (TimeStep <= cZero) &
       call CON_stop('+++ TimeStep MUST be positive +++')
  if (EndTimeRel > StartTimeRel + TimeStep*49999) &
       call CON_stop('+++ Output file cannot have more than 50000 lines +++')


  oUnit=6
  write(oUnit,'(a)') '#COOR'
  write(oUnit,'(a)') trim(CoordSys)
  write(oUnit,*)
  write(oUnit,'(a)') '#DATE'
  do i = 1, 7
     write(oUnit,*) iStartTime(i)
  end do
  write(oUnit,*)
  write(oUnit,'(a)') '#TIMELOOP'
  write(oUnit,'(3f12.3)') StartTimeRel, EndTimeRel, TimeStep
  write(oUnit,*)
  write(oUnit,'(a)') '#START'

  call time_int_to_real(iStartTime, StartTime)

  CurrentTime = StartTime + StartTimeRel
  EndTime = StartTime + EndTimeRel

  do while(CurrentTime <= EndTime)
     call time_real_to_int(CurrentTime, iCT)
     call CON_recalc(iCT(1), iCT(2), iCT(3), iCT(4), iCT(5), iCT(6))
     select case(CoordSys)
     case('HAE')
        write(oUnit,'(7i4,3f12.2)') iCT, Sun2EarthDirHae_D*SunEMBDistance*R
     case('HGI')
        write(oUnit,'(7i4,3f12.2)') iCT, Sun2EarthDirHgi_D*SunEMBDistance*R
     case('HGR')
        write(oUnit,'(7i4,3f12.2)') iCT,  matmul(HgrHgi_DD,Sun2EarthDirHgi_D)*SunEMBDistance*R
     end select
     CurrentTime = CurrentTime + TimeStep
  end do

  stop

9000 continue
  call CON_stop('+++ Error opening file "'//trim(NameFile)//'" +++')
9010 continue
  call CON_stop('+++ Error in date format +++')
9020 continue
  call CON_stop('+++ Error in timeloop format +++')

end program earth_traj
