program Interface

  use ModEIE_Interface

  implicit none

  character (len=100) :: inFileName
  integer             :: iError
  real*8              :: rTime
  integer, dimension(7) :: iTime_i
  real, dimension(4,2) :: templat, tempmlt, temppot

  iDebugLevel = 100

  write(6,*) 'Enter file name :'
  read(5,'(A100)') inFileName

  iError = 0

  call AMIE_SetFileName(inFileName)

  call readAMIEOutput(iError)

  call AMIE_GetnLats(IEi_HavenLats)
  call AMIE_GetnMLTs(IEi_HavenMLTs)
  IEi_HavenBLKs = 2

  if (iDebugLevel > 1) then
     write(*,*) "IEi_HavenBLKs : ", IEi_HavenBLKs
     write(*,*) "IEi_HavenLats : ", IEi_HavenLats
     write(*,*) "IEi_HavenMLTs : ", IEi_HavenMLTs
  endif

  allocate(IEr3_HaveLats(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs), &
       stat=iError)
  if (iError /= 0) then
     write(*,*) "Error in allocating array IEr3_HaveLats in Interface"
     call CON_stop('ERROR in Interface.f90')
  endif

  allocate(IEr3_HaveMlts(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs), &
       stat=iError)
  if (iError /= 0) then
     write(*,*) "Error in allocating array IEr3_HaveMlts in Interface"
     call CON_stop('ERROR in Interface.f90')
  endif

  allocate(IEr3_HavePotential(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs), &
       stat=iError)
  if (iError /= 0) then
     write(*,*) "Error in allocating array IEr3_HavePotential in Interface"
     call CON_stop('ERROR in Interface.f90')
  endif

  allocate(IEr3_HaveEFlux(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs), &
       stat=iError)
  if (iError /= 0) then
     write(*,*) "Error in allocating array IEr3_HaveEFlux in Interface"
     call CON_stop('ERROR in Interface.f90')
  endif

  allocate(IEr3_HaveAveE(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs), &
       stat=iError)
  if (iError /= 0) then
     write(*,*) "Error in allocating array IEr3_HaveAveE in Interface"
     call CON_stop('ERROR in Interface.f90')
  endif

  call AMIE_GetLats(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs,&
       IEr3_HaveLats,iError)

  call AMIE_GetMLTs(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs,&
       IEr3_HaveMLTs,iError)

  itime_i(1) = 1998
  itime_i(2) = 05
  itime_i(3) = 01
  itime_i(4) = 12
  itime_i(5) = 00
  itime_i(6) = 0
  itime_i(7) = 0
  call time_int_to_real(itime_i,rtime)

  call AMIE_GetPotential(rtime, IE_Interpolate_, &
       IEi_HavenMlts, IEi_HavenLats, IEi_HavenBLKs, IEr3_HavePotential, iError)

  call AMIE_GetAveE(rtime, IE_Closest_, &
       IEi_HavenMlts, IEi_HavenLats, IEi_HavenBLKs, IEr3_HaveAveE, iError)

  call AMIE_GetEFlux(rtime, IE_Closest_, &
       IEi_HavenMlts, IEi_HavenLats, IEi_HavenBLKs, IEr3_HaveEFlux, iError)

  call IO_SetnMLTs(4)
  call IO_SetnLats(2)

  templat(:,1) = -60.0
  templat(:,2) = 70.0
  tempmlt(1,:) =  0.0
  tempmlt(2,:) =  6.0
  tempmlt(3,:) = 12.0
  tempmlt(4,:) = 18.0

  call IO_SetGrid(tempmlt, templat, iError)

  call IO_GetPotential(temppot,iError)

  write(*,*) templat(:,2)
  write(*,*) tempmlt(:,2)
  write(*,*) temppot(:,2)

  call IE_End

end program Interface
