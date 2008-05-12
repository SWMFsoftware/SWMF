subroutine readAMIEoutput(iBLK, iError)

  use ModAMIE_Interface

  implicit none

  integer, intent(out) :: iError
  integer, intent(in)  :: iBLK

  integer :: iTime
  integer :: nfields
  integer :: ntemp, iyr, imo, ida, ihr, imi
  integer :: i,j, iField, iPot_, iAveE_, iEFlux_
  real*4  :: swv,bx,by,bz,aei,ae,au,al,dsti,dst,hpi,sjh,pot
  real*8  :: rtime
  integer, dimension(7) :: itime_i

  real*4, allocatable, dimension(:,:,:) :: AllData

  integer, parameter :: nFieldsMax = 100
  character (len=30), dimension(nFieldsMax) :: Fields

  logical :: IsBinary

  iError = 0
  open(11, file=AMIE_FileName, status='old',form='UNFORMATTED',iostat=iError)
  if (iError.ne.0) then
     write(*,*) "Error opening file:", AMIE_FileName
     stop
  endif
  AMIE_nLats = 0
  IsBinary = .true.

  read(11,iostat=iError) AMIE_nlats,AMIE_nmlts,AMIE_ntimes
  if ((iError.ne.0).or.(AMIE_nlats.gt.100)) then
     write(*,*) "Error reading variables AMIE_nlats, AMIE_nmlts, AMIE_ntimes"
     IsBinary = .false.
  endif
  close(11)

  if (IsBinary) then
     open(11, file=AMIE_FileName, status='old',form='UNFORMATTED')
     read(11) AMIE_nlats,AMIE_nmlts,AMIE_ntimes
  else
     open(11, file=AMIE_FileName, status='old')
     read(11,*) AMIE_nlats,AMIE_nmlts,AMIE_ntimes
  endif

  if (allocated(AMIE_Lats)) deallocate(AMIE_Lats)
  allocate(AMIE_Lats(AMIE_nLats), stat=iError)
  if (iError /= 0) then
     write(*,*) "Error in allocating array AMIE_Lats in "
     stop
  endif

  if (allocated(AMIE_Mlts)) deallocate(AMIE_Mlts)
  allocate(AMIE_Mlts(AMIE_nMlts), stat=iError)
  if (iError /= 0) then
     write(*,*) "Error in allocating array Mlts in "
     stop
  endif

  if (.not.allocated(AMIE_Potential)) then

     allocate(AMIE_Potential(AMIE_nMlts,AMIE_nLats,AMIE_nTimes,2), stat=iError)
     if (iError /= 0) then
        write(*,*) "Error in allocating array AMIE_Potential in "
        stop
     endif

     allocate(AMIE_EFlux(AMIE_nMlts,AMIE_nLats,AMIE_nTimes,2), stat=iError)
     if (iError /= 0) then
        write(*,*) "Error in allocating array AMIE_EFlux in "
        stop
     endif

     allocate(AMIE_AveE(AMIE_nMlts,AMIE_nLats,AMIE_nTimes,2), stat=iError)
     if (iError /= 0) then
        write(*,*) "Error in allocating array AMIE_AveE in "
        stop
     endif

     allocate(AMIE_Value(AMIE_nMlts,AMIE_nLats,AMIE_nTimes,2), stat=iError)
     if (iError /= 0) then
        write(*,*) "Error in allocating array AMIE_Value in "
        stop
     endif

     allocate(AMIE_Time(AMIE_nTimes,2), stat=iError)
     if (iError /= 0) then
        write(*,*) "Error in allocating array AMIETimes in "
        stop
     endif

  endif

  if (IsBinary) then
     read(11) (AMIE_Lats(i),i=1,AMIE_nLats)
     read(11) (AMIE_Mlts(i),i=1,AMIE_nMlts)
     read(11) nFields
  else
     read(11,*) (AMIE_Lats(i),i=1,AMIE_nLats)
     read(11,*) (AMIE_Mlts(i),i=1,AMIE_nMlts)
     read(11,*) nFields
  endif

  AMIE_Lats = 90.0 - AMIE_Lats

  if (nFields > nFieldsMax) then
     write(*,*) "Maximum number of fields in AMIE is ",nFieldsMax
     stop
  endif

  allocate(AllData(AMIE_nMlts,AMIE_nLats,nFields), stat=iError)
  if (iError /= 0) then
     write(*,*) "Error in allocating array AllData in "
     stop
  endif

  AMIE_iDebugLevel = 2

  do iField=1,nfields
     if (IsBinary) then
        read(11) Fields(iField)
     else
        read(11,'(a)') Fields(iField)
     endif

     if (AMIE_iDebugLevel > 1) write(*,*) Fields(iField)

     if ((index(Fields(iField),"Potential") > 0).and. &
         (index(Fields(iField),"odel") < 1)) then
        iPot_ = iField
        if (AMIE_iDebugLevel > 1) write(*,*) "<--- Potential Found", iPot_
     endif

     if ((index(Fields(iField),"Mean Energy") > 0) .and. &
         (index(Fields(iField),"odel") < 1)) then
        iAveE_ = iField
        if (AMIE_iDebugLevel > 1) write(*,*) "<--- Mean Energy Found", iAveE_
     endif

     if ((index(Fields(iField),"Energy Flux") > 0) .and. &
         (index(Fields(iField),"odel") < 1)) then
        iEFlux_ = iField
        if (AMIE_iDebugLevel > 1) write(*,*) "<--- Energy Flux Found", iEFlux_
     endif

  enddo

  do iTime=1,AMIE_ntimes

     if (IsBinary) then

        read(11) ntemp,iyr,imo,ida,ihr,imi
        read(11) swv,bx,by,bz,aei,ae,au,al,dsti,dst,hpi,sjh,pot

        do iField=1,nfields
           read(11) ((AllData(j,i,iField),j=1,AMIE_nMlts),i=1,AMIE_nLats)
        enddo

     else

        read(11,*) ntemp,iyr,imo,ida,ihr,imi
        read(11,*) swv,bx,by,bz,aei,ae,au,al,dsti,dst,hpi,sjh,pot

        do iField=1,nfields
           read(11,*) ((AllData(j,i,iField),j=1,AMIE_nMlts),i=1,AMIE_nLats)
        enddo

     endif

     itime_i(1) = iyr
     itime_i(2) = imo
     itime_i(3) = ida
     itime_i(4) = ihr
     itime_i(5) = imi
     itime_i(6) = 0
     itime_i(7) = 0
     call time_int_to_real(itime_i,rtime)
     AMIE_Time(iTime,iBLK) = ihr*3600.0 + imi*60.0

     ! We need Potential to be in Volts
     !         AveE to be in keV
     !         EFlux to be in W/m2

     AMIE_Potential(:,:,iTime,iBLK) = AllData(:,:,iPot_)
     AMIE_AveE(:,:,iTime,iBLK)      = AllData(:,:,iAveE_)
     ! Need to convert from erg/cm2/s to W/m2
     AMIE_EFlux(:,:,iTime,iBLK)     = AllData(:,:,iEFlux_) * 1.0e-7 * 100.0 * 100.0

  enddo

  write(*,*) "min max AMIE_Time : ", minval(AMIE_Time), maxval(AMIE_Time)

  close(11)

  deallocate(AllData, stat=iError)

end subroutine readAMIEoutput
