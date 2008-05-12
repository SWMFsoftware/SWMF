subroutine IE_Initialize(iOutputError)

  use ModErrors
  use ModIE_Interface
  use ModFiles

  implicit none

  integer, intent(out) :: iOutputError
  character (len=100)  :: inFileName
  integer              :: iError

  integer, parameter  :: South_ = 1
  integer, parameter  :: North_ = 2

  logical :: IsFound_EFieldModel

  iError = 0
  iOutputError = 0

  IsFound_EFieldModel = .false.

  call set_error_codes

  !\
  ! --------------------------------------------------------------------
  ! Electric Field Models
  ! --------------------------------------------------------------------
  !/

  if (iDebugLevel > 1) &
       write(*,*) "==> Efield Model : ",IE_NameOfEFieldModel

  if (iDebugLevel > 1) &
       write(*,*) "==> Model Directory : ",IE_NameOfModelDir

  if (index(IE_NameOfEFieldModel,'zero') > 0) then
     IsFound_EFieldModel = .true.
  endif

  if (index(IE_NameOfEFieldModel,'weimer96') > 0) then
     IsFound_EFieldModel = .true.
     call merge_str(IE_NameOfModelDir, weimer96_file)
     open(LunEField_,file=weimer96_file,status='old', iostat = iError)
     if (iError /= 0) then
        write(6,*) 'Error opening file :',weimer96_file
        iOutputError = ecFileNotFound_
     endif
     call ReadCoef96(LunEField_)
  endif

  if (index(IE_NameOfEFieldModel,'weimer01') > 0) then
     IsFound_EFieldModel = .true.
     call merge_str(IE_NameOfModelDir, weimer01_file)
     open(LunEField_,file=weimer01_file,status='old',&
          form='unformatted', iostat = iError)
     if (iError /= 0) then
        write(6,*) 'Error opening file :',weimer01_file
        iOutputError = ecFileNotFound_
     endif
     call ReadCoef01(LunEField_)
  endif

!  if (index(IE_NameOfEFieldModel,'samie') > 0) then
!     IsFound_EFieldModel = .true.
!     call merge_str(IE_NameOfModelDir, stat_amie_file)
!     open(LunEField_,file=stat_amie_file,status='old', iostat = iError)
!     if (iError /= 0) then
!        write(6,*) 'Error opening file :',stat_amie_file
!        iOutputError = ecFileNotFound_
!     endif
!     call read_amies(LunEField_)
!     close(LunEField_)
!  endif

  if (index(IE_NameOfEFieldModel,'millstone_hpi') > 0) then
     IsFound_EFieldModel = .true.
     call merge_str(IE_NameOfModelDir, millstone_hill_i_file)
     open(LunEField_,file=millstone_hill_i_file,status='old', iostat = iError)
     if (iError /= 0) then
        write(6,*) 'Error opening file :',millstone_hill_i_file
        iOutputError = ecFileNotFound_
     endif
     call mhinit(1, LunEField_, 1, iDebugLevel)
  endif

  if (index(IE_NameOfEFieldModel,'millstone_imf') > 0) then
     IsFound_EFieldModel = .true.
     call merge_str(IE_NameOfModelDir, millstone_hill_s_file)
     open(LunEField_,file=millstone_hill_s_file,status='old', iostat = iError)
     if (iError /= 0) then
        write(6,*) 'Error opening file :',millstone_hill_s_file
        iOutputError = ecFileNotFound_
     endif
     call mhinit(2, LunEField_, 1, iDebugLevel)
  endif

  if (index(IE_NameOfEFieldModel,'hmr89') > 0) then
     IsFound_EFieldModel = .true.
     call merge_str(IE_NameOfModelDir, hepner_maynard_file)
     open(LunEField_,file=hepner_maynard_file,status='old', iostat = iError)
     if (iError /= 0) then
        write(6,*) 'Error opening file :',hepner_maynard_file
        iOutputError = ecFileNotFound_
     endif
     call gethmr(LunEField_)
  endif

  if (index(IE_NameOfEFieldModel,'izmem') > 0) then
     IsFound_EFieldModel = .true.
     call merge_str(IE_NameOfModelDir, izmem_file)
     open(LunEField_,file=izmem_file,status='old', iostat = iError)
     if (iError /= 0) then
        write(6,*) 'Error opening file :',izmem_file
        iOutputError = ecFileNotFound_
     endif
     call izinit(LunEField_)
  endif

  !\
  ! --------------------------------------------------------------------
  ! Conductance Models
  ! --------------------------------------------------------------------
  !/

  if (iDebugLevel > 1) &
       write(*,*) "==> Conductance Model : ",IE_NameOfAuroralModel

  if (index(IE_NameOfAuroralModel,'ihp') > 0) call read_conductance_model
  if (index(IE_NameOfAuroralModel,'pem') > 0) call read_conductance_model

  if (index(IE_NameOfEFieldModel,'amie') > 0) then

     UseGridBasedIE = .true.

     IsFound_EFieldModel = .true.

     call AMIE_SetFileName(AMIEFileNorth)
     call readAMIEOutput(North_, iError)

     call AMIE_SetFileName(AMIEFileSouth)
     call readAMIEOutput(South_, iError)

     call AMIE_GetnLats(IEi_HavenLats)
     call AMIE_GetnMLTs(IEi_HavenMLTs)
     IEi_HavenBLKs = 2

     if (iDebugLevel > 1) then
        write(*,*) "=> IEi_HavenBLKs : ", IEi_HavenBLKs
        write(*,*) "=> IEi_HavenLats : ", IEi_HavenLats
        write(*,*) "=> IEi_HavenMLTs : ", IEi_HavenMLTs
     endif

     allocate(IEr3_HaveLats(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs), &
          stat=iError)
     if (iError /= 0) then
        write(*,*) "Error in allocating array IEr3_HaveLats in Interface"
        stop
     endif

     allocate(IEr3_HaveMlts(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs), &
          stat=iError)
     if (iError /= 0) then
        write(*,*) "Error in allocating array IEr3_HaveMlts in Interface"
        stop
     endif

     allocate(IEr3_HavePotential(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs), &
          stat=iError)
     if (iError /= 0) then
        write(*,*) "Error in allocating array IEr3_HavePotential in Interface"
        stop
     endif

     allocate(IEr3_HaveEFlux(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs), &
          stat=iError)
     if (iError /= 0) then
        write(*,*) "Error in allocating array IEr3_HaveEFlux in Interface"
        stop
     endif

     allocate(IEr3_HaveAveE(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs), &
          stat=iError)
     if (iError /= 0) then
        write(*,*) "Error in allocating array IEr3_HaveAveE in Interface"
        stop
     endif

     call AMIE_GetLats(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs,&
          IEr3_HaveLats,iError)

     call AMIE_GetMLTs(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs,&
          IEr3_HaveMLTs,iError)

     return

  endif

  if (.not.IsFound_EFieldModel) then
     iOutputError = ecEFieldModelNotFound_
  endif

end subroutine IE_Initialize
