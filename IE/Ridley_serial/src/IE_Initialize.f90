subroutine initialize_ie_ua_buffers(iOutputError)

  use ModErrors
  use ModIE_Interface
  use IE_ModMain, only: iDebugLevel
  use ModFiles
  use ModConst
  use CON_coupler, ONLY : Grid_C, ncells_decomposition_d, IE_

  implicit none

  integer, intent(out) :: iOutputError
  character (len=100)  :: inFileName
  integer              :: iError

  integer, parameter  :: South_ = 1
  integer, parameter  :: North_ = 2

  logical :: IsFound_EFieldModel
  logical :: IsFound_AuroralModel

  integer :: i, j, k, nCells_D(2)

  logical :: IsInitialized = .false.

  !------------------------------------------------------------------

  IE_NameOfEFieldModel = 'SPS'
  UseGridBasedIE = .true.

  if (IsInitialized) return

  IsInitialized = .true.

  iError = 0
  iOutputError = 0

  IsFound_EFieldModel  = .false.
  IsFound_AuroralModel = .false.

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

  if (index(IE_NameOfAuroralModel,'ihp') > 0) then
       call read_conductance_model(iOutputError)
       IsFound_AuroralModel = .true.
    endif
  if (index(IE_NameOfAuroralModel,'pem') > 0) then
     call read_conductance_model(iOutputError)
     IsFound_AuroralModel = .true.
  endif

  !\
  ! This is for the serial potential solver
  !/

  if (index(IE_NameOfEFieldModel,'SPS') > 0) then

     IsFound_EFieldModel = .true.
     IsFound_AuroralModel = .true.

     IEi_HavenBLKs = 2 !!! IEt_Grid % Descriptor % iRootMapDim_D(1)

     nCells_D      = ncells_decomposition_d(IE_) + 1
     IEi_HavenLats = nCells_D(1)
     IEi_HavenMlts = nCells_D(2)

     call allocate_arrays

     !\
     ! This is sort of screwed up, since other IE models are mlt,lat
     ! based, where the SPS is theta,phi based.  So, we have to convert
     ! between theta and lat and phi and mlt.  Also, we have to switch the
     ! order of the variables.
     !/

     do i = 1, IEi_HavenLats
        do j = 1, IEi_HavenMlts

           IEr3_HaveLats(j,i,1) = &
                90.0 - Grid_C(IE_) % Coord1_I(i) * 180.0 / cPi

           IEr3_HaveLats(j,i,2) = &
                - Grid_C(IE_) % Coord1_I(IEi_HavenLats-i+1) * 180.0 / cPi

           IEr3_HaveMLTs(j,i,:) = &
                mod(12 + Grid_C(IE_) % Coord2_I(j) * 24.0 / (cTwo * cPi),24.0)
           do k=1,2
              if (IEr3_HaveMLTs(j,i,k) < 0.1)  IEr3_HaveMLTs(j,i,k)=0.0
              if (IEr3_HaveMLTs(j,i,k) > 23.9) IEr3_HaveMLTs(j,i,k)=23.9
           enddo

        enddo
     enddo

     return

  endif

  if (index(IE_NameOfEFieldModel,'amie') > 0) then

     IsFound_EFieldModel = .true.
     IsFound_AuroralModel = .true.

     call AMIE_SetFileName(AMIEFileNorth)
     call readAMIEOutput(North_, iError)

     call AMIE_SetFileName(AMIEFileSouth)
     call readAMIEOutput(South_, iError)

     call AMIE_GetnLats(IEi_HavenLats)
     call AMIE_GetnMLTs(IEi_HavenMLTs)
     IEi_HavenBLKs = 2

     call allocate_arrays

     call AMIE_GetLats(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs,&
          IEr3_HaveLats,iError)

     call AMIE_GetMLTs(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs,&
          IEr3_HaveMLTs,iError)

     return

  endif

  if (.not.IsFound_EFieldModel) then
     iOutputError = ecEFieldModelNotFound_
     return
  endif

  if (.not.IsFound_AuroralModel) then
     iOutputError = ecAuroralModelNotFound_
     return
  endif

contains

  subroutine allocate_arrays

    if (iDebugLevel > 1) then
       write(*,*) "=> IEi_HavenBLKs : ", IEi_HavenBLKs
       write(*,*) "=> IEi_HavenLats : ", IEi_HavenLats
       write(*,*) "=> IEi_HavenMLTs : ", IEi_HavenMLTs
    endif

    allocate(IEr3_HaveLats(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs), &
         stat=iError)
    if (iError /= 0) then
       call CON_stop("Error in allocating IEr3_HaveLats in Interface")
    endif

    allocate(IEr3_HaveMlts(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs), &
         stat=iError)
    if (iError /= 0) then
       call CON_stop("Error in allocating IEr3_HaveMlts in Interface")
    endif

    allocate(IEr3_HavePotential(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs), &
         stat=iError)
    if (iError /= 0) then
       call CON_stop("Error in allocating IEr3_HavePotential in Interface")
    endif

    allocate(IEr3_HaveEFlux(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs), &
         stat=iError)
    if (iError /= 0) then
       call CON_stop("Error in allocating IEr3_HaveEFlux in Interface")
    endif

    allocate(IEr3_HaveAveE(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs), &
         stat=iError)
    if (iError /= 0) then
       call CON_stop("Error in allocating IEr3_HaveAveE in Interface")
    endif

  end subroutine allocate_arrays

end subroutine initialize_ie_ua_buffers
