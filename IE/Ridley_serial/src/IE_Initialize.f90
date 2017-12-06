!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine initialize_ie_ua_buffers(iOutputError)

  use ModErrors
  use ModIE_Interface
  use IE_ModMain, only: iDebugLevel
  use ModFiles
  use ModConst, only: cRadToDeg, cTwoPi
  use CON_coupler, ONLY : Grid_C, ncell_id, IE_

  implicit none

  integer, intent(out) :: iOutputError

  logical :: IsFound_EFieldModel
  logical :: IsFound_AuroralModel

  integer :: i, j, k, nCells_D(2)

  logical :: IsInitialized = .false.

  !------------------------------------------------------------------

  IE_NameOfEFieldModel = 'SPS'
  UseGridBasedIE = .true.

  if (IsInitialized) return

  IsInitialized = .true.

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

  IsFound_EFieldModel = .true.
  IsFound_AuroralModel = .true.

  IEi_HavenBLKs = 2 !!! IEt_Grid % Descriptor % iRootMapDim_D(1)

  nCells_D      = ncell_id(IE_) + 1
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
             90.0 - Grid_C(IE_) % Coord1_I(i) * cRadToDeg

        IEr3_HaveLats(j,i,2) = &
             - Grid_C(IE_) % Coord1_I(IEi_HavenLats-i+1) * cRadToDeg

        IEr3_HaveMLTs(j,i,:) = &
             mod(12 + Grid_C(IE_) % Coord2_I(j) * 24.0 /cTwoPi, 24.0)
        do k=1,2
           if (IEr3_HaveMLTs(j,i,k) < 0.1)  IEr3_HaveMLTs(j,i,k)=0.0
           if (IEr3_HaveMLTs(j,i,k) > 23.9) IEr3_HaveMLTs(j,i,k)=23.9
        enddo

     enddo
  enddo

contains

  subroutine allocate_arrays

    if (iDebugLevel > 1) then
       write(*,*) "=> IEi_HavenBLKs : ", IEi_HavenBLKs
       write(*,*) "=> IEi_HavenLats : ", IEi_HavenLats
       write(*,*) "=> IEi_HavenMLTs : ", IEi_HavenMLTs
    endif

    allocate(IEr3_HaveLats(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs))
    allocate(IEr3_HaveMlts(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs))
    allocate(IEr3_HavePotential(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs))
    allocate(IEr3_HaveEFlux(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs))
    allocate(IEr3_HaveAveE(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs))

  end subroutine allocate_arrays

end subroutine initialize_ie_ua_buffers
