!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!BOP =========================================================================
!ROUTINE: XY_get_for_pm - get XY data for sending it to PM
!INTERFACE:
subroutine XY_get_for_pm(Buffer_II, iSize, jSize, NameVar)

  !USES:
  use XY_ModField, ONLY: XY_calc_field, XY_Field_II, XY_FieldUnit

  implicit none

  !INPUT ARGUMENTS:
  integer, intent(in)           :: iSize,jSize            ! Size of buffer
  character (len=*),intent(in)  :: NameVar                ! Name of data

  !OUTPUT ARGUMENTS:
  real, intent(out)             :: Buffer_II(iSize,jSize) ! Data buffer

  !LOCAL VARIABLES:
  character (len=*),parameter :: NameSub='XY_get_for_pm'

  !DESCRIPTION:
  ! This subroutine shows an example for getting data for another component.
  ! This subroutine should be in the srcInterface directory of the
  ! specific XY component version. The subroutine is called by the
  ! couple\_xy\_pm method of the CON\_couple\_pm\_xy module in 
  ! CON/Interface/src.
  !
  ! Calculate field for the positions needed by the PM component.
  ! We assume here that the name of the variable, its dimension,
  ! and the grid descriptor of the PM component are 
  ! sufficient for extracting the data.
  !EOP
  !---------------------------------------------------------------------------
  !BOC
  ! Allocate and set XY_Field_II for the positions needed by the PM component
  call XY_calc_field(NameVar, iSize, jSize)

  ! Put appropriate part of the XY_Field_II array into the buffer
  ! Convert to SI units
  select case(NameVar)
  case('FieldNorth')
     Buffer_II = XY_FieldUnit * XY_Field_II(1:iSize,:)
  case('FieldSouth')
     Buffer_II = XY_FieldUnit * XY_Field_II(iSize+1:2*iSize)
  case default
     call CON_stop(NameSub//' invalid NameVar='//NameVar)
  end select
  !EOC
end subroutine XY_get_for_pm

!BOP =========================================================================
!ROUTINE: PM_put_from_xy - put into PM the data received from XY
!INTERFACE:
subroutine PM_put_from_xy(Buffer_II, iSize, jSize, NameVar)

  !USES:
  use PM_ModProc, ONLY: PM_iProc, PM_nProc
  use PM_ModMain, ONLY: PM_FieldNorth_II, PM_FieldSouth_II, PM_FieldUnit

  implicit none
  !INPUT ARGUMENTS:
  integer,          intent(in) :: iSize, jSize            ! Size of buffer
  real,             intent(in) :: Buffer_II(iSize, jSize) ! Data buffer
  character(len=*), intent(in) :: NameVar                 ! Name of data

  !LOCAL VARIABLES:
  character (len=*), parameter :: NameSub = 'PM_put_from_xy'

  !DESCRIPTION:
  ! This subroutine shows an example for receiving data from another component.
  ! This subroutine should be in the srcInterface directory of the 
  ! specific PM component version. The subroutine is called by the
  ! couple\_xy\_pm method of the CON\_couple\_pm\_xy module in
  ! CON/Interface/src.
  ! 
  ! The data is assumed to be on a spherical grid of size iSize by jSize.
  ! Based on the NameVar argument the data corresponds either to the 
  ! northern or to the southern hemisphere. The northern hemisphere is
  ! always solved with processor 0, while the southern hemisphere with 
  ! processor PM\_nProc-1 where PM\_nProc is the number of processors (1 or 2)
  ! used by component PM.
  !EOP
  !---------------------------------------------------------------------------
  !BOC
  select case(NameVar)
  case('FieldNorth')
     if (PM_iProc /= 0) RETURN
     PM_FieldNorth_II = Buffer_II / PM_FieldUnit
  case('FieldSouth')
     if (PM_iProc /= PM_nProc-1) RETURN
     PM_FieldSouth_II = Buffer_II / PM_FieldUnit
  case default
     call CON_stop(NameSub//' invalid NameVar='//NameVar)
  end select
  !EOC
end subroutine PM_put_from_xy

!BOP =========================================================================
!MODULE: CON_couple_pm_xy - couple PM and XY components
!
!DESCRIPTION:
! This is an example for a module for the data exchange between two
! components. This module should be in the CON/Interface/src directory.
! The order of the component ID-s in the name of the module is alphabetical.
!
! This example provides methods for coupling the XY component to the PM 
! component (one way only). A single field is passed on the spherical 
! grid of the PM component. The communication assumes that the XY component 
! collects all data on its root processor, while the grid of the PM component
! consists of two blocks (two hemispheres) and it runs on 1 or 2 processors.
!
!INTERFACE:
module CON_couple_pm_xy

  !USES:
  use CON_coupler ! provides all the CON methods needed for the coupling

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:
  public :: couple_pm_xy_init   ! initialize coupling
  public :: couple_xy_pm        ! couple XY to PM

  !LOCAL VARIABLES:
  logical, save :: UseMe        ! logical for participating processors
  integer, save :: iSize, jSize ! size of the 2D spherical structured PM grid

  !EOP
contains

  !BOP =======================================================================
  !IROUTINE: couple_pm_xy_init - initialize XY-PM couplings
  !INTERFACE:
  subroutine couple_pm_xy_init

    !LOCAL VARIABLES:
    logical :: IsInitialized = .false.  ! logical to store initialization
    integer :: nCells_D(2)              ! temporary array for grid size

    !DESCRIPTION:
    ! This subroutine initializes the data exchange between the XY and PM
    ! components. The subroutine must be called by the couple\_all\_init
    ! method of the CON\_couple\_all module.
    ! This particular initialization stores the PM grid size iSiza and jSize
    ! and sets the logical UseMe to .true. for participating PE-s.
    !EOP
    !------------------------------------------------------------------------
    !BOC
    if(IsInitialized) RETURN
    IsInitialized = .true.

    ! Get the block size iSize and jSize from the PM grid descriptor 
    ! using the ncells_decomposition_d method of CON_coupler
    nCells_D = ncells_decomposition_d(PM_)
    iSize = nCells_D(1); jSize = nCells_D(2)

    ! Set UseMe to .true. for the participating PE-s
    ! using the is_proc() function of CON_coupler
    UseMe = is_proc(PM_) .or. is_proc(XY_)
    !EOC
  end subroutine couple_pm_xy_init

  !BOP =======================================================================
  !IROUTINE: couple_xy_pm - couple XY to PM
  !INTERFACE:
  subroutine couple_xy_pm(tSimulation)

    !INPUT ARGUMENTS:
    real, intent(in) :: tSimulation  ! simulation time

    !LOCAL VARIABLES:
    character (len=*), parameter :: NameSub='couple_xy_pm'

    character (len=*), parameter, dimension(2) :: &   ! names of variables
         NameVar_B = (/ 'FieldNorth', 'FieldSouth' /) ! for both blocks

    real, dimension(:,:), allocatable :: Buffer_II    ! buffer for 2D field

    integer :: nSize                      ! MPI message size
    integer :: iStatus_I(MPI_STATUS_SIZE) ! MPI status variable
    integer :: iError                     ! MPI error code

    integer :: iBlock                     ! block index
    integer :: iProcTo                    ! rank of the receiving processor

    !DESCRIPTION:
    ! This subroutine is called by the couple\_two\_comp method
    ! of the CON\_couple\_all module.
    ! This particular coupler send data Field from XY to PM.
    !EOP
    !-------------------------------------------------------------------------
    !BOC
    ! Exclude PEs which are not involved
    if(.not.UseMe) RETURN

    ! Do Northern and then Southern hemispheres
    do iBlock = 1, 2

       ! Get the rank of the PM processor for this block
       ! using the pe_decomposition() method of CON_coupler
       iProcTo = pe_decomposition(PM_, iBlock)

       ! Allocate buffers for the variables both in XY and PM
       allocate(Buffer_II(iSize, jSize), stat=iError)
       call check_allocate(iError, NameSub//": "//NameVar_B(iBlock))

       ! Calculate Field on XY.
       ! The result will be on the root processor of XY.
       if( is_proc(XY_) ) &
            call XY_get_for_pm(Buffer_II, iSize, jSize, NameVar_B(iBlock))

       ! Transfer variables from the root processor of XY to PM.
       ! Identify the root processor of PM with the i_proc0() method
       ! and the receiving processor with the i_proc() method of CON_coupler.
       if( iProcTo /= i_proc0(XY_))then
          nSize = iSize*jSize
          if(is_proc0(XY_)) &
               call MPI_send(Buffer_II, nSize, MPI_REAL, iProcTo,&
               1, i_comm(), iError)
          if(i_proc() == iProcTo) &
               call MPI_recv(Buffer_II, nSize, MPI_REAL, i_proc0(XY_),&
               1,i_comm(), iStatus_I, iError)
       end if

       ! Put variables into PM
       if( i_proc() == iProcTo )&
            call PM_put_from_xy(Buffer_II, iSize, jSize, NameVar_B(iBlock))

       ! Deallocate buffer to save memory
       deallocate(Buffer_II)
    end do

    !EOC
  end subroutine couple_xy_pm

end module CON_couple_pm_xy
