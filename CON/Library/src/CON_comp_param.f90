!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module CON_comp_param

  ! parameters used by the SWMF Registry and query functions for them
  !
  ! This module contains the constants specific to the SWMF components,
  ! such as names, ID-s, number of components, etc.
  ! It also contains constants such as named indexes and file names
  ! used by CON_world and CON_comp_info.
  !
  ! The following Physics Components can occur in SWMF
  ! (listed in alphabetical order):
  !
  !    CZ Convection Zone
  !    EE Eruptive Event
  !    GM Global Magnetosphere
  !    IE Ionospheric Electrodynamics
  !    IH Inner Heliosphere
  !    IM Inner Magnetosphere
  !    OH Outer Heliosphere
  !    PC Particle-in-Cell
  !    PS Plasmasphere
  !    PT Particle Tracker
  !    PW Polar Wind
  !    RB Radiation Belts
  !    SC Solar Corona
  !    SP Solar Energetic Particles
  !    UA Upper Atmosphere

  implicit none

  !PUBLIC DATA MEMBERS:

  integer, parameter :: MaxComp   = 15 ! maximum number of components
  integer, parameter :: lNameComp =  2 ! length of component names

  ! Convert component index to component name
  character(len=lNameComp), parameter :: NameComp_I(MaxComp) = [ &
       "EE", "GM", "IE", "IH", "IM", "OH", "PC", "PS", "PT", "PW", &
       "RB", "SC", "SP", "UA", "CZ"]

  ! Named indexes for the components
  integer, parameter :: &
       EE_=1, GM_=2, IE_=3, IH_=4 , IM_=5, OH_=6, PC_=7, PS_=8, PT_=9, PW_=10,&
       RB_=11, SC_=12 ,SP_=13, UA_=14, CZ_=15

  ! Length of the version name of the component
  integer, parameter :: lNameVersion=34

  ! Named indexes of the MPI parameters
  integer, parameter :: ProcZero_  = 1 ! the spokesman for the group
  integer, parameter :: ProcLast_  = 2 ! the upper limit of the group rank
  integer, parameter :: ProcStride_= 3 ! the group stride
  integer, parameter :: Proc_      = 4 ! the processor index in the group
  integer, parameter :: nProc_     = 5 ! the number of processors of the group
  integer, parameter :: Comm_      = 6 ! the component global communicator
  integer, parameter :: Group_     = 7 ! the component MPI group

  integer, parameter :: nMpiParam = 7  ! number of MPI parameters

  ! Name of the file containing the processor map
  character(len=9)   :: NameMapFile = "PARAM.in"

  public :: is_valid_comp_name ! return true if component name is valid
  public :: i_comp_name        ! return index for component name

  ! revision history:
  !
  !  June    2003 - O. Volberg <volov@umich.edu> - initial version
  !  July 12 2003 - G. Toth    <gtoth@umich.edu> - rewrite
  !

  character(len=*),parameter,private :: NameMod = 'CON_comp_param'

contains
  !============================================================================
  logical function is_valid_comp_name(NameComp)

    character(len=*), intent(in) :: NameComp
    integer :: iComp
    character(len=*), parameter:: NameSub = 'is_valid_comp_name'
    !--------------------------------------------------------------------------
    iComp=i_comp_name(NameComp)
    is_valid_comp_name = iComp/=0

  end function is_valid_comp_name
  !============================================================================
  integer function i_comp_name(Name)

    character(len=*), intent(in) :: Name
    integer :: iComp

    ! This could be made more efficient by using select case(Name) construct

    character(len=*), parameter:: NameSub = 'i_comp_name'
    !--------------------------------------------------------------------------
    do iComp=1, MaxComp
       if(Name == NameComp_I(iComp)) then
          i_comp_name = iComp
          RETURN
       end if
    end do

    i_comp_name=0

  end function i_comp_name
  !============================================================================
  subroutine check_i_comp(iComp,NameCaller)

    use ModUtilities, ONLY: CON_stop

    integer, intent(in) :: iComp
    character (len=*), intent(in) :: NameCaller
    !--------------------------------------------------------------------------
    if( iComp<1 .or. iComp>MaxComp) then
       write(*,'(a,i3,a,i3)')NameCaller//' SWMF_ERROR iComp ',iComp,&
            ' is out of range 1 ..',MaxComp
       call CON_stop('Error in the caller method')
    end if

  end subroutine check_i_comp
  !============================================================================
end module CON_comp_param
!==============================================================================
