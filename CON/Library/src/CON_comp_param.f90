!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!               Space Weather Modeling Framework (SWMF)                !
!    Center for Space Environment Modeling, The University of Michigan !
!-----------------------------------------------------------------------
!BOP -------------------------------------------------------------------------
!
! !MODULE: CON_comp_param - parameters used by the SWMF Registry 
!          and query functions for them
! 
! !DESCRIPTION: 
! This module contains the constants specific to the SWMF components,
! such as names, ID-s, number of components, etc.
! It also contains constants such as named indexes and file names
! used by CON\_world and CON\_comp\_info.
!
! The following Physics Components can occur in SWMF
! (listed in alphabetical order): 
!
!\begin{itemize}
!\item[EE] Eruptive Event
!\item[GM] Global Magnetosphere
!\item[IE] Ionospheric Electrodynamics
!\item[IH] Inner Heliosphere
!\item[IM] Inner Magnetosphere
!\item[LA] Lower Atmosphere
!\item[OH] Outer Heliosphere
!\item[PS] Plasmasphere
!\item[PW] Polar Wind
!\item[RB] Radiation Belts
!\item[RT] Radiative Transfer
!\item[SC] Solar Corona
!\item[LC] Lower Solar Corona
!\item[SP] Solar Energetic Particles
!\item[UA] Upper Atmosphere
!\end{itemize}

!INTERFACE:	
!
module CON_comp_param
  !
  implicit none

  !PUBLIC DATA MEMBERS:

  integer, parameter :: MaxComp   = 15 ! maximum number of components
  integer, parameter :: lNameComp =  2 ! length of component names

  ! Convert component index to component name
  character(len=lNameComp), parameter :: NameComp_I(MaxComp) = (/ & 
       "EE", "GM", "IE", "IH", "IM", "LA", "OH", "PS", "PW", &
       "RB", "RT", "SP", "SC", "LC", "UA"/)

  ! Named indexes for the components
  integer, parameter :: &
       EE_=1, GM_=2, IE_=3, IH_=4 , IM_=5 ,LA_=6, OH_=7, PS_=8, PW_=9, &
       RB_=10 ,RT_=11, SP_=12 ,SC_=13 ,LC_=14, UA_=15

  ! Length of the version name of the component
  integer, parameter :: lNameVersion=40 

  ! Named indexes of the MPI parameters
  integer, parameter :: ProcZero_  = 1 ! the spokesman for the group
  integer, parameter :: ProcLast_  = 2 ! the upper limit of the group rank
  integer, parameter :: ProcStride_= 3 ! the group stride
  integer, parameter :: Proc_      = 4 ! the processor index in the group
  integer, parameter :: nProc_     = 5 ! the number of processors of the group 
  integer, parameter :: Comm_      = 6 ! the component global communicator
  integer, parameter :: Group_     = 7 ! the component MPI group

  integer, parameter :: nMpiParam = 7  ! number of MPI parameters

  ! Name of the processor map file
  character (len=*), parameter  :: NameMapFile = "LAYOUT.in"

  !
  !PUBLIC MEMBER FUNCTIONS:
  !
  public :: is_valid_comp_name ! return true if component name is valid
  public :: i_comp_name        ! return index for component name

  !REVISION HISTORY: 
  !
  !  June    2003 - O. Volberg <volov@umich.edu> - initial version
  !  July 12 2003 - G. Toth    <gtoth@umich.edu> - rewrite
  !  Sep  16 2009 - R. Oran    <oran@umich.edu>  - add LC component
  !
  !EOP ------------------------------------------------------------------------

  character(len=*),parameter,private :: NameMod = 'CON_comp_param'

contains

  !==================================================================
  logical function is_valid_comp_name(NameComp)
    character(len=*), parameter :: NameSub=NameMod//'::is_valid_comp_name'
    character(len=*), intent(in) :: NameComp
!    integer, intent (out), optional :: iCompOut
    integer :: iComp
    !---------------------------------------------------------------
    iComp=i_comp_name(NameComp)
    is_valid_comp_name = iComp/=0
!    if(present(iCompOut)) iCompOut=iComp

  end function is_valid_comp_name
  !==================================================================
  integer function i_comp_name(Name)
    character(len=*), parameter  :: NameSub=NameMod//'::i_comp_name'
    character(len=*), intent(in) :: Name
    integer :: iComp

    ! This could be made more efficient by using select case(Name) construct

    do iComp=1, MaxComp
       if(Name == NameComp_I(iComp)) then
          i_comp_name = iComp
          RETURN
       end if
    end do

    i_comp_name=0

  end function i_comp_name
  !===========================================================================
  subroutine check_i_comp(iComp,NameCaller)
    integer, intent(in) :: iComp
    character (len=*), intent(in) :: NameCaller
    !-------------------------------------------------------------------------
    if( iComp<1 .or. iComp>MaxComp) then
       write(*,'(a,i3,a,i3)')NameCaller//' SWMF_ERROR iComp ',iComp,&
            ' is out of range 1 ..',MaxComp
       call CON_stop('Error in the caller method')
    end if
    
  end subroutine check_i_comp
  !===========================================================================

end module CON_comp_param
