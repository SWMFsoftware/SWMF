!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module CON_main

  ! Main code of the Space Weather Modeling Framework performing a simulation

  use CON_world
  use CON_comp_param
  use CON_wrapper
  use CON_io, ONLY : SaveRestart, save_restart, IsRestartSaved

  use CON_time
  use CON_variables, ONLY: IsStandAlone, iErrorSwmf, lVerbose, DnTiming
  use CON_session
  use ModUtilities, ONLY: remove_file, touch_file
  use ModPlanetConst, ONLY: init_planet_const
  use CON_planet,     ONLY: set_planet_defaults
  use ModMpi, ONLY: MPI_WTIME, MPI_allreduce, MPI_INTEGER, MPI_MIN

  implicit none

  private ! except

  public :: initialize         ! initialize SWMF
  public :: finalize           ! finalize   SWMF

  ! revision history:
  !
  ! This module is a result of a continuous transformation from the main
  ! program of BATSRUS (developed at the University of Michigan)
  ! into what it is now. Probably not a single line is left untouched
  ! from the original code, but it originates from there.
  !
  ! The transformations were carried out by O.Volberg and G.Toth.
  !
  ! To allow an external program to call the SWMF, the main program
  ! is broken into subroutines which are called from swmf.f90 or
  ! swmf\_interface.f90. This module contains the subroutines
  ! initialize anf finalize.
  !

  character (len=*), parameter :: NameSub = 'CON_main'

  ! Local variable definitions.
  integer :: lComp, iComp

  logical :: DoTest, DoTestMe

contains
  !============================================================================
  subroutine initialize(iComm)
    integer, intent(in), optional :: iComm ! the MPI communicator for the SWMF

    ! The optional argument MPI communicator argument is present only
    ! when the SWMF is not running in stand alone mode.
    ! This subroutine executes the following major steps shown in
    ! pseudo F90 code
    ! \begin{itemize}
    ! \item Initialize the framework:                         \begin{verbatim}
    !       call world_init                                   \end{verbatim}
    ! \item Read the layout of the components from LAYOUT.in: \begin{verbatim}
    !       call world_setup                                  \end{verbatim}
    ! \item Obtain version information from and provide MPI information to the
    !       components:                                       \begin{verbatim}
    !       do iComp = 1, nComp
    !           call set_param_comp(iComp,"VERSION")
    !           call set_param_comp(iComp,"MPI")
    !       end do                                            \end{verbatim}
    ! \item Show the registered and unregistered components:  \begin{verbatim}
    !       call show_all_comp                                \end{verbatim}
    ! \end{itemize}
    ! The actual code is somewhat longer, since the main code also
    ! deals with timing, deleting the SWMF.STOP and SWMF.SUCCESS files
    ! at the beginning of the run, and initializing the planet information.
    ! There is also some verbose information printed.
    !--------------------------------------------------------------------------

    ! Initialize control component (MPI)
    call world_init(iComm)

    ! Set IsStandAlone variable:
    ! if the communicator is externally given,
    ! the SWMF is not running in stand alone mode.
    IsStandAlone = .not. present(iComm)

    ! Delete SWMF.SUCCESS, SWMF.DONE, SWMF.STOP and SWMF.KILL files if found
    if(is_proc0())then
       call remove_file('SWMF.SUCCESS') ! code successfully stopeed (no crash)
       call remove_file('SWMF.DONE')    ! code is done with the whole run
       call remove_file('SWMF.STOP')    ! stop the code when it listens
       call remove_file('SWMF.KILL')    ! kill the code ASAP
    end if

    if(is_proc0())then
       include 'show_git_info.h'
    end if

    ! Read component information from LAYOUT.in
    call world_setup

    ! Initialize CPU timing
    CpuTimeStart = MPI_WTIME()

    ! Read and store version name and number of registered components
    ! initialize the MPI parameters for the registered components
    do lComp = 1, n_comp()
       iComp=i_comp(lComp)
       call set_param_comp(iComp,"VERSION")
       if(.not.use_comp(iComp))then
          if(is_proc0())then
             write(*,'(a)')'CON_main SWMF_ERROR registered component '// &
                  NameComp_I(iComp)//' is OFF!'
             write(*,'(a)')&
                  'Compile in a working component or remove it from '//&
                  NameMapFile
          end if
          ! stop in a clean fashion (without abort)
          call world_clean
          iErrorSwmf = 3
          RETURN
       end if

       call set_param_comp(iComp,'MPI')     ! Initialize MPI parameters
       call set_param_comp(iComp,'STDOUT')  ! Set prefix string for STDOUT
    end do

    ! Show framework and component information
    if(is_proc0()) call show_all_comp

    ! Check for illegal overlap of components with shared source code
    call check_overlap_comp
    if(iErrorSwmf /= 0) RETURN

    ! Initialize CON_time
    call init_time

    ! Initialize the planetary constant library and
    ! set the default planet (Earth)
    call init_planet_const
    call set_planet_defaults

  end subroutine initialize
  !============================================================================

  subroutine finalize
    ! This subroutine executes the following major steps shown in
    ! pseudo F90 code
    ! begin{itemize}
    ! \item Save final restart files if required:             \begin{verbatim}
    !    if(SaveRestart % DoThis) call save_restart           \end{verbatim}
    ! \item Finalize components:                              \begin{verbatim}
    !    do iComp=1,nComp
    !        call finalize_comp(iComp,tSimulation)
    !    end do                                               \end{verbatim}
    ! \item Finish the execution:                             \begin{verbatim}
    !    call world_clean                                     \end{verbatim}
    ! end{itemize}
    ! The actual code is longer: verbose information is printed, timing
    ! report is shown and the SWMF.SUCCESS and SWMF.DONE files are written.

    use ModConst, ONLY: &
         cSecondPerYear, cSecondPerDay, cSecondPerHour, cSecondPerMinute

    !--------------------------------------------------------------------------
    if(is_proc0())then
       if(lVerbose>=0)then
          write(*,*)
          write(*,'(a)')'    Finished Numerical Simulation'
          write(*,'(a)')'    -----------------------------'
          if (DoTimeAccurate)then
             write(*, '(a,es13.5,a)') &
               '   Simulated time = ', tSimulation, ' s '
             if(tSimulation > 1e6*cSecondPerYear) then
                write(*, '(a,es13.5,a)') &
                     tSimulation/cSecondPerYear
             elseif(tSimulation > cSecondPerYear) then
                write(*, '(a,f13.6,a)') '    Simulated time = ', &
                     tSimulation/cSecondPerYear, ' years'
             elseif(tSimulation > cSecondPerDay) then
                write(*, '(a,f13.6,a)') '    Simulated time = ', &
                     tSimulation/cSecondPerDay, ' days'
             elseif(tSimulation > cSecondPerHour) then
                write(*, '(a,f13.6,a)') '    Simulated time = ', &
                     tSimulation/cSecondPerHour, ' hours'
             elseif(tSimulation > cSecondPerMinute) then
                write(*, '(a,f13.6,a)') '    Simulated time = ', &
                     tSimulation/cSecondPerMinute, ' mins'
             end if
          end if
       end if
    end if
    if (DnTiming > -2) call timing_report

    if(SaveRestart % DoThis .and. .not. IsRestartSaved)then
       if(UseEndTime .and. .not. IsForcedStop)then
          ! Save TimeEnd to RESTART.out for future simulation
          if(is_proc(CON_)) call save_restart(DoFinalize=.true.)
       else
          if(is_proc(CON_)) call save_restart
       end if
    end if
    do lComp = 1,n_comp()
       iComp = i_comp(lComp)
       call finalize_comp(iComp,tSimulation)
    end do

    if(is_proc0().and.lVerbose>0)then
       write(*,*)
       write(*,'(a)')'    Finished Finalizing SWMF'
       write(*,'(a)')'    ------------------------'
    end if

    call timing_stop('SWMF')

    if(DnTiming > -3)call timing_report_total

    ! Signal normal completion by writing an empty SWMF.SUCCESS file
    if(is_proc0()) call touch_file('SWMF.SUCCESS')
    if(is_proc0() .and. .not. IsForcedStop) call touch_file('SWMF.DONE')

    ! Stop running
    call world_clean

  end subroutine finalize
  !============================================================================
  subroutine show_all_comp

    use CON_comp_param, ONLY: lNameVersion

    ! Show the version information and layout for all registered components.
    ! Show the version information for all working but unregistered components.
    !
    ! revision history:
    ! 08/2003 G.Toth <gtoth@umich.edu> - initial version and impovements

    integer, parameter :: lWidth = 68

    logical                      :: IsOn, IsFound
    character (LEN=lNameVersion) :: NameVersion

    integer :: lComp, iComp, iProc0Comp, nProcComp, iStrideComp, nThread

    !--------------------------------------------------------------------------
    NameVersion  = 'SWMF by Univ. of Michigan'

    write(*,'(a)')'#'//repeat('=',lWidth)//'#'
    write(*,'(2a)') &
         '# ID  Version                             ',&
         'nproc proc0 stride nthread #'
    write(*,'(a)')'#'//repeat('-',lWidth)//'#'
    write(*,'(a,a,4i6,a6)') &
         '# CON ', NameVersion, n_proc(), 0, 1, 1,'#'
    write(*,'(a)')'#'//repeat('-',lWidth)//'#'

    ! Show registered components
    do lComp = 1, n_comp()
       iComp = i_comp(lComp)

       call get_comp_info(iComp, NameVersion=NameVersion, &
            iProcZero=iProc0Comp, nProc=nProcComp, iProcStride=iStrideComp,&
            nThread=nThread)

       write(*,'(a,4i6,a6)')&
            '# '//NameComp_I(iComp)//'  '//NameVersion, &
            nProcComp, iProc0Comp, iStrideComp, nThread,'#'
    end do

    ! Show unregistered but compiled (ON) components
    IsFound = .false.
    do iComp = 1, MaxComp
       if(use_comp(iComp)) CYCLE ! registered component
       call get_version_comp(iComp, IsOn, NameVersion)
       if(.not.IsOn) CYCLE       ! empty component
       if(.not.IsFound)then
          write(*,'(a)')'#'//repeat('-',lWidth)//'#'
          IsFound = .true.
       end if
       write(*,'(a)') '# '//NameComp_I(iComp)//'  '//NameVersion// &
            '        not registered       #'
    end do
    write(*,'(a)')'#'//repeat('=', lWidth)//'#'

  end subroutine show_all_comp
  !============================================================================

  subroutine check_overlap_comp
    ! One may use the same source code for two non-overlapping components.
    ! For example GM and IH may use the same MHD source code. This is indicated
    ! by their NameVersion being the same. When the version names are the same
    ! no overlap is allowed. This suboutine stops with an appropriate error
    ! message if such an overlap occurs.
    !
    ! revision history:
    ! 09/10/03 - G.Toth <gtoth@umich.edu> - initial version and prolog

    integer :: lComp, iComp, lComp2, iComp2, iProc, iProcMin, iError
    character (LEN=lNameVersion) :: NameVersion, NameVersion2

    !--------------------------------------------------------------------------
    do lComp=2,n_comp()
       iComp = i_comp(lComp)
       call get_comp_info(iComp, NameVersion=NameVersion)

       do lComp2 = 1, lComp - 1
          iComp2 = i_comp(lComp2)
          call get_comp_info(iComp2, NameVersion=NameVersion2)

          ! Store processor rank for illegal overlap
          if(NameVersion==NameVersion2 .and. &
               is_proc(iComp) .and. is_proc(iComp2)) then
             iProc = i_proc()
          else
             iProc = n_proc()+1
          end if

          ! Let all processors now if there was
          call MPI_allreduce(iProc,iProcMin,1,MPI_INTEGER,MPI_MIN,&
               i_comm(),iError)

          if( iProcMin <= n_proc() )then
             if(is_proc0())then
                write(*,'(a)')'CON_main SWMF_ERROR: Components '// &
                     NameComp_I(iComp)//' and '//NameComp_I(iComp2)
                write(*,'(a)')'CON_main SWMF_ERROR:'//&
                     ' have the same version '//trim(NameVersion)
                write(*,'(a,i3,a)')'CON_main SWMF_ERROR:'//&
                     ' and both components use processor ',iProcMin,' !'
                write(*,'(a)')'CON_main SWMF_ERROR:'// &
                     ' Remove overlap or use different/renamed versions!'
             end if
             call world_clean
             iErrorSwmf = 4
          end if
       end do
    end do

  end subroutine check_overlap_comp
  !============================================================================

end module CON_main
!==============================================================================

