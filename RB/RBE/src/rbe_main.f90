!******************************************************************************
!
!                               rbe_main.f90
! 
!  The code is created to be integrated into the Space Weather
!  Modeling Framework.
!
!  Created on 21 March 2006 by Mei-Ching Fok, Code 612.2, NASA GSFC.
!
!  Radiation Belt-Ring Current Forecasting Model, Version 02.
!
!  Contact: Mei-Ching Fok at mei-ching.h.fok@nasa.gov, 301-286-1083.
!
!  This code is designed to model radiation belt electrons and ions with 
!  energy ranges 10 keV - 4 MeV for e-, 10 keV - 1 MeV for H+.
!
!  A program calculates temporal evolutions of radiation belt particle fluxes,
!  considering drift, losscone loss, radial, pitch-angle and energy diffusions.
!
!  Model limit:
!       r:    r < rb        rb is set at 12 RE
!      LT:    0 - 24 
!       M:    corresponding to 10 keV - 4 MeV for e-, 10 keV - 1 MeV for H+.
!       K:    corresponding to sine of equatorial pitch angle from 0 - 1.
!
!  Magnetic field model: Tsyganenko 96, 04 model (t96_01.f or t04_s.f) or MHD.
!
!  Electric field model: Weimer 2k model (w2k.f, w2k.dat) or MHD.
!
!  Plasmasphere model: Dan Ober's model (pbo_2.f)
!
!  Input files: rbe_swmf.dat
!               storm.SWIMF
!               storm.symH 
!               w2k.dat
!               rbe_*.fin for initial run
!               outname_*_c.f2 for continuous run
!
!  Output files: outname_*.fls (* can be e or h for electrons or protons)
!******************************************************************************

program rbe

  use rbe_cread2, ONLY: nstept, nstep, IsStandAlone, iConvect
  use rbe_time,   ONLY: istep,t
  use ModMpi,ONLY: MPI_COMM_WORLD
  use ModReadParam
  use ModPrerunField, ONLY: UsePrerun, read_prerun, read_prerun_IE
  implicit none
  !---------------------------------------------------------------------------
  
  IsStandAlone=.true.

  ! Initial setup for the rbe model
  call read_file('PARAM.in',MPI_COMM_WORLD)
  call read_init('  ',iSessionIn=1,iLineIn=0)
  call RB_set_parameters('READ')
  if (usePrerun) call read_prerun(t)
  if (usePrerun .and. iConvect==2) call read_prerun_IE(t)
  call readInputData
  call timing_active(.true.)
  call timing_step(0)
  call timing_start('RBE')
  call timing_start('rbe_init')
  call rbe_init
  call timing_stop('rbe_init')
  call timing_report_total
  call timing_reset('#all',3)

  ! start the calculation, the time loop
  do istep = 1, nstept+nstep 
     call timing_step(istep)
     call timing_start('rbe_run')
     if (usePrerun) call read_prerun(t)
     if (usePrerun .and. iConvect==2) call read_prerun_IE(t)
     call rbe_run
     call timing_stop('rbe_run')
  end do

  call timing_stop('RBE')
  call timing_report

end program rbe
!============================================================================
subroutine CON_stop(String)

  character(len=*), intent(in) :: String
  write(*,*)'ERROR in RBE:',String
  stop

end subroutine CON_stop
!============================================================================
subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe
  DoTest   = .false.
  DoTestMe = .false.
end subroutine CON_set_do_test
