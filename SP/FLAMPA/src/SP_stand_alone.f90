!============================================================================!
!============================================================================!
! This module is a template that can be used to convert the MAIN program     !
! of a model, implemented in F90/95, so that it complies with the SWMF       !
! requirements and retain its ability to work as a STAND-ALONE module.       !
!============================================================================!
program CON_stand_alone
  use SP_ModMain
  implicit none
  !--------------------------------------------------------------------------!
  interface
     subroutine SP_diffusive_shock(&                  
          TypeAction,              &                 
          TimeToFinish,            &               
          nToFinish)                    
       character(LEN=*),intent(in):: TypeAction
       real,intent(in),optional   :: TimeToFinish
       integer,intent(in),optional:: nToFinish
     end subroutine SP_diffusive_shock
  end interface
  include 'stdout.h'
  !--------------------------------------------------------------------------!
  integer:: iIter,iX
  character(LEN=15):: NameFile
  !--------------------------------------------------------------------------!
  call set_stdout('SP: ')      !Set the string to be printed on screen.   !
  DoTest=.false.
  if(DoTest)write(iStdOut,'(a)')prefix//'Do the test only'
  if(DoTest)then
     !------------------------------- INIT ----------------------------------!
     call SP_diffusive_shock("INIT")
     call SP_allocate
     !-------------------------------- RUN ----------------------------------!
     ! The following is a test problem (DoTest=.true.):                      !
     ! 1D shock wave with unit speed propagates through half of the spatial  !
     ! interval (0,nX). The shock's compression ratio is set to 4.0.         !
     !-----------------------------------------------------------------------!
     do iIter=1,nX/2
        write(iStdOut,*)prefix,"iIter = ",iIter
        Rho_I(iIter) = cTwo+cHalf
        if (iIter>1) then
           Rho_I(iIter-1) = cFour
           do iX=1,iIter-1
              X_DI(1,iX) = cQuarter*real(iX+3*iIter)
           end do
        end if
        call SP_diffusive_shock("RUN",real(iIter))
     end do
  else
     !------------------------------- INIT ----------------------------------!
     EInjection=cOne
     !     BOverDeltaB2=cThree
     call SP_diffusive_shock("INIT")
     iFileToRead=114
     write(NameFile,'(a,i4.4,a)')'ihdata_',iFileToRead,'.dat'
     call read_ihdata_for_sp(NameFile,1,5)
     RhoOld_I=RhoSmooth_I
     RhoOld_I(iShock+1-nint(cOne/DsResolution):iShock)=&
          maxval(RhoOld_I(iShock+1-nint(cOne/DsResolution):iShock))
     RhoOld_I(iShock+1:iShock+nint(cOne/DsResolution))=&
          minval(RhoOld_I(iShock+1:iShock+nint(cOne/DsResolution)))
     SP_Time=TimeToRead
     DiffCoeffMin=1.0e+05*Rsun*DsResolution !m^2/s
     !-------------------------------- RUN ----------------------------------!
     do iFileToRead=115,493
        write(iStdOut,*)prefix,"iIter = ",iFileToRead-114
        write(NameFile,'(a,i4.4,a)')'ihdata_',iFileToRead,'.dat'
        call read_ihdata_for_sp(NameFile,1,5)
        call SP_diffusive_shock("RUN",TimeToRead)
     end do
  endif
  !--------------------------------- FINALIZE -------------------------------!
  call SP_diffusive_shock("FINALIZE")
  !------------------------------------ DONE --------------------------------!
end program CON_stand_alone
!============================================================================!
