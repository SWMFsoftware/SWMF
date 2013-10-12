!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
! File name: dgcpm_output_010.f90
!
! Contains: output routines for DGCPM
!	WRESULT
!
! Last Modified: December 2006, Mike Liemohn
!
! **********************************************************************
!				WRESULT
!       Routine prints all the results at time T after injection
!***********************************************************************
	SUBROUTINE WRESULT(IFIR)

	use ModIoDGCPM
        use ModTimeDGCPM
        use ModMainDGCPM

	implicit none

	integer i,j,IFIR
	INTEGER NTC
	CHARACTER*5 ST1,ST3
	CHARACTER*5 SUF
	CHARACTER*2 ST2,SUF2
	CHARACTER*1 SUF1
 	CHARACTER*80 filename
	SAVE NTC


!......... Dynamic Data Output
!.......Define parts of the output file names
	IF (IFIR.EQ.1) THEN
	  NTC=INT(NINT(TIME/TINT))
	ELSE
	  NTC=NTC+1
	END IF

        write(suf,'(i5.5)') ntc

	ST1=NAME
        
        call getdensity(vthetacells,nthetacells,vphicells,nphicells,   &
              dendgcpm)
 
!........ Static File Selector Writing
      filename=cOutputDir//ST1//'_dgcpm_static.dat'

      if(WriteStatic)then
        SELECT CASE (MagneticType)
            CASE ('DIPOLE')
                IF (OutputType.ne.'OLD') then
                    open(unit=12, file=filename, form='formatted')
                    write(12,*) 'NTHETA NPHI THETA PHI X Y OC VOL'
                    write(12,*) nthetacells, nphicells
                    write(12,*) 90.0-vthetacells
                    write(12,*) vphicells
                    write(12,*) mgridx
                    write(12,*) mgridy
                    write(12,*) mgridoc
                    write(12,*) mgridvol
                    close(unit = 12)
                ENDIF
             CASE ('T96')
                write(*,*) 'NOT IMPLIMENTED - NO OUTPUT'
             CASE ('SHUE')
                write(*,*) 'NOT IMPLIMENTED - NO OUTPUT'
          END SELECT
          WriteStatic = .false.
        endif

        filename=cOutputDir//ST1//'_dgcpm_dynamic_'//SUF//'.dat'

!........... Type Selection


      if(WriteDynamic)then
        open(unit=10, file=filename, form = 'formatted')
        
        SELECT CASE (OutputType)
            CASE ('SHORT')
                write(10,*) 'TIME NTHETA NPHI N POT'
                write(10,*) CurrentTime
                write(10,*) nthetacells, nphicells
                write(10,*) mgridden
                write(10,*) mgridpot
            CASE ('VELOCITY')
                write(10,*) 'TIME NTHETA NPHI X Y DEN POT VR VP'
                write(10,*) CurrentTime
                write(10,*) nthetacells, nphicells
                write(10,*) mgridx
                write(10,*) mgridy
                write(10,*) mgridden
                write(10,*) mgridpot
                write(10,*) mgridvr
                write(10,*) mgridvp
            CASE ('POTENTIAL')
                write(10,*) 'TIME NTHETA NPHI THETA PHI X Y POT CORO'
                write(10,*) CurrentTime
                write(10,*) nthetacells, nphicells
                write(10,*) 90.0 - vthetacells
                write(10,*) vphicells
                write(10,*) mgridx
                write(10,*) mgridy
                write(10,*) mgridpot
                write(10,*) mgridcoro
            CASE ('FLOWS')
                write(10,*) 'TIME NTHETA NPHI THETA PHI X Y OC VOL'
                write(10,*) CurrentTime
                write(10,*) nthetacells, nphicells
                write(10,*) 90.0 - vthetacells
                write(10,*) vphicells
                write(10,*) mgridx
                write(10,*) mgridy
                write(10,*) mgridfluxr
                write(10,*) mgridfluxa
            CASE ('OLD')
                write(10,*) nthetacells, nphicells
                write(10,*) 90.0-vthetacells
                write(10,*) vphicells
                write(10,*) mgridden
                write(10,*) mgridx
                write(10,*) mgridy
                write(10,*) mgridoc
                write(10,*) mgridpot
                write(10,*) mgridcoro
                write(10,*) mgridvr
                write(10,*) mgridvp
                write(10,*) mgridsource
                write(10,*) mgridfluxr
                write(10,*) mgridfluxa
                write(10,*) mgridn
                write(10,*) mgridvol
                write(10,*) CurrentTime
        END SELECT
    
        close(unit = 10)
        endif
  
      if(WriteRestart)then
        filename=cOutputDir//ST1//'_restart.dat'
        open(unit=15, file=filename, form = 'formatted')

            write(15,*) nthetacells, nphicells
            write(10,*) vphicells
            write(10,*) mgridden
            write(10,*) mgridx
            write(10,*) mgridy
            write(10,*) mgridoc
            write(10,*) mgridpot
            write(10,*) mgridvr
            write(10,*) mgridvp
            write(10,*) mgridn
            write(10,*) mgridvol

        close(unit = 15)
      endif
  
      RETURN
      END
