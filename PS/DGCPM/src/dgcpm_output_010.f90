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
!	IRES(7)  'pla'	Thermal plasma densities
!***********************************************************************
	SUBROUTINE WRESULT(IFIR)

	use ModIoDGCPM
	use ModHeidiDGCPM

	implicit none

	integer i,j,IFIR
	INTEGER NTC
	CHARACTER*5 ST1,ST3
	CHARACTER*3 SUF
	CHARACTER*2 ST2,SUF2
	CHARACTER*1 SUF1
 	CHARACTER*80 filename
	SAVE NTC

!.......Define parts of the output file names
	IF (IFIR.EQ.1) THEN
	  NTC=INT(NINT(TIME/TINT))
	ELSE
	  NTC=NTC+1
	END IF

        write(suf,'(i3.3)') ntc

	ST1=NAME

!.......Write the plasmaspheric thermal densities (IRES(7), 'pla' & 'dgcpm')
	print *, 'Printing plasmasphere'
!	  First create Dan's output file for his plotting software
	  filename=cOutputDir//ST1//'_dgcpm_'//SUF//'.dat'
	call getdensity(vthetacells,nthetacells,vphicells,nphicells,   &
              dendgcpm)
!	call saveit(vthetacells,nthetacells,vphicells,nphicells,   &
!              dendgcpm,gridx,gridy,gridoc,filename)

        write(*,*) "Saving plasmasphere File : ",filename

	call saveplasmasphere(filename)
!	  Next create an output like we have made before
!!  Deleted printing DGCPM results on the RAM/HEIDI grid

      RETURN
      END
!
! End of subroutine WRESULT
!
