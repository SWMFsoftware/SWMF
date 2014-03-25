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
        use ModIoUnit, ONLY: UNITTMP_

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
                    open(unit=UNITTMP_, file=filename, form='formatted')
                    write(UNITTMP_,*) 'NTHETA NPHI THETA PHI X Y OC VOL'
                    write(UNITTMP_,*) nthetacells, nphicells
                    write(UNITTMP_,*) 90.0-vthetacells
                    write(UNITTMP_,*) vphicells
                    write(UNITTMP_,*) mgridx
                    write(UNITTMP_,*) mgridy
                    write(UNITTMP_,*) mgridoc
                    write(UNITTMP_,*) mgridvol
                    close(unit = UNITTMP_)
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
        open(unit=UNITTMP_, file=filename, form = 'formatted')
        
        SELECT CASE (OutputType)
            CASE ('SHORT')
                write(UNITTMP_,*) 'TIME NTHETA NPHI N POT'
                write(UNITTMP_,*) CurrentTime
                write(UNITTMP_,*) nthetacells, nphicells
                write(UNITTMP_,*) mgridden
                write(UNITTMP_,*) mgridpot
            CASE ('VELOCITY')
                write(UNITTMP_,*) 'TIME NTHETA NPHI X Y DEN POT VR VP'
                write(UNITTMP_,*) CurrentTime
                write(UNITTMP_,*) nthetacells, nphicells
                write(UNITTMP_,*) mgridx
                write(UNITTMP_,*) mgridy
                write(UNITTMP_,*) mgridden
                write(UNITTMP_,*) mgridpot
                write(UNITTMP_,*) mgridvr
                write(UNITTMP_,*) mgridvp
            CASE ('POTENTIAL')
                write(UNITTMP_,*) 'TIME NTHETA NPHI THETA PHI X Y POT CORO'
                write(UNITTMP_,*) CurrentTime
                write(UNITTMP_,*) nthetacells, nphicells
                write(UNITTMP_,*) 90.0 - vthetacells
                write(UNITTMP_,*) vphicells
                write(UNITTMP_,*) mgridx
                write(UNITTMP_,*) mgridy
                write(UNITTMP_,*) mgridpot
                write(UNITTMP_,*) mgridcoro
            CASE ('FLOWS')
                write(UNITTMP_,*) 'TIME NTHETA NPHI THETA PHI X Y OC VOL'
                write(UNITTMP_,*) CurrentTime
                write(UNITTMP_,*) nthetacells, nphicells
                write(UNITTMP_,*) 90.0 - vthetacells
                write(UNITTMP_,*) vphicells
                write(UNITTMP_,*) mgridx
                write(UNITTMP_,*) mgridy
                write(UNITTMP_,*) mgridfluxr
                write(UNITTMP_,*) mgridfluxa
            CASE ('OLD')
                write(UNITTMP_,*) nthetacells, nphicells
                write(UNITTMP_,*) 90.0-vthetacells
                write(UNITTMP_,*) vphicells
                write(UNITTMP_,*) mgridden
                write(UNITTMP_,*) mgridx
                write(UNITTMP_,*) mgridy
                write(UNITTMP_,*) mgridoc
                write(UNITTMP_,*) mgridpot
                write(UNITTMP_,*) mgridcoro
                write(UNITTMP_,*) mgridvr
                write(UNITTMP_,*) mgridvp
                write(UNITTMP_,*) mgridsource
                write(UNITTMP_,*) mgridfluxr
                write(UNITTMP_,*) mgridfluxa
                write(UNITTMP_,*) mgridn
                write(UNITTMP_,*) mgridvol
                write(UNITTMP_,*) CurrentTime
        END SELECT
    
        close(unit = UNITTMP_)
        endif
  
      if(WriteRestart)then
        filename=cOutputDir//ST1//'_restart.dat'
        open(unit=UNITTMP_, file=filename, form = 'formatted')

            write(UNITTMP_,*) nthetacells, nphicells
            write(UNITTMP_,*) vphicells
            write(UNITTMP_,*) mgridden
            write(UNITTMP_,*) mgridx
            write(UNITTMP_,*) mgridy
            write(UNITTMP_,*) mgridoc
            write(UNITTMP_,*) mgridpot
            write(UNITTMP_,*) mgridvr
            write(UNITTMP_,*) mgridvp
            write(UNITTMP_,*) mgridn
            write(UNITTMP_,*) mgridvol

        close(unit = UNITTMP_)
      endif
  
      RETURN
      END

!=============================================================================
subroutine write_lslice
  
  use ModIoDGCPM,     ONLY: cOutputDir, iUnitSlice
  use ModMainDGCPM,   ONLY: vrcells, mgridden, vmltcells
  use ModSizeDGCPM,   ONLY: nthetacells, nphicells
  use ModTimeDGCPM,   ONLY: CurrentTime
  use ModIoUnit,      ONLY: io_unit_new
  use ModTimeConvert, ONLY: TimeType, time_real_to_int

  implicit none

  real          :: radiusSlice = 6.6 ! Radius of slice.
  real          :: density(nphicells) = 0
  integer       :: i, j
  type(TimeType):: TimeNow
  integer, save :: iL=-1
  logical, save :: IsFirstWrite=.true.
  character(len=12)  :: StrFmt1
  character(len=24)  :: StrFmt2
  character(len=100) :: NameFile

  ! Testing params:
  logical :: DoTest, DoTestMe
  character(len=*), parameter :: NameSub='write_lslice'
  !------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  if(IsFirstWrite)then
     ! Create file name and open.
     i = floor(radiusSlice)
     j = floor(100.0*radiusSlice-100.0*i)
     write(NameFile, "(a,'/slice_',i1.1,'_',i2.2,'.txt')")trim(cOutputDir), i, j
     iUnitSlice = io_unit_new()
     open(iUnitSlice, FILE=NameFile, STATUS='REPLACE')

     ! Write header.
     write(StrFmt1, '(a,i3,a)') '(a, ', nphicells,  'f7.3)'
     write(iUnitSlice, '(a,f7.4)')'Radial Slice at L=', radiusSlice
     write(iUnitSlice, StrFmt1) 'MLT: ', vmltcells

     ! Find location of grid L that is within radius of slice.
     do i=2, nthetacells
        if( (radiusSlice-vrcells(i)) < 0.0 ) then
           iL=i-1
           exit
        end if
     end do
     
     ! Stop code on invalid radii.
     if(iL<1) call CON_stop(NameSub// &
          ': Selected radius outside of code boundary')
     
     IsFirstWrite = .false.
  end if

  ! Linearly interpolate results to selected L.
  do j=1, nphicells
     density(j) = (mgridden(iL+1,j)-mgridden(iL,j)) / &
          (vrcells(iL+1)-vrcells(iL)) * &
          (radiusSlice-vrcells(iL)) + mgridden(iL,j)
  end do

  ! Use TimeType to get a well-formatted date/time string.
  TimeNow%Time = CurrentTime
  call time_real_to_int(TimeNow)

  ! Write data to file.
  write(StrFmt2, '(a,i3.3,a)') '(i5,5i3,1x,i3.3,', nphicells, 'f9.3)'
  write(iUnitSlice, StrFmt2) &
       TimeNow%iYear, TimeNow%iMonth, TimeNow%iDay, &
       TimeNow%iHour, TimeNow%iMinute, TimeNow%iSecond, &
       floor(TimeNow%FracSecond*1000.0), density/100.0**3

end subroutine write_lslice

!=============================================================================
subroutine write_mltslice

  use ModIoDGCPM,     ONLY: cOutputDir, nMltSlice, iUnitMlt
  use ModTimeDGCPM,   ONLY: CurrentTime
  use ModSizeDGCPM,   ONLY: nrcells, nphicells
  use ModMainDGCPM,   ONLY: vrcells, mgridden, vmltcells
  use ModTimeDGCPM,   ONLY: CurrentTime, StartTime
  use ModIoUnit,      ONLY: io_unit_new
  use ModTimeConvert, ONLY: TimeType, time_real_to_int

  implicit none

  ! Saved variables: file names, MLT indexes, etc.
  logical, save :: IsInitiated = .false.
  integer, save :: iMltIndex(nPhiCells) ! Up to one slice per MLT.

  ! Temporary vars:
  integer :: i
  real    :: MltNow = 0.0
  character(len=1000) :: StringHeader
  character(len=12)   :: StrFmt1
  character(len=24)   :: StrFmt2
  character(len=100)  :: NameMltFiles(nMltSlice)
  type(TimeType)      :: TimeNow

  ! Testing params:
  logical :: DoTest, DoTestMe
  character(len=*), parameter :: NameSub='write_mltslice'
  !------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  ! Initiate files, write headers.
  if(.not. IsInitiated) then
     ! Slices should align with grid.
     if(mod(nPhiCells, nMltSlice)>0)then
        write(*,*)'DGCPM ERROR: nPhiCells not multiple of nMltSlice.'
        write(*,*)'             Please change PARAM::#MLTSLICE to fix this.'
        write(*,*)'             nPhiCells, nMltSlice = ', nPhiCells, nMltSlice
        call CON_stop(NameSub//'Bad nMltSlice choice.')
     end if

     ! Create space for file units.
     allocate(iUnitMlt(nMltSlice))

     ! Create common header line.
     write(StrFmt1, '(a,i3,a)') '(a, ', nrcells, 'f7.3)'
     write(StringHeader,StrFmt1) 'L = ', vrcells

     ! Find slice indices, get file names, write headers.
     do i=1, nMltSlice
        iMltIndex(i) = (nPhiCells/nMltSlice) * (i-1) + 1
        MltNow = vmltcells(iMltIndex(i))
        write(NameMltFiles(i), "(a,'/MLT_',i2.2,'_',i2.2,'_t',i10.10,'.txt')") &
             cOutputDir, floor(MltNow), floor(100.0*(MltNow-floor(MltNow))), &
             floor(CurrentTime-StartTime)
        iUnitMlt(i) = io_unit_new()
        if(DoTestMe) write(*,*)'DGCPM: Opening file ', NameMltFiles(i)
        open(iUnitMlt(i), FILE=NameMltFiles(i), STATUS='REPLACE')
        write(iUnitMlt(i), '(a, f5.2)') 'Density (cm-3) at MLT = ', MltNow
        write(iUnitMlt(i), '(a)') trim(StringHeader)             
     end do
     IsInitiated = .true.
  end if

  ! Use TimeType to get a well-formatted date/time string.
  TimeNow%Time = CurrentTime
  call time_real_to_int(TimeNow)

  ! Write data to file.
  write(StrFmt2, '(a,i3.3,a)') '(i5,5i3,1x,i3.3,', nrcells, 'f9.3)'
  do i=1, nMltSlice
     write(iUnitMlt(i), StrFmt2) &
          TimeNow%iYear, TimeNow%iMonth, TimeNow%iDay, &
          TimeNow%iHour, TimeNow%iMinute, TimeNow%iSecond, &
          floor(TimeNow%FracSecond*1000.0), mgridden(:,iMltIndex(i))/100.0**3
  end do

end subroutine write_mltslice

!=============================================================================
