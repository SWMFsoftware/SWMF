!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
! File Name: dgcpm_test.f90
!
! Contains : Subroutines used to test the core functionality of
! the DGCPM model. Currently each test routine writes its own output
! file, but in the future they may be linked into a single output
! file that can be processed by IDL code automatically.
!
! Last Modified : January 2012, Aron Dodger

! *************************************************************************
!                                   TestFilling
! This subroutine's only use is to test the Filling subroutine and
! produce an output array that is the max filling and loss values as
! a result of the filling subroutine. This is done via a simple replacement
! hack for the mgridn and mgridden variables, and uses the mgridsource
! variable for its results. Values calculated are for closed fieldlines only.
! *************************************************************************

       subroutine testfilling(delt)

     use ModSizeDGCPM
     use ModMainDGCPM
     use ModTimeDGCPM
     use ModFunctionsDGCPM

     use ModIoUnit, ONLY: UnitTMP_
     
        implicit none

        real OldN(nrcells,nphicells), OldDen(nrcells,nphicells)
        real OldOC(nrcells,nphicells)
        real MaxFilling(nrcells,nphicells)
        real delt
        integer i, j

! Backup N, Den, and OC values. Sets all fieldlines to closed.

        do i=1,nrcells
            do j=1, nphicells
                OldN(i,j) = mgridn(i,j)
                OldDen(i,j) = mgridden(i,j)
                OldOC(i,j) = mgridoc(i,j)
            enddo
        enddo

! Begin Loop 1 - Dayside Fill calculation 
        do i=1, nrcells
            do j=1, nphicells
                mgridoc(i,j) = 1.
                mgridn(i,j) = 1.
            enddo
        enddo

        call filling(delt)

        do i=1, nrcells
            do j=1, nphicells
                if ((vphicells(j).GE.90.0).AND.(vphicells(j).LE.270.0)) then
                    MaxFilling(i,j) = mgridsource(i,j)
                endif
            enddo
        enddo

! Begin Loop 2 - Nightside Loss Calculation
    
        do i=1, nrcells
            do j=1, nphicells
                mgridn(i,j) = saturation(vrcells(i)) * mgridvol(i,j)
            enddo
        enddo

        call filling(delt)

        do i=1, nrcells
            do j=1, nphicells
                if ((vphicells(j).LT.90.).OR.(vphicells(j).GT.270.)) then
                    MaxFilling(i,j) = mgridsource(i,j)
                endif
            enddo
        enddo

! Write output file

        open(unit=UnitTmp_, file='Filling_Test.dat', form='formatted')
        write(UnitTmp_,'(F5.2)') Delt
        write(UnitTmp_,'(I3)') nrcells 
        write(UnitTmp_,'(I3)') nphicells
        write(UnitTmp_,*) vrcells
        write(UnitTmp_,*) mgridx
        write(UnitTmp_,*) mgridy
        write(UnitTmp_,*) MaxFilling
        close(unit = UnitTmp_)

! Fix mgridn and mgridden

        do i=1, nrcells
            do j=1, nphicells
                mgridn(i,j) = OldN(i,j)
                mgridden(i,j) = OldDen(i,j)
                mgridoc(i,j) = OldOC(i,j)
            enddo
        enddo

! End of Subroutine

        return
        End
