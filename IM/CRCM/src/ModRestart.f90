Module ModCrcmRestart
  implicit none
  
  SAVE

  logical :: IsRestart = .false.
contains
  !============================================================================

  subroutine crcm_read_restart
    use ModCrcmPlanet,ONLY: nspec
    use ModCrcmGrid,  ONLY: np,nt,nm,nk
    use ModCrcm,      ONLY: f2, phot, Pressure_IC, FAC_C
    use ModIoUnit,    ONLY: UnitTmp_
    use ModCrcmGrid,  ONLY: iProc,nProc,iComm
    use ModMpi

    integer :: iError,iSend
    !--------------------------------------------------------------------------
!    if (nProc>1) then
!       ! When nProc>1, proc0 reads and then bcasts restart infor
!       if (iProc == 0 ) then
!          open(unit=UnitTmp_,file='IM/restartIN/data.restart',status='old',form='unformatted')
!          read(UnitTmp_) f2             
!          read(UnitTmp_) phot             
!          read(UnitTmp_) Pressure_IC           
!          read(UnitTmp_) FAC_C             
!          close(UnitTmp_)
!       endif
!       !Bcast from proc 0 to all procs
!       call MPI_bcast(f2,nspec*np*nt*nm*nk,MPI_REAL,0,iComm,iError)
!       call MPI_bcast(phot,nspec*np*nt,MPI_REAL,0,iComm,iError)
!       call MPI_bcast(Pressure_IC,nspec*np*nt,MPI_REAL,0,iComm,iError)
!       call MPI_bcast(FAC_C,np*nt,MPI_REAL,0,iComm,iError)
!    else
       !when only 1 proc is used then just read restart info
       open(unit=UnitTmp_,file='IM/restartIN/data.restart',status='old',form='unformatted')
       read(UnitTmp_) f2             
       read(UnitTmp_) phot             
       read(UnitTmp_) Pressure_IC           
       read(UnitTmp_) FAC_C             
       close(UnitTmp_)
!    endif


  end subroutine crcm_read_restart
  

  !============================================================================
  subroutine crcm_write_restart
    use ModCrcmPlanet,ONLY: nspec
    use ModCrcmGrid,  ONLY: np,nt,nm,nk,MinLonPar,MaxLonPar
    use ModCrcm,      ONLY: f2,time, phot, Pressure_IC, FAC_C
    use ModIoUnit,    ONLY: UnitTmp_
    use ModCrcmGrid,  ONLY: iProc,nProc,iComm,nLonPar,nLonPar_P,nLonBefore_P
    use ModMpi

    integer ::iSendCount, iSpecies, iM, iK, iError
    real :: BufferSend_C(np,nt),BufferRecv_C(np,nt)
    integer,allocatable :: iRecieveCount_P(:),iDisplacement_P(:)
    !--------------------------------------------------------------------------

    ! When nProc>1 gather to proc 0 for writing.
    if (nProc>1) then
       if (.not.allocated(iRecieveCount_P)) &
            allocate(iRecieveCount_P(nProc), iDisplacement_P(nProc))       
       iSendCount = np*nLonPar
       iRecieveCount_P=np*nLonPar_P
       iDisplacement_P = np*nLonBefore_P
       !Gather f2
       do  iSpecies=1,nspec
          do iM=1,nm
             do iK=1,nK
                BufferSend_C(:,:)=f2(iSpecies,:,:,iM,iK)
                call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar),iSendCount, &
                     MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P, &
                     MPI_REAL, 0, iComm, iError)
                if (iProc==0) f2(iSpecies,:,:,iM,iK)=BufferRecv_C(:,:)
             enddo
          enddo
       enddo
       
       !gather FAC
       call MPI_GATHERV(FAC_C(:,MinLonPar:MaxLonPar), iSendCount, MPI_REAL, &
            FAC_C, iRecieveCount_P, iDisplacement_P, MPI_REAL, &
            0, iComm, iError)

       !gather pressure and phot
       do iSpecies=1,nspec
          BufferSend_C(:,:)=0.0
          BufferSend_C(:,:)=Pressure_IC(iSpecies,:,:)
          call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
               MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P,MPI_REAL, &
               0, iComm, iError)
          if (iProc==0) Pressure_IC(iSpecies,:,:)=BufferRecv_C(:,:)
          
          BufferSend_C(:,:)=0.0
          BufferSend_C(:,:)=phot(iSpecies,:,:)
          call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
               MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P,MPI_REAL, &
               0, iComm, iError)
          if (iProc==0) phot(iSpecies,:,:)=BufferRecv_C(:,:)
       enddo
    endif

    if (nProc>1) then
       !if nProc>1 only write out restart from proc 0
       if(iProc==0) then
          open(unit=UnitTmp_,file='IM/restartOUT/data.restart',form='unformatted')
          write(UnitTmp_) f2
          write(UnitTmp_) phot
          write(UnitTmp_) Pressure_IC
          write(UnitTmp_) FAC_C
          close(UnitTmp_)

          open(unit=UnitTmp_,file='IM/restartOUT/restart.H')
          write(UnitTmp_,'(a)') '#TIMESIMULATION'
          write(UnitTmp_,'(es15.8,a25)') time,'tSimulation'
          close(UnitTmp_)
       endif
    else
       !When nProc=1, write out restart
       open(unit=UnitTmp_,file='IM/restartOUT/data.restart',form='unformatted')
       write(UnitTmp_) f2
       write(UnitTmp_) phot
       write(UnitTmp_) Pressure_IC
       write(UnitTmp_) FAC_C
       close(UnitTmp_)
       
       open(unit=UnitTmp_,file='IM/restartOUT/restart.H')
       write(UnitTmp_,'(a)') '#TIMESIMULATION'
       write(UnitTmp_,'(es15.8,a25)') time,'tSimulation'
       close(UnitTmp_)
    endif
  end subroutine crcm_write_restart

end Module ModCrcmRestart
