Module ModCrcmRestart
  implicit none
  
  SAVE

  logical :: IsRestart = .false.
contains
  !============================================================================

  subroutine crcm_read_restart
    use ModCrcmPlanet,ONLY: nspec
    use ModCrcmGrid,  ONLY: np,nt,nm,nk
    use ModCrcm,      ONLY: f2, phot, Pressure_IC, FAC_C, Ppar_IC, Bmin_C
    use ModFieldTrace,ONLY: iba
    use ModGmCrcm,    ONLY: Den_IC
    use ModIoUnit,    ONLY: UnitTmp_
    use ModCrcmGrid,  ONLY: iProc,nProc,iComm
    use ModMpi

    integer :: iError,iSend
    !--------------------------------------------------------------------------
    !When nProc>1, proc0 reads and then bcasts restart infor
    !when only 1 proc is used then just read restart info
    if(iProc == 0)then
       open(unit=UnitTmp_,file='IM/restartIN/data.restart',status='old',form='unformatted')
       read(UnitTmp_) f2  
       read(UnitTmp_) Den_IC
       read(UnitTmp_) phot
       read(UnitTmp_) Pressure_IC           
       read(UnitTmp_) FAC_C             
       read(UnitTmp_) iba
       read(UnitTmp_) Ppar_IC
       read(UnitTmp_) Bmin_C
       close(UnitTmp_)
    end if

    if(nProc>1)then
       !Bcast from proc 0 to all procs
       call MPI_bcast(f2,nspec*np*nt*nm*nk, MPI_REAL, 0, iComm, iError)
       call MPI_bcast(Den_IC, nspec*np*nt, MPI_REAL, 0, iComm, iError)
       call MPI_bcast(phot, nspec*np*nt, MPI_REAL, 0, iComm, iError)
       call MPI_bcast(Pressure_IC, nspec*np*nt, MPI_REAL, 0, iComm, iError)
       call MPI_bcast(FAC_C, np*nt, MPI_REAL, 0, iComm, iError)
       call MPI_bcast(iba, nt, MPI_INTEGER, 0, iComm, iError)
       call MPI_bcast(Ppar_IC, nspec*np*nt, MPI_REAL, 0, iComm, iError)
       call MPI_bcast(Bmin_C, np*nt, MPI_REAL, 0, iComm, iError)
    endif

  end subroutine crcm_read_restart
  

  !============================================================================
  subroutine crcm_write_restart
    use ModCrcmPlanet,ONLY: nspec
    use ModCrcmGrid,  ONLY: np,nt,nm,nk,MinLonPar,MaxLonPar
    use ModCrcm,      ONLY: f2,time, phot, Pressure_IC, FAC_C, Ppar_IC, Bmin_C
    use ModFieldTrace,ONLY: iba    
    use ModGmCrcm,    ONLY: Den_IC
    use ModIoUnit,    ONLY: UnitTmp_
    use ModCrcmGrid,  ONLY: iProc,nProc,iComm,nLonPar,nLonPar_P,nLonBefore_P
    use ModMpi

    integer ::iSendCount, iSpecies, iM, iK, iError
    real :: BufferSend_C(np,nt),BufferRecv_C(np,nt)
    integer :: BufferSend_I(nt),BufferRecv_I(nt)
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
       BufferSend_C(:,:)=FAC_C(:,:)
       call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, MPI_REAL, &
            BufferRecv_C, iRecieveCount_P, iDisplacement_P, MPI_REAL, &
            0, iComm, iError)
       if (iProc==0) FAC_C(:,:)=BufferRecv_C(:,:)

       !gather density and pressures
       do iSpecies=1,nspec
          BufferSend_C(:,:)=Den_IC(iSpecies,:,:)
          call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
               MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P,MPI_REAL, &
               0, iComm, iError)
          if (iProc==0) Den_IC(iSpecies,:,:)=BufferRecv_C(:,:)

          BufferSend_C(:,:)=Pressure_IC(iSpecies,:,:)
          call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
               MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P,MPI_REAL, &
               0, iComm, iError)
          if (iProc==0) Pressure_IC(iSpecies,:,:)=BufferRecv_C(:,:)
          
          BufferSend_C(:,:)=phot(iSpecies,:,:)
          call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
               MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P,MPI_REAL, &
               0, iComm, iError)
          if (iProc==0) phot(iSpecies,:,:)=BufferRecv_C(:,:)

          BufferSend_C(:,:)=Ppar_IC(iSpecies,:,:)
          call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
               MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P,MPI_REAL, &
               0, iComm, iError)
          if (iProc==0) Ppar_IC(iSpecies,:,:)=BufferRecv_C(:,:)
       enddo

       !gather iba
       BufferSend_I(:)=iba(:)
       call MPI_GATHERV(BufferSend_I(MinLonPar:MaxLonPar), nLonPar, MPI_INTEGER, &
            BufferRecv_I, nLonPar_P, nLonBefore_P, MPI_INTEGER, 0, iComm, iError)
       if (iProc==0) iba(:)=BufferRecv_I(:)

       !gather Bmin
       BufferSend_C(:,:)=Bmin_C(:,:)
       call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
            MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P, MPI_REAL, &
            0, iComm, iError)
       if (iProc==0) Bmin_C(:,:)=BufferRecv_C(:,:)
    endif

    if(iProc==0) then
       open(unit=UnitTmp_,file='IM/restartOUT/data.restart',form='unformatted')
       write(UnitTmp_) f2
       write(UnitTmp_) Den_IC
       write(UnitTmp_) phot
       write(UnitTmp_) Pressure_IC
       write(UnitTmp_) FAC_C
       write(UnitTmp_) iba                
       write(UnitTmp_) Ppar_IC
       write(UnitTmp_) Bmin_C
       close(UnitTmp_)

       open(unit=UnitTmp_,file='IM/restartOUT/restart.H')
       write(UnitTmp_,'(a)') '#TIMESIMULATION'
       write(UnitTmp_,'(es15.8,a25)') time,'tSimulation'
       close(UnitTmp_)
    endif
  end subroutine crcm_write_restart

end Module ModCrcmRestart
