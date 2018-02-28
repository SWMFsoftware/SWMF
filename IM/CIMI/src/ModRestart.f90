Module ModCimiRestart
  implicit none
  
  SAVE

  logical :: IsRestart = .false.
  real    :: DtSaveRestart=-1.0

  character(len=100), public:: NameRestartInDir="IM/restartIN/"
  character(len=100), public:: NameRestartOutDir="IM/restartOUT/"
contains
  !============================================================================

  subroutine cimi_read_restart
    use ModCimiPlanet,ONLY: nspec
    use ModCimiGrid,  ONLY: np,nt,nm,nk,neng,d4Element_C
    use ModCimi,      ONLY: f2, phot, Pressure_IC, PressurePar_IC, FAC_C, &
         Ppar_IC, Bmin_C, &
         eTimeAccumult_ICI, eChangeOperator_VICI, eChangeGlobal, &
         pTimeAccumult_ICI, pChangeOperator_VICI, &
         driftin, driftout, rbsumGlobal, rcsumGlobal, nOperator
    use ModCimiTrace, ONLY: iba
    use ModGmCimi,    ONLY: Den_IC
    use ModIoUnit,    ONLY: UnitTmp_
    use ModUtilities, ONLY: open_file, close_file
    use ModCimiGrid,  ONLY: iProc,nProc,iComm
    use ModMpi

    integer :: iError

    character(len=*), parameter:: NameSub = 'cimi_read_restart'    
    !--------------------------------------------------------------------------
    !When nProc>1, proc0 reads and then bcasts restart infor
    !when only 1 proc is used then just read restart info
    if(iProc == 0)then
       call open_file(file=trim(NameRestartInDir)//'data.restart',&
            status='old',form='unformatted', NameCaller=NameSub)

       read(UnitTmp_) f2  
       read(UnitTmp_) Den_IC
       read(UnitTmp_) phot
       read(UnitTmp_) Pressure_IC           
       read(UnitTmp_) PressurePar_IC           
       read(UnitTmp_) FAC_C             
       read(UnitTmp_) iba
       read(UnitTmp_) Ppar_IC
       read(UnitTmp_) Bmin_C
       read(UnitTmp_) eTimeAccumult_ICI
       read(UnitTmp_) eChangeOperator_VICI
       read(UnitTmp_) eChangeGlobal
       read(UnitTmp_) pTimeAccumult_ICI
       read(UnitTmp_) pChangeOperator_VICI
       read(UnitTmp_) rbsumGlobal
       read(UnitTmp_) rcsumGlobal
       read(UnitTmp_) driftin
       read(UnitTmp_) driftout
       call close_file
    end if
    
    if(nProc>1)then
       !Bcast from proc 0 to all procs
       call MPI_bcast(f2,nspec*np*nt*nm*nk, MPI_REAL, 0, iComm, iError)
       call MPI_bcast(Den_IC, nspec*np*nt, MPI_REAL, 0, iComm, iError)
       call MPI_bcast(phot, nspec*np*nt, MPI_REAL, 0, iComm, iError)
       call MPI_bcast(Pressure_IC, nspec*np*nt, MPI_REAL, 0, iComm, iError)
       call MPI_bcast(PressurePar_IC, nspec*np*nt, MPI_REAL, 0, iComm, iError)
       call MPI_bcast(FAC_C, np*nt, MPI_REAL, 0, iComm, iError)
       call MPI_bcast(iba, nt, MPI_INTEGER, 0, iComm, iError)
       call MPI_bcast(Ppar_IC, nspec*np*nt, MPI_REAL, 0, iComm, iError)
       call MPI_bcast(Bmin_C, np*nt, MPI_REAL, 0, iComm, iError)
       call MPI_bcast(eTimeAccumult_ICI, nspec*np*nt*(neng+2), &
            		MPI_REAL, 0, iComm, iError)
       call MPI_bcast(eChangeOperator_VICI, nspec*np*nt*nOperator*(neng+2), &
            		MPI_REAL, 0, iComm, iError)
       call MPI_bcast(eChangeGlobal, nspec * nOperator, &
            		MPI_REAL, 0, iComm, iError)
       call MPI_bcast(pTimeAccumult_ICI, nspec*np*nt*(neng+2), &
            		MPI_REAL, 0, iComm, iError)
       call MPI_bcast(pChangeOperator_VICI, nspec*np*nt*nOperator*(neng+2), &
            		MPI_REAL, 0, iComm, iError)
       call MPI_bcast(rbsumGlobal, nspec, MPI_REAL, 0, iComm, iError)
       call MPI_bcast(rcsumGlobal, nspec, MPI_REAL, 0, iComm, iError)
       call MPI_bcast(driftin, nspec, MPI_REAL, 0, iComm, iError)
       call MPI_bcast(driftout, nspec, MPI_REAL, 0, iComm, iError)
    endif
    
  end subroutine cimi_read_restart

  !============================================================================
  subroutine cimi_write_restart
    use ModCimiPlanet,ONLY: nspec
    use ModCimiGrid,  ONLY: np,nt,nm,nk,neng,MinLonPar,MaxLonPar
    use ModCimi,      ONLY: f2,time, phot, Pressure_IC, PressurePar_IC, &
         FAC_C, Ppar_IC, Bmin_C, &
         eTimeAccumult_ICI, eChangeOperator_VICI, eChangeGlobal, &
         pTimeAccumult_ICI, pChangeOperator_VICI, &
         driftin, driftout, rbsumGlobal, rcsumGlobal, nOperator
    use ModCimiTrace,ONLY: iba    
    use ModGmCimi,    ONLY: Den_IC
    use ModIoUnit,    ONLY: UnitTmp_
    use ModUtilities, ONLY: open_file, close_file
    use ModCimiGrid,  ONLY: iProc,nProc,iComm,nLonPar,nLonPar_P,nLonBefore_P
    use ModImSat,     ONLY: DoWriteSats,nImSats,NameSat_I,SatLoc_3I
    use ModMpi

    integer ::iSendCount, iSpecies, iM, iK, iError, iOperator, iEnergy
    real :: BufferSend_C(np,nt),BufferRecv_C(np,nt)
    integer :: BufferSend_I(nt),BufferRecv_I(nt)
    integer,allocatable :: iRecieveCount_P(:),iDisplacement_P(:)

    ! Initialize necessary variables for saving restart.sat file
    integer :: iSat,iRow

    character(len=*), parameter:: NameSub = 'cimi_write_restart'
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
                call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), &
                     iSendCount, &
                     MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P, &
                     MPI_REAL, 0, iComm, iError)
                if (iProc==0) f2(iSpecies,:,:,iM,iK)=BufferRecv_C(:,:)
             enddo
          enddo
       enddo
       
       !gather FAC
       BufferSend_C(:,:)=FAC_C(:,:)
       call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
            MPI_REAL, BufferRecv_C, iRecieveCount_P, iDisplacement_P, &
            MPI_REAL, 0, iComm, iError)
       if (iProc==0) FAC_C(:,:)=BufferRecv_C(:,:)

       !gather density and pressures
       do iSpecies=1,nspec
          BufferSend_C(:,:)=Den_IC(iSpecies,:,:)
          call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
               MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P, &
               MPI_REAL, 0, iComm, iError)
          if (iProc==0) Den_IC(iSpecies,:,:)=BufferRecv_C(:,:)

          BufferSend_C(:,:)=Pressure_IC(iSpecies,:,:)
          call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
               MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P, &
               MPI_REAL, 0, iComm, iError)
          if (iProc==0) Pressure_IC(iSpecies,:,:)=BufferRecv_C(:,:)

          BufferSend_C(:,:)=PressurePar_IC(iSpecies,:,:)
          call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
               MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P, &
               MPI_REAL, 0, iComm, iError)
          if (iProc==0) PressurePar_IC(iSpecies,:,:)=BufferRecv_C(:,:)
          
          BufferSend_C(:,:)=phot(iSpecies,:,:)
          call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
               MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P, &
               MPI_REAL, 0, iComm, iError)
          if (iProc==0) phot(iSpecies,:,:)=BufferRecv_C(:,:)

          BufferSend_C(:,:)=Ppar_IC(iSpecies,:,:)
          call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
               MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P,&
               MPI_REAL, 0, iComm, iError)
          if (iProc==0) Ppar_IC(iSpecies,:,:)=BufferRecv_C(:,:)

          do iEnergy = 1, neng + 2
                
             !gather eTimeAccumult_ICI
             BufferSend_C(:,:)=&
                  eTimeAccumult_ICI(iSpecies,:,:,iEnergy)
             call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
                  MPI_REAL, BufferRecv_C, iRecieveCount_P, iDisplacement_P, &
                  MPI_REAL, 0, iComm, iError)
             if (iProc==0) &
                  eTimeAccumult_ICI(iSpecies,:,:,iEnergy)=BufferRecv_C(:,:)

             !gather pTimeAccumult_ICI
             BufferSend_C(:,:)=&
                  pTimeAccumult_ICI(iSpecies,:,:,iEnergy)
             call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
                  MPI_REAL, BufferRecv_C, iRecieveCount_P, iDisplacement_P, &
                  MPI_REAL, 0, iComm, iError)
             if (iProc==0) &
                  pTimeAccumult_ICI(iSpecies,:,:,iEnergy)=BufferRecv_C(:,:)

             do iOperator = 1, nOperator

                !gather eChangeOperator_VICI
                BufferSend_C(:,:)=&
                     eChangeOperator_VICI(iSpecies,:,:,iEnergy,iOperator)
                call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), &
                     iSendCount, &
                     MPI_REAL, BufferRecv_C, iRecieveCount_P, iDisplacement_P, &
                     MPI_REAL, 0, iComm, iError)
                if (iProc==0) &
                     eChangeOperator_VICI(iSpecies,:,:,iEnergy,iOperator)=&
                     	   BufferRecv_C(:,:)
                
                !gather pChangeOperator_VICI
                BufferSend_C(:,:)=&
                     pChangeOperator_VICI(iSpecies,:,:,iEnergy,iOperator)
                call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), &
                     iSendCount, &
                     MPI_REAL, BufferRecv_C, iRecieveCount_P, iDisplacement_P, &
                     MPI_REAL, 0, iComm, iError)
                if (iProc==0) &
                     pChangeOperator_VICI(iSpecies,:,:,iEnergy,iOperator)=&
                     	   BufferRecv_C(:,:)
                
             enddo ! end iOperator loop
          enddo ! end iEnergy loop
       enddo ! end iSpecies loop
       
       !gather iba
       BufferSend_I(:)=iba(:)
       call MPI_GATHERV(BufferSend_I(MinLonPar:MaxLonPar), nLonPar, &
            MPI_INTEGER, BufferRecv_I, nLonPar_P, nLonBefore_P, &
            MPI_INTEGER, 0, iComm, iError)
       if (iProc==0) iba(:)=BufferRecv_I(:)

       !gather Bmin
       BufferSend_C(:,:)=Bmin_C(:,:)
       call MPI_GATHERV(BufferSend_C(:,MinLonPar:MaxLonPar), iSendCount, &
            MPI_REAL, BufferRecv_C,iRecieveCount_P, iDisplacement_P, &
            MPI_REAL, 0, iComm, iError)
       if (iProc==0) Bmin_C(:,:)=BufferRecv_C(:,:)
    endif

    if(iProc==0) then
       call open_file(file=trim(NameRestartOutDir)//'data.restart', &
            form='unformatted', NameCaller=NameSub)
       write(UnitTmp_) f2
       write(UnitTmp_) Den_IC
       write(UnitTmp_) phot
       write(UnitTmp_) Pressure_IC
       write(UnitTmp_) PressurePar_IC
       write(UnitTmp_) FAC_C
       write(UnitTmp_) iba                
       write(UnitTmp_) Ppar_IC
       write(UnitTmp_) Bmin_C
       write(UnitTmp_) eTimeAccumult_ICI
       write(UnitTmp_) eChangeOperator_VICI
       write(UnitTmp_) eChangeGlobal
       write(UnitTmp_) pTimeAccumult_ICI
       write(UnitTmp_) pChangeOperator_VICI
       write(UnitTmp_) rbsumGlobal
       write(UnitTmp_) rcsumGlobal
       write(UnitTmp_) driftin
       write(UnitTmp_) driftout
       call close_file

       call open_file(file=trim(NameRestartOutDir)//'restart.H', &
            NameCaller=NameSub)
       write(UnitTmp_,'(a)') '#TIMESIMULATION'
       write(UnitTmp_,'(es15.8,a25)') time,'tSimulation'
       write(UnitTmp_,*)
       write(UnitTmp_,'(a)') '#RESTART'
       write(UnitTmp_,'(a)') 'T                       DoRestart'
       write(UnitTmp_,'(l1,a45)') DoWriteSats, 'DoReadRestartSatellite'
       call close_file

       if (DoWriteSats) then
          call open_file(file=trim(NameRestartOutDir)//'restart.sat', &
               form="unformatted", NameCaller=NameSub) 
          write(UnitTmp_) nImSats
       
          do iSat=1,nImSats
             
             write(UnitTmp_) NameSat_I(iSat)
             
          end do
       
          do iSat=1,nImSats
          
             do iRow=1,2
             
                write(UnitTmp_) SatLoc_3I(1:4,iRow,iSat)
                
             end do
             
          end do
          
          call close_file

       end if
       
    endif
  end subroutine cimi_write_restart

end Module ModCimiRestart
