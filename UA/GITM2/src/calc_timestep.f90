
!\
! ------------------------------------------------------------
! calc_timestep
! ------------------------------------------------------------
!/

subroutine calc_timestep_horizontal

  use ModGITM
  use ModInputs
  use ModTime
  use ModMpi
  implicit none

  logical :: DoWarning = .true.

  real :: DtCond, DtLocal, DtHorizontal, DtEnd

  real :: cSound_H(0:nLons+1,0:nLats+1)

  integer :: iBlock, iError, iLat, iLon, iAlt

  call report("calc_timestep_horizontal",2)

  DtLocal = 1.0e32

  do iBlock = 1, nBlocks

     do iAlt = 1,nAlts

        ! Calculate maximum propagation speeds for the horizontal directions
        cSound_H = sqrt(Gamma * Temperature(0:nLons+1,0:nLats+1,iAlt,iBlock))

        cMax_GDB(:,:,iAlt,iNorth_,iBlock) = &
             abs(Velocity(0:nLons+1,0:nLats+1,iAlt,iNorth_,iBlock)) + cSound_H

        cMax_GDB(:,:,iAlt,iEast_,iBlock) = &
             abs(Velocity(0:nLons+1,0:nLats+1,iAlt,iEast_,iBlock)) + cSound_H

        ! Find stability limit on the time step
        do iLat = 1,nLats
           do iLon = 1,nLons
           
              DtLocal = min(DtLocal, Cfl / ( &
                   cMax_GDB(iLon, iLat, iAlt, iEast_,  iBlock) / &
                   dLonDist_GB(iLat,iAlt,iBlock) + &
                   cMax_GDB(iLon, iLat, iAlt, iNorth_, iBlock) / &
                   dLatDist_GB(iLat,iAlt,iBlock)))

              if (UseIonAdvection) &
                   DtLocal = min(DtLocal, Cfl / ( &
                   (abs(IVelocity(iLon, iLat, iAlt, iEast_,  iBlock))+ &
                   cSound_H(iLon, iLat))/ &
                   dLonDist_GB(iLat,iAlt,iBlock) + &
                   (abs(IVelocity(iLon, iLat, iAlt, iNorth_, iBlock))+ &
                   cSound_H(iLon, iLat))/ &
                   dLatDist_GB(iLat,iAlt,iBlock))/2.0)

           end do
        end do

     end do

  end do

  DtEnd = EndTime - CurrentTime
  DtLocal = min(DtLocal, DtEnd)

  call MPI_AllREDUCE(DtLocal, DtHorizontal, 1, MPI_REAL, MPI_MIN, &
       iCommGITM, iError)

  if (iDebugLevel > 2) &
       write(*,*) "===> DtVertical, DtHorizontal : ", Dt, DtHorizontal

  Dt = min(Dt, DtHorizontal)

end subroutine calc_timestep_horizontal

!==========================================================================

subroutine calc_timestep_vertical

  use ModGITM
  use ModInputs
  use ModTime
  use ModMpi
  implicit none

  logical :: DoWarning = .true.

  real :: DtCond, DtLocal, DtVertical, DtEnd

  real :: cm(1:nAlts)

  integer :: iBlock, iError, iLat, iLon, iAlt
  integer :: iLatSave, iLonSave

  call report("calc_timestep_vertical",2)

  DtLocal = 1.0e32

  iLatSave = 0
  iLonSave = 0

  do iBlock = 1, nBlocks

     call calc_rates(iBlock)

     do iLon = 1,nLons
        do iLat = 1,nLats
!           cMax_GDB(iLon,iLat,:,iUp_,iBlock) = &
!                abs(Velocity(iLon,iLat,0:nAlts+1,iUp_,iBlock)) + &
!                sqrt(Gamma * Temperature(iLon,iLat,0:nAlts+1,iBlock))

           do iAlt = 0, nAlts+1
              cMax_GDB(iLon,iLat,iAlt,iUp_,iBlock) = &
                   maxval(abs(VerticalVelocity(iLon,iLat,iAlt,:,iBlock))) + &
                   sqrt(Gamma * Temperature(iLon,iLat,iAlt,iBlock))
           enddo
!           write(*,*) "neutral : ", &
!                maxval(abs(Velocity(iLon,iLat,0:nAlts+1,iUp_,iBlock))), &
!                maxval(cMax_GDB(iLon,iLat,1:nAlts,iUp_,iBlock)), &
!                maxval(abs(Velocity(iLon,iLat,0:nAlts+1,iUp_,iBlock)) + &
!                sqrt(Gamma * Temperature(iLon,iLat,0:nAlts+1,iBlock)))
!
!           do iAlt=0,nAlts
!              write(*,*) "cmax : ",iLon,iLat,&
!                   cMax_GDB(iLon,iLat,iAlt,iUp_,iBlock), &
!                   abs(Velocity(iLon,iLat,iAlt,iUp_,iBlock)), &
!                   sqrt(Gamma * Temperature(iLon,iLat,iAlt,iBlock))
!           enddo

           DtLocal = min(DtLocal, &
                Cfl / &
                maxval(cMax_GDB(iLon,iLat,1:nAlts,iUp_,iBlock)/dAlt(1:nAlts)))

!write(*,*) "dt : ", dtlocal, maxval(cMax_GDB(iLon,iLat,1:nAlts,iUp_,iBlock)), &
!     minval(Temperature(iLon,iLat,0:nAlts+1,iBlock)),&
!     maxval(Temperature(iLon,iLat,0:nAlts+1,iBlock)), &
!     maxval(abs(Velocity(iLon,iLat,0:nAlts+1,iUp_,iBlock)))
!
!do iAlt =1, nAlts
!   write(*,*) Temperature(iLon,iLat,iAlt,iBlock), &
!        Velocity(iLon,iLat,iAlt,iUp_,iBlock)
!enddo

           if (UseIonAdvection) then

              cm = abs(IVelocity(iLon,iLat,1:nAlts,iUp_,iBlock)) + &
                sqrt(Gamma * Temperature(iLon,iLat,1:nAlts,iBlock))

!              write(*,*) "ion : ", &
!                   maxval(abs(IVelocity(iLon,iLat,1:nAlts,iUp_,iBlock))), &
!                   maxval(cm(1:nAlts))

              DtLocal = min(DtLocal, &
                   Cfl / &
                   maxval(cm/dAlt(1:nAlts)))

!write(*,*) "dt ion : ", dtlocal, maxval(cM)

           endif

!           if (UseConduction) then
!              DtCond = 0.25*minval(cp(iLon,iLat,1:nAlts,iBlock)) / &
!                   maxval((KappaTemp(iLon,iLat,1:nAlts,iBlock)) / &
!                   (Rho(iLon,iLat,1:nAlts,iBlock))/dAlt(1:nAlts)**2)
!              if (DtCond < DtLocal) then
!                 iLatSave = iLat
!                 iLonSave = iLon
!              endif
!              if(DoWarning .and. DtCond < DtLocal .and. iProc == 0)then
!                 write(*,*)'Reduced time step for conduction from Dt=',&
!                      Dt,' to DtCond=',DtCond
!                 DoWarning = .false.
!              endif
!              DtLocal = min(DtLocal, DtCond)
!           endif

        enddo
     enddo

  enddo

  DtEnd = EndTime - CurrentTime
  dTLocal = min(DtLocal, DtEnd)

  call MPI_AllREDUCE(dTLocal, DtVertical, 1, MPI_REAL, MPI_MIN, &
       iCommGITM, iError)

!  if (DtVertical == DtLocal .and. iLatSave > 0) then
!     write(*,*) "Conduction : ",iLatSave, iLonSave
!  endif

  Dt = min(Dt, DtVertical)

  if (iDebugLevel > 2) &
       write(*,*) "===> DtVertical : ", Dt

  if (dt < cfl/100.0) then
     write(*,*) "Dt too slow!!!", dt
     call stop_gitm("Stopping in calc_timestep_vertical")
  endif

end subroutine calc_timestep_vertical

