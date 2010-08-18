
subroutine write_restart(dir)

  use ModGITM
  use ModInputs
  use ModTime
  use ModSphereInterface, only: iStartBLK

  implicit none

  character (len=*), intent(in) :: dir
  character (len=5) :: cBlock

  integer :: iBlock, iSpecies, i

  call report("write_restart",1)

  if (iProc == 0) then

     open(unit=iRestartUnit_, file=dir//"header.rst", status="unknown")

     write(iRestartUnit_, *) ""
     write(iRestartUnit_, '(a)') "#ISTEP"
     write(iRestartUnit_, *) iStep

     write(iRestartUnit_, *) ""
     write(iRestartUnit_, '(a)') "#TSIMULATION"
     write(iRestartUnit_, *) tSimulation

     write(iRestartUnit_, *) ""
     write(iRestartUnit_, '(a)') "#TIMESTART"
     do i=1,6
        write(iRestartUnit_, *) iStartTime(i)
     enddo

     write(iRestartUnit_, *) ""
     write(iRestartUnit_, '(a)') "#SPHERE"
     write(iRestartUnit_, *) Is1D
     write(iRestartUnit_, *) IsFullSphere
     write(iRestartUnit_, *) LatStart
     write(iRestartUnit_, *) LatEnd
     write(iRestartUnit_, *) LonStart

     write(iRestartUnit_, *) ""
     write(iRestartUnit_, '(a)') "#END"
     write(iRestartUnit_, *) ""

     close(iRestartUnit_)

  endif

  do iBlock = 1, nBlocks

     write(cBlock,'(a1,i4.4)') "b",iBlock+iStartBLK
     open(unit=iRestartUnit_, file=dir//"/"//cBlock//".rst", &
          status="unknown", form="unformatted")

     write(iRestartUnit_) Longitude(:,iBlock)
     write(iRestartUnit_) Latitude(:,iBlock)
     if(UseTopography)then
        write(iRestartUnit_) Altitude_GB(:,:,:,iBlock)
     else
        write(iRestartUnit_) Altitude_GB(1,1,:,iBlock)
     end if

     do iSpecies=1,nSpeciesTotal
        write(iRestartUnit_) NDensityS(:,:,:,iSpecies,iBlock)
     enddo

     do iSpecies=1,nIons
        write(iRestartUnit_) IDensityS(:,:,:,iSpecies,iBlock)
     enddo

     write(iRestartUnit_)  Temperature(:,:,:,iBlock)
     write(iRestartUnit_) ITemperature(:,:,:,iBlock)
     write(iRestartUnit_) eTemperature(:,:,:,iBlock)

     write(iRestartUnit_)  Velocity(:,:,:,:,iBlock)
     write(iRestartUnit_) IVelocity(:,:,:,:,iBlock)
     write(iRestartUnit_) VerticalVelocity(:,:,:,:,iBlock)

     close(iRestartUnit_)

  enddo

end subroutine write_restart

!=============================================================================

subroutine read_restart(dir)

  use ModGITM
  use ModInputs
  use ModTime
  use ModSphereInterface, only: iStartBLK

  implicit none

  character (len=*), intent(in) :: dir
  character (len=5) :: cBlock
  integer :: iBlock, i, iSpecies, iAlt
  !---------------------------------------------------------------------------
  call report("read_restart",1)

  do iBlock = 1, nBlocks

     write(cBlock,'(a1,i4.4)') "b",iBlock+iStartBLK
     if (iDebugLevel > 2) write(*,*) "===> Reading block ",cBlock

     open(unit=iRestartUnit_, file=dir//"/"//cBlock//".rst", &
          status="old", form="unformatted")

     if (iDebugLevel > 4) write(*,*) "=====> Reading Longitude"
     read(iRestartUnit_) Longitude(:,iBlock)
     if (iDebugLevel > 4) write(*,*) "=====> Reading Latitude"
     read(iRestartUnit_) Latitude(:,iBlock)
     if (iDebugLevel > 4) write(*,*) "=====> Reading Altitude"

     if(UseTopography)then
        read(iRestartUnit_) Altitude_GB(:,:,:,iBlock)
     else
        read(iRestartUnit_) Altitude_GB(1,1,:,iBlock)
        do iAlt = -1, nAlts+2
           Altitude_GB(:,:,iAlt,iBlock) = Altitude_GB(1,1,iAlt,iBlock)
        end do
     end if
!     read(iRestartUnit_) Altitude_GB(:,:,:,iBlock)

     do iSpecies=1,nSpeciesTotal
        if (iDebugLevel > 3) &
             write(*,*) "====> Reading Species",iSpecies, Mass(iSpecies)
        read(iRestartUnit_) NDensityS(:,:,:,iSpecies,iBlock)
     enddo

     Rho(:,:,:,iBlock) = 0.0
     NDensity(:,:,:,iBlock) = 0.0
     do iSpecies=1,nSpecies
        Rho(:,:,:,iBlock) = Rho(:,:,:,iBlock) + &
          Mass(iSpecies)*NDensityS(:,:,:,iSpecies,iBlock)
        NDensity(:,:,:,iBlock) = NDensity(:,:,:,iBlock) + &
          NDensityS(:,:,:,iSpecies,iBlock)
     enddo

     do iSpecies=1,nIons
        if (iDebugLevel > 4) &
             write(*,*) "=====> Reading Ion Species",iSpecies, Mass(iSpecies)
        read(iRestartUnit_) IDensityS(:,:,:,iSpecies,iBlock)
     enddo

     if (iDebugLevel > 4) write(*,*) "=====> Reading Temperature"
     read(iRestartUnit_)  Temperature(:,:,:,iBlock)
     
     if (isMars) then 
        SurfaceTemp = 0.0
        SurfaceTemp(:,:,iBlock) = Temperature(1:nLons,:1:nLats,0,iBlock)
     endif

     if (iDebugLevel > 4) write(*,*) "=====> Reading ITemperature"
     read(iRestartUnit_) ITemperature(:,:,:,iBlock)
     if (iDebugLevel > 4) write(*,*) "=====> Reading eTemperature"
     read(iRestartUnit_) eTemperature(:,:,:,iBlock)

     if (iDebugLevel > 4) write(*,*) "=====> Reading Velocity"
     read(iRestartUnit_)  Velocity(:,:,:,:,iBlock)
     if (iDebugLevel > 4) write(*,*) "=====> Reading IVelocity"
     read(iRestartUnit_) IVelocity(:,:,:,:,iBlock)
     if (iDebugLevel > 4) write(*,*) "=====> Reading VerticalVelocity"
     read(iRestartUnit_) VerticalVelocity(:,:,:,:,iBlock)

     close(iRestartUnit_)

  enddo

  if (.not. Is1D) call exchange_messages_sphere

end subroutine read_restart
