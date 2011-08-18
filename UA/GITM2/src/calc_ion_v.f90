subroutine calc_ion_v(iBlock)

  use ModGITM
  use ModInputs
  use ModConstants

  implicit none

  integer, intent(in) :: iBlock
  integer :: iLon, iLat, iAlt
  integer :: imax, jmax, kmax, iError, iDir
  real    :: maxi, MaxVParallel

  real, dimension(1:nLons,1:nLats,1:nAlts) ::           &
                  B02, ForceDotB, Nie, RhoNu, IRho, &
                  VIParallel, VNParallel, gDotB, gpDotB, UDotB

  real, dimension(1:nLons, 1:nLats, 1:nAlts, 3) ::           &
                  PressureGradient, Force, BLocal, &
                  ForceCrossB, ForcePerp

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2):: Pressure_G
  !---------------------------------------------------------------------------

  call report("Ion Forcing Terms",1)
  call start_timing("Ion Forcing")

  IVelocity(:,:,:,:,iBlock) = 0.0

  if (iDebugLevel > 4) write(*,*) "=====> pressure gradient", iproc

  Pressure_G = IPressure(:,:,:,iBlock)+ePressure(:,:,:,iBlock)
  call UAM_Gradient(Pressure_G, PressureGradient, iBlock)

  PressureGradient(:,:,nAlts,iUp_) = PressureGradient(:,:,nAlts-1,iUp_)

  if (Is1D) then
     PressureGradient(:,:,:,iEast_) = 0.0
     PressureGradient(:,:,:,iNorth_) = 0.0
  endif

  Force = 0.0

  IRho = IDensityS(1:nLons,1:nLats,1:nAlts,ie_,iBlock) * &
       MeanIonMass(1:nLons,1:nLats,1:nAlts)

  if (UseIonPressureGradient) Force = Force - PressureGradient

  if (UseIonGravity) then
     do iAlt = 1, nAlts
        Force(:,:,iAlt,iUp_) = Force(:,:,iAlt,iUp_) + &
             IRho(:,:,iAlt) * Gravity_GB(1:nLons,1:nLats,iAlt,iBlock)
     enddo
  endif

  Nie = IDensityS(1:nLons,1:nLats,1:nAlts,ie_,iBlock) * Element_Charge

  BLocal = B0(1:nLons,1:nLats,1:nAlts,1:3,iBlock)
  B02 = B0(1:nLons,1:nLats,1:nAlts,iMag_,iBlock)**2

  if (UseExB) then
     do iDir = 1, 3
        Force(:,:,:,iDir) = Force(:,:,:,iDir) + &
             Nie * EField(1:nLons,1:nLats,1:nAlts,iDir)
     enddo
  endif

  RhoNu = IRho * Collisions(1:nLons,1:nLats,1:nAlts,iVIN_)

  if (UseNeutralDrag) then
     do iDir = 1, 3
        Force(:,:,:,iDir) = Force(:,:,:,iDir) + &
             RhoNu * Velocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock)
     enddo
  endif

  ForceDotB = sum(Force * BLocal, dim=4)

  do iDir = 1, 3
     ForcePerp(:,:,:,iDir) = Force(:,:,:,iDir) - &
          Force(:,:,:,iDir) * B0(1:nLons,1:nLats,1:nAlts,iDir,iBlock) / &
          B0(1:nLons,1:nLats,1:nAlts,iMag_,iBlock)
  enddo

  VIParallel = 0.0
  VNParallel = 0.0

  if (maxval(blocal) == 0) then

     IVelocity(1:nLons,1:nLats,1:nAlts,iUp_,iBlock) = &
          Velocity(1:nLons,1:nLats,1:nAlts,iUp_,iBlock) + &
          (Gravity_GB(1:nLons, 1:nLats, 1:nAlts, iBlock) - &
          (PressureGradient(1:nLons,1:nLats,1:nAlts,iUp_) / IRho) / &
          Collisions(1:nLons,1:nLats,1:nAlts,iVIN_))

     IVelocity(1:nLons,1:nLats,1:nAlts,iEast_,iBlock) = &
          Velocity(1:nLons,1:nLats,1:nAlts,iEast_,iBlock) - &
          (PressureGradient(1:nLons,1:nLats,1:nAlts,iEast_) / IRho) / &
          Collisions(1:nLons,1:nLats,1:nAlts,iVIN_)

     IVelocity(1:nLons,1:nLats,1:nAlts,iNorth_,iBlock) = &
          Velocity(1:nLons,1:nLats,1:nAlts,iNorth_,iBlock) - &
          (PressureGradient(1:nLons,1:nLats,1:nAlts,iNorth_) / IRho) / &
          Collisions(1:nLons,1:nLats,1:nAlts,iVIN_)
         
  else

  UDotB = sum(Velocity(1:nLons,1:nLats,1:nAlts,:,iBlock) * BLocal, dim=4)/ &
       B0(1:nLons,1:nLats,1:nAlts,iMag_,iBlock)
  gpDotB = sum(PressureGradient(1:nLons,1:nLats,1:nAlts,:) * &
       BLocal, dim=4) / B0(1:nLons,1:nLats,1:nAlts,iMag_,iBlock)

  do iLon = 1,nLons
     do iLat = 1,nLats
        gDotB(iLon,iLat,:) = Gravity_GB(iLon, iLat, 1:nAlts, iBlock) &
             * BLocal(iLon,iLat,1:nAlts,iUp_) &
             /     B0(iLon,iLat,1:nAlts,iMag_,iBlock)
     enddo
  enddo

  VIParallel = UDotB + &
       ( gDotB - gpDotB / IRho) / Collisions(1:nLons,1:nLats,1:nAlts,iVIN_)

!  write(*,*) VIParallel(1,1,49)

!  do iDir = 1, 3
!
!     VIParallel = VIParallel + &
!          Nie**2/RhoNu * ForceDotB * BLocal(:,:,:,iDir) &
!          / (RhoNu**2 + Nie**2 * B02)
!
!     write(*,*) "Force: ", iDir, ForceDotB(1,1,40), VIParallel(1,1,40),&
!          BLocal(1,1,40,iDir)/(RhoNu(1,1,40)**2 + Nie(1,1,40)**2 * B02(1,1,40))
!
!     if (UseNeutralDrag) then
!        VNParallel = VNParallel + &
!             Velocity(1:nLons,1:nLats,1:nAlts,iDir,iBlock) * &
!             BLocal(:,:,:,iDir) / B0(1:nLons,1:nLats,1:nAlts,iMag_,iBlock)
!     endif
!  enddo

  ! Let's limit the Parallel Flow to something reasonable...

  MaxVParallel = 100.0

  if (UseNeutralDrag) then
     VIParallel = min( UDotB + MaxVParallel, VIParallel)
     VIParallel = max( UDotB - MaxVParallel, VIParallel)
  else
     VIParallel = min( MaxVParallel, VIParallel)
     VIParallel = max(-MaxVParallel, VIParallel)
  endif

  ForceCrossB(:,:,:,iEast_) = &
       Force(:,:,:,iNorth_) * BLocal(:,:,:,iUp_) - &
       Force(:,:,:,iUp_)    * BLocal(:,:,:,iNorth_)

  ForceCrossB(:,:,:,iNorth_) = &
       Force(:,:,:,iUp_)    * BLocal(:,:,:,iEast_) - &
       Force(:,:,:,iEast_)  * BLocal(:,:,:,iUp_)

  ForceCrossB(:,:,:,iUp_)    = &
       Force(:,:,:,iEast_)  * BLocal(:,:,:,iNorth_) - &
       Force(:,:,:,iNorth_) * BLocal(:,:,:,iEast_)

  do iDir = 1, 3
     IVelocity(1:nLons,1:nLats,1:nAlts,iDir, iBlock) = &
          VIParallel*BLocal(:,:,:,iDir)/&
          B0(1:nLons,1:nLats,1:nAlts,iMag_,iBlock) + &
          ( RhoNu * ForcePerp(:,:,:,iDir) &
          + Nie * ForceCrossB(:,:,:,iDir) &
          ) / (RhoNu**2 + Nie**2 * B02)
  enddo

endif

  IVelocity(:,:,:,:,iBlock) = min( 3000.0, IVelocity(:,:,:,:,iBlock))
  IVelocity(:,:,:,:,iBlock) = max(-3000.0, IVelocity(:,:,:,:,iBlock))

  call end_timing("Ion Forcing")

  if (iDebugLevel > 4) write(*,*) "=====> done with calc_ion_v", iproc

end subroutine calc_ion_v
