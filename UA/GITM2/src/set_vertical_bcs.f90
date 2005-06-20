!\
! ------------------------------------------------------------
! set_boundary
! ------------------------------------------------------------
!/

subroutine set_vertical_bcs(LogRho,LogNS,Vel_GD,Temp, LogINS, iVel, VertVel)

  ! Fill in ghost cells at the top and bottom

  use ModSizeGitm, only: nAlts
  use ModPlanet, only: nSpecies, nIonsAdvect, Mass, nIons, IsEarth
  use ModGITM, only: gravity, TempUnit, dAlt, iEast_, iNorth_, iUp_, Altitude
  use ModInputs
  use ModConstants
  use ModTime, only: UTime, iJulianDay
  use ModVertical, only: Lat, Lon
  implicit none

  real, intent(inout) :: &
       LogRho(-1:nAlts+2), &
       LogNS(-1:nAlts+2,nSpecies), &
       LogINS(-1:nAlts+2,nIonsAdvect), &
       Vel_GD(-1:nAlts+2,3), &
       IVel(-1:nAlts+2,3), &
       Temp(-1:nAlts+2), &
       VertVel(-1:nAlts+2,nSpecies)

  integer :: iSpecies, iAlt
  real    :: InvScaleHeightS, InvScaleHgt, Alt, Lst, Ap = 4.0, dn
  logical :: IsFirstTime = .true., UseMsisBCs = .false.

  integer, dimension(25) :: sw

  if (IsEarth) UseMsisBCs = UseMsis

  if (IsFirstTime .and. UseMsisBCs) then
     call meter6(.true.)
     sw = 1
     IsFirstTime = .true.
  endif

  if (UseMsisBCs) then
     do iAlt = -1, 0
        Alt = Altitude(iAlt)/1000.0
        Lst = mod(UTime/3600.0+Lon/15.0,24.0)
        call msis_bcs(iJulianDay,UTime,Alt,Lat,Lon,Lst, &
             F107A,F107,AP,LogNS(iAlt,:), Temp(iAlt), &
             LogRho(iAlt))
     enddo
  else

     !  ! Bottom
     !
     !  do iSpecies=1,nIonsAdvect
     !     LogINS(0,iSpecies) = LogINS(1,iSpecies) - &
     !          dAlt(1)* Gravity(1)/Temp(1)
     !     LogINS(-1,iSpecies) = LogINS(1,iSpecies) - &
     !          2*dAlt(1)* Gravity(1)/Temp(1)
     !  enddo
     !
     !  ! Bottom
     !  ! Density and temperature are fixed
     !
     !  Temp(0)  = max(Temp(1),TempMin/TempUnit)
     !  Temp(-1) = max(Temp(1),TempMin/TempUnit)

  endif

  ! Slipping wall condition
  Vel_GD(-1:0,iEast_)  = 0.0
  Vel_GD(-1:0,iNorth_) = 0.0
!  Vel_GD(-1:0,iEast_)  = Vel_GD(1,iEast_)
!  Vel_GD(-1:0,iNorth_) = Vel_GD(1,iNorth_)
 ! Vel_GD( 0,iUp_)      = - Vel_GD(1,iUp_)
 ! Vel_GD(-1,iUp_)      = - Vel_GD(2,iUp_)

 Vel_GD( 0,iUp_)      = 0.0
 Vel_GD(-1,iUp_)      = 0.0

!!$ Vel_GD( 0,iUp_)      = Vel_GD( 1,iUp_)
!!$ Vel_GD(-1,iUp_)      = Vel_GD( 1,iUp_)
!!$
!!$ VertVel(0,:)  = VertVel(1,:)
!!$ VertVel(-1,:) = VertVel(1,:)

 VertVel(0,:)  = 0
 VertVel(-1,:) = 0


 IVel( 0,iUp_)      = 0.0
 IVel(-1,iUp_)      = 0.0

  ! Top

  Vel_GD(nAlts+1:nAlts+2,iEast_)  = Vel_GD(nAlts,iEast_)
  Vel_GD(nAlts+1:nAlts+2,iNorth_) = Vel_GD(nAlts,iNorth_)

  IVel(nAlts+1:nAlts+2,iEast_)  = IVel(nAlts,iEast_)
  IVel(nAlts+1:nAlts+2,iNorth_) = IVel(nAlts,iNorth_)

!  Vel_GD(nAlts+1,iUp_) = -Vel_GD(nAlts  ,iUp_)
!  Vel_GD(nAlts+2,iUp_) = -Vel_GD(nAlts-1,iUp_)

  if(Vel_GD(nAlts,iUp_)>0.)then

     Vel_GD(nAlts+1:nAlts+2,iUp_) = Vel_GD(nAlts,iUp_)
     IVel(nAlts+1:nAlts+2,iUp_)   = IVel(nAlts,iUp_)

     VertVel(nAlts+1,:) = VertVel(nAlts,:)
     VertVel(nAlts+2,:) = VertVel(nAlts,:)

  else
     ! Vel_GD(nAlts+1:nAlts+2,iUp_) = 0.0 ! -Vel(nAlts)
     Vel_GD(nAlts+1,iUp_) = -Vel_GD(nAlts,iUp_)
     Vel_GD(nAlts+2,iUp_) = -Vel_GD(nAlts-1,iUp_)

     VertVel(nAlts+1,:) = -VertVel(nAlts,:)
     VertVel(nAlts+2,:) = -VertVel(nAlts-1,:)

  endif

!  if (UseConduction) then
!     Temp(nAlts+1) = max(Temp(nAlts), TempMin/TempUnit(1,1,1))
!     Temp(nAlts+2) = max(Temp(nAlts), TempMin/TempUnit(1,1,1))
!  else
!     Temp(nAlts+1) = TempMax/TempUnit(1,1,nalts+1)
!     Temp(nAlts+2) = TempMax/TempUnit(1,1,nalts+2)
!  endif


!!!! CHANGE !!!!
  Temp(nAlts+1) = Temp(nAlts)
  Temp(nAlts+2) = Temp(nAlts-1)

  do iAlt = nAlts+1, nAlts+2
     InvScaleHgt  =  &
          -(Gravity(iAlt-1)+Gravity(iAlt)) / 2 / Temp(iAlt)
     LogRho(iAlt) = &
          LogRho(iAlt-1)-(Altitude(iAlt)-Altitude(iAlt-1))*InvScaleHgt
  enddo

  do iSpecies=1,nIonsAdvect
     dn = (LogINS(nAlts,iSpecies) - LogINS(nAlts-1,iSpecies))
     if (dn > 0.0) dn = -dn*1.0
!     LogINS(nAlts+1,iSpecies) = max( &
!          LogINS(nAlts,iSpecies) + dn, LogINS(nAlts,iSpecies) / 2.0)
!     LogINS(nAlts+2,iSpecies) = max( &
!          LogINS(nAlts+1,iSpecies) + dn, LogINS(nAlts+1,iSpecies) / 2.0)
     LogINS(nAlts+1,iSpecies) = LogINS(nAlts,iSpecies) + dn
     LogINS(nAlts+2,iSpecies) = LogINS(nAlts+1,iSpecies) + dn
  enddo

  do iSpecies=1,nSpecies
     do iAlt = nAlts+1, nAlts+2
        InvScaleHeightS = -Gravity(iAlt) * &
             Mass(iSpecies) / (Temp(iAlt)*Boltzmanns_Constant)
!!!! CHANGE !!!!
!             Mass(iSpecies) / (Temp(iAlt)*TempUnit(1,1,iAlt)*Boltzmanns_Constant)
        LogNS(iAlt,iSpecies) = &
             LogNS(iAlt-1,iSpecies) &
             -(Altitude(iAlt)-Altitude(iAlt-1))*InvScaleHeightS
        if (LogNS(nAlts+1,iSpecies) > 75.0 .or. &
             LogNS(nAlts+2,iSpecies) > 75.0) then
           write(*,*) "======> bcs : ", iSpecies, 1.0e-3/InvScaleHeightS, &
                Gravity(nAlts), Mass(iSpecies), Temp(nAlts), &
                LogNS(nAlts,iSpecies), LogNS(nAlts+1,iSpecies), &
                dAlt(nAlts), LogNS(nAlts+2,iSpecies)
        endif
     enddo
  enddo

end subroutine set_vertical_bcs

