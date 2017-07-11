!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine chapman_integrals(iBlock)

  use ModGITM
  use ModEUV
  use ModPlanet
  use ModConstants
  use ModInputs, only : iDebugLevel
  implicit none

  integer, intent(in) :: iBlock

  integer :: iLon, iLat, iAlt, iSpecies, iiAlt
  real :: ScaleHeightS

  real, dimension(nLons, nLats, nAlts, nSpecies) :: &
       integrals, xp, erfcy

  real, dimension(nSpecies) :: MaxCh

  real :: a,b,c,d,f,g, y

  Chapman(:,:,:,:,iBlock) = 1.0e26
  Integrals = 0.0

  call report("chapman",2)
  call start_timing("chapman_integrals")

  !
  ! This is all from Smith and Smith, JGR 1972, vol. 77, page 3592
  ! "Numerical evaluation of chapman's grazing incidence integral ch(X,x)"
  ! Xp is supposed to be R/H
  ! JMB Update: 05/2017.  Corrected a small error in the y-calc for 
  ! erfc(y)
  !

  a = 1.06069630
  b = 0.55643831
  c = 1.06198960
  d = 1.72456090
  f = 0.56498823
  g = 0.06651874

  do iSpecies = 1, nSpecies

     MaxCh(iSpecies) = 0.0

     do iLon = 1, nLons
        do iLat = 1, nLats

           do iAlt = nAlts, 1, -1

              if (iAlt < nAlts) then
                 Integrals(iLon,iLat,iAlt,iSpecies) = &
                      Integrals(iLon,iLat,iAlt+1,iSpecies) + &
                      NDensityS(iLon,iLat,iAlt,iSpecies,iBlock)&
                      *dAlt_GB(iLon,iLat,iAlt,iBlock)
              else
                 ScaleHeightS = &
                      Temperature(iLon,iLat,iAlt,iBlock) &
                      * TempUnit(iLon,iLat,iAlt) * Boltzmanns_Constant &
                      / (-Gravity_GB(iLon,iLat,iAlt,iBlock) * Mass(iSpecies))
                 Integrals(iLon,iLat,iAlt,iSpecies) = &
                      NDensityS(iLon,iLat,iAlt,iSpecies, iBlock) * ScaleHeightS
              endif
              ScaleHeightS = &
                   Temperature(iLon,iLat,iAlt,iBlock) &
                   *TempUnit(iLon,iLat,iAlt) * Boltzmanns_Constant &
                   / (-Gravity_GB(iLon,iLat,iAlt,iBlock) * Mass(iSpecies))

!              xp(iLon,iLat,iAlt,iSpecies) = &
!                   sqrt(0.5 * pi * &
!                   RadialDistance_GB(iLon, iLat, iAlt, iBlock) / ScaleHeightS)
!              y = xp(iLon,iLat,iAlt,iSpecies) * abs(cos(SZA(iLon,iLat,iBlock)))

              ! JMB Update
              xp(iLon,iLat,iAlt,iSpecies) = &
                   RadialDistance_GB(iLon, iLat, iAlt, iBlock) / ScaleHeightS
              ! Eqn (10) Smith & Smith
              y = sqrt(0.5*xp(iLon,iLat,iAlt,iSpecies)) * &
                  abs(cos(SZA(iLon,iLat,iBlock)))
!              SpeciesChapY(iLon,iLat,iAlt,iSpecies) = y

              ! Eqn (12) Smith and Smith
              if (y < 8) then
                 erfcy(iLon,iLat,iAlt,iSpecies) = (a+b*y) / (c+d*y+y**2)
              else
                 erfcy(iLon,iLat,iAlt,iSpecies) = f / (g + y)
              endif
           enddo

        enddo
     enddo

     do iLon = 1, nLons
        do iLat = 1, nLats
           do iAlt = 1, nAlts

              if (SZA(iLon,iLat,iBlock)<TwoPi/4 .or. &
                   SZA(iLon,iLat,iBlock)>3*TwoPi/4) then


!                 chapman(iLon,iLat,iAlt,iSpecies,iBlock) = &
!                      Integrals(iLon,iLat,iAlt,iSpecies) * &
!                      xp(iLon,iLat,iAlt,iSpecies) * &
!                      erfcy(iLon,iLat,iAlt,iSpecies)

                 ! JMB: 05/04/2017
                 ! Eqn (13) Smith & Smith
                 chapman(iLon,iLat,iAlt,iSpecies,iBlock) = &
                      Integrals(iLon,iLat,iAlt,iSpecies) * &
                      sqrt(0.5*pi*xp(iLon,iLat,iAlt,iSpecies)) * &
                      erfcy(iLon,iLat,iAlt,iSpecies)

                 if (chapman(iLon,iLat,iAlt,iSpecies,iBlock) > &
                      MaxCh(iSpecies)) then
                    MaxCh(iSpecies) = chapman(iLon,iLat,iAlt,iSpecies,iBlock)
                 endif

              else

                ! Determine if the grazing angle (X = 90 deg) lies
                ! within our modeling domain.
                ! This y is a projection of the r-vector at point P
                ! onto the terminator/grazing angle
                 y = RadialDistance_GB(iLon, iLat, iAlt, iBlock) &
                      *abs(cos(SZA(iLon,iLat,iBlock)-pi/2))

                 if (y > RadialDistance_GB(iLon, iLat, 2, iBlock)) then

                    iiAlt = iAlt
                    do while (RadialDistance_GB(iLon, iLat, iiAlt-1, iBlock)>y)
                       iiAlt = iiAlt - 1
                    enddo

                    ! This is approximate.  We should do linear interpolation
                    ! here, but the integral is going to be huge anyways...

                    ScaleHeightS = &
                         Temperature(iLon,iLat,iiAlt,iBlock) &
                         *TempUnit(iLon,iLat,iiAlt) * Boltzmanns_Constant &
                         / (-Gravity_GB(iLon,iLat,iiAlt,iBlock)*Mass(iSpecies))

                    ! xp is evaluated at the tangent point, which is where
                    ! at a radial distance of r*cos(sza)


!                    chapman(iLon,iLat,iAlt,iSpecies,iBlock) = &
!                         xp(iLon,iLat,iiAlt,iSpecies) * &
!                         ScaleHeightS * ( &
!                         2 * NDensityS(iLon,iLat,iiAlt,iSpecies,iBlock) - &
!                         NDensityS(iLon,iLat,iAlt,iSpecies,iBlock) * &
!                         erfcy(iLon,iLat,iAlt,iSpecies))

                    ! JMB Update.  Equation (19) Smith & Smith
                    chapman(iLon,iLat,iAlt,iSpecies,iBlock) = &
                         sqrt(0.5*pi*xp(iLon,iLat,iiAlt,iSpecies)) * &
                         ScaleHeightS * ( &
                         2 * NDensityS(iLon,iLat,iiAlt,iSpecies,iBlock) - &
                         NDensityS(iLon,iLat,iAlt,iSpecies,iBlock) * &
                         erfcy(iLon,iLat,iAlt,iSpecies))



                    if (chapman(iLon,iLat,iAlt,iSpecies,iBlock) > &
                         MaxCh(iSpecies)) then
                       MaxCh(iSpecies)=chapman(iLon,iLat,iAlt,iSpecies,iBlock)
                    endif

                 ! you're within the "shadow" and so make champan huge
                 else

                    chapman(iLon, iLat, iAlt, iSpecies, iBlock) = 1.0e26

                 endif

              endif

           enddo
        enddo
     enddo

     if (MaxCh(iSpecies) == 0.0) &
          MaxCh(iSpecies) = maxval(chapman(:,:,:, iSpecies,iBlock))
     if (iDebugLevel > 3) &
          write(*,*) "====> MaxChapman : ", iSpecies, MaxCh(iSpecies)

  enddo

  call end_timing("chapman_integrals")

end subroutine chapman_integrals
