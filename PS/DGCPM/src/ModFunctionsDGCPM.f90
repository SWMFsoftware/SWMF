! File Name: dgcpm_functions.f90
!
! Contains : All major calculation functions originally contained
! within pbo.f. Moved from this file to be accessable by outside 
! f90 routines, and in an attempt to move the code to the f90 standard.
!
! Last Modified : January 2012, Aron Dodger
        Module ModFunctionsDGCPM

        Contains
!cccccccccccccccccccccccccc
!cc function saturation ccc
!cccccccccccccccccccccccccc

      real function saturation(l)

! Carpenter and Anderson's saturation density in units of
! particles / m**3 
! Carpenter and Anderson, JGR, p. 1097, 1992.

! input: l in re
      real l

! output: saturation in particles / m**3

      saturation = (1.0e6) * 10.0**((-0.3145*l)+3.9043)

      return
      end function saturation

!cccccccccccccccccccccc
!cc function trough ccc
!cccccccccccccccccccccc

      real function trough(l)

! Carpenter and Anderson's trough density in units of
! particles / m**3 
! Carpenter and Anderson, JGR, p. 1097, 1992.

! input: l in re
      real l

! output: trough in particles / m**3

      trough = (1.0e6) * 0.5 * ((10.0/l)**4.0)

      return
      end function trough

!ccccccccccccccccccccccccccccccccc
!cc function DipoleFluxTubeVol ccc
!ccccccccccccccccccccccccccccccccc

      real function dipoleFluxTubeVol(l)

! calculates the unit volume of a dipole magnetic field flux tube
! (the volume in m**3 per unit of magnetic flux(weber))

! 1 tesla = newton/(ampere-meter) or (volt-sec)/meter**2
! 1 weber = Tesla-m**2 or joule/ampere or volt-sec

! input: l in re
      real l

! output: dipoleFluxTubevol in m/Tesla or m**3/weber

      real pi, re, mu, m

      pi = 3.14159          ! rad
      re = 6.378e6          ! radius of Earth in meters
      mu = 4.0*pi*1.0e-7    ! newtons/amps**2
      m = 8.05e22           ! amps*meter**2

      dipoleFluxTubeVol = ((4.0*pi)/(mu*m)) * (32.0/35.0) * (l**4) * & 
         sqrt(1.0-(1.0/l)) * (1.0+(1.0/(2.0*l))+(3.0/(8.0*l*l)) + &
         (5.0/(16.0*l*l*l))) * (re**4.0)

      return
      end function dipoleFluxTubeVol

      End Module ModFunctionsDGCPM

