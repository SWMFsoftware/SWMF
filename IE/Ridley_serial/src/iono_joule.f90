!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine ionosphere_jouleheating_ionflux(iBlock, ETh, EPs, SigmaP, &
     Joule, IonNumFlux)

  !\
  ! Joule heating is determined by SigmaP * E^2
  !/

  use ModIonosphere
  use CON_world
  use IE_ModMain
  use CON_planet_field
  use ModCoordTransform, ONLY: sph_to_xyz

  implicit none

  integer, intent(in) :: iBlock

  real, dimension(1:IONO_nTheta, 1:IONO_nPsi) :: &
        SigmaP, ETh, EPs, IonNumFlux_tmp
  real, intent(out), dimension(1:IONO_nTheta, 1:IONO_nPsi) ::Joule, IonNumFlux

  integer            :: i, j, iHemisphere
  real, dimension(3) :: B_D, bIono_D, XyzIono_D, Xyz_tmp  
  real               :: bIono, B, ratioOH
  ! FAST observation height = 4000km
  real, parameter    :: height_fast = 4.0e6  

  !\
  ! Joule heating is assumed to be equal to Poynting Flux, 
  ! since according to the Poynting Theorem, divergence of S = J*E,
  ! which is integral(S dA) = integral(J*E dV),
  ! but J in the formula is current density with a unit of A/m^3.
  ! In this routine, J is height-integrated (A/m^2), therefore, 
  ! J * E = S in this code.
  !/

	Joule(:,:) =0.0 
	IonNumFlux(:,:) = 0.0
  ! Joule heating (or Poynting flux) at ionosphere altitude
  Joule = SigmaP * (ETh**2+EPs**2) 

  !\
  ! Ion number flux at ionosphere altitude, based on Strangeway et al(2005):
  ! At an altitude of 4000km(where the FAST orbit), f = 2.142e7 * (S)^1.265, 
  ! where unit of f is 1/(cm^2 * s), and unit of S is mW/m^2
  !/

  ! Mapping to the ionosphere altitude along dipole magnetic field lines
  select case(iBlock)
  case(1)
     ! North
     do i = 1, IONO_nTheta
        do j = 1, IONO_nPsi
           call sph_to_xyz(Radius, IONO_NORTH_Theta(i,j), IONO_NORTH_Psi(i,j), XyzIono_D)
           call get_planet_field(Time_Simulation, XyzIono_D, 'SMG', bIono_D)
           bIono = sqrt(sum(bIono_D**2))

           call map_planet_field(Time_Simulation, XyzIono_D, 'SMG', &
                (height_fast + IONO_Radius),Xyz_tmp, iHemisphere)
           call get_planet_field(Time_Simulation, Xyz_tmp, 'SMG', B_D)
           b = sqrt(sum(B_D**2))
           
           ! Flux at ionosphere altitude by mapping both S and f down to the ionosphere 
           ! from 4000km, caution in the units (flux in /m^2/s, jouleheating in W/m^2)
           ! The total ion number flux = No * Vo
           IonNumFlux(i,j) = 2.142e7 * (Joule(i,j)*1.0e3)**1.265 * (b/bIono)**0.265 * 1.0e4

           ! Next, set up the ratio of Op/Hp as a function of Joule Heating. 
           ! A hypothesis about this is: r = a*log(b*Joule)
           ! A threshold of Joule is needed to have the Op outflow, which is about 2mW/m^2
           ! A high joule heating results in a ratio value approaching to 0.5, which is around 85mW/m^2
           ! b*J0 = 10 when r=0, and J0=2mW/m^2    ==> b = 10/J0
           ! a*log(b*Jm)=0.45=r, when Jm=85mW/m^2  ==> a = 0.45/log(b*Jm)
           
!           if (Joule(i,j)*1.0e3 < 2.0) then 
!              ratioOH = 0 
!           else 
!              ratioOH = 0.45/(log10(1./2.*85.)) * log10(1./2.*Joule(i,j)*1.0e3)
!           endif
!           ! Assume the velocity of Hp is 10 times larger than Op
!           IonNumFlux(i,j,1) = IonNumFlux_tmp(i,j) /(1 + ratioOH * 0.1)
!           IonNumFlux(i,j,2) = IonNumFlux(i,j,1) * ratioOH * 0.1
!
!           call map_planet_field(Time_Simulation, XyzIono_D, 'SMG NORM', &
!                rBody, Xyz_tmp, iHemisphere)
!           call get_planet_field(Time_Simulation, Xyz_tmp, 'SMG NORM', B_D)
!           b_rbody = sqrt(sum(B_D**2))
!
!           ! Flux at rBody (inner boundary of GM)
!           gm_inner_IonNumFlux(i,j) = IonNumFlux(i,j) * b_rbody/bIono
!
        end do
     end do

  case(2)
     ! South
     do i = 1, IONO_nTheta
        do j = 1, IONO_nPsi
           call sph_to_xyz(Radius, IONO_SOUTH_Theta(i,j), IONO_SOUTH_Psi(i,j), XyzIono_D)
           call get_planet_field(Time_Simulation, XyzIono_D, 'SMG', bIono_D)
           bIono = sqrt(sum(bIono_D**2))

           call map_planet_field(Time_Simulation, XyzIono_D, 'SMG', &
                height_fast + IONO_Radius, Xyz_tmp, iHemisphere)
           call get_planet_field(Time_Simulation, Xyz_tmp, 'SMG', B_D)
           b = sqrt(sum(B_D**2))
           
           ! Flux at ionosphere altitude by mapping both S and f down from 4000km
           IonNumFlux(i,j) = 2.142e7 * (Joule(i,j)*1.0e03)**1.265 * (b/bIono)**0.265*1.0e04

!           ratioOH = 0.45/(log10(10./2.)) * log10(10./2.*Joule(i,j)*1.0e3)
!
!           ! Assume the velocity of Hp is 10 times larger than Op
!           IonNumFlux(i,j,1) = IonNumFlux_tmp(i,j) /(1 + ratioOH * 0.1)
!           IonNumFlux(i,j,2) = IonNumFlux(i,j,1) * ratioOH * 0.1
!
!           
!           call map_planet_field(Time_Simulation, XyzIono_D, 'SMG NORM', &
!                rBody, Xyz_tmp, iHemisphere)
!           call get_planet_field(Time_Simulation, Xyz_tmp, 'SMG NORM', B_D)
!           b_rbody = sqrt(sum(B_D**2))
!
!           ! Flux at rBody (inner boundary of GM)
!           gm_inner_IonNumFlux(i,j) = IonNumFlux(i,j) * b_rbody/bIono
!
        end do
     end do

  end select
end subroutine ionosphere_jouleheating_ionflux
