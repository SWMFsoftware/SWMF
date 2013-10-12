!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!******************************************************************************
! Use Newton's method to find chemical equilibrium for lower boundary
!******************************************************************************

subroutine calc_chemical_equilibrium(DensityHp, DensityH3p)
  use ModCommonVariables
  implicit none
  real, intent(out) :: DensityHp, DensityH3p
  real :: OldDensityHp, OldDensityH3p, NewDensityHp, NewDensityH3p
  real :: DensityH2, DensityH, DensityH2O, DensityCH4,Temperature 
  real :: alpha1, alpha2, beta1, beta2, beta3
  real :: eta0, eta1, eta2, eta3
  real :: zeta0,zeta1,zeta2
  real :: InitialGuess,kr1T,kr2T
  integer :: i
  
  call MODATM(1400.0e5,DensityH2, DensityH, DensityH2O, DensityCH4,Temperature)
  

  kr1T=kr1*(Temperature**(-0.7))
  kr2T=kr2*(Temperature**(-0.5))


  alpha1= jp1*DensityH2 + jp3*DensityH + jp4*DensityH2O
  alpha2= -kc2*DensityH2**2. - kc3*DensityCH4 - kc8*DensityH2O &
          - kc9(1)*DensityH2
  
  beta1 = kc1*(jp2/kc1)*DensityH2
  beta2 = kc2*DensityH2**2.+ kc9(1)*DensityH2
  beta3 =-kc6*DensityCH4-kc7*DensityH2O
  
  eta0  =-kr2T*alpha1**2.
  eta1  =beta3*kr1T*alpha1 - 2.*kr2T*alpha1*alpha2
  eta2  =beta1*kr1T**2. + beta3*kr1T*alpha2 + kr2T*kr1T*alpha1 - kr2T*alpha2**2.
  eta3  =beta2*kr1T**2.-beta3*kr1T**2.+kr2T*kr1T*alpha2 

  
  InitialGuess=1500.
  
  OldDensityHp=InitialGuess
  
  do i=1,50
     NewDensityHp=OldDensityHp &
          -(eta3*OldDensityHp**3.+ eta2*OldDensityHp**2. &
          + eta1*OldDensityHp + eta0) &
          / (3.*eta3*OldDensityHp**2.+ 2.*eta2*OldDensityHp &
          + eta1)
     OldDensityHp=NewDensityHp
!     write(*,*)NewDensityHp
  enddo
  
  DensityHp=NewDensityHp

!  DensityH3p=(kr1T*DensityHp**2.- alpha1 - alpha2*DensityHp)/(-kr1T*DensityHp)

  zeta0=beta1+beta2*DensityHp
  zeta1=beta3-kr2T*DensityHp
  zeta2=-kr2T

  DensityH3p=(-zeta1-(zeta1**2.0-4.0*zeta2*zeta0)**0.5)/(2.0*zeta2)

  return
end subroutine calc_chemical_equilibrium
