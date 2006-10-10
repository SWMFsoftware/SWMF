!******************************************************************************
! Use Newton's method to find chemical equilibrium for lower boundary
!******************************************************************************

subroutine calc_chemical_equilibrium(DensityHp, DensityH3p)
  implicit none
  real, intent(out) :: DensityHp, DensityH3p
  real :: OldDensityHp, OldDensityH3p, NewDensityHp, NewDensityH3p
  real :: DensityH2, DensityH, DensityH2O, DensityCH4,Temperature 
  real :: alpha1, alpha2, beta1, beta2, beta3
  real :: eta0, eta1, eta2, eta3
  real :: zeta0,zeta1,zeta2
  real :: jp1, jp2, jp3, jp4, kc1, kc2, kc3, kc6, kc7, kc8, kr1, kr2
  real :: InitialGuess
  integer :: i
  
  call MODATM(1400.0e5,DensityH2, DensityH, DensityH2O, DensityCH4,Temperature)
  

  jp1=9.5e-11
  jp2=5.4e-10
  jp3=7.3e-10
  jp4=1.3e-10
  kc1=2.0e-9
  kc2=3.2e-29
  kc3=4.15e-9
  kc6=2.4e-9
  kc7=5.3e-9
  kc8=8.2e-9
  kr1=2.0e-12
  kr2=4.6e-6*Temperature**(-0.65)
  

  alpha1= jp1*DensityH2 + jp3*DensityH + jp4*DensityH2O
  alpha2= -kc2*DensityH2**2. - kc3*DensityCH4 - kc8*DensityH2O
  
  beta1 = kc1*(jp2/kc1)*DensityH2
  beta2 =kc2*DensityH2**2.
  beta3 =-kc6*DensityCH4-kc7*DensityH2O
  
  eta0  =-kr2*alpha1**2.
  eta1  =beta3*kr1*alpha1 - 2.*kr2*alpha1*alpha2
  eta2  =beta1*kr1**2. + beta3*kr1*alpha2 + kr2*kr1*alpha1 - kr2*alpha2**2.
  eta3  =beta2*kr1**2.-beta3*kr1**2.+kr2*kr1*alpha2 

  
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

!  DensityH3p=(kr1*DensityHp**2.- alpha1 - alpha2*DensityHp)/(-kr1*DensityHp)

  zeta0=beta1+beta2*DensityHp
  zeta1=beta3-kr2*DensityHp
  zeta2=-kr2

  DensityH3p=(-zeta1-(zeta1**2.0-4.0*zeta2*zeta0)**0.5)/(2.0*zeta2)

  return
end subroutine calc_chemical_equilibrium
