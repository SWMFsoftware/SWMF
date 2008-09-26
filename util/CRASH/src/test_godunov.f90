!^CFG COPYRIGHT UM
!=====================Test for the Godunov scheme with the eos==========!!!!
!
! Thermodynamical variables and other notations
!         \rho, Rho - the mass density
!         {\cal E}, E - internal energy of the unit of mass
!         e, i - electron, ion
!        V, vol - volume or volumetric
!        \left(\frac{\partial...}{\partial...}\right)_V - thermodynamical
!             derivative at constant volume
!        T_{e,i}, Te,Ti - electron and ion temperature
!        iMaterial - integer variable, a signature of the material:
!        iMaterial=0 - xenon
!        iMaterial=1 - beryllium  
program test_Godunov
  use ModExactRS
  use ModAtomicMass
  use ModEos
  implicit none
  integer,parameter::nX=1000,nDim=1,nVar=3,iMaterial=0
  real::Cons_VC(nVar,1:nX),Prime_VG(nVar,0:nX+1),Flux_VF(nVar,1:nX+1)
  real::CMax_F(1:nX+1),Gamma_G(0:nX+1),Energy0_G(0:nX+1),InternalEnergy_G(0:nX+1)
  real::GammaMax
  real,parameter::cLength=5.0e-3    !5 mm
  real::Dt,Time
  real,parameter::DX=cLength/nX
  integer::iX
  real,parameter::RhoInit=100.0 !Approximately 30 times the normal density
  real,parameter::UInit = 1.50*RhoInit*(8.31e+3/122.2)*3.0e3 !3000 K 
  real,parameter::UPiston = 3.0e+4 !30 km/s
  real,parameter::CFL=0.8
  real,parameter::TimeOut=5.0e-7 !500 ns
!^CFG COPYRIGHT UM
!=====================Equation Of State (EOS)===========================!!!!
!
! Thermodynamical variables and other notations
!         \rho, Rho - the mass density
!         {\cal E}, E - internal energy of the unit of mass
!         e, i - electron, ion
!        V, vol - volume or volumetric
!        \left(\frac{\partial...}{\partial...}\right)_V - thermodynamical
!             derivative at constant volume
!        T_{e,i}, Te,Ti - electron and ion temperature
!        iMaterial - integer variable, a signature of the material:
!        iMaterial=0 - xenon
!        iMaterial=1 - beryllium         

  Dt=0.0; Time =0.0 
  !Initial condition, in the frame of reference comoving with a piston
  do iX=1,nX
     Cons_VC(1,iX)=RhoInit
     Cons_VC(2,iX)=-UPiston*RhoInit
     Cons_VC(3,iX)= UInit +0.50*RhoInit*UPiston**2
  end do
     
  open(24,file='Godunov_scheme_for_Xe',status='replace')
  write(24,'(a)') 'X [mm] Rho [100 kg/m^3]   u [10,000 m/s]  P [MBar]  Gamma'
  do 
     write(*,*)'iStep,StartTime=',Time
     !Get primitives:
     do iX=1,nX
        Prime_VG(1,iX)=Cons_VC(1,iX)
        Prime_VG(2,iX)=Cons_VC(2,iX)/Prime_VG(1,iX)
        InternalEnergy_G(iX)=Cons_VC(3,iX)-0.50*Prime_VG(1,iX)*Prime_VG(2,iX)**2
        call eos(UDensityTotal=InternalEnergy_G(iX),& !Input energy density[J/m^3]
               Rho=Prime_VG(1,iX),            & !Input mass density,[kg/m^3] 
               iMaterial=iMaterial,           & !Input: sort of material
               PTotalOut=Prime_VG(3,iX),      & !Output,pressure [Pa]
               GammaOut=Gamma_G(iX),          & !Output,polytropic index
               Energy0Out=Energy0_G(iX))     !Output   (E-P/(\gamma-1))/\rho
     end do
     if(Time>TimeOut)then
        do iX=1,nX
           write(24,'(5(E14.6,2X))')DX*(iX-0.5)/1.0e-3,&
           Prime_VG(1,iX)/RhoInit,&
           Prime_VG(2,iX)/1.0e4,&
           Prime_VG(3,iX)/1.0e11,&
           Gamma_G(iX)
        end do
        exit
     end if
     !Fix ends
     Prime_VG(:,0)   =Prime_VG(:, 1); Gamma_G(0)   =Gamma_G(1)
     Prime_VG(:,nX+1)=Prime_VG(:,nX); Gamma_G(nX+1)=Gamma_G(nX)
     Energy0_G(0)   =Energy0_G(1); Energy0_G( nX+1)   =Energy0_G(nX) 
     InternalEnergy_G(0)   =InternalEnergy_G(1)
     InternalEnergy_G( nX+1)   =InternalEnergy_G(nX) 
     !Reflection at the left end
     Prime_VG(2,0) = - Prime_VG(2,0)

     
     !Get fluxes
     do iX=1,nX+1
        GammaMax=5.0/3.0
        call get_godunov_flux(1,(/1.0/),3,Prime_VG(:,iX-1),Prime_VG(:,iX),&
             Flux_VF(:,iX),CMax_F(iX),GammaMax,GammaMax,&
             (InternalEnergy_G(iX-1)-Prime_VG(3,iX-1)/(GammaMax-1.0))/Prime_VG(1,iX-1),&
             (InternalEnergy_G(iX)-Prime_VG(3,iX)/(GammaMax-1.0))/Prime_VG(1,iX))
     end do
     Cons_VC(:,1:nX)=Cons_VC(:,1:nX)+&
          (Dt/DX)*(Flux_VF(:,1:nX)-Flux_VF(:,2:nX+1))
     Dt=min(CFL*DX/maxval(CMax_F(1:nX+1)),1.001*(TimeOut-Time))

     Time = Time + Dt
  end do
  close(24)
end program test_Godunov
!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  implicit none
  character (len=*), intent(in) :: StringError
end subroutine CON_stop
!============================================================================
subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe
end subroutine CON_set_do_test


