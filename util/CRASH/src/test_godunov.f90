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
  integer,parameter::nStep=500,nX=1000,nDim=1,nVar=3,iMaterial=0
  real::Cons_VC(nVar,1:nX),Prime_VG(nVar,0:nX+1),Flux_VF(nVar,1:nX+1)
  real::CMax_F(1:nX+1),Gamma_G(0:nX+1)
  real,parameter::cLength=5.0e-3    !5 mm
  real::Dt,Time
  real,parameter::DX=cLength/nX
  integer::iStep,iX
  real,parameter::RhoInit=100.0 !Approximately 30 times the normal density
  real,parameter::UInit = 1.50*RhoInit*(8.31e+3/122.2)*3.0e3 !3000 K 
  real,parameter::UPiston = 3.0e+4 !30 km/s
  real,parameter::CFL=0.9           !
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
  write(24,'(a)') 'X [m] Rho [kg/m^3]   u [m/s]  P [Pa]  Gamma'
  do iStep=1,nStep
     write(*,*)'iStep,StartTime=',iStep,Time
     !Get primitives:
     do iX=1,nX
        Prime_VG(1,iX)=Cons_VC(1,iX)
        Prime_VG(2,iX)=Cons_VC(2,iX)/Prime_VG(1,iX)
        Prime_VG(3,iX)=Cons_VC(3,iX)-0.50*Prime_VG(1,iX)*Prime_VG(2,iX)**2
        call eos(UDensityTotal=Prime_VG(3,iX),& !Input energy density[J/m^3]
               Rho=Prime_VG(1,iX),            & !Input mass density,[kg/m^3] 
               iMaterial=iMaterial,           & !Input: sort of material
               PTotalOut=Prime_VG(3,iX),      & !Output,pressure [Pa]
               GammaOut=Gamma_G(iX))            !Output,polytropic index
     end do
     if(iStep==nStep)then
        do iX=1,nX
           write(24,'(5E13.6)')DX*iX,Prime_VG(:,iX),Gamma_G(iX)
        end do
        exit
     end if
     !Fix ends
     Prime_VG(:,0)   =Prime_VG(:, 1); Gamma_G(0)   =Gamma_G(1)
     Prime_VG(:,nX+1)=Prime_VG(:,nX); Gamma_G(nX+1)=Gamma_G(nX)
     
     !Reflection at the left end
     Prime_VG(2,0) = -Prime_VG(2,0)

     
     !Get fluxes
     do iX=1,nX+1
        call get_godunov_flux(1,(/1.0/),3,Prime_VG(:,iX-1),Prime_VG(:,iX),&
             Flux_VF(:,iX),CMax_F(iX),Gamma_G(iX-1),Gamma_G(iX))
!        if(iX==1)write(*,*)Flux_VF(:,1)
     end do
     Cons_VC(:,1:nX)=Cons_VC(:,1:nX)+(Dt/DX)*(Flux_VF(:,1:nX)-Flux_VF(:,2:nX+1))
     Time = Time + Dt
     call set_dt
  end do
  close(24)
contains
  subroutine set_dt
    Dt=CFL*DX/maxval(CMax_F(1:nX+1))
  end subroutine set_dt
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


