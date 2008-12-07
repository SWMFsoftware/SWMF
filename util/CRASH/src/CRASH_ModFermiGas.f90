!^CFG COPYRIGHT UM

module ModFermiGas
  implicit none
  SAVE
  PRIVATE !Except

  ! Logarithm of the electron statistical  
  ! weight=exp(-\mu/T)
  real, public:: LogGe 

  ! At LogGe >= LogGeMinBoltzmann the electrons are treated as the  
  !Boltzmann gas
  real, public:: LogGeMinBoltzmann = 2.0

  ! At LogGe >= LogGeMinFermi the effects of Fermi statistics are  
  ! accounted for
  real, public, parameter:: LogGeMinFermi = 0.0

  !The module does nothing if this logical is false  
  logical,public :: UseFermiGas = .false.
  
  !Correcting coefficients in thermodynamic functions  
  real,public :: rMinus = 1.0, rPlus = 1.0  

  integer,parameter :: nStep = 20
  integer,parameter :: NuEqMinus12_ = 1, NuEq12_ = 2, NuEq32_ = 3
  
  real :: FermiFunctionTable_II(0:nStep, NuEqMinus12_ : NuEq32_)
  
  public:: init_fermi_function, iterate_ge 
  public:: read_fermi_gas_param
  
  ! public :: test_fermi_function !Uncomment for testing

contains
  !==============================
  !The usage in the user_read_param:
  !case('#FERMIGASEFFECT')
  !call read_fermi_gas_param
  !In PARAM.in:
  !#FERMIGASEFFECT
  !T
  !4.0
  subroutine read_fermi_gas_param
    use ModReadParam
    !----------------
    call read_var('UseFermiGas',UseFermiGas)
    call read_var('LogGeMinBoltzmann',LogGeMinBoltzmann)
  end subroutine read_fermi_gas_param
  !=====================================
  !The Fermi functions are defined to be
  !Fe_{\nu}(g_e)=(1/\Gamma(\nu+1))\int_0^\infty{x^\nu dx/(g_e \exp(x)+1)}
  !At abs(g_e)>=1 the integral can be taken by developing into a convergent 
  !series in powers of 1/g_e. At g_e=1 the integral may be also expressed 
  !in terms of the Riemann zeta-function, which is used in testing
  !==============================
  !Initialaze the calculation of Fermi functions for \nu+1=1/2, 3/2, 5/2
  subroutine init_fermi_function
    integer:: iStep,iNu
    real :: Ge, GeN, SumTerm, RealN
    real, parameter :: MuPlus1_I(3) = (/ 0.5, 1.50, 2.50 /)
    !------------------------------------------------------
    !Fill in the lookup table
    do iStep=0,nStep-1
       Ge = -exp(&
           -(LogGeMinFermi*(nStep-iStep)+iStep*LogGeMinBoltzmann)/nStep  )
       do iNu = NuEqMinus12_ ,  NuEq32_
          GeN = 1.0; SumTerm = 1.0 ; RealN = 1.0
          FermiFunctionTable_II(iStep, iNu) = 1.0 
          do while(abs(SumTerm) > 3.0e-4)
             RealN = RealN + 1.0
             GeN = GeN * Ge
             SumTerm = GeN/RealN** MuPlus1_I(iNu)
             FermiFunctionTable_II(iStep, iNu) = &
                  FermiFunctionTable_II(iStep, iNu) + SumTerm
          end do
       end do
    end do 
   
    !To achieve a continuity at G_e_Min_Boltzmann, put the functions  to be 1 
    FermiFunctionTable_II(nStep, :) = 1.0

  end subroutine init_fermi_function
  !========================================================
  subroutine iterate_ge(zAvr,Delta2I,LogGe1,Diff)
    real, intent(in) :: zAvr    !<i>
    real, intent(in) :: Delta2I !<i^2>-<i>^2
    real, intent(in) :: LogGe1  !log(G_{e1}

    real, intent(out):: Diff    !<i>-G_{e1}Fe_{1/2}(g_e)

    integer :: iStep, iStep1
    real :: Residual
    real::FermiFunction_I(NuEqMinus12_:NuEq32_)
    !------------------------


    Residual = (LogGe - LogGeMinFermi)/&
         (LogGeMinBoltzmann - LogGeMinFermi)*nStep

    Residual = min(max(Residual,0.0),real(nStep))
    iStep    = int(Residual)
    iStep1   = min(iStep+1,nStep)
   
    Residual = Residual - real(iStep)
    FermiFunction_I = FermiFunctionTable_II(iStep,:)*(1.0-Residual)+&
         FermiFunctionTable_II(iStep1,:)*Residual
    


    rMinus = FermiFunction_I(NuEqMinus12_)/FermiFunction_I(NuEq12_)
    rPlus  = FermiFunction_I(NuEq32_) / FermiFunction_I(NuEq12_)

    Diff = (zAvr - exp(LogGe1-LogGe) * FermiFunction_I(NuEq12_))/&
         (Delta2I/zAvr + rMinus)
    LogGe = LogGe - Diff/zAvr

  end subroutine iterate_ge
  !========================================================
  subroutine test_fermi_function
    integer :: iStep
    !---------------
    LogGeMinBoltzmann = 2.0*log(10.)

    write(*,*)'zeta-function multiplied by (1-2*0.5^\nu) = ',&
         -1.46035*(1.0 -2.0*sqrt(0.5)),&
         2.612*(1.0-sqrt(0.5)),1.341*(1.0-0.5*sqrt(0.5))
    call init_fermi_function
    write(*,*)'The computed values at g_e=1 are:', &
         FermiFunctionTable_II(0 , :)

    write(*,*)'log10g_e  g_e*Fe1/2   rMinus rPlus'
    do iStep = 0, nStep
       write(*,*) (LogGeMinFermi*(nStep-iStep)+iStep*LogGeMinBoltzmann)/nStep/&
            log(10.0), FermiFunctionTable_II(iStep,2),&
            FermiFunctionTable_II(iStep,1)/FermiFunctionTable_II(iStep,2),&
            FermiFunctionTable_II(iStep,3)/FermiFunctionTable_II(iStep,2)
    end do
  end subroutine test_fermi_function
  
end module ModFermiGas
