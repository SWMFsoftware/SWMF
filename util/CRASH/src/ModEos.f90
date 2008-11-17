!^CFG COPYRIGHT UM
module ModEos

  ! Equation Of State (EOS)
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
      
  use ModPolyimide
  use ModStatSum
  use ModAtomicMass

  implicit none
  PRIVATE !Except
  include 'CRASH_definitions.h'
  
  public::eos, UsePreviousTe, Xe_,Be_,Plastic_

  !Eos function for CRASH
  !Usage : call eos(iMaterial,Rho,(TeIn=...| ETotalIn=...| PTotalIn=...),Parameter1Out=,Parameter2Out=...)
  !Or:     call eos(RhoToARatio_I,(TeIn=...| ETotalIn=...| PTotalIn=...),Parameter1Out=,Parameter2Out=...)
  !where the parameters iMaterial and Rho [kg/m^3] are the input parameters for CRASH material iMaterial
  !of the density of Rho; 
  !Or RhoToARatio_I is the array of atomic concentrations of the materials in their mixture
  !
  !Energetic input parameter is the (electron) temperature in K, or the internal energy density [J/m^3] 
  !(the total of electron and ion energies, assuming their temperatures to be equal), or the total 
  !presure [Pa]. In more than one input energetic parameter is given, the first one in the sequence of 
  !(Te,ETotal,PTotal) is used. If none is provided, the code stops 
  !
  !Output parameters: TeOut,  PTotal,ETotal,Gamma( Gamma*PTotal/Rho being the adiabatic compressibility,
  !or, the same, speed of sound squared) 

  interface eos
     module procedure eos_material
     module procedure eos_mixture
  end interface
contains
  !============================================================================

  subroutine eos_material(iMaterial,Rho,&
       TeIn,ETotalIn,pTotalIn,TeOut,ETotalOut,pTotalOut,GammaOut,IsError)
    
    integer, intent(in):: iMaterial     ! sort of material
    real,    intent(in):: Rho           ! mass density [kg/m^3]

    real,    optional, intent(in)  :: TeIn          ! temperature SI[K]
    real,    optional, intent(in)  :: ETotalIn      ! internal energy density [J/m^3]
    real,    optional, intent(in)  :: pTotalIn      ! pressure, SI [Pa]

    real,    optional, intent(out) :: TeOut         ! temperature SI[K]
    real,    optional, intent(out) :: pTotalOut     ! pressure, SI [Pa]
    real,    optional, intent(out) :: ETotalOut     ! internal energy density [J/m^3]
    real,    optional, intent(out) :: GammaOut      ! polytropic index
    logical, optional, intent(out) :: IsError
    
    real   :: NAtomic
    !-------------------------------
    if(iMaterial==2)then
       call set_mixture(nPolyimide, nZPolyimide_I, CPolyimide_I)

       !Get the atomic concentration
       NAtomic = Rho / ( cAtomicMass * cAPolyimide )
    else
       call set_element(nZ_I(iMaterial))

       ! Get the atomic concentration
       NAtomic=Rho/(cAtomicMass*cAtomicMass_I(nZ_I(iMaterial)))
    end if
    call eos_generic(NAtomic,&
         TeIn,ETotalIn,pTotalIn,TeOut,ETotalOut,pTotalOut,GammaOut,IsError)
  end subroutine eos_material
  !============================================================================


   subroutine eos_mixture(RhoToARatio_I,&
       TeIn,ETotalIn,pTotalIn,TeOut,ETotalOut,pTotalOut,GammaOut,IsError)

     real,intent(in) :: RhoToARatio_I( 0:2 )  !Input mass density, SI [kg/m^3] divided by A
     
     real,    optional, intent(in)  :: TeIn          ! temperature SI[K]
     real,    optional, intent(in)  :: ETotalIn      ! internal energy density [J/m^3]
     real,    optional, intent(in)  :: pTotalIn      ! pressure, SI [Pa]

     real,    optional, intent(out) :: TeOut         ! temperature SI[K]
     real,    optional, intent(out) :: pTotalOut     ! pressure, SI [Pa]
     real,    optional, intent(out) :: ETotalOut     ! internal energy density [J/m^3]
     real,    optional, intent(out) :: GammaOut      ! polytropic index
     logical, optional, intent(out) :: IsError
    

     real :: RhoToATotal, NAtomic
    
     integer, parameter :: nAll = 1 + 1 + nPolyimide   
 
     integer, parameter :: nZAll_I(nAll) = (/54 , &  !Xe
                                              4 , &  !Be
                                              6 , &  !C
                                              1 , &  !H
                                              7 , &  !N
                                              8 /)   !O
    real :: ConcentrationAll_I(nAll)
    !-------------------------------!
    RhoToATotal = sum( RhoToARatio_I ) 
    
    !Relative atomic concentrations of Xe, Be and polyimide:
    ConcentrationAll_I( 1:3 ) = RhoToARatio_I / RhoToATotal
    
    !Specify concentrations for C, H, N, O

    ConcentrationAll_I( 3:6 ) = ConcentrationAll_I( 3 ) * CPolyimide_I

    call set_mixture(nAll, nZAll_I, ConcentrationAll_I)

    !Get the atomic concentration
    NAtomic = RhoToATotal / cAtomicMass 

    call eos_generic(NAtomic,&
         TeIn,ETotalIn,pTotalIn,TeOut,ETotalOut,pTotalOut,GammaOut,IsError)
  end subroutine eos_mixture
  !=============================================================================

  subroutine  eos_generic(NAtomic,&
       TeIn,ETotalIn,pTotalIn,TeOut,ETotalOut,pTotalOut,GammaOut,IsError)
    
    
    real,               intent(in) :: NAtomic       ! Atomic concentration,[1/m^3]
     
    real,    optional, intent(in)  :: TeIn          ! temperature SI[K]
    real,    optional, intent(in)  :: ETotalIn      ! internal energy density [J/m^3]
    real,    optional, intent(in)  :: pTotalIn      ! pressure, SI [Pa]
    
    real,    optional, intent(out) :: TeOut         ! temperature SI[K]
    real,    optional, intent(out) :: pTotalOut     ! pressure, SI [Pa]
    real,    optional, intent(out) :: ETotalOut     ! internal energy density [J/m^3]
    real,    optional, intent(out) :: GammaOut      ! polytropic index
    logical, optional, intent(out) :: IsError
    
    real :: EPerAtom, pPerAtom,TeInEV  !Both in eV

    character (len=*), parameter:: NameSub='ModEos::eos'
    !----------------------------------------------------------------------!
    if(present(TeIn))then

       TeInEV = TeIn * cKToEV
       call set_ionization_equilibrium(TeInEV, NAtomic, IsError )

    elseif(present(ETotalIn))then

       ! Get an energy per the atomic cell, express in eV
       EPerAtom = ETotalIn/ (cEV * NAtomic)

       ! Find temperature from dentity and internal energy
       call set_temperature(EPerAtom, Natomic, IsError)
    
    elseif(present(pTotalIn))then
       ! Divide pressure by Na , express in eV
       pPerAtom = pTotalIn / (cEV * Natomic)

       !Find temperature from dentity and pressure
       call pressure_to_temperature(pPerAtom, Natomic, IsError)
    else
       call CON_stop(NameSub//'None of Te|ETotal|pTotal is among input parameters')
    end if

  

    if(present(TeOut))     TeOut     = Te*cEVToK
    if(present(ETotalOut)) ETotalOut = NAtomic*cEV*internal_energy()
    if(present(PTotalOut)) pTotalOut = pressure()
    if(present(GammaOut))  call get_gamma(GammaOut=GammaOut)
  !     if(present(Energy0Out))&
  !          Energy0Out = (UDensityTotal - pressure()/(GammaOut-1.0))/Rho

  !  end if

  end subroutine eos_generic
end module ModEos
