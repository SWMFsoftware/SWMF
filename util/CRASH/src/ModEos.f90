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

  use ModStatSum
  use ModAtomicMass

  implicit none

  include 'CRASH_definitions.h'

contains
  !============================================================================
  subroutine eos(uDensityTotal, Rho, iMaterial, &
       TeOut, pTotalOut, GammaOut, Energy0Out, IsError)

    ! Equation of state

    real,    intent(in):: uDensityTotal ! total energy density [J/m^3]
    real,    intent(in):: Rho           ! mass density [kg/m^3] 
    integer, intent(in):: iMaterial     ! sort of material

    real,    optional, intent(out) :: TeOut         ! temperature SI[K]
    real,    optional, intent(out) :: pTotalOut     ! pressure, SI [Pa]
    real,    optional, intent(out) :: GammaOut      ! polytropic index
    real,    optional, intent(out) :: Energy0Out    ! (E-P/(\gamma-1))/\rho
    logical, optional, intent(out) :: IsError

    real   :: NAtomic, UPerAtom
    logical:: IsDegenerated

    character (len=*), parameter:: NameSub='ModEos::eos'
    !----------------------------------------------------------------------!
    if(iMaterial == 2) then

       if( present(GammaOut) .or. present(Energy0Out) ) &
            call CON_stop(NameSub// &
            ': thermodynamic derivatives for polyimide '// &
            'have not been implemented yet')

       call eos_polyimide(UDensityTotal, Rho, TeOut, PTotalOut, &
            IsError=IsError)

       RETURN
    end if

    call set_element(nZ_I(iMaterial))

    ! Get the atomic concentration
    NAtomic=Rho/(cAtomicMass*cAtomicMass_I(nZ_I(iMaterial)))

    ! Get an energy per the atomic cell, express in eV
    UPerAtom = UDensityTotal / (cEV * NAtomic)

    ! Find temperature from dentity and internal energy
    call set_temperature(uPerAtom, Natomic, IsDegenerated)

    if(present(IsError)) IsError = IsDegenerated
    if(IsDegenerated)then
       if(present(IsError))RETURN
       write(*,*) NameSub,' iMaterial =', iMaterial
       write(*,*) NameSub,' UDensityTotal, UPerAtom=', UPerAtom, UDensityTotal
       write(*,*) NameSub,' Rho, NAtomic =', Rho, NAtomic 
       call CON_stop(NameSub//': no EOS for Fermi degenerated state')
    end if

    if(present(TeOut))     TeOut     = Te*cEVToK
    if(present(PTotalOut)) pTotalOut = pressure()
    if(present(GammaOut))then
       call get_thermodyn_derivatives(GammaOut=GammaOut)
       if(present(Energy0Out))&
            Energy0Out = (UDensityTotal - pressure()/(GammaOut-1.0))/Rho
    end if

  end subroutine eos

  !=======================================================================!

  subroutine pressure_to_eint(pTotal, Rho, iMaterial, &
       TeOut, GammaOut, IsError)

    ! Calculate internal energy from pressure

    real,              intent(in) :: PTotal           ! pressure  [J/m^3=Pa]
    real,              intent(in) :: Rho              ! mass density [kg/m^3] 
    integer,           intent(in) :: iMaterial        ! sort of material
    real,    optional, intent(out):: TeOut            ! electron temp. [K]
    real,    optional, intent(out):: GammaOut         ! polytropic index
    logical, optional, intent(out):: IsError          ! true if error occured

    real    :: Natomic, pToNaRatio
    logical :: IsDegenerated

    character (len=*), parameter:: NameSub='ModEos::pressure_to_eint'
    !------------------------------------------------------------------------
    call set_element(nZ_I(iMaterial))

    ! Get the atomic concentration
    Natomic = Rho / (cAtomicMass*cAtomicMass_I(nZ_I(iMaterial)))

    ! Divide pressure by Na , express in eV
    pToNaRatio = pTotal / (cEV * Natomic)

    !Find temperature from dentity and internal energy
    call pressure_to_temperature(pToNaRatio, Natomic, IsDegenerated)

    if(present(IsError)) IsError = IsDegenerated
    if(IsDegenerated)then
       if(present(IsError))RETURN
       write(*,*) NameSub,' pTotal, Rho, iMaterial =', pTotal, Rho, iMaterial
       write(*,*) NameSub,' pToNaRatio, Natomic =', pToNaRatio, Natomic
       call CON_stop(NameSub//': no EOS for Fermi degenerated state')
    end if

    if(present(TeOut))    TeOut = Te*cEvToK

    if(present(GammaOut)) call get_thermodyn_derivatives(GammaOut=GammaOut)

    ! ???
    !if(present(UDensityTotalOut)) &
    !    UDensityTotalOut = NAtomic*cEV*internal_energy()

  end subroutine pressure_to_eint

  !====================================================================

  subroutine eos_polyimide(UDensityTotal, Rho, &
       TeOut, PTotalOut, GammaOut, Energy0Out, IsError)

    ! Equation of state for polyimide

    use ModStatSumMix

    real, intent(in):: UDensityTotal ! total energy density [J/m^3]
    real, intent(in):: Rho           ! mass density [kg/m^3] 

    real,    optional, intent(out) :: TeOut         ! temperature SI[K]
    real,    optional, intent(out) :: PTotalOut     ! pressure, SI [Pa]
    real,    optional, intent(out) :: GammaOut      ! polytropic index
    real,    optional, intent(out) :: Energy0Out    ! (E-P/(\gamma-1))/\rho
    logical, optional, intent(out) :: IsError

    integer,parameter :: nPolyimide = 4

    real, parameter :: cPolyimide_I(nPolyimide) = &
         (/22.0, 10.0, 2.0, 5.0/)/(22.0 + 10.0 + 2.0 +5.0)

    integer, parameter :: nZPolyimide_I(nPolyimide) = (/6, 1, 7, 8/)

    real, parameter :: cAPolyimide = &
         (cAtomicMass_I(6) * 22.0 + &
         cAtomicMass_I(1) * 10.0 + &
         cAtomicMass_I(7) *  2.0 + &
         cAtomicMass_I(8) *  5.0) / &
         (22.0 + 10.0 + 2.0 +5.0)

    real   :: NAtomic, UPerAtom
    logical:: IsDegenerated

    character(len=*), parameter:: NameSub = "ModEos::eos_polyimide"
    !----------------------------------------------------------------------!
    call set_mixture(nPolyimide, nZPolyimide_I, CPolyimide_I)

    !Get the atomic concentration
    NAtomic = Rho / ( cAtomicMass * cAPolyimide )

    !Get an energy per the atomic cell, express in eV
    UPerAtom = UDensityTotal / ( cEV * NAtomic )

    !Find temperature from dentity and internal energy
    call set_temperature_in_mix( UPerAtom, NAtomic, IsDegenerated )

    if(present(IsError)) IsError = IsDegenerated
    if(IsDegenerated)then
       if(present(IsError))RETURN
       write(*,*) NameSub,' UDensityTotal, UPerAtom=', UPerAtom, UDensityTotal
       write(*,*) NameSub,' Rho, NAtomic =', Rho, NAtomic 
       call CON_stop(NameSub//': no EOS for Fermi degenerated state')
    end if

    if(present(TeOut))TeOut = TeMix*cEVToK
    if(present(PTotalOut))PTotalOut = pressure_mix()

  end subroutine eos_polyimide

end module ModEos
