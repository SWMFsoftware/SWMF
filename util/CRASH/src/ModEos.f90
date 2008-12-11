!^CFG COPYRIGHT UM
module ModEos

  ! Equation Of State (EOS) for ionized plasma
  !
  ! Thermodynamic variables and other notations
  !
  !        Rho - the mass density
  !        E - internal energy of the unit of mass
  !        e, i - electron, ion
  !        V, vol - volume or volumetric
  !        \left(\frac{\partial...}{\partial...}\right)_V - thermodynamical
  !             derivative at constant volume
  !        Te, Ti - electron and ion temperature
  !        iMaterial - integer index of the material:
  !        iMaterial=0 - xenon
  !        iMaterial=1 - beryllium   
  !        iMaterial=2 - plastic
  !        iMaterial=90 - plasma with eos E=aT^4/4; p=aT^4/12; C_V=aT^3 
  !
  ! In the initial CRASH treatment of materials,
  !
  ! 1. We can treat any given material as having a single average Z.
  !
  ! 2. Mixtures can be treated either as an average material or as mixtures.
  !
  ! 3. If mixtures are treated as mixtures, then collisional rates should be 
  !    calculated using the "effective Z", which is the average of Z squared 
  !    divided by the average Z.
  !
  ! 4. The average Z can be determined from the Saha equation, but must not 
  !    exceed the nuclear charge of the material in question.
  !
  ! 5. We do not need to account for electron degeneracy in the initial model.
  !
  ! 6. In our regime of interest, the electrons behave as an ideal gas in an 
  !    ion-sphere environment within which Coulomb interactions do affect the 
  !    electron pressure and internal energy. The electron pressure and 
  !    internal energy are best calculated using equations 3.47 through 3.50 
  !    in R. P. Drake, High Energy Density Physics
  !
  ! 7. The ion pressure is the ideal gas pressure. The ion internal energy 
  !    includes the particle energy of random motion and the energy of 
  !    ionization. The model of eqs. 3.74 through 3.76 in the mentioned book 
  !    is acceptable. Alternatively, a more complex model using actual 
  !    ionization energies would be acceptable.
  !
  ! 8. The materials that matter are
  !    - Beryllium
  !    - Xenon
  !    - Polyimide (C_22 H_10 N_2 O_5)
  !

  use ModPolyimide
  use ModStatSum
  use ModAtomicMass
  use ModESimT4
  use CRASH_ModFermiGas, ONLY: UseFermiGas, LogGeMinBoltzmann

  implicit none

  private !Except

  integer, public, parameter:: Xe_=0      ! Xenon
  integer, public, parameter:: Be_=1      ! Beryllium
  integer, public, parameter:: Plastic_=2 ! Polyimide (C_22 H_10 N_2 O_5)

  public:: UsePreviousTe ! inherited from ModStatSum
  public:: eos, read_eos_parameters
  interface eos
     module procedure eos_material
     module procedure eos_mixture
  end interface

  ! Local variables

  integer, parameter:: Test_=90 ! test material with powerlaw EOS function
  integer, parameter :: nZ_I(Xe_:Plastic_)=(/54, 4, 6 /)

  ! integer, parameter:: C_=3, H_=4, N_=5, O_=6 !Composition elements

contains

  !============================================================================
  subroutine read_eos_parameters

    ! The usage in the user_read_param:
    ! #EOS
    ! T                     UseFermiGas
    ! 4.0                   LogGeMinBoltzmann)

    use ModReadParam, ONLY: read_var
    !-----------------------------------------------------------------------
    ! For now. But it should/could read other things
    call read_var('UseFermiGas',       UseFermiGas)
    call read_var('LogGeMinBoltzmann', LogGeMinBoltzmann)

  end subroutine read_eos_parameters
  !============================================================================

  subroutine eos_material(iMaterial,Rho,&
       TeIn,eTotalIn,pTotalIn,&
       TeOut,eTotalOut,pTotalOut,GammaOut,CvTotalOut,IsError)

    ! Eos function for single material

    integer, intent(in):: iMaterial     ! index of material
    real,    intent(in):: Rho           ! mass density [kg/m^3]

    ! One of the following three energetic input parameters must be present
    real,    optional, intent(in)  :: TeIn       ! temperature SI[K]
    real,    optional, intent(in)  :: eTotalIn   ! internal energy density
    real,    optional, intent(in)  :: pTotalIn   ! pressure

    ! One or more of the output parameters can be present
    real,    optional, intent(out) :: TeOut      ! temperature
    real,    optional, intent(out) :: pTotalOut  ! pressure
    real,    optional, intent(out) :: eTotalOut  ! internal energy density
    real,    optional, intent(out) :: GammaOut   ! polytropic index
    real,    optional, intent(out) :: CvTotalOut ! specific heat / unit volume
    logical, optional, intent(out) :: IsError

    real   :: Natomic
    !-------------------------------------------------------------------------
    if(iMaterial == Test_)then
       call eos_esimt4(TeIn, eTotalIn, pTotalIn, &
            TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut)
       RETURN
    elseif(iMaterial == Plastic_)then
       call set_mixture(nPolyimide, nZPolyimide_I, CPolyimide_I)

       !Get the atomic concentration
       Natomic = Rho / ( cAtomicMass * cAPolyimide )
    else
       call set_element(nZ_I(iMaterial))

       ! Get the atomic concentration
       Natomic=Rho/(cAtomicMass*cAtomicMass_I(nZ_I(iMaterial)))
    end if
    call eos_generic(Natomic, TeIn, eTotalIn, pTotalIn, &
         TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut,IsError)

  end subroutine eos_material

  !============================================================================

  subroutine eos_mixture(RhoToARatio_I,&
       TeIn, eTotalIn, pTotalIn, &
       TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut, IsError)

    ! Eos function for mixed material
    real, intent(in) :: RhoToARatio_I(Xe_:Plastic_) ! Mass densities/A

    ! One of the following three energetic input parameters must be present
    real,    optional, intent(in)  :: TeIn          ! temperature
    real,    optional, intent(in)  :: eTotalIn      ! internal energy density
    real,    optional, intent(in)  :: pTotalIn      ! pressure

    ! One or more of the output parameters can be present
    real,    optional, intent(out) :: TeOut         ! temperature
    real,    optional, intent(out) :: pTotalOut     ! pressure
    real,    optional, intent(out) :: eTotalOut     ! internal energy density 
    real,    optional, intent(out) :: GammaOut      ! polytropic index
    real,    optional, intent(out) :: CvTotalOut    ! specific heat per volume
    logical, optional, intent(out) :: IsError       ! error flag

    real :: RhoToATotal, Natomic

    integer, parameter :: nAll = 1 + 1 + nPolyimide   

    integer, parameter :: nZAll_I(nAll) = (/54 , &  !Xe
         4 , &  !Be
         6 , &  !C
         1 , &  !H
         7 , &  !N
         8 /)   !O
    real :: ConcentrationAll_I(nAll)
    !-------------------------------------------------------------------------
    RhoToATotal = sum( RhoToARatio_I ) 

    !Relative atomic concentrations of Xe, Be and polyimide:
    ConcentrationAll_I( 1:3 ) = RhoToARatio_I / RhoToATotal

    !Specify concentrations for C, H, N, O

    ConcentrationAll_I( 3:6 ) = ConcentrationAll_I( 3 ) * CPolyimide_I

    call set_mixture(nAll, nZAll_I, ConcentrationAll_I)

    !Get the atomic concentration
    Natomic = RhoToATotal / cAtomicMass 

    call eos_generic(Natomic, TeIn, eTotalIn, pTotalIn, &
         TeOut, eTotalOut, PtotalOut, GammaOut, CvTotalOut, IsError)

  end subroutine eos_mixture

  !============================================================================

  subroutine eos_generic(Natomic, TeIn, eTotalIn, pTotalIn,&
       TeOut, eTotalOut, pTotalOut, GammaOut, CvTotalOut, IsError)

    real,              intent(in)  :: Natomic    ! Atomic concentration
    real,    optional, intent(in)  :: TeIn       ! temperature
    real,    optional, intent(in)  :: eTotalIn   ! internal energy density
    real,    optional, intent(in)  :: PtotalIn   ! pressure

    real,    optional, intent(out) :: TeOut      ! temperature
    real,    optional, intent(out) :: pTotalOut  ! pressure
    real,    optional, intent(out) :: eTotalOut  ! internal energy density 
    real,    optional, intent(out) :: GammaOut   ! polytropic index
    real,    optional, intent(out) :: CvTotalOut ! Specific heat per volume
    logical, optional, intent(out) :: IsError

    real :: ePerAtom, pPerAtom,TeInEV  !Both in eV

    character (len=*), parameter:: NameSub='ModEos::eos'
    !----------------------------------------------------------------------!
    if(present(TeIn))then

       TeInEV = TeIn * cKToEV
       call set_ionization_equilibrium(TeInEV, Natomic, IsError )

    elseif(present(eTotalIn))then

       ! Get an energy per the atomic cell, express in eV
       ePerAtom = eTotalIn/ (cEV * Natomic)

       ! Find temperature from dentity and internal energy
       call set_temperature(ePerAtom, Natomic, IsError)

    elseif(present(pTotalIn))then
       ! Divide pressure by Na , express in eV
       pPerAtom = pTotalIn / (cEV * Natomic)

       !Find temperature from dentity and pressure
       call pressure_to_temperature(pPerAtom, Natomic, IsError)
    else
       call CON_stop(NameSub// &
            'None of Te, eTotal, or pTotal is among the input parameters')
    end if

    if(present(TeOut))      TeOut     = Te*cEVToK
    if(present(eTotalOut))  eTotalOut = Natomic*cEV*internal_energy()
    if(present(PTotalOut))  pTotalOut = pressure()
    if(present(GammaOut))   call get_gamma(GammaOut=GammaOut)
    if(present(CvTotalOut)) CvTotalOut = (Natomic*cBoltzmann)*heat_capacity()

  end subroutine eos_generic

end module ModEos
