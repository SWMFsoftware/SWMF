!^CFG COPYRIGHT UM
module CRASH_ModTransport
  use CRASH_ModStatSumMix,ONLY: &
       Z2PerA, Z2, ZAv, CoulombLog, Na, Te
  !Calculation of the transport coefficients 
  !in dense plasmas: electron heat conduction,
  !Electron-ion temperature relaxation etc.

  implicit none
contains
  !=================================================================
  real function nu_ei()
    !The conventional auxiliary characterstic of the
    !electron-ion collisions
    use ModConst,ONLY: cSigmaThomson, cElectronMass, &
                       cLightSpeed, cEV
    use ModNumConst,ONLY: cPi
   
    !The rest mass energy
    real, parameter:: MeC2InEV = &
         cElectronMass * cLightSpeed * cLightSpeed /cEV
    real :: TePerMeC2
    !-------------------------------
    !The expression to eveluate:
    !\nu_{\rm eff} = \frac{4}{3}\sqrt{\frac{2\pi}{m_e}}*
    !                *\frac{e^4*Na*CoulombLog}{Te^{3/2}}
    
    !We re-rite it in the following manner:
    ! \nu_{\rm eff} = 
    !  =\frac{8\pi}{3}\left(\frac{e^2}{m_ec^2}\right)^2 * 
    ! (which is the Thomson cross-section)
    ! ...*(m_ec^2 [eV]/T_e[eV])^2 *\sqrt{T_e/(2\pi m_e}
    !The last term equals c\sqrt{T_e[eV]/{2\pi* m_ec^2[eV]}
    !With this account:
    
    TePerMeC2 = Te/MeC2InEV !Dimensionless

    nu_ei = cSigmaThomson * cLightSpeed * Na * CoulombLog / &
         ( TePerMeC2 * sqrt(2.0 * cPi * TePerMeC2 ) ) ! s^{-1}
    
  end function nu_ei
  !================
  real function te_ti_relaxation()
    use ModConst,ONLY: cElectronMass, cAtomicMass 
    
    !The Landau formula, may be written only for T_i in our case:
    !d(3/2)T_i/dt =  TeTiRelaxation * (T_e - T_i)
    !so that the relaxation rate for the ion ENERGY is expressed in terms
    !of the ion and electron TEMPERATURES.
    !Therefore, TeTiRelaxation = (3/2) * (2m_e/M_A)*Z*<Z^2/A>*\nu_{ei}:

    real,parameter :: cMassRatio = 3.0 * cElectronMass/ cAtomicMass
    !-
    te_ti_relaxation =  cMassRatio * ZAv * Z2PerA *  nu_ei() !s^{-1}

    !The practical application:
    !The ralaxation per the unit of volume:
    !dE_i/dt = TeTiRelaxation * Na * (T_e-T_i)* 
    !* {1 cEV or k_{Bolttmann}}, 
    !depending on the unit of temperature.
    !The dimensionality of dE_i/dt as calculated in this way is [J/(m^3s)]
    !The energy transfer to electrons: dE_e/dt = - dE_i/dt
    !
    !The equations for temperature:
    !
    !dT_i/dt = (2/3) TeTiRelaxation * (T_e - T_i), however
    !dT_e/dt = - TeTiRelaxation * (T_e - T_i) * (k_B*N_a)/(C_{Ve}
    !
    !Only in a fully ionized plasma of a very high temperature the 
    !latter is equal to TeTiRelaxation*(T_i-T_e)*2/(3Z), which is the
    !expression used by Landau.
  end function te_ti_relaxation
  !===========================
  real function electron_heat_conductivity()
    use ModConst, ONLY: cBoltzmann, cEV, cElectronMass
    !----------------------------------
    !The dependence on Z as well as the numerical coefficient
    !are to be modified:
    
    electron_heat_conductivity = cBoltzmann * (Te * cEV/cElectronMass) * &
         (Na * ZAv) /(Z2 * nu_ei())  ![J/(M*s*K)]
  end function electron_heat_conductivity

end module CRASH_ModTransport
