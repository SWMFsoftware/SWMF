!^CFG COPYRIGHT UM
module HEDPlasma
  use ModAdvance,Te_GB     =>B0xCell_BLK,& !Electron temperature, K
                 EIon_GB   =>B0yCell_BLK,& !IONIZATION energy per mass
                 Gamma_GB  =>B0zCell_BLK,& !Effective polytropic index
                 Cv_CB     =>NormB0_CB,  & !Specific heat
                 ESource_CB=>DivB0_CB ,  & !Energy source
                 MSource_DCB=>CurlB0_DCB   !Momentum Source
  implicit none
  SAVE
end module HEDPlasma
