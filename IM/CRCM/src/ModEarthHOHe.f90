Module ModCrcmPlanet
  implicit none

  real :: re_m, dipmom
  real :: Hiono = 120 ! ionosphere altitude in km
  
  !Species Information: Number of Ion Species (Earth=H+,O+,He+)
  integer,parameter :: nspec=4  ! number of species

  ! a0,a1,... are coef. of a polynomial which defines the 
  ! exponent of the charge exchange cross section
  real,parameter,dimension(nspec-1) :: a0_I = (/-18.767,-18.987,-20.789/)
  real,parameter,dimension(nspec-1) :: a1_I = (/-0.11017,-0.10613,0.92316/)
  real,parameter,dimension(nspec-1) :: a2_I = (/-3.8173e-2,-5.4841e-3,-0.68017/)
  real,parameter,dimension(nspec-1) :: a3_I = (/-0.1232,-1.6262e-2,0.66153/)
  real,parameter,dimension(nspec-1) :: a4_I = (/-5.0488e-2,-7.0554e-3,-0.20998/)

  !set species Mass in amu
  real,parameter,dimension(nspec)   :: amu_I = (/1.0,16.0,4.0,5.4462e-4/)


  !set density and temp factor
  real, dimension(nspec) :: dFactor_I =(/0.84,0.15,0.01,1.0/)
  real, dimension(nspec) :: tFactor_I =(/1.0,1.0,1.0,0.128205/)

  !set plot parameters
  character(len=300), parameter :: & 
       NamePlotVar='x y P[nP] HpP[nP] OpP[nP] HepP[nP] eP[nP] ' &
       //'Ppar[nP] HpPpar[nP] OpPpar[nP] HepPpar[nP] ePpar[nP] '&
       //'Phot[nP] HpPhot[nP] OpPhot[nP] HepPhot[nP] ePhot[nP] ' &
       //'Pparhot[nP] HpPparhot[nP] OpPparhot[nP] HepPparhot[nP] ePparhot[nP] ' &
       //'N[/m3] HpN[/m3] OpN[/m3] HepN[/m3] eN[/m3] ' &
       //'Beq[T] Vol[m3/Wb] Pot[Volts] FAC[Amp/m2] g rbody'

  integer, dimension(nspec+1) :: iPplot_I    =(/1,2,3,4,5/)
  integer, dimension(nspec+1) :: iPparplot_I =(/6,7,8,9,10/) 
  integer, dimension(nspec+1) :: iPhotplot_I =(/11,12,13,14,15/) 
  integer, dimension(nspec+1) :: iPparhotplot_I =(/16,17,18,19,20/)
  integer, dimension(nspec+1) :: iNplot_I    =(/21,22,23,24,25/)
  integer, parameter          :: Beq_=26,Vol_=27,Pot_=28, FAC_=29, nVar=29
end Module ModCrcmPlanet
