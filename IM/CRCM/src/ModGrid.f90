Module ModCrcmGrid
  implicit none
  
  ! define dimensions of CRCM grids
  integer,parameter :: np1=51,nt1=48,neng1=12,npit1=12,nspec1=1  
  integer,parameter :: nm=35,nk=28 ! dimension of CRCM magnetic moment and K

  integer,parameter :: np=51    ! dimension of the CRCM latitude grid
  integer,parameter :: nt=48    ! dimension of the CRCM local-time grid
  integer,parameter :: neng=12  ! dimension of the CRCM energy grid
  integer,parameter :: npit=12  ! dimension of the CRCM pitch-angle grid
  integer,parameter :: nspec=1  ! number of species
  

  real :: xlat(np),xmlt(nt),phi(nt1),dlat(np1),energy(neng),sinAo(npit)
end Module ModCrcmGrid


