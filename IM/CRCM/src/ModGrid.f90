Module ModCrcmGrid
  implicit none
  
  ! define dimensions of CRCM grids
  integer,parameter :: np1=51,nt1=48,neng1=12,npit1=12!,nspec1=1  
  integer,parameter :: nm=35,nk=28 ! dimension of CRCM magnetic moment and K

  integer,parameter :: np=51    ! dimension of the CRCM latitude grid
  integer,parameter :: nt=48    ! dimension of the CRCM local-time grid
  integer,parameter :: neng=12  ! dimension of the CRCM energy grid
  integer,parameter :: npit=12  ! dimension of the CRCM pitch-angle grid

  ! These have to be initialized so that IM_set_grid does not fail on non-IM PEs
  real:: xlat(np) = 0.0, phi(nt1)=0.0

  real :: xlatr(np), xmlt(nt), dlat(np1), energy(neng), sinAo(npit)

end Module ModCrcmGrid


