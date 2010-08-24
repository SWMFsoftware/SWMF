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

  real,parameter, dimension(0:52):: xlat_data=(/&
       11.812,13.777,15.742,17.705,19.665,21.622,23.576,25.527,27.473, &
       29.414,31.350,33.279,35.200,37.112,39.012,40.897,42.763,44.604, &
       46.409,48.163,49.837,51.382,52.725,53.823,54.720,55.488,56.175, &
       56.812,57.413,57.990,58.547,59.090,59.622,60.144,60.659,61.168, &
       61.671,62.170,62.666,63.159,63.649,64.137,64.624,65.109,65.593, &
       66.077,66.560,67.043,67.526,68.009,68.492,68.975,69.458/)


  real :: xlatr(np), xmlt(nt), dlat(np1), energy(neng), sinAo(npit)

end Module ModCrcmGrid


