
Module ModConductance

  integer :: longmx, latmx, ndx
  real    :: steplat

  integer, parameter :: latmdx = 40
  integer, parameter :: lonmdx = 30
  integer, parameter :: mndx   = 10
  integer, parameter :: nConductanceSolutions = 4

  real, dimension(mndx) :: halmin, pedmin, avk50, efx50

  real, dimension(0:lonmdx,0:latmdx,mndx) :: &
       halar, pedar, avkar, efxar

  real, dimension(nConductanceSolutions) :: ConductanceBackground

  real, dimension(4) :: bkgc, bkgcer
  real               :: flxmax, dbymax, dbzmax

  integer, parameter :: pedersen_ = 1
  integer, parameter :: hall_     = 2
  integer, parameter :: eflux_    = 3
  integer, parameter :: avee_     = 4

end Module ModConductance
