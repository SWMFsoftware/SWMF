!^CFG COPYRIGHT UM
!^CFG FILE PROJECTION

module ModProject

  implicit none
  save
  !\
  ! Parameters for projection scheme:
  !
  ! proj_method determines the iterative scheme to be used
  ! if proj_method='cg'       then Conjugate Gradient method is used
  !                           (this is only good for symmetric matrices!)
  !    proj_method='bicgstab' then BiConjugate Gradient method is used
  !
  ! proj_typestop determines the stopping condition for the iterations.
  ! if proj_typestop='rel' then sum((div B)**2) is reduced by a factor 
  !    proj_divbcoeff
  ! if proj_typestop='max' then |div B| is kept below a limit determeined by 
  !    proj_divbcoeff and
  !    proj_divbconst 
  !    if proj_divbcoeff>0 then 
  !       |div B|< proj_divbcoeff * divbmax_0  
  !          where divbmax_0 is the max(|div B|) after the first time step
  !    otherwise
  !       |div B|< proj_divbconst
  ! 
  ! proj_matvecmax is an upper limit on matrix-vector products for iterations
  !
  ! proj_boundtype determines the boundary conditions for the Poisson equation
  ! if proj_boundtype='zero' : phi=0 in ghost cells outside of comput. domain
  !
  !/

  character (len=10):: proj_method   ='cg        '
  character (len=3)::  proj_typestop ='rel'
  real ::              proj_divbcoeff=0.1
  real ::              proj_divbconst=0.0
  integer ::           proj_matvecmax=50
  character (len=10):: proj_boundtype='zero'

  ! Counter for matrix vector multiplications, and for errors of solver

  integer :: nmatvectotal, poissonerror

  ! Minimum value for divB (a small number)

  real, parameter :: divbmin=1E-10

  ! Maximum value for divB (absolute or relative)
  real :: DivbMax

end module ModProject
