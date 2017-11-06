Module ModSAMI
  ! allocatable total matrices
  
  ! Total 
  !
  !      real denit(nz,nf,nlt,nion)
  
! Total 
  
  real, dimension(:,:,:,:), allocatable :: denit,dennt
  !      real denit(nz,nf,nlt,nion),dennt(nz,nf,nlt,nneut)
  real, dimension(:,:,:,:), allocatable :: vsit, sumvsit
!      real vsit(nz,nf,nlt,nion),sumvsit(nz,nf,nlt,nion)
  real, dimension(:,:,:,:), allocatable :: tit
  !      real tit(nz,nf,nlt,nion)
  
  real, dimension(:,:,:), allocatable :: ut, vt, vpit, net
  !      real ut(nz,nf,nlt),vt(nz,nf,nlt),vpit(nz,nf,nlt),net(nz,nf,nlt)
  real, dimension(:,:,:), allocatable :: tet, tnt
  !      real tet(nz,nf,nlt),tnt(nz,nf,nlt)
  
  !     height integrated pedersen/hall conductivities
  real, dimension(:,:,:), allocatable :: u1t, u2t, u3t, u4t
  !      real u1t(nz,nf,nlt),u2t(nz,nf,nlt),u3t(nz,nf,nlt),u4t(nz,nf,nlt)

  real, dimension(:,:,:), allocatable :: vnqt, vnpt, vnphit
  !      real vnqt(nz,nf,nlt),vnpt(nz,nf,nlt),vnphit(nz,nf,nlt)
  
  real, dimension(:,:,:), allocatable :: jpt, jphit
  !      real jpt(nz,nf,nlt),jphit(nz,nf,nlt)

  real, dimension(:,:,:), allocatable :: u1pt, u2st, u3ht
  
  real, dimension(:,:,:), allocatable :: sigmapict,sigmahict
  
  real, dimension(:,:,:), allocatable :: sigmapt,sigmaht
  
! Output matrices for restart
  
!  real :: deniout(nz,nf,nl,nion,numwork),  &
!          tiout(nz,nf,nl,nion,numwork),    &
!          vsiout(nz,nf,nl,nion,numwork) 
!  real teout(nz,nf,nl,numwork) 
!  real*8 dphi(nnx+1,nnyt)

  real,allocatable:: deniout(:,:,:,:,:),  &
          tiout(:,:,:,:,:),  &
          vsiout(:,:,:,:,:)
  real,allocatable :: teout(:,:,:,:)
  real*8,allocatable ::dphi(:,:)

! make phi arrays allocatable
  real, allocatable ::  phi(:,:),phialt(:,:),philon(:,:),p_crit(:)

! Left and right taskid
  integer :: left,right
  ! mpi info
  integer :: iComm,nProc,iProc
end Module ModSAMI
