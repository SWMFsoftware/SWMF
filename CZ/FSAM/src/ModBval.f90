module ModBval

  implicit none

  private

  public :: bvalv
  public :: bvals
  public :: bvalp
  public :: bvalemf1
  public :: bvalemf2
  public :: bvalemf3
  public :: bvaleta
  public :: bvalb

contains
  !-----------------------------------------------------------------------
  subroutine message_pass_nodes(idtop, idbot, nI1, nJ1, nK1, nI2, nJ2, nK2, &
       sendtop, sendbot, recvtop, recvbot)
    use ModMpi
    implicit none

    integer, intent(in)  :: idtop, idbot, nI1, nJ1, nK1, nI2, nJ2, nK2
    real,    intent(in)  :: sendtop(nI1,nJ1,nK1), sendbot(nI2,nJ2,nK2)
    real,    intent(out) :: recvtop(nI2,nJ2,nK2), recvbot(nI1,nJ1,nK1)

    integer :: ireq(4), istatus_arr(MPI_STATUS_SIZE,4), ierr
    integer :: ncount_top, ncount_bot
    real, allocatable, dimension(:) :: sendtopbuf, sendbotbuf, &
         recvtopbuf, recvbotbuf 
    !-----------------------------------------------------------------------
    ncount_top = nI1*nJ1*nK1
    ncount_bot = nI2*nJ2*nK2
    allocate(sendtopbuf(1:ncount_top), sendbotbuf(1:ncount_bot))
    allocate(recvbotbuf(1:ncount_top), recvtopbuf(1:ncount_bot))
    sendtopbuf = reshape(sendtop,(/ncount_top/))
    sendbotbuf = reshape(sendbot,(/ncount_bot/))

    call MPI_IRECV(recvtopbuf, ncount_bot, MPI_DOUBLE_PRECISION, idtop, 1, &
         MPI_COMM_WORLD, ireq(2), ierr)
    call MPI_IRECV(recvbotbuf, ncount_top, MPI_DOUBLE_PRECISION, idbot, 2, &
         MPI_COMM_WORLD, ireq(4), ierr)
    call MPI_ISEND(sendbotbuf, ncount_bot, MPI_DOUBLE_PRECISION, idbot, 1, &
         MPI_COMM_WORLD, ireq(1), ierr)
    call MPI_ISEND(sendtopbuf, ncount_top, MPI_DOUBLE_PRECISION, idtop, 2, &
         MPI_COMM_WORLD, ireq(3), ierr)
    call MPI_WAITALL(4,ireq,istatus_arr,ierr)

    recvtop = reshape(recvtopbuf, (/nI2, nJ2, nK2/))
    recvbot = reshape(recvbotbuf,(/nI1, nJ1, nK1/))

    deallocate(sendtopbuf, sendbotbuf, recvbotbuf, recvtopbuf)

  end subroutine message_pass_nodes

subroutine bvalv
  use ModPar
  use ModGrid
  use ModField,     ONLY: v1, v2, v3
  use ModBoundary,  ONLY: noib
  implicit none

  integer idtop, idbot
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  !-------------------------------------------------------
  !------------------ i - boundary ----------------------
  idtop = myid1 + 1 + nproc1*myid2
  idbot = myid1 - 1 + nproc1*myid2
  if(myid1 .eq. 0) idbot = nproc1 - 1 + nproc1*myid2
  if(myid1 .eq. nproc1-1) idtop = 0 + nproc1*myid2
  ! v1
  allocate(recvtop(is:isp2,js:je,ks:ke),recvbot(iem1:ie,js:je,ks:ke))
  call message_pass_nodes(idtop, idbot, &
    2, je-js+1, ke-ks+1, 3, je-js+1, ke-ks+1, &
    v1(iem1:ie,js:je,ks:ke), v1(is:isp2,js:je,ks:ke), recvtop, recvbot)
  v1(iep1:iep3,js:je,ks:ke) = recvtop
  v1(ism2:ism1,js:je,ks:ke) = recvbot
  deallocate(recvtop, recvbot)
  ! v2
  allocate(recvtop(is:isp1,js:je,ks:ke),recvbot(iem1:ie,js:je,ks:ke))
  call message_pass_nodes(idtop, idbot, &
    2, je-js+1, ke-ks+1, 2, je-js+1, ke-ks+1, &
    v2(iem1:ie,js:je,ks:ke), v2(is:isp1,js:je,ks:ke), recvtop, recvbot)
  v2(iep1:iep2,js:je,ks:ke) = recvtop
  v2(ism2:ism1,js:je,ks:ke) = recvbot
  deallocate(recvtop, recvbot)
  ! v3
  allocate(recvtop(is:isp1,js:je,ks:ke),recvbot(iem1:ie,js:je,ks:ke))
  call message_pass_nodes(idtop, idbot, &
       2, je-js+1, ke-ks+1, 2, je-js+1, ke-ks+1, &
       v3(iem1:ie,js:je,ks:ke), v3(is:isp1,js:je,ks:ke), recvtop, recvbot)
  v3(iep1:iep2,js:je,ks:ke) = recvtop
  v3(ism2:ism1,js:je,ks:ke) = recvbot
  deallocate(recvtop, recvbot)
  ! Inner i boundary.
  if(myid1 .eq. 0) then
     v1(is  ,js:je,ks:ke) = 0.D0
     v1(ism1,js:je,ks:ke) = - v1(isp1,js:je,ks:ke)
     v1(ism2,js:je,ks:ke) = - v1(isp2,js:je,ks:ke)
     v2(ism1,js:je,ks:ke) = g2bi(is)  /g2bi(ism1)*v2(is  ,js:je,ks:ke)
     v2(ism2,js:je,ks:ke) = g2bi(isp1)/g2bi(ism2)*v2(isp1,js:je,ks:ke)
     v3(ism1,js:je,ks:ke) = g31bi(is)  /g31bi(ism1)*v3(is,  js:je,ks:ke)
     v3(ism2,js:je,ks:ke) = g31bi(isp1)/g31bi(ism2)*v3(isp1,js:je,ks:ke)
  endif
  ! Outer i boundary.
  if(myid1 .eq. nproc1-1) then
     if (noib .eq. 2) then
        v1(iep1,js:je,ks:ke) = v1(ie,  js:je,ks:ke)
        v1(iep2,js:je,ks:ke) = v1(iep1,js:je,ks:ke)
        v1(iep3,js:je,ks:ke) = v1(iep2,js:je,ks:ke)
        v2(iep1,js:je,ks:ke) = v2(ie,  js:je,ks:ke)
        v2(iep2,js:je,ks:ke) = v2(iep1,js:je,ks:ke)
        v3(iep1,js:je,ks:ke) = v3(ie,  js:je,ks:ke)
        v3(iep2,js:je,ks:ke) = v3(iep1,js:je,ks:ke)
     else
        v1(iep1,js:je,ks:ke) = 0.D0
        v1(iep2,js:je,ks:ke) = - v1(ie,  js:je,ks:ke)
        v1(iep3,js:je,ks:ke) = - v1(iem1,js:je,ks:ke)
        v2(iep1,js:je,ks:ke) = g2bi(ie)  /g2bi(iep1)*v2(ie  ,js:je,ks:ke)
        v2(iep2,js:je,ks:ke) = g2bi(iem1)/g2bi(iep2)*v2(iem1,js:je,ks:ke)
        v3(iep1,js:je,ks:ke) = g31bi(ie)  /g31bi(iep1)*v3(ie,  js:je,ks:ke)
        v3(iep2,js:je,ks:ke) = g31bi(iem1)/g31bi(iep2)*v3(iem1,js:je,ks:ke)
     endif
  endif
  !------------------ j - boundary -----------------
  idtop = myid1 + (myid2 + 1)*nproc1
  idbot = myid1 + (myid2 - 1)*nproc1
  if(myid2 .eq. 0) idbot = myid1 + (nproc2 - 1)*nproc1
  if(myid2 .eq. nproc2-1) idtop = myid1 + 0*nproc1
  ! v1
  allocate(recvtop(ism2:iep3,js:jsp1,ks:ke),recvbot(ism2:iep3,jem1:je,ks:ke))
  call message_pass_nodes(idtop, idbot, &
    iep3-ism2+1, 2, ke-ks+1, iep3-ism2+1, 2, ke-ks+1, &
    v1(ism2:iep3,jem1:je,ks:ke), v1(ism2:iep3,js:jsp1,ks:ke), &
    recvtop, recvbot)
  v1(ism2:iep3,jep1:jep2,ks:ke) = recvtop
  v1(ism2:iep3,jsm2:jsm1,ks:ke) = recvbot
  deallocate(recvtop,recvbot)
  ! v2
  allocate(recvtop(ism2:iep2,js:jsp2,ks:ke),recvbot(ism2:iep2,jem1:je,ks:ke))
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, 2, ke-ks+1, iep2-ism2+1, 3, ke-ks+1, &
    v2(ism2:iep2,jem1:je,ks:ke), v2(ism2:iep2,js:jsp2,ks:ke), &
    recvtop, recvbot)
  v2(ism2:iep2,jep1:jep3,ks:ke) = recvtop
  v2(ism2:iep2,jsm2:jsm1,ks:ke) = recvbot
  deallocate(recvtop,recvbot)
  ! v3
  allocate(recvtop(ism2:iep2,js:jsp1,ks:ke),recvbot(ism2:iep2,jem1:je,ks:ke))
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, 2, ke-ks+1, iep2-ism2+1, 2, ke-ks+1, &
    v3(ism2:iep2,jem1:je,ks:ke), v3(ism2:iep2,js:jsp1,ks:ke), &
    recvtop, recvbot)
  v3(ism2:iep2,jep1:jep2,ks:ke) = recvtop
  v3(ism2:iep2,jsm2:jsm1,ks:ke) = recvbot
  deallocate(recvtop,recvbot)
  ! Inner j boundary
  if(myid2 .eq. 0) then
     v1(ism2:iep3,jsm1,ks:ke) = v1(ism2:iep3,js  ,ks:ke)
     v1(ism2:iep3,jsm2,ks:ke) = v1(ism2:iep3,jsp1,ks:ke)
     v2(ism2:iep2,js  ,ks:ke) = 0.D0
     v2(ism2:iep2,jsm1,ks:ke) = - v2(ism2:iep2,jsp1,ks:ke)
     v2(ism2:iep2,jsm2,ks:ke) = - v2(ism2:iep2,jsp2,ks:ke)
     v3(ism2:iep2,jsm1,ks:ke) = g32bi(js)  /g32bi(jsm1)*v3(ism2:iep2,js  ,ks:ke)
     v3(ism2:iep2,jsm2,ks:ke) = g32bi(jsp1)/g32bi(jsm2)*v3(ism2:iep2,jsp1,ks:ke)
  endif
  ! Outer j boundary
  if(myid2 .eq. nproc2-1) then
     v1(ism2:iep3,jep1,ks:ke) = v1(ism2:iep3,je  ,ks:ke)
     v1(ism2:iep3,jep2,ks:ke) = v1(ism2:iep3,jem1,ks:ke)
     v2(ism2:iep2,jep1,ks:ke) = 0.D0
     v2(ism2:iep2,jep2,ks:ke) = - v2(ism2:iep2,je  ,ks:ke)
     v2(ism2:iep2,jep3,ks:ke) = - v2(ism2:iep2,jem1,ks:ke)
     v3(ism2:iep2,jep1,ks:ke) = g32bi(je)  /g32bi(jep1)*v3(ism2:iep2,je  ,ks:ke)
     v3(ism2:iep2,jep2,ks:ke) = g32bi(jem1)/g32bi(jep2)*v3(ism2:iep2,jem1,ks:ke)
  endif

  !------------------ k - boundary -----------------
  ! Inner k boundary
  v1(ism2:iep3,jsm2:jep2,ksm1) = v1(ism2:iep3,jsm2:jep2,ke  )
  v1(ism2:iep3,jsm2:jep2,ksm2) = v1(ism2:iep3,jsm2:jep2,kem1)
  v2(ism2:iep2,jsm2:jep3,ksm1) = v2(ism2:iep2,jsm2:jep3,ke  )
  v2(ism2:iep2,jsm2:jep3,ksm2) = v2(ism2:iep2,jsm2:jep3,kem1)
  v3(ism2:iep2,jsm2:jep2,ksm1) = v3(ism2:iep2,jsm2:jep2,ke  )
  v3(ism2:iep2,jsm2:jep2,ksm2) = v3(ism2:iep2,jsm2:jep2,kem1)
  ! Outer k boundary
  v1(ism2:iep3,jsm2:jep2,kep1) = v1(ism2:iep3,jsm2:jep2,ks  )
  v1(ism2:iep3,jsm2:jep2,kep2) = v1(ism2:iep3,jsm2:jep2,ksp1)
  v2(ism2:iep2,jsm2:jep3,kep1) = v2(ism2:iep2,jsm2:jep3,ks  )
  v2(ism2:iep2,jsm2:jep3,kep2) = v2(ism2:iep2,jsm2:jep3,ksp1)
  v3(ism2:iep2,jsm2:jep2,kep1) = v3(ism2:iep2,jsm2:jep2,ks  )
  v3(ism2:iep2,jsm2:jep2,kep2) = v3(ism2:iep2,jsm2:jep2,ksp1)
  v3(ism2:iep2,jsm2:jep2,kep3) = v3(ism2:iep2,jsm2:jep2,ksp2)

end subroutine bvalv

subroutine bvals
  use ModPar
  use ModGrid
  use ModField,        ONLY: s
  use ModBoundary,     ONLY: niib, noib, siib1
  use ModGetQtys
  implicit none

  integer :: idtop, idbot
  real  :: smeanbot, smeantop
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  !---------------------------------------------------------
  !------------------ i - boundary -----------------
  idtop = myid1 + 1 + myid2*nproc1
  idbot = myid1 - 1 + myid2*nproc1
  if(myid1 .eq. 0) idbot = nproc1 - 1 + nproc1*myid2
  if(myid1 .eq. nproc1-1) idtop = 0 + nproc1*myid2
  allocate(recvtop(iep1:iep2,js:je,ks:ke), recvbot(ism2:ism1,js:je,ks:ke))
  call message_pass_nodes(idtop, idbot, &
    2, je-js+1, ke-ks+1, 2, je-js+1, ke-ks+1, &
    s(iem1:ie,js:je,ks:ke), s(is:isp1,js:je,ks:ke), recvtop, recvbot)
  s(iep1:iep2,js:je,ks:ke) = recvtop
  s(ism2:ism1,js:je,ks:ke) = recvbot
  deallocate(recvtop, recvbot)
  ! Inner i boundary
  if (niib .ne. 3) call get_mean_entropy(smeanbot, smeantop)
  if(myid1 .eq. 0) then
     if(niib .eq. 3) then
        s(ism1,js:je,ks:ke) = 2.D0*siib1(js:je,ks:ke) - s(  is,js:je,ks:ke)
        s(ism2,js:je,ks:ke) = 2.D0*siib1(js:je,ks:ke) - s(isp1,js:je,ks:ke)
     else
        s(  is,js:je,ks:ke) = smeanbot + siib1(js:je,ks:ke)
        s(ism1,js:je,ks:ke) = s(  is,js:je,ks:ke)
        s(ism2,js:je,ks:ke) = s(isp1,js:je,ks:ke)
     endif
  endif
  ! Outer i boundary
  if(myid1 .eq. nproc1-1) then
     if(noib .eq. 2) then
        s(iep1,js:je,ks:ke) = s(  ie,js:je,ks:ke)
        s(iep2,js:je,ks:ke) = s(iep1,js:je,ks:ke)
     elseif(noib .eq. 3 .or. noib .eq. 5) then
        s(iep1,js:je,ks:ke) = - s(  ie,js:je,ks:ke)
        s(iep2,js:je,ks:ke) = - s(iem1,js:je,ks:ke)
     else
        s(iep1,js:je,ks:ke) = s(  ie,js:je,ks:ke)
        s(iep2,js:je,ks:ke) = s(iem1,js:je,ks:ke)
     endif
  endif

  !------------------ j - boundary -----------------
  idtop = myid1 + (myid2 + 1)*nproc1
  idbot = myid1 + (myid2 - 1)*nproc1
  if(myid2 .eq. 0) idbot = myid1 + (nproc2 - 1)*nproc1
  if(myid2 .eq. nproc2-1) idtop = myid1 + 0*nproc1
  allocate(recvtop(ism2:iep2,jep1:jep2,ks:ke),recvbot(ism2:iep2,jsm2:jsm1,ks:ke))
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, 2, ke-ks+1, iep2-ism2+1, 2, ke-ks+1, &
    s(ism2:iep2,jem1:je,ks:ke), s(ism2:iep2,js:jsp1,ks:ke), recvtop, recvbot)
  s(ism2:iep2,jep1:jep2,ks:ke) = recvtop
  s(ism2:iep2,jsm2:jsm1,ks:ke) = recvbot
  deallocate(recvtop, recvbot)
  ! Inner j boundary
  if(myid2 .eq. 0) then
     s(ism2:iep2,jsm1,ks:ke) = s(ism2:iep2,  js,ks:ke)
     s(ism2:iep2,jsm2,ks:ke) = s(ism2:iep2,jsp1,ks:ke)
  endif
  ! Outer j boundary
  if(myid2 .eq. nproc2-1) then
     s(ism2:iep2,jep1,ks:ke) = s(ism2:iep2,  je,ks:ke)
     s(ism2:iep2,jep2,ks:ke) = s(ism2:iep2,jem1,ks:ke)
  endif
  !------------------ k - boundary -----------------
  ! Inner k boundary
  s(ism2:iep2,jsm2:jep2,ksm1) = s(ism2:iep2,jsm2:jep2,  ke)
  s(ism2:iep2,jsm2:jep2,ksm2) = s(ism2:iep2,jsm2:jep2,kem1)
  ! Outer k boundary
  s(ism2:iep2,jsm2:jep2,kep1) = s(ism2:iep2,jsm2:jep2,  ks)
  s(ism2:iep2,jsm2:jep2,kep2) = s(ism2:iep2,jsm2:jep2,ksp1)

end subroutine bvals

subroutine bvaldels
! specify boundary of dels2 such that dels2(ism2:iep2,js:jep1,ks:ke) are
! defined
  use ModPar
  use ModGrid
  use ModDel,         ONLY: dels2
  implicit none

  integer :: idtop, idbot
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  !-------------------------------------------------
  !------------------ i - boundary -----------------
  idtop = myid1 + 1 + myid2*nproc1
  idbot = myid1 - 1 + myid2*nproc1
  if(myid1==0) idbot = nproc1 - 1 + nproc1*myid2
  if(myid1==nproc1-1) idtop = 0 + nproc1*myid2
  allocate(recvtop(iep1:iep2,js:jep1,ks:ke), recvbot(ism2:ism1,js:jep1,ks:ke))
  call message_pass_nodes(idtop, idbot, &
    2, jep1-js+1, ke-ks+1, 2, jep1-js+1, ke-ks+1, &
    dels2(iem1:ie,js:jep1,ks:ke), dels2(is:isp1,js:jep1,ks:ke), &
    recvtop, recvbot)
  dels2(iep1:iep2,js:jep1,ks:ke) = recvtop
  dels2(ism2:ism1,js:jep1,ks:ke) = recvbot
  deallocate(recvtop, recvbot)
  ! Inner i boundary.
  if(myid1==0) then
     dels2(ism1,js:jep1,ks:ke) = dels2(is  ,js:jep1,ks:ke)
     dels2(ism2,js:jep1,ks:ke) = dels2(isp1,js:jep1,ks:ke)
  endif
  ! Outer i boundary.
  if(myid1==nproc1-1) then
     dels2(iep1,js:jep1,ks:ke) = dels2(ie  ,js:jep1,ks:ke)
     dels2(iep2,js:jep1,ks:ke) = dels2(iem1,js:jep1,ks:ke)
  endif
  !------------------ j - boundary -----------------
  ! Outer j boundary
  if((myid2==nproc2-1) .and. &
    (abs(x2a(jep1)-pi/2.D0)<1.D-6 .or. abs(x2a(jep1)-pi)<1.D-6)) &
    dels2(ism2:iep2,jep1,ks:ke) = 0.D0
  ! Inner j boundary
  if((myid2==0) .and. (abs(x2a(js))<1.D-6)) &
    dels2(ism2:iep2,js,ks:ke) = 0.D0

end subroutine bvaldels

subroutine bvalp
  use ModPar
  use ModGrid
  use ModField,     ONLY: p
  use ModDel
  use ModBoundary
  use ModBack,     ONLY: fact
  implicit none

  integer :: i, k
  integer :: idtop, idbot
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  !-------------------------------------------------
  !------------------ i - boundary -----------------
  idtop = myid1 + 1 + myid2*nproc1
  idbot = myid1 - 1 + myid2*nproc1
  if(myid1==0) idbot = nproc1 - 1 + nproc1*myid2
  if(myid1==nproc1-1) idtop = 0 + nproc1*myid2
  allocate(recvtop(iep1:iep2,js:je,ks:ke), recvbot(ism2:ism1,js:je,ks:ke))
  call message_pass_nodes(idtop, idbot, &
    2, je-js+1, ke-ks+1, 2, je-js+1, ke-ks+1, &
       p(iem1:ie,js:je,ks:ke), p(is:isp1,js:je,ks:ke), recvtop, recvbot)
  p(iep1:iep2,js:je,ks:ke) = recvtop
  p(ism2:ism1,js:je,ks:ke) = recvbot
  deallocate(recvtop, recvbot)
  ! Inner i boundary
  if(myid1==0) then
     p(ism1,js:je,ks:ke) = p(is  ,js:je,ks:ke)*prbt &
      + dels1(is,js:je,ks:ke)*srbt
     p(ism2,js:je,ks:ke) = p(ism1,js:je,ks:ke)
  endif
  ! Outer i boundary
  if(myid1==nproc1-1) then
     p(iep1,js:je,ks:ke) = p(ie  ,js:je,ks:ke)*prtp &
      + dels1(iep1,js:je,ks:ke)*srtp
     p(iep2,js:je,ks:ke) = p(iep1,js:je,ks:ke)
  endif
  !------------------ j - boundary -----------------
  idtop = myid1 + (myid2 + 1)*nproc1
  idbot = myid1 + (myid2 - 1)*nproc1
  if(myid2==0) idbot = myid1 + (nproc2 - 1)*nproc1
  if(myid2==nproc2-1) idtop = myid1 + 0*nproc1
  allocate(recvtop(ism2:iep2,jep1:jep2,ks:ke),recvbot(ism2:iep2,jsm2:jsm1,ks:ke))
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, 2, ke-ks+1, iep2-ism2+1, 2, ke-ks+1, &
    p(ism2:iep2,jem1:je,ks:ke), p(ism2:iep2,js:jsp1,ks:ke), recvtop, recvbot)
  p(ism2:iep2,jep1:jep2,ks:ke) = recvtop
  p(ism2:iep2,jsm2:jsm1,ks:ke) = recvbot
  deallocate(recvtop, recvbot)
  ! Inner j boundary
  call bvaldels
  if(myid2==0) then
       do i=ism2,iep2
         do k=ks,ke
             p(i,jsm1,k) = p(i,js  ,k) &
               -fact(i+myid1*(in-5))*g2b(i)*dx2b(js) &
               *dels2(i,js,k)
             p(i,jsm2,k) = p(i,jsm1,k)
         enddo
       enddo
  endif
  ! Outer j boundary
  if(myid2==nproc2-1) then
       do i=ism2,iep2
         do k=ks,ke
             p(i,jep1,k) = p(i,je  ,k) &
               +fact(i+myid1*(in-5))*g2b(i)*dx2b(jep1) &
               *dels2(i,jep1,k)
             p(i,jep2,k) = p(i,jep1,k)
         enddo
       enddo
  endif
  !------------------ k - boundary -----------------
  ! Inner k boundary
  p(ism2:iep2,jsm2:jep2,ksm1) = p(ism2:iep2,jsm2:jep2,  ke)
  p(ism2:iep2,jsm2:jep2,ksm2) = p(ism2:iep2,jsm2:jep2,kem1)
  ! Outer k boundary
  p(ism2:iep2,jsm2:jep2,kep1) = p(ism2:iep2,jsm2:jep2,  ks)
  p(ism2:iep2,jsm2:jep2,kep2) = p(ism2:iep2,jsm2:jep2,ksp1)

end subroutine bvalp

subroutine bvalemf1(emf1)
  use ModPar,      ONLY: myid1, myid2, nproc1, nproc2, in, jn, kn
  use ModGrid
  use ModBoundary,  ONLY: noib
  implicit none

  real, intent(inout), dimension(in,jn,kn) :: emf1

  integer idtop, idbot
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  !-------------------------------------------------------
  !------------------ i - boundary ----------------------
  idtop = myid1 + 1 + nproc1*myid2
  idbot = myid1 - 1 + nproc1*myid2
  if(myid1 .eq. 0) idbot = nproc1 - 1 + nproc1*myid2
  if(myid1 .eq. nproc1-1) idtop = 0 + nproc1*myid2
  allocate(recvtop(iep1:iep2,js:jep1,ks:kep1), &
           recvbot(ism2:ism1,js:jep1,ks:kep1))
  call message_pass_nodes(idtop, idbot, &
    2, jep1-js+1, kep1-ks+1, 2, jep1-js+1, kep1-ks+1, &
    emf1(iem1:ie,js:jep1,ks:kep1), emf1(is:isp1,js:jep1,ks:kep1), &
    recvtop, recvbot)
  emf1(iep1:iep2,js:jep1,ks:kep1) = recvtop
  emf1(ism2:ism1,js:jep1,ks:kep1) = recvbot
  deallocate(recvtop, recvbot)
  ! Inner i boundary.
  if(myid1 .eq. 0) then
     emf1(ism1,js:jep1,ks:kep1) = emf1(is  ,js:jep1,ks:kep1)
     emf1(ism2,js:jep1,ks:kep1) = emf1(isp1,js:jep1,ks:kep1)
  endif
  ! Outer i boundary.
  if(myid1 .eq. nproc1-1) then
     if(noib .eq. 2) then
        emf1(iep1,js:jep1,ks:kep1) = emf1(ie  ,js:jep1,ks:kep1)
        emf1(iep2,js:jep1,ks:kep1) = emf1(ie  ,js:jep1,ks:kep1)
     elseif (noib .eq. 5) then
        emf1(iep1,js:jep1,ks:kep1) = - emf1(ie  ,js:jep1,ks:kep1)
        emf1(iep2,js:jep1,ks:kep1) = - emf1(iem1,js:jep1,ks:kep1)
     else
        emf1(iep1,js:jep1,ks:kep1) = emf1(ie  ,js:jep1,ks:kep1)
        emf1(iep2,js:jep1,ks:kep1) = emf1(iem1,js:jep1,ks:kep1)
     endif
  endif
  !------------------ j - boundary -----------------
  idtop = myid1 + (myid2 + 1)*nproc1
  idbot = myid1 + (myid2 - 1)*nproc1
  if(myid2 .eq. 0) idbot = myid1 + (nproc2 - 1)*nproc1
  if(myid2 .eq. nproc2-1) idtop = myid1 + 0*nproc1
  allocate(recvtop(ism2:iep2,jep2:jep3,ks:kep1), &
           recvbot(ism2:iep2,jsm2:jsm1,ks:kep1))
  call message_pass_nodes(idtop, idbot, &
    iep2-ism2+1, 2, kep1-ks+1, iep2-ism2+1, 2, kep1-ks+1, &
    emf1(ism2:iep2,jem1:je,ks:kep1),emf1(ism2:iep2,jsp1:jsp2,ks:kep1), &
       recvtop,recvbot)
  emf1(ism2:iep2,jep2:jep3,ks:kep1) = recvtop
  emf1(ism2:iep2,jsm2:jsm1,ks:kep1) = recvbot
  deallocate(recvtop,recvbot)
  ! Inner j boundary.
  if(myid2 .eq. 0) then
     emf1(ism2:iep2,js  ,ks:kep1) = 0.D0
     emf1(ism2:iep2,jsm1,ks:kep1) = - emf1(ism2:iep2,jsp1,ks:kep1)
     emf1(ism2:iep2,jsm2,ks:kep1) = - emf1(ism2:iep2,jsp2,ks:kep1)
  endif
  ! Outer j boundary.
  if(myid2 .eq. nproc2-1) then
     emf1(ism2:iep2,jep1,ks:kep1) = 0.D0
     emf1(ism2:iep2,jep2,ks:kep1) = - emf1(ism2:iep2,je  ,ks:kep1)
     emf1(ism2:iep2,jep3,ks:kep1) = - emf1(ism2:iep2,jem1,ks:kep1)
  endif
  !------------------ k - boundary -----------------
  ! Inner k boundary
  emf1(ism2:iep2,jsm2:jep3,ksm1) = emf1(ism2:iep2,jsm2:jep3,ke  )
  emf1(ism2:iep2,jsm2:jep3,ksm2) = emf1(ism2:iep2,jsm2:jep3,kem1)
  ! Outer k boundary.
  emf1(ism2:iep2,jsm2:jep3,kep2) = emf1(ism2:iep2,jsm2:jep3,ksp1)
  emf1(ism2:iep2,jsm2:jep3,kep3) = emf1(ism2:iep2,jsm2:jep3,ksp2)

  return
end subroutine bvalemf1

subroutine bvalemf2(emf2)
  use ModPar,      ONLY: myid1, myid2, nproc1, nproc2, in, jn, kn
  use ModGrid
  use ModBoundary,  ONLY: noib
  implicit none

  real, intent(inout), dimension(in,jn,kn) :: emf2

  integer idtop, idbot
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  !-----------------------------------------------------------------------------
  !------------------ j - boundary ----------------- 
  idtop = myid1 + (myid2 + 1)*nproc1
  idbot = myid1 + (myid2 - 1)*nproc1
  if(myid2 .eq. 0) idbot = myid1 + (nproc2 - 1)*nproc1
  if(myid2 .eq. nproc2-1) idtop = myid1 + 0*nproc1
  allocate(recvtop(is:iep1,jep1:jep2,ks:kep1), &
           recvbot(is:iep1,jsm2:jsm1,ks:kep1))
  call message_pass_nodes(idtop, idbot, &
    iep1-is+1, 2, kep1-ks+1, iep1-is+1, 2, kep1-ks+1, &
    emf2(is:iep1,jem1:je, ks:kep1),emf2(is:iep1,js:jsp1, ks:kep1), &
    recvtop,recvbot)
    emf2(is:iep1,jep1:jep2,ks:kep1) = recvtop
    emf2(is:iep1,jsm2:jsm1,ks:kep1) = recvbot
  deallocate(recvtop,recvbot)
  if(myid2 .eq. 0) then
     emf2(is:iep1,jsm1,ks:kep1) = emf2(is:iep1,js,ks:kep1)
     emf2(is:iep1,jsm2,ks:kep1) = emf2(is:iep1,jsp1,ks:kep1)
  endif
  ! Outer j boundary.
  if(myid2 .eq. nproc2-1) then
     emf2(is:iep1,jep1,ks:kep1) = emf2(is:iep1,je,ks:kep1)
     emf2(is:iep1,jep2,ks:kep1) = emf2(is:iep1,jem1,ks:kep1)
  endif
  !------------------ k - boundary -----------------
  emf2(is:iep1,jsm2:jep2,ksm1) = emf2(is:iep1,jsm2:jep2,ke)
  emf2(is:iep1,jsm2:jep2,ksm2) = emf2(is:iep1,jsm2:jep2,kem1)
  emf2(is:iep1,jsm2:jep2,kep2) = emf2(is:iep1,jsm2:jep2,ksp1)
  emf2(is:iep1,jsm2:jep2,kep3) = emf2(is:iep1,jsm2:jep2,ksp2)
  !-------------------------------------------------------
  !------------------ i - boundary ----------------------
  idtop = myid1 + 1 + myid2*nproc1
  idbot = myid1 - 1 + myid2*nproc1
  if(myid1 .eq. 0) idbot = nproc1 - 1 + nproc1*myid2
  if(myid1 .eq. nproc1-1) idtop = 0 + nproc1*myid2
  allocate(recvtop(iep2:iep3,jsm2:jep2,ksm2:kep3), &
           recvbot(ism2:ism1,jsm2:jep2,ksm2:kep3))
  call message_pass_nodes(idtop, idbot, &
    2, jep2-jsm2+1, kep3-ksm2+1, 2, jep2-jsm2+1, kep3-ksm2+1, &
    emf2(iem1:ie,jsm2:jep2,ksm2:kep3), emf2(isp1:isp2,jsm2:jep2,ksm2:kep3), &
    recvtop,recvbot)
  emf2(iep2:iep3,jsm2:jep2,ksm2:kep3) = recvtop
  emf2(ism2:ism1,jsm2:jep2,ksm2:kep3) = recvbot
  deallocate(recvtop, recvbot)
  ! Inner i boundary.
  if(myid1 .eq. 0) then
     emf2(is,jsm2:jep2,ksm2:kep3) = 0.D0
     emf2(ism1,jsm2:jep2,ksm2:kep3) = - emf2(isp1,jsm2:jep2,ksm2:kep3)
     emf2(ism2,jsm2:jep2,ksm2:kep3) = - emf2(isp2,jsm2:jep2,ksm2:kep3)
  endif
  ! Outer i boundary.
  if(myid1 .eq. nproc1-1) then
     if(noib .eq. 2) then
        emf2(iep2,jsm2:jep2,ksm2:kep3) = &
          emf2(iep1,jsm2:jep2,ksm2:kep3)*2.D0 &
          - emf2(ie,jsm2:jep2,ksm2:kep3)
        emf2(iep3,jsm2:jep2,ksm2:kep3) = &
          emf2(iep2,jsm2:jep2,ksm2:kep3)*2.D0 &
          - emf2(iep1,jsm2:jep2,ksm2:kep3)
     elseif (noib .eq. 5) then
        emf2(iep2,jsm2:jep2,ksm2:kep3) = emf2(ie,jsm2:jep2,ksm2:kep3)
        emf2(iep3,jsm2:jep2,ksm2:kep3) = emf2(iem1,jsm2:jep2,ksm2:kep3)
     else
        emf2(iep1,jsm2:jep2,ksm2:kep3) = 0.D0
        emf2(iep2,jsm2:jep2,ksm2:kep3) = - emf2(ie,jsm2:jep2,ksm2:kep3)
        emf2(iep3,jsm2:jep2,ksm2:kep3) = - emf2(iem1,jsm2:jep2,ksm2:kep3)
     endif
  endif

  return
end subroutine bvalemf2

subroutine bvalemf3(emf3)
  use ModPar,      ONLY: myid1, myid2, nproc1, nproc2, in, jn, kn
  use ModGrid
  use ModBoundary,  ONLY: noib
  implicit none

  real, intent(inout), dimension(in,jn,kn) :: emf3

  integer idtop, idbot
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  !-----------------------------------------------------------------------------
  !------------------ k - boundary -----------------
  emf3(is:iep1,js:jep1,ksm1) = emf3(is:iep1,js:jep1,ke)
  emf3(is:iep1,js:jep1,ksm2) = emf3(is:iep1,js:jep1,kem1)
  emf3(is:iep1,js:jep1,kep1) = emf3(is:iep1,js:jep1,ks)
  emf3(is:iep1,js:jep1,kep2) = emf3(is:iep1,js:jep1,ksp1)
  !------------------ i - boundary -----------------
  idtop = myid1 + 1 + myid2*nproc1
  idbot = myid1 - 1 + myid2*nproc1
  if(myid1 .eq. 0) idbot = nproc1 - 1 + nproc1*myid2
  if(myid1 .eq. nproc1-1) idtop = 0 + nproc1*myid2
  allocate(recvtop(iep2:iep3,js:jep1,ksm2:kep2), &
           recvbot(ism2:ism1,js:jep1,ksm2:kep2))
  call message_pass_nodes(idtop, idbot, &
    2, jep1-js+1, kep2-ksm2+1, 2, jep1-js+1, kep2-ksm2+1, &
    emf3(iem1:ie,js:jep1,ksm2:kep2), emf3(isp1:isp2,js:jep1,ksm2:kep2), &
    recvtop,recvbot)
  emf3(iep2:iep3,js:jep1,ksm2:kep2) = recvtop
  emf3(ism2:ism1,js:jep1,ksm2:kep2) = recvbot
  deallocate(recvtop, recvbot)
  ! Inner i boundary
  if(myid1 .eq. 0) then
     emf3(is,  js:jep1,ksm2:kep2) = 0.D0
     emf3(ism1,js:jep1,ksm2:kep2) = - emf3(isp1,js:jep1,ksm2:kep2)
     emf3(ism2,js:jep1,ksm2:kep2) = - emf3(isp2,js:jep1,ksm2:kep2)
  endif
  ! Outer i boundary
  if(myid1 .eq. nproc1-1) then
     if(noib .eq. 2) then
        emf3(iep2,js:jep1,ksm2:kep2) = &
          emf3(iep1,js:jep1,ksm2:kep2)*2.D0 &
          - emf3(ie,js:jep1,ksm2:kep2)
        emf3(iep3,js:jep1,ksm2:kep2) = &
          emf3(iep2,js:jep1,ksm2:kep2)*2.D0 &
          - emf3(iep1,js:jep1,ksm2:kep2)
     elseif (noib .eq. 5) then
        emf3(iep2,js:jep1,ksm2:kep2) = emf3(ie,js:jep1,ksm2:kep2)
        emf3(iep3,js:jep1,ksm2:kep2) = emf3(iem1,js:jep1,ksm2:kep2)
     else
        emf3(iep1,js:jep1,ksm2:kep2) = 0.D0
        emf3(iep2,js:jep1,ksm2:kep2) = - emf3(ie,js:jep1,ksm2:kep2)
        emf3(iep3,js:jep1,ksm2:kep2) = - emf3(iem1,js:jep1,ksm2:kep2)
     end if
  endif
  !------------------ j - boundary -----------------
  idtop = myid1 + (myid2 + 1)*nproc1
  idbot = myid1 + (myid2 - 1)*nproc1
  if(myid2 .eq. 0) idbot = myid1 + (nproc2 - 1)*nproc1
  if(myid2 .eq. nproc2-1) idtop = myid1 + 0*nproc1
  allocate(recvtop(ism2:iep3,jep2:jep3,ksm2:kep2), &
           recvbot(ism2:iep3,jsm2:jsm1,ksm2:kep2))
  call message_pass_nodes(idtop, idbot, &
    iep3-ism2+1, 2, kep2-ksm2+1, iep3-ism2+1, 2, kep2-ksm2+1, &
    emf3(ism2:iep3,jem1:je,ksm2:kep2),emf3(ism2:iep3,jsp1:jsp2,ksm2:kep2), &
    recvtop,recvbot)
  emf3(ism2:iep3,jep2:jep3,ksm2:kep2) = recvtop
  emf3(ism2:iep3,jsm2:jsm1,ksm2:kep2) = recvbot
  deallocate(recvtop,recvbot)
  ! Inner j boundary.   
  if(myid2 .eq. 0) then
     emf3(ism2:iep3,js  ,ksm2:kep2) = 0.D0
     emf3(ism2:iep3,jsm1,ksm2:kep2) = - emf3(ism2:iep3,jsp1,ksm2:kep2)
     emf3(ism2:iep3,jsm2,ksm2:kep2) = - emf3(ism2:iep3,jsp2,ksm2:kep2)
  endif
  ! Outer j boundary.                      
  if(myid2 .eq. nproc2-1) then
     emf3(ism2:iep3,jep1,ksm2:kep2) = 0.D0
     emf3(ism2:iep3,jep2,ksm2:kep2) = - emf3(ism2:iep3,je,ksm2:kep2)
     emf3(ism2:iep3,jep3,ksm2:kep2) = - emf3(ism2:iep3,jem1,ksm2:kep2)
  endif

  return
end subroutine bvalemf3

subroutine bvaleta(eta)
  use ModPar,          ONLY: myid1, myid2, nproc1, nproc2, in, jn, kn
  use ModGrid
  use ModBoundary,     ONLY: noib
  implicit none

  real, intent(inout), dimension(in,jn,kn) :: eta

  integer :: idtop, idbot
  real, allocatable :: recvtop(:,:,:), recvbot(:,:,:)
  !-------------------------------------------------
  !------------------ i - boundary -----------------
  idtop = myid1 + 1 + myid2*nproc1
  idbot = myid1 - 1 + myid2*nproc1
  if(myid1==0) idbot = nproc1 - 1 + nproc1*myid2
  if(myid1==nproc1-1) idtop = 0 + nproc1*myid2
  allocate(recvtop(iep1:iep2,js:je,ks:ke), recvbot(ism2:ism1,js:je,ks:ke))
  call message_pass_nodes(idtop, idbot,&
    2, je-js+1, ke-ks+1, 2, je-js+1, ke-ks+1, &
    eta(iem1:ie,js:je,ks:ke), eta(is:isp1,js:je,ks:ke), recvtop, recvbot)
  eta(iep1:iep2,js:je,ks:ke) = recvtop
  eta(ism2:ism1,js:je,ks:ke) = recvbot
  deallocate(recvtop, recvbot)
  ! Inner i boundary.
  if(myid1==0) then
     eta(ism1,js:je,ks:ke) = eta(is  ,js:je,ks:ke)
     eta(ism2,js:je,ks:ke) = eta(isp1,js:je,ks:ke)
  endif
  ! Outer i boundary.
  if(myid1==nproc1-1) then
     if(noib==2)then
        eta(iep1,js:je,ks:ke) = eta(ie  ,js:je,ks:ke)
        eta(iep2,js:je,ks:ke) = eta(iep1,js:je,ks:ke)
     else
        eta(iep1,js:je,ks:ke) = eta(ie  ,js:je,ks:ke)
        eta(iep2,js:je,ks:ke) = eta(iem1,js:je,ks:ke)
     endif
  endif
  !------------------ j - boundary -----------------
  idtop = myid1 + (myid2 + 1)*nproc1
  idbot = myid1 + (myid2 - 1)*nproc1
  if(myid2 .eq. 0) idbot = myid1 + (nproc2 - 1)*nproc1
  if(myid2 .eq. nproc2-1) idtop = myid1 + 0*nproc1
  allocate(recvtop(ism2:iep2,jep1:jep2,ks:ke),recvbot(ism2:iep2,jsm2:jsm1,ks:ke))
  call message_pass_nodes(idtop, idbot,&
    iep2-ism2+1, 2, ke-ks+1, iep2-ism2+1, 2, ke-ks+1,&
    eta(ism2:iep2,jem1:je,ks:ke), eta(ism2:iep2,js:jsp1,ks:ke),&
    recvtop, recvbot)
  eta(ism2:iep2,jep1:jep2,ks:ke) = recvtop
  eta(ism2:iep2,jsm2:jsm1,ks:ke) = recvbot
  deallocate(recvtop, recvbot)
  ! Inner j boundary
  if(myid2 .eq. 0) then
     eta(ism2:iep2,jsm1,ks:ke) = eta(ism2:iep2,  js,ks:ke)
     eta(ism2:iep2,jsm2,ks:ke) = eta(ism2:iep2,jsp1,ks:ke)
  endif
  ! Outer j boundary
  if(myid2 .eq. nproc2-1) then
     eta(ism2:iep2,jep1,ks:ke) = eta(ism2:iep2,  je,ks:ke)
     eta(ism2:iep2,jep2,ks:ke) = eta(ism2:iep2,jem1,ks:ke)
  endif
  !------------------ k - boundary -----------------
  ! Inner k boundary
  eta(ism2:iep2,jsm2:jep2,ksm1) = eta(ism2:iep2,jsm2:jep2,  ke)
  eta(ism2:iep2,jsm2:jep2,ksm2) = eta(ism2:iep2,jsm2:jep2,kem1)
  ! Outer k boundary
  eta(ism2:iep2,jsm2:jep2,kep1) = eta(ism2:iep2,jsm2:jep2,  ks)
  eta(ism2:iep2,jsm2:jep2,kep2) = eta(ism2:iep2,jsm2:jep2,ksp1)

end subroutine bvaleta

subroutine bvalb
  use ModPar,        ONLY: myid1, nproc1
  use ModGrid
  use ModField,      ONLY: b1, b2, b3
  use ModBoundary,   ONLY: noib
  implicit none
  integer i, j, k
!----------------------------------------------
!      special boundary conditions for B
!
!      Outer i boundary.
!
       if(myid1 .eq. nproc1-1) then

         if(noib .eq. 5) then
           do k=ksm2,kep2
             do j=jsm2,jep2
               b1(iep2,j,k)=b1(ie,j,k)
               b2(iep1,j,k)=-b2(ie,j,k)
               b2(iep2,j,k)=-b2(iem1,j,k)
               b3(iep1,j,k)=-b3(ie,j,k)
               b3(iep2,j,k)=-b3(iem1,j,k)
             enddo
           enddo
         endif

       endif

end subroutine bvalb

end module ModBval
