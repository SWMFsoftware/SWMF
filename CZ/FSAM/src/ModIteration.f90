module ModIteration

  implicit none

  private

  public :: nudt

  public :: nudt_cond

  public :: pred_corr_step

  public :: conductstep

contains

  !===================================================================================

  subroutine nudt
    use ModPar,    ONLY: huge, tiny, myid1, in
    use ModGrid,   ONLY: dx1a, dx2a, dx3a, g2b, g32b, is, ie, js, je, ks, ke
    use ModBack,   ONLY: d0, ovrevar, ovrmvar, ovbvfq
    use ModField,  ONLY: b1, b2, b3, v1, v2, v3
    use ModSundry, ONLY: dt, courno
    use ModMpi
    use ModFSAM
    implicit none

    integer :: ierr, i, j, k
    real    :: dtnewlc, dtnew
    real    :: q1, va, vv1, vv2, vv3, cmax, dr1, dr2, dr3, drmin
    !--------------------------------------------------------------------------------

    dtnewlc = huge
    do k=ks,ke; do j=js,je; do i=is,ie
       va = sqrt((abs(b1(i,j,k)) + abs(b1(i+1,j,k)))**2 + &
            (abs(b2(i,j,k)) + abs(b2(i,j+1,k)))**2 + &
            (abs(b3(i,j,k)) + abs(b3(i,j,k+1)))**2)*0.5d0/sqrt(d0(i+myid1*(in-5)))
       vv1 = 0.5d0*(abs(v1(i,j,k)) + abs(v1(i+1,j,k)))
       vv2 = 0.5d0*(abs(v2(i,j,k)) + abs(v2(i,j+1,k)))
       vv3 = 0.5d0*(abs(v3(i,j,k)) + abs(v3(i,j,k+1)))
       cmax = sqrt((vv1+va)**2 + (vv2+va)**2 + (vv3+va)**2)
       dr1 = dx1a(i)
       dr2 = dx2a(j)*g2b(i)
       dr3 = dx3a(k)*g2b(i)*g32b(j)
       drmin = min(dr1, dr2, dr3)
       q1 = max(2.d0*ovrevar(i+myid1*(in-5)), 2.d0*ovrmvar(i+myid1*(in-5)))
       dtnewlc = min(dtnewlc, drmin/(cmax + 2.d0*q1/drmin + tiny))
    enddo; enddo; enddo

    call MPI_ALLREDUCE(dtnewlc, dtnew, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
         iComm, ierr)
    dt = min(courno*dtnew, ovbvfq*1.d-2)

  end subroutine nudt

  !===================================================================================

  subroutine nudt_cond
    use ModPar,     ONLY: huge, tiny, myid1, in
    use ModGrid,    ONLY: dx1a, dx2a, dx3a, g2b, g32b, is, ie, js, je, ks, ke
    use ModBack,    ONLY: ovrth
    use ModSundry,  ONLY: dtcond, courno
    use ModMpi
    use ModFSAM
    implicit none

    integer :: ierr, i, j, k
    real    :: dtnewlc, dtnew, q1, dr1, dr2, dr3, drmin
    !-------------------------------------------------------------------------------

    dtnewlc = huge
    do k=ks,ke; do j=js,je; do i=is,ie
       dr1 = dx1a(i)
       dr2 = dx2a(j)*g2b(i)
       dr3 = dx3a(k)*g2b(i)*g32b(j)
       drmin = min(dr1, dr2, dr3)
       q1 = ovrth(i+myid1*(in-5))
       dtnewlc = min(dtnewlc,0.25D0*drmin**2/(q1+tiny))
    enddo; enddo; enddo

    call MPI_ALLREDUCE(dtnewlc, dtnew, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
         iComm, ierr)
    dtcond = courno*dtnew

  end subroutine nudt_cond

  !===================================================================================

  subroutine pred_corr_step
    use ModGrid,    ONLY: ism2,iep2,jsm2,jep2,ksm2,kep2
    use ModField,   ONLY: v1, v2, v3, b1, b2, b3, s, &
         v1_prev, v2_prev, v3_prev, b1_prev, b2_prev, b3_prev, s_prev
    use ModRHS,     ONLY: p_init
    implicit none

    integer :: i, j, k
    real    :: tFactor
    !--------------------------------------------------------------------------------

    ! copy dependent variables
    do k=ksm2,kep2; do j=jsm2,jep2; do i=ism2,iep2
       s_prev(i,j,k)  = s(i,j,k)
       v1_prev(i,j,k) = v1(i,j,k)
       v2_prev(i,j,k) = v2(i,j,k)
       v3_prev(i,j,k) = v3(i,j,k)
       b1_prev(i,j,k) = b1(i,j,k)
       b2_prev(i,j,k) = b2(i,j,k)
       b3_prev(i,j,k) = b3(i,j,k)
    enddo; enddo; enddo

    tFactor = 0.5D0
    call update_states(tFactor)
    call p_init
    tFactor = 1.D0
    call update_states(tFactor)

  end subroutine pred_corr_step

  !===================================================================================

  subroutine conductstep(ntcond)
    use ModPar
    use ModGrid
    use ModField,      ONLY: s, s_prev
    use ModSundry,     ONLY: dt, time, ctime, dtcond
    use ModDel
    use ModRHS,        ONLY: p_init, heattran
    use ModBval,       ONLY: bvals
    implicit none

    integer, intent(out) :: ntcond
    integer :: i, j, k, iflg
    real    :: delt, dtc
    ! paramters for testing, write ModDel module data
    character(len=4) :: idcpu
    logical :: DoWriteDel = .false.
    !---------------------------------------------------------------------------------

    ntcond = 0
    delt   = 0.D0
    iflg   = -1
    do while (iflg .ne. 1)
       dtc=dtcond
       if (dtc .ge. dt-delt) then
          dtc=dt-delt
          iflg=1
       endif
       do k=ksm2,kep2; do j=jsm2,jep2; do i=ism2,iep2
          s_prev(i,j,k) = s(i,j,k)
          dels(i,j,k)   = 0.D0
       enddo; enddo; enddo
       call heattran
       do k=ks,ke; do j=js,je; do i=is,ie
          s(i,j,k) = s_prev(i,j,k) + dels(i,j,k)*0.5D0*dtc
       enddo; enddo; enddo
       ctime = time + delt + 0.5D0*dtc
       call bvals

       do k=ksm2,kep2; do j=jsm2,jep2; do i=ism2,iep2
          dels (i,j,k) = 0.D0
       enddo; enddo; enddo
       call heattran
       do k=ks,ke; do j=js,je; do i=is,ie
          s(i,j,k) = s_prev(i,j,k) + dels(i,j,k)*dtc
       enddo; enddo; enddo
       ctime = time + delt + dtc
       call bvals

       delt   = delt + dtc
       ntcond = ntcond + 1
    enddo

    call p_init

    if (DoWriteDel) then
       write(idcpu,'(i4.4)') myid
       open(unit=17,file='emf.cpu'//idcpu//'.cond',form='unformatted', &
            access='stream')
       write(17) in, jn, kn
       write(17) (((emf1s(i,j,k),i=ism2,iep2),j=jsm2,jep3),k=ksm2,kep3)
       write(17) (((emf2s(i,j,k),i=ism2,iep3),j=jsm2,jep2),k=ksm2,kep3)
       write(17) (((emf3s(i,j,k),i=ism2,iep3),j=jsm2,jep3),k=ksm2,kep2)
       close(17)
       open(unit=17,file='dels.cpu'//idcpu,form='unformatted', &
            access='stream')
       write(17) in, jn, kn
       write(17) (((dels1(i,j,k),i=is,iep1),j=js,je),k=ks,ke)
       write(17) (((dels2(i,j,k),i=is,ie),j=js,jep1),k=ks,ke)
       write(17) (((dels3(i,j,k),i=is,ie),j=js,je),k=ks,kep1)
       write(17) (((dels(i,j,k),i=is,ie),j=js,je),k=ks,ke)
       close(17)
    endif

  end subroutine conductstep

  !===================================================================================

  subroutine update_states(tFactor)
    use ModPar
    use ModGrid
    use ModField
    use ModBack,    ONLY: fact, d0
    use ModSundry,  ONLY: dt, time, ctime
    use ModDel,     ONLY: emf1s, emf2s, emf3s, dels1, dels2, dels3, dels
    use ModRHS,     ONLY: hsmoc, diffemf, resistive_emf
    use ModBval,    ONLY: bvalb, bvalv, bvals
    implicit none

    real, intent(in) :: tFactor
    integer :: i, ip1, j, jp1, k, kp1
    real    :: q1, q2
    !--------------------------------------------------------------------------------

    call hsmoc
    call diffemf
    call resistive_emf

    !-----------------------------> Update b1 <-----------------------------
    do k=ksm2,kep2; do j=jsm2,jep2; do i=ism2,iep2
       b1(i,j,k) = b1_prev(i,j,k) + tFactor*dt* &
            ( emf3s(i,j+1,k)*dx3a(k)*g31a(i)*g32a(j+1) &
            - emf3s(i,j,k)*dx3a(k)*g31a(i)*g32a(j) &
            - emf2s(i,j,k+1)*dx2a(j)*g2a(i) &
            + emf2s(i,j,k)*dx2a(j)*g2a(i) ) &
            *g2ai(i)*g31ai(i)*g32bi(j)*dx2ai(j)*dx3ai(k)
    enddo; enddo; enddo

    !-----------------------------> Update b2 <-----------------------------
    do k=ksm2,kep2; do j=jsm2,jep2; do i=ism2,iep2
       b2(i,j,k) = b2_prev(i,j,k) + tFactor*dt* &
            ( emf1s(i,j,k+1)*dx1a(i) - emf1s(i,j,k)*dx1a(i) &
            - emf3s(i+1,j,k)*dx3a(k)*g31a(i+1)*g32a(j) &
            + emf3s(i,j,k)*dx3a(k)*g31a(i)*g32a(j) ) &
            *g31bi(i)*g32ai(j)*dx1ai(i)*dx3ai(k)
    enddo; enddo; enddo

    !-----------------------------> Update b3 <-----------------------------
    do k=ksm2,kep2; do j=jsm2,jep2; do i=ism2,iep2
       b3(i,j,k) = b3_prev(i,j,k) + tFactor*dt* &
            ( emf2s(i+1,j,k)*dx2a(j)*g2a(i+1) - emf2s(i,j,k)*dx2a(j)*g2a(i) &
            - emf1s(i,j+1,k)*dx1a(i) + emf1s(i,j,k)*dx1a(i) ) &
            *g2bi(i)*dx1ai(i)*dx2ai(j)
    enddo; enddo; enddo

    call bvalb

    ! update velocity and entropy
    do k=ks,ke; do j=js,je; do i=is,ie
       q1 = 0.5D0*(1.D0/fact(myid1*(in-5)+i) + 1.D0/fact(myid1*(in-5)+i-1))
       q2 = 1.D0/fact(myid1*(in-5)+i)
       v1(i,j,k) = v1_prev(i,j,k) + ( dels1(i,j,k) &
            - (p(i,j,k) - p(i-1,j,k))*dx1bi(i)*q1) &
            *tFactor*dt/(0.5D0*(d0(myid1*(in-5)+i) + d0(myid1*(in-5)+i-1)))
       v2(i,j,k) = v2_prev(i,j,k) + ( dels2(i,j,k) &
            - (p(i,j,k) - p(i,j-1,k))*dx2bi(j)*g2bi(i)*q2) &
            *tFactor*dt/d0(myid1*(in-5)+i)
       v3(i,j,k) = v3_prev(i,j,k) + ( dels3(i,j,k) &
            - (p(i,j,k)-p(i,j,k-1))*dx3bi(k)*g31bi(i)*g32bi(j)*q2) &
            *tFactor*dt/d0(myid1*(in-5)+i)
       s(i,j,k) = s_prev(i,j,k) + dels(i,j,k)*tFactor*dt
    enddo; enddo; enddo
    ctime = time + tFactor*dt
    call bvalv
    call bvals

  end subroutine update_states

end module ModIteration
