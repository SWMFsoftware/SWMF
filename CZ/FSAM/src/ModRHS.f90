module ModRHS

  implicit none

  private

  public :: p_init

  public :: hsmoc

  public :: diffemf

  public :: resistive_emf

  public :: heattran

contains

  !=====================================================================================

  subroutine p_init
    use ModGetP,     ONLY: getp
    implicit none

    !----------------------------------------------------------------------------------

    call stv

    call lorentz

    call advection

    call viscous

    call radheat

    call getp

  end subroutine p_init

  !=====================================================================================

  subroutine stv
    use ModPar,     ONLY: myid1, in, jn, kn
    use ModGrid
    use ModBack,    ONLY: d0, gamma, gacc_str
    use ModField,   ONLY: v2, v3, s
    use ModDel,     ONLY: dels1, dels2, dels3
    use ModCoeff,   ONLY: ovrb
    implicit none

    integer :: i, j, k
    real    :: rho, gammaloc, gaccloc
    real, allocatable, dimension(:,:,:) :: f1, f2
    !-----------------------------------------------------------------------------------

    ! compute rotational pseudo-force and coriolis force in cell center
    allocate(f1(1:in,1:jn,1:kn), f2(1:in,1:jn,1:kn))
    do k=ksm1,kep1; do j=jsm1,jep1; do i=ism1,iep1
       f1(i,j,k) = 0.25d0*(v2(i,j,k) + v2(i,j+1,k))**2*g2bi(i)*dg2ad1(i) + &
            0.25d0*(v3(i,j,k) + v3(i,j,k+1))**2*g31bi(i)*dg31ad1(i) + &
            ovrb*(v3(i,j,k) + v3(i,j,k+1))*sin(x2b(j))
       f2(i,j,k) = 0.25d0*(v3(i,j,k) + v3(i,j,k+1))**2*g32bi(j)*g2bi(i)* &
            dg32ad2(j) + ovrb*(v3(i,j,k) + v3(i,j,k+1))*cos(x2b(j))
    enddo; enddo; enddo

    do k=ks,ke
       do j=js,je; do i=is,iep1
          ! Add source terms to "dels1"
          rho = 0.5d0*(d0(i-1+myid1*(in-5)) + d0(i+myid1*(in-5)))
          gammaloc = 0.5d0*(gamma(i-1+myid1*(in-5)) + gamma(i+myid1*(in-5)))
          gaccloc = 0.5d0*(gacc_str(i-1+myid1*(in-5)) + gacc_str(i+myid1*(in-5)))
          dels1(i,j,k) = (1.d0 - 1.d0/gammaloc)*rho*gaccloc* &
               (s(i-1,j,k) + s(i,j,k))*0.5d0 + &
               0.5d0*(f1(i-1,j,k) + f1(i,j,k))*rho
       enddo; enddo
       do j = js, jep1; do i = is, ie
          ! Add source terms to "dels2"
          rho = d0(i+myid1*(in-5))
          dels2(i,j,k)  = 0.5d0*(f2(i,j-1,k) + f2(i,j,k))*rho
       end do;  end do
    enddo
    ! Add source terms to "dels3"
    dels3(is:ie,js:je,ks:kep1) = 0.d0

    deallocate(f1, f2)

  end subroutine stv

  !=====================================================================================
  subroutine lorentz
    use ModPar,     ONLY: ijkn, in, jn, kn
    use ModGrid
    use ModField,   ONLY: b1, b2, b3
    use ModDel,     ONLY: dels1, dels2, dels3
    implicit none

    integer :: i, j, k
    real    :: q1
    real, dimension(1:ijkn) :: scra1, scra2, scra3
    real, dimension(1:in,1:jn,1:kn) :: b12, b23, b13
    !----------------------------------------------------------------------------------

    do k=ks,kep1; do i=is,iep1; do j=js,jep1
       b12(i,j,k) = 0.25d0*(b1(i,j-1,k) + b1(i,j,k))*(b2(i,j,k) + b2(i-1,j,k))
       b13(i,j,k) = 0.25d0*(b3(i,j,k) + b3(i-1,j,k))*(b1(i,j,k-1) + b1(i,j,k))
       b23(i,j,k) = 0.25d0*(b3(i,j,k) + b3(i,j-1,k))*(b2(i,j,k-1) + b2(i,j,k))
    end do; enddo; enddo

    ! 1-force
    do k=ks,ke; do j=js,je; do i=is,iep1
       q1 = g31ai(i)*g32bi(j)
       dels1(i,j,k) = dels1(i,j,k) + (g32a(j+1)*b12(i,j+1,k) - &
            g32a(j)*b12(i,j,k))*g2ai(i)*g32bi(j)*dx2ai(j) + &
            (b13(i,j,k+1) - b13(i,j,k))*q1*dx3ai(k)
    end do; end do; end do
    do k=ks,ke; do j=js,je
       do i=ism1,iep1
          scra1(i) = g2b(i)**2*(b2(i,j,k)**2 + b2(i,j+1,k)**2)
          scra2(i) = g31b(i)**2*(b3(i,j,k)**2 + b3(i,j,k+1)**2)
          scra3(i) = (g2b(i)*g31b(i)*(b1(i+1,j,k) + b1(i,j,k)))**2*0.5d0
       enddo
       do i=is,iep1
          dels1(i,j,k) = dels1(i,j,k) - 0.25d0*dx1bi(i)*( &
               g2ai(i)**2*(scra1(i) - scra1(i-1)) + &
               g31ai(i)**2*(scra2(i) -scra2(i-1)) - &
               (g2ai(i)*g31ai(i))**2*(scra3(i) - scra3(i-1)))
       enddo
    enddo; enddo

    ! 2-force
    do k=ks,ke; do j=js,jep1; do i=is,ie
       q1 = g31bi(i)*g32ai(j)
       dels2(i,j,k) = dels2(i,j,k) + (b23(i,j,k+1) - b23(i,j,k))*q1*dx3ai(k) + &
            (g2a(i+1)*g31a(i+1)*b12(i+1,j,k)*g2a(i+1) - &
            g2a(i)*g2a(i)*g31a(i)*b12(i,j,k))*g2bi(i)*g2bi(i)*g31bi(i)*dx1ai(i)
    end do; end do; end do
    do k=ks,ke; do i=is,ie
       do j=jsm1,jep1
          scra1(j) = g32b(j)**2*(b3(i,j,k)**2 + b3(i,j,k+1)**2)
          scra2(j) = b1(i,j,k)**2 + b1(i+1,j,k)**2
          scra3(j) = (g32b(j)*(b2(i,j,k) + b2(i,j+1,k)))**2*0.5d0
       enddo
       do j=js,jep1
          dels2(i,j,k) = dels2(i,j,k) - 0.25d0*dx2bi(j)*g2bi(i)*( &
               g32ai(j)**2*(scra1(j) - scra1(j-1)) + &
               (scra2(j) - scra2(j-1)) - &
               g32ai(j)**2*(scra3(j) - scra3(j-1)))
       enddo
    enddo; enddo

    ! 3-force
    do k=ks,kep1; do j=js,je; do i=is,ie
       dels3(i,j,k) = dels3(i,j,k) + (g2a(i+1)*g31a(i+1)*b13(i+1,j,k)*g31a(i+1) - &
            g31a(i)*g2a(i)*g31a(i)*b13(i,j,k))*g2bi(i)*g31bi(i)*g31bi(i)*dx1ai(i) + &
            (g32a(j+1)*b23(i,j+1,k)*g32a(j+1) - g32a(j)*g32a(j)*b23(i,j,k)) &
            *g32bi(j)*g32bi(j)*g2bi(i)*dx2ai(j)
    end do; end do; end do
    do j=js,je; do i=is,ie
       do k=ksm1,kep1
          scra1(k) = b1(i,j,k)**2 + b1(i+1,j,k)**2
          scra2(k) = b2(i,j,k)**2 + b2(i,j+1,k)**2
          scra3(k) = (b3(i,j,k) + b3(i,j,k+1))**2*0.5d0
       enddo
       do k=ks,kep1
          dels3(i,j,k) = dels3(i,j,k) - 0.25d0*dx3bi(k)*g31bi(i)*g32bi(j)*( &
               scra1(k) - scra1(k-1) + scra2(k) - scra2(k-1) - scra3(k)+scra3(k-1))
       enddo
    enddo; enddo
  end subroutine lorentz

  !=======================================================================================
 subroutine advection
    use ModPar,     ONLY: ijkn, myid1, tiny, in, jn, kn, ijkn, myid
    use ModGrid
    use ModField,   ONLY: b1, b2, b3, v1, v2, v3, s
    use ModBack,    ONLY: d0, temp0, s0, sr0, gamma, gacc_str
    use ModDel,     ONLY: dels1, dels2, dels3, dels
    use ModSundry,  ONLY: nlf_v, alpha_v, nlf_s, alpha_s
    use ModCoeff,   ONLY: ovrb
    use ModInterp,  ONLY: xtvd, xtvd_a
    implicit none

    integer :: i, j, k, nlf
    real, dimension(1:ijkn) :: qint_v1, qint_v2g2i, qint_v3g3i, qint_s, &
         qlr_v1, qlr_v2g2i, qlr_v3g3i, qlr_s, dq_v1, dq_v2g2i, dq_v3g3i, dq_s
    real, dimension(1:ijkn) :: qint, dq, qlr, dqint
    real, dimension(1:ijkn) :: flux
    real :: alpha, q1, qq, rr, cmax
    real, dimension(1:in,1:jn,1:kn) :: cmaxc, cmax12, cmax13, cmax23, cmax1, cmax2, cmax3, &
         scratchc, scratch12, scratch13, scratch23
    real, dimension(1:ijkn) :: vtmp
    real, dimension(1:in) :: d0m, temp0m
    !---------------------------------------------------------------------------------

    nlf   = nlf_v
    alpha = alpha_v
    do i=ism1,iep1
       temp0m(i) = 0.5d0*(temp0(i-1+myid1*(in-5)) + temp0(i+myid1*(in-5)))
       d0m(i)    = 0.5d0*(d0(i-1+myid1*(in-5)) + d0(i+myid1*(in-5)))
    enddo

    do k=ksm1,kep1; do j=jsm1,jep1; do i=ism1,iep1
       cmaxc(i,j,k) = sqrt((0.5d0*(b1(i,j,k) + b1(i+1,j,k)))**2 + &
               (0.5d0*(b2(i,j,k) + b2(i,j+1,k)))**2 + &
               (0.5d0*(b3(i,j,k) + b3(i,j,k+1)))**2)/sqrt(d0(i+myid1*(in-5)))
       cmax12(i,j,k) = sqrt((0.5d0*(b1(i,j-1,k) + b1(i,j,k)))**2 + &
            (0.5d0*(b2(i-1,j,k) + b2(i,j,k)))**2 + &
            (0.125d0*(b3(i,j,k) + b3(i,j,k+1) + b3(i-1,j,k) + b3(i-1,j,k+1) + &
            b3(i-1,j-1,k) + b3(i-1,j-1,k+1) + b3(i,j-1,k) + b3(i,j-1,k+1)))**2)/ &
            sqrt(d0m(i))
       cmax13(i,j,k) = sqrt((0.5d0*(b1(i,j,k-1) + b1(i,j,k)))**2 + &
            (0.5d0*(b3(i-1,j,k) + b3(i,j,k)))**2 + &
            (0.125d0*(b2(i,j,k) + b2(i,j+1,k) + b2(i-1,j,k) + b2(i-1,j+1,k) + &
            b2(i-1,j,k-1) + b2(i-1,j+1,k-1) + b2(i,j,k-1) + b2(i,j+1,k-1)))**2)/ &
            sqrt(d0m(i))
       cmax23(i,j,k) = sqrt((0.5d0*(b2(i,j,k-1) + b2(i,j,k)))**2 + &
            (0.5d0*(b3(i,j-1,k) + b3(i,j,k)))**2 + &
            (0.125d0*(b1(i,j,k) + b1(i+1,j,k) + b1(i,j-1,k) + b1(i+1,j-1,k) + &
            b1(i,j-1,k-1) + b1(i+1,j-1,k-1) + b1(i,j,k-1) + b1(i+1,j,k-1)))**2)/ &
            sqrt(d0(i+myid1*(in-5)))
       cmax1(i,j,k) = sqrt(b1(i,j,k)**2 + &
            (0.25d0*(b2(i,j,k) + b2(i,j+1,k) + b2(i-1,j,k) + b2(i-1,j+1,k)))**2 + &
            (0.25d0*(b3(i,j,k) + b3(i,j,k+1) + b3(i-1,j,k) + b3(i-1,j,k+1)))**2)/ &
            sqrt(d0m(i))
       cmax2(i,j,k) = sqrt(b2(i,j,k)**2 + &
            (0.25d0*(b3(i,j,k) + b3(i,j,k+1) + b3(i,j-1,k) + b3(i,j-1,k+1)))**2 + &
            (0.25d0*(b1(i,j,k) + b1(i+1,j,k) + b1(i,j-1,k) + b1(i+1,j-1,k)))**2)/ &
            sqrt(d0(i+myid1*(in-5)))
       cmax3(i,j,k) = sqrt(b3(i,j,k)**2 + &
            (0.25d0*(b1(i,j,k) + b1(i+1,j,k) + b1(i,j,k-1) + b1(i+1,j,k-1)))**2 + &
            (0.25d0*(b2(i,j,k) + b2(i,j+1,k) + b2(i,j,k-1) + b2(i,j+1,k-1)))**2)/ &
            sqrt(d0(i+myid1*(in-5)))
    enddo; enddo; enddo

    ! advection for dels1
    do k=ks,ke; do j=js,je
       qint(ism2:iep3) = v1(ism2:iep3,j,k)
       call xtvd_a(in, ijkn, dx1a, dx1ai, qint, dq, qlr, dqint)
       vtmp(ism1:iep1) = 0.5d0*(v1(ism1:iep1,j,k) + v1(is:iep2,j,k))
       do i=ism1,iep1
          cmax = cmaxc(i,j,k)
          rr = abs(dq(i))/(abs(dqint(i)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          q1 = 0.5d0*(abs(vtmp(i)) + cmax*qq)* &
               d0(i+myid1*(in-5))*dq(i)
          flux(i) = (d0(i+myid1*(in-5))*vtmp(i)*qlr(i) - &
               q1)*g2b(i)*g31b(i)
          scratchc(i,j,k) = q1*dqint(i)*dx1ai(i) 
       enddo
       do i=is,iep1
          dels1(i,j,k) = dels1(i,j,k) - (flux(i) - flux(i-1))*dx1bi(i)*g2ai(i)*g31ai(i)
       enddo
    enddo; enddo
    do k=ks,ke; do i=is,iep1
       qint(jsm2:jep2) = v1(i,jsm2:jep2,k)
       call xtvd(jn, ijkn, dx2a, dx2bi, qint, dq, qlr, dqint)
       vtmp(1:jn) = 0.5d0*(v2(i-1,:,k) + v2(i,:,k))
       do j=js,jep1
          cmax = cmax12(i,j,k)
          rr = abs(dq(j))/(abs(dqint(j)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          q1 = 0.5d0*(abs(vtmp(j)) + cmax*qq)*d0m(i)*dq(j)
          flux(j) = (d0m(i)*vtmp(j)*qlr(j) - q1)*g32a(j)
          scratch12(i,j,k) = q1*dqint(j)*g2ai(i)*dx2bi(j) 
       enddo
       do j=js,je
          dels1(i,j,k) = dels1(i,j,k) - (flux(j+1) - flux(j))*dx2ai(j)*g2ai(i)*g32bi(j)
       enddo
    enddo; enddo
    do j=js,je; do i=is,iep1
       qint(ksm2:kep2) = v1(i,j,ksm2:kep2)
       call xtvd(kn, ijkn, dx3a, dx3bi, qint, dq, qlr, dqint)
       vtmp(1:kn) = 0.5d0*(v3(i-1,j,:) + v3(i,j,:))
       do k=ks,kep1
          cmax = cmax13(i,j,k)
          rr = abs(dq(k))/(abs(dqint(k)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          q1 = 0.5d0*(abs(vtmp(k)) + cmax*qq)*d0m(i)*dq(k)
          flux(k) = d0m(i)*vtmp(k)*qlr(k) - q1
          scratch13(i,j,k) = q1*dqint(k)*g31ai(i)*g32bi(j)*dx3bi(k)
       enddo
       do k=ks,ke
          dels1(i,j,k) = dels1(i,j,k) - (flux(k+1) - flux(k))*dx3ai(k)*g31ai(i)*g32bi(j)
       enddo
    enddo; enddo

    ! advection for dels2
    do k=ks,ke; do j=js,jep1
       qint(ism2:iep2) = v2(ism2:iep2,j,k)*g2bi(ism2:iep2)
       call xtvd(in, ijkn, dx1a, dx1bi, qint, dq, qlr, dqint)
       vtmp(1:in) = 0.5d0*(v1(:,j-1,k) + v1(:,j,k))
       do i=is,iep1
          cmax = cmax12(i,j,k)
          rr = abs(dq(i))/(abs(dqint(i)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          q1 = 0.5d0*(abs(vtmp(i)) + cmax*qq)*d0m(i)*dq(i)
          flux(i) = (d0m(i)*vtmp(i)*qlr(i) - q1) &
               *g2a(i)*g31a(i)*(g2a(i))**2
          scratch12(i,j,k) = scratch12(i,j,k) + q1*(g2a(i))**2*dqint(i)*dx1bi(i)
           enddo
       do i=is,ie
          dels2(i,j,k) = dels2(i,j,k) - (flux(i+1) - flux(i))* &
               dx1ai(i)*g2bi(i)*g31bi(i)*g2bi(i)
       enddo
    enddo; enddo
    do k=ks,ke; do i=is,ie
       qint(jsm2:jep3) = v2(i,jsm2:jep3,k)*g2bi(i)
       call xtvd_a(jn, ijkn, dx2a, dx2ai, qint, dq, qlr, dqint)
       vtmp(jsm1:jep1) = 0.5d0*(v2(i,jsm1:jep1,k) + v2(i,js:jep2,k))
       do j=jsm1,jep1
          cmax = cmaxc(i,j,k)
          rr = abs(dq(j))/(abs(dqint(j)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          q1 = 0.5d0*(abs(vtmp(j)) + cmax*qq)* &
               d0(i+myid1*(in-5))*dq(j)
          flux(j) = (d0(i+myid1*(in-5))*vtmp(j)*qlr(j) - q1) &
               *g32b(j)*(g2b(i))**2
          scratchc(i,j,k) = scratchc(i,j,k) + q1*(g2b(i))**2*dqint(j)*g2bi(i)*dx2ai(j)
       enddo
       do j=js,jep1
          dels2(i,j,k) = dels2(i,j,k) - (flux(j) - flux(j-1))* &
               dx2bi(j)*g2bi(i)*g32ai(j)*g2bi(i)
       enddo
    enddo; enddo
    do j=js,jep1; do i=is,ie
       qint(ksm2:kep2) = v2(i,j,ksm2:kep2)*g2bi(i)
       call xtvd(kn, ijkn, dx3a, dx3bi, qint, dq, qlr, dqint)
       vtmp(1:kn) = 0.5d0*(v3(i,j-1,:) + v3(i,j,:))
       do k=ks,kep1
          cmax = cmax23(i,j,k)
          rr = abs(dq(k))/(abs(dqint(k)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          q1 = 0.5d0*(abs(vtmp(k)) + cmax*qq)* &
               d0(i+myid1*(in-5))*dq(k)
          flux(k) = (d0(i+myid1*(in-5))*vtmp(k)*qlr(k) - q1) &
               *(g2b(i))**2
          scratch23(i,j,k) = q1*(g2b(i))**2*dqint(k)*g31bi(i)*g32ai(j)*dx3bi(k) 
       enddo
       do k=ks,ke
          dels2(i,j,k) = dels2(i,j,k) - (flux(k+1) - flux(k))* &
               dx3ai(k)*g31bi(i)*g32ai(j)*g2bi(i)
       enddo
    enddo; enddo

    ! advection for dels3   
    do k=ks,kep1; do j=js,je
       qint(ism2:iep2) = v3(ism2:iep2,j,k)*g31bi(ism2:iep2)*g32bi(j)
       call xtvd(in, ijkn, dx1a, dx1bi, qint, dq, qlr, dqint)
       vtmp(1:in) = 0.5d0*(v1(:,j,k-1) + v1(:,j,k))
       do i = is,iep1
          cmax = cmax13(i,j,k)
          rr = abs(dq(i))/(abs(dqint(i)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          q1 = 0.5d0*(abs(vtmp(i)) + cmax*qq)*d0m(i)*dq(i)
          flux(i) = (d0m(i)*vtmp(i)*(qlr(i) + ovrb) - q1)&
               *g2a(i)*g31a(i)*(g31a(i)*g32b(j))**2
          scratch13(i,j,k) = scratch13(i,j,k) + q1*(g31a(i)*g32b(j))**2*dqint(i)*dx1bi(i)
       enddo
       do i=is,ie
          dels3(i,j,k) = dels3(i,j,k) - (flux(i+1) - flux(i))* &
               dx1ai(i)*g2bi(i)*g31bi(i)*g31bi(i)*g32bi(j)
       enddo
    enddo; enddo
    do k=ks,kep1; do i=is,ie
       qint(jsm2:jep2) = v3(i,jsm2:jep2,k)*g31bi(i)*g32bi(jsm2:jep2)
       vtmp(1:jn) = 0.5d0*(v2(i,:,k-1) + v2(i,:,k))
       call xtvd(jn, ijkn, dx2a, dx2bi, qint, dq, qlr, dqint)
       do j=js,jep1
          cmax = cmax23(i,j,k)
          rr = abs(dq(j))/(abs(dqint(j)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          q1 = 0.5d0*(abs(vtmp(j)) + cmax*qq)* &
               d0(i+myid1*(in-5))*dq(j)
          flux(j) = (d0(i+myid1*(in-5))*vtmp(j)* &
               (qlr(j) + ovrb)- q1)*g32a(j)*(g31b(i)*g32a(j))**2
          scratch23(i,j,k) = scratch23(i,j,k) + q1*(g31b(i)*g32a(j))**2*dqint(j)*g2bi(i)*dx2bi(j)
       enddo
       do j=js,je
          dels3(i,j,k) = dels3(i,j,k) - (flux(j+1) - flux(j))* &
               dx2ai(j)*g2bi(i)*g32bi(j)*g31bi(i)*g32bi(j)
       enddo
    enddo; enddo
    do j=js,je; do i=is,ie
       qint(ksm2:kep3) = v3(i,j,ksm2:kep3)*g31bi(i)*g32bi(j)
       call xtvd_a(kn, ijkn, dx3a, dx3ai, qint, dq, qlr, dqint)
       vtmp(ksm1:kep1) = 0.5d0*(v3(i,j,ksm1:kep1) + v3(i,j,ks:kep2))
       do k=ksm1,kep1
          cmax = cmaxc(i,j,k)
          rr = abs(dq(k))/(abs(dqint(k)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          q1 = 0.5d0*(abs(vtmp(k)) + cmax*qq)* &
               d0(i+myid1*(in-5))*dq(k)
          flux(k) = (d0(i+myid1*(in-5))*vtmp(k)*&
               (qlr(k) + ovrb) - q1)*(g31b(i)*g32b(j))**2
          scratchc(i,j,k) = scratchc(i,j,k) + &
               q1*(g31b(i)*g32b(j))**2*dqint(k)*g31bi(i)*g32bi(j)*dx3ai(k)
       enddo
       do k=ks,kep1
          dels3(i,j,k) = dels3(i,j,k) - (flux(k) - flux(k-1))* &
               dx3bi(k)*g31bi(i)*g32bi(j)*g31bi(i)*g32bi(j)
       enddo
    enddo; enddo
    ! add heating due to slope limited diffusion 
    do k=ks,ke; do j=js,je; do i=is,ie
       dels(i,j,k) = 1.d0/d0(i+myid1*(in-5))/temp0(i+myid1*(in-5))* &
            (scratchc(i,j,k) + &
            (scratch13(i,j,k)+scratch13(i,j,k+1)+scratch13(i+1,j,k)+scratch13(i+1,j,k+1))*0.25d0 + &
            (scratch23(i,j,k)+scratch23(i,j,k+1)+scratch23(i,j+1,k)+scratch23(i,j+1,k+1))*0.25d0 + &
            (scratch12(i,j,k)+scratch12(i+1,j,k)+scratch12(i,j+1,k)+scratch12(i+1,j+1,k))*0.25d0)
    enddo; enddo; enddo

    ! advection for dels
    nlf=nlf_s
    alpha=alpha_s
    do k=ks,ke; do j=js,je
       qint(ism2:iep2) = s(ism2:iep2,j,k) + s0(ism2+myid1*(in-5):iep2+myid1*(in-5)) &
         + sr0(ism2+myid1*(in-5):iep2+myid1*(in-5))
       call xtvd(in, ijkn, dx1a, dx1bi, qint, dq, qlr, dqint)
       do i = is,iep1
          cmax = cmax1(i,j,k)
          rr = abs(dq(i))/(abs(dqint(i)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          flux(i) = d0m(i)*temp0m(i)*&
               (v1(i,j,k)*qlr(i) - 0.5d0*(abs(v1(i,j,k)) + cmax*qq) &
               *dq(i))*g2a(i)*g31a(i)
       enddo
       do i=is,ie
          dels(i,j,k) = dels(i,j,k) - (flux(i+1) - flux(i))*dx1ai(i)*g2bi(i)*g31bi(i)/ &
               d0(i+myid1*(in-5))/temp0(i+myid1*(in-5)) + &
               1.d0/temp0(i+myid1*(in-5))*qint(i)*0.5d0*(v1(i,j,k) + v1(i+1,j,k))* &
               ( (temp0(i+1+myid1*(in-5))-temp0(i+myid1*(in-5))) &
               *dx1bi(i+1) + &
               (temp0(i+myid1*(in-5))-temp0(i-1+myid1*(in-5))) &
               *dx1bi(i) )*0.5D0
       enddo
    enddo; enddo
    do k=ks,ke; do i=is,ie
       qint(jsm2:jep2) = s(i,jsm2:jep2,k) + s0(i+myid1*(in-5)) &
         + sr0(i+myid1*(in-5))
       call xtvd(jn, ijkn, dx2a, dx2bi, qint, dq, qlr, dqint)
       do j=js,jep1
          cmax = cmax2(i,j,k)
          rr = abs(dq(j))/(abs(dqint(j)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          flux(j) = (v2(i,j,k)*qlr(j) - 0.5d0*(abs(v2(i,j,k)) + cmax*qq)*dq(j))*g32a(j)
       enddo
       do j=js,je
          dels(i,j,k) = dels(i,j,k) - (flux(j+1) - flux(j))*dx2ai(j)*g2bi(i)*g32bi(j)
       enddo
    enddo; enddo
    do j=js,je; do i=is,ie
       qint(ksm2:kep2) = s(i,j,ksm2:kep2) + s0(i+myid1*(in-5)) &
         + sr0(i+myid1*(in-5))
       call xtvd(kn, ijkn, dx3a, dx3bi, qint, dq, qlr, dqint)
       do k=ks,kep1
          cmax = cmax3(i,j,k)
          rr = abs(dq(k))/(abs(dqint(k)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          flux(k) = v3(i,j,k)*qlr(k) - 0.5d0*(abs(v3(i,j,k)) + cmax*qq)*dq(k)
       enddo
       do k=ks,ke
          dels(i,j,k) = dels(i,j,k) - (flux(k+1) - flux(k))*dx3ai(k)*g31bi(i)*g32bi(j)
       enddo
    enddo; enddo

  end subroutine advection

  !=====================================================================================

  subroutine viscous
    use ModPar,      ONLY: inmax, in, jn, kn, myid1
    use ModGrid
    use ModBack,     ONLY: d0, temp0, ovrevar
    use ModField,    ONLY: v1, v2, v3
    use ModDel,      ONLY: dels1, dels2, dels3, dels
    implicit none

    integer :: i, j, k

    real :: d0ovrevar(inmax), ovrevarovtemp0(inmax)
    real, allocatable, dimension(:,:,:) :: rdivv_flx, scratch1, scratch2, scratch3
    !-------------------------------------------------------------------------------

    allocate(rdivv_flx(in,jn,kn))
    allocate(scratch1(in,jn,kn))
    allocate(scratch2(in,jn,kn))
    allocate(scratch3(in,jn,kn))
    ! define d0ovrevar and ovrevarovtemp0                                                 
    d0ovrevar = d0*ovrevar
    ovrevarovtemp0 = ovrevar/temp0

    ! terms involving div(v)                                                              
    do k=ksm1,kep1; do j=jsm1,jep1; do i=ism1,iep1
       rdivv_flx(i,j,k) = (g2a(i+1)*g31a(i+1)*v1(i+1,j,k) &
            - g2a(i)*g31a(i)*v1(i,j,k))*dx1ai(i)*g2bi(i)*g31bi(i) &
            + (g32a(j+1)*v2(i,j+1,k) - g32a(j)*v2(i,j,k) ) &
            *dx2ai(j)*g2bi(i)*g32bi(j) + (v3(i,j,k+1) - v3(i,j,k)) &
            *dx3ai(k)*g31bi(i)*g32bi(j)
       rdivv_flx(i,j,k) = rdivv_flx(i,j,k)*d0ovrevar(i+myid1*(in-5))
    enddo; enddo; enddo
    do k=ks,ke; do j=js,je; do i=is,iep1
       dels1(i,j,k) = dels1(i,j,k) - 2.d0/3.d0*(rdivv_flx(i,j,k) &
            - rdivv_flx(i-1,j,k))*dx1bi(i)
    enddo; enddo; enddo
    do k=ks,ke; do j=js,jep1; do i=is,ie
       dels2(i,j,k) = dels2(i,j,k) - 2.d0/3.d0*(rdivv_flx(i,j,k) &
            - rdivv_flx(i,j-1,k))*g2bi(i)*dx2bi(j)
    enddo; enddo; enddo
    do k=ks,kep1; do j=js,je; do i=is,ie
       dels3(i,j,k) = dels3(i,j,k) - 2.d0/3.d0*(rdivv_flx(i,j,k) &
            - rdivv_flx(i,j,k-1))*g31bi(i)*g32bi(j)*dx3bi(k)
    enddo; enddo; enddo
    ! s11
    do k=ksm1,kep1; do j=jsm1,jep1; do i=ism1,iep1
       scratch1(i,j,k) = 2.d0*(v1(i+1,j,k) - v1(i,j,k))*dx1ai(i)
       rdivv_flx(i,j,k) = scratch1(i,j,k)*d0ovrevar(i+myid1*(in-5))
    enddo; enddo; enddo
    do k=ks,ke; do j=js,je; do i=is,iep1
       dels1(i,j,k) = dels1(i,j,k) + (g2b(i)*g31b(i)*rdivv_flx(i,j,k) &
            - g2b(i-1)*g31b(i-1)*rdivv_flx(i-1,j,k)) &
            *dx1bi(i)*g2ai(i)*g31ai(i)
    enddo; enddo; enddo
    ! s12                                                                                 
    do k=ks,ke; do j=js,jep1; do i=is,iep1
       scratch2(i,j,k) = (g2ai(i)*(v1(i,j,k)-v1(i,j-1,k))*dx2bi(j) &
            +(g2bi(i)*v2(i,j,k)-g2bi(i-1)*v2(i-1,j,k))*dx1bi(i)*g2a(i))
       rdivv_flx(i,j,k) = scratch2(i,j,k)*0.5d0* &
            (d0ovrevar(i-1+myid1*(in-5)) + d0ovrevar(i+myid1*(in-5)))
    enddo; enddo; enddo
    do k=ks,ke; do j=js,je; do i=is,iep1
       dels1(i,j,k) = dels1(i,j,k) + (g32a(j+1)*rdivv_flx(i,j+1,k) - &
            g32a(j)*rdivv_flx(i,j,k))*dx2ai(j)*g2ai(i)*g32bi(j)
    enddo; enddo; enddo
    do k=ks,ke; do j=js,jep1; do i=is,ie
       dels2(i,j,k) = dels2(i,j,k) + (g2a(i+1)*g31a(i+1)*rdivv_flx(i+1,j,k) &
            - g2a(i)*g31a(i)*rdivv_flx(i,j,k))*dx1ai(i)*g2bi(i)*g31bi(i) &
            + 0.5d0*(rdivv_flx(i,j,k)+rdivv_flx(i+1,j,k))*g2bi(i)*dg2ad1(i)
    enddo; enddo; enddo
    ! s13                                                                                 
    do k=ks,kep1; do j=js,je; do i=is,iep1
       scratch3(i,j,k) = ((g31bi(i)*v3(i,j,k) - g31bi(i-1)*v3(i-1,j,k)) &
            *dx1bi(i)*g31a(i) + g31ai(i)*g32bi(j)* &
            (v1(i,j,k) - v1(i,j,k-1))*dx3bi(k))
       rdivv_flx(i,j,k) = scratch3(i,j,k)*0.5d0* &
            (d0ovrevar(i-1+myid1*(in-5)) + d0ovrevar(i+myid1*(in-5)))
    enddo; enddo; enddo
    do k=ks,ke; do j=js,je; do i=is,iep1
       dels1(i,j,k) = dels1(i,j,k) + (rdivv_flx(i,j,k+1) - rdivv_flx (i,j,k)) &
            *dx3ai(k)*g31ai(i)*g32bi(j)
    enddo; enddo; enddo
    do k=ks,kep1; do j=js,je; do i=is,ie
       dels3(i,j,k) = dels3(i,j,k) + (g2a(i+1)*g31a(i+1) &
            *g31a(i+1)*g32b(j)*rdivv_flx(i+1,j,k) &
            - g2a(i)*g31a(i)*g31a(i)*g32b(j)*rdivv_flx(i,j,k)) &
            *dx1ai(i)*g2bi(i)*g31bi(i)*g31bi(i)*g32bi(j)
    enddo; enddo; enddo
    do k=ks,ke; do j=js,je; do i=is,ie
       dels(i,j,k) = dels(i,j,k) + ovrevarovtemp0(i+myid1*(in-5)) &
            *((0.25d0*(scratch2(i,j,k) + scratch2(i+1,j,k) &
            + scratch2(i,j+1,k) + scratch2(i+1,j+1,k)))**2 &
            +(0.25*(scratch3(i,j,k) + scratch3(i,j,k+1) &
            +scratch3(i+1,j,k) + scratch3(i+1,j,k+1)))**2 )
    enddo; enddo; enddo
    ! s22                                                                                 
    do k=ksm1,kep1; do j=jsm1,jep1; do i=ism1,iep1
       scratch2(i,j,k)=(2.d0*g2bi(i)*(v2(i,j+1,k)-v2(i,j,k))*dx2ai(j) &
            +(v1(i,j,k)+v1(i+1,j,k))*g2bi(i)*dg2ad1(i))
       rdivv_flx(i,j,k)=scratch2(i,j,k)*d0ovrevar(i+myid1*(in-5))
    enddo; enddo; enddo
    do k=ks,ke; do j=js,jep1; do i=is,ie
       dels2(i,j,k) = dels2(i,j,k) &
            + (g32b(j)*rdivv_flx(i,j,k)-g32b(j-1)*rdivv_flx(i,j-1,k)) &
            * dx2bi(j) * g2bi(i) * g32ai(j)
    enddo; enddo; enddo
    do k=ks,ke; do j=js,je; do i=is,iep1
       dels1(i,j,k) = dels1(i,j,k) - 0.5d0*(rdivv_flx(i,j,k)+rdivv_flx(i-1,j,k)) &
            * g2ai(i) * dg2bd1(i)
    enddo; enddo; enddo
    do k=ks,ke; do j=js,je; do i=is,ie
       dels(i,j,k) = dels(i,j,k) + ovrevarovtemp0(i+myid1*(in-5))/6.d0 &
            * (scratch1(i,j,k)-scratch2(i,j,k))**2
    enddo; enddo; enddo
    ! s23                                                                                 
    do k=ks,kep1; do j=js,jep1; do i=is,ie
       scratch3(i,j,k) = (g31bi(i)*g32ai(j)*(v2(i,j,k) - v2(i,j,k-1))*dx3bi(k) + &
            (g32bi(j)*v3(i,j,k) - g32bi(j-1)*v3(i,j-1,k))*dx2bi(j)*g2bi(i)*g32a(j))
       rdivv_flx(i,j,k) = scratch3(i,j,k)*d0ovrevar(i+myid1*(in-5))
    enddo; enddo; enddo
    do k=ks,ke; do j=js,jep1; do i=is,ie
       dels2(i,j,k) = dels2(i,j,k) + (rdivv_flx(i,j,k+1) - rdivv_flx(i,j,k)) &
            * dx3ai(k) * g31bi(i) * g32ai(j)
    enddo; enddo; enddo
    do k=ks,kep1; do j=js,je; do i=is,ie
       dels3(i,j,k) = dels3(i,j,k) + (g32a(j+1)*g31b(i)*g32a(j+1)*rdivv_flx(i,j+1,k) &
            - g32a(j)*g31b(i)*g32a(j)*rdivv_flx(i,j,k)) &
            * g2bi(i) * g32bi(j) * dx2ai(j) * g31bi(i) * g32bi(j)
    enddo; enddo; enddo
    do k=ks,ke; do j=js,je; do i=is,ie
       dels(i,j,k) = dels(i,j,k) + ovrevarovtemp0(i+myid1*(in-5)) &
            *(0.25d0*(scratch3(i,j,k) + scratch3(i,j+1,k) + &
            scratch3(i,j,k+1) + scratch3(i,j+1,k+1)))**2
    enddo; enddo; enddo
    ! s33                                                                                 
    do k=ksm1,kep1; do j=jsm1,jep1; do i=ism1,iep1
       scratch3(i,j,k) = (2.d0*g31bi(i)*g32bi(j)*(v3(i,j,k+1) - v3(i,j,k))*dx3ai(k) + &
            (v1(i,j,k) + v1(i+1,j,k))*g31bi(i)*dg31ad1(i) + &
            (v2(i,j,k) + v2(i,j+1,k))*g2bi(i)*g32bi(j)*dg32ad2(j))
       rdivv_flx(i,j,k) = scratch3(i,j,k)*d0ovrevar(i+myid1*(in-5))
    enddo; enddo; enddo
    do k=ks,ke; do j=js,je; do i=is,iep1
       dels1(i,j,k) = dels1(i,j,k) &
            - 0.5d0*(rdivv_flx(i,j,k)+rdivv_flx(i-1,j,k)) * g31ai(i) * dg31bd1(i)
    enddo; enddo; enddo
    do k=ks,ke; do j=js,jep1; do i=is,ie
       dels2(i,j,k) = dels2(i,j,k) - &
            0.5d0*(rdivv_flx(i,j,k)+rdivv_flx(i,j-1,k))*g2bi(i)*g32ai(j)*dg32bd2(j)
    enddo; enddo; enddo
    do k=ks,kep1; do j=js,je; do i=is,ie
       dels3(i,j,k) = dels3(i,j,k) + &
            (rdivv_flx(i,j,k) - rdivv_flx(i,j,k-1))*dx3bi(k)*g31bi(i)*g32bi(j)
    enddo; enddo; enddo
    do k=ks,ke; do j=js,je; do i=is,ie
       dels(i,j,k) = dels(i,j,k) + &
            ovrevarovtemp0(i+myid1*(in-5))/6.d0*((scratch2(i,j,k) - scratch3(i,j,k))**2 &
            + (scratch3(i,j,k) - scratch1(i,j,k))**2)
    enddo; enddo; enddo

    deallocate(rdivv_flx, scratch1, scratch2, scratch3)

  end subroutine viscous

  !=====================================================================================

  subroutine radheat
    use ModPar,   ONLY: myid1, in
    use ModGrid,  ONLY: is, ie, js, je, ks, ke
    use ModBack,  ONLY: heatrad, d0, temp0
    use ModDel,   ONLY: dels
    implicit none

    integer :: i
    !---------------------------------------------------------------------------------    

    do i = is, ie
       dels(i,js:je,ks:ke) = dels(i,js:je,ks:ke) + heatrad(i+myid1*(in-5)) &
            /d0(i+myid1*(in-5))/temp0(i+myid1*(in-5))
    enddo

  end subroutine radheat

  !=====================================================================================
  
  subroutine hsmoc
    use ModPar
    use ModSundry,   ONLY: dt
    use ModGrid
    use ModBack,     ONLY: d0, ovrmvar, temp0
    use ModField,    ONLY: b1, b2, b3, v1, v2, v3
    use ModDel,      ONLY: emf1s, emf2s, emf3s, dels
    use ModBval
    use ModInterp,   ONLY: xint1d, xzc1d,xint1d1
    use ModWork
    use ModSundry,  ONLY: dtmdi2
    use ModBval,    ONLY: bvaleta
    use ModInterp,  ONLY: xdel
    use ModMpi
    use ModFSAM,    ONLY: iComm
    use ModBoundary

    implicit none

    integer :: i, j, k
    real    :: absb, sgnp, sgnm, q2, srdi, srd, aave, gfact, vsnm1, bsnm1

    real, dimension(in)   :: srd1, srd2, srd3
    real, dimension(ijkn) :: bt, vt, vfl, bint, vint, bave, &
         vave, vchp, vchm, vtmp, btmp, vpch, vmch, bpch, bmch, &
         srdp, srdm, srdpi, srdmi
    real, dimension(jn,kn) :: v3intj, b3intj, v2intk, b2intk
    real, dimension(kn,in) :: v1intk, b1intk, v3inti, b3inti
    real, dimension(in,jn) :: v2inti, b2inti, v1intj, b1intj

    integer :: ierr
    real :: dr1, dr2, dr3, drmin, current, dtmdi2tmp, dtmdi2lc
    real :: dbint(ijkn), dbint1
    real, dimension(:,:,:),allocatable :: emf1, emf2, emf3, eta
    real :: bsnp
    real, dimension(:,:,:),allocatable :: curl1, curl2, curl3, bsnp1

    ! paramters for testing, write emf1-3s if DoWriteEmf
    character(len=4) :: idcpu
    logical :: DoWriteEmf = .false., DoWriteDels=.false.
    !-----------------------------------------------------------------------------------
    do i = ism1,iep1
       srd1(i) = sqrt(0.5d0*(d0(i+myid1*(in-5)) + d0(i-1+myid1*(in-5))))
       srd2(i) = sqrt(d0(i+myid1*(in-5)))
       srd3(i) = srd2(i)
    enddo

    allocate(emf1(in,jn,kn), emf2(in,jn,kn), emf3(in,jn,kn))
    allocate(eta(in,jn,kn), bsnp1(in,jn,kn))
    allocate(curl1(in,jn,kn),curl2(in,jn,kn),curl3(in,jn,kn))

    emf1 = 0.d0
    curl1 = 0.d0
    ! emf1
    do i=is,ie
       ! Compute upwinded b3 and v3 in the 2-direction and use them to
       ! compute the wave speeds for the characteristic cones for the MoC
       do k=ks,kep1
          vt(jsm2:jep2) = v3(i,jsm2:jep2,k)
          bt(jsm2:jep2) = b3(i,jsm2:jep2,k)
          do j=jsm2,jep2
             vfl(j) = 0.5d0*(v2(i,j,k) + v2(i,j,k-1))
          enddo
          if(x2a(js) .eq. 0.d0) then
             vt(jsm1) = - vt(js)
             vt(jsm2) = - vt(jsp1)
             bt(jsm1) = - bt(js)
             bt(jsm2) = - bt(jsp1)
          endif
          if(x2a(jep1) .eq. pi) then
             vt(jep1) = - vt(je)
             vt(jep2) = - vt(jem1)
             bt(jep1) = - bt(je)
             bt(jep2) = - bt(jem1)
          endif
          gfact = dt*g2bi(i)
          call xint1d(jn, bt, vt, vfl, gfact, dx2bi, dx2a, bint, vint)
          v3intj(js:jep1,k) = vint(js:jep1)
          b3intj(js:jep1,k) = bint(js:jep1)
          ! diffemf
          bint(jsm2:jep2) = bt(jsm2:jep2)*g32b(jsm2:jep2)
          call xdel(jn, bint, dbint, dx2a, dx2bi)
          do j=js,jep1
             emf1(i,j,k) = - dbint(j)*g2bi(i)*dx2bi(j)*g32ai(j)
             curl1(i,j,k) =  (bint(j) - bint(j-1))*g2bi(i)*dx2bi(j)*g32ai(j)
          enddo
       enddo
       ! Select an effective density and determine the characteristic velocity
       ! using upwinded values for each characteristic in the 3-direction
       do j=js,jep1
          vave(ks:kep1) = v3intj(j,ks:kep1)
          bave(ks:kep1) = b3intj(j,ks:kep1)
          srd = srd2(i)
          srdi = 1.d0/srd
          do k=ks,kep1
             absb    = abs(bave(k))
             vchp(k) = vave(k) - absb*srdi
             vchm(k) = vave(k) + absb*srdi
          enddo
          vtmp(ksm2:kep2) = v2(i,j,ksm2:kep2)
          btmp(ksm2:kep2) = b2(i,j,ksm2:kep2)
          gfact = dt*g31bi(i)*g32ai(j)
          call xzc1d(kn, btmp, vtmp, vchp, vchm, gfact, dx3bi, dx3a, &
               bpch, bmch, vpch, vmch)
          do k=ks,kep1
             q2    = sign(1.d0, bave(k))
             emf1s(i,j,k) = 0.5d0*((vpch(k)*srd + vmch(k)*srd + q2*(bpch(k) - bmch(k))) &
                  /(2.d0*srd)*bave(k) - &
                  ((bpch(k)*srdi + bmch (k)*srdi + q2*(vpch(k) - vmch(k))) &
                  /(2.d0*srdi)*vave(k)))
          enddo
       enddo
       ! Compute upwinded b2 and v2 in the 3-direction and use them to  
       ! compute the wave speeds for the chracteristic cones for the MoC
       do j=js,jep1
          vt(ksm2:kep2) = v2(i,j,ksm2:kep2)
          bt(ksm2:kep2) = b2(i,j,ksm2:kep2)
          do  k=ksm2,kep2
             vfl(k) = 0.5d0*(v3(i,j,k) + v3(i,j-1,k))
          enddo
          gfact = dt*g31bi(i)*g32ai(j)
          !call xint1d(kn, bt, vt, vfl, gfact, dx3bi, dx3a, bint, vint)
          call xint1d1(kn, bt, vt, vfl, gfact, dx3bi, dx3a, bint, vint,dbint)
          v2intk(j,ks:kep1) = vint(ks:kep1)
          b2intk(j,ks:kep1) = bint(ks:kep1)
          ! diffemf
          !bint(ksm2:kep2) = bt(ksm2:kep2)
          !call xdel(kn, bint, dbint, dx3a, dx3bi)
          do k=ks,kep1
             emf1(i,j,k) = emf1(i,j,k) + &
                  dbint(k)*g31bi(i)*g32ai(j)*dx3bi(k)
             curl1(i,j,k) = curl1(i,j,k) - (bt(k) - bt(k-1))*g31bi(i)*g32ai(j)*dx3bi(k)
          enddo

       enddo
       ! Select an effective density and determine the characteristic velocity
       ! using upwinded values for each characteristic in the 2-direction
       srd = srd3(i)
       srdi = 1.d0/srd
       do  k=ks,kep1
          vave(js:jep1) = v2intk(js:jep1,k)
          bave(js:jep1) = b2intk(js:jep1,k)
          do j=js,jep1
             absb     = abs(bave(j))
             vchp(j)  = vave(j) - absb*srdi
             vchm(j)  = vave(j) + absb*srdi
          enddo
          vtmp(jsm2:jep2) = v3(i,jsm2:jep2,k)
          btmp(jsm2:jep2) = b3(i,jsm2:jep2,k)
          if (x2a(js) .eq. 0.0) then
             vtmp(jsm1) = - vtmp(js)
             vtmp(jsm2) = - vtmp(jsp1)
             btmp(jsm1) = - btmp(js)
             btmp(jsm2) = - btmp(jsp1)
          endif
          if (x2a(jep1) .eq. pi) then
             vtmp(jep1) = - vtmp(je)
             vtmp(jep2) = - vtmp(jem1)
             btmp(jep1) = - btmp(je)
             btmp(jep2) = - btmp(jem1)
          endif
          gfact = dt*g2bi(i)
          call xzc1d(jn, btmp, vtmp, vchp, vchm, gfact, dx2bi, dx2a, &
               bpch, bmch, vpch, vmch)
          do j=js,jep1
             q2    = sign(1.d0, bave(j))
             vsnm1 = (vpch(j)*srd + vmch(j)*srd + q2*(bpch(j) - bmch(j))) &
                  /(2.d0*srd)*bave(j)
             bsnm1 = (bpch(j)*srdi + bmch(j)*srdi + q2*(vpch(j) - vmch(j))) &
                  /(2.d0*srdi)*vave(j)
             emf1s(i,j,k) = emf1s(i,j,k) + 0.5d0*(bsnm1-vsnm1)
             numemf1(i,j,k) = emf1s(i,j,k) &
                  - 0.25D0*(v2(i,j,k)+v2(i,j,k-1))*(b3(i,j,k)+b3(i,j-1,k)) &
                  + 0.25D0*(b2(i,j,k)+b2(i,j,k-1))*(v3(i,j,k)+v3(i,j-1,k))
          enddo
       enddo
    enddo


    emf2 = 0.d0
    curl2 = 0.d0
    ! emf2
    do j=js,je
       ! Compute upwinded b1 and v1 in the 3-direction and use them to
       ! compute the wave speeds for the characteristic cones for the MoC
       do i=is,iep1
          vt(ksm2:kep2) = v1(i,j,ksm2:kep2)
          bt(ksm2:kep2) = b1(i,j,ksm2:kep2)
          do k=ksm2,kep2
             vfl(k) = 0.5*(v3(i,j,k) + v3(i-1,j,k))
          enddo
          gfact = dt*g31ai(i)*g32bi(j)
          !call xint1d(kn, bt, vt, vfl, gfact, dx3bi, dx3a, bint, vint)
          call xint1d1(kn, bt, vt, vfl, gfact, dx3bi, dx3a, bint, vint,dbint)
          v1intk(ks:kep1,i) = vint(ks:kep1)
          b1intk(ks:kep1,i) = bint(ks:kep1)
          !diffemf
          !bint(ksm2:kep2) = bt(ksm2:kep2)
          !call xdel(kn, bint, dbint, dx3a, dx3bi)
          do k=ks,kep1
             emf2(i,j,k) = - dbint(k)*g31ai(i)*g32bi(j)*dx3bi(k)
             curl2(i,j,k) = (bt(k) - bt(k-1))*g31ai(i)*g32bi(j)*dx3bi(k)
          enddo
       enddo
       ! Select an effective density and determine the characteristic velocity  
       ! using upwinded values for each characteristic in the 1-direction
       do k=ks,kep1
          vave(is:iep1) = v1intk(k,is:iep1)
          bave(is:iep1) = b1intk(k,is:iep1)
          do i=is,iep1
             absb    = abs(bave(i))
             aave    = 0.5d0*absb*(srd3(i) + srd3(i-1))/(srd3(i)*srd3(i-1))
             sgnp    = sign(0.5d0, vave(i) - aave)
             sgnm    = sign(0.5d0, vave(i) + aave)
             srdp(i) = (0.5d0 + sgnp)*srd3(i-1) + (0.5d0 - sgnp)*srd3(i)
             srdm(i) = (0.5d0 + sgnm)*srd3(i-1) + (0.5d0 - sgnm)*srd3(i)
             srdpi(i) = 1.d0/srdp(i)
             srdmi(i) = 1.d0/srdm(i)
             vchp(i) = vave(i) - absb*srdpi(i)
             vchm(i) = vave(i) + absb*srdmi(i)
          enddo
          vtmp(ism2:iep2) = v3(ism2:iep2,j,k)
          btmp(ism2:iep2) = b3(ism2:iep2,j,k)
          gfact = dt
          call xzc1d(in, btmp, vtmp, vchp, vchm, gfact, dx1bi, dx1a, &
               bpch, bmch, vpch, vmch, DoBc = .true.)
          do  i=is,iep1
             q2    = sign(1.d0, bave(i))
             emf2s(i,j,k) = 0.5d0*((vpch(i)*srdp(i) + vmch(i)*srdm (i) + &
                  q2*(bpch(i) - bmch(i)))/(srdp(i) + srdm(i))*bave(i) - &
                  ((bpch(i)*srdpi(i) + bmch (i)*srdmi(i) + &
                  q2*(vpch(i) - vmch(i)))/(srdpi(i) + srdmi(i))*vave(i)))
          enddo
       enddo
       ! Compute upwinded b3 and v3 in the 1-direction and use them to     
       ! compute the wave speeds for the chracteristic cones for the MoC
       do k=ks,kep1
          vt(ism2:iep2) = v3(ism2:iep2,j,k)
          bt(ism2:iep2) = b3(ism2:iep2,j,k)
          do i=ism2,iep2
             vfl(i) = 0.5d0*(v1(i,j,k) + v1(i,j,k-1))
          enddo
          gfact = dt
          call xint1d(in, bt, vt, vfl, gfact, dx1bi, dx1a, bint, vint)
          v3inti(k,is:iep1) = vint(is:iep1)
          b3inti(k,is:iep1) = bint(is:iep1)
          !diffemf
          bint(ism2:iep2) = bt(ism2:iep2)*g31b(ism2:iep2)
          call xdel(in, bint, dbint, dx1a, dx1bi)
          do i=is,iep1
             curl2(i,j,k) = curl2(i,j,k) - (bint(i) - bint(i-1))*g31ai(i)*dx1bi(i)
             emf2(i,j,k) = emf2(i,j,k) + dbint(i)*g31ai(i)*dx1bi(i)
          enddo
       enddo
       ! Select an effective density and determine the characteristic velocity 
       ! using upwinded values for each characteristic in the 3-direction.
       do i=is,iep1
          srd = srd1(i)
          srdi = 1.d0/srd
          vave(ks:kep1) = v3inti(ks:kep1,i)
          bave(ks:kep1) = b3inti(ks:kep1,i)
          do k=ks,kep1
             absb     = abs(bave(k))
             vchp(k) = vave(k) - absb*srdi
             vchm(k) = vave(k) + absb*srdi
          enddo
          vtmp(ksm2:kep2) = v1(i,j,ksm2:kep2)
          btmp(ksm2:kep2) = b1(i,j,ksm2:kep2)
          gfact = dt*g31ai(i)*g32bi(j)
          call xzc1d(kn, btmp, vtmp, vchp, vchm, gfact, dx3bi, dx3a, &
               bpch, bmch, vpch, vmch)
          do k=ks,kep1
             q2    = sign(1.d0, bave(k))
             vsnm1 = (vpch(k)*srd + vmch(k)*srd + q2*(bpch(k) - bmch(k))) &
                  /(2.d0*srd)*bave(k)
             bsnm1 = (bpch(k)*srdi + bmch(k)*srdi + q2*(vpch(k) - vmch(k))) &
                  /(2.d0*srdi)*vave(k)
             emf2s(i,j,k) = emf2s(i,j,k) + 0.5d0*(bsnm1-vsnm1)
             numemf2(i,j,k) = emf2s(i,j,k) &
                  - 0.25D0*(v3(i,j,k)+v3(i-1,j,k))*(b1(i,j,k)+b1(i,j,k-1)) &
                  + 0.25D0*(b3(i,j,k)+b3(i-1,j,k))*(v1(i,j,k)+v1(i,j,k-1))
          enddo
       enddo
    enddo


    emf3 = 0.d0
    curl3 = 0.d0
    ! emf3
    do k=ks,ke
       ! Compute upwinded b2 and v2 in the 1-direction and use them to  
       ! compute the wave speeds for the chracteristic cones for MoC
       do j=js,jep1
          vt(ism2:iep2) = v2(ism2:iep2,j,k)
          bt(ism2:iep2) = b2(ism2:iep2,j,k)
          do i=ism2,iep2
             vfl(i) = 0.5d0*(v1(i,j,k) + v1(i,j-1,k))
          enddo
          gfact = dt
          call xint1d(in, bt, vt, vfl, gfact, dx1bi, dx1a, bint, vint)
          v2inti(is:iep1,j) = vint(is:iep1)
          b2inti(is:iep1,j) = bint(is:iep1)
          !diffemf
          bint(ism2:iep2) = bt(ism2:iep2)*g2b(ism2:iep2)
          call xdel(in, bint, dbint, dx1a, dx1bi)
          do i=is,iep1
             emf3(i,j,k) = - dbint(i)*dx1bi(i)*g2ai(i)
             curl3(i,j,k) = (bint(i) - bint(i-1))*dx1bi(i)*g2ai(i)
          enddo
       enddo
       ! Select an effective density and determine the characteristic velocity 
       ! using upwinded values for each characteristic in the 2-direction
       do i=is,iep1
          srd  = srd1(i)
          srdi = 1.d0/srd
          vave(js:jep1) = v2inti(i,js:jep1)
          bave(js:jep1) = b2inti(i,js:jep1)
          do j=js,jep1
             absb     = abs(bave(j))
             vchp(j) = vave(j) - absb*srdi
             vchm(j) = vave(j) + absb*srdi
          enddo
          vtmp(jsm2:jep2) = v1(i,jsm2:jep2,k)
          btmp(jsm2:jep2) = b1(i,jsm2:jep2,k)
          gfact = dt*g2ai(i)
          call xzc1d(jn, btmp, vtmp, vchp, vchm, gfact, dx2bi, dx2a, &
               bpch, bmch, vpch, vmch)
          do j=js,jep1
             q2    = sign(1.d0, bave(j))
             emf3s(i,j,k) = 0.5d0*((vpch(j)*srd + vmch(j)*srd + q2*(bpch(j) - bmch(j))) &
                  /(2.d0*srd)*bave(j) - &
                  ((bpch(j)*srdi + bmch(j)*srdi + q2*(vpch(j) - vmch(j))) &
                  /(2.d0*srdi)*vave(j)))
          enddo
       enddo
       ! Select an effective density and determine the characteristic velocity
       ! using upwinded values for each characteristic in the 1-direction
       do i=is,iep1
          vt(jsm2:jep2) = v1(i,jsm2:jep2,k)
          bt(jsm2:jep2) = b1(i,jsm2:jep2,k)
          do j=jsm2,jep2
             vfl(j) = 0.5*(v2(i,j,k) + v2(i-1,j,k))
          enddo
          gfact = dt*g2ai(i)
          !call xint1d(jn, bt, vt, vfl, gfact, dx2bi, dx2a, bint, vint)
          call xint1d1(jn, bt, vt, vfl, gfact, dx2bi, dx2a, bint, vint,dbint)
          b1intj(i,js:jep1) = bint(js:jep1)
          v1intj(i,js:jep1) = vint(js:jep1)
          !diffemf
          !bint(jsm2:jep2) = bt(jsm2:jep2)
          !call xdel(jn, bint, dbint, dx2a, dx2bi)
          do j=js,jep1
             emf3(i,j,k) = emf3(i,j,k) + dbint(j)*g2ai(i)*dx2bi(j)
             curl3(i,j,k) = curl3(i,j,k)-(bt(j)-bt(j-1))*g2ai(i)*dx2bi(j)
          enddo
       enddo
       ! Select an effective density and determine the characteristic velocity 
       ! using upwinded values for each characteristic in the 1-direction
       do j=js,jep1
          do i=is,iep1
             vave(i) = v1intj(i,j)
             bave(i) = b1intj(i,j)
             absb    = abs(bave(i))
             aave    = 0.5d0*absb*(srd2(i) + srd2(i-1))/(srd2(i)*srd2(i-1))
             sgnp    = sign(0.5d0, vave(i) - aave)
             sgnm    = sign(0.5d0, vave(i) + aave)
             srdp(i) = (0.5d0 + sgnp)*srd2(i-1) + (0.5d0 - sgnp)*srd2(i)
             srdm(i) = (0.5d0 + sgnm)*srd2(i-1) + (0.5d0 - sgnm)*srd2(i)
             srdpi(i) = 1.d0/srdp(i)
             srdmi(i) = 1.d0/srdm(i)
             vchp(i) = vave(i) - absb*srdpi(i)
             vchm(i) = vave(i) + absb*srdmi(i)
          enddo
          vtmp(ism2:iep2) = v2(ism2:iep2,j,k)
          btmp(ism2:iep2) = b2(ism2:iep2,j,k)
          gfact = dt
          call xzc1d(in, btmp, vtmp, vchp, vchm, gfact, dx1bi, dx1a, &
               bpch, bmch, vpch, vmch, DoBc = .true.)
          do  i=is,iep1
             q2     = sign(1.d0, bave(i))
             vsnm1 = (vpch(i)*srdp(i) + vmch(i)*srdm(i) + q2*(bpch(i) - bmch(i))) &
                  /(srdp(i) + srdm(i))*bave(i)
             bsnm1 = (bpch(i)*srdpi(i) + bmch(i)*srdmi(i) + q2*(vpch(i) - vmch(i))) &
                  /(srdpi(i) + srdmi(i))*vave(i)
             emf3s(i,j,k) = emf3s(i,j,k) + 0.5d0*(bsnm1-vsnm1)
             numemf3(i,j,k) = emf3s(i,j,k) &
                  - 0.25D0*(v1(i,j,k)+v1(i,j-1,k))*(b2(i,j,k)+b2(i-1,j,k)) &
                  + 0.25D0*(b1(i,j,k)+b1(i,j-1,k))*(v2(i,j,k)+v2(i-1,j,k))
          enddo
       enddo
    enddo

    if (DoWriteEmf) then
       write(idcpu,'(i4.4)') myid
       open(unit=17,file='emf.cpu'//idcpu//'.hsmoc',form='unformatted', &
            access='stream')
       write(17) in, jn, kn
       write(17) (((emf1s(i,j,k),i=ism2,iep2),j=jsm2,jep3),k=ksm2,kep3)
       write(17) (((emf2s(i,j,k),i=ism2,iep3),j=jsm2,jep2),k=ksm2,kep3)
       write(17) (((emf3s(i,j,k),i=ism2,iep3),j=jsm2,jep3),k=ksm2,kep2)
       close(17)
    endif

    ! calculate eta
    dtmdi2lc = 0.d0
    do k=ks,ke; do j=js,je; do i=is,ie
       current = sqrt(0.25d0*(emf1(i,j,k)**2 + emf1(i,j+1,k)**2 + &
            emf1(i,j,k+1)**2 + emf1(i,j+1,k+1)**2) + &
            0.25d0*(emf2(i,j,k)**2 + emf2(i+1,j,k)**2 + &
            emf2(i,j,k+1)**2 + emf2(i+1,j,k+1)**2) + &
            0.25d0*(emf3(i,j,k)**2 + emf3(i+1,j,k)**2 + &
            emf3(i,j+1,k)**2 + emf3(i+1,j+1,k)**2))
       dr1 = dx1a(i)
       dr2 = dx2a(j)*g2b(i)
       dr3 = dx3a(k)*g31b(i)*g32b(j)
       drmin = min(dr1, dr2, dr3)
       eta(i,j,k) = 0.5d0*drmin**2*current/sqrt(d0(i+myid1*(in-5)))
       dtmdi2tmp = (4.d0*eta(i,j,k)/drmin**2)**2
       dtmdi2lc = max(dtmdi2lc, dtmdi2tmp)
    enddo; end do; end do
    call bvaleta(eta)
    call MPI_ALLREDUCE(dtmdi2lc, dtmdi2, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
         iComm, ierr)

    ! emf1
    bsnp1 = 0.d0
    do k=ks,kep1; do j=js,jep1; do i=is,ie
       bsnp1(i,j,k) = 0.25d0*(eta(i,j,k) + eta(i,j-1,k) + eta(i,j-1,k-1) + eta(i,j,k-1))
    enddo; enddo; enddo
    emf1 = emf1*bsnp1
    emf1s = emf1s + emf1
    numemf1 = numemf1 + emf1
    if(myid2.eq.0) numemf1(:,js,:) = 0.d0
    if(myid2.eq.nproc2-1) numemf1(:,jep1,:) = 0.d0

    ! emf2
    bsnp1 = 0.d0
    do k=ks,kep1; do j=js,je; do i=is,iep1
       bsnp1(i,j,k) = 0.25d0*(eta(i,j,k) + eta(i-1,j,k) + eta(i-1,j,k-1) + eta(i,j,k-1))
    enddo; enddo; enddo
    emf2 = emf2*bsnp1
    emf2s = emf2s + emf2
    numemf2 = numemf2 + emf2
    if(myid1.eq.0) numemf2(is,:,:) = 0.d0
    if(myid1.eq.nproc1-1.and.noib.ne.2.and.noib.ne.5) numemf2(iep1,:,:)=0.d0
    
    ! emf3
    bsnp1 = 0.d0
    do k=ks,ke; do j=js,jep1; do i=is,iep1
       bsnp1(i,j,k) = 0.25d0*(eta(i,j,k) + eta(i-1,j,k) + eta(i-1,j-1,k) + eta(i,j-1,k))
    enddo; enddo; enddo
    emf3 = emf3*bsnp1
    emf3s = emf3s + emf3
    numemf3 = numemf3 + emf3
    if(myid1.eq.0) numemf3(is,:,:) = 0.d0
    if(myid1.eq.nproc1-1.and.noib.ne.2.and.noib.ne.5) numemf3(iep1,:,:) = 0.d0
    if(myid2.eq.0) numemf3(:,js,:) = 0.d0
    if(myid2.eq.nproc2-1) numemf3(:,jep1,:) = 0.d0

    do i=is,iep1
       emf1(i,:,:)  = curl1(i,:,:)*(-ovrmvar(i+myid1*(in-5)))
       emf2(i,:,:)  = curl2(i,:,:)*(- 0.5d0*(ovrmvar(i+myid1*(in-5)) & 
            + ovrmvar(i-1+myid1*(in-5))))
       emf3(i,:,:) = curl3(i,:,:)*(-0.5d0*(ovrmvar(i+myid1*(in-5)) &
            + ovrmvar(i-1+myid1*(in-5)))) 
    enddo

    emf1s = emf1s + emf1
    call bvalemf1(emf1s)
    if(myid2.eq.0) emf1(:,js,:) = 0.d0
    if(myid2.eq.nproc2-1) emf1(:,jep1,:) = 0.d0

    emf2s = emf2s + emf2
    call bvalemf2(emf2s)
    if(myid1.eq.0) emf2(is,:,:) = 0.d0
    if(myid1.eq.nproc1-1.and.noib.ne.2.and.noib.ne.5) emf2(iep1,:,:)=0.d0

    emf3s = emf3s + emf3
    call bvalemf3(emf3s)
    if(myid1.eq.0) emf3(is,:,:) = 0.d0
    if(myid1.eq.nproc1-1.and.noib.ne.2.and.noib.ne.5) emf3(iep1,:,:) = 0.d0
    if(myid2.eq.0) emf3(:,js,:) = 0.d0
    if(myid2.eq.nproc2-1) emf3(:,jep1,:) = 0.d0

    ! compute resistive heating 
    do k=ks,ke; do j=js,je; do i=is,ie
       dels(i,j,k) = dels(i,j,k) - &
            ((curl1(i,j,k) + curl1(i,j+1,k) + curl1(i,j,k+1) + curl1(i,j+1,k+1)) &
            *(numemf1(i,j,k) + numemf1(i,j+1,k) + numemf1(i,j,k+1) + numemf1(i,j+1,k+1)) + &
            (curl2(i,j,k) + curl2(i+1,j,k) + curl2(i,j,k+1) + curl2(i+1,j,k+1)) &
            *(numemf2(i,j,k) + numemf2(i+1,j,k) + numemf2(i,j,k+1) + numemf2(i+1,j,k+1)) + &
            (curl3(i,j,k) + curl3(i+1,j,k) + curl3(i,j+1,k) + curl3(i+1,j+1,k)) &
            *(numemf3(i,j,k) + numemf3(i+1,j,k) + numemf3(i,j+1,k) + numemf3(i+1,j+1,k))) &
            /16.d0/d0(i+myid1*(in-5)) &
            /temp0(i+myid1*(in-5)) + &
            !dels(i,j,k) = dels(i,j,k) + &
            ((emf1(i,j,k) + emf1(i,j+1,k) + emf1(i,j,k+1) + emf1(i,j+1,k+1))**2 + &
            (emf2(i,j,k) + emf2(i+1,j,k) + emf2(i,j,k+1) + emf2(i+1,j,k+1))**2 + &
            (emf3(i,j,k) + emf3(i+1,j,k) + emf3(i,j+1,k) + emf3(i+1,j+1,k))**2) &
            /16.d0/(ovrmvar(i+myid1*(in-5)) + tiny)/d0(i+myid1*(in-5)) &
            /temp0(i+myid1*(in-5))
    enddo; enddo; enddo

    if (DoWriteEmf) then
       write(idcpu,'(i4.4)') myid
       open(unit=17,file='emf.cpu'//idcpu//'.resis',form='unformatted', &
            access='stream')
       write(17) in, jn, kn
       write(17) (((emf1(i,j,k),i=ism2,iep2),j=jsm2,jep3),k=ksm2,kep3)
       write(17) (((emf2(i,j,k),i=ism2,iep3),j=jsm2,jep2),k=ksm2,kep3)
       write(17) (((emf3(i,j,k),i=ism2,iep3),j=jsm2,jep3),k=ksm2,kep2)
       close(17)
    endif
    deallocate(emf1, emf2, emf3, curl1, curl2, curl3, eta, bsnp1)

  end subroutine hsmoc
  subroutine diffemf
  end subroutine diffemf
  subroutine resistive_emf
  end subroutine resistive_emf
  
  !===================================================================================

  subroutine heattran
    use Modpar,    ONLY: ijkn, in, jn, kn, myid1
    use ModGrid
    use ModField,  ONLY: b1, b2, b3, s
    use ModBack,   ONLY: ovrth, d0, temp0, s0, bquench_k
    use ModDel,    ONLY: dels
    implicit none

    integer :: i, j, k
    real    :: flux(ijkn), Bmag
    real, allocatable, dimension(:,:,:) :: scratch2
    !---------------------------------------------------------------------------------

    allocate(scratch2(in,jn,kn))

    if (bquench_k > 0.d0) then
       do k=ksm1,kep1; do j=jsm1,jep1; do i=ism1,iep1
          Bmag = sqrt((0.5d0*(b1(i,j,k) + b1(i+1,j,k)))**2 + &
               (0.5d0*(b2(i,j,k) + b2(i,j+1,k)))**2 + (0.5d0*(b3(i,j,k) + b3(i,j,k+1)))**2)
          scratch2(i,j,k) = ovrth(i+myid1*(in-5))/(1.d0 + (Bmag/bquench_k)**2)
       enddo; enddo; enddo
    else
       do k=ksm1,kep1; do j=jsm1,jep1; do i=ism1,iep1
          scratch2(i,j,k) = ovrth(i+myid1*(in-5))
       enddo; enddo; enddo
    endif

    do k=ks,ke; do j=js,je
       do i=is,iep1
          flux(i) = 0.5D0*(scratch2(i-1,j,k) + scratch2(i,j,k)) &
               *0.5D0*(d0(i-1+myid1*(in-5)) + d0(i+myid1*(in-5))) &
               *0.5D0*(temp0(i-1+myid1*(in-5)) + temp0(i+myid1*(in-5))) &
               *(s(i,j,k) + s0(i+myid1*(in-5)) - s(i-1,j,k) - s0(i-1+myid1*(in-5))) &
               *dx1bi(i)*g2a(i)*g31a(i)
       enddo
       do i=is,ie
          dels(i,j,k) = dels(i,j,k) + ((flux(i+1) - flux(i)) &
               *dx1ai(i)*g2bi(i)*g31bi(i))/d0(i+myid1*(in-5))/temp0(i+myid1*(in-5))
       enddo
    enddo; enddo

    do k=ks,ke; do i=is,ie
       do j=js,jep1
          flux(j) = 0.5D0*(scratch2(i,j-1,k) + scratch2(i,j,k)) &
               *(s(i,j,k) - s(i,j-1,k))*dx2bi(j)*g2bi(i)*g32a(j)
       enddo
       do j=js,je
          dels(i,j,k) = dels(i,j,k) + (flux(j+1) - flux(j))*dx2ai(j)*g2bi(i)*g32bi(j)
       enddo
    enddo; enddo

    do j=js,je; do i=is,ie
       do k=ks,kep1
          flux(k) = 0.5D0*(scratch2(i,j,k-1) + scratch2(i,j,k)) &
               *(s(i,j,k) - s(i,j,k-1))*dx3bi(k)*g31bi(i)*g32bi(j)
       enddo
       do k=ks,ke
          dels(i,j,k) = dels(i,j,k) + (flux(k+1) - flux(k))*dx3ai(k)*g31bi(i)*g32bi(j)
       enddo
    enddo; enddo

    deallocate(scratch2)

  end subroutine heattran

end module ModRHS
