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
    real, dimension(1:ijkn) :: bave, bstar
    real, allocatable :: scratch1(:,:,:), scratch2(:,:,:), scratch3(:,:,:)
    !----------------------------------------------------------------------------------

    allocate(scratch1(in,jn,kn), scratch2(in,jn,kn), scratch3(in,jn,kn))
    
    ! 1-force
    ! By following the Alfven velocity in the 2-direction, evaluate
    ! "bstar" from "b1" to estimate the 2-transverse Lorentz force
    do k=ks,ke; do i=is,iep1
       do j=js,jep1
          bave (j) = 0.5d0*(b2(i,  j,k) + b2(i-1,j,k))
          bstar(j) = 0.5d0*(b1(i,j-1,k) + b1(  i,j,k))
       end do
       do j=js,je
          dels1(i,j,k) = dels1(i,j,k) + (g32a(j+1)*bave(j+1)*bstar(j+1) - &
               g32a(j)*bave(j)*bstar(j))*g2ai(i)*g32bi(j)*dx2ai(j)
       end do
    end do; end do
    ! By following the Alfven velocity in the 3-direction, evaluate                         
    ! "bstar" from "b1" to estimate the 3-transverse Lorentz force.                         
    do j=js,je; do i=is,iep1
       do k=ks,kep1
          bave (k) = 0.5d0*(b3(i,j,  k) + b3(i-1,j,k))
          bstar(k) = 0.5d0*(b1(i,j,k-1) + b1(  i,j,k))
       end do
       q1 = g31ai(i)*g32bi(j)
       do k=ks,ke
          dels1(i,j,k) = dels1(i,j,k) + (bave(k+1)*bstar(k+1) - bave(k)*bstar(k))* &
               q1*dx3ai(k)
       end do
    end do; end do
    do k=ks,kep1; do j=js,jep1; do i=is,iep1
       scratch1(i,j,k) = (g2b(i)*b2(i,j,k) + g2b(i-1)*b2(i-1,j,k))* &
            (g2b(i)*b2(i,j,k) - g2b(i-1)*b2(i-1,j,k))*g2ai(i)*g2ai(i)
       scratch2(i,j,k) = (g31b(i)*b3(i,j,k) + g31b(i-1)*b3(i-1,j,k))* &
            (g31b(i)*b3(i,j,k) - g31b(i-1)*b3(i-1,j,k))*g31ai(i)*g31ai(i)
       scratch3(i,j,k) = ((g2b(i)*g31b(i)*(b1(i+1,j,k) + b1(i,j,k)))**2 - &
            (g2b(i-1)*g31b(i-1)*(b1(i,j,k) + b1(i-1,j,k)))**2)*0.5d0* &
            (g2ai(i)*g31ai(i))**2
    enddo; enddo; enddo
    do k=ks,ke; do j=js,je; do i=is,iep1
       dels1(i,j,k)  = dels1(i,j,k) - (scratch1(i,j,k) + scratch1(i,j+1,k) + &
            scratch2(i,j,k) + scratch2(i,j,k+1) - scratch3(i,j,k))*dx1bi(i)*0.25d0
    enddo; enddo; enddo

    ! 2-force
    ! By following the Alfven velocity in the 3-direction, evaluate
    ! "bstar" from "b2" to estimate the 3-transverse Lorentz force
    do j=js,jep1; do i=is,ie
       do k=ks,kep1
          bave (k) = 0.5d0*(b3(i,j,  k) + b3(i,j-1,k))
          bstar(k) = 0.5d0*(b2(i,j,k-1) + b2(i,  j,k))
       end do
       q1 = g31bi(i)*g32ai(j)
       do k=ks,ke
          dels2(i,j,k) = dels2(i,j,k) + (bave(k+1)*bstar(k+1) - bave(k)*bstar(k))* &
               q1*dx3ai(k)
       end do
    end do; end do
    ! By following the Alfven velocity in the 1-direction, evaluate                         
    ! "bstar" from "b2" to estimate the 1-transverse Lorentz force.                         
    do k=ks,ke; do j=js,jep1
       do i=is,iep1
          bave (i) = 0.5d0*(b1(  i,j,k) + b1(i,j-1,k))
          bstar(i) = 0.5d0*(b2(i-1,j,k) + b2(i,  j,k))*g2a(i)
       end do
       do i=is,ie
          dels2(i,j,k) = dels2(i,j,k) + (g2a(i+1)*g31a(i+1)*bave(i+1)*bstar(i+1) - &
               g2a(i)*g31a(i)*bave(i)*bstar(i))*g2bi(i)*g2bi(i)*g31bi(i)*dx1ai(i)
       end do
    end do; end do
    do k=ks,kep1; do j=js,jep1; do i=is,iep1
       scratch1(i,j,k) = (g32b(j)*b3(i,j,k) + g32b(j-1)*b3(i,j-1,k))* &
            (g32b(j)*b3(i,j,k) - g32b(j-1)*b3(i,j-1,k))* g32ai(j)*g32ai(j)
       scratch2(i,j,k) = (b1(i,j,k) + b1(i,j-1,k))*(b1(i,j,k) - b1(i,j-1,k))
       scratch3(i,j,k) = ((g32b(j)*(b2(i,j,k) + b2(i,j+1,k)))**2 - &
            (g32b(j-1)*(b2(i,j,k)+b2(i,j-1,k)))**2)*0.5*g32ai(j)*g32ai(j)
    enddo; enddo; enddo
    do k=ks,ke; do j=js,jep1; do i=is,ie
       dels2(i,j,k)  = dels2(i,j,k) - (scratch1(i,j,k) + scratch1(i,j,k+1) + &
            scratch2(i,j,k) + scratch2(i+1,j,k) - scratch3(i,j,k))*dx2bi(j)*g2bi(i)*0.25d0
    enddo; enddo; enddo
    
    ! 3-force                                                                               
    do k=ks,kep1
       ! By following the Alfven velocity in the 1-direction, evaluate                      
       !"bstar" from "b3" to estimate the 1-transverse Lorentz force.                       
       do j=js,je
          do i=is,iep1
             bave (i) = 0.5d0*(b1(  i,j,k) + b1(i,j,k-1))
             bstar(i) = 0.5d0*(b3(i-1,j,k) + b3(i,j,  k))*g31a(i)
          end do
          do i=is,ie
             dels3(i,j,k) = dels3(i,j,k) + (g2a(i+1)*g31a(i+1)*bave(i+1)*bstar(i+1) - &
                  g2a(i)*g31a(i)*bave(i)*bstar(i))*g2bi(i)*g31bi(i)*g31bi(i)*dx1ai(i)
          end do
       end do
       ! By following the Alfven velocity in the 2-direction, evaluate                      
       ! "bstar" from "b3" to estimate the 2-transverse Lorentz force.                      
       do i=is,ie
          do j=js,jep1
             bave (j) = 0.5d0*(b2(i,  j,k) + b2(i,j,k-1))
             bstar(j) = 0.5d0*(b3(i,j-1,k) + b3(i,j,  k))*g32a(j)
          end do
          do j=js,je
             dels3(i,j,k) = dels3(i,j,k) + (g32a(j+1)*bave(j+1)*bstar(j+1) - &
                  g32a(j)*bave(j)*bstar(j))*g32bi(j)*g32bi(j)*g2bi(i)*dx2ai(j)
          end do
       end do
    end do
    do k=ks,kep1; do j=js,jep1; do i=is,iep1
       scratch1(i,j,k) = (b1(i,j,k) + b1(i,j,k-1))*(b1(i,j,k) - b1(i,j,k-1))
       scratch2(i,j,k) = (b2(i,j,k) + b2(i,j,k-1))*(b2(i,j,k) - b2(i,j,k-1))
       scratch3(i,j,k) = ((b3(i,j,k) + b3(i,j,k+1))**2 - (b3(i,j,k)+b3(i,j,k-1))**2) &
            *0.5d0
    enddo; enddo; enddo
    do k=ks,kep1; do j=js,je; do i=is,ie
       dels3(i,j,k)  = dels3(i,j,k) - (scratch1(i,j,k) + scratch1(i+1,j,k) + &
            scratch2(i,j,k) + scratch2(i,j+1,k) - scratch3(i,j,k))*dx3bi(k)* &
            g31bi(i)*g32bi(j)*0.25d0
    enddo; enddo; enddo

    deallocate(scratch1, scratch2, scratch3)
  
  end subroutine lorentz
  
  !=======================================================================================
  
  subroutine advection
    use ModPar,     ONLY: ijkn, myid1, tiny, in, jn, kn, ijkn, myid
    use ModGrid
    use ModField,   ONLY: b1, b2, b3, v1, v2, v3, s
    use ModBack,    ONLY: d0, temp0, s0, gamma, gacc_str
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
    real, allocatable, dimension(:,:,:) :: scratch1, scratch2, scratch3
    !---------------------------------------------------------------------------------
    
    allocate(scratch1(in,jn,kn), scratch2(in,jn,kn), scratch3(in,jn,kn))
    nlf   = nlf_v
    alpha = alpha_v
    
    ! advection for dels1
    do k=ks,ke; do j=js,je
       qint(ism2:iep3) = v1(ism2:iep3,j,k)
       call xtvd_a(in, ijkn, dx1a, dx1ai, qint, dq, qlr)
       do i=ism1,iep1
          cmax = sqrt((0.5d0*(b1(i,j,k) + b1(i+1,j,k)))**2 + &
               (0.5d0*(b2(i,j,k) + b2(i,j+1,k)))**2 + &
               (0.5d0*(b3(i,j,k) + b3(i,j,k+1)))**2)/sqrt(d0(i+myid1*(in-5)))
          rr = abs(dq(i))/(abs(v1(i+1,j,k) - v1(i,j,k)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          q1 = 0.5d0*(abs(0.5d0*(v1(i,j,k) + v1(i+1,j,k))) + cmax*qq)* &
               d0(i+myid1*(in-5))*dq(i)
          flux(i) = (d0(i+myid1*(in-5))*0.5d0*(v1(i,j,k) + v1(i+1,j,k))*qlr(i) - &
               q1)*g2b(i)*g31b(i)
          scratch1(i,j,k) = q1*(v1(i+1,j,k) - v1(i,j,k))*dx1ai(i)
       enddo
       do i=is,iep1
          dels1(i,j,k) = dels1(i,j,k) - (flux(i) - flux(i-1))*dx1bi(i)*g2ai(i)*g31ai(i)
       enddo
    enddo; enddo
    do k=ks,ke; do i=is,iep1
       qint(jsm2:jep2) = v1(i,jsm2:jep2,k)
       call xtvd(jn, ijkn, dx2a, dx2bi, qint, dq, qlr)
       do j=js,jep1
          cmax = sqrt((0.5d0*(b1(i,j-1,k) + b1(i,j,k)))**2 + &
               (0.5d0*(b2(i-1,j,k) + b2(i,j,k)))**2 + &
               (0.125d0*(b3(i,j,k) + b3(i,j,k+1) + b3(i-1,j,k) + b3(i-1,j,k+1) + &
               b3(i-1,j-1,k) + b3(i-1,j-1,k+1) + b3(i,j-1,k) + b3(i,j-1,k+1)))**2)/ &
               sqrt(0.5d0*(d0(i-1+myid1*(in-5)) + d0(i+myid1*(in-5))))
          rr = abs(dq(j))/(abs(v1(i,j,k) - v1(i,j-1,k)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          q1 = 0.5d0*(abs(0.5d0*(v2(i-1,j,k) + v2(i,j,k))) + cmax*qq)* &
               0.5d0*(d0(i-1+myid1*(in-5)) + d0(i+myid1*(in-5)))*dq(j)
          flux(j) = (0.5d0*(d0(i-1+myid1*(in-5)) + d0(i+myid1*(in-5)))* &
               0.5d0*(v2(i,j,k) + v2(i-1,j,k))*qlr(j) - q1)*g32a(j)
          scratch2(i,j,k) = q1*(v1(i,j,k) - v1(i,j-1,k))*g2ai(i)*dx2bi(j)
       enddo
       do j=js,je
          dels1(i,j,k) = dels1(i,j,k) - (flux(j+1) - flux(j))*dx2ai(j)*g2ai(i)*g32bi(j)
       enddo
    enddo; enddo
    do j=js,je; do i=is,iep1
       qint(ksm2:kep2) = v1(i,j,ksm2:kep2)
       call xtvd(kn, ijkn, dx3a, dx3bi, qint, dq, qlr)
       do k=ks,kep1
          cmax = sqrt((0.5d0*(b1(i,j,k-1) + b1(i,j,k)))**2 + &
               (0.5d0*(b3(i-1,j,k) + b3(i,j,k)))**2 + &
               (0.125d0*(b2(i,j,k) + b2(i,j+1,k) + b2(i-1,j,k) + b2(i-1,j+1,k) + &
               b2(i-1,j,k-1) + b2(i-1,j+1,k-1) + b2(i,j,k-1) + b2(i,j+1,k-1)))**2)/ &
               sqrt(0.5d0*(d0(i-1+myid1*(in-5)) + d0(i+myid1*(in-5))))
          rr = abs(dq(k))/(abs(v1(i,j,k) - v1(i,j,k-1)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          q1 = 0.5d0*(abs(0.5d0*(v3(i-1,j,k) + v3(i,j,k))) + cmax*qq)* &
               0.5d0*(d0(i-1+myid1*(in-5)) + d0(i+myid1*(in-5)))*dq(k)
          flux(k) = 0.5d0*(d0(i-1+myid1*(in-5)) + d0(i+myid1*(in-5)))*0.5d0* &
               (v3(i-1,j,k) + v3(i,j,k))*qlr(k) - q1
          scratch3(i,j,k) = q1*(v1(i,j,k)-v1(i,j,k-1))*g31ai(i)*g32bi(j)*dx3bi(k)
       enddo
       do k=ks,ke
          dels1(i,j,k) = dels1(i,j,k) - (flux(k+1) - flux(k))*dx3ai(k)*g31ai(i)*g32bi(j)
       enddo
    enddo; enddo
    ! add heating due to slope limited diffusion
    do k=ks,ke; do j=js,je; do i=is,ie
       dels(i,j,k) = 1.0d0/d0(i+myid1*(in-5))/temp0(i+myid1*(in-5))* &
            (scratch1(i,j,k) + (scratch2(i,j,k) + scratch2(i+1,j,k) + &
            scratch2(i,j+1,k) + scratch2(i+1,j+1,k))*0.25d0 + &
            (scratch3(i,j,k) + scratch3(i+1,j,k) + scratch3(i,j,k+1) + &
            scratch3(i+1,j,k+1))*0.25d0)
    enddo; enddo; enddo

    ! advection for dels2
     do k=ks,ke; do j=js,jep1
       qint(ism2:iep2) = v2(ism2:iep2,j,k)*g2bi(ism2:iep2)
       call xtvd(in, ijkn, dx1a, dx1bi, qint, dq, qlr)
       do i=is,iep1
          cmax = sqrt((0.5d0*(b1(i,j-1,k) + b1(i,j,k)))**2 + &
               (0.5d0*(b2(i-1,j,k) + b2(i,j,k)))**2 + &
               (0.125d0*(b3(i,j,k) + b3(i,j,k+1) + b3(i,j-1,k) + b3(i,j-1,k+1) + &
               b3(i-1,j-1,k) + b3(i-1,j-1,k+1) + b3(i-1,j,k) + b3(i-1,j,k+1)))**2)/ &
               sqrt(0.5d0*(d0(i-1+myid1*(in-5)) + d0(i+myid1*(in-5))))
          rr = abs(dq(i))/(abs(qint(i) - qint(i-1)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          q1 = 0.5d0*(abs(0.5*(v1(i,j-1,k) + v1(i,j,k))) + cmax*qq) &
               *0.5d0*(d0(i-1+myid1*(in-5)) + d0(i+myid1*(in-5)))*dq(i)
          flux(i) = (0.5d0*(d0(i-1+myid1*(in-5)) + d0(i+myid1*(in-5))) &
               *0.5d0*(v1(i,j,k) + v1(i,j-1,k))*qlr(i) - q1) &
               *g2a(i)*g31a(i)*(g2a(i))**2
          scratch1(i,j,k) = q1*(g2a(i))**2*(qint(i) - qint(i-1))*dx1bi(i)
       enddo
       do i=is,ie
          dels2(i,j,k) = dels2(i,j,k) - (flux(i+1) - flux(i))* &
               dx1ai(i)*g2bi(i)*g31bi(i)*g2bi(i)
       enddo
    enddo; enddo
    do k=ks,ke; do i=is,ie
       qint(jsm2:jep3) = v2(i,jsm2:jep3,k)*g2bi(i)
       call xtvd_a(jn, ijkn, dx2a, dx2ai, qint, dq, qlr)
       do j=jsm1,jep1
          cmax = sqrt((0.5d0*(b2(i,j,k) + b2(i,j+1,k)))**2 + &
               (0.5d0*(b3(i,j,k) + b3(i,j,k+1)))**2 + &
               (0.5d0*(b1(i,j,k) + b1(i+1,j,k)))**2)/sqrt(d0(i+myid1*(in-5)))
          rr = abs(dq(j))/(abs(qint(j+1) - qint(j)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          q1 = 0.5d0*(abs(0.5d0*(v2(i,j,k) + v2(i,j+1,k))) + cmax*qq)* &
               d0(i+myid1*(in-5))*dq(j)
          flux(j) = (d0(i+myid1*(in-5))*0.5d0*(v2(i,j,k) + v2(i,j+1,k))*qlr(j) - q1) &
               *g32b(j)*(g2b(i))**2
          scratch2(i,j,k) = q1*(g2b(i))**2*(qint(j+1) - qint(j))*g2bi(i)*dx2ai(j)
       enddo
       do j=js,jep1
          dels2(i,j,k) = dels2(i,j,k) - (flux(j) - flux(j-1))* &
               dx2bi(j)*g2bi(i)*g32ai(j)*g2bi(i)
       enddo
    enddo; enddo
    do j=js,jep1; do i=is,ie
       qint(ksm2:kep2) = v2(i,j,ksm2:kep2)*g2bi(i)
       call xtvd(kn, ijkn, dx3a, dx3bi, qint, dq, qlr)
       do k=ks,kep1
          cmax = sqrt((0.5d0*(b2(i,j,k-1) + b2(i,j,k)))**2 + &
               (0.5d0*(b3(i,j-1,k) + b3(i,j,k)))**2 + &
               (0.125d0*(b1(i,j,k) + b1(i+1,j,k) + b1(i,j-1,k) + b1(i+1,j-1,k) + &
               b1(i,j-1,k-1) + b1(i+1,j-1,k-1) + b1(i,j,k-1) + b1(i+1,j,k-1)))**2)/ &
               sqrt(d0(i+myid1*(in-5)))
          rr = abs(dq(k))/(abs(qint(k) - qint(k-1)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          q1 = 0.5d0*(abs(0.5d0*(v3(i,j-1,k) + v3(i,j,k))) + cmax*qq)* &
               d0(i+myid1*(in-5))*dq(k)
          flux(k) = (d0(i+myid1*(in-5))*0.5d0*(v3(i,j-1,k) + v3(i,j,k))*qlr(k) - q1) &
               *(g2b(i))**2
          scratch3(i,j,k) = q1*(g2b(i))**2*(qint(k) - qint(k-1)) &
               *g31bi(i)*g32ai(j)*dx3bi(k)
       enddo
       do k=ks,ke
          dels2(i,j,k) = dels2(i,j,k) - (flux(k+1) - flux(k))* &
               dx3ai(k)*g31bi(i)*g32ai(j)*g2bi(i)
       enddo
    enddo; enddo
    ! add heating due to slope limited diffusion                                            
    do k=ks,ke; do j=js,je; do i=is,ie
       dels(i,j,k) = dels(i,j,k) + 1.d0/d0(i+myid1*(in-5))/temp0(i+myid1*(in-5))* &
            ((scratch1(i,j,k) + scratch1(i,j+1,k) + scratch1(i+1,j,k) + &
            scratch1(i+1,j+1,k))*0.25d0 + scratch2(i,j,k) + &
            (scratch3(i,j,k) + scratch3(i,j+1,k) + scratch3(i,j,k+1) + &
            scratch3(i,j+1,k+1))*0.25d0)
    enddo; enddo; enddo

    ! advection for dels3   
    do k=ks,kep1; do j=js,je
       qint(ism2:iep2) = v3(ism2:iep2,j,k)*g31bi(ism2:iep2)*g32bi(j)
       call xtvd(in, ijkn, dx1a, dx1bi, qint, dq, qlr)
       do i = is,iep1
          cmax = sqrt((0.5d0*(b1(i,j,k-1) + b1(i,j,k)))**2 + &
               (0.5d0*(b3(i-1,j,k) + b3(i,j,k)))**2 + &
               (0.125d0*(b2(i,j,k) + b2(i,j+1,k) + b2(i,j,k-1) + b2(i,j+1,k-1) + &
               b2(i-1,j,k-1) + b2(i-1,j+1,k-1) + b2(i-1,j,k) + b2(i-1,j+1,k)))**2)/ &
               sqrt(0.5*(d0(i-1+myid1*(in-5)) + d0(i+myid1*(in-5))))
          rr = abs(dq(i))/(abs(qint(i) - qint(i-1)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          q1 = 0.5d0*(abs(0.5d0*(v1(i,j,k-1) + v1(i,j,k))) + cmax*qq)* &
               0.5d0*(d0(i-1+myid1*(in-5)) + d0(i+myid1*(in-5)))*dq(i)
          flux(i) = (0.5d0*(d0(i-1+myid1*(in-5)) + d0(i+myid1*(in-5)))* &
               0.5d0*(v1(i,j,k) + v1(i,j,k-1))*(qlr(i) + ovrb) - q1)&
               *g2a(i)*g31a(i)*(g31a(i)*g32b(j))**2
          scratch1(i,j,k) = q1*(g31a(i)*g32b(j))**2*(qint(i) - qint(i-1))*dx1bi(i)
       enddo
       do i=is,ie
          dels3(i,j,k) = dels3(i,j,k) - (flux(i+1) - flux(i))* &
               dx1ai(i)*g2bi(i)*g31bi(i)*g31bi(i)*g32bi(j)
       enddo
    enddo; enddo
    do k=ks,kep1; do i=is,ie
       qint(jsm2:jep2) = v3(i,jsm2:jep2,k)*g31bi(i)*g32bi(jsm2:jep2)
       call xtvd(jn, ijkn, dx2a, dx2bi, qint, dq, qlr)
       do j=js,jep1
          cmax = sqrt((0.5d0*(b3(i,j-1,k) + b3(i,j,k)))**2 + &
               (0.5d0*(b2(i,j,k-1) + b2(i,j,k)))**2 + &
               (0.125d0*(b1(i,j,k) + b1(i+1,j,k) + b1(i,j,k-1) + b1(i+1,j,k-1) + &
               b1(i,j-1,k-1) + b1(i+1,j-1,k-1) + b1(i,j-1,k) + b1(i+1,j-1,k)))**2)/ &
               sqrt(d0(i+myid1*(in-5)))
          rr = abs(dq(j))/(abs(qint(j) - qint(j-1)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          q1 = 0.5d0*(abs(0.5d0*(v2(i,j,k-1) + v2(i,j,k))) + cmax*qq)* &
               d0(i+myid1*(in-5))*dq(j)
          flux(j) = (d0(i+myid1*(in-5))*0.5d0*(v2(i,j,k-1) + v2(i,j,k))* &
               (qlr(j) + ovrb)- q1)*g32a(j)*(g31b(i)*g32a(j))**2
          scratch2(i,j,k) = q1*(g31b(i)*g32a(j))**2*(qint(j) - qint(j-1))*g2bi(i)*dx2bi(j)
       enddo
       do j=js,je
          dels3(i,j,k) = dels3(i,j,k) - (flux(j+1) - flux(j))* &
               dx2ai(j)*g2bi(i)*g32bi(j)*g31bi(i)*g32bi(j)
       enddo
    enddo; enddo
    do j=js,je; do i=is,ie
       qint(ksm2:kep3) = v3(i,j,ksm2:kep3)*g31bi(i)*g32bi(j)
       call xtvd_a(kn, ijkn, dx3a, dx3ai, qint, dq, qlr)
       do k=ksm1,kep1
          cmax = sqrt((0.5d0*(b3(i,j,k) + b3(i,j,k+1)))**2 + &
               (0.5d0*(b1(i,j,k) + b1(i+1,j,k)))**2 + &
               (0.5d0*(b2(i,j,k) + b2(i,j+1,k)))**2)/ &
               sqrt(d0(i+myid1*(in-5)))
          rr = abs(dq(k))/(abs(qint(k+1) - qint(k)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          q1 = 0.5d0*(abs(0.5d0*(v3(i,j,k) + v3(i,j,k+1))) + cmax*qq)* &
               d0(i+myid1*(in-5))*dq(k)
          flux(k) = (d0(i+myid1*(in-5))*0.5d0*(v3(i,j,k) + v3(i,j,k+1))*&
               (qlr(k) + ovrb) - q1)*(g31b(i)*g32b(j))**2
          scratch3(i,j,k) = q1*(g31b(i)*g32b(j))**2*(qint(k+1) - qint(k)) &
               *g31bi(i)*g32bi(j)*dx3ai(k)
       enddo
       do k=ks,kep1
          dels3(i,j,k) = dels3(i,j,k) - (flux(k) - flux(k-1))* &
               dx3bi(k)*g31bi(i)*g32bi(j)*g31bi(i)*g32bi(j)
       enddo
    enddo; enddo
    ! add heating due to slope limited diffusion 
    do k=ks,ke; do j=js,je; do i=is,ie
       dels(i,j,k) = dels(i,j,k) + 1.d0/d0(i+myid1*(in-5))/temp0(i+myid1*(in-5))* &
            ((scratch1(i,j,k) + scratch1(i,j,k+1) + scratch1(i+1,j,k) + &
            scratch1(i+1,j,k+1))*0.25d0 + (scratch2(i,j,k) + scratch2(i,j,k+1) + &
            scratch2(i,j+1,k) + scratch2(i,j+1,k+1))*0.25d0 + scratch3(i,j,k))
    enddo; enddo; enddo

    ! advection for dels
    nlf=nlf_s
    alpha=alpha_s
    do k=ks,ke; do j=js,je
       qint(ism2:iep2) = s(ism2:iep2,j,k) + s0(ism2+myid1*(in-5):iep2+myid1*(in-5))
       call xtvd(in, ijkn, dx1a, dx1bi, qint, dq, qlr)
       do i = is,iep1
          cmax = sqrt(b1(i,j,k)**2 + &
               (0.25d0*(b2(i,j,k) + b2(i,j+1,k) + b2(i-1,j,k) + b2(i-1,j+1,k)))**2 + &
               (0.25d0*(b3(i,j,k) + b3(i,j,k+1) + b3(i-1,j,k) + b3(i-1,j,k+1)))**2)/ &
               sqrt(0.5d0*(d0(i-1+myid1*(in-5)) + d0(i+myid1*(in-5))))
          rr = abs(dq(i))/(abs(qint(i) - qint(i-1)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          flux(i) = (0.5d0*(d0(i-1+myid1*(in-5)) + d0(i+myid1*(in-5)))* &
               0.5d0*(temp0(i-1+myid1*(in-5)) + temp0(i+myid1*(in-5)))* &
               v1(i,j,k)*qlr(i) - 0.5d0*(abs(v1(i,j,k)) + cmax*qq)*&
               0.5d0*(d0(i-1+myid1*(in-5)) + d0(i+myid1*(in-5)))* &
               0.5d0*(temp0(i-1+myid1*(in-5)) + temp0(i+myid1*(in-5)))* &
               dq(i))*g2a(i)*g31a(i)
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
       qint(jsm2:jep2) = s(i,jsm2:jep2,k) + s0(i+myid1*(in-5))
       call xtvd(jn, ijkn, dx2a, dx2bi, qint, dq, qlr)
       do j=js,jep1
          cmax = sqrt(b2(i,j,k)**2 + &
               (0.25d0*(b3(i,j,k) + b3(i,j,k+1) + b3(i,j-1,k) + b3(i,j-1,k+1)))**2 + &
               (0.25d0*(b1(i,j,k) + b1(i+1,j,k) + b1(i,j-1,k) + b1(i+1,j-1,k)))**2)/ &
               sqrt(d0(i+myid1*(in-5)))
          rr = abs(dq(j))/(abs(s(i,j,k) - s(i,j-1,k)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          flux(j) = (v2(i,j,k)*qlr(j) - 0.5d0*(abs(v2(i,j,k)) + cmax*qq)*dq(j))*g32a(j)
       enddo
       do j=js,je
          dels(i,j,k) = dels(i,j,k) - (flux(j+1) - flux(j))*dx2ai(j)*g2bi(i)*g32bi(j)
       enddo
    enddo; enddo
    do j=js,je; do i=is,ie
       qint(ksm2:kep2) = s(i,j,ksm2:kep2) + s0(i+myid1*(in-5))
       call xtvd(kn, ijkn, dx3a, dx3bi, qint, dq, qlr)
       do k=ks,kep1
          cmax = sqrt(b3(i,j,k)**2 + &
               (0.25d0*(b1(i,j,k) + b1(i+1,j,k) + b1(i,j,k-1) + b1(i+1,j,k-1)))**2 + &
               (0.25d0*(b2(i,j,k) + b2(i,j+1,k) + b2(i,j,k-1) + b2(i,j+1,k-1)))**2)/ &
               sqrt(d0(i+myid1*(in-5)))
          rr = abs(dq(k))/(abs(s(i,j,k) - s(i,j,k-1)) + tiny)
          qq = (min(1.d0, alpha*rr))**nlf
          flux(k) = v3(i,j,k)*qlr(k) - 0.5d0*(abs(v3(i,j,k)) + cmax*qq)*dq(k)
       enddo
       do k=ks,ke
          dels(i,j,k) = dels(i,j,k) - (flux(k+1) - flux(k))*dx3ai(k)*g31bi(i)*g32bi(j)
       enddo
    enddo; enddo
    
    deallocate(scratch1, scratch2, scratch3)
    
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
    use ModPar,      ONLY: in, jn, kn, ijkn, myid1, pi, myid
    use ModSundry,   ONLY: dt
    use ModGrid
    use ModBack,     ONLY: d0
    use ModField,    ONLY: b1, b2, b3, v1, v2, v3
    use ModDel,      ONLY: emf1s, emf2s, emf3s
    use ModBval,     ONLY: bvalemf1, bvalemf2, bvalemf3
    use ModInterp,   ONLY: xint1d, xzc1d
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
    
    real, allocatable :: vsnp1(:,:,:), bsnp1(:,:,:)
    ! paramters for testing, write emf1-3s if DoWriteEmf
    character(len=4) :: idcpu
    logical :: DoWriteEmf = .false.
    !-----------------------------------------------------------------------------------
    
    allocate(vsnp1(in,jn,kn), bsnp1(in,jn,kn))
    
    do i = ism1,iep1
       srd1(i) = sqrt(0.5d0*(d0(i+myid1*(in-5)) + d0(i-1+myid1*(in-5))))
       srd2(i) = sqrt(d0(i+myid1*(in-5)))
       srd3(i) = srd2(i)
    enddo
    
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
             vsnp1(i,j,k) = (vpch(k)*srd + vmch(k)*srd + q2*(bpch(k) - bmch(k))) &
                  /(2.d0*srd)*bave(k)
             bsnp1(i,j,k) = (bpch(k)*srdi + bmch (k)*srdi + q2*(vpch(k) - vmch(k))) &
                  /(2.d0*srdi)*vave(k)
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
          call xint1d(kn, bt, vt, vfl, gfact, dx3bi, dx3a, bint, vint)
          v2intk(j,ks:kep1) = vint(ks:kep1)
          b2intk(j,ks:kep1) = bint(ks:kep1)
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
             vsnp1(i,j,k) = 0.5d0*(vsnp1(i,j,k) + bsnm1)
             bsnp1(i,j,k) = 0.5d0*(vsnm1 + bsnp1(i,j,k))
          enddo
       enddo
    enddo
    emf1s = vsnp1 - bsnp1
    call bvalemf1(emf1s)
    
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
          call xint1d(kn, bt, vt, vfl, gfact, dx3bi, dx3a, bint, vint)
          v1intk(ks:kep1,i) = vint(ks:kep1)
          b1intk(ks:kep1,i) = bint(ks:kep1)
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
             vsnp1(i,j,k) = (vpch(i)*srdp(i) + vmch(i)*srdm (i) + &
                  q2*(bpch(i) - bmch(i)))/(srdp(i) + srdm(i))*bave(i)
             bsnp1(i,j,k) = (bpch(i)*srdpi(i) + bmch (i)*srdmi(i) + &
                  q2*(vpch(i) - vmch(i)))/(srdpi(i) + srdmi(i))*vave(i)
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
             vsnp1(i,j,k) = 0.5d0*(vsnp1(i,j,k) + bsnm1)
             bsnp1(i,j,k) = 0.5d0*(vsnm1 + bsnp1(i,j,k))
          enddo
       enddo
    enddo
    emf2s = vsnp1 - bsnp1
    call bvalemf2(emf2s)
    
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
             vsnp1(i,j,k) = (vpch(j)*srd + vmch(j)*srd + q2*(bpch(j) - bmch(j))) &
                  /(2.d0*srd)*bave(j)
             bsnp1(i,j,k) = (bpch(j)*srdi + bmch(j)*srdi + q2*(vpch(j) - vmch(j))) &
                  /(2.d0*srdi)*vave(j)
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
          call xint1d(jn, bt, vt, vfl, gfact, dx2bi, dx2a, bint, vint)
          b1intj(i,js:jep1) = bint(js:jep1)
          v1intj(i,js:jep1) = vint(js:jep1)
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
             vsnp1(i,j,k) = 0.5d0*(vsnp1(i,j,k) + bsnm1)
             bsnp1(i,j,k) = 0.5d0*(vsnm1 + bsnp1(i,j,k))
          enddo
       enddo
    enddo
    emf3s = vsnp1 - bsnp1
    call bvalemf3(emf3s)
    
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
    
    deallocate(vsnp1, bsnp1)
    
  end subroutine hsmoc
  
  !======================================================================================= 
  
  subroutine diffemf
    use ModPar,     ONLY: in, jn, kn, ijkn, pi, myid1, myid
    use ModGrid
    use ModField,   ONLY: b1, b2, b3
    use ModBack,    ONLY: d0
    use ModSundry,  ONLY: dtmdi2
    use ModDel,     ONLY: emf1s, emf2s, emf3s
    use ModBval,    ONLY: bvaleta, bvalemf1, bvalemf2, bvalemf3
    use ModInterp,  ONLY: xdel
    use ModMpi
    use ModFSAM,    ONLY: iComm
    implicit none
    
    integer :: i, j, k, ierr
    real :: dr1, dr2, dr3, drmin, current, dtmdi2tmp, dtmdi2lc
    real :: bint(ijkn), dbint(ijkn)
    real, allocatable, dimension(:,:,:) :: emf1, emf2, emf3, eta, vsnp1, bsnp1
    
    ! paramters for testing, write emf1-3s if DoWriteEmf
    character(len=4) :: idcpu
    logical :: DoWriteEmf = .false.
    !-----------------------------------------------------------------------------------
    
    allocate(emf1(in,jn,kn), emf2(in,jn,kn), emf3(in,jn,kn))
    allocate(eta(in,jn,kn), vsnp1(in,jn,kn), bsnp1(in,jn,kn))

    ! emf1
    do i=is,ie
       ! del2(b3 h3 dx3)
       do k=ks,kep1
          bint(jsm2:jep2) = b3(i,jsm2:jep2,k)*g31b(i)*g32b(jsm2:jep2)*dx3b(k)
          if(x2a(js) .eq. 0.0) then
             bint(jsm1) = - bint(js)
             bint(jsm2) = - bint(jsp1)
          endif
          if(x2a(jep1) .eq. pi) then
             bint(jep1) = - bint(je)
             bint(jep2) = - bint(jem1)
          endif
          call xdel(jn, bint, dbint, dx2a, dx2bi)
          do j=js,jep1
             vsnp1(i,j,k) = - dbint(j)*g2bi(i)*dx2bi(j)*g31bi(i)*g32ai(j)*dx3bi(k)
          enddo
       enddo
       ! del3(b2 h2 dx2) 
       do j=js,jep1
          bint(ksm2:kep2) = b2(i,j,ksm2:kep2)*g2b(i)*dx2b(j)
          call xdel(kn, bint, dbint, dx3a, dx3bi)
          do k=ks,kep1
             bsnp1(i,j,k) = -dbint(k)*g2bi(i)*dx2bi(j)*g31bi(i)*g32ai(j)*dx3bi(k)
          enddo
       enddo
    enddo
    emf1(is:ie,js:jep1,ks:kep1) = vsnp1(is:ie,js:jep1,ks:kep1) - &
         bsnp1(is:ie,js:jep1,ks:kep1)
    
    !emf2
    do j=js,je
       ! del3(b1 h1 dx1)
       do i=is,iep1
          bint(ksm2:kep2) = b1(i,j,ksm2:kep2)*dx1b(i)
          call xdel(kn, bint, dbint, dx3a, dx3bi)
          do k=ks,kep1
             vsnp1(i,j,k) = - dbint(k)*g31ai(i)*g32bi(j)*dx3bi(k)*dx1bi(i)
          enddo
       enddo
       ! del1(b3 h3 dx3)
       do k=ks,kep1
          bint(ism2:iep2) = b3(ism2:iep2,j,k)*g31b(ism2:iep2)*g32b(j)*dx3b(k)
          call xdel(in, bint, dbint, dx1a, dx1bi)
          do i=is,iep1
             bsnp1(i,j,k) = - dbint(i)*g31ai(i)*g32bi(j)*dx3bi(k)*dx1bi(i)
          enddo
       enddo
    enddo
    emf2(is:iep1,js:je,ks:kep1) = vsnp1(is:iep1,js:je,ks:kep1) - &
         bsnp1(is:iep1,js:je,ks:kep1)
    
    ! emf3
    do k=ks,ke
       ! del1(b2 h2 dx2)
       do j=js,jep1
          bint(ism2:iep2) = b2(ism2:iep2,j,k)*g2b(ism2:iep2)*dx2b(j)
          call xdel(in, bint, dbint, dx1a, dx1bi)
          do i=is,iep1
             vsnp1(i,j,k) = - dbint(i)*dx1bi(i)*g2ai(i)*dx2bi(j)
          enddo
       enddo
       ! del2(b1 h1 dx1)
       do i=is,iep1
          bint(jsm2:jep2) = b1(i,jsm2:jep2,k)*dx1b(i)
          call xdel(jn, bint, dbint, dx2a, dx2bi)
          do j=js,jep1
             bsnp1(i,j,k) = -dbint(j)*dx1bi(i)*g2ai(i)*dx2bi(j)
          enddo
       enddo
    enddo
    emf3(is:iep1,js:jep1,ks:ke) = vsnp1(is:iep1,js:jep1,ks:ke) - &
         bsnp1(is:iep1,js:jep1,ks:ke)
    
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
    do k=ks,kep1; do j=js,jep1; do i=is,ie
       bsnp1(i,j,k) = 0.25d0*(eta(i,j,k) + eta(i,j-1,k) + eta(i,j-1,k-1) + eta(i,j,k-1))
    enddo; enddo; enddo
    do i=is,ie
       do k=ks,kep1
          bint(jsm2:jep2) = b3(i,jsm2:jep2,k)*g31b(i)*g32b(jsm2:jep2)*dx3b(k)
          if(x2a(js) .eq. 0.d0) then
             bint(jsm1) = - bint(js)
             bint(jsm2) = - bint(jsp1)
          endif
          if(x2a(jep1) .eq. pi) then
             bint(jep1) = - bint(je)
             bint(jep2) = - bint(jem1)
          endif
          call xdel(jn, bint, dbint, dx2a, dx2bi)
          do j=js,jep1
             vsnp1(i,j,k) = dbint(j)*( - bsnp1(i,j,k))*g2bi(i)*dx2bi(j) &
                  *g31bi(i)*g32ai(j)*dx3bi(k)
          enddo
       enddo
       do j=js,jep1
          bint(ksm2:kep2) = b2(i,j,ksm2:kep2)*g2b(i)*dx2b(j)
          call xdel(kn, bint, dbint, dx3a, dx3bi)
          do k=ks,kep1
             bsnp1(i,j,k) = dbint(k)*( - bsnp1(i,j,k))*g2bi(i)*dx2bi(j) &
                  *g31bi(i)*g32ai(j)*dx3bi(k)
          enddo
       enddo
    enddo
    emf1 = vsnp1 - bsnp1
    call bvalemf1(emf1)
    emf1s(1:in-1,1:jn,1:kn) = emf1s(1:in-1,1:jn,1:kn) + emf1(1:in-1,1:jn,1:kn)
    
    ! emf2
    do k=ks,kep1; do j=js,je; do i=is,iep1
       bsnp1(i,j,k) = 0.25d0*(eta(i,j,k) + eta(i-1,j,k) + eta(i-1,j,k-1) + eta(i,j,k-1))
    enddo; enddo; enddo
    do j=js,je
       do i=is,iep1
          bint(ksm2:kep2) = b1(i,j,ksm2:kep2)*dx1b(i)
          call xdel(kn, bint, dbint, dx3a, dx3bi)
          do k=ks,kep1
             vsnp1(i,j,k) = dbint(k)*( - bsnp1(i,j,k))*g31ai(i)*g32bi(j)*dx3bi(k) &
                  *dx1bi(i)
          enddo
       enddo
       do k=ks,kep1
          bint(ism2:iep2) = b3(ism2:iep2,j,k)*g31b(ism2:iep2)*g32b(j)*dx3b(k)
          call xdel(in, bint, dbint, dx1a, dx1bi)
          do i=is,iep1
             bsnp1(i,j,k) = dbint(i)*( - bsnp1(i,j,k))*g31ai(i)*g32bi(j)*dx3bi(k) &
                  *dx1bi(i)
          enddo
       enddo
    enddo
    emf2 = vsnp1 - bsnp1
    call bvalemf2(emf2)
    emf2s(1:in,1:jn-1,1:kn) = emf2s(1:in,1:jn-1,1:kn) + emf2(1:in,1:jn-1,1:kn)
    
    ! emf3
    do k=ks,ke; do j=js,jep1; do i=is,iep1
       bsnp1(i,j,k) = 0.25d0*(eta(i,j,k) + eta(i-1,j,k) + eta(i-1,j-1,k) + eta(i,j-1,k))
    enddo; enddo; enddo
    do k=ks,ke
       do j=js,jep1
          bint(ism2:iep2) = b2(ism2:iep2,j,k)*g2b(ism2:iep2)*dx2b(j)
          call xdel(in, bint, dbint, dx1a, dx1bi)
          do i=is,iep1
             vsnp1(i,j,k) = dbint(i)*( - bsnp1(i,j,k))*dx1bi(i)*g2ai(i)*dx2bi(j)
          enddo
       enddo
       do i=is,iep1
          bint(jsm2:jep2) = b1(i,jsm2:jep2,k)*dx1b(i)
          call xdel(jn, bint, dbint, dx2a, dx2bi)
          do j=js,jep1
             bsnp1(i,j,k) = dbint(j)*( - bsnp1(i,j,k))*dx1bi(i)*g2ai(i)*dx2bi(j)
          enddo
       enddo
    enddo
    emf3 = vsnp1 - bsnp1
    call bvalemf3(emf3)
    emf3s(1:in,1:jn,1:kn-1) = emf3s(1:in,1:jn,1:kn-1) + emf3(1:in,1:jn,1:kn-1)
    
    if(DoWriteEmf)then
       write(idcpu,'(i4.4)') myid
       open(unit=17,file='emf.cpu'//idcpu//'.diffemf',form='unformatted', &
            access='stream')
       write(17) in, jn, kn
       write(17) (((emf1(i,j,k),i=ism2,iep2),j=jsm2,jep3),k=ksm2,kep3)
       write(17) (((emf2(i,j,k),i=ism2,iep3),j=jsm2,jep2),k=ksm2,kep3)
       write(17) (((emf3(i,j,k),i=ism2,iep3),j=jsm2,jep3),k=ksm2,kep2)
       close(17)
    endif
    
    deallocate(emf1, emf2, emf3, eta, vsnp1, bsnp1)
    
  end subroutine diffemf
  
  !======================================================================================= 

  subroutine resistive_emf
    use ModPar,    ONLY: myid1, pi, tiny, ijkn, in, jn, kn, myid
    use ModGrid
    use ModField,  ONLY: b1, b2, b3
    use ModBack,   ONLY: ovrmvar, d0, temp0
    use ModDel
    use ModBval,   ONLY: bvalemf1, bvalemf2, bvalemf3
    implicit none
    
    integer :: i, j, k
    real    :: bint(ijkn), dbint
    real, allocatable, dimension(:,:,:) :: emf1, emf2, emf3, vsnp1, bsnp1
    
    ! paramters for testing, write emf1-3s if DoWriteEmf           
    ! paramters for testing, write dels1-3s and dels if DoWriteDels           
    character(len=4) :: idcpu
    logical :: DoWriteEmf = .false.
    logical :: DoWriteDels = .false.
    !----------------------------------------------------------------------------------
    
    allocate(vsnp1(in,jn,kn), bsnp1(in,jn,kn))
    allocate(emf1(in,jn,kn),emf2(in,jn,kn),emf3(in,jn,kn))
    
    ! emf1
    do i=is,ie
       ! calculate del2(b3 h3 dx3)
       do k=ks,kep1
          bint(jsm2:jep2) = b3(i,jsm2:jep2,k)*g31b(i)*g32b(jsm2:jep2)*dx3b(k)
          if(x2a(js) .eq. 0.d0) then
             bint(jsm1) = - bint(js)
             bint(jsm2) = - bint(jsp1)
          endif
          if(x2a(jep1) .eq. pi) then
             bint(jep1) = - bint(je)
             bint(jep2) = - bint(jem1)
          endif
          do j=js,jep1
             dbint = bint(j) - bint(j-1)
             vsnp1(i,j,k) = dbint*(-ovrmvar(i+myid1*(in-5)))*g2bi(i)*dx2bi(j) &
                  *g31bi(i)*g32ai(j)*dx3bi(k)
          enddo
       enddo
       ! calculate del3(b2 h2 dx2)
       do j=js,jep1
          bint(ksm2:kep2) = b2(i,j,ksm2:kep2)*g2b(i)*dx2b(j)
          do k=ks,kep1
             dbint = bint(k) - bint(k-1)
             bsnp1(i,j,k) = dbint*(-ovrmvar(i+myid1*(in-5)))*g2bi(i)*dx2bi(j) &
                  *g31bi(i)*g32ai(j)*dx3bi(k)
          enddo
       enddo
    enddo
    emf1 = vsnp1 - bsnp1
    call bvalemf1(emf1)
    emf1s(1:in-1,1:jn,1:kn) = emf1s(1:in-1,1:jn,1:kn) + emf1(1:in-1,1:jn,1:kn)
    
    ! emf2
    do j=js,je
       ! calculate del3(b1 h1 dx1) 
       do i=is,iep1
          bint(ksm2:kep2) = b1(i,j,ksm2:kep2)*dx1b(i)
          do k=ks,kep1
             dbint = bint(k) - bint(k-1)
             vsnp1(i,j,k) = dbint*(- 0.5d0*(ovrmvar(i+myid1*(in-5)) & 
                  + ovrmvar(i-1+myid1*(in-5))))*g31ai(i)*g32bi(j)*dx3bi(k)*dx1bi(i)
          enddo
       enddo
       ! calculate del1(b3 h3 dx3) 
       do k=ks,kep1
          bint(ism2:iep2) = b3(ism2:iep2,j,k)*g31b(ism2:iep2)*g32b(j)*dx3b(k)
          do i=is,iep1
             dbint = bint(i) - bint(i-1)
             bsnp1(i,j,k) = dbint*(-0.5d0*(ovrmvar(i+myid1*(in-5)) &
                  + ovrmvar(i-1+myid1*(in-5))))*g31ai(i)*g32bi(j)*dx3bi(k)*dx1bi(i)
          enddo
       enddo
    enddo
    emf2 = vsnp1 - bsnp1
    call bvalemf2(emf2)
    emf2s(1:in,1:jn-1,1:kn) = emf2s(1:in,1:jn-1,1:kn) + emf2(1:in,1:jn-1,1:kn)
    
    ! emf3
    do k=ks,ke
       ! calculate del1(b2 h2 dx2)
       do j=js,jep1
          bint(ism2:iep2) = b2(ism2:iep2,j,k)*g2b(ism2:iep2)*dx2b(j)
          do i=is,iep1
             dbint = bint(i) - bint(i-1)
             vsnp1(i,j,k) = dbint*(-0.5d0*(ovrmvar(i+myid1*(in-5)) &
                  + ovrmvar(i-1+myid1*(in-5))))*dx1bi(i)*g2ai(i)*dx2bi(j)
          enddo
       enddo
       ! calculate del2(b1 h1 dx1)
       do i=is,iep1
          bint(jsm2:jep2) = b1(i,jsm2:jep2,k)*dx1b(i)
          do j=js,jep1
             dbint = bint(j)-bint(j-1)
             bsnp1(i,j,k) = dbint*(-0.5d0*(ovrmvar(i+myid1*(in-5))&
                  + ovrmvar(i-1+myid1*(in-5))))*dx1bi(i)*g2ai(i)*dx2bi(j)
          enddo
       enddo
    enddo
    emf3 = vsnp1 - bsnp1
    call bvalemf3(emf3)
    emf3s(1:in,1:jn,1:kn-1) = emf3s(1:in,1:jn,1:kn-1) + emf3(1:in,1:jn,1:kn-1)
    
    ! compute resistive heating 
    do k=ks,ke; do j=js,je; do i=is,ie
       dels(i,j,k) = dels(i,j,k) + &
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
    
    if (DoWriteDels) then
       write(idcpu,'(i4.4)') myid
       open(unit=17,file='dels.cpu'//idcpu,form='unformatted', &
            access='stream')
       write(17) in, jn, kn
       write(17) (((dels1(i,j,k),i=is,iep1),j=js,je),k=ks,ke)
       write(17) (((dels2(i,j,k),i=is,ie),j=js,jep1),k=ks,ke)
       write(17) (((dels3(i,j,k),i=is,ie),j=js,je),k=ks,kep1)
       write(17) (((dels(i,j,k),i=is,ie),j=js,je),k=ks,ke)
       close(17)
    endif
    
    deallocate(vsnp1, bsnp1, emf1, emf2, emf3)
    
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
