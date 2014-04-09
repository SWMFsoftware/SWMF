module ModGetP

  implicit none

  private

  public :: getp

contains

  subroutine getp
    use ModControl
    use ModPar
    use ModGrid
    use ModField,     ONLY: p
    use ModDel
    use ModBack,      ONLY: fact, ovdxxb
    use ModBlk
    use ModBoundary
    use ModBval,      ONLY: bvalp
    use ModFSAM,      ONLY: ErrPSolv
    use ModMpi

    ! local blktri variables
    integer iflg_blk
    integer kk,ll,nww_ch,ierr_blk
    integer, parameter :: np_blk=1,mp_blk=1, &
         n_blk=inmax-5,m_blk=jnmax-5, &
         idimy_blk=jnmax-5,nww=50000
    real an(n_blk),bn(n_blk),cn(n_blk)
    real am(m_blk),bm(m_blk),cm(m_blk)
    real yy(m_blk,n_blk),ww(nww)
    save ww, an, bn, cn

    ! local fft variables
    integer, parameter :: nrlpts=knmax-5
    integer kint, kp
    real phk
    real rldat(nrlpts)
    real wsave(2*nrlpts+15)
    save wsave
    real qq(in-5,jn-5,nrlpts)
    real srhsb(jn-5,nrlpts),srhst(jn-5,nrlpts)
    real srhsl(in-5,nrlpts),srhsr(in-5,nrlpts)
    real tqq(inmax-5,jnmax-5,kns)
    real tsrhsb(jnmax-5,kns),tsrhst(jnmax-5,kns)
    real tsrhsl(inmax-5,kns),tsrhsr(inmax-5,kns)

    logical, save :: firstcall = .TRUE.
    real :: errorslc(2), errors(2)

    integer i,j,k,ip1,jp1,kp1,it,ierr

    real divdels(in,jn,kn)
    real rsd(in,jn,kn)

    real rnormlc,rnorm
    real pmeanlc(inmax),pmeanarr(inmax)
    real s1meanlc(inmax),s1meanarr(inmax)

    !---------------------------------------------------
    !
    !      initialize fft and blktri
    !
    if (firstcall) then

       call rffti(nrlpts,wsave)
       kk=log(dble(n_blk))/log(dble(2))+1
       ll=2**(kk+1)
       nww_ch=(kk-2)*ll+kk+5+max(2*n_blk,6*m_blk)
       if(nww .lt. nww_ch) then
          write(6,*) 'nww is set to', nww
          write(6,*) 'nww should be',nww_ch
          call MPI_ABORT(MPI_COMM_WORLD, 1,ierr)
       endif

       do i=1,n_blk
          an(i)=ar(i)
          bn(i)=br(i)
          cn(i)=cr(i)
       enddo
       an(1)=0.D0
       bn(1)=bn(1)+ar(1)*prbt
       cn(n_blk)=0.D0
       bn(n_blk)=bn(n_blk)+cr(n_blk)*prtp
       iflg_blk=0
       !       write(6,*) 'initializing blktri'
       call blktri(iflg_blk,np_blk,n_blk,an,bn,cn,mp_blk, &
            m_blk,am,bm,cm,idimy_blk,yy,ierr_blk,ww)
       !        write(6,*) 'return from blktri'
       !        write(6,*) 'ierr_blk,ww(1)=',ierr_blk,ww(1)
       firstcall = .FALSE.

    endif
    !
    !      compute divergence
    !
    do k=ks,ke
       kp1 = k + kone
       do j=js,je
          jp1 = j + jone
          do i=is,ie
             ip1        = i + ione
             divdels(i,j,k) = &
                  ( g2a(ip1) * g31a(ip1) * dels1(ip1,j,k) &
                  - g2a(i  ) * g31a(i  ) * dels1(i  ,j,k) ) &
                  * g2bi(i) * g31bi(i) * dx1ai(i) &
                  + ( g32a(jp1) * dels2(i,jp1,k) &
                  - g32a(j  ) * dels2(i,j  ,k) ) &
                  * g2bi(i) * g32bi(j) * dx2ai(j) &
                  + ( dels3(i,j,kp1) - dels3(i,j,k) ) &
                  * g31bi(i) * g32bi(j) * dx3ai(k)
             divdels(i,j,k)=divdels(i,j,k) &
                  *g2b(i)**2*fact(myid1*(in-5)+i)
          enddo
       enddo
    enddo
    !
    !      remove the horizontally averaged component
    !
    rnormlc=0.D0
    if(myid1 .eq. 0) then
       do k=ks,ke
          do j=js,je
             rnormlc=rnormlc &
                  +g32b(j)*dx2a(j)*dx3a(k)
          enddo
       enddo
    endif
    call MPI_ALLREDUCE(rnormlc,rnorm,1, &
         MPI_DOUBLE_PRECISION,MPI_SUM, &
         MPI_COMM_WORLD,ierr)

    do i=1,inmax
       s1meanlc(i)=0.D0
    enddo

    do i=is,ie
       do k=ks,ke
          do j=js,je
             s1meanlc(i+myid1*(in-5))=s1meanlc(i+myid1*(in-5)) &
                  +dels1(i,j,k)*g32b(j)*dx2a(j)*dx3a(k)
          enddo
       enddo
    enddo
    if(myid1 .eq. nproc1-1) then
       do k=ks,ke
          do j=js,je
             s1meanlc(iep1+myid1*(in-5))=s1meanlc(iep1+myid1*(in-5)) &
                  +dels1(iep1,j,k)*g32b(j)*dx2a(j)*dx3a(k)
          enddo
       enddo
    endif
    call MPI_ALLREDUCE(s1meanlc(is:inmax-2),s1meanarr(is:inmax-2),inmax-4, &
         MPI_DOUBLE_PRECISION,MPI_SUM, &
         MPI_COMM_WORLD,ierr)
    do i=is,inmax-2
       s1meanarr(i)=s1meanarr(i)/rnorm
    enddo
    do k=ks,ke
       do j=js,je
          do i=is,ie
             ip1        = i + ione
             rsd(i,j,k)=divdels(i,j,k) &
                  - (g2a(ip1) * g31a(ip1) * s1meanarr(ip1+myid1*(in-5)) &
                  - g2a(i  ) * g31a(i  ) * s1meanarr(i+myid1*(in-5))) &
                  * dx1ai(i) * fact(myid1*(in-5)+i)
          enddo
       enddo
    enddo
    !
    !      blktri
    !
    !      fft right hand side
    do j=js,je
       do i=is,ie
          do k=ks,ke
             rldat(k-ks+1)=rsd(i,j,k)
          enddo
          call rfftf(nrlpts,rldat,wsave) 
          do k=1,nrlpts
             qq(i-is+1,j-js+1,k)=rldat(k)/dble(nrlpts)
          enddo
       enddo
       do k=ks,ke
          rldat(k-ks+1)=dels1(is,j,k)-s1meanarr(is+myid1*(in-5))
       enddo
       call rfftf(nrlpts,rldat,wsave)
       do k=1,nrlpts
          srhsb(j-js+1,k)=rldat(k)/dble(nrlpts)
       enddo
       do k=ks,ke
          rldat(k-ks+1)=dels1(iep1,j,k)-s1meanarr(iep1+myid1*(in-5))
       enddo
       call rfftf(nrlpts,rldat,wsave)
       do k=1,nrlpts
          srhst(j-js+1,k)=rldat(k)/dble(nrlpts)
       enddo
    enddo
    do i=is,ie
       do k=ks,ke
          rldat(k-ks+1)=-fact(i+myid1*(in-5))*g2b(i)*dx2b(js) &
               *dels2(i,js,k)
       enddo
       call rfftf(nrlpts,rldat,wsave)
       do k=1,nrlpts
          srhsl(i-is+1,k)=rldat(k)/dble(nrlpts)
       enddo
    enddo
    do i=is,ie
       do k=ks,ke
          rldat(k-ks+1)=fact(i+myid1*(in-5))*g2b(i)*dx2b(jep1) &
               *dels2(i,jep1,k)
       enddo
       call rfftf(nrlpts,rldat,wsave)
       do k=1,nrlpts
          srhsr(i-is+1,k)=rldat(k)/dble(nrlpts)
       enddo
    enddo
    !
    !      transpose
    call transposef(qq,srhsb,srhst,srhsl,srhsr, &
         tqq,tsrhsb,tsrhst,tsrhsl,tsrhsr)
    !      for each of the fourier components solve blktri
    do k=1,kns
       kp=mod(k+myid*kns-1,knmax-5)+1
       kint=kp/2
       phk=2.D0*pi*dble(kint)/dble(nrlpts)
       do j=1,m_blk
          am(j)=ath(j)
          bm(j)=bth(j) &
               -(2.D0-2.D0*cos(phk))*bthk(j)
          cm(j)=cth(j)
       enddo
       am(1)=0.D0
       bm(1)=bm(1)+ath(1)
       cm(m_blk)=0.D0
       bm(m_blk)=bm(m_blk)+cth(m_blk)

       do j=js,jnmax-3
          do i=is,inmax-3
             yy(j-js+1,i-is+1)=tqq(i-is+1,j-js+1,k)
          enddo
       enddo
       do j=js,jnmax-3
          yy(j-js+1,1)=yy(j-js+1,1) &
               -ar(1)*srbt*tsrhsb(j-js+1,k)
          yy(j-js+1,inmax-5)=yy(j-js+1,inmax-5) &
               -cr(inmax-5)*srtp*tsrhst(j-js+1,k)
       enddo
       do i=is,inmax-3
          yy(1,i-is+1)=yy(1,i-is+1) &
               -ath(1)*tsrhsl(i-is+1,k)
          yy(jnmax-5,i-is+1)=yy(jnmax-5,i-is+1) &
               -cth(jnmax-5)*tsrhsr(i-is+1,k)
       enddo
       iflg_blk=1
       !         write(6,*) 'blktri solve'
       call blktri(iflg_blk,np_blk,n_blk,an,bn,cn,mp_blk, &
            m_blk,am,bm,cm,idimy_blk,yy,ierr_blk,ww)
       !         write(6,*) 'return from blktri'
       !         write(6,*) 'kint,ierr_blk,ww(1)=',kint,ierr_blk,ww(1)
       !
       do j=js,jnmax-3
          do i=is,inmax-3
             tqq(i-is+1,j-js+1,k)=yy(j-js+1,i-is+1)
          enddo
       enddo
       !
    enddo
    !
    !      inverse transpose
    call rtranspose(tqq,qq)
    !      inverse fft
    !
    do j=js,je
       do i=is,ie
          do k=1,nrlpts
             rldat(k)=qq(i-is+1,j-js+1,k)
          enddo
          call rfftb(nrlpts,rldat,wsave)
          do k=1,nrlpts
             p(i,j,k-1+ks)=rldat(k)
          enddo
       enddo
    enddo
    !
    !      remove the horizontal averaged component
    !
    do i=1,inmax
       pmeanlc(i)=0.D0
    enddo

    do i=is,ie
       do k=ks,ke
          do j=js,je
             pmeanlc(i+myid1*(in-5))=pmeanlc(i+myid1*(in-5)) &
                  +p(i,j,k)*g32b(j)*dx2a(j)*dx3a(k)
          enddo
       enddo
    enddo
    call MPI_ALLREDUCE(pmeanlc(is:inmax-3),pmeanarr(is:inmax-3),inmax-5, &
         MPI_DOUBLE_PRECISION,MPI_SUM, &
         MPI_COMM_WORLD,ierr)
    do i=is,inmax-3
       pmeanarr(i)=pmeanarr(i)/rnorm
    enddo
    do k=ks,ke
       do j=js,je
          do i=is,ie
             p(i,j,k)=p(i,j,k) &
                  -pmeanarr(i+myid1*(in-5))
          enddo
       enddo
    enddo
    !
    !      compute and add the horizontally averaged component
    !
    pmeanarr(is-1)=srbt*s1meanarr(is)*0.5D0
    do i=is,inmax-2
       pmeanarr(i)=pmeanarr(i-1) &
            +2.D0/((1.D0/fact(i-1)+1.D0/fact(i))*ovdxxb(i))*s1meanarr(i)
    enddo
    do k=ks,ke
       do j=js,je
          do i=is,ie
             p(i,j,k)=p(i,j,k)+pmeanarr(i+myid1*(in-5))
          enddo
       enddo
    enddo
    !
    call bvalp

    call lapl(p,rsd)
    !
    !      ref error and error
    !
    if((ErrPSolv).and.(mod(itnow,itintv)==0))then
       call lapl(p,rsd)
       errorslc = 0.D0
       errorslc(1) = sum(abs(divdels(is:ie,js:je,ks:ke)))
       errorslc(2) = sum(abs(rsd(is:ie,js:je,ks:ke)-divdels(is:ie,js:je,ks:ke)))
       call MPI_ALLREDUCE(errorslc, errors,2, MPI_DOUBLE_PRECISION,MPI_SUM, &
            MPI_COMM_WORLD,ierr)
       if(myid==0) then
          write(6,'(2(a,e23.15))') 'ref error=', errors(1), ', error=', errors(2)
          call flush(6)
       endif
    endif

  end subroutine getp

  !----------------------------------------------------------
  subroutine lapl(p,lp)
    use ModPar
    use ModGrid
    use ModBack,     ONLY: fact
    implicit none

    real, intent(in) :: p(in,jn,kn)
    real, intent(out) :: lp(in,jn,kn)
    real sp1up,sp1dn,sp2up,sp2dn,sp3up,sp3dn
    integer i,j,k,ip1,jp1,kp1
    !
    do k=ks,ke
       kp1 = k + kone
       do j=js,je
          jp1 = j + jone
          do i=is,ie
             ip1        = i + ione
             sp1up = (p(i+1,j,k) - p(i,j,k)) &
                  *dx1bi(i+1) &
                  *0.5D0*(1.D0/fact(myid1*(in-5)+i+1) &
                  +1.D0/fact(myid1*(in-5)+i)) &
                  *g2a(i+1)*g31a(i+1)
             sp1dn = (p(i,j,k) - p(i-1,j,k)) &
                  *dx1bi(i) &
                  *0.5D0*(1.D0/fact(myid1*(in-5)+i) &
                  +1.D0/fact(myid1*(in-5)+i-1)) &
                  *g2a(i)*g31a(i)
             sp2up = (p(i,j+1,k) - p(i,j,k)) &
                  *dx2bi(j+1)*g32a(j+1)
             sp2dn = (p(i,j,k) - p(i,j-1,k)) &
                  *dx2bi(j)*g32a(j)
             sp3up = (p(i,j,k+1) - p(i,j,k))*dx3bi(k+1)
             sp3dn = (p(i,j,k) - p(i,j,k-1))*dx3bi(k)
             lp(i,j,k) = ( sp1up - sp1dn ) &
                  * fact(myid1*(in-5)+i) * dx1ai(i) &
                  + ( sp2up - sp2dn ) &
                  * g32bi(j) * dx2ai(j) &
                  + ( sp3up - sp3dn ) &
                  * dx3ai(k) * g32bi(j) * g32bi(j)
          enddo
       enddo
    enddo
    !
  end subroutine lapl
  !----------------------------------------------------------

  subroutine transposef(qq,srhsb,srhst,srhsl,srhsr,tqq,tsrhsb,tsrhst,tsrhsl,tsrhsr)
    use ModPar
    use ModMpi
    implicit none

    real, intent(in) :: qq(in-5,jn-5,knmax-5)
    real, intent(in), dimension(jn-5,knmax-5) :: srhsb, srhst
    real, intent(in), dimension(in-5,knmax-5) :: srhsl, srhsr
    real, intent(out) :: tqq(inmax-5,jnmax-5,kns)
    real, intent(out), dimension(jnmax-5,kns) :: tsrhsb, tsrhst
    real, intent(out), dimension(inmax-5,kns) :: tsrhsl, tsrhsr
    integer :: i, j, k, ip, isend, irecv, myid1_recv, myid2_recv, &
         i1, i2, j1, j2, k1, k2

    integer, parameter :: ncount=(in-5+2)*(jn-5+2)*kns
    integer :: ireq(2*(nproc-1)), &
         istatus_arr(MPI_STATUS_SIZE,2*(nproc-1)), ierr

    real, allocatable :: buf(:,:,:), buft(:,:,:), sendbf(:), recvbf(:)
    !-----------------------------------------------------------------------
    ! copy dianonal data which remain local
    i1 = myid1*(in-5) + 1
    i2 = myid1*(in-5) + in - 5
    j1 = myid2*(jn-5) + 1
    j2 = myid2*(jn-5) + jn - 5
    k1 = myid*kns + 1
    k2 = myid*kns + kns
    tqq(i1:i2,j1:j2,:)=qq(:,:,k1:k2)
    tsrhsb(j1:j2,:)    = srhsb(:,k1:k2)
    tsrhst(j1:j2,:)    = srhst(:,k1:k2)
    tsrhsl(i1:i2,:)    = srhsl(:,k1:k2)
    tsrhsr(i1:i2,:)    = srhsr(:,k1:k2)

    ! swap data with all other CPUS
    allocate(buf(0:in-5+1,0:jn-5+1,knmax-5))
    buf = 0.D0
    buf(1:in-5,1:jn-5,:) = qq
    ! if boundary processors, fill the top/bottom row or right/left column for boundary conditions
    if(myid1==0) buf(0,1:jn-5,:) = srhsb
    if(myid1==(nproc1-1)) buf(in-5+1,1:jn-5,:) = srhst
    if(myid2==0) buf(1:in-5,0,:) = srhsl
    if(myid2==(nproc2-1)) buf(1:in-5,jn-5+1,:) = srhsr

    allocate(sendbf(ncount*(nproc-1)), recvbf(ncount*(nproc-1)))
    do ip=1,nproc-1
       isend = mod(nproc-ip+myid, nproc) ! from myid-1 to myid+1
       sendbf((ip-1)*ncount+1:ip*ncount) = &
            reshape(buf(:,:,isend*kns+1:isend*kns+kns),(/ncount/))
    enddo
    deallocate(buf)

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    do ip=1,nproc-1
       isend = mod(nproc-ip+myid, nproc) ! from myid-1 to myid+1
       irecv = mod(myid+ip, nproc) ! from myid+1 to myid-1
       call MPI_IRECV(recvbf(1+(ip-1)*ncount), ncount, MPI_DOUBLE_PRECISION, irecv, &
            ip, MPI_COMM_WORLD, ireq(ip+nproc-1), ierr)
       call MPI_ISEND(sendbf(1+(ip-1)*ncount), ncount, MPI_DOUBLE_PRECISION, isend, &
            ip, MPI_COMM_WORLD, ireq(ip), ierr)
    enddo
    call MPI_WAITALL(2*(nproc-1),ireq,istatus_arr,ierr)

    allocate(buft(0:in-5+1,0:jn-5+1,1:kns))
    do ip=1,nproc-1
       irecv = mod(myid+ip, nproc)  ! from myid +1 to myid-1
       myid1_recv=mod(irecv,nproc1)
       myid2_recv=irecv/nproc1
       i1 = myid1_recv*(in-5) + 1
       i2 = myid1_recv*(in-5) + in - 5
       j1 = myid2_recv*(jn-5) + 1
       j2 = myid2_recv*(jn-5) + jn - 5
       buft = reshape(recvbf((ip-1)*ncount+1:ip*ncount),(/in-5+2,jn-5+2,kns/))
       tqq(i1:i2,j1:j2,1:kns) = buft(1:in-5,1:jn-5,:)
       if (myid1_recv==0) tsrhsb(j1:j2,:) = buft(0,1:jn-5,:)
       if (myid1_recv==(nproc1-1)) tsrhst(j1:j2,:) = buft(in-5+1,1:jn-5,:)
       if (myid2_recv==0) tsrhsl(i1:i2,:) = buft(1:in-5,0,:)
       if (myid2_recv==(nproc2-1)) tsrhsr(i1:i2,:) = buft(1:in-5,jn-5+1,:)
    enddo
    deallocate(sendbf, recvbf, buft)

  end subroutine transposef
  !-------------------------------------------------------------
  subroutine rtranspose(tqq,qq)
    use ModPar
    use ModMpi
    implicit none

    real, intent(out) :: qq(in-5,jn-5,knmax-5)
    real, intent(in) :: tqq(inmax-5,jnmax-5,kns)

    integer, parameter :: ncount=(in-5)*(jn-5)*kns
    real, allocatable, dimension(:) :: sendbf, recvbf
    integer :: i1, i2, j1, j2, k1, k2, ip, &
         isend,irecv,myid1_send,myid2_send
    integer :: ireq(2*(nproc-1)), &
         istatus_arr(MPI_STATUS_SIZE,2*(nproc-1)), ierr
    !-------------------------------------------------------------
    !
    ! copy dianonal data which remain local
    !
    i1 = myid1*(in-5) + 1
    i2 = myid1*(in-5) + in - 5
    j1 = myid2*(jn-5) + 1
    j2 = myid2*(jn-5) + jn - 5
    k1 = myid*kns + 1
    k2 = myid*kns + kns
    qq(:,:,k1:k2) = tqq(i1:i2,j1:j2,:)
    !
    ! swap data with all other CPUS
    !
    allocate(sendbf(ncount*(nproc-1)), recvbf(ncount*(nproc-1)))
    do ip=1,nproc-1
       isend=mod(nproc-ip+myid,nproc) ! from myid-1 to myid+1
       myid1_send=mod(isend,nproc1)
       myid2_send=isend/nproc1
       i1 = myid1_send*(in-5) + 1
       i2 = myid1_send*(in-5) + in - 5
       j1 = myid2_send*(jn-5) + 1
       j2 = myid2_send*(jn-5) + jn - 5
       k1 = myid*kns + 1
       k2 = myid*kns + kns
       sendbf((ip-1)*ncount+1:ip*ncount) = &
            reshape(tqq(i1:i2,j1:j2,1:kns),(/ncount/))
    enddo

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    do ip=1,nproc-1
       isend=mod(nproc-ip+myid,nproc)
       irecv=mod(myid+ip,nproc)
       call MPI_IRECV(recvbf(1+(ip-1)*ncount),ncount, &
            MPI_DOUBLE_PRECISION,irecv,ip, &
            MPI_COMM_WORLD,ireq(ip+nproc-1),ierr)
       call MPI_ISEND(sendbf(1+(ip-1)*ncount),ncount, &
            MPI_DOUBLE_PRECISION,isend,ip, &
            MPI_COMM_WORLD,ireq(ip),ierr)
    enddo
    call MPI_WAITALL(2*(nproc-1),ireq,istatus_arr,ierr)

    do ip=1,nproc-1
       irecv=mod(myid+ip,nproc) ! from myid+1 to myid-1
       qq(:,:,irecv*kns+1:irecv*kns+kns) = &
            reshape(recvbf((ip-1)*ncount+1:ip*ncount),(/in-5,jn-5,kns/))
    enddo
    deallocate(sendbf, recvbf)

  end subroutine rtranspose

end module ModGetP
