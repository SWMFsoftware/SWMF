!This code is a copyright protected software (c) 2002- University of Michigan

!=============================================================================
subroutine ray_pass

!  call ray_pass_new
  call ray_pass_old

end subroutine ray_pass
!=============================================================================
subroutine ray_pass_new

  use ModMain, ONLY : nBlock,Unused_B
  use ModParallel, ONLY : neiLEV
  use ModRaytrace, ONLY: RayFace
  use BATL_lib, ONLY: message_pass_node

  implicit none

  integer :: iBLK, iface

!  do i=1,3; do j=1,2

!!     call pass_and_max_nodes(.false.,rayface(i,j,:,:,:,:))
!!$     call pass_and_max_nodes(.true.,rayface(i,j,:,:,:,:))
!  end do; end do

  call message_pass_node(6, RayFace, 'max')

  do iBLK=1,nBlock
     if(Unused_B(iBLK))cycle
     do iface=1,6
        if(neiLEV(iface,iBLK)==1)call prolong_ray_after_pass(iface,iBLK)
     end do
  end do

end subroutine ray_pass_new

!===========================================================================
subroutine prolong_ray_after_pass(iface,iBLK)
  use ModRaytrace
  implicit none

  ! For faces that are shared with a coarser neighbor, interpolate 
  ! for all points which are not coinciding and where the ray is going out.
  !
  ! a at odd  j and even k requires interpolation in direction k
  ! b at even j and odd  k requires interpolation in direction j
  ! c at even j and even k requires interpolation in both directions

  !(1,5)           (5,5)
  !   O-- b --O-- b --O
  !   |   |   |   |   |
  !   |   |   |   |   |
  !   a - c - a - c - a
  !   |   |   |   |   |
  !   |   |   |   |   |
  !   O-- b --O-- b --O
  !   |   |   |   |   |
  !   |   |   |   |   |
  !   a - c - a - c - a
  !   |   |   |   |   |
  !   |   |   |   |   |
  !   O-- b --O-- b --O
  !(1,1)           (5,1)

  integer, intent(in) :: iface,iBLK
  integer :: iray
  integer :: j, k, nFaceJ, nFaceK
  integer, parameter :: nFaceMax=max(nI+1,nJ+1,nK+1)
  real    :: qrayface(3,2,nFaceMax,nFaceMax)
  integer :: qrayend_ind(2,nFaceMax,nFaceMax)

  ! Interpolation weights
  real, dimension(4), parameter:: weight4=0.25
  real, dimension(2), parameter:: weight2=0.5

  !-------------------------------------------------------------------------

  ! Extract qrayface and qrayend_ind for the appropriate face
  !NOTE: qrayend_ind assignment split to two lines to avoid reshaping compiler bug!
  select case(iface)
  case(1)
     nFaceJ=nJ+1; nFaceK=nK+1
     qrayface( :,:,1:nFaceJ,1:nFaceK)=rayface(   :,:,1,1:nJ+1,1:nK+1,iBLK)
     qrayend_ind(1,1:nFaceJ,1:nFaceK)=rayend_ind(1,1,1,1:nJ+1,1:nK+1,iBLK)
     qrayend_ind(2,1:nFaceJ,1:nFaceK)=rayend_ind(1,2,1,1:nJ+1,1:nK+1,iBLK)
  case(2)
     nFaceJ=nJ+1; nFaceK=nK+1
     qrayface( :,:,1:nFaceJ,1:nFaceK)=rayface(   :,:,nI+1,1:nJ+1,1:nK+1,iBLK)
     qrayend_ind(1,1:nFaceJ,1:nFaceK)=rayend_ind(1,1,nI+1,1:nJ+1,1:nK+1,iBLK)
     qrayend_ind(2,1:nFaceJ,1:nFaceK)=rayend_ind(1,2,nI+1,1:nJ+1,1:nK+1,iBLK)
  case(3)
     nFaceJ=nI+1; nFaceK=nK+1
     qrayface( :,:,1:nFaceJ,1:nFaceK)=rayface(   :,:,1:nI+1,1,1:nK+1,iBLK)
     qrayend_ind(1,1:nFaceJ,1:nFaceK)=rayend_ind(1,1,1:nI+1,1,1:nK+1,iBLK)
     qrayend_ind(2,1:nFaceJ,1:nFaceK)=rayend_ind(1,2,1:nI+1,1,1:nK+1,iBLK)
  case(4)
     nFaceJ=nI+1; nFaceK=nK+1
     qrayface( :,:,1:nFaceJ,1:nFaceK)=rayface(   :,:,1:nI+1,nJ+1,1:nK+1,iBLK)
     qrayend_ind(1,1:nFaceJ,1:nFaceK)=rayend_ind(1,1,1:nI+1,nJ+1,1:nK+1,iBLK)
     qrayend_ind(2,1:nFaceJ,1:nFaceK)=rayend_ind(1,2,1:nI+1,nJ+1,1:nK+1,iBLK)
  case(5)
     nFaceJ=nI+1; nFaceK=nJ+1
     qrayface( :,:,1:nFaceJ,1:nFaceK)=rayface(   :,:,1:nI+1,1:nJ+1,1,iBLK)
     qrayend_ind(1,1:nFaceJ,1:nFaceK)=rayend_ind(1,1,1:nI+1,1:nJ+1,1,iBLK)
     qrayend_ind(2,1:nFaceJ,1:nFaceK)=rayend_ind(1,2,1:nI+1,1:nJ+1,1,iBLK)
  case(6)
     nFaceJ=nI+1; nFaceK=nJ+1
     qrayface( :,:,1:nFaceJ,1:nFaceK)=rayface(   :,:,1:nI+1,1:nJ+1,nK+1,iBLK)
     qrayend_ind(1,1:nFaceJ,1:nFaceK)=rayend_ind(1,1,1:nI+1,1:nJ+1,nK+1,iBLK)
     qrayend_ind(2,1:nFaceJ,1:nFaceK)=rayend_ind(1,2,1:nI+1,1:nJ+1,nK+1,iBLK)
  case default
     call stop_mpi('Impossible value for iface in prolong_ray')
  end select

  do iray=1,2
     do k=1,nfaceK
        if(mod(k,2)==1)then
           do j=2,nfaceJ,2
              ! Case b: even j and odd k

              if(qrayend_ind(iray,j,k)/=ray_out_)CYCLE

              call rayface_interpolate(&
                   qrayface(:,iray,j-1:j+1:2,k),weight2,2,&
                   qrayface(:,iray,j,k))
           end do
        else
           do j=1,nJ+1,2
              ! Case a: odd j and even k

              if(qrayend_ind(iray,j,k)/=ray_out_)CYCLE

              call rayface_interpolate(&
                   qrayface(:,iray,j,k-1:k+1:2),weight2,2,&
                   qrayface(:,iray,j,k))
           end do
           do j=2,nJ,2
              ! Case c: even j and even k

              if(qrayend_ind(iray,j,k)/=ray_out_)CYCLE

              call rayface_interpolate(&
                   qrayface(:,iray,j-1:j+1:2,k-1:k+1:2),weight4,4,&
                   qrayface(:,iray,j,k))
           end do ! j
        end if ! mod(k,2)
     end do ! k
  end do ! iray

  ! Put back result into rayface
  select case(iface)
  case(1)
     rayface(:,:,     1,1:nJ+1,1:nK+1,iBLK)=qrayface(:,:,1:nFaceJ,1:nFaceK)
  case(2)
     rayface(:,:,  nI+1,1:nJ+1,1:nK+1,iBLK)=qrayface(:,:,1:nFaceJ,1:nFaceK)
  case(3)
     rayface(:,:,1:nI+1,     1,1:nK+1,iBLK)=qrayface(:,:,1:nFaceJ,1:nFaceK)
  case(4)
     rayface(:,:,1:nI+1,  nJ+1,1:nK+1,iBLK)=qrayface(:,:,1:nFaceJ,1:nFaceK)
  case(5)
     rayface(:,:,1:nI+1,1:nJ+1,     1,iBLK)=qrayface(:,:,1:nFaceJ,1:nFaceK)
  case(6)
     rayface(:,:,1:nI+1,1:nJ+1,  nK+1,iBLK)=qrayface(:,:,1:nFaceJ,1:nFaceK)
  end select

end subroutine prolong_ray_after_pass

!==================================================================================
subroutine ray_pass_old

  ! Exchange and update rayface values between blocks direction by direction

  ! Notation: _o out        (cells to be sent for equal blocks) 
  !           _g get        (cells to be received)
  !           _r restricted (to be sent to a coarser block)
  !           _s subface    (one quarter of a face)

  use ModProcMH
  use ModMain, ONLY : nblockMax,okdebug,Unused_B,optimize_message_pass
  use BATL_lib, ONLY: iNode_B, iTree_IA, Coord0_
  use ModRaytrace
  use ModParallel, ONLY : NOBLK,neiLEV,neiBLK,neiPE
  use ModMpi
  implicit none

  ! Local variables

  ! idir=1,2,3 correspond to east-west, south-north, bot-top.
  integer :: idir, isweep

  ! Face (east..top), side (1 for east,south,bot, 2 for others)
  integer :: iface, otherface, iside

  ! number of subfaces (1 or 4), subface (1..nsubface) and child (1..8) index
  integer ::  nsubface, isubface

  ! Array ranges for outgoing, incoming, restricted and subfaces
  integer :: imin_o,imax_o,jmin_o,jmax_o,kmin_o,kmax_o
  integer :: imin_g,imax_g,jmin_g,jmax_g,kmin_g,kmax_g
  integer :: imin_r,imax_r,jmin_r,jmax_r,kmin_r,kmax_r
  integer :: imin_s,imax_s,jmin_s,jmax_s,kmin_s,kmax_s

  ! Block index (1..nBLK)
  integer :: iBLK

  ! Descriptors for neighbor
  integer :: neiP,neiB,neiL

  ! MPI variables
  integer :: itag, request, number_receive_requests, receive_requests(nBLK*6)
  integer :: status(MPI_STATUS_SIZE, nBLK*6)  

  ! Maximum size of the RESTRICTED rayface layer to be received
  ! for the 6 ray variables (3 coord*2 ray dir.)
  integer, parameter :: maxsize_r= 6*max((nI/2+1)*(nJ/2+1),(nI/2+1)*(nK/2+1),(nJ/2+1)*(nK/2+1))

  ! Receive buffer to hold 4 incoming RESTRICTED rayface values
  ! for all blocks and for both sides
  real, dimension(maxsize_r,4,nBLK,2) :: buffer

  ! Actual size of messages: full, restricted/sparse and actual face
  integer :: isize, isize_r, isize1

  ! Equal and restricted values to be sent are stored in these buffers
  real, dimension(:,:,:,:,:), allocatable :: eq_buf, re_buf

  logical :: oktest, oktest_me

  integer :: iError

  !---------------------------------------------------------------------------

  call set_oktest('ray_pass',oktest, oktest_me)
  if(oktest)write(*,*)'ray_pass me=',iProc

  do idir=1,3
     select case(optimize_message_pass)
     case('face')
        ! Send messages face by face
        call ray_pass_faces(2*idir-1,2*idir-1,.true.,.true.,.true.)
        call ray_pass_faces(2*idir  ,2*idir  ,.true.,.true.,.true.)
     case('min')
        ! Send messages face by face and kind by kind
        do isweep=2*idir-1,2*idir
           ! Send equal
           call ray_pass_faces(isweep,isweep,.true.,.false.,.false.)
           ! Send restricted
           call ray_pass_faces(isweep,isweep,.false.,.true.,.false.)
           ! Send prolonged
           call ray_pass_faces(isweep,isweep,.false.,.false.,.true.)
        end do
     case default
        ! Send messages for both faces
        call ray_pass_faces(2*idir-1,2*idir,.true.,.true.,.true.)
     end select
  end do ! idir

  if(oktest_me)write(*,*)'ray_pass starting prolongation'

  do iBLK=1,nBlockMax
     if(Unused_B(iBLK))cycle

     do iface=1,6
        if(neiLEV(iface,iBLK)==1)call prolong_ray
     end do
  end do

  if(oktest_me)write(*,*)'ray_pass finished'

contains

  subroutine ray_pass_faces(&
       ifacemin,ifacemax,do_equal,do_restricted,do_prolonged)

    integer, intent(in):: ifacemin,ifacemax
    logical, intent(in):: do_equal,do_restricted,do_prolonged

    ! BATL related
    integer:: iNode, iDim, iSideFace
    !------------------------------------------------------------------------

    if(oktest)write(*,*)&
         'ray_pass_faces:me,ifacemin,ifacemax,do_eq,do_re,do_pr=',&
         iProc,ifacemin,ifacemax,do_equal,do_restricted,do_prolonged

    ! Debug
    if(okdebug)buffer  =0.00

    number_receive_requests = 0
    receive_requests = MPI_REQUEST_NULL

    do iface=ifacemin,ifacemax

       ! Set index ranges for the face
       call setranges_ray

       if(okdebug.and.oktest)then
          write(*,*)&
               'setranges_ray for receive done: me,iface,isize,isize_r',&
               iProc, iface, isize, isize_r
          write(*,*)'_o=',imin_o,imax_o,jmin_o,jmax_o,kmin_o,kmax_o
          write(*,*)'_g=',imin_g,imax_g,jmin_g,jmax_g,kmin_g,kmax_g
          write(*,*)'_r=',imin_r,imax_r,jmin_r,jmax_r,kmin_r,kmax_r
       end if

       do iBLK = 1,nBlockMax
          if(Unused_B(iBLK))CYCLE
          ! Post non-blocking receive for opposite face of neighbor block
          neiL=neiLEV(otherface,iBLK)
          select case(neiL)
          case(0)
             if(.not.do_equal)cycle
             nsubface=1
             isize1=isize
          case(1)
             if(.not.do_prolonged)cycle
             nsubface=1
             isize1=isize_r
          case(-1)
             if(.not.do_restricted)cycle
             nsubface=4
             isize1=isize_r
          case(NOBLK)
             ! Do nothing
             CYCLE
          case default
             write(*,*)'me,iBLK,otherface,neiL=',&
                  iProc,iBLK,otherface,neiL
             call stop_mpi(&
                  'Error in message pass: Invalid value for neiLEV')
          end select

          if(okdebug.and.oktest_me)write(*,*)&
               'receive: me,neiL,nsubface,isize1',&
               iProc,neiL,nsubface,isize1

          do isubface=1,nsubface
             neiP=neiPE(isubface,otherface,iBLK)
             if(neiP/=iProc)then
                ! Remote receive
                itag = 100*iBLK+10*iface+isubface
                if(oktest.and.okdebug)write(*,*)&
                     'Remote receive, me,itag,neiL,neiP=',&
                     iProc,itag,neiL,neiP

                call MPI_irecv(buffer(1,isubface,iBLK,iside),&
                     isize1,MPI_REAL,neiP,itag,iComm,request,iError)
                number_receive_requests = number_receive_requests + 1
                receive_requests(number_receive_requests) = request
             end if
          end do ! isubface
       end do ! iBLK
    end do ! iface

    !\
    ! Wait for all receive commands to be posted for all processors
    !/
    call barrier_mpi

    if(oktest)write(*,*)'receives posted: me=',iProc

    !\
    ! Send blocking messages with Rsend (ready to receive)
    !/
    do iface=ifacemin,ifacemax

       ! Set index ranges for the face
       call setranges_ray

       if(okdebug.and.oktest)write(*,*)&
            'setranges_ray for send done: me, iface=',iProc, iface

       if(do_equal)&
            allocate(eq_buf(3,2,imin_o:imax_o,jmin_o:jmax_o,kmin_o:kmax_o))
       if(do_restricted.or.do_prolonged)&
            allocate(re_buf(3,2,imin_r:imax_r,jmin_r:jmax_r,kmin_r:kmax_r))

       if(okdebug.and.oktest)write(*,*)'allocation done, me,iface=',&
            iProc,iface

       do iBLK=1,nBlockMax
          if(Unused_B(iBLK))CYCLE
          neiL=neiLEV(iface,iBLK)

          if(okdebug.and.oktest)write(*,*)&
               'sending: me, iface,iBLK,neiL=',iProc,iface,iBLK,neiL
          select case(neiL)
          case(0)
             if(.not.do_equal)cycle

             neiP=neiPE(1,iface,iBLK)
             neiB=neiBLK(1,iface,iBLK)
             if(neiP==iProc)then
                ! Local copy
                if(okdebug.and.oktest)write(*,*)&
                     'local equal copy: me,iface,iBLK=',iProc,iface,iBLK

                ! Debug
                !write(*,*)'rayface(_o,iBLK)=',&
                !  rayface(:,:,imin_o:imax_o,jmin_o:jmax_o,kmin_o:kmax_o,iBLK)
                !write(*,*)'before: rayface(_g,neiB)=',&
                !  rayface(:,:,imin_g:imax_g,jmin_g:jmax_g,kmin_g:kmax_g,neiB)

                rayface(:,:,imin_g:imax_g,jmin_g:jmax_g,kmin_g:kmax_g,neiB)=&
                     max(&
                     rayface(:,:,imin_g:imax_g,jmin_g:jmax_g,kmin_g:kmax_g,neiB),&
                     rayface(:,:,imin_o:imax_o,jmin_o:jmax_o,kmin_o:kmax_o,iBLK))

                ! Debug
                !write(*,*)'after: rayface(_g,neiB)=',&
                !  rayface(:,:,imin_g:imax_g,jmin_g:jmax_g,kmin_g:kmax_g,neiB)

             else
                ! Remote send
                itag = 100*neiB+10*iface+1
                if(oktest.and.okdebug)write(*,*)&
                     'Remote equal send, me,itag,neiP=',iProc,itag,neiP

                eq_buf=&
                     rayface(:,:,imin_o:imax_o,jmin_o:jmax_o,kmin_o:kmax_o,iBLK)

                call MPI_Rsend(eq_buf,&
                     isize,MPI_REAL,neiP,itag,iComm,iError)
             end if
          case(1)
             if(.not.do_restricted)cycle

             ! Restrict rayface in _o range into _r
             re_buf=&
                  rayface(:,:,imin_o:imax_o:2,jmin_o:jmax_o:2,kmin_o:kmax_o:2,iBLK)

             neiP=neiPE(1,iface,iBLK)
             neiB=neiBLK(1,iface,iBLK)
             ! Subface index =1,2,3, or 4 with respect to the coarse neighbor
             
             ! iSubFace = iSubFace_IA(iFace,iNode_B(iBlk))
             iNode = iNode_B(iBlk)
             iSubFace = 0
             do iDim = 1, nDim
                iSideFace = modulo(iTree_IA(Coord0_+iDim,iNode) - 1,2)

                if(iDim == (iFace+1)/2) CYCLE

                if(iSubFace == 0) then
                   iSubFace = iSideFace + 1
                else
                   iSubFace = iSubFace + 2*iSideFace
                end if

             end do

             ! Swap subface 2 and 3 for iFace = 1..4 for BATSRUS tradition...
             if(iFace <= 4 .and. iSubFace >= 2 .and. iSubFace <= 3) &
                  iSubFace = 5 - iSubFace

             if(neiP==iProc)then
                ! Local copy into appropriate subface
                call setsubrange_ray(.false.)
                if(okdebug.and.oktest)write(*,*)&
                     'local restricted copy: me,iface,iBLK,_s=',&
                     iProc,iface,iBLK,&
                     imin_s,imax_s,jmin_s,jmax_s,kmin_s,kmax_s

                rayface(:,:,imin_s:imax_s,jmin_s:jmax_s,kmin_s:kmax_s,neiB)=&
                     max(re_buf,&
                     rayface(:,:,imin_s:imax_s,jmin_s:jmax_s,kmin_s:kmax_s,neiB))
             else
                ! Remote send
                itag = 100*neiB+10*iface+isubface
                if(oktest.and.okdebug)write(*,*)&
                     'Remote restricted send, me,iface,itag=',&
                     iProc,iface,itag
                call MPI_Rsend(re_buf,isize_r,&
                     MPI_REAL,neiP,itag,iComm,iError)
             end if
          case(-1)
             if(.not.do_prolonged)cycle

             do isubface=1,4
                neiP=neiPE(isubface,iface,iBLK)
                neiB=neiBLK(isubface,iface,iBLK)

                call setsubrange_ray(.true.)

                if(neiP==iProc)then
                   ! Local copy of appropriate subface

                   if(okdebug.and.oktest)write(*,*)&
                        'local prolonged copy: me,isubface,iface,iBLK,_s=',&
                        iProc,isubface,iface,iBLK,&
                        imin_s,imax_s,jmin_s,jmax_s,kmin_s,kmax_s

                   rayface(:,:,imin_g:imax_g:2,&
                               jmin_g:jmax_g:2,&
                               kmin_g:kmax_g:2,neiB)=max(&
                   rayface(:,:,imin_g:imax_g:2,&
                               jmin_g:jmax_g:2,&
                               kmin_g:kmax_g:2,neiB),&
                   rayface(:,:,imin_s:imax_s,jmin_s:jmax_s,kmin_s:kmax_s,iBLK))
                else
                   ! Remote send
                   re_buf=&
                        rayface(:,:,imin_s:imax_s,jmin_s:jmax_s,kmin_s:kmax_s,iBLK)

                   itag = 100*neiB+10*iface+1
                   if(oktest.and.okdebug)write(*,*)&
                        'Remote prolong send, me,iface,itag=',&
                        iProc,iface,itag

                   call MPI_Rsend(re_buf,isize_r,&
                        MPI_REAL,neiP,itag,iComm,iError)
                end if
             end do ! isubface

          case(NOBLK)
             ! There is no neighbor, do nothing
             CYCLE
          case default
             write(*,*)'me,iBLK,iface,neiL=',&
                  iProc,iBLK,iface,neiL
             call stop_mpi('Error in message pass: Invalid value for neiLEV')
          end select ! neiL
       end do ! iBLK

       if(do_equal)deallocate(eq_buf)
       if(do_restricted.or.do_prolonged)deallocate(re_buf)

       if(oktest_me)write(*,*)'messages sent, me, iface=',iProc,iface
    end do ! iface


    !\
    ! WAIT FOR ALL MESSAGES TO BE RECEIVED
    !/
    if (number_receive_requests > 0) &
         call MPI_waitall(number_receive_requests,receive_requests,status,iError)

    if(oktest_me)write(*,*)'messages received, me, idir=',iProc, idir

    ! Copy ghost cells received from non-local neigbors
    ! and stored in the buffer into sol_BLK

    do iface=ifacemin,ifacemax

       ! Set index ranges for the face
       call setranges_ray

       if(okdebug.and.oktest)write(*,*)&
            'setranges_ray for buf2ray done: me, iface=',iProc, iface

       do iBLK = 1,nBlockMax
          if(Unused_B(iBLK))CYCLE
          select case(neiLEV(otherface,iBLK))
          case(0)
             if(okdebug.and.oktest)&
                  write(*,*)'buf2rayface: me, iBLK=',iProc,iBLK
             if(do_equal.and.neiPE(1,otherface,iBLK)/=iProc)&
                  call buf2rayface(buffer(1,1,iBLK,iside),&
                  imin_g,imax_g,jmin_g,jmax_g,kmin_g,kmax_g)
          case(1)
             if(okdebug.and.oktest)&
                  write(*,*)'buf2sparserayface: me, iBLK=',iProc,iBLK
             if(do_prolonged.and.neiPE(1,otherface,iBLK)/=iProc)&
                  call buf2sparserayface(buffer(1,1,iBLK,iside),&
                  imin_r,imax_r,jmin_r,jmax_r,kmin_r,kmax_r)
          case(-1)
             if(do_restricted)then
                do isubface=1,4
                   if(okdebug.and.oktest)&
                        write(*,*)'buf2subrayface: me, isubface, iBLK=',&
                        iProc,isubface,iBLK
                   if(neiPE(isubface,otherface,iBLK)/=iProc)&
                        call buf2subrayface(&
                        buffer(1,isubface,iBLK,iside),&
                        imin_r,imax_r,jmin_r,jmax_r,kmin_r,kmax_r)
                end do
             end if
          end select ! neiL
       end do ! iBLK
    end do ! iface

    if(oktest)write(*,*)'ray_pass_faces finished: me, ifacemin, ifacemax=',&
         iProc,ifacemin,ifacemax


  end subroutine ray_pass_faces

  !===========================================================================

  subroutine setranges_ray

    ! Set ranges orthogonal to idir based on the value of idir

    if(idir/=1)then
       imin_g=1;  imax_g=nI+1
       imin_o=1;  imax_o=nI+1
       imin_r=1;  imax_r=nI/2+1
    end if

    if(idir/=2)then
       jmin_g=1;  jmax_g=nJ+1
       jmin_o=1;  jmax_o=nJ+1
       jmin_r=1;  jmax_r=nJ/2+1
    endif

    if(idir/=3)then
       kmin_g=1;  kmax_g=nK+1
       kmin_o=1;  kmax_o=nK+1
       kmin_r=1;  kmax_r=nK/2+1
    end if

    ! Set ranges in direction of idir based on the value of iface

    select case(iface)
    case(1)
       otherface=2; iside=1
       imin_o=1;    imax_o=1
       imin_g=nI+1; imax_g=nI+1
       imin_r=1;    imax_r=1
    case(3)
       otherface=4; iside=1
       jmin_o=1;    jmax_o=1
       jmin_g=nJ+1; jmax_g=nJ+1
       jmin_r=1;    jmax_r=1 
    case(5)
       otherface=6; iside=1
       kmin_o=1;    kmax_o=1
       kmin_g=nK+1; kmax_g=nK+1
       kmin_r=1;    kmax_r=1
    case(2)
       otherface=1; iside=2
       imin_o=nI+1; imax_o=nI+1
       imin_g=1;    imax_g=1
       imin_r=1;    imax_r=1
    case(4)
       otherface=3; iside=2
       jmin_o=nJ+1; jmax_o=nJ+1
       jmin_g=1;    jmax_g=1
       jmin_r=1;    jmax_r=1
    case(6)
       otherface=5; iside=2
       kmin_o=nK+1; kmax_o=nK+1
       kmin_g=1;    kmax_g=1
       kmin_r=1;    kmax_r=1
    end select

    ! Size of full and restricted cell layers
    isize  =6*(imax_g-imin_g+1)*(jmax_g-jmin_g+1)*(kmax_g-kmin_g+1)
    isize_r=6*(imax_r-imin_r+1)*(jmax_r-jmin_r+1)*(kmax_r-kmin_r+1)

  end subroutine setranges_ray

  !===========================================================================

  subroutine setsubrange_ray(oksend)

    logical, intent(in) :: oksend

    ! Select appropriate quarter of ghost cell layer

    select case(iface)
    case(1,2)
       if(oksend)then
          imin_s=imin_o; imax_s=imax_o
       else
          imin_s=imin_g; imax_s=imax_g
       end if

       select case(isubface)
          ! Beware, case(2) and case(3) are swapped
       case(1)
          jmin_s=jmin_r; jmax_s=jmax_r; 
          kmin_s=kmin_r; kmax_s=kmax_r
       case(3)
          jmin_s=jmin_r+nJ/2; jmax_s=jmax_r+nJ/2; 
          kmin_s=kmin_r; kmax_s=kmax_r
       case(2)
          jmin_s=jmin_r; jmax_s=jmax_r; 
          kmin_s=kmin_r+nK/2; kmax_s=kmax_r+nK/2;
       case(4)
          jmin_s=jmin_r+nJ/2; jmax_s=jmax_r+nJ/2; 
          kmin_s=kmin_r+nK/2; kmax_s=kmax_r+nK/2;
       end select
    case(3,4)
       if(oksend)then
          jmin_s=jmin_o; jmax_s=jmax_o
       else
          jmin_s=jmin_g; jmax_s=jmax_g
       end if
       select case(isubface)
          ! Beware, case(2) and case(3) are swapped
       case(1)
          imin_s=imin_r;      imax_s=imax_r; 
          kmin_s=kmin_r;      kmax_s=kmax_r
       case(3)
          imin_s=imin_r+nI/2; imax_s=imax_r+nI/2; 
          kmin_s=kmin_r;      kmax_s=kmax_r
       case(2)
          imin_s=imin_r;      imax_s=imax_r; 
          kmin_s=kmin_r+nK/2; kmax_s=kmax_r+nK/2;
       case(4)
          imin_s=imin_r+nI/2; imax_s=imax_r+nI/2; 
          kmin_s=kmin_r+nK/2; kmax_s=kmax_r+nK/2;
       end select
    case(5,6)
       if(oksend)then
          kmin_s=kmin_o; kmax_s=kmax_o
       else
          kmin_s=kmin_g; kmax_s=kmax_g
       end if
       select case(isubface)
          ! Beware, case(2) and case(3) are not swapped
       case(1)
          imin_s=imin_r;      imax_s=imax_r; 
          jmin_s=jmin_r;      jmax_s=jmax_r
       case(2)
          imin_s=imin_r+nI/2; imax_s=imax_r+nI/2; 
          jmin_s=jmin_r;      jmax_s=jmax_r
       case(3)
          imin_s=imin_r;      imax_s=imax_r; 
          jmin_s=jmin_r+nJ/2; jmax_s=jmax_r+nJ/2;
       case(4)
          imin_s=imin_r+nI/2; imax_s=imax_r+nI/2; 
          jmin_s=jmin_r+nJ/2; jmax_s=jmax_r+nJ/2;
       end select
    end select

  end subroutine setsubrange_ray

  !===========================================================================
  subroutine buf2rayface(&
       buf,imin,imax,jmin,jmax,kmin,kmax)

    integer, intent(in) :: imin,imax,jmin,jmax,kmin,kmax
    real, dimension(3,2,imin:imax,jmin:jmax,kmin:kmax),intent(inout) :: buf
    !-------------------------------------------------------------------------

    ! Take maximum of rayface and buf (more positive values are more real)
    ! for the full face

    rayface(:,:,imin_g:imax_g,jmin_g:jmax_g,kmin_g:kmax_g,iBLK)=max(buf,&
         rayface(:,:,imin_g:imax_g,jmin_g:jmax_g,kmin_g:kmax_g,iBLK))

  end subroutine buf2rayface

  !===========================================================================
  subroutine buf2sparserayface(buf,imin,imax,jmin,jmax,kmin,kmax)

    integer, intent(in) :: imin,imax,jmin,jmax,kmin,kmax
    real, dimension(3,2,imin:imax,jmin:jmax,kmin:kmax),intent(inout) :: buf
    !-------------------------------------------------------------------------

    ! Take maximum of rayface and buf (more positive values are more real)
    ! for a factor of 2 coarser grid

    rayface(:,:,imin_g:imax_g:2,jmin_g:jmax_g:2,kmin_g:kmax_g:2,iBLK)=max(buf,&
         rayface(:,:,imin_g:imax_g:2,jmin_g:jmax_g:2,kmin_g:kmax_g:2,iBLK))

  end subroutine buf2sparserayface

  !===========================================================================
  subroutine buf2subrayface(buf,imin,imax,jmin,jmax,kmin,kmax)

    integer, intent(in) :: imin,imax,jmin,jmax,kmin,kmax
    real, dimension(3,2,imin:imax,jmin:jmax,kmin:kmax),intent(inout) :: buf
    !-------------------------------------------------------------------------

    ! Set subface range to write into

    call setsubrange_ray(.false.)

    ! Take maximum of rayface and buf (more positive values are more real)
    ! for the appropriate subface

    rayface(:,:,imin_s:imax_s,jmin_s:jmax_s,kmin_s:kmax_s,iBLK)=max(buf,&
         rayface(:,:,imin_s:imax_s,jmin_s:jmax_s,kmin_s:kmax_s,iBLK))

  end subroutine buf2subrayface

  !===========================================================================
  subroutine prolong_ray

    ! For faces that are shared with a coarser neighbor, interpolate 
    ! for all points which are not coinciding and where the ray is going out.
    !
    ! a at odd  j and even k requires interpolation in direction k
    ! b at even j and odd  k requires interpolation in direction j
    ! c at even j and even k requires interpolation in both directions

    !(1,5)           (5,5)
    !   O-- b --O-- b --O
    !   |   |   |   |   |
    !   |   |   |   |   |
    !   a - c - a - c - a
    !   |   |   |   |   |
    !   |   |   |   |   |
    !   O-- b --O-- b --O
    !   |   |   |   |   |
    !   |   |   |   |   |
    !   a - c - a - c - a
    !   |   |   |   |   |
    !   |   |   |   |   |
    !   O-- b --O-- b --O
    !(1,1)           (5,1)

    integer :: iray
    integer :: j, k, nFaceJ, nFaceK
    integer, parameter :: nFaceMax=max(nI+1,nJ+1,nK+1)
    real    :: qrayface(3,2,nFaceMax,nFaceMax)
    integer :: qrayend_ind(2,nFaceMax,nFaceMax)

    ! Interpolation weights
    real, dimension(4), parameter:: weight4=0.25
    real, dimension(2), parameter:: weight2=0.5

    !-------------------------------------------------------------------------

    if(oktest)write(*,*)'Prolong_ray, me, iBLK, iface=',iProc, iBLK, iface

    ! Extract qrayface and qrayend_ind for the appropriate face
    !NOTE: qrayend_ind assignment split to two lines to avoid reshaping compiler bug!
    select case(iface)
    case(1)
       nFaceJ=nJ+1; nFaceK=nK+1
       qrayface( :,:,1:nFaceJ,1:nFaceK)=rayface(   :,:,1,1:nJ+1,1:nK+1,iBLK)
       qrayend_ind(1,1:nFaceJ,1:nFaceK)=rayend_ind(1,1,1,1:nJ+1,1:nK+1,iBLK)
       qrayend_ind(2,1:nFaceJ,1:nFaceK)=rayend_ind(1,2,1,1:nJ+1,1:nK+1,iBLK)
    case(2)
       nFaceJ=nJ+1; nFaceK=nK+1
       qrayface( :,:,1:nFaceJ,1:nFaceK)=rayface(   :,:,nI+1,1:nJ+1,1:nK+1,iBLK)
       qrayend_ind(1,1:nFaceJ,1:nFaceK)=rayend_ind(1,1,nI+1,1:nJ+1,1:nK+1,iBLK)
       qrayend_ind(2,1:nFaceJ,1:nFaceK)=rayend_ind(1,2,nI+1,1:nJ+1,1:nK+1,iBLK)
    case(3)
       nFaceJ=nI+1; nFaceK=nK+1
       qrayface( :,:,1:nFaceJ,1:nFaceK)=rayface(   :,:,1:nI+1,1,1:nK+1,iBLK)
       qrayend_ind(1,1:nFaceJ,1:nFaceK)=rayend_ind(1,1,1:nI+1,1,1:nK+1,iBLK)
       qrayend_ind(2,1:nFaceJ,1:nFaceK)=rayend_ind(1,2,1:nI+1,1,1:nK+1,iBLK)
    case(4)
       nFaceJ=nI+1; nFaceK=nK+1
       qrayface( :,:,1:nFaceJ,1:nFaceK)=rayface(   :,:,1:nI+1,nJ+1,1:nK+1,iBLK)
       qrayend_ind(1,1:nFaceJ,1:nFaceK)=rayend_ind(1,1,1:nI+1,nJ+1,1:nK+1,iBLK)
       qrayend_ind(2,1:nFaceJ,1:nFaceK)=rayend_ind(1,2,1:nI+1,nJ+1,1:nK+1,iBLK)
    case(5)
       nFaceJ=nI+1; nFaceK=nJ+1
       qrayface( :,:,1:nFaceJ,1:nFaceK)=rayface(   :,:,1:nI+1,1:nJ+1,1,iBLK)
       qrayend_ind(1,1:nFaceJ,1:nFaceK)=rayend_ind(1,1,1:nI+1,1:nJ+1,1,iBLK)
       qrayend_ind(2,1:nFaceJ,1:nFaceK)=rayend_ind(1,2,1:nI+1,1:nJ+1,1,iBLK)
    case(6)
       nFaceJ=nI+1; nFaceK=nJ+1
       qrayface( :,:,1:nFaceJ,1:nFaceK)=rayface(   :,:,1:nI+1,1:nJ+1,nK+1,iBLK)
       qrayend_ind(1,1:nFaceJ,1:nFaceK)=rayend_ind(1,1,1:nI+1,1:nJ+1,nK+1,iBLK)
       qrayend_ind(2,1:nFaceJ,1:nFaceK)=rayend_ind(1,2,1:nI+1,1:nJ+1,nK+1,iBLK)
    case default
       call stop_mpi('Impossible value for iface in prolong_ray')
    end select

    do iray=1,2
       do k=1,nfaceK
          if(mod(k,2)==1)then
             do j=2,nfaceJ,2
                ! Case b: even j and odd k

                if(qrayend_ind(iray,j,k)/=ray_out_)CYCLE

                call rayface_interpolate(&
                     qrayface(:,iray,j-1:j+1:2,k),weight2,2,&
                     qrayface(:,iray,j,k))
             end do
          else
             do j=1,nJ+1,2
                ! Case a: odd j and even k

                if(qrayend_ind(iray,j,k)/=ray_out_)CYCLE

                call rayface_interpolate(&
                     qrayface(:,iray,j,k-1:k+1:2),weight2,2,&
                     qrayface(:,iray,j,k))
             end do
             do j=2,nJ,2
                ! Case c: even j and even k

                if(qrayend_ind(iray,j,k)/=ray_out_)CYCLE

                call rayface_interpolate(&
                     qrayface(:,iray,j-1:j+1:2,k-1:k+1:2),weight4,4,&
                     qrayface(:,iray,j,k))
             end do ! j
          end if ! mod(k,2)
       end do ! k
    end do ! iray

    ! Put back result into rayface
    select case(iface)
    case(1)
       rayface(:,:,     1,1:nJ+1,1:nK+1,iBLK)=qrayface(:,:,1:nFaceJ,1:nFaceK)
    case(2)
       rayface(:,:,  nI+1,1:nJ+1,1:nK+1,iBLK)=qrayface(:,:,1:nFaceJ,1:nFaceK)
    case(3)
       rayface(:,:,1:nI+1,     1,1:nK+1,iBLK)=qrayface(:,:,1:nFaceJ,1:nFaceK)
    case(4)
       rayface(:,:,1:nI+1,  nJ+1,1:nK+1,iBLK)=qrayface(:,:,1:nFaceJ,1:nFaceK)
    case(5)
       rayface(:,:,1:nI+1,1:nJ+1,     1,iBLK)=qrayface(:,:,1:nFaceJ,1:nFaceK)
    case(6)
       rayface(:,:,1:nI+1,1:nJ+1,  nK+1,iBLK)=qrayface(:,:,1:nFaceJ,1:nFaceK)
    end select

  end subroutine prolong_ray
  !===========================================================================

end subroutine ray_pass_old
