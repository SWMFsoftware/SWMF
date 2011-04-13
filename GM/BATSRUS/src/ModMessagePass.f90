!^CFG COPYRIGHT UM

module ModMessagePass

  implicit none

contains

  subroutine message_pass_dir(&
       idirmin,idirmax,width,sendcorners,prolongorder,nVar,sol_BLK,&
       Sol_VGB,restrictface,DoTwoCoarseLayers)

    ! Exchange messages for array sol_BLK, ..., sol9_BLK or for Sol_VGB
    ! for directions idir=idirmin..idirmax using a ghost cell layer of 
    ! thickness width. 
    ! idir=1,2,3 correspond to east-west, south-north, bot-top.
    ! The actual number of variables to pass is nVar with 1<=nVar<=9.
    ! If sendcorners=.true., we send ghostcells already updated,
    ! in effect, the edge and corner information gets propagated
    ! For prolongorder=1 the prolongation operator for coarse to fine messages
    ! is only first order accurate. For prolongorder=2 there are various 2nd
    ! order prolongation operators implemented depending on the global variable
    !
    !   prolong_type='central'  (gradient based on central differencing)
    !                'central2' (fully 2nd order using one sided if needed)
    !                'lr'       (left and right gradients are calculated)
    !                'lr2'      (fully 2nd order using opposite side if needed) 
    !                'lr3'      (2nd order assuming that equal & restr. are done
    !                            for all dir., so sendcorners must be false.)
    !                'minmod'   (minmod limited gradient is calculated)
    !                'mindod2'  (fully 2nd order using unlimited if needed)
    !
    ! If restrictface=.true., average only the 4 fine cells next to the face
    ! If DoTwoCoarseLayers=.true., prolong the 2nd coarse cell into the 
    !                              2nd layer of fine cells for TVD res change.
    !                              Only for Sol_VGB and prolongorder=1.

    ! Notation: _o original   (for equal blocks) 
    !           _g ghost      (for equal blocks)
    !           R  originals to be restricted
    !           _r restricted (to be sent to a coarser block)
    !           P  originals to be prolonged
    !           _p prolonged  (to be sent to a refined block)
    !           _s subface    (one quarter of a face)

    use ModProcMH
    use ModMain, ONLY: nI,nJ,nK,nBLK,prolong_type,nBlockMax,okdebug,unusedBLK,&
         east_,west_,south_,north_,bot_,top_, optimize_message_pass, BlkTest,&
         UseBatl
    use BATL_lib, ONLY: message_pass_cell

    use ModAMR, ONLY : unusedBlock_BP,child2subface
    use ModParallel, ONLY : NOBLK, neiLEV,neiPE,neiBLK, BLKneighborCHILD
    use ModVarIndexes,ONLY:MaxVarState=>nVar
    use ModMpi
    implicit none

    ! Maximum number of variables to be sent and maximum width
    integer, parameter :: maxvar=max(9,MaxVarState), maxwidth=2

    ! Parameters

    integer, intent(in) :: idirmin,idirmax,width
    logical, intent(in) :: sendcorners
    integer, intent(in) :: prolongorder,nvar

    ! Optional parameters

    real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), optional, intent(inout) :: &
         sol_BLK

    real, dimension(nVar,-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), &
         optional, intent(inout) :: Sol_VGB

    logical, optional, intent(in) :: restrictface, DoTwoCoarseLayers

    ! Local variables

    integer :: iError

    logical :: qrestrictface

    integer :: isweep, idir, iface, otherface, iside, isubface
    integer :: ibuf, width1, width2, iBLK
    integer :: imin_g,imax_g,jmin_g,jmax_g,kmin_g,kmax_g
    integer :: imin_o,imax_o,jmin_o,jmax_o,kmin_o,kmax_o
    integer :: imin_r,imax_r,jmin_r,jmax_r,kmin_r,kmax_r
    integer :: iminR,imaxR,jminR,jmaxR,kminR,kmaxR
    integer :: imin_p,imax_p,jmin_p,jmax_p,kmin_p,kmax_p
    integer :: iminP,imaxP,jminP,jmaxP,kminP,kmaxP
    integer :: imin_s,imax_s,jmin_s,jmax_s,kmin_s,kmax_s
    integer       :: neiP,neiB,neiL,isize,isize_r,ichild,itag,request
    integer       :: number_receive_requests=0
    integer, save :: receive_requests(nBLK*6)
    integer       :: status(MPI_STATUS_SIZE, nBLK*6)  

    ! Maximum size of the ghost cell layer to be sent. 
    integer, parameter :: maxsize=maxwidth*maxvar* &
         max((nI+4)*(nJ+4),(nI+4)*nK,(nJ+4)*nK)

    ! Buffer to hold incoming layers of ghost cell values.
    ! For all blocks and all faces
    real, dimension(:,:), allocatable :: buffer

    ! Equal, restricted and prolonged values are stored in these arrays
    real, dimension(:,:,:,:), allocatable :: eq_buf, re_buf, pr_buf, &
         pr_buf_s, avrg_state

    ! Small allocatable arrays for prolongation
    real, dimension(:,:,:), allocatable :: qsol,avrg,avrg1,gradx,grady,gradz,&
         gradxl,gradxr,gradyl,gradyr,gradzl,gradzr,pr_buf1

    ! Number of coarse layers to be prolonged: 
    ! 1 for 1st order, 2 for TVD res change.
    ! nCoarse1 = nCoarseLayer-1
    integer :: nCoarseLayer, nCoarse1

    logical :: oktest, oktest_me

    !--------------------------------------------------------------------------

    if(.not.allocated(buffer)) allocate(buffer(nBLK*maxsize,2))

    nCoarseLayer=1
    if(present(DoTwoCoarseLayers))then
       if(DoTwoCoarseLayers) nCoarseLayer=2
    end if

    if(UseBatl)then
       if(.not.present(Sol_VGB)) &
            call stop_mpi('message_pass_dir without Sol_VGB, and UseBatl')

       call message_pass_cell(nVar,Sol_VGB,width,prolongorder,nCoarseLayer,&
            sendcorners, restrictface)

       RETURN
    end if


    call set_oktest('message_pass_dir',oktest, oktest_me)

    if(present(restrictface))then
       qrestrictface=restrictface
    else
       qrestrictface=.false.
    endif

    nCoarse1 = nCoarseLayer - 1

    if(oktest)write(*,*)&
         'message_pass_dir me,idirmin,max,width,sendcorners,prolorder,nvar,',&
         'Sol_BLK,Sol_VGB,restrictface,nCoarseLayer=',&
         iProc,idirmin,idirmax,width,sendcorners,prolongorder,nvar,&
         present(Sol_BLK),present(Sol_VGB),qrestrictface,nCoarseLayer

    if(present(Sol_VGB).eqv.present(Sol_BLK))&
         call stop_mpi('message_pass_dir must be called with '// &
         'exactly one of the Sol_VGB and Sol_BLK parameters')

    if(oktest_me)write(*,*)'prolong_type=',prolong_type,&
         ' optimize_message_pass=',optimize_message_pass

    if(nvar>maxvar.or.nvar<1)&
         call stop_mpi('message_pass_dir can only be used with 1 <= nvar <= 9')

    if(idirmin<1 .or. idirmax>3 .or. idirmin>idirmax)call stop_mpi(&
         'Incorrect value for 1<=idirmin<=idirmax<=3 in message_pass_dir!')

    if(width/=1.and.width/=2)call stop_mpi(&
         'width must be 1 or 2 in message_pass_dir!')

    if(sendcorners.and.width==1)call stop_mpi(&
         'width must be 2 if sendcorners=T in message_pass_dir!')

    if(nCoarseLayer==2) then
       if(.not.present(Sol_VGB)) call stop_mpi('ERROR in message_pass_dir: '// &
            'DoTwoCoarseLayer=2 requires Sol_VGB to be present')
       if(prolongorder/=1) call stop_mpi('ERROR in message_pass_dir: '// &
            'DoTwoCoarseLayer=2 requires prolongorder=1')
    end if

    ! These will be useful for index limits
    width1=width-1; width2=2*width-1


    if(prolongorder==1)then
       select case(optimize_message_pass)
       case('face')
          ! Send messages face by face
          do isweep=2*idirmin-1,2*idirmax
             call msg_pass_faces(isweep,isweep,.true.,.true.,.true.)
          end do
       case('min')
          ! Send messages face by face and kind by kind
          do isweep=2*idirmin-1,2*idirmax
             ! Send equal
             call msg_pass_faces(isweep,isweep,.true.,.false.,.false.)
             ! Send restricted
             call msg_pass_faces(isweep,isweep,.false.,.true.,.false.)
             ! Send prolonged
             call msg_pass_faces(isweep,isweep,.false.,.false.,.true.)
          end do
       case default
          ! Send messages direction by direction
          do isweep=idirmin,idirmax
             call msg_pass_faces(2*isweep-1,2*isweep,.true.,.true.,.true.)
          end do
       end select
    else
       ! Prolongorder is 2, prolonged must be done after equal and coarse
       select case(optimize_message_pass)
       case('face')
          if(sendcorners)then
             ! Send messages face by face, but do prolonged messages after
             ! equal and restricted are done for both faces in that direction
             do isweep=idirmin,idirmax
                ! Do equal and restricted for the 2 faces
                call msg_pass_faces(2*isweep-1,2*isweep-1,.true.,.true.,.false.)
                call msg_pass_faces(2*isweep  ,2*isweep,  .true.,.true.,.false.)
                ! Do prolonged for the 2 faces
                call msg_pass_faces(2*isweep-1,2*isweep-1,.false.,.false.,.true.)
                call msg_pass_faces(2*isweep,  2*isweep,  .false.,.false.,.true.)
             end do
          else
             ! Do all restricted and equal face by face
             do isweep=2*idirmin-1,2*idirmax
                call msg_pass_faces(isweep,isweep,.true. ,.true. ,.false.)
             end do
             ! Do prolonged for all the faces
             do isweep=2*idirmin-1,2*idirmax
                call msg_pass_faces(isweep,isweep,.false.,.false.,.true. )
             end do
          end if
       case('min')
          if(sendcorners)then
             ! Send messages face by face and kind by kind, but do prolonged 
             ! messages later when equal and restricted are done for both faces
             ! in that direction
             do isweep=idirmin,idirmax
                ! Do equal then restricted for the two faces
                call msg_pass_faces(2*isweep-1,2*isweep-1,.true.,.false.,.false.)
                call msg_pass_faces(2*isweep-1,2*isweep-1,.false.,.true.,.false.)
                call msg_pass_faces(2*isweep,  2*isweep,  .true.,.false.,.false.)
                call msg_pass_faces(2*isweep,  2*isweep,  .false.,.true.,.false.)
                ! Do prolonged for the two faces
                call msg_pass_faces(2*isweep-1,2*isweep-1,.false.,.false.,.true.)
                call msg_pass_faces(2*isweep,  2*isweep,  .false.,.false.,.true.)
             end do
          else
             ! Send equal and restricted messages face by face and kind by kind
             do isweep=2*idirmin-1,2*idirmax
                call msg_pass_faces(isweep,isweep,.true. ,.false.,.false.)
                call msg_pass_faces(isweep,isweep,.false.,.true. ,.false.)
             end do
             ! Do prolonged for all the faces
             do isweep=2*idirmin-1,2*idirmax
                call msg_pass_faces(isweep,isweep,.false.,.false.,.true. )
             end do
          end if
       case default
          if(sendcorners)then
             !Send messages direction by direction but do prolonged later
             do isweep=idirmin,idirmax
                ! Do equal and coarse
                call msg_pass_faces(2*isweep-1,2*isweep,.true. ,.true. ,.false.)
                ! Do prolonged
                call msg_pass_faces(2*isweep-1,2*isweep,.false.,.false.,.true. )
             end do
          else
             !Send messages for all directions but do prolonged later
             do isweep=idirmin,idirmax
                ! Do equal and restricted
                call msg_pass_faces(2*isweep-1,2*isweep,.true. ,.true. ,.false.)
             end do
             do isweep=idirmin,idirmax
                ! Do prolonged
                call msg_pass_faces(2*isweep-1,2*isweep,.false.,.false.,.true. )
             end do
          end if
       end select
    end if ! prolongorder

  contains

    subroutine msg_pass_faces(&
         ifacemin,ifacemax,do_equal,do_restricted,do_prolonged)

      integer:: ifacemin,ifacemax
      logical:: do_equal,do_restricted,do_prolonged
      !------------------------------------------------------------------------

      ! Debug
      if(okdebug)buffer=0.00

      number_receive_requests = 0
      receive_requests = MPI_REQUEST_NULL

      do iface=ifacemin,ifacemax

         ! Set index ranges for the face
         call setranges

         ! Size of ghost cell layer including all the variables
         isize=(imax_g-imin_g+1)*(jmax_g-jmin_g+1)*(kmax_g-kmin_g+1)*nvar

         ! Size of restricted ghost cell layer
         isize_r=isize/4

         do iBLK = 1,nBlockMax
            if(unusedBLK(iBLK)) CYCLE
            ! Post non-blocking receive for opposite face of neighbor block
            select case(neiLEV(otherface,iBLK))
            case(0)
               if(do_equal)call remote_receive(1,isize)
            case(-1)
               if(do_restricted)call remote_receive(4,isize_r)
            case(1)
               if(do_prolonged)call remote_receive(1,isize)
            case(NOBLK)
               ! Do nothing
            case default
               write(*,*)'me,iBLK,otherface,neiLEV=',&
                    iProc,iBLK,otherface,neiLEV(otherface,iBLK)
               call stop_mpi(&
                    'Error in message pass: Invalid value for neiLEV')
            end select
         end do ! iBLK
      end do ! iface

      !\
      ! Wait for all receive commands to be posted for all processors
      !/
      call barrier_mpi

      !\
      ! Send blocking messages with Rsend (ready to receive)
      !/
      do iface=ifacemin,ifacemax

         ! Set index ranges for the face
         call setranges

         if(present(Sol_VGB))then
            if(do_equal)&
                 allocate(eq_buf(nVar,imin_o:imax_o,jmin_o:jmax_o,kmin_o:kmax_o))
            if(do_restricted)&
                 allocate(re_buf(nVar,imin_r:imax_r,jmin_r:jmax_r,kmin_r:kmax_r))
            if(do_prolonged)&
                 allocate(&
                 pr_buf(nVar,imin_p:imax_p,jmin_p:jmax_p,kmin_p:kmax_p),&
                 pr_buf_s(nVar,imin_o:imax_o,jmin_o:jmax_o,kmin_o:kmax_o),&
                 avrg_state(nVar,iminP:imaxP,jminP:jmaxP,kminP:kmaxP))
         else
            if(do_equal.or.do_prolonged)&
                 allocate(eq_buf(imin_o:imax_o,jmin_o:jmax_o,kmin_o:kmax_o,nVar))
            if(do_restricted)&
                 allocate(re_buf(imin_r:imax_r,jmin_r:jmax_r,kmin_r:kmax_r,nVar))
            if(do_prolonged)&
                 allocate(&
                 pr_buf(imin_p:imax_p,jmin_p:jmax_p,kmin_p:kmax_p,nVar),&
                 pr_buf_s(imin_o:imax_o,jmin_o:jmax_o,kmin_o:kmax_o,nVar),&
                 avrg1(iminP:imaxP,jminP:jmaxP,kminP:kmaxP))
         end if

         if(do_prolonged.and.prolongorder>1)then
            allocate(&
                 qsol(iminP-1:imaxP+1,jminP-1:jmaxP+1,kminP-1:kmaxP+1),&
                 avrg(iminP:imaxP,jminP:jmaxP,kminP:kmaxP),&
                 pr_buf1(imin_p:imax_p,jmin_p:jmax_p,kmin_p:kmax_p))
            if(prolong_type(1:7)=='central'.or.prolong_type(1:6)=='minmod')&
                 allocate(&
                 gradx(iminP:imaxP,jminP:jmaxP,kminP:kmaxP),&
                 grady(iminP:imaxP,jminP:jmaxP,kminP:kmaxP),&
                 gradz(iminP:imaxP,jminP:jmaxP,kminP:kmaxP))
            if(prolong_type(1:2)=='lr'.or.prolong_type(1:6)=='minmod')&
                 allocate(&
                 gradxl(iminP:imaxP,jminP:jmaxP,kminP:kmaxP),&
                 gradyl(iminP:imaxP,jminP:jmaxP,kminP:kmaxP),&
                 gradzl(iminP:imaxP,jminP:jmaxP,kminP:kmaxP),&
                 gradxr(iminP:imaxP,jminP:jmaxP,kminP:kmaxP),&
                 gradyr(iminP:imaxP,jminP:jmaxP,kminP:kmaxP),&
                 gradzr(iminP:imaxP,jminP:jmaxP,kminP:kmaxP))
         end if

         do iBLK=1,nBlockMax
            if(unusedBLK(iBLK))CYCLE
            neiL=neiLEV(iface,iBLK)
            select case(neiL)
            case(0)
               if(do_equal)call send_equal
            case(1)
               if(do_restricted)call send_restricted
            case(-1)
               if(do_prolonged)call send_prolonged
            case(NOBLK)
               ! There is no neighbor, do nothing
            case default
               write(*,*)'me,iBLK,iface,neiLEV=',&
                    iProc,iBLK,iface,neiLEV(iface,iBLK)
               call stop_mpi(&
                    'Error in message pass: Invalid value for neiLEV')
            end select
         end do ! iBLK

         if(do_equal)deallocate(eq_buf)
         if(do_restricted)deallocate(re_buf)
         if(do_prolonged)then
            deallocate(pr_buf,pr_buf_s)
            if(present(Sol_VGB))then
               deallocate(avrg_state)
            else
               deallocate(avrg1)
            end if
            if(prolongorder>1)then
               deallocate(qsol,avrg,pr_buf1)
               if(prolong_type(1:7)=='central'.or.prolong_type(1:6)=='minmod')&
                    deallocate(gradx,grady,gradz)
               if(prolong_type(1:2)=='lr'.or.prolong_type(1:6)=='minmod')&
                    deallocate(gradxl,gradyl,gradzl,gradxr,gradyr,gradzr)
            end if
         end if

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
         call setranges

         do iBLK = 1,nBlockMax
            if(unusedBLK(iBLK))CYCLE
            ibuf=(iBLK-1)*isize+1
            select case(neiLEV(otherface,iBLK))
            case(0)
               neiP=neiPE(1,otherface,iBLK)
               if(do_equal .and. neiP/=iProc .and. &
                    .not.unusedBlock_BP(neiBLK(1,otherface,iBLK),neiP))then
                  if(present(Sol_VGB))then
                     call buffer2face_state(buffer(ibuf,iside),&
                          imin_g,imax_g,jmin_g,jmax_g,kmin_g,kmax_g,nvar)
                  else
                     call buffer2face(buffer(ibuf,iside),&
                          imin_g,imax_g,jmin_g,jmax_g,kmin_g,kmax_g,nvar)
                  end if
               end if
            case(-1)
               if(do_restricted)then
                  if(present(Sol_VGB))then
                     call buffer2subfaces_state(buffer(ibuf,iside),&
                          imin_r,imax_r,jmin_r,jmax_r,kmin_r,kmax_r,nvar)
                  else
                     call buffer2subfaces(buffer(ibuf,iside),&
                          imin_r,imax_r,jmin_r,jmax_r,kmin_r,kmax_r,nvar)
                  end if
               end if
            case(1)
               neiP=neiPE(1,otherface,iBLK)
               if(do_prolonged .and. neiP/=iProc .and. &
                    .not.unusedBlock_BP(neiBLK(1,otherface,iBLK),neiP))then
                  if(present(Sol_VGB))then
                     call buffer2face_state(buffer(ibuf,iside),&
                          imin_g,imax_g,jmin_g,jmax_g,kmin_g,kmax_g,nvar)
                  else
                     call buffer2face(buffer(ibuf,iside),&
                          imin_g,imax_g,jmin_g,jmax_g,kmin_g,kmax_g,nvar)
                  end if
               end if
            end select
         end do ! iBLK
      end do ! iface

      if(oktest)write(*,*)'msg_pass_faces finished: me, ifacemin, ifacemax=',&
           iProc,ifacemin,ifacemax

    end subroutine msg_pass_faces

    !==========================================================================

    subroutine setranges

      ! Set ranges orthogonal to idir based on the value of idir
      ! Include ghostcells for directions<idir if sendcorners is true.

      idir=(iface+1)/2

      if(idir>1.and.sendcorners)then
         imin_g=-1; imax_g=nI+2  ; imin_o=-1; imax_o=nI+2;
         imin_r= 0; imax_r=nI/2+1; iminR =-1; imaxR =nI+2;
         imin_p=-1; imax_p=nI*2+2; iminP = 0; imaxP =nI+1;
      else if(idir/=1)then
         imin_g=1;  imax_g=nI;     imin_o=1;  imax_o=nI
         imin_r=1;  imax_r=nI/2;   iminR =1;  imaxR =nI;
         imin_p=1;  imax_p=nI*2;   iminP =1;  imaxP =nI;
      end if

      if(idir>2.and.sendcorners)then
         jmin_g=-1; jmax_g=nJ+2  ; jmin_o=-1; jmax_o=nJ+2;
         jmin_r= 0; jmax_r=nJ/2+1; jminR =-1; jmaxR =nJ+2;
         jmin_p=-1; jmax_p=nJ*2+2; jminP = 0; jmaxP =nJ+1;
      else if(idir/=2)then
         jmin_g=1;  jmax_g=nJ;     jmin_o=1;  jmax_o=nJ;
         jmin_r=1;  jmax_r=nJ/2;   jminR =1;  jmaxR =nJ;
         jmin_p=1;  jmax_p=nJ*2;   jminP =1;  jmaxP =nJ;
      endif

      if(idir/=3)then
         kmin_g=1;  kmax_g=nK;     kmin_o=1;  kmax_o=nK;
         kmin_r=1;  kmax_r=nK/2;   kminR =1;  kmaxR =nK;
         kmin_p=1;  kmax_p=nK*2;   kminP =1;  kmaxP =nK;
      end if

      ! Set ranges in direction of idir based on the value of iface

      select case(iface)
      case(1)
         otherface=2
         imin_o=1;    imax_o=width;
         imin_g=nI+1; imax_g=nI+width;
         iminR =1;    imaxR=2*width;
         imin_r=1;    imax_r=width;
         iminP =1;    imaxP=nCoarseLayer;
         imin_p=1;    imax_p=2;
      case(3)
         otherface=4
         jmin_o=1;    jmax_o=width;
         jmin_g=nJ+1; jmax_g=nJ+width; 
         jminR =1;    jmaxR=2*width;
         jmin_r=1;    jmax_r=width; 
         jminP =1;    jmaxP=nCoarseLayer;
         jmin_p=1;    jmax_p=2;
      case(5)
         otherface=6
         kmin_o=1;    kmax_o=width;  
         kmin_g=nK+1; kmax_g=nK+width;
         kminR =1;    kmaxR =2*width; 
         kmin_r=1;    kmax_r=width;
         kminP =1;    kmaxP=nCoarseLayer;
         kmin_p=1;    kmax_p=2;
      case(2)
         otherface=1
         imin_o=nI-width1;   imax_o=nI;
         imin_g=-width1;     imax_g=0; 
         iminR=nI-width2;    imaxR =nI;
         imin_r=nI/2-width1; imax_r=nI/2;
         iminP =nI-nCoarse1; imaxP =nI;
         imin_p=nI*2-1;      imax_p=nI*2;
      case(4)
         otherface=3
         jmin_o=nJ-width1;   jmax_o=nJ;
         jmin_g=-width1;     jmax_g=0; 
         jminR =nJ-width2;   jmaxR =nJ;
         jmin_r=nJ/2-width1; jmax_r=nJ/2;
         jminP =nJ-nCoarse1; jmaxP =nJ;
         jmin_p=nJ*2-1;      jmax_p=nJ*2;
      case(6)
         otherface=5
         kmin_o=nK-width1;   kmax_o=nK;
         kmin_g=-width1;     kmax_g=0; 
         kminR=nK-width2;    kmaxR=nK;
         kmin_r=nK/2-width1; kmax_r=nK/2;
         kminP=nK-nCoarse1;  kmaxP=nK;
         kmin_p=nK*2-1;      kmax_p=nK*2;
      end select

      iside=(iface-otherface+3)/2

    end subroutine setranges

    !==========================================================================

    subroutine setsubrange_g

      ! Select appropriate quarter of ghost cell layer

      select case(iface)
      case(1,2)
         imin_s=imin_g; imax_s=imax_g
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
         jmin_s=jmin_g; jmax_s=jmax_g
         select case(isubface)
            ! Beware, case(2) and case(3) are swapped
         case(1)
            imin_s=imin_r; imax_s=imax_r; 
            kmin_s=kmin_r; kmax_s=kmax_r
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
         kmin_s=kmin_g; kmax_s=kmax_g
         select case(isubface)
            ! Beware, case(2) and case(3) are not swapped
         case(1)
            imin_s=imin_r; imax_s=imax_r; 
            jmin_s=jmin_r; jmax_s=jmax_r
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

    end subroutine setsubrange_g

    !==========================================================================

    subroutine setsubrange_p

      ! Select appropriate quarter of prolonged originals
      ! Select the appropriate layer in the orthogonal direction for width==1

      select case(iface)
      case(1,2)
         if(iface==1)then
            imin_s=1;           imax_s=width;
         else
            imin_s=nI*2-width1; imax_s=nI*2;
         end if
         select case(isubface)
            ! Beware, case(2) and case(3) are not swapped
         case(1)
            jmin_s=jmin_o; jmax_s=jmax_o; 
            kmin_s=kmin_o; kmax_s=kmax_o
         case(3)
            jmin_s=jmin_o+nJ; jmax_s=jmax_o+nJ; 
            kmin_s=kmin_o; kmax_s=kmax_o
         case(2)
            jmin_s=jmin_o; jmax_s=jmax_o; 
            kmin_s=kmin_o+nK; kmax_s=kmax_o+nK;
         case(4)
            jmin_s=jmin_o+nJ; jmax_s=jmax_o+nJ; 
            kmin_s=kmin_o+nK; kmax_s=kmax_o+nK;
         end select
      case(3,4)
         if(iface==3)then
            jmin_s=1;           jmax_s=width;
         else
            jmin_s=nJ*2-width1; jmax_s=nJ*2;
         end if
         select case(isubface)
            ! Beware, case(2) and case(3) are swapped
         case(1)
            imin_s=imin_o; imax_s=imax_o; 
            kmin_s=kmin_o; kmax_s=kmax_o
         case(3)
            imin_s=imin_o+nI; imax_s=imax_o+nI; 
            kmin_s=kmin_o;    kmax_s=kmax_o
         case(2)
            imin_s=imin_o;    imax_s=imax_o; 
            kmin_s=kmin_o+nK; kmax_s=kmax_o+nK;
         case(4)
            imin_s=imin_o+nI; imax_s=imax_o+nI; 
            kmin_s=kmin_o+nK; kmax_s=kmax_o+nK;
         end select
      case(5,6)
         if(iface==5)then
            kmin_s=1;           kmax_s=width;
         else
            kmin_s=nK*2-width1; kmax_s=nK*2;
         end if
         select case(isubface)
            ! Beware, case(2) and case(3) are not swapped
         case(1)
            imin_s=imin_o; imax_s=imax_o; 
            jmin_s=jmin_o; jmax_s=jmax_o
         case(2)
            imin_s=imin_o+nI; imax_s=imax_o+nI; 
            jmin_s=jmin_o;    jmax_s=jmax_o
         case(3)
            imin_s=imin_o;    imax_s=imax_o; 
            jmin_s=jmin_o+nJ; jmax_s=jmax_o+nJ;
         case(4)
            imin_s=imin_o+nI; imax_s=imax_o+nI; 
            jmin_s=jmin_o+nJ; jmax_s=jmax_o+nJ;
         end select
      end select

    end subroutine setsubrange_p

    !==========================================================================

    subroutine restrict(ivar,sol_BLK)

      ! Take average of 8 small cells in sol_BLK and put it into re_buf
      integer, intent(in) :: ivar
      real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(in) :: sol_BLK

      re_buf(:,:,:,ivar)=0.125*(&
           sol_BLK(iminR  :imaxR:2,jminR  :jmaxR:2,kminR  :kmaxR:2,iBLK)+&
           sol_BLK(iminR+1:imaxR:2,jminR  :jmaxR:2,kminR  :kmaxR:2,iBLK)+&
           sol_BLK(iminR  :imaxR:2,jminR+1:jmaxR:2,kminR  :kmaxR:2,iBLK)+&
           sol_BLK(iminR  :imaxR:2,jminR  :jmaxR:2,kminR+1:kmaxR:2,iBLK)+&
           sol_BLK(iminR+1:imaxR:2,jminR+1:jmaxR:2,kminR  :kmaxR:2,iBLK)+&
           sol_BLK(iminR+1:imaxR:2,jminR  :jmaxR:2,kminR+1:kmaxR:2,iBLK)+&
           sol_BLK(iminR  :imaxR:2,jminR+1:jmaxR:2,kminR+1:kmaxR:2,iBLK)+&
           sol_BLK(iminR+1:imaxR:2,jminR+1:jmaxR:2,kminR+1:kmaxR:2,iBLK))

      ! redo the cells next to resolution changes if required by restrictface
      if(qrestrictface)call restrict_face(ivar,sol_BLK)

    end subroutine restrict

    !==========================================================================

    subroutine restrict_face(ivar,sol_BLK)

      ! Take average of 4 fine cells in sol_BLK next to the resolution change
      ! and put the result into re_buf

      integer, intent(in) :: ivar
      real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(in) :: sol_BLK

      select case(iface)
      case(east_)
         re_buf(imin_r,:,:,ivar)=0.25*(&
              sol_BLK(iminR,jminR  :jmaxR:2,kminR  :kmaxR:2,iBLK)+&
              sol_BLK(iminR,jminR+1:jmaxR:2,kminR  :kmaxR:2,iBLK)+&
              sol_BLK(iminR,jminR  :jmaxR:2,kminR+1:kmaxR:2,iBLK)+&
              sol_BLK(iminR,jminR+1:jmaxR:2,kminR+1:kmaxR:2,iBLK))
      case(west_)
         re_buf(imax_r,:,:,ivar)=0.25*(&
              sol_BLK(imaxR,jminR  :jmaxR:2,kminR  :kmaxR:2,iBLK)+&
              sol_BLK(imaxR,jminR+1:jmaxR:2,kminR  :kmaxR:2,iBLK)+&
              sol_BLK(imaxR,jminR  :jmaxR:2,kminR+1:kmaxR:2,iBLK)+&
              sol_BLK(imaxR,jminR+1:jmaxR:2,kminR+1:kmaxR:2,iBLK))
      case(south_)
         re_buf(:,jmin_r,:,ivar)=0.25*(&
              sol_BLK(iminR  :imaxR:2,jminR,kminR  :kmaxR:2,iBLK)+&
              sol_BLK(iminR+1:imaxR:2,jminR,kminR  :kmaxR:2,iBLK)+&
              sol_BLK(iminR  :imaxR:2,jminR,kminR+1:kmaxR:2,iBLK)+&
              sol_BLK(iminR+1:imaxR:2,jminR,kminR+1:kmaxR:2,iBLK))
      case(north_)
         re_buf(:,jmax_r,:,ivar)=0.25*(&
              sol_BLK(iminR  :imaxR:2,jmaxR,kminR  :kmaxR:2,iBLK)+&
              sol_BLK(iminR+1:imaxR:2,jmaxR,kminR  :kmaxR:2,iBLK)+&
              sol_BLK(iminR  :imaxR:2,jmaxR,kminR+1:kmaxR:2,iBLK)+&
              sol_BLK(iminR+1:imaxR:2,jmaxR,kminR+1:kmaxR:2,iBLK))
      case(bot_)
         re_buf(:,:,kmin_r,ivar)=0.25*(&
              sol_BLK(iminR  :imaxR:2,jminR  :jmaxR:2,kminR,iBLK)+&
              sol_BLK(iminR+1:imaxR:2,jminR  :jmaxR:2,kminR,iBLK)+&
              sol_BLK(iminR  :imaxR:2,jminR+1:jmaxR:2,kminR,iBLK)+&
              sol_BLK(iminR+1:imaxR:2,jminR+1:jmaxR:2,kminR,iBLK))
      case(top_)
         re_buf(:,:,kmax_r,ivar)=0.25*(&
              sol_BLK(iminR  :imaxR:2,jminR  :jmaxR:2,kmaxR,iBLK)+&
              sol_BLK(iminR+1:imaxR:2,jminR  :jmaxR:2,kmaxR,iBLK)+&
              sol_BLK(iminR  :imaxR:2,jminR+1:jmaxR:2,kmaxR,iBLK)+&
              sol_BLK(iminR+1:imaxR:2,jminR+1:jmaxR:2,kmaxR,iBLK))
      end select
    end subroutine restrict_face

    !==========================================================================

    subroutine restrict_state

      ! Take average of 8 small cells in Sol_VGB and put it into re_buf

      re_buf(:,:,:,:)=0.125*(&
           Sol_VGB(:,iminR  :imaxR:2,jminR  :jmaxR:2,kminR  :kmaxR:2,iBLK)+&
           Sol_VGB(:,iminR+1:imaxR:2,jminR  :jmaxR:2,kminR  :kmaxR:2,iBLK)+&
           Sol_VGB(:,iminR  :imaxR:2,jminR+1:jmaxR:2,kminR  :kmaxR:2,iBLK)+&
           Sol_VGB(:,iminR  :imaxR:2,jminR  :jmaxR:2,kminR+1:kmaxR:2,iBLK)+&
           Sol_VGB(:,iminR+1:imaxR:2,jminR+1:jmaxR:2,kminR  :kmaxR:2,iBLK)+&
           Sol_VGB(:,iminR+1:imaxR:2,jminR  :jmaxR:2,kminR+1:kmaxR:2,iBLK)+&
           Sol_VGB(:,iminR  :imaxR:2,jminR+1:jmaxR:2,kminR+1:kmaxR:2,iBLK)+&
           Sol_VGB(:,iminR+1:imaxR:2,jminR+1:jmaxR:2,kminR+1:kmaxR:2,iBLK))

      ! redo the cells next to resolution changes if required by restrictface
      if(qrestrictface)call restrict_face_state

    end subroutine restrict_state

    !==========================================================================

    subroutine restrict_face_state

      ! Take average of 4 fine cells in sol_BLK next to the resolution change
      ! and put the result into re_buf

      select case(iface)
      case(east_)
         re_buf(:,imin_r,:,:)=0.25*(&
              Sol_VGB(:,iminR,jminR  :jmaxR:2,kminR  :kmaxR:2,iBLK)+&
              Sol_VGB(:,iminR,jminR+1:jmaxR:2,kminR  :kmaxR:2,iBLK)+&
              Sol_VGB(:,iminR,jminR  :jmaxR:2,kminR+1:kmaxR:2,iBLK)+&
              Sol_VGB(:,iminR,jminR+1:jmaxR:2,kminR+1:kmaxR:2,iBLK))
      case(west_)
         re_buf(:,imax_r,:,:)=0.25*(&
              Sol_VGB(:,imaxR,jminR  :jmaxR:2,kminR  :kmaxR:2,iBLK)+&
              Sol_VGB(:,imaxR,jminR+1:jmaxR:2,kminR  :kmaxR:2,iBLK)+&
              Sol_VGB(:,imaxR,jminR  :jmaxR:2,kminR+1:kmaxR:2,iBLK)+&
              Sol_VGB(:,imaxR,jminR+1:jmaxR:2,kminR+1:kmaxR:2,iBLK))
      case(south_)
         re_buf(:,:,jmin_r,:)=0.25*(&
              Sol_VGB(:,iminR  :imaxR:2,jminR,kminR  :kmaxR:2,iBLK)+&
              Sol_VGB(:,iminR+1:imaxR:2,jminR,kminR  :kmaxR:2,iBLK)+&
              Sol_VGB(:,iminR  :imaxR:2,jminR,kminR+1:kmaxR:2,iBLK)+&
              Sol_VGB(:,iminR+1:imaxR:2,jminR,kminR+1:kmaxR:2,iBLK))
      case(north_)
         re_buf(:,:,jmax_r,:)=0.25*(&
              Sol_VGB(:,iminR  :imaxR:2,jmaxR,kminR  :kmaxR:2,iBLK)+&
              Sol_VGB(:,iminR+1:imaxR:2,jmaxR,kminR  :kmaxR:2,iBLK)+&
              Sol_VGB(:,iminR  :imaxR:2,jmaxR,kminR+1:kmaxR:2,iBLK)+&
              Sol_VGB(:,iminR+1:imaxR:2,jmaxR,kminR+1:kmaxR:2,iBLK))
      case(bot_)
         re_buf(:,:,:,kmin_r)=0.25*(&
              Sol_VGB(:,iminR  :imaxR:2,jminR  :jmaxR:2,kminR,iBLK)+&
              Sol_VGB(:,iminR+1:imaxR:2,jminR  :jmaxR:2,kminR,iBLK)+&
              Sol_VGB(:,iminR  :imaxR:2,jminR+1:jmaxR:2,kminR,iBLK)+&
              Sol_VGB(:,iminR+1:imaxR:2,jminR+1:jmaxR:2,kminR,iBLK))
      case(top_)
         re_buf(:,:,:,kmax_r)=0.25*(&
              Sol_VGB(:,iminR  :imaxR:2,jminR  :jmaxR:2,kmaxR,iBLK)+&
              Sol_VGB(:,iminR+1:imaxR:2,jminR  :jmaxR:2,kmaxR,iBLK)+&
              Sol_VGB(:,iminR  :imaxR:2,jminR+1:jmaxR:2,kmaxR,iBLK)+&
              Sol_VGB(:,iminR+1:imaxR:2,jminR+1:jmaxR:2,kmaxR,iBLK))
      end select
    end subroutine restrict_face_state

    !==========================================================================

    subroutine prolong1(ivar,sol_BLK)

      ! First order prolongation of sol_BLK(imin,...,iBLK) into pr_buf

      integer, intent(in) :: ivar
      real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), intent(in) :: sol_BLK

      !------------------------------------------------------------------------
      avrg1=sol_BLK(iminP:imaxP,jminP:jmaxP,kminP:kmaxP,iBLK)

      pr_buf(imin_p  :imax_p:2,jmin_p  :jmax_p:2,kmin_p  :kmax_p:2,ivar)=avrg1
      pr_buf(imin_p+1:imax_p:2,jmin_p  :jmax_p:2,kmin_p  :kmax_p:2,ivar)=avrg1
      pr_buf(imin_p  :imax_p:2,jmin_p+1:jmax_p:2,kmin_p  :kmax_p:2,ivar)=avrg1
      pr_buf(imin_p  :imax_p:2,jmin_p  :jmax_p:2,kmin_p+1:kmax_p:2,ivar)=avrg1
      pr_buf(imin_p+1:imax_p:2,jmin_p+1:jmax_p:2,kmin_p  :kmax_p:2,ivar)=avrg1
      pr_buf(imin_p+1:imax_p:2,jmin_p  :jmax_p:2,kmin_p+1:kmax_p:2,ivar)=avrg1
      pr_buf(imin_p  :imax_p:2,jmin_p+1:jmax_p:2,kmin_p+1:kmax_p:2,ivar)=avrg1
      pr_buf(imin_p+1:imax_p:2,jmin_p+1:jmax_p:2,kmin_p+1:kmax_p:2,ivar)=avrg1

    end subroutine prolong1

    !==========================================================================

    subroutine prolong2(iVar)

      ! Second order prolongation of Sol_BLK..Sol9_BLK or Sol_VGB into pr_buf

      integer, intent(in) :: ivar

      !------------------------------------------------------------------------
      ! Obtain qsol with an extra layer of ghost cells for calculating gradients
      if(present(Sol_VGB))then
         qsol(iminP-1:imaxP+1,jminP-1:jmaxP+1,kminP-1:kmaxP+1)=&
              Sol_VGB(iVar,iminP-1:imaxP+1,jminP-1:jmaxP+1,kminP-1:kmaxP+1,iBLK)
      else
         qsol(iminP-1:imaxP+1,jminP-1:jmaxP+1,kminP-1:kmaxP+1)=&
              Sol_BLK(iminP-1:imaxP+1,jminP-1:jmaxP+1,kminP-1:kmaxP+1,iBLK)
      end if

      ! The part without the extra ghost cells
      avrg=qsol(iminP:imaxP,jminP:jmaxP,kminP:kmaxP)

      ! Although we may have ghost cell values in directions/=idir, we don't 
      ! use those so that one directional gradient can be used without sending
      ! ghost cells in any other direction.

      select case(prolong_type)
      case('central','lr','minmod')
         ! Use first order extrapolation in the other directions
         if(idir/=1)then
            qsol(iminP-1,jminP:jmaxP,kminP:kmaxP)=&
                 qsol(iminP,jminP:jmaxP,kminP:kmaxP)
            qsol(imaxP+1,jminP:jmaxP,kminP:kmaxP)=&
                 qsol(imaxP,jminP:jmaxP,kminP:kmaxP)
         end if
         if(idir/=2)then
            qsol(iminP:imaxP,jminP-1,kminP:kmaxP)=&
                 qsol(iminP:imaxP,jminP,kminP:kmaxP)
            qsol(iminP:imaxP,jmaxP+1,kminP:kmaxP)=&
                 qsol(iminP:imaxP,jmaxP,kminP:kmaxP)
         end if
         if(idir/=3)then
            qsol(iminP:imaxP,jminP:jmaxP,kminP-1)=&
                 qsol(iminP:imaxP,jminP:jmaxP,kminP)
            qsol(iminP:imaxP,jminP:jmaxP,kmaxP+1)=&
                 qsol(iminP:imaxP,jminP:jmaxP,kmaxP)
         endif
      case('central2','lr2','minmod2','lr3')
         if(prolong_type/='lr3' .or. sendcorners)then
            ! Use second order extrapolation in the other directions
            if(idir/=1)then
               qsol(iminP-1,jminP:jmaxP,kminP:kmaxP)=&
                    2*qsol(iminP  ,jminP:jmaxP,kminP:kmaxP)&
                    -qsol(iminP+1,jminP:jmaxP,kminP:kmaxP)

               qsol(imaxP+1,jminP:jmaxP,kminP:kmaxP)=&
                    2*qsol(imaxP  ,jminP:jmaxP,kminP:kmaxP)&
                    -qsol(imaxP-1,jminP:jmaxP,kminP:kmaxP)
            end if
            if(idir/=2)then
               qsol(iminP:imaxP,jminP-1,kminP:kmaxP)=&
                    2*qsol(iminP:imaxP,jminP  ,kminP:kmaxP)&
                    -qsol(iminP:imaxP,jminP+1,kminP:kmaxP)

               qsol(iminP:imaxP,jmaxP+1,kminP:kmaxP)=&
                    2*qsol(iminP:imaxP,jmaxP  ,kminP:kmaxP)&
                    -qsol(iminP:imaxP,jmaxP-1,kminP:kmaxP)
            end if
            if(idir/=3)then
               qsol(iminP:imaxP,jminP:jmaxP,kminP-1)=&
                    2*qsol(iminP:imaxP,jminP:jmaxP,kminP  )&
                    -qsol(iminP:imaxP,jminP:jmaxP,kminP+1)

               qsol(iminP:imaxP,jminP:jmaxP,kmaxP+1)=&
                    2*qsol(iminP:imaxP,jminP:jmaxP,kmaxP  )&
                    -qsol(iminP:imaxP,jminP:jmaxP,kmaxP-1)
            endif
         end if
      case default
         call stop_mpi(&
              'Error 0 in prolong, unknown prolong_type='//prolong_type)
      end select ! prolong_type

      ! Calculate the gradients
      select case(prolong_type)
      case('central','central2')
         gradx=0.125*(&
              qsol(iminP+1:imaxP+1,jminP:jmaxP,kminP:kmaxP)-&
              qsol(iminP-1:imaxP-1,jminP:jmaxP,kminP:kmaxP))
         grady=0.125*(&
              qsol(iminP:imaxP,jminP+1:jmaxP+1,kminP:kmaxP)-&
              qsol(iminP:imaxP,jminP-1:jmaxP-1,kminP:kmaxP))
         gradz=0.125*(&
              qsol(iminP:imaxP,jminP:jmaxP,kminP+1:kmaxP+1)-&
              qsol(iminP:imaxP,jminP:jmaxP,kminP-1:kmaxP-1))
      case('lr','minmod','lr2','minmod2','lr3')
         gradxl=0.25*(&
              qsol(iminP  :imaxP  ,jminP:jmaxP,kminP:kmaxP)-&
              qsol(iminP-1:imaxP-1,jminP:jmaxP,kminP:kmaxP))
         gradxr=0.25*(&
              qsol(iminP+1:imaxP+1,jminP:jmaxP,kminP:kmaxP)-&
              qsol(iminP  :imaxP  ,jminP:jmaxP,kminP:kmaxP))
         gradyl=0.25*(&
              qsol(iminP:imaxP,jminP  :jmaxP  ,kminP:kmaxP)-&
              qsol(iminP:imaxP,jminP-1:jmaxP-1,kminP:kmaxP))
         gradyr=0.25*(&
              qsol(iminP:imaxP,jminP+1:jmaxP+1,kminP:kmaxP)-&
              qsol(iminP:imaxP,jminP  :jmaxP  ,kminP:kmaxP))
         gradzl=0.25*(&
              qsol(iminP:imaxP,jminP:jmaxP,kminP  :kmaxP  )-&
              qsol(iminP:imaxP,jminP:jmaxP,kminP-1:kmaxP-1))
         gradzr=0.25*(&
              qsol(iminP:imaxP,jminP:jmaxP,kminP+1:kmaxP+1)-&
              qsol(iminP:imaxP,jminP:jmaxP,kminP  :kmaxP  ))
         if(prolong_type=='minmod'.or.prolong_type=='minmod2')then
            gradx=sign(1.,gradxl)*max(0.,min(abs(gradxl),sign(1.,gradxl)*gradxr))
            grady=sign(1.,gradyl)*max(0.,min(abs(gradyl),sign(1.,gradyl)*gradyr))
            gradz=sign(1.,gradzl)*max(0.,min(abs(gradzl),sign(1.,gradzl)*gradzr))
         end if
      case default
         call stop_mpi('Error 1 in prolong, unknown prolong_type='//prolong_type)
      end select

      ! Apply gradients
      select case(prolong_type)
      case('central','central2','minmod','minmod2')
         pr_buf1(imin_p  :imax_p:2,jmin_p  :jmax_p:2,kmin_p  :kmax_p:2)=&
              avrg-gradx-grady-gradz
         pr_buf1(imin_p+1:imax_p:2,jmin_p  :jmax_p:2,kmin_p  :kmax_p:2)=&
              avrg+gradx-grady-gradz
         pr_buf1(imin_p  :imax_p:2,jmin_p+1:jmax_p:2,kmin_p  :kmax_p:2)=&
              avrg-gradx+grady-gradz
         pr_buf1(imin_p  :imax_p:2,jmin_p  :jmax_p:2,kmin_p+1:kmax_p:2)=&
              avrg-gradx-grady+gradz
         pr_buf1(imin_p+1:imax_p:2,jmin_p+1:jmax_p:2,kmin_p  :kmax_p:2)=&
              avrg+gradx+grady-gradz
         pr_buf1(imin_p+1:imax_p:2,jmin_p  :jmax_p:2,kmin_p+1:kmax_p:2)=&
              avrg+gradx-grady+gradz
         pr_buf1(imin_p  :imax_p:2,jmin_p+1:jmax_p:2,kmin_p+1:kmax_p:2)=&
              avrg-gradx+grady+gradz
         pr_buf1(imin_p+1:imax_p:2,jmin_p+1:jmax_p:2,kmin_p+1:kmax_p:2)=&
              avrg+gradx+grady+gradz
      case('lr','lr2','lr3')
         pr_buf1(imin_p  :imax_p:2,jmin_p  :jmax_p:2,kmin_p  :kmax_p:2)=&
              avrg-gradxl-gradyl-gradzl
         pr_buf1(imin_p+1:imax_p:2,jmin_p  :jmax_p:2,kmin_p  :kmax_p:2)=&
              avrg+gradxr-gradyl-gradzl
         pr_buf1(imin_p  :imax_p:2,jmin_p+1:jmax_p:2,kmin_p  :kmax_p:2)=&
              avrg-gradxl+gradyr-gradzl
         pr_buf1(imin_p  :imax_p:2,jmin_p  :jmax_p:2,kmin_p+1:kmax_p:2)=&
              avrg-gradxl-gradyl+gradzr
         pr_buf1(imin_p+1:imax_p:2,jmin_p+1:jmax_p:2,kmin_p  :kmax_p:2)=&
              avrg+gradxr+gradyr-gradzl
         pr_buf1(imin_p+1:imax_p:2,jmin_p  :jmax_p:2,kmin_p+1:kmax_p:2)=&
              avrg+gradxr-gradyl+gradzr
         pr_buf1(imin_p  :imax_p:2,jmin_p+1:jmax_p:2,kmin_p+1:kmax_p:2)=&
              avrg-gradxl+gradyr+gradzr
         pr_buf1(imin_p+1:imax_p:2,jmin_p+1:jmax_p:2,kmin_p+1:kmax_p:2)=&
              avrg+gradxr+gradyr+gradzr
      case default
         call stop_mpi('Error 2 in prolong, unknown prolong_type='//prolong_type)
      end select

      if(present(Sol_VGB))then
         pr_buf(iVar,:,:,:)=pr_buf1
      else
         pr_buf(:,:,:,iVar)=pr_buf1
      end if

    end subroutine prolong2

    !==========================================================================

    subroutine prolong1_state

      ! First order prolongation of Sol_VGB into pr_buf

      !------------------------------------------------------------------------
      ! First order prolongation
      avrg_state=Sol_VGB(:,iminP:imaxP,jminP:jmaxP,kminP:kmaxP,iBLK)

      pr_buf(:,imin_p  :imax_p:2,jmin_p  :jmax_p:2,kmin_p  :kmax_p:2)&
           = avrg_state
      pr_buf(:,imin_p+1:imax_p:2,jmin_p  :jmax_p:2,kmin_p  :kmax_p:2)&
           = avrg_state
      pr_buf(:,imin_p  :imax_p:2,jmin_p+1:jmax_p:2,kmin_p  :kmax_p:2)&
           = avrg_state
      pr_buf(:,imin_p  :imax_p:2,jmin_p  :jmax_p:2,kmin_p+1:kmax_p:2)&
           = avrg_state
      pr_buf(:,imin_p+1:imax_p:2,jmin_p+1:jmax_p:2,kmin_p  :kmax_p:2)&
           = avrg_state
      pr_buf(:,imin_p+1:imax_p:2,jmin_p  :jmax_p:2,kmin_p+1:kmax_p:2)&
           = avrg_state
      pr_buf(:,imin_p  :imax_p:2,jmin_p+1:jmax_p:2,kmin_p+1:kmax_p:2)&
           = avrg_state
      pr_buf(:,imin_p+1:imax_p:2,jmin_p+1:jmax_p:2,kmin_p+1:kmax_p:2)&
           = avrg_state

    end subroutine prolong1_state

    !==========================================================================

    subroutine prolong2_state

      ! Prolongation for 2nd order TVD resolution change.
      ! Put 2 coarse layers of Sol_VGB into pr_buf

      !------------------------------------------------------------------------
      ! First order prolongation
      avrg_state=Sol_VGB(:,iminP:imaxP,jminP:jmaxP,kminP:kmaxP,iBLK)

      select case(iFace)
      case(east_, west_)
         pr_buf(:,imin_p:imax_p,jmin_p  :jmax_p:2,kmin_p  :kmax_p:2) = avrg_state
         pr_buf(:,imin_p:imax_p,jmin_p+1:jmax_p:2,kmin_p  :kmax_p:2) = avrg_state
         pr_buf(:,imin_p:imax_p,jmin_p  :jmax_p:2,kmin_p+1:kmax_p:2) = avrg_state
         pr_buf(:,imin_p:imax_p,jmin_p+1:jmax_p:2,kmin_p+1:kmax_p:2) = avrg_state
      case(south_, north_)
         pr_buf(:,imin_p  :imax_p:2,jmin_p:jmax_p,kmin_p  :kmax_p:2) = avrg_state
         pr_buf(:,imin_p+1:imax_p:2,jmin_p:jmax_p,kmin_p  :kmax_p:2) = avrg_state
         pr_buf(:,imin_p  :imax_p:2,jmin_p:jmax_p,kmin_p+1:kmax_p:2) = avrg_state
         pr_buf(:,imin_p+1:imax_p:2,jmin_p:jmax_p,kmin_p+1:kmax_p:2) = avrg_state
      case(bot_, top_)
         pr_buf(:,imin_p  :imax_p:2,jmin_p  :jmax_p:2,kmin_p:kmax_p) = avrg_state
         pr_buf(:,imin_p+1:imax_p:2,jmin_p  :jmax_p:2,kmin_p:kmax_p) = avrg_state
         pr_buf(:,imin_p  :imax_p:2,jmin_p+1:jmax_p:2,kmin_p:kmax_p) = avrg_state
         pr_buf(:,imin_p+1:imax_p:2,jmin_p+1:jmax_p:2,kmin_p:kmax_p) = avrg_state
      end select

    end subroutine prolong2_state

    !==========================================================================
    subroutine send_equal

      ! Same level

      neiB=neiBLK(1,iface,iBLK)
      neiP=neiPE(1,iface,iBLK)
      if(unusedBlock_BP(neiB,neiP)) RETURN

      if(neiP==iProc)then
         ! Local copy
         if(present(Sol_VGB))then
            Sol_VGB(:,imin_g:imax_g,jmin_g:jmax_g,kmin_g:kmax_g,neiB)= &
                 Sol_VGB(:,imin_o:imax_o,jmin_o:jmax_o,kmin_o:kmax_o,iBLK)
         else
            sol_BLK(imin_g:imax_g,jmin_g:jmax_g,kmin_g:kmax_g,neiB)= &
                 sol_BLK(imin_o:imax_o,jmin_o:jmax_o,kmin_o:kmax_o,iBLK)
         endif
      else
         ! Remote send
         itag = 100*neiB+10*iface+1
         if(oktest.and.okdebug)write(*,*)'Remote equal send, me,iface,itag=',&
              iProc,iface,itag

         if(present(Sol_VGB))then
            eq_buf = Sol_VGB(:,imin_o:imax_o,jmin_o:jmax_o,kmin_o:kmax_o,iBLK)
         else
            eq_buf(:,:,:,1)=&
                 sol_BLK(imin_o:imax_o,jmin_o:jmax_o,kmin_o:kmax_o,iBLK)
         end if
         call MPI_Rsend(eq_buf,isize, MPI_REAL, neiP, itag, iComm, iError)
      end if
    end subroutine send_equal

    !==========================================================================
    subroutine send_restricted

      ! Neighbor is coarser

      neiB=neiBLK(1,iface,iBLK)
      neiP=neiPE(1,iface,iBLK)
      if (unusedBlock_BP(neiB,neiP)) RETURN

      ! Restrict the R range of cells into _r

      if(present(Sol_VGB))then
         call restrict_state
      else
         call restrict(1,sol_BLK)
      end if

      ! Subface index =1,2,3, or 4 with respect to the coarse neighbor
      ichild=BLKneighborCHILD(0,0,0,1,iBLK)
      if(ichild<1.or.ichild>8)then
         write(*,*)'me,iBLK,iface,neiL,neiP,neiB,ichild:',&
              iProc,iBLK,iface,neiL,neiP,neiB,ichild
         call stop_mpi('error in message_pass: ichild incorrect')
      end if
      isubface=child2subface(ichild,iface)
      if(isubface==0)then
         write(*,*)'me,iBLK,iface,neiL,neiP,neiB,ichild:',&
              iProc,iBLK,iface,neiL,neiP,neiB,ichild
         call stop_mpi('error in message_pass: isubface=0')
      end if
      if(neiP==iProc)then
         ! Local copy into appropriate subface
         call setsubrange_g
         if(present(Sol_VGB))then
            Sol_VGB(:,imin_s:imax_s,jmin_s:jmax_s,kmin_s:kmax_s,neiB)=re_buf
         else
            sol_BLK(imin_s:imax_s,jmin_s:jmax_s,kmin_s:kmax_s,neiB)=&
                 re_buf(:,:,:,1)
         end if
      else
         ! Remote send
         itag = 100*neiB+10*iface+isubface
         if(oktest.and.okdebug)write(*,*)&
              'Remote restricted send, me,iface,itag=',&
              iProc,iface,itag
         call MPI_Rsend(re_buf, isize_r, MPI_REAL, neiP, itag,iComm,iError)
      end if
    end subroutine send_restricted

    !==========================================================================
    subroutine send_prolonged

      integer :: iVar

      ! Neighbor is finer, so prolong the _o range of cells into _p
      if(prolongorder > 1)then
         ! Unoptimized 2nd order prolongation
         do iVar=1,nVar
            call prolong2(iVar)
         end do
      else if(present(Sol_VGB))then
         if(nCoarseLayer == 1)then
            ! Optimized 1st order prolongation
            call prolong1_state
         else
            ! Prolong 2 coarse layers for 2nd order TVD resolution change
            call prolong2_state
         end if
      else
         ! Variable by variable 1st order prolongation
         call prolong1(1,sol_BLK)
      end if

      do isubface=1,4

         neiB=neiBLK(isubface,iface,iBLK)
         neiP=neiPE(isubface,iface,iBLK)
         if(unusedBlock_BP(neiB,neiP)) CYCLE

         call setsubrange_p
         if(neiP==iProc)then
            ! Local copy of appropriate subface
            if(present(Sol_VGB))then
               Sol_VGB(:,imin_g:imax_g,jmin_g:jmax_g,kmin_g:kmax_g,neiB)=&
                    pr_buf(:,imin_s:imax_s,jmin_s:jmax_s,kmin_s:kmax_s)
            else
               sol_BLK(imin_g:imax_g,jmin_g:jmax_g,kmin_g:kmax_g,neiB)=&
                    pr_buf(imin_s:imax_s,jmin_s:jmax_s,kmin_s:kmax_s,1)
            end if
         else
            ! Remote send
            itag = 100*neiB+10*iface+1
            if(oktest.and.okdebug)write(*,*)&
                 'Remote prolong send, me,iface,itag=',&
                 iProc,iface,itag

            if(present(Sol_VGB))then
               pr_buf_s=pr_buf(1:nvar,imin_s:imax_s,jmin_s:jmax_s,kmin_s:kmax_s)
            else
               pr_buf_s=pr_buf(imin_s:imax_s,jmin_s:jmax_s,kmin_s:kmax_s,1:nvar)
            end if
            call MPI_Rsend(pr_buf_s,isize,MPI_REAL,neiP,itag,iComm,iError)
         end if
      end do

    end subroutine send_prolonged

    !==========================================================================
    subroutine remote_receive(nsubface, isize1)

      ! Receive nsubface subfaces of size isize1

      integer, intent(in) :: nsubface, isize1
      !------------------------------------------------------------------------

      do isubface=1,nsubface

         neiP=neiPE(isubface,otherface,iBLK)
         if(neiP==iProc) CYCLE

         neiB=neiBLK(isubface,otherface,iBLK)
         if(unusedBlock_BP(neiB,neiP)) CYCLE

         ! Remote receive
         itag = 100*iBLK+10*iface+isubface
         if(oktest.and.okdebug)write(*,*)&
              'Remote recieve, me,iface,itag,neiL=',&
              iProc,iface,itag,neiLEV(otherface,iBLK)

         ibuf = (iBLK-1)*isize+(isubface-1)*isize_r+1
         call MPI_irecv(buffer(ibuf,iside), isize1, MPI_REAL, &
              neiP, itag, iComm, request, iError)
         number_receive_requests = number_receive_requests + 1
         receive_requests(number_receive_requests) = request
      end do

    end subroutine remote_receive

    !==========================================================================
    subroutine buffer2face(&
         buf,imin_g,imax_g,jmin_g,jmax_g,kmin_g,kmax_g,nvar)

      integer, intent(in) :: imin_g,imax_g,jmin_g,jmax_g,kmin_g,kmax_g,nvar
      real, dimension(imin_g:imax_g,jmin_g:jmax_g,kmin_g:kmax_g,nvar),&
           intent(inout) :: buf
      !------------------------------------------------------------------------

      ! Read buffer into 4 index variables
      sol_BLK(imin_g:imax_g,jmin_g:jmax_g,kmin_g:kmax_g,iBLK)=buf(:,:,:,1)

    end subroutine buffer2face

    !==========================================================================
    subroutine buffer2face_state(&
         buf,imin_g,imax_g,jmin_g,jmax_g,kmin_g,kmax_g,nvar)

      integer, intent(in) :: imin_g,imax_g,jmin_g,jmax_g,kmin_g,kmax_g,nvar
      real, dimension(nVar,imin_g:imax_g,jmin_g:jmax_g,kmin_g:kmax_g),&
           intent(inout) :: buf
      !------------------------------------------------------------------------

      ! Read buffer into Sol_VGB

      Sol_VGB(:,imin_g:imax_g,jmin_g:jmax_g,kmin_g:kmax_g,iBLK)=buf

    end subroutine buffer2face_state

    !==========================================================================
    subroutine buffer2subfaces(&
         buf,imin_r,imax_r,jmin_r,jmax_r,kmin_r,kmax_r,nvar)

      integer, intent(in) :: imin_r,imax_r,jmin_r,jmax_r,kmin_r,kmax_r,nvar
      real, dimension(imin_r:imax_r,jmin_r:jmax_r,kmin_r:kmax_r,nvar,4),&
           intent(inout) :: buf
      !------------------------------------------------------------------------

      ! Loop over 4 subfaces to read restricted values
      do isubface=1,4
         neiP = neiPE(isubface,otherface,iBLK)
         if(neiP==iProc) CYCLE
         if (unusedBlock_BP(neiBLK(isubface,otherface,iBLK),neiP)) CYCLE
         call setsubrange_g

         sol_BLK(imin_s:imax_s,jmin_s:jmax_s,kmin_s:kmax_s,iBLK)=&
              buf(:,:,:,1,isubface)
      end do ! isubface

    end subroutine buffer2subfaces

    !==========================================================================
    subroutine buffer2subfaces_state(&
         buf,imin_r,imax_r,jmin_r,jmax_r,kmin_r,kmax_r,nvar)

      integer, intent(in) :: imin_r,imax_r,jmin_r,jmax_r,kmin_r,kmax_r,nvar
      real, dimension(nVar,imin_r:imax_r,jmin_r:jmax_r,kmin_r:kmax_r,4),&
           intent(inout) :: buf
      !------------------------------------------------------------------------

      ! Loop over 4 subfaces to read restricted values
      do isubface=1,4
         neiP = neiPE(isubface,otherface,iBLK)
         if(neiP==iProc) CYCLE
         if (unusedBlock_BP(neiBLK(isubface,otherface,iBLK),neiP)) CYCLE
         call setsubrange_g

         Sol_VGB(:,imin_s:imax_s,jmin_s:jmax_s,kmin_s:kmax_s,iBLK)=&
              buf(:,:,:,:,isubface)
      end do ! isubface

    end subroutine buffer2subfaces_state
    !==========================================================================

  end subroutine message_pass_dir

end module ModMessagePass
