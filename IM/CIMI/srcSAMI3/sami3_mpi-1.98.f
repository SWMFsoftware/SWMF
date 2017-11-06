


!     *******************************************
!     *******************************************
 
!                  SAMI3_MPI-1.98

!     *******************************************
!     *******************************************
!
!      program main
      subroutine sami_init
      use ModMpi
      use ModSAMI
      use ModCoupleCimi,ONLY:DoCoupleCimi,set_sami_grid_for_mod

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc' 

c     Some local variables

!       real fism(linesuv)
      ! these look to be temporary variables for passing MPI 
       real denitmp(nz,nf,nl), titmp(nz,nf,nl), vsitmp(nz,nf,nl)

!       real deni_mnptmp(nz,nion),ti_mnptmp(nz,nion),te_mnptmp(nz)

       real denntmp(nz,nf,nl)

!!       real deni_inp(nz,nf,nion),ti_inp(nz,nf,nion),
!!     .        vsi_inp(nz,nf,nion),denn_inp(nz,nf,nion)
!!       real te_inp(nz,nf)
!!       real te_inp(nz,nf),tn_inp(nz,nf)
!!!       real phi(nnx,nny),phialt(nnx,nny),philon(nnx,nny),p_crit(nnx-1)


!!!       logical tflag,ttflag

!!!c allocatable total matrices
!!!
!!!! Total 
!!!c
!!!!      real denit(nz,nf,nlt,nion)
!!!
!!!! Total 
!!!
!!!       real, dimension(:,:,:,:), allocatable :: denit,dennt
!!!!      real denit(nz,nf,nlt,nion),dennt(nz,nf,nlt,nneut)
!!!       real, dimension(:,:,:,:), allocatable :: vsit, sumvsit
!!!!      real vsit(nz,nf,nlt,nion),sumvsit(nz,nf,nlt,nion)
!!!       real, dimension(:,:,:,:), allocatable :: tit
!!!!      real tit(nz,nf,nlt,nion)
!!!
!!!       real, dimension(:,:,:), allocatable :: ut, vt, vpit, net
!!!!      real ut(nz,nf,nlt),vt(nz,nf,nlt),vpit(nz,nf,nlt),net(nz,nf,nlt)
!!!       real, dimension(:,:,:), allocatable :: tet, tnt
!!!!      real tet(nz,nf,nlt),tnt(nz,nf,nlt)
!!!
!!!!     height integrated pedersen/hall conductivities
!!!       real, dimension(:,:,:), allocatable :: u1t, u2t, u3t, u4t
!!!!      real u1t(nz,nf,nlt),u2t(nz,nf,nlt),u3t(nz,nf,nlt),u4t(nz,nf,nlt)
!!!
!!!       real, dimension(:,:,:), allocatable :: vnqt, vnpt, vnphit
!!!!      real vnqt(nz,nf,nlt),vnpt(nz,nf,nlt),vnphit(nz,nf,nlt)
!!!
!!!       real, dimension(:,:,:), allocatable :: jpt, jphit
!!!!      real jpt(nz,nf,nlt),jphit(nz,nf,nlt)
!!!
!!!       real, dimension(:,:,:), allocatable :: u1pt, u2st, u3ht
!!!
!!!       real, dimension(:,:,:), allocatable :: sigmapict,sigmahict
!!!
!!!       real, dimension(:,:,:), allocatable :: sigmapt,sigmaht


!!!c Output matrices for restart
!!!
!!!       real deniout(nz,nf,nl,nion,numwork),  
!!!     &      tiout(nz,nf,nl,nion,numwork),  
!!!     &      vsiout(nz,nf,nl,nion,numwork) 
!!!       real teout(nz,nf,nl,numwork) 
!!!       real*8 dphi(nnx+1,nnyt)
 

c     Begin MPI stuff

!      include 'mpif.h'
      integer status(MPI_STATUS_SIZE)
      
      integer :: nError
C ************************ initializations ***********************************
C Find out how many tasks are in this partition and what my task id is.  Then
C define the number of worker tasks and the array partition size as chunksize.
C Note:  For this example, the MP_PROCS environment variable should be set
C to an odd number...to insure even distribution of the array to numtasks-1
C worker tasks.
C *****************************************************************************

!!!      call mpi_init(ierr)
!!!      call mpi_comm_rank(iComm, taskid, ierr)
!!!      call mpi_comm_size(iComm, numtasks, ierr)
      numtasks = nProc
      taskid = iProc
      write(*,*)'taskid =',taskid

      allocate(phi(nnx,nny),phialt(nnx,nny),philon(nnx,nny),p_crit(nnx-1))



      numworkers = numtasks-1

!output matrices for restart
      allocate(deniout(nz,nf,nl,nion,numwork))
      allocate(tiout(nz,nf,nl,nion,numwork))
      allocate(vsiout(nz,nf,nl,nion,numwork))
      allocate(teout(nz,nf,nl,numwork))
      allocate(dphi(nnx+1,nnyt))


c Check to see if the number of processors selected agrees with 
c the number of divisions in params

      if(taskid .eq. 0) then
      if(numwork .ne. numworkers) then
         print *, ' numworkers is ',numworkers
         print *, ' numwork (in param3_mpi) is',numwork
         print *, ' in order for the code to work correctly '
         print *, ' these two numbers must be the same '
         print *, ' Either set np = numwork +1 and rerun or '
         print *, ' change numwork and recompile '
         call mpi_abort(iComm, nError, ierr)
         call mpi_finalize(ierr)

      endif
      endif
c
c Determine what is left (down) and right (up)
c Here we assume that taskid=0 is the Master and does nothing but
c deal with handling the data
c

!      call getenva("HOSTNAME",strng)

!      call system ("echo $HOSTNAME")
      print *, ' taskid ',taskid

      if(taskid .eq. numtasks -1) then
         right = 1
      else
         right = taskid +1
      endif
      if(taskid .eq. 1) then
         left = numtasks -1
      else
         left = taskid -1
      endif

! open daily fism file

!      if ( taskid .eq. 0 ) then
!        open(unit=402,file='fism_daily.inp ',form='unformatted')
!        read(402) fism
!        close(402)
!      endif

! open input files
! only need these on the master

      if(taskid .eq. 0) then
         open ( unit=10, file='sami3_mpi-1.98.namelist'  )
         open ( unit=20, file='deni-init.inp'        )
         open ( unit=30, file='ichem.inp'            )
         open ( unit=50, file='phabsdt_euvac.inp'    )  ! euvac
         open ( unit=60, file='phiondt_euvac.inp'    )  ! euvac
         open ( unit=61, file='phionnt.inp'          )
         open ( unit=65, file='euvflux_euvac.inp'    )  ! euvac
         open ( unit=66, file='thetant.inp'          )
         open ( unit=67, file='zaltnt.inp'           )
      endif

      call initial(phialt,philon,p_crit)
      
c We are out of initial now.  
c We will overwrite the values of
c deni, vsi, ti, te if this is a restart (restart = true)

      if(restart) then
         if(taskid .eq. 0) then
            print *,'doing restart'
            open ( unit=210, file='time.rst', form='formatted' )
            open ( unit=211, file='deni.rst', form='unformatted' )
            open ( unit=212, file='vsi.rst', form='unformatted' )
            open ( unit=213, file='ti.rst', form='unformatted' )
            open ( unit=214, file='te.rst', form='unformatted' )

            read(210,*) hrinit
            read(211) deniout
            read(212) vsiout
            read(213) tiout
            read(214) teout

            close (210)
            close (211)
            close (212)
            close (213)
            close (214)

            do iwrk = 1,numworkers
               do nntmp = 1,nion
                  do ktmp = 1,nl
                     do jtmp = 1,nf
                        do itmp = 1,nz
                           denitmp(itmp,jtmp,ktmp)
     .                          =  deniout(itmp,jtmp,ktmp,nntmp,iwrk) 
                           titmp(itmp,jtmp,ktmp)
     .                          =  tiout(itmp,jtmp,ktmp,nntmp,iwrk) 
                           vsitmp(itmp,jtmp,ktmp)
     .                          =  vsiout(itmp,jtmp,ktmp,nntmp,iwrk) 
                        enddo
                     enddo
                  enddo

                  call mpi_send(denitmp, nz*nf*nl, MPI_REAL, iwrk, 9, 
     .                 iComm, ierr)
                  call mpi_send(titmp, nz*nf*nl, MPI_REAL, iwrk, 9, 
     .                 iComm, ierr)
                  call mpi_send(vsitmp, nz*nf*nl, MPI_REAL, iwrk, 9, 
     .                 iComm, ierr)
               enddo

               do ktmp = 1,nl
                  do jtmp = 1,nf
                     do itmp = 1,nz
                        te(itmp,jtmp,ktmp)
     .                       =  teout(itmp,jtmp,ktmp,iwrk) 
                     enddo
                  enddo
               enddo
               call mpi_send(te, nz*nf*nl, MPI_REAL, iwrk, 9, 
     .              iComm, ierr)

            enddo
         endif

c Now let's get those files

         if(taskid .gt. 0 .and. taskid .le. numworkers) then

            do nntmp = 1,nion
               call mpi_recv(denitmp, nz*nf*nl, MPI_REAL, 0, 9, 
     .              iComm, status, ierr)
               call mpi_recv(titmp, nz*nf*nl, MPI_REAL, 0, 9, 
     .              iComm, status, ierr)
               call mpi_recv(vsitmp, nz*nf*nl, MPI_REAL, 0, 9, 
     .              iComm, status, ierr)
               do ktmp = 1,nl
                  do jtmp = 1,nf
                     do itmp = 1,nz
                        deni(itmp,jtmp,ktmp,nntmp)
     .                       =  denitmp(itmp,jtmp,ktmp) 
                        ti(itmp,jtmp,ktmp,nntmp)
     .                       =  titmp(itmp,jtmp,ktmp) 
                        vsi(itmp,jtmp,ktmp,nntmp)
     .                       =  vsitmp(itmp,jtmp,ktmp) 
                     enddo
                  enddo
               enddo
            enddo
            call mpi_recv(te, nz*nf*nl, MPI_REAL, 0, 9, 
     .           iComm, status, ierr)

         endif

c tell the workers the starting time
c this call has to be seen by the master and workers

        call mpi_bcast(hrinit,1,MPI_REAL,0,iComm,ierr)

      endif

! open output files

      if(taskid .eq. 0) then
         if ( fmtout ) then
            call open_f
         else
            call open_u
         endif
      endif 

      if(taskid .eq. 0) then
         close (10)
         close (20)
         close (30)
         close (50)
         close (60)
         close (61)
         close (65)
         close (66)
         close (67)
         close (68)
      endif

      ! when coupling to cimi make sure coupling module knows sami's grid
      if (DoCoupleCimi) then
         call set_sami_grid_for_mod(nzp1,nfp1,nlt,nnx,nny,blatpt,blonpt)
      endif
      end !end of sami_init
!===============================================================================

      subroutine sami_run(DtAdvance)
      use ModMpi
      use ModSAMI
      use ModCoupleCimi,ONLY:DoCoupleCimi,PotCimiOnSamiGrid_C
     &     ,iStartTime_I
      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc' 
      
      ! input representing how much sami_run advances
      real, intent(in) :: DtAdvance 

c     Some local variables

!       real fism(linesuv)


       real denitmp(nz,nf,nl), titmp(nz,nf,nl), vsitmp(nz,nf,nl)

       real deni_mnptmp(nz,nion),ti_mnptmp(nz,nion),te_mnptmp(nz)

       real denntmp(nz,nf,nl)

!!       real deni_inp(nz,nf,nion),ti_inp(nz,nf,nion),
!!     .        vsi_inp(nz,nf,nion),denn_inp(nz,nf,nion)
!!       real te_inp(nz,nf)
!!       real te_inp(nz,nf),tn_inp(nz,nf)
       real :: DeltaTime, DtStop
       real, save :: DtCourant
       real, allocatable :: DtCourant_I(:)
       real, parameter :: StopTolerance = 1e-5
       logical, save :: IsFirstCall = .true.,IsFirstCallMaster = .true.
       logical, save :: IsFirstCallWorker = .true.
       logical tflag,ttflag,lprnt

       integer status(MPI_STATUS_SIZE)
       
       integer :: nError


******************** master task *******************************************
       if (IsFirstCall) then
          DtCourant = Dt0
          IsFirstCall = .false.
       endif
       Dt = DtCourant
       DeltaTime = 0.0 
       DtStop    = DtAdvance
       
       if(DoCoupleCimi) then
          !temporarily remove this statement
!     force hrinit to be consistant with CIMI start time
!          hrinit = real(iStartTime_I(4))+real(iStartTime_I(5))/60.0
!     &         +real(iStartTime_I(6))/3600.0
       endif
       
       if(taskid .eq. 0) then
          if(IsFirstCallMaster) then
!     allocate the total matrices only on master
             
             allocate (denit(nz,nf,nlt,nion),dennt(nz,nf,nlt,nneut))
             allocate  (vsit(nz,nf,nlt,nion),sumvsit(nz,nf,nlt,nion))
             allocate  (tit(nz,nf,nlt,nion))
             allocate (ut(nz,nf,nlt),vt(nz,nf,nlt),vpit(nz,nf,nlt),
     .            net(nz,nf,nlt))
             allocate (tet(nz,nf,nlt),tnt(nz,nf,nlt))
             allocate (u1t(nz,nf,nlt),u2t(nz,nf,nlt),u3t(nz,nf,nlt),
     .            u4t(nz,nf,nlt))
             allocate (vnqt(nz,nf,nlt),vnpt(nz,nf,nlt),vnphit(nz,nf,nlt))
             allocate (jpt(nz,nf,nlt),jphit(nz,nf,nlt))
             allocate (u1pt(nz,nf,nlt),u2st(nz,nf,nlt),u3ht(nz,nf,nlt))
             allocate (sigmapict(nz,nf,nlt),sigmahict(nz,nf,nlt))
             allocate (sigmapt(nz,nf,nlt),sigmaht(nz,nf,nlt))
             allocate (DtCourant_I(numworkers))
             hrut    = hrinit
             timemax = hrmax * sphr
             istep   = 0
             tprnt   = 0.
             tneut   = 0.
             time    = 0.
             ntm     = 0
!     ieuv    = 1
             print *,'hrinit = ',hrinit
             print *,'max',max(dt0/3600.,dthr)
             ntmmax  = min(maxstep, int( hrmax / max(dt0/3600.,dthr) ))
             print *,'ntmmax ',ntmmax
             IsFirstCallMaster=.false.
          endif


         ifintot = numworkers
         ifintot1 = numworkers
         ifintot2 = numworkers

         tflag  = .true.
         icnt10 =  0

!         do while (      time .le. timemax  
!     .        .and. istep .le. maxstep     )
         
         do while ( tflag )
             !if(taskid == 0) write(*,*) 'Starting timeloop', deltatime,dt
            do iwrk = 1,numworkers
               call mpi_iprobe(iwrk,10,iComm,flagit10,
     .              status,ierr)
               if (flagit10) then 
                  icnt10 = icnt10 + 1
                  call mpi_recv(xxx,1,MPI_REAL,iwrk,10,
     .                 iComm,status,ierr)
               endif
               if (icnt10 .eq. numworkers) tflag=.false.
            enddo

C Now wait to receive back the results from each worker task
            
            do  iwrk = 1, numworkers
               source = iwrk
               dest = source

               call mpi_iprobe(source, 2, 
     .              iComm, flagit, status, ierr)
               
               if(flagit .and. ifintot2 .gt. 0) then
                  
                  call mpi_recv(hipcp, nf*nl, MPI_REAL, iwrk, 2, 
     .                 iComm, status, ierr)
                  call mpi_recv(hihcm, nf*nl, MPI_REAL, iwrk, 2, 
     .                 iComm, status, ierr)
                  call mpi_recv(hipcphi, nf*nl, MPI_REAL, iwrk, 2, 
     .                 iComm, status, ierr)
                  call mpi_recv(hidphig, nf*nl, MPI_REAL, iwrk, 2, 
     .                 iComm, status, ierr)
                  call mpi_recv(hidpg, nf*nl, MPI_REAL, iwrk, 2, 
     .                 iComm, status, ierr)
                  call mpi_recv(hidphiv, nf*nl, MPI_REAL, iwrk, 2, 
     .                 iComm, status, ierr)
                  call mpi_recv(hidpv, nf*nl, MPI_REAL, iwrk, 2, 
     .                 iComm, status, ierr)
                  call mpi_recv(hipc, nf*nl, MPI_REAL, iwrk, 2, 
     .                 iComm, status, ierr)
                  call mpi_recv(hihc, nf*nl, MPI_REAL, iwrk, 2, 
     .                 iComm, status, ierr)
                  call mpi_recv(hidv, nf*nl, MPI_REAL, iwrk, 2, 
     .                 iComm, status, ierr)
                  
                  call mpi_recv(hrut, 1, MPI_REAL, iwrk, 2, 
     .                 iComm, status, ierr)
                  
                  do k = 2,nl-1
                     kk = (iwrk-1)*(nl -2) + k - 1
                     if(kk .eq. 0) kk = nlt
                     if(kk .eq. nltp1) kk = 1
                     do j = 1,nf
                        hipcpt(j,kk)   = hipcp(j,k)
                        hihcmt(j,kk)   = hihcm(j,k)
                        hipcphit(j,kk) = hipcphi(j,k)
                        hidphigt(j,kk) = hidphig(j,k)
                        hidpgt(j,kk)   = hidpg(j,k)
                        hidphivt(j,kk) = hidphiv(j,k)
                        hidpvt(j,kk)   = hidpv(j,k)
                        hipct(j,kk)    = hipc(j,k)
                        hihct(j,kk)    = hihc(j,k)
                        hidvt(j,kk)    = hidv(j,k)
                     enddo
                  enddo
                  
                  
                  ifintot2 = ifintot2 - 1
                  
               endif
               
               if ( ifintot2 .eq. 0 ) then
                  ifintot2 = numworkers
                  call potpphi(phi,phialt,philon,dphi,hrut,p_crit)
                  do jwrk = 1,numworkers
                     call mpi_send(phi,nnx*nny,MPI_REAL,jwrk,3,
     .                    iComm,ierr)
                     call mpi_send(phialt,nnx*nny,MPI_REAL,jwrk,3,
     .                    iComm,ierr)
                     call mpi_send(philon,nnx*nny,MPI_REAL,jwrk,3,
     .                    iComm,ierr)
                     call mpi_send(p_crit,nnx-1,MPI_REAL,jwrk,3,
     .                    iComm,ierr)
                  enddo
!     print *,'out of sends in potential'
               endif
               
               call mpi_iprobe(source, 0, 
     .              iComm, flagit, status, ierr)
               
               if(flagit .and. ifintot .gt. 0) then
                  
!     This is just for outputting the data
!     only sent as often as data dumps are requested
                  
                  call mpi_recv(time, 1, MPI_REAL, iwrk, 0, 
     .                 iComm, status, ierr)
                  call mpi_recv(hrut, 1, MPI_REAL, iwrk, 0, 
     .                 iComm, status, ierr)
                  call mpi_recv(istep, 1, MPI_INTEGER, iwrk, 0, 
     .                 iComm, status, ierr)
                  do nntmp = 1,nion
                     call mpi_recv(denitmp, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                    iComm, status, ierr)
                     call mpi_recv(denntmp, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                    iComm, status, ierr)
                     call mpi_recv(titmp, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                    iComm, status, ierr)
                     call mpi_recv(vsitmp, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                    iComm, status, ierr)
                     do itmp = 1,nz
                        do jtmp = 1,nf
                           do ktmp = 1,nl
                              deni(itmp,jtmp,ktmp,nntmp)
     .                             =  denitmp(itmp,jtmp,ktmp) 
                              denn(itmp,jtmp,ktmp,nntmp)
     .                             =  denntmp(itmp,jtmp,ktmp) 
                              ti(itmp,jtmp,ktmp,nntmp)
     .                             =  titmp(itmp,jtmp,ktmp) 
                              vsi(itmp,jtmp,ktmp,nntmp)
     .                             =  vsitmp(itmp,jtmp,ktmp) 
                           enddo
                        enddo
                     enddo
                  enddo
                  call mpi_recv(te, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                 iComm, status, ierr)
                  call mpi_recv(u1p, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                 iComm, status, ierr)
                  call mpi_recv(u2s, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                 iComm, status, ierr)
                  call mpi_recv(u3h, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                 iComm, status, ierr)
                  call mpi_recv(u1, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                 iComm, status, ierr)
                  call mpi_recv(u2, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                 iComm, status, ierr)
                  call mpi_recv(u3, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                 iComm, status, ierr)
                  call mpi_recv(u4, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                 iComm, status, ierr)
                  call mpi_recv(sigmap, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                 iComm, status, ierr)
                  call mpi_recv(sigmah, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                 iComm, status, ierr)
                  call mpi_recv(sigmapic, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                 iComm, status, ierr)
                  call mpi_recv(sigmahic, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                 iComm, status, ierr)
 !                 call mpi_recv(hipcp, nf*nl, MPI_REAL, iwrk, 0, 
 !    .                 iComm, status, ierr)
 !                 call mpi_recv(hipcphi, nf*nl, MPI_REAL, iwrk, 0, 
 !    .                 iComm, status, ierr)
 !                 call mpi_recv(hihcm, nf*nl, MPI_REAL, iwrk, 0, 
 !    .                 iComm, status, ierr)
 !                 call mpi_recv(hidpv, nf*nl, MPI_REAL, iwrk, 0, 
 !    .                 iComm, status, ierr)
 !                 call mpi_recv(hidphiv, nf*nl, MPI_REAL, iwrk, 0, 
 !    .                 iComm, status, ierr)
 !                 call mpi_recv(hidpg, nf*nl, MPI_REAL, iwrk, 0, 
 !    .                 iComm, status, ierr)
 !                 call mpi_recv(hidphig, nf*nl, MPI_REAL, iwrk, 0, 
 !    .                 iComm, status, ierr)
                  call mpi_recv(vnq, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                 iComm, status, ierr)
                  call mpi_recv(vnp, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                 iComm, status, ierr)
                  call mpi_recv(vnphi, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                 iComm, status, ierr)
                  call mpi_recv(jp, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                 iComm, status, ierr)
                  call mpi_recv(jphi, nz*nf*nl, MPI_REAL, iwrk, 0, 
     .                 iComm, status, ierr)
!                  call mpi_recv(hipc, nf*nl, MPI_REAL, iwrk, 0, 
!     .                 iComm, status, ierr)
!                  call mpi_recv(hihc, nf*nl, MPI_REAL, iwrk, 0, 
!     .                 iComm, status, ierr)
!                  call mpi_recv(hidv, nf*nl, MPI_REAL, iwrk, 0, 
!     .                 iComm, status, ierr)

! Put the submatrices into the correct matrix

                  do nn = 1,nion
                     do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) +k -1
                        if(kk .eq. 0) kk = nlt
                        if(kk .eq. nltp1) kk = 1
                        do j = 1,nf
                           do i = 1,nz
                              denit(i,j,kk,nn) = deni(i,j,k,nn)
                              dennt(i,j,kk,nn) = denn(i,j,k,nn)
                              tit(i,j,kk,nn) = ti(i,j,k,nn)
                              vsit(i,j,kk,nn) = vsi(i,j,k,nn)
                           enddo
                        enddo
                     enddo
                     
! Put the submatrices into the total matrix for restart
                     
                     do k = 1,nl
                        do j = 1,nf
                           do i = 1,nz
                              deniout(i,j,k,nn,iwrk) = deni(i,j,k,nn)
                              tiout(i,j,k,nn,iwrk) = ti(i,j,k,nn)
                              vsiout(i,j,k,nn,iwrk) = vsi(i,j,k,nn)
                           enddo
                        enddo
                     enddo
                  enddo  ! for nion loop
                  
                  do k = 1,nl
                     do j = 1,nf
                        do i = 1,nz
                           teout(i,j,k,iwrk) = te(i,j,k)
                        enddo
                     enddo
                  enddo  
                  
                  do k = 2,nl-1
                     kk = (iwrk-1)*(nl -2) +k -1
                     if(kk .eq. 0) kk = nlt
                     if(kk .eq. nltp1) kk = 1
                     do j = 1,nf
                        do i = 1,nz
                           tet(i,j,kk)       = te(i,j,k)
                           u1pt(i,j,kk)      = u1p(i,j,k)
                           u2st(i,j,kk)      = u2s(i,j,k)
                           u3ht(i,j,kk)      = u3h(i,j,k)
                           u1t(i,j,kk)       = u1(i,j,k)
                           u2t(i,j,kk)       = u2(i,j,k)
                           u3t(i,j,kk)       = u3(i,j,k)
                           u4t(i,j,kk)       = u4(i,j,k)
                           sigmapt(i,j,kk)   = sigmap(i,j,k)
                           sigmaht(i,j,kk)   = sigmah(i,j,k)
                           sigmapict(i,j,kk) = sigmapic(i,j,k)
                           sigmahict(i,j,kk) = sigmahic(i,j,k)
                           vnqt(i,j,kk)      = vnq(i,j,k)
                           vnpt(i,j,kk)      = vnp(i,j,k)
                           vnphit(i,j,kk)    = vnphi(i,j,k)
                           jpt(i,j,kk)       = jp(i,j,k)
                           jphit(i,j,kk)     = jphi(i,j,k)
                        enddo
                     enddo
                  enddo
                  
!                  do k = 2,nl-1
!                     kk = (iwrk-1)*(nl -2) +k -1
!                     if(kk .eq. 0) kk = nlt
!                     if(kk .eq. nltp1) kk = 1
!                     do j = 1,nf
!                       hipcpt(j,kk)   = hipcp(j,k)
!                       hipcphit(j,kk) = hipcphi(j,k)
!                       hihcmt(j,kk)   = hihcm(j,k)
!                       hidpvt(j,kk)   = hidpv(j,k)
!                       hidphivt(j,kk) = hidphiv(j,k)
!                       hidpgt(j,kk)   = hidpg(j,k)
!                       hidphigt(j,kk) = hidphig(j,k)
!                       hipct(j,kk)    = hipc(j,k)
!                       hihct(j,kk)    = hihc(j,k)
!                       hidvt(j,kk)    = hidv(j,k)
!                     enddo
!                  enddo

                  ifintot = ifintot -1                  
                  
               endif

               call mpi_iprobe(source, 1, 
     .              iComm, flagit1, status, ierr)
               
               if(flagit1 .and. ifintot1 .gt. 0) then
                  call mpi_recv(dtmp, 1, MPI_REAL, iwrk, 1, 
     .                 iComm, status, ierr)
                  DtCourant_I(iwrk) = dtmp

                  call mpi_recv(DeltaTime, 1, MPI_REAL, iwrk, 1, 
     .                 iComm, status, ierr)
                  
                  call mpi_recv(time, 1, MPI_REAL, iwrk, 1, 
     .                 iComm, status, ierr)

                  call mpi_recv(istep, 1, MPI_INTEGER, iwrk, 1, 
     .                 iComm, status, ierr)

        call mpi_recv(deni_mnptmp,nz*nion,MPI_REAL,iwrk,1,
     .          iComm, status, ierr)
        call mpi_recv(ti_mnptmp,nz*nion,MPI_REAL,iwrk,1,
     .          iComm, status, ierr)
        call mpi_recv(te_mnptmp,nz,MPI_REAL,iwrk,1,
     .          iComm, status, ierr)

              if ( ifintot1 .eq. numworkers ) then
                 do ni = nion1,nion2
                   do i = 1,nz
                     deni_mnp(i,ni) = 0.
                     ti_mnp(i,ni)   = 0.
                   enddo
                 enddo

                 do i = 1,nz
                   te_mnp(i)     = 0.
                 enddo
              endif

            do ni = nion1,nion2
              do i = 1,nz
                deni_mnp(i,ni) = deni_mnp(i,ni) + 
     .                           deni_mnptmp(i,ni)/numworkers
                ti_mnp(i,ni)   = ti_mnp(i,ni) + 
     .                           ti_mnptmp(i,ni)/numworkers
              enddo
            enddo

            do i = 1,nz
              te_mnp(i)     = te_mnp(i) + te_mnptmp(i)/numworkers
            enddo

            ifintot1 = ifintot1 - 1

               endif

            enddo               ! end worker loop

! if we are here, we should have gathered up all the data

            if(ifintot .eq. 0) then
               ifintot = numworkers
               ntm = ntm + 1
!               call output ( hrut,ntm,istep,phi )
              do k = 1,nlt
                do j = 1,nf-1
                  do i = 1,nz
!             u1(i,j,k) =  eph(i,j,k)
!             u2(i,j,k) =  delph(i,j,k)
!             u3(i,j,k) =  phihp(j,k)
                     if(DoCoupleCimi)then
                        u4t(i,j,k) =   PotCimiOnSamiGrid_C(k,j)
                     else
                        u4t(i,j,k) =   0.0
                     endif
                 enddo
               enddo
             enddo
               
 
               call output ( hrut,ntm,istep,phi,denit,dennt,vsit,
     &                       sumvsit,tit,ut,vt,vpit,net,tet,tnt,u1t,
     &                     u2t,u3t,u4t,vnqt,vnpt,vnphit,jpt,jphit,
     &                     u1pt,u2st,u3ht,sigmapict,sigmahict,
     &                     sigmapt,sigmaht )
                         
! write the restart files and close those files

               open ( unit=210, file='time.rst', form='formatted' )
               open ( unit=211, file='deni.rst', form='unformatted' )
               open ( unit=212, file='vsi.rst', form='unformatted' )
               open ( unit=213, file='ti.rst', form='unformatted' )
               open ( unit=214, file='te.rst', form='unformatted' )
               open(2322,file='dphi.rst',form='unformatted')

               write(210,*) hrut
               write(211) deniout
               write(212) vsiout
               write(213) tiout
               write(214) teout
               write(2322) dphi

               close (210)
               close (211)
               close (212)
               close (213)
               close (214)
               close(2322)

!               do iwrk = 1,numworkers

!                  call mpi_send(ntm,1,MPI_INTEGER,iwrk,18,
!     .                 iComm,status,ierr)

!                  call mpi_send(ntmmax,1,MPI_INTEGER,iwrk,18,
!     .                 iComm,status,ierr)

!               enddo

!!!               if ( ntm .ge. ntmmax ) tflag = .false.



            endif

c Need to fix up dt calculation

            if(ifintot1 .eq. 0) then
               ifintot1 = numworkers
               ! get the DtCourant and dt and send
               DtCourant = minval(DtCourant_I)
               DtStop = DtAdvance - DeltaTime
               dt = min(DtCourant, DtStop)
               if (dt.eq.0) dt=DtCourant
               
               print *,'master sending dt',dt
               do  iwrk = 1,numworkers
                  call mpi_send(dt, 1, MPI_REAL, iwrk, 1, 
     .                 iComm, ierr)
                  call mpi_send(DtCourant, 1, MPI_REAL, iwrk, 1, 
     .                 iComm, ierr)
                  call mpi_send(deni_mnp,nz*nion,MPI_REAL,iwrk,1,
     .                 iComm, ierr)
                  call mpi_send(ti_mnp,nz*nion,MPI_REAL,iwrk,1,
     .                 iComm, ierr)
                  call mpi_send(te_mnp,nz,MPI_REAL,iwrk,1,
     .                 iComm, ierr)
!     print *,'master sending dt/mnp'
               enddo
            endif
            
            
         enddo                  ! end while (tflag)
      
         print *, 'MASTER: All Done with sami_run!' 


!         call mpi_abort(iComm, errorcode, ierr)
!         call mpi_finalize(ierr)

      
      endif

c******************* end master task ***************************************

c******************** worker task *******************************************

      if(taskid .gt. 0) then

! field line loop: actual run
         if (IsFirstCallWorker) then
            hrut    = hrinit
            timemax = hrmax * sphr
            istep   = 0
            tprnt   = 0.
            tneut   = 0.
            time    = 0.
            ntm     = 0
            ntmmax  = min(maxstep, int(hrmax / max(dt0/3600.,dthr)))
            IsFirstCallWorker = .false.
         endif
         ttflag  = .true.
!         print *,'iwrk ',taskid,ntmmax   
         
! initialize neutrals
! neutral density, temperature, and neutral wind
! already done in initialization

         if (restart) then
           do nll = 2,nl-1
              call neutambt (hrinit,nll) 
           enddo
         endif

!         do while (      istep .le. maxstep 
!     .          .and. time  .le. timemax+dt  )

         do while ( ttflag )

!            if (istep > maxstep) 
!     &           call con_stop(
!     &           'SAMI3 ERROR: stuck in time loop beyond maxstep')
!           print *,'start worker loop',istep,ttflag,taskid

! parallel transport

! Below is  nll = 2,nl-1 because of guard cells

            do nll = 2,nl-1
               do nfl = 1,nf
!                  call mag_zenith (hrut,nfl,nll)
                  call zenith (hrut,nfl,nll)
!                  call zenith (hrinit,nfl,nll)
                  call transprt (nfl,nll,p_crit)
               enddo         
            enddo

! Do data exchanges between guard cells

c buffer and send to the LEFT

            do k = 1,nion
               do j = 1,nf
                  do i = 1,nz
                     tl1s(i,j,k) = deni(i,j,2,k)
                  enddo
               enddo
            enddo
            do k = nion+1,nion+nion
               do j = 1,nf
                  do i = 1,nz
                     tl1s(i,j,k) = ti(i,j,2,k-nion)
                  enddo
               enddo
            enddo
            k = nion + nion + 1
            do j = 1,nf
               do i = 1,nz
                  tl1s(i,j,k) = te(i,j,2)
               enddo
            enddo

            call mpi_sendrecv(tl1s, (nion+nion+1)*nz*nf, MPI_REAL, 
     .           left, 0, tl1r, (nion+nion+1)*nz*nf, MPI_REAL, 
     .           right, 0, iComm, status, ierr)

c Now everybody receives

            do k = 1,nion
               do j = 1,nf
                  do i = 1,nz
                     deni(i,j,nl,k) = tl1r(i,j,k)
                  enddo
               enddo
            enddo
            do k = nion+1,nion+nion
               do j = 1,nf
                  do i = 1,nz
                     ti(i,j,nl,k-nion) = tl1r(i,j,k)
                  enddo
               enddo
            enddo
            k = nion + nion + 1
            do j = 1,nf
               do i = 1,nz
                  te(i,j,nl) = tl1r(i,j,k)
               enddo
            enddo

c Buffer and send to the RIGHT

            do k = 1,nion
               do j = 1,nf
                  do i = 1,nz
                     tr1s(i,j,k) = deni(i,j,nl-1,k)
                  enddo
               enddo
            enddo
            do k = nion+1,nion+nion
               do j = 1,nf
                  do i = 1,nz
                     tr1s(i,j,k) = ti(i,j,nl-1,k-nion)
                  enddo
               enddo
            enddo
            k = nion + nion + 1
            do j = 1,nf
               do i = 1,nz
                  tr1s(i,j,k) = te(i,j,nl-1)
               enddo
            enddo

            call mpi_sendrecv(tr1s, (nion+nion+1)*nz*nf, MPI_REAL, 
     .           right, 0, tr1r, (nion +nion +1)*nz*nf, MPI_REAL, 
     .           left, 0, iComm, status, ierr)

            do k = 1,nion
               do j = 1,nf
                  do i = 1,nz
                     deni(i,j,1,k) = tr1r(i,j,k)
                  enddo
               enddo
            enddo
            do k = nion+1,nion+nion
               do j = 1,nf
                  do i = 1,nz
                     ti(i,j,1,k-nion) = tr1r(i,j,k)
                  enddo
               enddo
            enddo
            k = nion + nion + 1
            do j = 1,nf
               do i = 1,nz
                  te(i,j,1) = tr1r(i,j,k)
               enddo
            enddo

! We are now finished exchanging guard cell data

! Sending hipcp and hidphig to master to calculate the
! potential

           call mpi_send(hipcp, nf*nl, MPI_REAL, 0, 2, 
     .          iComm, ierr)
           call mpi_send(hihcm, nf*nl, MPI_REAL, 0, 2, 
     .          iComm, ierr)
           call mpi_send(hipcphi, nf*nl, MPI_REAL, 0, 2, 
     .          iComm, ierr)
           call mpi_send(hidphig, nf*nl, MPI_REAL, 0, 2, 
     .          iComm, ierr)
           call mpi_send(hidpg, nf*nl, MPI_REAL, 0, 2, 
     .          iComm, ierr)
           call mpi_send(hidphiv, nf*nl, MPI_REAL, 0, 2, 
     .          iComm, ierr)
           call mpi_send(hidpv, nf*nl, MPI_REAL, 0, 2, 
     .          iComm, ierr)
           call mpi_send(hipc, nf*nl, MPI_REAL, 0, 2, 
     .          iComm, ierr)
           call mpi_send(hihc, nf*nl, MPI_REAL, 0, 2, 
     .          iComm, ierr)
           call mpi_send(hidv, nf*nl, MPI_REAL, 0, 2, 
     .          iComm, ierr)

! send master current time: hrut

           call mpi_send(hrut, 1, MPI_REAL, 0, 2, 
     .          iComm, ierr)


! now get the potential from master


        call mpi_recv(phi, nnx*nny, MPI_REAL, 0, 3, 
     .          iComm, status, ierr)
        call mpi_recv(phialt, nnx*nny, MPI_REAL, 0, 3, 
     .          iComm, status, ierr)
        call mpi_recv(philon, nnx*nny, MPI_REAL, 0, 3, 
     .          iComm, status, ierr)
        call mpi_recv(p_crit,nnx-1, MPI_REAL, 0, 3, 
     .          iComm, status, ierr)

!            print *,taskid,' just got phi'

! perpendicular transport

!        if (DoCoupleCimi) then
!           call exb(hrut,phi+PotCimiOnSamiGrid_C,phialt,philon) 
!        else
           call exb(hrut,phi,phialt,philon) 
!        endif
        call neut(hrut)
        call courant 

!       average magnetic pole grid values (deni,Ti,Te)

        j0           = nf

        do ni = nion1,nion2
          do i = 1,nz
            deni_mnp0 = 0.
            ti_mnp0   = 0.
            do k = 2,nl-1
              if ( alts (i,j0,k) .lt. alt_crit_avg) then
                deni_mnp0      = deni_mnp0 + deni(i,j0,k,ni)
                ti_mnp0        = ti_mnp0 + ti(i,j0,k,ni)
              endif
            enddo
            deni_mnp(i,ni) = deni_mnp0 / float(nl-2)
            ti_mnp(i,ni)   = ti_mnp0 / float(nl-2)
          enddo
        enddo

        do i = 1,nz
          te_mnp0 = 0.
          do k = 2,nl-1
            if ( alts (i,j0,k) .lt. alt_crit_avg) then
              te_mnp0     = te_mnp0 + te(i,j0,k)
            endif
          enddo
          te_mnp(i)   = te_mnp0 / float(nl-2)
        enddo


c send local dt

!         print *,'send local dt',taskid

        call mpi_send(dt, 1, MPI_REAL, 0, 1, 
     .          iComm, ierr)
        call mpi_send(DeltaTime, 1, MPI_REAL, 0, 1, 
     .          iComm, ierr)
        call mpi_send(time, 1, MPI_REAL, 0, 1, 
     .          iComm, ierr)
        call mpi_send(istep, 1, MPI_INTEGER, 0, 1, 
     .          iComm, ierr)

        call mpi_send(deni_mnp,nz*nion,MPI_REAL,0,1,
     .          iComm, ierr)
        call mpi_send(ti_mnp,nz*nion,MPI_REAL,0,1,
     .          iComm, ierr)
        call mpi_send(te_mnp,nz,MPI_REAL,0,1,
     .          iComm, ierr)

c get global dt

!         print *,'now receive dt',taskid

        call mpi_recv(dt, 1, MPI_REAL, 0, 1, 
     .          iComm, status, ierr)
        call mpi_recv(DtCourant, 1, MPI_REAL, 0, 1, 
     .          iComm, status, ierr)

!        call mpi_recv(flux,linesuv , MPI_REAL, 0, 1, 
!     .          iComm, status, ierr)

        call mpi_recv(deni_mnp,nz*nion,MPI_REAL,0,1,
     .          iComm, status, ierr)
        call mpi_recv(ti_mnp,nz*nion,MPI_REAL,0,1,
     .          iComm, status, ierr)
        call mpi_recv(te_mnp,nz,MPI_REAL,0,1,
     .          iComm, status, ierr)

!         print *,'done receiving dt/mnp',taskid,dt,hrut
     
!         do ni = nion1,nion2
!           do i = 1,nz
!             print *,'taskid ...',deni_mnp(i,ni),deni(i,nf,0,ni)
!           enddo
!         enddo

! update neutrals

        if( tneut .ge. 0.25 ) then
!         print *,'calling neutambt',taskid
           do nll = 2,nl-1
              call neutambt (hrut,nll) 
!              call neutambt (hrinit,nll) 
           enddo
           tneut = 0.
           if ( hrut .lt. hrpr+hrinit ) then
             print *,'No output yet: hr = ',hrut
           endif
        endif

! time/step advancement

        istep  = istep + 1
        time   = time  + dt
        hrut   = time / sphr + hrinit
        DeltaTime = DeltaTime+dt
        tneut  = tneut + dt / sphr
!            print *,'time/dt/hrut',time,dt,hrut,taskid


! output data
        ! check if it is time to print
        lprnt = floor((time+1.0e-5)/(dthr*3600.)) /= 
     &       floor((time+1.0e-5-dt)/(dthr*3600.))
        

        if ( lprnt .and. hrut .ge. hrpr+hrinit) then

!         print *,'sending output to master',taskid
! We no longer call output from here, but send data to the MASTER
! The four things we want to send are  deni, ti, vsi, te

           ntm    = ntm + 1

           call mpi_send(time, 1, MPI_REAL, 0, 0, 
     .          iComm, ierr)
           call mpi_send(hrut, 1, MPI_REAL, 0, 0, 
     .          iComm, ierr)
           call mpi_send(istep, 1, MPI_INTEGER, 0, 0, 
     .          iComm, ierr)
           do nntmp = 1,nion
              do itmp = 1,nz
                 do jtmp = 1,nf
                    do ktmp = 1,nl
                       denitmp(itmp,jtmp,ktmp) 
     .                      = deni(itmp,jtmp,ktmp,nntmp)
                       denntmp(itmp,jtmp,ktmp) 
     .                      = denn(itmp,jtmp,ktmp,nntmp)
                       titmp(itmp,jtmp,ktmp) 
     .                      = ti(itmp,jtmp,ktmp,nntmp)
                       vsitmp(itmp,jtmp,ktmp) 
     .                      = vsi(itmp,jtmp,ktmp,nntmp)
                    enddo
                 enddo
              enddo
              call mpi_send(denitmp, nz*nf*nl, MPI_REAL, 0, 0, 
     .             iComm, ierr)
              call mpi_send(denntmp, nz*nf*nl, MPI_REAL, 0, 0, 
     .             iComm, ierr)
              call mpi_send(titmp, nz*nf*nl, MPI_REAL, 0, 0, 
     .             iComm, ierr)
              call mpi_send(vsitmp, nz*nf*nl, MPI_REAL, 0, 0, 
     .             iComm, ierr)
           enddo
           call mpi_send(te, nz*nf*nl, MPI_REAL, 0, 0, 
     .          iComm, ierr)
           call mpi_send(u1p, nz*nf*nl, MPI_REAL, 0, 0, 
     .          iComm,  ierr)
           call mpi_send(u2s, nz*nf*nl, MPI_REAL, 0, 0, 
     .          iComm,  ierr)
           call mpi_send(u3h, nz*nf*nl, MPI_REAL, 0, 0, 
     .          iComm,  ierr)
           call mpi_send(u1, nz*nf*nl, MPI_REAL, 0, 0, 
     .          iComm,  ierr)
           call mpi_send(u2, nz*nf*nl, MPI_REAL, 0, 0, 
     .          iComm,  ierr)
           call mpi_send(u3, nz*nf*nl, MPI_REAL, 0, 0, 
     .          iComm,  ierr)
           call mpi_send(u4, nz*nf*nl, MPI_REAL, 0, 0, 
     .          iComm,  ierr)
           call mpi_send(sigmap, nz*nf*nl, MPI_REAL, 0, 0, 
     .          iComm, ierr)
           call mpi_send(sigmah, nz*nf*nl, MPI_REAL, 0, 0, 
     .          iComm, ierr)
           call mpi_send(sigmapic, nz*nf*nl, MPI_REAL, 0, 0, 
     .          iComm, ierr)
           call mpi_send(sigmahic, nz*nf*nl, MPI_REAL, 0, 0, 
     .          iComm, ierr)
!           call mpi_send(hipcp, nf*nl, MPI_REAL, 0, 0, 
!     .          iComm, ierr)
!           call mpi_send(hipcphi, nf*nl, MPI_REAL, 0, 0, 
!     .          iComm, ierr)
!           call mpi_send(hihcm, nf*nl, MPI_REAL, 0, 0, 
!     .          iComm, ierr)
!           call mpi_send(hidpv, nf*nl, MPI_REAL, 0, 0, 
!     .          iComm, ierr)
!           call mpi_send(hidphiv, nf*nl, MPI_REAL, 0, 0, 
!     .          iComm, ierr)
!           call mpi_send(hidpg, nf*nl, MPI_REAL, 0, 0, 
!     .          iComm, ierr)
!           call mpi_send(hidphig, nf*nl, MPI_REAL, 0, 0, 
!     .          iComm, ierr)
           call mpi_send(vnq, nz*nf*nl, MPI_REAL, 0, 0, 
     .          iComm, ierr)
           call mpi_send(vnp, nz*nf*nl, MPI_REAL, 0, 0, 
     .          iComm, ierr)
           call mpi_send(vnphi, nz*nf*nl, MPI_REAL, 0, 0, 
     .          iComm, ierr)
           call mpi_send(jp, nz*nf*nl, MPI_REAL, 0, 0, 
     .          iComm, ierr)
           call mpi_send(jphi, nz*nf*nl, MPI_REAL, 0, 0, 
     .          iComm, ierr)
!           call mpi_send(hipc, nf*nl, MPI_REAL, 0, 0, 
!     .          iComm, ierr)
!           call mpi_send(hihc, nf*nl, MPI_REAL, 0, 0, 
!     .          iComm, ierr)
!           call mpi_send(hidv, nf*nl, MPI_REAL, 0, 0, 
!     .          iComm, ierr)
 
          tprnt   = 0.

!                  call mpi_recv(ntm,1,MPI_INTEGER,0,18,
!     .                 MPI_COM_WORLD,status,ierr)

!                  call mpi_recv(ntmmax,1,MPI_INTEGER,0,18,
!     .                 iComm,status,ierr)
          
!             print *,'in worker',ntm,ntmmax,taskid

!           if ( ntm .ge. ntmmax ) ttflag = .false.

          lprnt   = .false.

       endif
       
       if ( DtAdvance-DeltaTime < StopTolerance ) ttflag = .false.
 
            
!!! End of sami time advance unit
        enddo    ! end time loop

      endif ! end worker task
!!! Finalizing
      xxx = 1.
      call mpi_send(xxx, 1, MPI_REAL, 0, 10, 
     &    iComm, ierr)

      call MPI_BARRIER(iComm,ierr)

!      call mpi_finalize(ierr)
!      print *,'done finalizing,taskid',taskid
!      print *,'Finished sami_run for task',taskid
      
      return
      end

!==============================================================================
      subroutine sami_finalize
! close files
      
      close (10)
      close (20)
      close (40)
      close (70)
      close (71)
      close (72)
      close (73)
      close (74)
      close (75)
      close (78)
      close (79)
      close (80)
      close (90)
      close (91)
      close (92)
      close (93)
      close (94)
      close (95)
      close (96)
      close (97)
      close (98)
      close (81)
      close (82)
      close (83)
      close (84)
      close (85)
      close (86)
      close (87)
      close (88)
      close (711)
      close (712)
      close (713)
      close (714)
      close (715)
      close (1712)
      close (1713)
      close (1714)
      close (1715)
      close (569)
      close (716)
      close (717)
      close (1718)
      close (811)
      close (812)
      close (813)
      close (814)
      close (815)
      close (816)
      close (817)
      close (911)
      close (912)
      close (913)
      close (914)
      close (915)
      close (916)
      close (917)
      close (384)
      close (385)
      close (386)
      
      close (196)
      close (197)
      close (198)
      
      close (201)
      close (202)
      
      close (491)
      close (492)
      close (493)
      close (494)
      close (495)
      close (496)
      close (497)
      close (498)
      

      end

*******************************************
*******************************************

!            initial

*******************************************
*******************************************

      subroutine initial(phialt,philon,p_crit)
      use ModMpi
      UseModSAMI, ONLY: iComm

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc' 

      real, dimension(:,:,:), allocatable :: altstmp,glatstmp,glonstmp
      real, dimension(:,:,:), allocatable :: altst,glatst,glonst
      real, dimension(:,:,:), allocatable :: baltst,blatst,blonst
      real, dimension(:,:,:), allocatable :: xst,yst,zst
      real, dimension(:,:,:), allocatable :: altptmp,blatptmp,blonptmp
!      real, dimension(:,:,:), allocatable :: xalttmp,xblattmp,xblontmp
      real, dimension(:,:,:), allocatable :: baltpt
!      real, dimension(:,:,:), allocatable :: baltpt,blonpt
!      real, dimension(:,:,:), allocatable :: baltpt,blatpt,blonpt
      real, dimension(:,:,:), allocatable :: vpsnxt,vpsnyt,vpsnzt
      real, dimension(:,:,:), allocatable :: vhsnxt,vhsnyt,vhsnzt
      real, dimension(:,:,:), allocatable :: xpt,ypt,zpt
      real, dimension(:,:,:), allocatable :: bdirsxt,bdirsyt,bdirszt
      real, dimension(:,:,:), allocatable :: gsthetaxt
      real, dimension(:,:,:), allocatable :: gsthetayt
      real, dimension(:,:,:), allocatable :: gsthetazt
      real, dimension(:,:,:), allocatable :: gsphixt
      real, dimension(:,:,:), allocatable :: gsphiyt
      real, dimension(:,:,:), allocatable :: gsphizt
      real, dimension(:,:,:), allocatable :: gsrxt
      real, dimension(:,:,:), allocatable :: gsryt
      real, dimension(:,:,:), allocatable :: gsrzt
      real, dimension(:,:,:), allocatable :: xrgt
      real, dimension(:,:,:), allocatable :: xthgt
      real, dimension(:,:,:), allocatable :: xphigt
      character(len=14),parameter :: plotdir='IM/plotsSAMI3/'
!      real fism(linesuv)

      real f1026(nz,nf,nl,91),f584(nz,nf,nl,91),
     .     f304 (nz,nf,nl,91),f1216(nz,nf,nl,91)

!      include 'mpif.h'
      integer status(MPI_STATUS_SIZE)

      real zi(29),denii(29,7)
      real phionr(linesuv,5)

      real fluxdat(linesuv,2) ! original euvac stuff

c Some local variables

!      real altstmp(nz,nf,nl), glatstmp(nz,nf,nl), glonstmp(nz,nf,nl)
!      real altst(nz,nf,nlt), glatst(nz,nf,nlt), glonst(nz,nf,nlt)
!      real baltst(nz,nf,nlt), blatst(nz,nf,nlt), blonst(nz,nf,nlt)
!      real xst(nz,nf,nlt), yst(nz,nf,nlt), zst(nz,nf,nlt)
!      real altptmp(nzp1,nfp1,nlp1), blatptmp(nzp1,nfp1,nlp1), 
!     .                              blonptmp(nzp1,nfp1,nlp1)
!      real baltpt(nzp1,nfp1,nlt), blatpt(nzp1,nfp1,nlt), 
!     .                            blonpt(nzp1,nfp1,nlt)

      real phialt(nnx,nny),philon(nnx,nny),p_crit(nnx-1)

      namelist / go / fmtout,maxstep,hrmax,dthr,hrpr,dt0,
     .                grad_in,glat_in,glon_in,
     .                fejer,
     .                rmin,rmax,
     .                altmin,
     .                fbar,f10p7,ap,
     .                year,day,mmass,
     .                nion1,nion2,hrinit,tvn0,tvexb0,ver,veh,vw,
     .                 gams1,gams1m,gamp1,nz1,
     .                 gams2,gams2m,gamp2,nz2,
     .                 gams3,gams3m,gamp3,nz3,
     .                 gams4,gams4m,gamp4,nz4,
     .                 r_min1,r_max1,
     .                 r_max2,alt_crit_avg,
     .                 blat_max3,blat_max4,
     .                snn,stn,denmin,alt_crit,cqe,plat,plon,
     .                dellon,psmooth,hall,restart,
     .                storm_ti,storm_tf,vexb_max,
     .                lmadala,lcr,lvs,lweimer,decay_time,pcrit,
     .                lhwm93,lhwm14,
     .                 vsi0,delta_vsi0,anu_drag0

      if ( taskid .eq. 0 ) then
        allocate (altstmp(nz,nf,nl), 
     .            glatstmp(nz,nf,nl), 
     .            glonstmp(nz,nf,nl))
        allocate (altst(nz,nf,nlt), 
     .            glatst(nz,nf,nlt), 
     .            glonst(nz,nf,nlt))
        allocate (baltst(nz,nf,nlt), 
     .            blatst(nz,nf,nlt), 
     .            blonst(nz,nf,nlt))
        allocate (xst(nz,nf,nlt), 
     .            yst(nz,nf,nlt), 
     .            zst(nz,nf,nlt))
        allocate (altptmp(nzp1,nfp1,nlp1), 
     .            blatptmp(nzp1,nfp1,nlp1), 
     .            blonptmp(nzp1,nfp1,nlp1))
!        allocate (xalttmp(nzp1,nf,nl), 
!     .            xblattmp(nzp1,nf,nl), 
!     .            xblontmp(nzp1,nf,nl))
        allocate  (baltpt(nzp1,nfp1,nlt))
!     .            blatpt(nzp1,nfp1,nlt), 
!     .            blonpt(nzp1,nfp1,nlt))
        allocate (vpsnxt(nz,nf,nlt), 
     .            vpsnyt(nz,nf,nlt), 
     .            vpsnzt(nz,nf,nlt))
        allocate (vhsnxt(nz,nf,nlt), 
     .            vhsnyt(nz,nf,nlt), 
     .            vhsnzt(nz,nf,nlt))
        allocate (xpt(nzp1,nfp1,nlt), 
     .            ypt(nzp1,nfp1,nlt), 
     .            zpt(nzp1,nfp1,nlt))  
        allocate (bdirsxt(nz,nf,nlt), 
     .            bdirsyt(nz,nf,nlt), 
     .            bdirszt(nz,nf,nlt))
        allocate (gsthetaxt(nz,nf,nlt), 
     .            gsthetayt(nz,nf,nlt), 
     .            gsthetazt(nz,nf,nlt))
        allocate (gsphixt(nz,nf,nlt), 
     .            gsphiyt(nz,nf,nlt), 
     .            gsphizt(nz,nf,nlt))
        allocate (gsrxt(nz,nf,nlt), 
     .            gsryt(nz,nf,nlt), 
     .            gsrzt(nz,nf,nlt))
        allocate (xrgt(nz,nf,nlt), 
     .            xthgt(nz,nf,nlt), 
     .            xphigt(nz,nf,nlt))
      endif

! read in parameters and initial ion density data 
         print *,'read namelist',taskid
      if(taskid .eq. 0) then
         read(10,go)
      endif

c send the namelist data to all the other processors

      call mpi_bcast(fmtout,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(maxstep,1,MPI_INTEGER,0,iComm,ierr)
      call mpi_bcast(hrmax,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(dthr,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(hrpr,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(dt0,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(grad_in,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(glat_in,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(glon_in,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(fejer,1,MPI_LOGICAL,0,iComm,ierr)
      call mpi_bcast(rmin,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(rmax,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(altmin,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(fbar,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(f10p7,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(ap,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(year,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(day,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(mmass,1,MPI_INTEGER,0,iComm,ierr)
      call mpi_bcast(nion1,1,MPI_INTEGER,0,iComm,ierr)
      call mpi_bcast(nion2,1,MPI_INTEGER,0,iComm,ierr)
      call mpi_bcast(tvn0,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(tvexb0,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(ver,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(veh,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(vw,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(gams1,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(gams1m,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(gamp1,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(nz1,1,MPI_INTEGER,0,iComm,ierr)
      call mpi_bcast(gams2,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(gams2m,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(gamp2,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(nz2,1,MPI_INTEGER,0,iComm,ierr)
      call mpi_bcast(gams3,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(gams3m,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(gamp3,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(nz3,1,MPI_INTEGER,0,iComm,ierr)
      call mpi_bcast(gams4,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(gams4m,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(gamp4,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(nz4,1,MPI_INTEGER,0,iComm,ierr)
      call mpi_bcast(r_min1,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(r_max1,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(r_max2,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(alt_crit_avg,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(blat_max3,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(blat_max4,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(snn,7,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(stn,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(denmin,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(alt_crit,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(cqe,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(plat,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(plon,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(dellon,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(psmooth,1,MPI_INTEGER,0,iComm,ierr)
      call mpi_bcast(hall,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(storm_ti,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(storm_tf,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(vexb_max,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(restart,1,MPI_LOGICAL,0,iComm,ierr)
      call mpi_bcast(lmadala,1,MPI_LOGICAL,0,iComm,ierr)
      call mpi_bcast(lcr,1,MPI_LOGICAL,0,iComm,ierr)
      call mpi_bcast(lvs,1,MPI_LOGICAL,0,iComm,ierr)
      call mpi_bcast(lweimer,1,MPI_LOGICAL,0,iComm,ierr)
      call mpi_bcast(decay_time,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(pcrit,1,MPI_REAL,0,iComm,ierr)
      if (.not. restart)
     .  call mpi_bcast(hrinit,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(lhwm93,1,MPI_LOGICAL,0,iComm,ierr)
      call mpi_bcast(lhwm14,1,MPI_LOGICAL,0,iComm,ierr)
      call mpi_bcast(vsi0,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(delta_vsi0,1,MPI_REAL,0,iComm,ierr)
      call mpi_bcast(anu_drag0,1,MPI_REAL,0,iComm,ierr)
!            print *,'done with bcast', taskid
      dt = dt0

      ami(pthp)  = 1.
      ami(pthep) = 4.
      ami(ptnp)  = 14.
      ami(ptop)  = 16.
      ami(ptn2p) = 28.
      ami(ptnop) = 30.
      ami(pto2p) = 32.

      amn(pth)  = 1.
      amn(pthe) = 4.
      amn(ptn)  = 14.
      amn(pto)  = 16.
      amn(ptn2) = 28.
      amn(ptno) = 30.
      amn(pto2) = 32.

      alpha0(pth)  = 0.67
      alpha0(pthe) = 0.21
      alpha0(ptn)  = 1.10
      alpha0(pto)  = 0.79
      alpha0(ptn2) = 1.76
      alpha0(ptno) = 1.74
      alpha0(pto2) = 1.59

      do i = 1,7
        aap(i) = ap
      enddo

! read in initial density data

      if(taskid .eq. 0) then
         do i = 1,29
            read(20,102) zi(i),(denii(i,j),j=1,7)
 102        format(1x,f7.1,1p7e8.1)
         enddo
      endif
!       print *,'bcast zi', taskid
      call mpi_bcast(zi,29,
     &       MPI_REAL,0,iComm,ierr)
      call mpi_bcast(denii,29*7,
     &       MPI_REAL,0,iComm,ierr)

! read in chemistry data
! in format statement 104 need to 'hardwire' nneut (= 7)

      if(taskid .eq. 0) then
         do k = 1,nchem
            read(30,103) (ichem(k,j),j=1,3)
 103        format(3i3)
         enddo
      endif
!       print *,'bcast ichem', taskid
      call mpi_bcast(ichem,nchem*3,
     &       MPI_INTEGER,0,iComm,ierr)
!       print *,'done bcast ichem', taskid
! generate the mesh data by everybody but the Master
!      write(*,*) 'start with blonp0a_mpi',taskid
      if (taskid .eq. 0) call blonp0a
!      write(*,*) 'done with blonp0a_mpi',taskid
      if(taskid .gt. 0) then 
!         write(*,*) 'start with grid3_mpi',taskid
         call grid3_mpi
!         write(*,*) 'done with grid3_mpi',taskid
      endif
C Now wait to receive back the results from each worker task

      if(taskid .eq. 0) then
         
         ifintot = numworkers

         do while( ifintot .gt. 0)

         do  iwrk = 1, numworkers
            source = iwrk
            dest = source
            
            call mpi_iprobe(source, 0, 
     .           iComm, flagit, status, ierr)
               
            if(flagit .and. ifintot .gt. 0) then

!  The three things we want to receive are  altpt blatpt blonpt

               call mpi_recv(altptmp, nzp1*nfp1*nlp1, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(blatptmp, nzp1*nfp1*nlp1, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(blonptmp, nzp1*nfp1*nlp1, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)

!  Put the submatrices into the correct matrix

               do k = 2,nl-1
                  kk = (iwrk-1)*(nl -2) +  k - 1
                  if(kk .eq. 0) kk = nlt
                  if(kk .eq. nltp1) kk = 1
                  do j = 1,nfp1
                     do i = 1,nzp1
                        baltpt(i,j,kk) = altptmp(i,j,k)
                        blatpt(i,j,kk) = blatptmp(i,j,k)
                        blonpt(i,j,kk) = blonptmp(i,j,k)
                     enddo
                  enddo
               enddo

!  The three things we want to receive are  xpt ypt zpt

               call mpi_recv(altptmp, nzp1*nfp1*nlp1, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(blatptmp, nzp1*nfp1*nlp1, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(blonptmp, nzp1*nfp1*nlp1, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)

!  Put the submatrices into the correct matrix

               do k = 2,nl-1
                  kk = (iwrk-1)*(nl -2) + k - 1
                  if(kk .eq. 0) kk = nlt
                  if(kk .eq. nltp1) kk = 1
                  do j = 1,nfp1
                     do i = 1,nzp1
                        xpt(i,j,kk) = altptmp(i,j,k)
                        ypt(i,j,kk) = blatptmp(i,j,k)
                        zpt(i,j,kk) = blonptmp(i,j,k)
                     enddo
                  enddo
               enddo

! We want to receive pp

               call mpi_recv(altptmp, nzp1*nfp1*nlp1, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)

! Put the submatrices into the correct matrix

               do k = 2,nl-1
                  kk = (iwrk-1)*(nl -2) +k -1
                  if(kk .eq. 0) kk = nlt
                  if(kk .eq. nltp1) kk = 1
                  do j = 1,nfp1
                     do i = 1,nzp1
                        ppt(i,j,kk)   = altptmp(i,j,k)
                     enddo
                  enddo
               enddo

! The three things we want to receive are  altst glatst glonst

               call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)

! Put the submatrices into the correct matrix

               do k = 2,nl-1
                  kk = (iwrk-1)*(nl -2) +k -1
                  if(kk .eq. 0) kk = nlt
                  if(kk .eq. nltp1) kk = 1
                  do j = 1,nf
                     do i = 1,nz
                        altst(i,j,kk)  = altstmp(i,j,k)
                        glatst(i,j,kk) = glatstmp(i,j,k)
                        glonst(i,j,kk) = glonstmp(i,j,k)
                     enddo
                  enddo
               enddo

! The three things we want to receive are  baltst blatst blonst

               call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)

! Put the submatrices into the correct matrix

               do k = 2,nl-1
                  kk = (iwrk-1)*(nl -2) + k - 1
                  if(kk .eq. 0) kk = nlt
                  if(kk .eq. nltp1) kk = 1
                  do j = 1,nf
                     do i = 1,nz
                        baltst(i,j,kk) = altstmp(i,j,k)
                        blatst(i,j,kk) = glatstmp(i,j,k)
                        blonst(i,j,kk) = glonstmp(i,j,k)
                     enddo
                  enddo
               enddo

! The three things we want to receive are  xst yst zst

               call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)

! Put the submatrices into the correct matrix

               do k = 2,nl-1
                  kk = (iwrk-1)*(nl -2) +k -1
                  if(kk .eq. 0) kk = nlt
                  if(kk .eq. nltp1) kk = 1
                  do j = 1,nf
                     do i = 1,nz
                        xst(i,j,kk) = altstmp(i,j,k)
                        yst(i,j,kk) = glatstmp(i,j,k)
                        zst(i,j,kk) = glonstmp(i,j,k)
                     enddo
                  enddo
               enddo

! The three things we want to receive are  xrg,xthg,xphig

               call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)

! Put the submatrices into the correct matrix

               do k = 2,nl-1
                  kk = (iwrk-1)*(nl -2) +k -1
                  if(kk .eq. 0) kk = nlt
                  if(kk .eq. nltp1) kk = 1
                  do j = 1,nf
                     do i = 1,nz
                        xrgt(i,j,kk)   = altstmp(i,j,k)
                        xthgt(i,j,kk)  = glatstmp(i,j,k)
                        xphigt(i,j,kk) = glonstmp(i,j,k)
                     enddo
                  enddo
               enddo

!  The three things we want to receive are vpsnx vpsny vpsnz

               call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)

!  Put the submatrices into the correct matrix

               do k = 2,nl-1
                  kk = (iwrk-1)*(nl -2) + k - 1
                  if(kk .eq. 0) kk = nlt
                  if(kk .eq. nltp1) kk = 1
                  do j = 1,nf
                     do i = 1,nz
                        vpsnxt(i,j,kk) = altstmp(i,j,k)
                        vpsnyt(i,j,kk) = glatstmp(i,j,k)
                        vpsnzt(i,j,kk) = glonstmp(i,j,k)
                     enddo
                  enddo
               enddo

!  The three things we want to receive are vhsnx vhsny vhsnz

               call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)

!  Put the submatrices into the correct matrix

               do k = 2,nl-1
                  kk = (iwrk-1)*(nl -2) + k - 1
                  if(kk .eq. 0) kk = nlt
                  if(kk .eq. nltp1) kk = 1
                  do j = 1,nf
                     do i = 1,nz
                        vhsnxt(i,j,kk) = altstmp(i,j,k)
                        vhsnyt(i,j,kk) = glatstmp(i,j,k)
                        vhsnzt(i,j,kk) = glonstmp(i,j,k)
                     enddo
                  enddo
               enddo

!  The three things we want to receive are bdirsx bdirsy bdirsz

               call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)

!  Put the submatrices into the correct matrix

               do k = 2,nl-1
                  kk = (iwrk-1)*(nl -2) + k - 1
                  if(kk .eq. 0) kk = nlt
                  if(kk .eq. nltp1) kk = 1
                  do j = 1,nf
                     do i = 1,nz
                        bdirsxt(i,j,kk) = altstmp(i,j,k)
                        bdirsyt(i,j,kk) = glatstmp(i,j,k)
                        bdirszt(i,j,kk) = glonstmp(i,j,k)
                     enddo
                  enddo
               enddo

!               print *,'main',bdirsxt(1,nf,1),bdirsxt(1,nf-1,1)

!  The three things we want to receive are gsthetax/y/z

               call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)

!  Put the submatrices into the correct matrix

               do k = 2,nl-1
                  kk = (iwrk-1)*(nl -2) + k - 1
                  if(kk .eq. 0) kk = nlt
                  if(kk .eq. nltp1) kk = 1
                  do j = 1,nf
                     do i = 1,nz
                        gsthetaxt(i,j,kk) = altstmp(i,j,k)
                        gsthetayt(i,j,kk) = glatstmp(i,j,k)
                        gsthetazt(i,j,kk) = glonstmp(i,j,k)
                     enddo
                  enddo
               enddo

!  The three things we want to receive are gsphix/y/z

               call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)

!  Put the submatrices into the correct matrix

               do k = 2,nl-1
                  kk = (iwrk-1)*(nl -2) + k - 1
                  if(kk .eq. 0) kk = nlt
                  if(kk .eq. nltp1) kk = 1
                  do j = 1,nf
                     do i = 1,nz
                        gsphixt(i,j,kk) = altstmp(i,j,k)
                        gsphiyt(i,j,kk) = glatstmp(i,j,k)
                        gsphizt(i,j,kk) = glonstmp(i,j,k)
                     enddo
                  enddo
               enddo


!  The three things we want to receive are gsrx/y/z

               call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)
               call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, 
     .              iwrk, 0, iComm, status, ierr)

!  Put the submatrices into the correct matrix

               do k = 2,nl-1
                  kk = (iwrk-1)*(nl -2) + k - 1
                  if(kk .eq. 0) kk = nlt
                  if(kk .eq. nltp1) kk = 1
                  do j = 1,nf
                     do i = 1,nz
                        gsrxt(i,j,kk) = altstmp(i,j,k)
                        gsryt(i,j,kk) = glatstmp(i,j,k)
                        gsrzt(i,j,kk) = glonstmp(i,j,k)
                     enddo
                  enddo
               enddo


                  ifintot = ifintot -1
                  
            endif

           enddo                  ! end worker loop

         enddo

! if we are here, we should have gathered up all the data

! output grid data

         if ( fmtout ) then
            close(69)
            close(76)
            close(77)
            open ( unit=69, file=plotdir//'zaltf.dat',form='formatted' )
            open ( unit=76, file=plotdir//'glatf.dat',form='formatted' )
            open ( unit=77, file=plotdir//'glonf.dat',form='formatted' )
            write(69,100) altst
            write(76,100) glatst
            write(77,100) glonst
            close(69)
            close(76)
            close(77)
         else
            open ( unit=69, file=plotdir//'zaltu.dat',form='unformatted' )
            open ( unit=76, file=plotdir//'glatu.dat',form='unformatted' )
            open ( unit=77, file=plotdir//'glonu.dat',form='unformatted' )
            write(69) altst
            write(76) glatst
            write(77) glonst
            close(69)
            close(76)
            close(77)

            open (144,file=plotdir//'glons_cg.rst',form='unformatted')
            write(144) glonst
            close(144)


            open ( unit=69, file=plotdir//'baltu.dat',form='unformatted' )
            open ( unit=76, file=plotdir//'blatu.dat',form='unformatted' )
            open ( unit=77, file=plotdir//'blonu.dat',form='unformatted' )
            write(69) baltst
            write(76) blatst
            write(77) blonst
            close(69)
            close(76)
            close(77)

            open (144,file=plotdir//'blons_cg.rst',form='unformatted')
            write(144) blonst
            close(144)

            open ( unit=69, file=plotdir//'xsu.dat',form='unformatted' )
            open ( unit=76, file=plotdir//'ysu.dat',form='unformatted' )
            open ( unit=77, file=plotdir//'zsu.dat',form='unformatted' )
            write(69) xst
            write(76) yst
            write(77) zst
            close(69)
            close(76)
            close(77)

            open(unit=169,file=plotdir//'baltpu.dat',form='unformatted')
            open(unit=176,file=plotdir//'blatpu.dat',form='unformatted')
            open(unit=177,file=plotdir//'blonpu.dat',form='unformatted')
            write(169) baltpt
            write(176) blatpt
            write(177) blonpt
            close(169)
            close(176)
            close(177)

            open (144,file=plotdir//'blonp_cg.rst',form='unformatted')
            write(144) blonpt
            close(144)

!  define phialt and philon (alt and lon of phi) and p_crit

               do j = 1,nny
                 do i = 1,nnx-1
                   jj = j + 1
                   phialt(i,j) = baltpt(nzp1/2,jj,i) + re
                   philon(i,j) = blonpt(nzp1/2,jj,i)
                   philon(i,j) = blonp0t(i+1)
                 enddo
               enddo

               do j = 1,nny
                   i = nnx
                   jj = j + 1
                   phialt(i,j) = baltpt(nzp1/2,jj,i-1) + re
                   phialt(i,j) = baltpt(nzp1/2,jj,i-1) + re +
     .                    (baltpt(nzp1/2,jj,i-1)-baltpt(nzp1/2,jj,i-2))
                   philon(i,j) = 360.
               enddo

               do i = 1,nnx-1
                 p_crit(i) = 6.
               enddo

!               do i=1,nnx
!                 print *,i,philon(i,25)
!               enddo

            open(unit=169,file=plotdir//'xpu.dat' ,form='unformatted' )
            open(unit=176,file=plotdir//'ypu.dat' ,form='unformatted' )
            open(unit=177,file=plotdir//'zpu.dat' ,form='unformatted' )
            write(169) xpt
            write(176) ypt
            write(177) zpt
            close(169)
            close(176)
            close(177)

            open(unit=169,file=plotdir//'xrgu.dat' ,form='unformatted' )
            open(unit=176,file=plotdir//'xthgu.dat',form='unformatted' )
            open(unit=177,file=plotdir//'xphigu.dat',form='unformatted')
            write(169) xrgt
            write(176) xthgt
            write(177) xphigt
            close(169)
            close(176)
            close(177)

            open(unit=169,file=plotdir//'vpsnxu.dat',form='unformatted')
            open(unit=176,file=plotdir//'vpsnyu.dat',form='unformatted')
            open(unit=177,file=plotdir//'vpsnzu.dat',form='unformatted')
            write(169) vpsnxt
            write(176) vpsnyt
            write(177) vpsnzt
            close(169)
            close(176)
            close(177)

            open(unit=169,file=plotdir//'vhsnxu.dat',form='unformatted')
            open(unit=176,file=plotdir//'vhsnyu.dat',form='unformatted')
            open(unit=177,file=plotdir//'vhsnzu.dat',form='unformatted')
            write(169) vhsnxt
            write(176) vhsnyt
            write(177) vhsnzt
            close(169)
            close(176)
            close(177)

            open(unit=169,file=plotdir//'bdirsxu.dat',form='unformatted')
            open(unit=176,file=plotdir//'bdirsyu.dat',form='unformatted')
            open(unit=177,file=plotdir//'bdirszu.dat',form='unformatted')
            write(169) bdirsxt
            write(176) bdirsyt
            write(177) bdirszt
            close(169)
            close(176)
            close(177)

         open(unit=169,file=plotdir//'gsthetaxu.dat',form='unformatted')
         open(unit=176,file=plotdir//'gsthetayu.dat',form='unformatted')
         open(unit=177,file=plotdir//'gsthetazu.dat',form='unformatted')
            write(169) gsthetaxt
            write(176) gsthetayt
            write(177) gsthetazt
            close(169)
            close(176)
            close(177)

            open(unit=169,file=plotdir//'gsphixu.dat',form='unformatted')
            open(unit=176,file=plotdir//'gsphiyu.dat',form='unformatted')
            open(unit=177,file=plotdir//'gsphizu.dat',form='unformatted')
            write(169) gsphixt
            write(176) gsphiyt
            write(177) gsphizt
            close(169)
            close(176)
            close(177)


            open(unit=169,file=plotdir//'gsrxu.dat',form='unformatted')
            open(unit=176,file=plotdir//'gsryu.dat',form='unformatted')
            open(unit=177,file=plotdir//'gsrzu.dat',form='unformatted')
            write(169) gsrxt
            write(176) gsryt
            write(177) gsrzt
            close(169)
            close(176)
            close(177)

         endif
      endif

 100  format (1x,1p10e16.6)

! The rest of the initialization is done also for everthing but the master

      if(taskid .gt. 0) then

! MS: chicrit is the zenith angle below which the Sun is visible.
! For points on the surface this is just pi/2, but at higher
! altitudes it is bigger.

        do k = 1,nl
          do j = 1,nf
            do i = 1,nz
              coschicrit(i,j,k) = cos(pie - 
     .                     asin( 1./ (1. + alts(i,j,k)/re) ))
            enddo
          enddo
        enddo


! put deni on mesh via linear interpolation
! and put on lower limit

! initialize all ions 

      j0 = 1
      do n = 1,nion
        do k = 1,nl
          do j = 1,nf
            do i = 1,nz
              jj = 1
              do while (  alts(i,j,k) .ge. zi(jj) .and. jj .le. 28 )
                j0 = jj
                jj = jj + 1
              enddo
              if ( n .eq. 1 ) nn = pthp
              if ( n .eq. 2 ) nn = pthep
              if ( n .eq. 3 ) nn = ptnp
              if ( n .eq. 4 ) nn = ptop
              if ( n .eq. 5 ) nn = ptn2p
              if ( n .eq. 6 ) nn = ptnop
              if ( n .eq. 7 ) nn = pto2p
              slope   = ( denii(j0+1,n) - denii(j0,n) ) 
     .                  / ( zi   (j0+1)   - zi   (j0) )
              deni(i,j,k,nn) = denii(j0,n) + 
     .                       ( alts(i,j,k) - zi(j0) ) * slope
              deni(i,j,k,nn) = amax1 ( deni(i,j,k,nn) , denmin )
!              deni(i,j,k,nn) = amax1 ( 10.*deni(i,j,k,nn) , denmin )
              if ( alts(i,j,k) .gt. zi(29) ) then
                if ( n .eq. 1 )  then
                  nn = pthp
                  deni(i,j,k,nn) = 
     .             amax1(denii(29,n)*zi(29)/alts(i,j,k),denmin)
                else
                  deni(i,j,k,nn) = denmin
                endif
              endif
            enddo
          enddo
        enddo
      enddo

!     initialize helium density = 10% hydrogen density

      do k = 1,nl
        do j = 1,nf
          do i = 1,nz
            deni(i,j,k,pthep) = 0.1 * deni(i,j,k,pthp)
          enddo
        enddo
      enddo

!      if ( taskid .eq. 1 ) then
!        do j = 1,nf
!          print *,'j',alts(nz/2,j,nl/2),deni(nz/2,j,k,pthp)
!        enddo
!      endif

!      do n = nion2+1,nion
!        do k = 1,nl
!          do j = 1,nf
!            do i = 1,nz
!              deni(i,j,k,n) = denmin
!            enddo
!          enddo
!        enddo
!      enddo

! print *,'done initializing deni',taskid

! initialize neutrals
! neutral density, temperature, and neutral wind

         if ( .not. restart ) then
           do nll = 2,nl-1
              call neutambt (hrinit,nll) 
           enddo
         endif

! electron and ion temperature initialization

         do k = nion1,nion2
            do n = 1,nl
               do j = 1,nf
                  do i = 1,nz
                     ti(i,j,n,k)    = tni(i,j,n)
                  enddo
               enddo
            enddo
         enddo

            do n = 1,nl
               do j = 1,nf
                  do i = 1,nz
                     te(i,j,n)      = tni(i,j,n)
                  enddo
               enddo
            enddo


!       average magnetic pole grid values (deni,Ti,Te)

        j0  = nf 

        do ni = nion1,nion2
          do i = 1,nz 
            deni_mnp0 = 0. 
            ti_mnp0   = 0. 
            do k = 2,nl-1
!              if ( alts (i,j0,k) .lt. alt_crit_avg) then 
                deni_mnp0      = deni_mnp0 + deni(i,j0,k,ni)
                ti_mnp0        = ti_mnp0 + ti(i,j0,k,ni)
!              endif
            enddo
            deni_mnp(i,ni) = deni_mnp0 / float(nl-2)
            ti_mnp(i,ni)   = ti_mnp0 / float(nl-2)
          enddo
        enddo

        do i = 1,nz 
          te_mnp0 = 0. 
          do k = 2,nl-1
!            if ( alts (i,j0,k) .lt. alt_crit_avg) then 
              te_mnp0     = te_mnp0 + te(i,j0,k)
!            endif
          enddo
          te_mnp(i)   = te_mnp0 / float(nl-2)
        enddo



! initialize ion velocity to zero 

         do nn = nion1,nion2
            do k = 1,nl
               do j = 1,nf
                  do i = 1,nz
                     vsi(i,j,k,nn)     = 0. 
                     sumvsi(i,j,k,nn)  = 0.
                  enddo
               enddo
            enddo
         enddo

      endif

c endif for taskid > 0 initialization

! read in photoabsorption rates

      if(taskid .eq. 0) then
         do i = 1,linesuv
            read (50,105) (sigabsdt(i,j), j=1,3)
 105        format (3f7.2) 
         enddo 
      endif

      if(taskid .eq. 0) then
         do j = 1,3 
            do i = 1,linesuv
               sigabsdt(i,j) = tm18 * sigabsdt(i,j) 
            enddo 
         enddo 
      endif
      call mpi_bcast(sigabsdt,linesuv*3,
     &     MPI_REAL,0,iComm,ierr)

! initialize photoionization rates to zero

      do j = 1,nneut
        do i = 1,linesuv
          sigidt(i,j)  = 0.
        enddo
        do i = 1,linesnt
          sigint(i,j)  = 0.
        enddo
      enddo

! read in daytime photoionization line data
! (only n, o, he, n_2, o_2)

      if(taskid .eq. 0) then
         do i = 1,linesuv
            read(60,106) (phionr(i,j), j=1,5)
            sigidt(i,ptn ) = phionr(i,1)
            sigidt(i,pto ) = phionr(i,2)
! can increase He+ photoproduction rate
! bailey and sellek used 2.5
! JK used 1.5
            sigidt(i,pthe) = phionr(i,3)
!
            sigidt(i,ptn2) = phionr(i,4)
            sigidt(i,pto2) = phionr(i,5)
         enddo
 106     format(5f7.2)
         
         do j = 1,nion
            do i = 1,linesuv
               sigidt(i,j) = tm18 * sigidt(i,j) 
            enddo 
         enddo 
         
! read in nighttime photoionization line data
! (only o, n_2, n0, o_2)

         do i = 1,linesnt
            read(61,106) (phionr(i,j), j=1,4)
            sigint(i,pto ) = phionr(i,1)
            sigint(i,ptn2) = phionr(i,2)
            sigint(i,ptno) = phionr(i,3)
            sigint(i,pto2) = phionr(i,4)
         enddo
         
         do j = 1,nion
            do i = 1,linesnt
               sigint(i,j) = tm18 * sigint(i,j) 
            enddo 
         enddo 
      endif
      call mpi_bcast(sigidt,linesuv*7,
     &     MPI_REAL,0,iComm,ierr)
      call mpi_bcast(sigint,linesnt*7,
     &     MPI_REAL,0,iComm,ierr)

!     below is altered from original
!     now just for flux spectra from harry warren

!      if(taskid .eq. 0) then
!         do i = 1,linesuv
!           flux(i) = fism(i)
!         enddo 
!       endif

 
! read in f74113, ai data and set euv flux 
! (from richards et al., jgr 99, 8981, 1994) 
 
       p  = 0.5 * ( f10p7 + fbar ) 
  
       if(taskid .eq. 0) then 
          do i = 1,linesuv 
             read (65,107) (fluxdat(i,j),j=1,2) 
             f74   = fluxdat(i,1) 
!             if ( i .eq. 1 ) f74 = 4.4 * f74
             ai    = fluxdat(i,2) 
             flux(i) = f74 * ( 1. + ai * ( p - 80.) ) * 1.e9 
          enddo  
       endif 
 107    format (f6.3,1pe11.4) 



      call mpi_bcast(flux,linesuv,
     &     MPI_REAL,0,iComm,ierr)


! read in angles for nighttime deposition fluxes

      if(taskid .eq. 0) then
         do i = 1,linesnt
            read(66,108) (thetant(i,j), j=1,4)
         enddo
      endif
 108  format (4f7.1)
      call mpi_bcast(thetant,linesnt*4,
     &     MPI_REAL,0,iComm,ierr)


! read in min/max altitude for nighttime deposition fluxes
!   zaltnt(i,1): zmin(i)
!   zaltnt(i,2): zmax(i)

      if(taskid .eq. 0) then
         do i = 1,linesnt
            read(67,108) (zaltnt(i,j), j=1,2)
         enddo
      endif

      call mpi_bcast(zaltnt,linesnt*2,
     &     MPI_REAL,0,iComm,ierr)

! Do this for everything but the master

      if (taskid .gt. 0) then

! call nighttime euv flux subroutines
! (lyman beta 1026, he i 584, he ii 304, lyman alpha 1216)

         do k = 1,nl
            do j = 1,nf
               call sf1026 ( f1026,1,j,k )
               call sf584  ( f584 ,2,j,k )
               call sf304  ( f304 ,3,j,k )
               call sf1216 ( f1216,4,j,k )
               do l = 1,91
                  do i = 1,nz
                     fluxnt(i,j,k,l,1) = f1026(i,j,k,l)            
                     fluxnt(i,j,k,l,2) = f584 (i,j,k,l)            
                     fluxnt(i,j,k,l,3) = f304 (i,j,k,l)            
                     fluxnt(i,j,k,l,4) = f1216(i,j,k,l)            
                  enddo
               enddo
            enddo
         enddo

! initialize e x b drift to 0

         do k = 1,nl
            do j = 1,nf
               do i = 1,nzp1
                  vexbs(i,j,k) = 0.
               enddo
            enddo
         enddo

         do k = 1,nl
            do j = 1,nfp1
               do i = 1,nz
                  vexbp(i,j,k) = 0.
               enddo
            enddo
         enddo

         do k = 1,nlp1
            do j = 1,nf
               do i = 1,nz
                  vexbh(i,j,k) = 0.
               enddo
            enddo
         enddo

! intialize diagnostic variables to 0

         do k = 1,nl
            do j = 1,nf
               do i = 1,nz
                  u1p(i,j,k) = 0.
                  u2s(i,j,k) = 0.
                  u3h(i,j,k) = 0.
                  u1(i,j,k) = 0.
                  u2(i,j,k) = 0.
                  u3(i,j,k) = 0.
                  u4(i,j,k) = 0.
                  u5(i,j,k) = 0.
               enddo
            enddo
         enddo

         do k = 1,nion
            do n = 1,nl
               do j = 1,nf
                  do i = 1,nz
                     t1(i,j,n,k) = 0.
                     t2(i,j,n,k) = 0.
                     t3(i,j,n,k) = 0.
                  enddo
               enddo
            enddo
         enddo
      endif

      if ( taskid .eq. 0 ) then
        deallocate (altstmp,glatstmp,glonstmp)
        deallocate (altst,glatst,glonst)
        deallocate (baltst,blatst,blonst)
        deallocate (xst,yst,zst)
        deallocate (altptmp,blatptmp,blonptmp)
!        deallocate (xalttmp,xblattmp,xblontmp)
!        deallocate (baltpt,blatpt,blonpt)  
!        deallocate (baltpt,blonpt)
        deallocate (baltpt)
        deallocate (vpsnxt,vpsnyt,vpsnzt)
        deallocate (vhsnxt,vhsnyt,vhsnzt)
        deallocate (xpt,ypt,zpt)
        deallocate (bdirsxt,bdirsyt,bdirszt)
        deallocate (gsthetaxt,gsthetayt,gsthetazt)
        deallocate (gsphixt,gsphiyt,gsphizt)
        deallocate (gsrxt,gsryt,gsrzt)
        deallocate (xrgt,xthgt,xphigt)
      endif

      print *,' finished initialization taskid = ',taskid

      return
      end

*******************************************
*******************************************

!            neutambt            

*******************************************
*******************************************


! calculate neutral densities and temperature
! from nrlmsise00

! no obtained from eq. (128) - bailey and balan (red book)

! neutral density and temperature

! input:
!    iyd - year and day as yyddd
!    sec - ut(sec)
!    alt - altitude(km) (greater than 85 km)
!    glat - geodetic latitude(deg)
!    glong - geodetic longitude(deg)
!    stl - local apparent solar time(hrs)
!    f107a - 3 month average of f10.7 flux
!    f107 - daily f10.7 flux for previous day
!    ap - magnetic index(daily) or when sw(9)=-1. :
!       - array containing:
!         (1) daily ap
!         (2) 3 hr ap index for current time
!         (3) 3 hr ap index for 3 hrs before current time
!         (4) 3 hr ap index for 6 hrs before current time
!         (5) 3 hr ap index for 9 hrs before current time
!         (6) average of eight 3 hr ap indicies from 12 to 33 hrs prior
!             to current time
!         (7) average of eight 3 hr ap indicies from 36 to 59 hrs prior
!             to current time
!    mass - mass number (only density for selected gas is
!             calculated.  mass 0 is temperature.  mass 48 for all.
! output: 
!    d(1) - he number density(cm-3)
!    d(2) - o number density(cm-3)
!    d(3) - n2 number density(cm-3)
!    d(4) - o2 number density(cm-3)
!    d(5) - ar number density(cm-3)
!    d(6) - total mass density(gm/cm3)
!    d(7) - h number density(cm-3)
!    d(8) - n number density(cm-3)
!    d(9) - anomalous O (see msis)
!    t(1) - exospheric temperature
!    t(2) - temperature at alt

! neutral winds

!    iyd - year and day as yyddd
!    sec - ut(sec)  (not important in lower atmosphere)
!    alt - altitude(km) 
!    glat - geodetic latitude(deg)
!    glong - geodetic longitude(deg)
!    stl - local apparent solar time(hrs)
!    f107a - 3 month average of f10.7 flux (use 150 in lower atmos.)
!    f107 - daily f10.7 flux for previous day ( " )
!    ap - two element array with
!         ap(1) = magnetic index(daily) (use 4 in lower atmos.)
!         ap(2)=current 3hr ap index (used only when sw(9)=-1.)
! note:  ut, local time, and longitude are used independently in the
!        model and are not of equal importance for every situation.  
!        for the most physically realistic calculation these three
!        variables should be consistent.
! output
!    w(1) = meridional (m/sec + northward)
!    w(2) = zonal (m/sec + eastward)


      subroutine neutambt (hr,nll)

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc' 

      real d(9),temp(2)
      real whm93(2),app(2)

      hruti = hr

      do j = 1,nf
        do i = 1,nz
          glonsij = glons(i,j,nll)
          if ( lcr )
!            glons0(i,j,nll) = glons(i,j,nll)+15.*hrinit
     .      glonsij = glons0(i,j,nll) -
     .                (hruti) * 15.
!     .                (hruti - hrinit) * 15.
            glonsij = mod(glonsij,360.)
          hrl   = mod(hruti + glonsij / 15.,24.)
          call msistim ( int(year),int(day),hrl,
     .                   glonsij,iyd,sec )
          call gtd7 ( iyd,sec,alts(i,j,nll)
     .         ,glats(i,j,nll),glonsij,
     .                hrl,fbar,f10p7,aap,mmass,d,temp )
          denni(i,j,nll,pth )  = snn(pth)  * d(7)
          denni(i,j,nll,pthe)  = snn(pthe) * d(1)
          denni(i,j,nll,ptn )  = snn(ptn)  * d(8)
          denni(i,j,nll,pto )  = snn(pto)  * d(2)
          denni(i,j,nll,ptn2)  = snn(ptn2) * d(3) + 1.e-30
          denni(i,j,nll,pto2)  = snn(pto2) * d(4) + 1.e-30
          tni(i,j,nll)         = stn * temp(2)
          denni(i,j,nll,ptno)  = 0.4 * exp( -3700. / tni(i,j,nll) ) 
     .                         * denni(i,j,nll,pto2) 
     .                         + 5.0e-7 * denni(i,j,nll,pto) 
        enddo
      enddo

      do j = 1,nf
        do i = 1,nz
          app(1)   = ap
          app(2)   = ap
          glonsij = glons(i,j,nll)
          if ( lcr )
     .      glonsij = glons0(i,j,nll) -
!     .                (hruti - hrinit) * 15.
     .                (hruti) * 15.
          glonsij = mod(glonsij,360.)
          hrl     = mod(hruti + glonsij / 15.,24.)
          call msistim ( int(year),int(day),hrl,glonsij,iyd,sec )
          if(lhwm93)  call gws5 ( iyd,sec,alts(i,j,nll),
     .                glats(i,j,nll),glonsij,
     .                hrl,fbar,f10p7,app,whm93        )
          if(lhwm14)  call hwm14 ( iyd,sec,alts(i,j,nll),
     .                glats(i,j,nll),glonsij,
     .                hrl,fbar,f10p7,app,whm93        )

          vi(i,j,nll)   = 100. * whm93(1) * tvn0 ! convert to cm/sec
          ui(i,j,nll)   = 100. * whm93(2) * tvn0 ! convert to cm/sec
          wi(i,j,nll)   = vw   * tvn0
        enddo
      enddo

      hrutf  = hr + .25

      do j = 1,nf
        do i = 1,nz
          glonsij = glons(i,j,nll)
          if ( lcr )
     .      glonsij = glons0(i,j,nll) -
!     .                (hrutf - hrinit) * 15.
     .                (hrutf) * 15.
          glonsij = mod(glonsij,360.)
          hrl   = mod(hrutf + glonsij / 15.,24.)
          call msistim ( int(year),int(day),hrl,
     .                   glonsij,iyd,sec )
          call gtd7 ( iyd,sec,alts(i,j,nll)
     .         ,glats(i,j,nll),glonsij,
     .                hrl,fbar,f10p7,aap,mmass,d,temp )
          dennf(i,j,nll,pth )  = snn(pth)  * d(7)
          dennf(i,j,nll,pthe)  = snn(pthe) * d(1)
          dennf(i,j,nll,ptn )  = snn(ptn)  * d(8)
          dennf(i,j,nll,pto )  = snn(pto)  * d(2)
          dennf(i,j,nll,ptn2)  = snn(ptn2) * d(3) + 1.e-30
          dennf(i,j,nll,pto2)  = snn(pto2) * d(4) + 1.e-30
          tnf(i,j,nll)         = stn * temp(2)
          dennf(i,j,nll,ptno)  = 0.4 * exp( -3700. / tnf(i,j,nll) ) 
     .                         * dennf(i,j,nll,pto2) 
     .                         + 5.0e-7 * dennf(i,j,nll,pto) 
        enddo
      enddo


      do j = 1,nf
        do i = 1,nz
          app(1)   = ap
          app(2)   = ap
          glonsij = glons(i,j,nll)
          if ( lcr )
     .      glonsij = glons0(i,j,nll) -
     .                (hrutf) * 15.
!     .                (hrutf - hrinit) * 15.

          glonsij = mod(glonsij,360.)
          hrl   = mod(hrutf + glonsij / 15.,24.)
          call msistim ( int(year),int(day),hrl,glonsij,iyd,sec )
          if(lhwm93)  call gws5 ( iyd,sec,alts(i,j,nll),
     .                glats(i,j,nll),glonsij,
     .                hrl,fbar,f10p7,app,whm93        )
          if(lhwm14)  call hwm14 ( iyd,sec,alts(i,j,nll),
     .                glats(i,j,nll),glonsij,
     .                hrl,fbar,f10p7,app,whm93        )

          vf(i,j,nll)   = 100. * whm93(1) * tvn0 
          uf(i,j,nll)   = 100. * whm93(2) * tvn0 
          wf(i,j,nll)   = vw   * tvn0

        enddo
      enddo

!     set density, temperature, velocity to current time

         do k = 1,nneut
            do n = 1,nl
               do j = 1,nf
                  do i = 1,nz
                     denn(i,j,n,k)  = denni(i,j,n,k)
                  enddo
               enddo
            enddo
         enddo

          do n = 1,nl
             do j = 1,nf
                do i = 1,nz
                   tn(i,j,n)      = tni(i,j,n)
                   u(i,j,n)       = ui(i,j,n)
                   v(i,j,n)       = vi(i,j,n)
                   w(i,j,n)       = wi(i,j,n)
!                   u1(i,j,n)      = u(i,j,n)
!                   u2(i,j,n)      = v(i,j,n)
                enddo
             enddo
          enddo


         return
         end

*******************************************
*******************************************

!            neut

*******************************************
*******************************************

         subroutine neut(hrut)

         include 'param3_mpi-1.98.inc'
         include 'com3_mpi-1.98.inc' 

         tfactor  = ( hrut - hruti ) / ( hrutf - hruti )
         tfactor1 = 1. - tfactor

         do k = 1,nneut
            do n = 1,nl
               do j = 1,nf
                  do i = 1,nz
                     denn(i,j,n,k) = denni(i,j,n,k) * tfactor1
     .                             + dennf(i,j,n,k) * tfactor
                  enddo
               enddo
            enddo
         enddo

            do n = 1,nl
               do j = 1,nf
                  do i = 1,nz
                     tn(i,j,n) =   tni(i,j,n) * tfactor1
     .                           + tnf(i,j,n) * tfactor
                     u(i,j,n)  =   ui(i,j,n)  * tfactor1
     .                           + uf(i,j,n)  * tfactor
                     v(i,j,n)  =   vi(i,j,n)  * tfactor1
     .                           + vf(i,j,n)  * tfactor
                     w(i,j,n)  =   wi(i,j,n)  * tfactor1
     .                           + wf(i,j,n)  * tfactor
                  enddo
               enddo
            enddo

          return
          end



*******************************************
*******************************************

!            transprt

*******************************************
*******************************************

      subroutine transprt (nfl,nll,p_crit)

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc'

      real prod(nz,nion),loss(nz,nion),lossr,
     .     phprodr(nz,nion),chrate(nz,nchem),
     .     chloss(nz,nion),chprod(nz,nion),relossr(nz,nion)
      real deni_old(nz,nion),te_old(nz),ti_old(nz,nion),vsi_old(nz,nion)
      real tvn(nz,nl)
      real nuin(nz,nion,nneut),
     .     nuij(nz,nion,nion),sumnuj(nz,nion)
      real vsin(nz,nion),vsidn(nz,nion),denin(nz,nion),cs(nz,nion)
      real ten(nz),tin(nz,nion)
      real p_crit(nnx-1)

! calculation of production and loss
!   phprodr: photo production rates
!   chrate:  chemical rates (ichem)
!   chloss:  chemical loss term
!   chprod:  chemical production term
!   relossr: recombination loss rates

! initialize tvn and gs and cfs (centrifugal force)

      do i = 1,nz
        tvn(i,nll) = 0.
        gs(i,nll)  = 0.
        cfs(i,nll)  = 0.
      enddo

      do i = 1,nz
        ne(i,nfl,nll)   = 1.
        te_old(i)       = te(i,nfl,nll)
        do j = nion1,nion2
          deni_old(i,j) = deni(i,nfl,nll,j)
          ne(i,nfl,nll) = ne(i,nfl,nll) + deni(i,nfl,nll,j)  
          ti_old(i,j)   = ti(i,nfl,nll,j)
          vsi_old(i,j)  = vsi(i,nfl,nll,j)
        enddo
      enddo

      call photprod ( phprodr,nfl,nll               ) ! calculates phprodr
      call chemrate ( chrate,nfl,nll                ) ! calculates chrate
      call chempl   ( chrate,chloss,chprod,nfl,nll  ) ! calcualtes chloss,chprod
      call recorate ( relossr,nfl,nll               ) ! calculates relossr

      do i = 1,nz 
        do j = nion1,nion2
          prod  (i,j) =  phprodr(i,j) * denn(i,nfl,nll,j)        
     .                   + chprod(i,j)
          lossr       =  relossr(i,j) * deni(i,nfl,nll,j) * 
     .                   ne(i,nfl,nll) 
     .                   + chloss(i,j)
          loss (i,j)  =  lossr / deni(i,nfl,nll,j)
        enddo

!     loss term for hydrogen and helium

        if ( alts(i,nfl,nll) .gt. pcrit*re ) then
          loss(i,pthp)  = loss(i,pthp)  + 1./decay_time 
          loss(i,pthep) = loss(i,pthep) + 1./decay_time 
! add loss term for O+ (JK)
!          loss(i,ptop) = loss(i,ptop) + 1./decay_time 
        endif

        gs(i,nll)   =  gzero * xrg(i,nfl,nll)
     .                 * ( re / (re + alts(i,nfl,nll)) ) ** 2

!        if (nll.eq.nl/2 .and. nfl.eq.44)
!     .   print *,'old s',i,gs(i,nll)

!        gs(i,nll)   = -gzero 
!     .                 * ( re / (re + alts(i,nfl,nll)) ) ** 2
!     .                 * ( gsrx(i,nfl,nll)*bdirsx(i,nfl,nll) +
!     .                     gsry(i,nfl,nll)*bdirsy(i,nfl,nll) +
!     .                     gsrz(i,nfl,nll)*bdirsz(i,nfl,nll)  )


!JK     centrifugal force (see notes 2012/01/04)
        fzero = 3.369
        clat = cos(pie*glats(i,nfl,nll)/180.0)
        slat = sin(pie*glats(i,nfl,nll)/180.0)
        cfs(i,nll)   =  -fzero * 
     .                 (clat*xrg(i,nfl,nll) + slat*xthg(i,nfl,nll))
     .                 * (re + alts(i,nfl,nll)) * clat / re 

!        if (nll.eq.nl/2 .and. nfl.eq.44)
!     .   print *,'new s',i,gs(i,nll)

!       approximation: not good for offset dipole
!       note: sign of gp expicitly accounted for
!       in derivation, i.e., g = -gp phat
!       so gp is positive here

!       gp(i,nfl,nll) = sqrt( gzero**2 - gs(i,nll)**2 ) 

!        if (nll.eq.nl/2 .and. nfl.eq.144)
!     .   print *,'old p',i,gp(i,nfl,nll)

!       should be good for offset dipole
!       note: sign of gp expicitly accounted for
!       in derivation, i.e., g = -gp phat
!       so gp is positive here (JH)

        gp(i,nfl,nll)   = gzero 
     .                 * ( re / (re + alts(i,nfl,nll)) ) ** 2
     .                 * ( gsrx(i,nfl,nll)*vpsnx(i,nfl,nll) +
     .                     gsry(i,nfl,nll)*vpsny(i,nfl,nll) +
     .                     gsrz(i,nfl,nll)*vpsnz(i,nfl,nll)  )


!        if (nll.eq.nl/2 .and. nfl.eq.144)
!     .   print *,'new p',i,gp(i,nfl,nll)

        vnq(i,nfl,nll) = v(i,nfl,nll) *
     .                   ( gsthetax(i,nfl,nll) * bdirsx(i,nfl,nll) +
     .                     gsthetay(i,nfl,nll) * bdirsy(i,nfl,nll) +
     .                     gsthetaz(i,nfl,nll) * bdirsz(i,nfl,nll)   ) +
     .                   u(i,nfl,nll) *
     .                   ( gsphix(i,nfl,nll) * bdirsx(i,nfl,nll) +
     .                     gsphiy(i,nfl,nll) * bdirsy(i,nfl,nll) +
     .                     gsphiz(i,nfl,nll) * bdirsz(i,nfl,nll)   )   +
     .                   w(i,nfl,nll) *
     .                   ( gsrx(i,nfl,nll) * bdirsx(i,nfl,nll) +
     .                     gsry(i,nfl,nll) * bdirsy(i,nfl,nll) +
     .                     gsrz(i,nfl,nll) * bdirsz(i,nfl,nll)   ) 

        vnp(i,nfl,nll) = v(i,nfl,nll) *
     .                   ( gsthetax(i,nfl,nll) * vpsnx(i,nfl,nll) +
     .                     gsthetay(i,nfl,nll) * vpsny(i,nfl,nll) +
     .                     gsthetaz(i,nfl,nll) * vpsnz(i,nfl,nll)   ) +
     .                   u(i,nfl,nll) *
     .                   ( gsphix(i,nfl,nll) * vpsnx(i,nfl,nll) +
     .                     gsphiy(i,nfl,nll) * vpsny(i,nfl,nll) +
     .                     gsphiz(i,nfl,nll) * vpsnz(i,nfl,nll)   )   +
     .                   w(i,nfl,nll) *
     .                   ( gsrx(i,nfl,nll) * vpsnx(i,nfl,nll) +
     .                     gsry(i,nfl,nll) * vpsny(i,nfl,nll) +
     .                     gsrz(i,nfl,nll) * vpsnz(i,nfl,nll)   ) 

        vnphi(i,nfl,nll) = v(i,nfl,nll) *
     .                   ( gsthetax(i,nfl,nll) * vhsnx(i,nfl,nll) +
     .                     gsthetay(i,nfl,nll) * vhsny(i,nfl,nll) +
     .                     gsthetaz(i,nfl,nll) * vhsnz(i,nfl,nll)   ) +
     .                   u(i,nfl,nll) *
     .                   ( gsphix(i,nfl,nll) * vhsnx(i,nfl,nll) +
     .                     gsphiy(i,nfl,nll) * vhsny(i,nfl,nll) +
     .                     gsphiz(i,nfl,nll) * vhsnz(i,nfl,nll)   )   +
     .                   w(i,nfl,nll) *
     .                   ( gsrx(i,nfl,nll) * vhsnx(i,nfl,nll) +
     .                     gsry(i,nfl,nll) * vhsny(i,nfl,nll) +
     .                     gsrz(i,nfl,nll) * vhsnz(i,nfl,nll)   ) 

!        u3(i,nfl,nll) = vnq(i,nfl,nll)
!        u4(i,nfl,nll) = vnp(i,nfl,nll)

        tvn(i,nll)    = vnq(i,nfl,nll) 
          
      enddo

      call update ( tvn,nuin,sumnuj,nuij,nfl,nll )

      do i = 1,nz
        do nni = nion1,nion2
          sumvsi(i,nfl,nll,nni) = 0.
          do nj = nion1,nion2
          sumvsi(i,nfl,nll,nni) =   sumvsi(i,nfl,nll,nni) +
     .                              nuij(i,nni,nj)*vsi(i,nfl,nll,nj)
          enddo
        enddo
      enddo

! define new arrays for velocity and density

      do ni = nion1,nion2
        do i = 1,nz
          vsin (i,ni) = vsi(i,nfl,nll,ni)
          vsidn(i,ni) = vsid(i,nfl,nll,ni)
          denin(i,ni) = deni(i,nfl,nll,ni)
        enddo
      enddo

! define sound velocity used in vsisolv

      do ni = nion1,nion2
        do i = 1,nz
          cfac     = 1.6667 * 8.6174e-5 * te(i,nfl,nll) / ami(ni)
          cs(i,ni) = 9.79e5 * sqrt(cfac)
        enddo
      enddo

! update variables

      do ni = nion1,nion2

        call vsisolv ( vsin(1,ni),vsidn(1,ni),vsi_old(1,ni),
     .                 sumnuj(1,ni),nfl,nll,cs(1,ni) )

! compensating filter

       call smoothz ( vsin(1,ni), 1 )

! put stuff back into velocity array

        do i = 1,nz
          vsi(i,nfl,nll,ni)  = vsin(i,ni)
          vsid(i,nfl,nll,ni) = vsidn(i,ni)
!          if ( alts(i,nfl,nll) .gt. 6.*re .and. ni .eq. pthp ) then
!            vsi(i,nfl,nll,ni)  = 0.
!            vsid(i,nfl,nll,ni) = 0.
!          endif
        enddo

        call densolv2 ( ni,denin(1,ni),
     .       prod(1,ni),loss(1,ni),deni_old(1,ni),nfl,nll )

! put stuff back into density array

        do i = 1,nz
          deni(i,nfl,nll,ni) = denin(i,ni)
        enddo

! put floor on density

        do i = 1,nz
          deni(i,nfl,nll,ni) = amax1 ( deni(i,nfl,nll,ni), denmin )
! below commented out (JK)
          if ( alts(i,nfl,nll) .gt. pcrit*re .and. ni .eq. pthp ) 
     .         deni(i,nfl,nll,ni) = amax1 ( deni(i,nfl,nll,ni), .1 )
          if ( alts(i,nfl,nll) .gt. pcrit*re .and. ni .eq. pthep ) 
     .         deni(i,nfl,nll,ni) = amax1 ( deni(i,nfl,nll,ni), .01 )

        enddo


      enddo

!JK   Basic plasmapause:  for L>__ and nz/2-_ < i < nz/2+_, reduce density
!     L = 4 is height 3*re = 19,110 km
!     L = 3 is height 2*re = 12,740 km
!      altloss = 4.0*re
!      framp = 0.99
!      nz_1 = nz/2 - 3
!      nz_2 = nz/2 + 3
!      do i=nz_1,nz_2
!        if (alts(i,nfl,nll) .gt. altloss) then
!          do ni=1,4
!            if (ni .eq. 1) ip = pthp
!            if (ni .eq. 2) ip = pthep
!            if (ni .eq. 3) ip = ptop
!            if (ni .eq. 4) ip = ptnp
!            deni(i,nfl,nll,ip)=framp*deni(i,nfl,nll,ip)
!            deni(i,nfl,nll,ip) = amax1 ( deni(i,nfl,nll,ip), denmin ) 
!          enddo
!        endif
!      enddo

! JH kill high altitude density 

!      framp = 0.999
!      nz_1 = nz/2 - 7
!      nz_2 = nz/2 + 7

!     epsden is decay rate of density (i.e., loss rate)

!!      epsden = 1.e-2
!!      framp  = (1. - epsden)
!!      pcrit = 6.

!      pcrit = 6.
!      nz_1 = nz/2 - 20
!      nz_2 = nz/2 + 20
!       kk = (taskid - 1) * (nl - 2) + (nll - 1)
!       if (alts(nz/2,nfl,nll)/re .ge. p_crit(kk) ) then

!!        do i=1,nz
!!          if (alts(i,nfl,nll)/re .ge. pcrit ) then
!!            do ni=1,4
!!              if (ni .eq. 1) ip  = pthp
!!              if (ni .eq. 2) ip  = pthep
!!              if (ni .eq. 3) ip  = ptop
!!              if (ni .eq. 4) ip  = ptnp
!!              deni(i,nfl,nll,ip) = framp*deni(i,nfl,nll,ip)
!!              deni(i,nfl,nll,ip) = amax1 ( deni(i,nfl,nll,ip), 10. )
!!            enddo
!!          endif
!!        enddo

! define new arrays for temperature

      do ni = nion1,nion2
        do i = 1,nz
          tin(i,ni)  = ti(i,nfl,nll,ni)
        enddo
      enddo

      do i = 1,nz
        ten(i)  = te(i,nfl,nll)
      enddo

! temperatures (with floors and warnings)

      tn0 = 200. ! floor on temperature

      call etemp  (ten,te_old,phprodr,nfl,nll)
      do i = 1,nz
        te(i,nfl,nll)  = amax1(tn(i,nfl,nll),ten(i))
!        te(i,nfl,nll)  = amax1(tn0,ten(i))
        te(i,nfl,nll)  = amin1(te(i,nfl,nll),1.e4)
        if ( te(i,nfl,nll) .lt. 0 ) then
          print *,' T(e) negative: i,nfl,nll taskid',i,nfl,nll,taskid
          stop
        endif
      enddo



      call htemp  (tin(1,pthp) ,ti_old(1,pthp) ,tvn,nuin,nfl,nll)
      do i = 1,nz
        ti(i,nfl,nll,pthp)  = amax1(tn(i,nfl,nll),tin(i,pthp))
!        ti(i,nfl,nll,pthp)  = amax1(tn0,tin(i,pthp))
        ti(i,nfl,nll,pthp)  = amin1(ti(i,nfl,nll,pthp),1.e4)
        if ( ti(i,nfl,nll,pthp) .lt. 0 ) then
          print *,' T(H) negative: i,nfl,nll',i,nfl,nll 
          stop
        endif
      enddo

      call hetemp (tin(1,pthep),ti_old(1,pthep),tvn,nuin,nfl,nll)
      do i = 1,nz
        ti(i,nfl,nll,pthep)  = amax1(tn(i,nfl,nll),tin(i,pthep))
!        ti(i,nfl,nll,pthep)  = amax1(tn0,tin(i,pthep))
        ti(i,nfl,nll,pthep)  = amin1(ti(i,nfl,nll,pthep),1.e4)
        if ( ti(i,nfl,nll,pthep) .lt. 0 ) then
          print *,' T(He) negative: i,nfl,nll',i,nfl,nll 
          stop
        endif
      enddo

      call otemp  (tin(1,ptop) ,ti_old(1,ptop) ,tvn,nuin,nfl,nll)
      do i = 1,nz
        ti(i,nfl,nll,ptop)  = amax1(tn(i,nfl,nll),tin(i,ptop))
!        ti(i,nfl,nll,ptop)  = amax1(tn0,tin(i,ptop))
        ti(i,nfl,nll,ptop)  = amin1(ti(i,nfl,nll,ptop),1.e4)
        if ( ti(i,nfl,nll,ptop) .lt. 0 ) then
          print *,' T(O) negative: i,nfl,nll',i,nfl,nll 
          stop
        endif
      enddo

      do i = 1,nz
        ti(i,nfl,nll,ptnp )    = ti(i,nfl,nll,ptop)
        ti(i,nfl,nll,ptn2p)    = ti(i,nfl,nll,ptop)
        ti(i,nfl,nll,ptnop)    = ti(i,nfl,nll,ptop)
        ti(i,nfl,nll,pto2p)    = ti(i,nfl,nll,ptop)
      enddo

      return
      end


*******************************************
*******************************************

!            photprod

*******************************************
*******************************************

! photoproduction rates 

      subroutine photprod ( phprodr,nfl,nll )

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc' 

      real phprodr(nz,nion),xmass(3)
      integer idx(3)
      real*8 ch1,atm_chapman

! scale height of neutral atmosphere

      hcof = 1.e-5 * bolt / ( gzero * amu * re ** 2 )

      nll1 = nll

!      do iz = 1,nz
!         coschi = cx(iz,nfl,nll)
!         do j = nion1,nion2
!         phprodr ( iz,j ) = 0.
!         enddo
!      enddo

!      return

      do iz = 1,nz
         coschi = cx(iz,nfl,nll)
         do j = nion1,nion2
            phprodr ( iz,j ) = 0.
         enddo

! only consider o, n2, o2 for absorption

         idx(1) = pto
         idx(2) = ptn2
         idx(3) = pto2

         rp    = alts(iz,nfl,nll1) + re
         rp2   = rp * rp  
         
         if ( coschi .ge. coschicrit(iz,nfl,nll1) ) then ! sun is up

! daytime deposition 

            do i = 1,3 
               hscale   = hcof * tn(iz,nfl,nll1) * rp2 / amn(idx(i))
               xscale   = rp / hscale
               ch1      = atm_chapman(xscale,rtod*acos(coschi))
               if ( ch1 .gt. 1.e22 ) ch1 = 1.e22
               xmass(i) = denn(iz,nfl,nll1,idx(i)) * hscale * ch1 * 1.e5
            enddo

            do l=1,linesuv
               exa =  xmass(1) * sigabsdt(l,1) 
     .              + xmass(2) * sigabsdt(l,2) 
     .              + xmass(3) * sigabsdt(l,3)
               if(exa .gt. 35.) exa = 35.
               flx = flux(l) * exp(-exa) 
               do j=nion1,nion2
                  phprodr(iz,j) = phprodr(iz,j) + sigidt(l,j) * flx
               enddo
            enddo

! photoelectron ionization

               pei_rate = exp(-(alts(iz,nfl,nll)-150.)/40.)
               if ( alts(iz,nfl,nll) .gt. 200. ) pei_rate = 0.2

               phprodr(iz,ptop)  = phprodr(iz,ptop)  * (1. + pei_rate)
               phprodr(iz,ptn2p) = phprodr(iz,ptn2p) * (1. + pei_rate)
               phprodr(iz,pto2p) = phprodr(iz,pto2p) * 1.2

! add nighttime ionization

! ignore for now

            ang    = acos ( coschi )
            itheta0 = int ( ang / po180 ) - 90
            itheta  = int ( amax1 ( float(itheta0), 1. ) )
            del     = ang/po180 - int(ang/po180)

            do l = 1,linesnt
               do j=nion1,nion2
                if (itheta0 .lt. 1) then
                  fluxntt = fluxnt(iz,nfl,nll1,itheta,l)
                else
                  fluxntt = fluxnt(iz,nfl,nll1,itheta,l) * (1.-del) +
     .                      fluxnt(iz,nfl,nll1,itheta+1,l) * del
                endif
                phprodr(iz,j) =   phprodr(iz,j) 
     .                          + sigint(l,j) * fluxntt
                  
               enddo
            enddo

         else                   ! sun is dowm

! nighttime deposition

! ignore for now

            ang    = acos ( coschi )
!            itheta = nint ( ang / po180 ) - 90
!            itheta = int ( amax1 ( float(itheta), 1. ) )

            itheta0 = int ( ang / po180 ) - 90
            itheta  = int ( amax1 ( float(itheta0), 1. ) )
            del     = ang/po180 - int(ang/po180)

            do l = 1,linesnt
               do j=nion1,nion2
                if (itheta0 .lt. 1) then
                  fluxntt = fluxnt(iz,nfl,nll1,itheta,l)
                else
                  fluxntt = fluxnt(iz,nfl,nll1,itheta,l) * (1.-del) +
     .                      fluxnt(iz,nfl,nll1,itheta+1,l) * del
                endif
               phprodr(iz,j) =   phprodr(iz,j) 
     .                          + sigint(l,j) * fluxntt
!                phprodr(iz,j) = 0.                  
               enddo
            enddo

         endif
      enddo

      return
      end



*******************************************
*******************************************

!            photprod

*******************************************
*******************************************

! photoproduction rates (not interpolated at night)

      subroutine photprod_old ( phprodr,nfl,nll )

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc' 

      real phprodr(nz,nion),xmass(3)
      integer idx(3)
      real*8 atm_chapman

! scale height of neutral atmosphere

      hcof = 1.e-5 * bolt / ( gzero * amu * re ** 2 )

      nll1 = nll

!      do iz = 1,nz
!         coschi = cx(iz,nfl,nll)
!         do j = nion1,nion2
!         phprodr ( iz,j ) = 0.
!         enddo
!      enddo

!      return

      do iz = 1,nz
         coschi = cx(iz,nfl,nll)
         do j = nion1,nion2
            phprodr ( iz,j ) = 0.
         enddo

! only consider o, n2, o2 for absorption

         idx(1) = pto
         idx(2) = ptn2
         idx(3) = pto2

         rp    = alts(iz,nfl,nll1) + re
         rp2   = rp * rp  
         
         if ( coschi .ge. coschicrit(iz,nfl,nll1) ) then ! sun is up

! daytime deposition 

            do i = 1,3 
               hscale   = hcof * tn(iz,nfl,nll1) * rp2 / amn(idx(i))
               xscale   = rp / hscale
               ch1      = atm_chapman(xscale,rtod*acos(coschi))
               if ( ch1 .gt. 1.e22 ) ch1 = 1.e22
               xmass(i) = denn(iz,nfl,nll1,idx(i)) * hscale * ch1 * 1.e5
            enddo

            do l=1,linesuv
               exa =  xmass(1) * sigabsdt(l,1) 
     .              + xmass(2) * sigabsdt(l,2) 
     .              + xmass(3) * sigabsdt(l,3)
               if(exa .gt. 35.) exa = 35.
               flx = flux(l) * exp(-exa) 
               do j=nion1,nion2
                  phprodr(iz,j) = phprodr(iz,j) + sigidt(l,j) * flx
               enddo
            enddo

! photoelectron ionization

               pei_rate = exp(-(alts(iz,nfl,nll)-150.)/40.)
               if ( alts(iz,nfl,nll) .gt. 200. ) pei_rate = 0.2

               phprodr(iz,ptop)  = phprodr(iz,ptop)  * (1. + pei_rate)
               phprodr(iz,ptn2p) = phprodr(iz,ptn2p) * (1. + pei_rate)
               phprodr(iz,pto2p) = phprodr(iz,pto2p) * 1.2

! add nighttime ionization

            ang    = acos ( coschi )
            itheta = nint ( ang / po180 ) - 90
            itheta = int ( amax1 ( float(itheta), 1. ) )
            do l = 1,linesnt
               do j=nion1,nion2
                  phprodr(iz,j) =   phprodr(iz,j) 
     .                 + sigint(l,j) * fluxnt(iz,nfl,nll1,itheta,l)
               enddo
            enddo

! nighttime deposition

         else                   ! sun is dowm

            ang    = acos ( coschi )
            itheta = nint ( ang / po180 ) - 90
            itheta = int ( amax1 ( float(itheta), 1. ) )
            do l = 1,linesnt
               do j=nion1,nion2
                  phprodr(iz,j) =   phprodr(iz,j) 
     .                 + sigint(l,j) * fluxnt(iz,nfl,nll1,itheta,l)
               enddo

            enddo

         endif
      enddo

      return
      end

*******************************************
*******************************************

!            chemrate

*******************************************
*******************************************

!     chemical producation and loss rates
!     bb: bailley and balan (red book, 1996)

      subroutine chemrate ( chrate,nfl,nll )

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc'

      real chrate(nz,nchem)

      do iz = 1,nz

      ti300o = ti(iz,nfl,nll,ptop) / 300.

! h+ + o --> o+ + h (bb)

      chrate (iz,1) = 2.2e-11 
     .            * sqrt( ti(iz,nfl,nll,pthp) )       

! he+ + n2 --> n2+ + he (bb)

      chrate (iz,2) = 3.5e-10                        

! he+ + n2 --> n+ + n + he (schunk)

      chrate (iz,3) = 8.5e-10                        

! he+ + o2 --> o+ + o + he (bb)

      chrate (iz,4) = 8.0e-10                        

! he+ + o2 --> o2+ + he

      chrate (iz,5) = 2.0e-10                        

! n+ + o2 --> no+ + o  (schunk)

      chrate (iz,6) = 2.0e-10                        

! n+ + o2 --> o2+ + n(2d) (schunk)

      chrate (iz,7) = 4.0e-10                        

! n+ + 0 --> o+ + n

      chrate (iz,8) = 1.0e-12                        

! n+ + no --> no+ + o (schunk)

      chrate (iz,9) = 2.0e-11                        

! o+ + h --> h+ + o   (bb)

      chrate(iz,10) = 2.5e-11 
     .             * sqrt( tn(iz,nfl,nll) )           
 
! o+ + n2 --> no+ + n (bb)

      chrate(iz,11) = 1.533e-12 -                    
     .             5.920e-13 * ti300o +
     .             8.600e-14 * ti300o ** 2

        if ( ti(iz,nfl,nll,ptop) .gt. 1700 ) 
     .    chrate(iz,11) = 2.730e-12 -
     .                    1.155e-12 * ti300o +
     .                    1.483e-13 * ti300o ** 2

! o+ + o2 --> o2+ + o

      chrate(iz,12) = 2.820e-11 -                    
     .             7.740e-12 * ti300o +
     .             1.073e-12 * ti300o ** 2 -
     .             5.170e-14 * ti300o ** 3 +
     .             9.650e-16 * ti300o ** 4

! o+ + no --> no+ + o

      chrate(iz,13) = 1.0e-12                        

! n2+ + o --> no+ + n(2d) (bb)

      chrate(iz,14) = 1.4e-10 / ti300o ** .44        

! n2+ + o2 --> o2+ + n2 (schunk)

      chrate(iz,15) = 5.0e-11 / sqrt( ti300o )       

! n2+ + o2 --> no+ + no

      chrate(iz,16) = 1.0e-14                        

! n2+ + no --> no+ + n2 (schunk)

      chrate(iz,17) = 3.3e-10                        

! o2+ + n --> no+ + o (schunk)

      chrate(iz,18) = 1.2e-10                        

! o2+ + n(2d) --> n+ + o2

      chrate(iz,19) = 2.5e-10                        

! o2+ + no --> no+ + o2 (bb)

      chrate(iz,20) = 4.4e-10                        

! o2+ + n2 --> no+ + no (schunk)

      chrate(iz,21) = 5.0e-16                        

      enddo

      return
      end

*******************************************
*******************************************

!            recorate

*******************************************
*******************************************

! recombination rates 
! bb: bailley and balan (red book, 1996)

      subroutine recorate ( relossr,nfl,nll )

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc' 

      real relossr(nz,nion)

      do iz = 1,nz

        te300 = te(iz,nfl,nll) / 300.

        relossr(iz,pthp)  = 4.43e-12 / te300 ** .7 
        relossr(iz,pthep) = relossr(iz,pthp)
        relossr(iz,ptnp)  = relossr(iz,pthp)
        relossr(iz,ptop)  = relossr(iz,pthp)
        relossr(iz,ptn2p) = 1.8e-7 / te300 ** 0.39     !   (schunk)
        relossr(iz,ptnop) = 4.2e-7 / te300 ** 0.85     !   (bb)
        relossr(iz,pto2p) = 1.6e-7 / te300 ** 0.55     !   (schunk)

      enddo

      return
      end


*******************************************
*******************************************

!          chempl

*******************************************
*******************************************
        
! chemical loss (chloss) and production (chprod)

! chrate: chemical reaction rates calculated in chemrate
! ichem: input data file showing loss, neutral, production
!        species for each reaction 

      subroutine chempl ( chrate,chloss,chprod,nfl,nll )

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc' 

      real chrate(nz,nchem),chloss(nz,nion),chprod(nz,nion)

      do i = nion1,nion2
        do iz = 1,nz
          chloss(iz,i)   = 0.
          chprod(iz,i)   = 0.
        enddo
      enddo

      do k = 1,nchem
        il = ichem(k,1) ! ion species (reacting) loss
        in = ichem(k,2) ! neutral species reacting
        ip = ichem(k,3) ! ion species produced
        do iz = 1,nz
           chem  = denn(iz,nfl,nll,in) * chrate(iz,k)
           tdeni = deni(iz,nfl,nll,il) * chem
           chloss(iz,il) = tdeni + chloss(iz,il)
           chprod(iz,ip) = tdeni + chprod(iz,ip)           
        enddo
      enddo

      return
      end


*******************************************
*******************************************

!            update

*******************************************
*******************************************

      subroutine update ( tvn,nuin,sumnuj,nuij,nfl,nll )

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc'

      real nuin(nz,nion,nneut),nuij(nz,nion,nion)
      real nuint(nz,nion)
      real sumnuj(nz,nion),nufacij,nufacin
      real tvn(nz,nl)
      real k0,mi

!!      do k = 1,nl
!!        print *,taskid,denn(51,10,3,k)
!!      enddo

! ion-neutral collision frequency

! initialize everything to 0

      do nn = 1,nneut
        do ni = nion1,nion2
          do iz = 1,nz
            nuin (iz,ni,nn) = 0.
            nuint(iz,ni)    = 0.
          enddo
        enddo
      enddo

! collision frequencies/factors

! hydrogen (H)

      ni = pthp
      do nn = 1,nneut
        do i = 1,nz
          if ( nn .eq. pto ) then
            teff    = ti(i,nfl,nll,ni) 
            fac     = ( 1.00 - .047 * alog10(teff) ) ** 2
            tfactor = sqrt(teff) * fac
            nuin(i,ni,nn)  = 6.61e-11 * denn(i,nfl,nll,nn) * tfactor
          else
            amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
            amimn   = amn(nn) / ( ami(ni) + amn(nn) )
            nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
            nuin(i,ni,nn) = nufacin * denn(i,nfl,nll,nn)
          endif
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo      

! helium (He)

      ni = pthep
      do nn = 1,nneut
        do i = 1,nz
          amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
          amimn   = amn(nn) / ( ami(ni) + amn(nn) )
          nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
          nuin(i,ni,nn) = nufacin * denn(i,nfl,nll,nn)
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo      

! nitrogen (N)

      ni = ptnp
      do nn = 1,nneut
        do i = 1,nz
          amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
          amimn   = amn(nn) / ( ami(ni) + amn(nn) )
          nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
          nuin(i,ni,nn) = nufacin * denn(i,nfl,nll,nn)
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo      

! oxygen (O)

      ni = ptop
      do nn = 1,nneut
        do i = 1,nz
          if ( nn .eq. pto ) then
            teff    = 0.5 * ( ti(i,nfl,nll,ni) + tn(i,nfl,nll) )
            fac     = ( 1.04 - .067 * alog10(teff) ) ** 2
            tfactor = sqrt(teff) * fac
            nuin(i,ni,nn)  = 4.45e-11 * denn(i,nfl,nll,nn) * tfactor
          else
            amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
            amimn   = amn(nn) / ( ami(ni) + amn(nn) )
            nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
            nuin(i,ni,nn) = nufacin * denn(i,nfl,nll,nn)
          endif
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo      

! nitrogen 2(N2)

      ni = ptn2p
      do nn = 1,nneut
        do i = 1,nz
          if ( nn .eq. ptn2 ) then
            teff    = 0.5 * ( ti(i,nfl,nll,ni) + tn(i,nfl,nll) )
            fac     = ( 1.00 - .069 * alog10(teff) ) ** 2
            tfactor = sqrt(teff) * fac
            nuin(i,ni,nn) = 5.14e-11 * denn(i,nfl,nll,nn) * tfactor
          else
            amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
            amimn   = amn(nn) / ( ami(ni) + amn(nn) )
            nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
            nuin(i,ni,nn) = nufacin * denn(i,nfl,nll,nn)
          endif
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo      

! nitrous oxide (N0)

      ni = ptnop
      do nn = 1,nneut
        do i = 1,nz
          amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
          amimn   = amn(nn) / ( ami(ni) + amn(nn) )
          nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
          nuin(i,ni,nn) = nufacin * denn(i,nfl,nll,nn)
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo      

! oxygen 2(O2)

      ni = pto2p
      do nn = 1,nneut
        do i = 1,nz
          if ( nn .eq. pto2 ) then
            teff    = 0.5 * ( ti(i,nfl,nll,ni) + tn(i,nfl,nll) )
            fac     = ( 1.00 - .073 * alog10(teff) ) ** 2
            tfactor = sqrt(teff) * fac
            nuin(i,ni,nn) = 2.59e-11 * denn(i,nfl,nll,nn) * tfactor
          else
            amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
            amimn   = amn(nn) / ( ami(ni) + amn(nn) )
            nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
            nuin(i,ni,nn) = nufacin * denn(i,nfl,nll,nn)
          endif
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo      

! ion-ion collision frequency

      do ni = nion1,nion2
        do i = 1,nz
          do nj = nion1,nion2
            if(ni .ne. nj) then
              alame1  = ( ami(ni) + ami(nj) ) * evtok /
     .                ( ami(ni)*ti(i,nfl,nll,nj) + 
     .                  ami(nj)*ti(i,nfl,nll,ni) ) 
              alame2  = deni(i,nfl,nll,ni) * evtok / ti(i,nfl,nll,ni) +
     .                  deni(i,nfl,nll,nj) * evtok / ti(i,nfl,nll,nj)
              if ( alame2 .lt. 0 ) then
              print *,'ni,i,nj,nfl,nll,tii,tij,alame1,alame2,nii,nij',
     .              ni,i,nj,nfl,nll,ti(i,nfl,nll,ni),ti(i,nfl,nll,nj),
     .              alame1,alame2,
     .              deni(i,nfl,nll,ni),deni(i,nfl,nll,nj)
                stop
              endif
              alame   = alame1 * sqrt(alame2)
              alam    = 23. - alog(alame)
              amufac  = (ami(nj)/ami(ni))/(ami(ni) +ami(nj))
              nufacij = 9.2e-2*alam*sqrt(amufac)
              nuij(i,ni,nj) =  nufacij * deni(i,nfl,nll,nj)
     .                        / sqrt( ti(i,nfl,nll,ni)**3 )
            else 
              nuij(i,ni,nj) = 0.
            endif
          enddo
        enddo
      enddo

!     add this for restart (sarah mcdonald)
!     it's used in term 1 below

      do i = 1,nz
        do nni = nion1,nion2
          sumvsi(i,nfl,nll,nni) = 0.
          do nj = nion1,nion2
          sumvsi(i,nfl,nll,nni) =   sumvsi(i,nfl,nll,nni) +
     .                              nuij(i,nni,nj)*vsi(i,nfl,nll,nj)
          enddo
        enddo
      enddo

! sumnuj: sum of ion-ion coll freq and nuin

      do ni = nion1,nion2
        do i = 1,nz
          sumnuj(i,ni) = 0.
          do nj = nion1,nion2
            sumnuj(i,ni) = sumnuj(i,ni) + nuij(i,ni,nj)
          enddo
          sumnuj(i,ni) = sumnuj(i,ni) + nuint(i,ni)
        enddo     
      enddo              

! update ne

      do i = 1,nz
      ne(i,nfl,nll) = 1.
        do ni = nion1,nion2
          ne(i,nfl,nll) = ne(i,nfl,nll) + deni(i,nfl,nll,ni)
        enddo
      enddo

! get a new value for vsid

      do i = 2,nz-1
        do ni = nion1,nion2
          mi    = amu * ami(ni)
          k0    = bolt / mi
          term1 = nuint(i,ni) * tvn(i,nll) + 
     .            sumvsi(i,nfl,nll,ni) + gs(i,nll) + cfs(i,nll)
          pip   = 0.5 * (   deni(i+1,nfl,nll,ni) * ti(i+1,nfl,nll,ni)
     .                    + deni(i,nfl,nll,ni)   * ti(i,nfl,nll,ni)   )
          pim   = 0.5 * (   deni(i,nfl,nll,ni)   * ti(i,nfl,nll,ni)
     .                    + deni(i-1,nfl,nll,ni) * ti(i-1,nfl,nll,ni) )
          denid = 
     .             (        deni(i-1,nfl,nll,ni) 
     .               + 4. * deni(i,nfl,nll,ni) 
     .               +      deni(i+1,nfl,nll,ni)  ) / 6.
          term2 =  - bms(i,nfl,nll) * k0 /  denid
     .             * ( pip - pim ) / d22s(i,nfl,nll)
          pep   = 0.5 * (   ne(i+1,nfl,nll) * te(i+1,nfl,nll)
     .                    + ne(i,nfl,nll)   * te(i,nfl,nll)   )
          pem   = 0.5 * (   ne(i,nfl,nll)   * te(i,nfl,nll)
     .                    + ne(i-1,nfl,nll) * te(i-1,nfl,nll) )
          dened = 
     .  ( ne(i-1,nfl,nll) + 4. * ne(i,nfl,nll) + ne(i+1,nfl,nll) ) / 6.
          term3 =  - bms(i,nfl,nll) * k0 /  dened
     .             * ( pep - pem ) / d22s(i,nfl,nll)

          vsid(i,nfl,nll,ni)  =  term1 + term2 + term3 

c$$$          if ( deni(i,nfl,nll,ni) .le. .005*ne(i,nfl,nll) ) 
c$$$     .     vsid(i,nfl,nll,ni) =   vsid(i,nfl,nll,ni) 
c$$$     .         * exp ( -.005*ne(i,nfl,nll)/deni(i,nfl,nll,ni) )

        enddo
      enddo 

! fix up end points for vsid

      do ni = nion1,nion2
        vsid (1,nfl,nll,ni)    = vsid (2,nfl,nll,ni)
        vsid (nz,nfl,nll,ni)   = vsid (nz-1,nfl,nll,ni)
      enddo

! calculate the electron-neutral collision frequency
! nuen = 5.4e-10*n_n*T_e^1/2 (kelley, the earth's ionosphere, p. 462)

      do i = 1,nz
        nuen(i,nfl,nll) = 0
        do nn = 1,nneut
          nuen(i,nfl,nll) = nuen(i,nfl,nll) + 5.4e-10 *
     .                      denn(i,nfl,nll,nn) * sqrt(te(i,nfl,nll))
        enddo
       enddo

!      do i = 1,nz
!        nuen(i,nfl,nll) = 3.3e-10 * sqrt(te(i,nfl,nll)) *
!     .   ( denn(i,nfl,nll,pto2) +
!     .     2. * denn(i,nfl,nll,pto) + denn(i,nfl,nll,ptno) )
!     .  + 2.33e-11 * denn(i,nfl,nll,ptn2) * te(i,nfl,nll)
!      enddo

! calculate pedersen and hall conductivities

      do i = 1,nz
        dene    = ne(i,nfl,nll)
        oce     = 1.76e7 * bmag * bms(i,nfl,nll)
        sige    = dene * charge * sol / ( bmag * bms(i,nfl,nll) )
        cole    = nuen(i,nfl,nll) / oce
        denome  = 1. + cole * cole
        sigpe   = sige * cole / denome
        sighe   = sige * cole * cole / denome
        sigpi   = 0.
        sighi   = 0.
        sighic  = 0.
        sigpic  = 0.
        do ni = nion1,nion2
          oci    = 9580. * bmag * bms(i,nfl,nll) / ami(ni)
          sigi   = deni(i,nfl,nll,ni) * charge * sol / 
     .             ( bmag * bms(i,nfl,nll) )
          coli   = nuint(i,ni) / oci
          denomi = 1. + coli * coli
          sigpi  = sigpi  + sigi * coli / denomi
          sigpic = sigpic + sigi * coli / denomi / oci
          sighi  = sighi  + sigi * coli * coli / denomi
          sighic = sighic + sigi / denomi / oci
        enddo
        sigmap(i,nfl,nll)   = sigpi + sigpe
        sigmah(i,nfl,nll)   = sighi - sighe
        sigmapic(i,nfl,nll) = sigpic
        sigmahic(i,nfl,nll) = sighic
        if (alts(i,nfl,nll) .ge. 1.e4) then
          sigmap(i,nfl,nll)   = 0.00
          sigmah(i,nfl,nll)   = 0.00
          sigmapic(i,nfl,nll) = 0.00
          sigmahic(i,nfl,nll) = 0.00
        endif
      enddo

      if ( .not. hall ) then
        do i=1,nz
          sigmah(i,nfl,nll)   = 0.
          sigmahic(i,nfl,nll) = 0.
        enddo
      endif

      hipcp(nfl,nll)    = 0. 
      hipcphi(nfl,nll)  = 0. 
      hihcm(nfl,nll)    = 0. 
      hidpv(nfl,nll)    = 0. 
      hidphiv(nfl,nll)  = 0. 
      hidpg(nfl,nll)    = 0. 
      hidphig(nfl,nll)  = 0. 

      hipc(nfl,nll)     = 0. 
      hihc(nfl,nll)     = 0. 
      hidv(nfl,nll)     = 0. 
     
      do i = 1,nz
        ang     = .5 * pie - blats(i,nfl,nll) * pie / 180.
        bang    = blats(nz,nfl,nll) * pie / 180.
        del     = sqrt ( 1. + 3. * cos(ang) * cos(ang) )
        b       = bmag * bms(i,nfl,nll)

! original hipcp multiplied by p
! now it is multiplied by .25 * tan(bang)

        hipcp(nfl,nll) = hipcp(nfl,nll) +
     .                   sigmap(i,nfl,nll) * del / bms(i,nfl,nll) *
!!     .                   dels(i,nfl,nll) * 
!!     .                   0.25 / tan(bang)
     .                   dels(i,nfl,nll) 

! original hipcphi divided by p
! now it is divided by tan(bang)

        hipcphi(nfl,nll) = hipcphi(nfl,nll) +
     .                     sigmap(i,nfl,nll) / del / bms(i,nfl,nll) *
!!     .                     dels(i,nfl,nll) * tan(bang)
     .                     dels(i,nfl,nll) 

! original hihcm
! now it is divided by 2

        hihcm(nfl,nll) = hihcm(nfl,nll) +
     .                   sigmah(i,nfl,nll) / bms(i,nfl,nll) *
!!     .                   dels(i,nfl,nll) * 0.5
     .                   dels(i,nfl,nll) 


! original fdpv
! now it is multiplied by 0.5

        fdpv = (bmag/sol) * ( sigmap(i,nfl,nll) * vnphi(i,nfl,nll) +
     .                        sigmah(i,nfl,nll) * vnp(i,nfl,nll)     )
!!     .                    * 0.5

        hidpv(nfl,nll) = hidpv(nfl,nll) +
     .                   brs(i,nfl,nll) * 1.e5  * sin(ang) *
     .                   fdpv * dels(i,nfl,nll)    

! original fdpg
! now it is multiplied by 0.5

        fdpg = (bmag/sol) * sigmapic(i,nfl,nll) * gp(i,nfl,nll)
!!     .                    * 0.5


        hidpg(nfl,nll) = hidpg(nfl,nll) +
     .                   brs(i,nfl,nll) * 1.e5  * sin(ang) *
     .                   fdpg * dels(i,nfl,nll)    

! fdphiv and fdphig now divide by
! sin^2(bang)*tan(bang)

!        angfac = tan(bang)/(cos(bang)*cos(bang))

        fdphiv = (bmag/sol) * ( -sigmap(i,nfl,nll) * vnp(i,nfl,nll) +
     .                           sigmah(i,nfl,nll) * vnphi(i,nfl,nll)  )
!     .                      * angfac

        hidphiv(nfl,nll) = hidphiv(nfl,nll) +
     .                     re * 1.e5 * ( sin(ang) ** 3 ) / del * 
     .                     fdphiv * dels(i,nfl,nll) 


        fdphig = (bmag/sol) * sigmahic(i,nfl,nll) * gp(i,nfl,nll)
!     .                      * angfac

        hidphig(nfl,nll) = hidphig(nfl,nll) +
     .                     re * 1.e5 * ( sin(ang) ** 3 ) / del * 
     .                     fdphig * dels(i,nfl,nll) 

!       integrated quantities for current

        hipc(nfl,nll) = hipc(nfl,nll) +
     .                  sigmap(i,nfl,nll) * del / re / sin(ang) ** 3 *
     .                  dels(i,nfl,nll) / 1.e5
        hihc(nfl,nll) = hihc(nfl,nll) +
     .                  sigmah(i,nfl,nll) / brs(i,nfl,nll) / sin(ang) *
     .                  dels(i,nfl,nll) / 1.e5
        hidv(nfl,nll) = hidv(nfl,nll) +
     .                  fdpv * dels(i,nfl,nll)    

      enddo

        hipcp(nfl,nll)   = hipcp(nfl,nll)   * ps(nz/2,nfl,nll)
        hipcphi(nfl,nll) = hipcphi(nfl,nll) / ps(nz/2,nfl,nll)


! calculate collisional ion velocity 
! not used; simply a diagnostic

!      do i = 1,nz
!        do ni = nion1,nion2
!          vsic(i,nfl,nll,ni) = vsid(i,nfl,nll,ni) / sumnuj(i,ni)
!        enddo
!      enddo

      return
      end

*******************************************
*******************************************

!            htemp

*******************************************
*******************************************

      subroutine htemp ( tti,tiold,tvn,nuin,nfl,nll )

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc'

      real tiold(nz),kapi(nz),s1i(nz),s2i(nz),s3i(nz),s4i(nz),s5i(nz)
      real tvn(nz,nll),nuin(nz,nion,nneut),s6i(nz),s7i(nz),tti(nz)
      real lambda
      real divvexb(nz)

      convfac = amu / bolt / 3.

      do i = 1,nz
        s1i(i)  = 0.
        s2i(i)  = 0.
        s3i(i)  = 0.
        s4i(i)  = 0.
        s5i(i)  = 0.
        s6i(i)  = 0.
        s7i(i)  = 0.
        kapi(i) = 0.
      enddo

      do i = 1,nz

! from schunk/nagy book

        lambda = 23. - 
     .           0.5*alog(deni(i,nfl,nll,pth)/
     .           (ti(i,nfl,nll,pth)/evtok)**3) 
         schunkfac = 0. 
         do ni = nion1,nion2 
           if (ni .ne. pth) then 
             schunkfac = schunkfac + 
     .                   deni(i,nfl,nll,ni)/deni(i,nfl,nll,pth) *  
     .                   sqrt(ami(ni)/(ami(ni)+ami(pth))**5) *  
     .          (3*ami(pth)**2 + 1.6*ami(pth)*ami(ni) + 1.3*ami(ni)**2) 
            endif 
          enddo 
        kapi(i) = 15.*3.1e4 * sqrt ( ti(i,nfl,nll,pth)**5 ) / 
     .            sqrt(ami(pth)) /  
     .            (1. + 1.75*schunkfac) / lambda 

        kapi(i)  = 0.6667 * kapi(i) * evtok 

! neutrals

        do nn = 1,nneut
          redmass = 
     .     ami(pthp) * amn(nn) / ( ami(pthp) + amn(nn) ) ** 2 
          s2i(i) = s2i(i) + 2. * nuin(i,pthp,nn) * redmass 
          s3i(i) = s3i(i)
     .      + convfac * amn(nn) 
     .                * abs ( vsi(i,nfl,nll,pthp) - tvn(i,nll) ) ** 2
     .      * 2. * nuin(i,pthp,nn) * redmass 
        enddo

        s1i(i) = s2i(i) * tn(i,nfl,nll) 

! electrons 

        s4i(i) = 7.7e-6 * ne(i,nfl,nll) / ami(pthp) 
     .                   / te(i,nfl,nll) / sqrt(te(i,nfl,nll))
     .                   * .66667 * evtok
        s5i(i) = s4i(i) * te(i,nfl,nll)

! other ions

        do ni = nion1,nion2
          if ( ni .ne. pthp ) then
            tfac    =    ti(i,nfl,nll,pthp) / ami(pthp) 
     .                +  ti(i,nfl,nll,ni) / ami(ni)
            xs6i    = 3.3e-4 * deni(i,nfl,nll,ni) / ami(pthp) / ami(ni)
     .                / tfac / sqrt(tfac) * .66667 * evtok
            xs7i    = xs6i * ti(i,nfl,nll,ni)
            s6i(i) = s6i(i) + xs6i
            s7i(i) = s7i(i) + xs7i 
          endif
        enddo

      enddo

! MS: Neglected term, divergence of ExB drift 
! Divergence of the ExB drift; requires equatorial drift 

      nzh = nz / 2 
      vexbeq = vexb(nzh,nfl,nll)
      do i = 1,nz 
        divvexb(i) = 6.*vexbeq / 
     .               (ps(i,nfl,nll)*re*1.e5) * 
     .               cos(blats(i,nfl,nll)*po180)**2 * 
     .               (1.+sin(blats(i,nfl,nll)*po180)**2) / 
     .               (1.+3.*sin(blats(i,nfl,nll)*po180)**2)**2 
        s2i(i) = s2i(i) - 0.3333 * divvexb(i)
      enddo 

      call tisolv(tti,tiold,kapi,s1i,s2i,s3i,s4i,s5i,s6i,s7i,pthp,
     .            nfl,nll)

      return
      end

*******************************************
*******************************************

!            hetemp

*******************************************
*******************************************

      subroutine hetemp ( tti,tiold,tvn,nuin,nfl,nll )

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc'

      real tiold(nz),kapi(nz),s1i(nz),s2i(nz),s3i(nz),s4i(nz),s5i(nz)
      real tvn(nz,nl),nuin(nz,nion,nneut),s6i(nz),s7i(nz),tti(nz)
      real lambda 
      real divvexb(nz)

      convfac = amu / bolt / 3.

      do i = 1,nz
        s1i(i)  = 0.
        s2i(i)  = 0.
        s3i(i)  = 0.
        s4i(i)  = 0.
        s5i(i)  = 0.
        s6i(i)  = 0.
        s7i(i)  = 0.
        kapi(i) = 0.
      enddo

      do i = 1,nz

! from schunk/nagy book

        lambda = 23. - 
     .           0.5*alog(deni(i,nfl,nll,pthe)/
     .           (ti(i,nfl,nll,pthe)/evtok)**3) 
         schunkfac = 0. 
         do ni = nion1,nion2 
           if (ni .ne. pthe) then 
             schunkfac = schunkfac + 
     .                   deni(i,nfl,nll,ni)/deni(i,nfl,nll,pthe) *  
     .                   sqrt(ami(ni)/(ami(ni)+ami(pthe))**5) *  
     .          (3*ami(pthe)**2 + 1.6*ami(pthe)*ami(ni) + 
     .             1.3*ami(ni)**2) 
            endif 
          enddo 
        kapi(i) = 15.*3.1e4 * sqrt ( ti(i,nfl,nll,pthe)**5 ) / 
     .            sqrt(ami(pthe)) /  
     .            (1. + 1.75*schunkfac) / lambda 

        kapi(i)  = 0.6667 * kapi(i) * evtok 

! neutrals

        do nn = 1,nneut
          redmass = 
     .     ami(pthep) * amn(nn) / ( ami(pthep) + amn(nn) ) ** 2 
          s2i(i) = s2i(i) + 2. * nuin(i,pthep,nn) * redmass 
          s3i(i) = s3i(i)
     .      + convfac * amn(nn) 
     .                * abs ( vsi(i,nfl,nll,pthep) - tvn(i,nll) ) ** 2 
     .      * 2. * nuin(i,pthep,nn) * redmass 
        enddo

        s1i(i) = s2i(i) * tn(i,nfl,nll) 

! electrons

        s4i(i) = 7.7e-6 * ne(i,nfl,nll) / ami(pthep) 
     .                   / te(i,nfl,nll) / sqrt(te(i,nfl,nll))
     .                   * .66667 * evtok
        s5i(i) = s4i(i) * te(i,nfl,nll)

! other ions

        do ni = nion1,nion2
          if ( ni .ne. pthep ) then
            tfac    =   ti(i,nfl,nll,pthep) / ami(pthep) 
     .                + ti(i,nfl,nll,ni) / ami(ni)
            xs6i    = 3.3e-4 * deni(i,nfl,nll,ni) / ami(pthep) / ami(ni)
     .                / tfac / sqrt(tfac) * .66667 * evtok
            xs7i    = xs6i * ti(i,nfl,nll,ni)
            s6i(i) = s6i(i) + xs6i
            s7i(i) = s7i(i) + xs7i 
          endif
        enddo

      enddo

! MS: Neglected term, divergence of ExB drift 
! Divergence of the ExB drift; requires equatorial drift 

      nzh = nz / 2 
      vexbeq = vexb(nzh,nfl,nll)
      do i = 1,nz 
        divvexb(i) = 6.*vexbeq / 
     .               (ps(i,nfl,nll)*re*1.e5) * 
     .               cos(blats(i,nfl,nll)*po180)**2 * 
     .               (1.+sin(blats(i,nfl,nll)*po180)**2) / 
     .               (1.+3.*sin(blats(i,nfl,nll)*po180)**2)**2 
        s2i(i) = s2i(i) - 0.3333 * divvexb(i)
      enddo 

      call tisolv(tti,tiold,kapi,s1i,s2i,s3i,s4i,s5i,s6i,s7i,pthep,
     .            nfl,nll)

      return
      end

*******************************************
*******************************************

!            otemp

*******************************************
*******************************************

      subroutine otemp ( tti,tiold,tvn,nuin,nfl,nll )

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc'

      real tiold(nz),kapi(nz),s1i(nz),s2i(nz),s3i(nz),s4i(nz),s5i(nz)
      real tvn(nz,nl),nuin(nz,nion,nneut),s6i(nz),s7i(nz),tti(nz)
      real lambda
      real divvexb(nz)

      convfac = amu / bolt / 3.

      do i = 1,nz
        s1i(i)  = 0.
        s2i(i)  = 0.
        s3i(i)  = 0.
        s4i(i)  = 0.
        s5i(i)  = 0.
        s6i(i)  = 0.
        s7i(i)  = 0.
        kapi(i) = 0.
      enddo

      do i = 1,nz

! from schunk/nagy book

        lambda = 23. - 
     .           0.5*alog(deni(i,nfl,nll,pto)/
     .           (ti(i,nfl,nll,pto)/evtok)**3) 
         schunkfac = 0. 
         do ni = nion1,nion2 
           if (ni .ne. pto) then 
             schunkfac = schunkfac + 
     .                   deni(i,nfl,nll,ni)/deni(i,nfl,nll,pto) *  
     .                   sqrt(ami(ni)/(ami(ni)+ami(pto))**5) *  
     .          (3*ami(pto)**2 + 1.6*ami(pto)*ami(ni) + 1.3*ami(ni)**2) 
            endif 
          enddo 
        kapi(i) = 15.*3.1e4 * sqrt ( ti(i,nfl,nll,pto)**5 ) / 
     .            sqrt(ami(pto)) /  
     .            (1. + 1.75*schunkfac) / lambda 

        kapi(i)  = 0.6667 * kapi(i) * evtok 

! neutrals

        do nn = 1,nneut
          redmass = 
     .     ami(ptop) * amn(nn) / ( ami(ptop) + amn(nn) ) ** 2 
          s2i(i) = s2i(i) + 2. * nuin(i,ptop,nn) * redmass 
          s3i(i) = s3i(i)
     .      + convfac * amn(nn) 
     .                * abs ( vsi(i,nfl,nll,ptop) - tvn(i,nll) )
     .      * 2. * nuin(i,ptop,nn) * redmass 
        enddo

        s1i(i) = s2i(i) * tn(i,nfl,nll) 

! electrons

        s4i(i) = 7.7e-6 * ne(i,nfl,nll) / ami(ptop) 
     .                   / te(i,nfl,nll) / sqrt(te(i,nfl,nll))
     .                   * .66667 * evtok
        s5i(i) = s4i(i) * te(i,nfl,nll)

! other ions

        do ni = nion1,nion2
          if ( ni .ne. ptop ) then
            tfac    =    ti(i,nfl,nll,ptop) / ami(ptop) 
     .                 + ti(i,nfl,nll,ni) / ami(ni)
            xs6i    = 3.3e-4 * deni(i,nfl,nll,ni) / ami(ptop) / ami(ni)
     .                / tfac / sqrt(tfac) * .66667 * evtok
            xs7i    = xs6i * ti(i,nfl,nll,ni)
            s6i(i) = s6i(i) + xs6i
            s7i(i) = s7i(i) + xs7i 
          endif
        enddo

      enddo

! MS: Neglected term, divergence of ExB drift 
! Divergence of the ExB drift; requires equatorial drift 

      nzh = nz / 2 
      vexbeq = vexb(nzh,nfl,nll)
      do i = 1,nz 
        divvexb(i) = 6.*vexbeq / 
     .               (ps(i,nfl,nll)*re*1.e5) * 
     .               cos(blats(i,nfl,nll)*po180)**2 * 
     .               (1.+sin(blats(i,nfl,nll)*po180)**2) / 
     .               (1.+3.*sin(blats(i,nfl,nll)*po180)**2)**2 
        s2i(i) = s2i(i) - 0.3333 * divvexb(i)
      enddo 

      call tisolv(tti,tiold,kapi,s1i,s2i,s3i,s4i,s5i,s6i,s7i,ptop,
     .            nfl,nll)

      return
      end



*******************************************
*******************************************

!            etemp

*******************************************
*******************************************

      subroutine etemp ( tte,te_old,phprodr,nfl,nll )

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc'


      real tte(nz),te_old(nz),kape(nz)
      real s1e(nz),s2e(nz),s3e(nz),s4e(nz),phprodr(nz,nion)
      real s5e(nz),qphe(nz),phprod(nz)
      real qen(nz,nneut)
      real ne300s,ne300n,n2300
      real ratio(nz)
      real divvexb(nz)

      do i = 1,nz
        s1e(i)  = 0.
        s2e(i)  = 0.
        s3e(i)  = 0.
        s4e(i)  = 0.
        kape(i) = 0.
        do ni = 1,nneut
          qen(i,ni) = 0.
        enddo
      enddo

      do i = 1,nz

        fac1 = denn(i,nfl,nll,pto)  * 1.1e-16  
     .          * ( 1. + 5.7e-4 * te(i,nfl,nll) )
        fac2 = denn(i,nfl,nll,ptn2) * 2.82e-17 
     .          * ( 1  - 1.2e-4 * te(i,nfl,nll) )* sqrt(te(i,nfl,nll))
        fac3 = denn(i,nfl,nll,pto2) * 2.2e-16  
     .         * ( 1. + 3.6e-2  * sqrt(te(i,nfl,nll)) )
        akpefac = fac1 + fac2 + fac3

        kape(i) = 7.7e5 * sqrt ( te(i,nfl,nll)**5 ) * 0.6667 * evtok 
     .      / ( 1. + 3.22e4 * ( te(i,nfl,nll)**2 / 
     .                          ne(i,nfl,nll) * akpefac) )


! neutrals (Tn - Te) term

! N2

! vibrational state from red book (p. 269) milward et al.

        qen(i,ptn2) = .6667 *  evtok * denn(i,nfl,nll,ptn2) *
     .                  ( 1.2e-19 * ( 1. - 1.2e-4 * te(i,nfl,nll) ) 
     .                            * te(i,nfl,nll) +
     .                    2.e-14 / sqrt(te(i,nfl,nll)) 
     .                    + 6.5e-22 * ( tn(i,nfl,nll) - 310 ) ** 2 *
     .                      exp(.0023*(te(i,nfl,nll) - tn(i,nfl,nll)))) 

! O2

        qen(i,pto2) = .6667 * evtok * denn(i,nfl,nll,pto2) *
     .                 ( 7.9e-19 * ( 1. + 3.6e-2 * sqrt(te(i,nfl,nll)))
     .                           *  sqrt(te(i,nfl,nll)) +
     .                   7.e-14 / sqrt(te(i,nfl,nll)) )

! O

        qen(i,pto) = .6667 * 7.2e-18 * evtok * denn(i,nfl,nll,pto) *
     .                  sqrt(te(i,nfl,nll))

! H

        qen(i,pth) = .6667 * 6.3e-16 * evtok * denn(i,nfl,nll,pth) *
     .                  ( 1. - 1.35e-4 * te(i,nfl,nll) ) * 
     .                  sqrt(te(i,nfl,nll))

        do nn = 1,nneut
          s2e(i) = s2e(i) + qen(i,nn)
        enddo

        s1e(i) = s2e(i) * tn(i,nfl,nll) 

! ions (Ti - Te) term

        do ni = nion1,nion2
          xs3e    = 7.7e-6 * deni(i,nfl,nll,ni) / ami(ni) 
     .                     / te(i,nfl,nll) / sqrt(te(i,nfl,nll))
     .                     * .66667 * evtok
          xs4e    = xs3e * ti(i,nfl,nll,ni)
          s3e(i) = s3e(i) + xs3e
          s4e(i) = s4e(i) + xs4e 
        enddo

      enddo

! photoelectron heating
! red book (millward et al. p. 269)

! calculate total ion photoproduction (= photoelectron)

      do i = 1,nz
        phprod(i)   = 0.
        do ni = nion1,nion2
          phprod(i) = phprodr(i,ni) * denn(i,nfl,nll,ni) + phprod(i)
        enddo
      enddo

! iz300s/iz300n are redefined here

      do i = 1,nz
        ratio(i) = ne(i,nfl,nll) / 
     .             (0.1*denn(i,nfl,nll,pto)+
     .              denn(i,nfl,nll,pto2)+denn(i,nfl,nll,ptn2)) 
      enddo

      i = 1 
      do while ( ratio(i) .le. 3.e-3 .and. i .lt. nz ) 
         iz300s(nfl,nll) = i 
         i         = i + 1 
      enddo 
 
      i = nz 
      do while ( ratio(i) .le. 3.e-3 .and. i .gt. 1 )  
         iz300n(nfl,nll) = i 
         i         = i - 1 
      enddo 

      if ( iz300s(nfl,nll) .gt. iz300n(nfl,nll) ) then

        do i = 1,nz
            xarg =   ne(i,nfl,nll) 
     .             / (        denn(i,nfl,nll,pto2) 
     .                 +      denn(i,nfl,nll,ptn2) 
     .                 + .1 * denn(i,nfl,nll,pto)   )
            x    = alog ( xarg )
            earg =     12.75 
     .               + 6.941 * x 
     .               + 1.166 * x ** 2 
     .               + 0.08034 * x ** 3
     .               + 0.001996 * x ** 4
            epsi = exp ( -earg )
            qphe(i) = epsi * phprod(i)
          enddo
        else
          do i = 1,iz300s(nfl,nll)
            xarg =   ne(i,nfl,nll) 
     .             / (        denn(i,nfl,nll,pto2) 
     .                 +      denn(i,nfl,nll,ptn2) 
     .                 + .1 * denn(i,nfl,nll,pto)   )
            x    = alog ( xarg )
            earg =     12.75 
     .               + 6.941 * x 
     .               + 1.166 * x ** 2 
     .               + 0.08034 * x ** 3
     .               + 0.001996 * x ** 4
            epsi = exp ( -earg )
            qphe(i) = epsi * phprod(i)
          enddo

! smooth things at '300' km 

        izs   = iz300s(nfl,nll)
        facts = (3.e-3-ratio(izs)) / 
     .          (ratio(izs+1)-ratio(izs))
        ne300s = ne(izs,nfl,nll) + (ne(izs+1,nfl,nll)-
     .                              ne(izs,nfl,nll)) * facts
        o2300 = denn(izs,nfl,nll,pto2) + 
     .         (denn(izs+1,nfl,nll,pto2)-
     .          denn(izs,nfl,nll,pto2)) * facts
        n2300 = denn(izs,nfl,nll,ptn2) + 
     .         (denn(izs+1,nfl,nll,ptn2)-
     .          denn(izs,nfl,nll,ptn2)) * facts
        o300 = denn(izs,nfl,nll,pto) + 
     .         (denn(izs+1,nfl,nll,pto)-denn(izs,nfl,nll,pto)) * facts
        phprod300 = phprod(izs) + 
     .        (phprod(izs+1)-phprod(izs)) * facts
        xarg300 = ne300s / ( o2300 + n2300 + 0.1*o300 )
        x300 = alog( xarg300)
        earg300 =     12.75 + 
     .        6.941 * x300 + 
     .        1.166 * x300 ** 2 + 
     .        0.08034 * x300 ** 3 + 
     .        0.001996 * x300 ** 4
        epsi300 = exp ( -earg300 )
        q0s = epsi300 * phprod300 / ne300s

        do i = iz300n(nfl,nll),nz
          xarg =   ne(i,nfl,nll) 
     .           / (       denn(i,nfl,nll,pto2) 
     .              +      denn(i,nfl,nll,ptn2) 
     .              + .1 * denn(i,nfl,nll,pto) )
          x    = alog ( xarg )
          earg =     12.75 
     .             + 6.941 * x 
     .             + 1.166 * x ** 2 
     .             + 0.08034 * x ** 3
     .             + 0.001996 * x ** 4
          epsi = exp ( -earg )
          qphe(i) = epsi * phprod(i)
        enddo

        izn   = iz300n(nfl,nll)
        factn = (3.e-3-ratio(izn)) / 
     .           (ratio(izn-1)-ratio(izn))
        ne300n = ne(izn,nfl,nll) + 
     .        (ne(izn-1,nfl,nll)-ne(izn,nfl,nll)) * factn
        o2300 = denn(izn,nfl,nll,pto2) + 
     .        (denn(izn-1,nfl,nll,pto2)-
     .         denn(izn,nfl,nll,pto2)) * factn
        n2300 = denn(izn,nfl,nll,ptn2) + 
     .        (denn(izn-1,nfl,nll,ptn2)-
     .         denn(izn,nfl,nll,ptn2)) * factn
        o300 = denn(izn,nfl,nll,pto) + 
     .        (denn(izn-1,nfl,nll,pto)-denn(izn,nfl,nll,pto)) * factn
        phprod300 = phprod(izn) + 
     .        (phprod(izn-1)-phprod(izn)) * factn
        xarg300 = ne300n / ( o2300 + n2300 + 0.1*o300 )
        x300 = alog( xarg300)
        earg300 =     12.75 + 
     .        6.941 * x300 + 
     .        1.166 * x300 ** 2 + 
     .        0.08034 * x300 ** 3 + 
     .        0.001996 * x300 ** 4
        epsi300 = exp ( -earg300 )
        q0n = epsi300 * phprod300 / ne300n

        xbms = bms(izs,nfl,nll) + 
     .         (bms(izs+1,nfl,nll)-bms(izs,nfl,nll)) * facts
        xbmn = bms(izn,nfl,nll) + 
     .         (bms(izn-1,nfl,nll)-bms(izn,nfl,nll)) * factn

        dels300s = dels(iz300s(nfl,nll),nfl,nll) * facts
        dels300n = dels(iz300n(nfl,nll)-1,nfl,nll) * factn

        ! MS: Old code used a wasteful way to calculate xn. 
        ! Cleaner version here. 
        xn = 0. 
        ! Set bottom integration bound to 300 km. 
        xn =   xn + 0.5 * ( ne(iz300n(nfl,nll)-1,nfl,nll) + ne300n ) * 
     .        (dels(iz300n(nfl,nll)-1,nfl,nll) - dels300n ) 
        do i =iz300n(nfl,nll)-2,iz300s(nfl,nll)+1,-1 
           xn = xn + 0.5 * ( ne(i,nfl,nll) + ne(i+1,nfl,nll) ) * 
     .                       dels(i,nfl,nll) 
        enddo 

        if ( q0s .lt. 0 .or. q0n .lt. 0 ) then
          print *,' q0s = ',q0s,' q0n = ',q0n,' nfl = ',nfl
        endif

! 1/22/00

! put in dels (arc length along field line)

        xs    = 0.
         do i = iz300s(nfl,nll)+1,iz300n(nfl,nll)-1 
           if (i .eq. iz300s(nfl,nll)+1) then 
             xs = xs + 0.5*( ne300s + ne(i,nfl,nll) ) * 
     .           (dels(iz300s(nfl,nll),nfl,nll) - dels300s) 
           else 
              xs = xs + 0.5 * ( ne(i,nfl,nll) + ne(i-1,nfl,nll) ) 
     .                         * dels(i-1,nfl,nll) 
              xn = xn - 0.5 * ( ne(i,nfl,nll) + ne(i-1,nfl,nll) ) 
     .                         * dels(i-1,nfl,nll) 
           endif 
 
           xints = cqe*xs  
           xintn = cqe*xn 
           xqs    = ne(i,nfl,nll) * q0s * bms(i,nfl,nll) / 
     .                              xbms * exp(-xints) 
           xqn    = ne(i,nfl,nll) * q0n * bms(i,nfl,nll) / 
     .                              xbmn * exp(-xintn) 
           qphe(i) = xqs + xqn 
        enddo 
      endif

      do i = 1,nz
        s5e(i) = 0.66667 * evtok * qphe(i) / ne(i,nfl,nll) !* .15
      enddo 

! MS: Neglected term, divergence of ExB drift 
! Divergence of the ExB drift; requires equatorial drift 

      nzh    = nz / 2 
      vexbeq = vexb(nzh,nfl,nll)
      do i = 1,nz 
        divvexb(i) = 6.*vexbeq / 
     .                 (ps(i,nfl,nll)*re*1.e5) * 
     .                 cos(blats(i,nfl,nll)*po180)**2 * 
     .                 (1.+sin(blats(i,nfl,nll)*po180)**2) / 
     .                 (1.+3.*sin(blats(i,nfl,nll)*po180)**2)**2 
        s2e(i) = s2e(i) - 0.3333 * divvexb(i)
      enddo 

      call tesolv(tte,te_old,kape,s1e,s2e,s3e,s4e,s5e,nfl,nll)

      return
      end


*******************************************
*******************************************

!            densolv2

*******************************************
*******************************************

      subroutine densolv2( ni,tdeni,prod,loss,oldion,nfl,nll )

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc' 

      real tdeni(nz)
      real oldion(nz), prod(nz), loss(nz)
      real a(nz), b(nz), c(nz), d(nz)

! initialize

      do j = 1,nz
        a(j) = 0.
        b(j) = 0.
        c(j) = 0.
        d(j) = 0.
      enddo


      do j = 2,nz-1

        ujm1  = vsi(j-1,nfl,nll,ni)/bms(j-1,nfl,nll) 
        uj    = vsi(j,nfl,nll,ni)  /bms(j,nfl,nll)
        ujp1  = vsi(j+1,nfl,nll,ni)/bms(j+1,nfl,nll)
        ur = .5*( uj +ujp1)
        ul = .5*( uj +ujm1)

        if (ur .ge. 0. .and. ul .ge. 0.) then
          a0 = -ul
          b0 =  ur
          c0 =  0.
        endif
        if (ur .le. 0. .and. ul .le. 0.) then
          a0 = 0.
          b0 = -ul
          c0 = ur
        endif
        if (ur .ge. 0. .and. ul .le. 0.) then
          a0 = 0.
          b0 = ur - ul
          c0 = 0.
        endif
        if (ur .le. 0. .and. ul .ge. 0.) then
          a0 = -ul
          b0 = 0.
          c0 = ur
        endif

        a(j) =  a0 * bms(j,nfl,nll) ** 2 / d22s(j,nfl,nll)
        b(j) = 1. / dt + loss(j) + 
     .         b0 * bms(j,nfl,nll) ** 2 / d22s(j,nfl,nll)
        c(j) = c0 * bms(j,nfl,nll) ** 2 / d22s(j,nfl,nll)
        d(j) = oldion(j) / dt + prod(j)

      enddo

! we will assume that they are determined by the production and loss
! at both ends of the field line

!     lower bc

      a(1) = 0.
      b(1) = 1.
      c(1) = 0.
      d(1) = 
     .  sqrt ( tdeni(1) * prod(1) / loss(1) ) + denmin

! upper bc 

      a(nz) = 0.
      b(nz) = 1.
      c(nz) = 0.
      d(nz) = 
     .  sqrt ( tdeni(nz) * prod(nz) / loss(nz) ) + denmin

      call rtryds ( a,b,c,d,tdeni,nz )
      
      return
      end



*******************************************
*******************************************

!            vsisolv

*******************************************
*******************************************

      subroutine vsisolv ( vii,vid,viold,snuj,nfl,nll,cs )

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc'

      dimension a(nz), b(nz), c(nz), d(nz)
      real vii(nz),vid(nz),viold(nz),snuj(nz),cs(nz)

! initialize

      do j = 1,nz
        a(j) = 0.
        b(j) = 0.
        c(j) = 0.
        d(j) = 0.
      enddo

      do j = 2,nz-1

        ujm1 = vii(j-1)
        uj   = vii(j)
        ujp1 = vii(j+1)
        ur = .25*( uj +ujp1)
        ul = .25*( uj +ujm1)

        if (ur .ge. 0. .and. ul .ge. 0.) then
          a0 = -ul
          b0 =  ur
          c0 =  0.
        endif
        if (ur .le. 0. .and. ul .le. 0.) then
          a0 = 0.
          b0 = -ul
          c0 = ur
        endif
        if (ur .ge. 0. .and. ul .le. 0.) then
          a0 = 0.
          b0 = ur - ul
          c0 = 0.
        endif
        if (ur .le. 0. .and. ul .ge. 0.) then
          a0 = -ul
          b0 = 0.
          c0 = ur
        endif

! anomalous drag to prevent ions from going supersonic

          delcs        = 0.1 * cs(j)
          alpha_drag_p = ( vii(j) - 0.9 * cs(j) )  / delcs
          alpha_drag_n = ( vii(j) + 0.9 * cs(j) )  / delcs
          anu_drag    = 0.5 * anu_drag0 * ( 1. + tanh(alpha_drag_p)) + 
     .                  0.5 * anu_drag0 * ( 1. - tanh(alpha_drag_n)) 

         
        a(j) = a0 / d22s(j,nfl,nll) * bms(j,nfl,nll)
        b(j) = 1/dt + snuj(j) + b0 / d22s(j,nfl,nll) *
     .                          bms(j,nfl,nll) + anu_drag
        c(j) = c0 / d22s(j,nfl,nll) * bms(j,nfl,nll)
        d(j) = viold(j)/dt + vid(j)

      enddo

! we will assume that the bc's are zero
! at both ends of the field line

! lower bc

      a(1) = 0.
      b(1) = 1.
      c(1) = 0.
      d(1) = 0.

! upper bc
 
      a(nz) = 0.
      b(nz) = 1.
      c(nz) = 0.
      d(nz) = 0.

      call rtryds(a,b,c,d,vii,nz)

      return
      end

*******************************************
*******************************************

!            tisolv

*******************************************
*******************************************

      subroutine tisolv(tti,tio,kap,s1,s2,s3,s4,s5,s6,s7,npt,nfl,nll)

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc'

      real a(nz),b(nz),c(nz),d(nz)
      real s1(nz),s2(nz),s3(nz),tti(nz),tio(nz),kap(nz)
      real s4(nz),s5(nz),s6(nz),s7(nz)

! initialize

      do j = 1,nz
        a(j) = 0.
        b(j) = 0.
        c(j) = 0.
        d(j) = 0.
      enddo


      do j = 2,nz-1
        ujm1 = bms(j-1,nfl,nll)*vsi(j-1,nfl,nll,npt)
        uj   = bms(j,nfl,nll)  *vsi(j,nfl,nll,npt)
        ujp1 = bms(j+1,nfl,nll)*vsi(j+1,nfl,nll,npt)
        ur = .5*( uj +ujp1)
        ul = .5*( uj +ujm1)

        if (ur .ge. 0. .and. ul .ge. 0.) then
          a0 = -ul
          b0 =  ur
          c0 =  0.
        endif
        if (ur .le. 0. .and. ul .le. 0.) then
          a0 = 0.
          b0 = -ul
          c0 = ur
        endif
        if (ur .ge. 0. .and. ul .le. 0.) then
          a0 = 0.
          b0 = ur - ul
          c0 = 0.
        endif
        if (ur .le. 0. .and. ul .ge. 0.) then
          a0 = -ul
          b0 = 0.
          c0 = ur
        endif

! anu_drag is velocity-dependent cooling term
! it limits ion cooling at high velocities due
! to 'compressional' heating of ions at magnetic apex

! not needed with anomalous drag term on ion velocity

c$$$          alpha_drag_p = ( vsi(j,nfl,nll,npt) - vsi0 )  / delta_vsi0  
c$$$          alpha_drag_n = ( vsi(j,nfl,nll,npt) + vsi0 )  / delta_vsi0 
c$$$          anu_drag    = 0.5 * anu_drag0 * ( 1. + tanh(alpha_drag_p)) + 
c$$$     .                  0.5 * anu_drag0 * ( 1. - tanh(alpha_drag_n)) 

        a(j) =     a0 / d22s(j,nfl,nll) 
     .         - ( bms(j,nfl,nll)**2 / deni(j,nfl,nll,npt) ) / 
     .                                 d22s(j,nfl,nll)
     .           *.5 * ( kap(j) + kap(j-1) ) / ds(j,nfl,nll)
         b(j) = 1. / dt + b0 / d22s(j,nfl,nll) 
     .         -.333333 * ( bms(j,nfl,nll) 
     .                     * (vsi(j+1,nfl,nll,npt) - 
     .                        vsi(j-1,nfl,nll,npt) ) 
     .                     + 5. * vsi(j,nfl,nll,npt) 
     .                          * (bms(j+1,nfl,nll) - 
     .                             bms(j-1,nfl,nll) ) )
     .         / d2s(j,nfl,nll) 
     .         +  ( bms(j,nfl,nll)**2 / deni(j,nfl,nll,npt) ) / 
     .                                  d22s(j,nfl,nll)
     .           *(.5* (kap(j+1) + kap(j) ) / ds(j+1,nfl,nll) 
     .         +.5 * (kap(j) + kap(j-1) ) / ds(j,nfl,nll)) 
     .         + s2(j) + s4(j) + s6(j) !+ anu_drag
         c(j) =     c0 / d22s(j,nfl,nll) 
     .         - ( bms(j,nfl,nll)**2 / deni(j,nfl,nll,npt) ) /
     .                                 d22s(j,nfl,nll)
     .           *.5 * (kap(j+1) + kap(j) ) / ds(j+1,nfl,nll)
        d(j) = tio(j)/dt + s1(j) + s3(j) + s5(j) + s7(j)

      enddo

! we will assume that the bc's are the neutral temperature
! at both ends of the field line

! lower bc

      a(1) = 0.
      b(1) = 1.
      c(1) = 0.
      d(1) = tn(1,nfl,nll)

! upper bc
 
      a(nz) = 0.
      b(nz) = 1.
      c(nz) = 0.
      d(nz) = tn(nz,nfl,nll)

      call rtryds ( a,b,c,d,tti,nz )


      return
      end

*******************************************
*******************************************

!            tesolv

*******************************************
*******************************************

      subroutine tesolv(tte,te_old,kap,s1,s2,s3,s4,s5,nfl,nll)

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc'

      dimension a(nz),b(nz),c(nz),d(nz)
      dimension s1(nz),s2(nz),s3(nz),s4(nz),s5(nz)
      real kap(nz),te_old(nz),tte(nz)

! initialize

      do j = 1,nz
        a(j) = 0.
        b(j) = 0.
        c(j) = 0.
        d(j) = 0.
      enddo

! note: ne used here is in a common block

      do j = 2,nz-1

        a(j) = - bms(j,nfl,nll)**2 / ne(j,nfl,nll) / d22s(j,nfl,nll)
     .         *.5 * ( kap(j) + kap(j-1) ) / ds(j,nfl,nll)
        b(j) = 1. / dt + bms(j,nfl,nll)**2 / ne(j,nfl,nll) / 
     .                                      d22s(j,nfl,nll)
     .        *(  .5 * (kap(j+1) + kap(j)   ) /ds(j+1,nfl,nll) 
     .           +.5 * (kap(j)   + kap(j-1) ) /ds(j,nfl,nll)   ) 
     .        + s2(j) + s3(j)
        c(j) = - bms(j,nfl,nll)**2 / ne(j,nfl,nll) /d22s(j,nfl,nll)
     .         *.5 * ( kap(j+1) + kap(j) )/ ds(j+1,nfl,nll)
        d(j) = te_old(j)/dt + s1(j) + s4(j) + s5(j)

       enddo

! we will assume that the bc's are the neutral temperature
! at both ends of the field line

! lower bc

      a(1) = 0.
      b(1) = 1.
      c(1) = 0.
      d(1) = tn(1,nfl,nll)

! upper bc
 
      a(nz) = 0.
      b(nz) = 1.
      c(nz) = 0.
      d(nz) = tn(nz,nfl,nll)

      call rtryds(a,b,c,d,tte,nz)

      return
      end

*******************************************
*******************************************

!            rtryds

*******************************************
*******************************************

      subroutine rtryds(a,b,c,d,x,n)
      
      include 'param3_mpi-1.98.inc'

! arrays a,b,c, and d may be used for stoage of alfa, beta and x
! in the actual call of this routine, but remember, whatever you
! use will be lost by the definition of of alfa and beta here.
! form,  a(k)*x(k-1) + b(k)*x(k) + c(k)*x(k+1) = d(k)

! i have modified the input sequence to the routine, but have left it
! otherwise intact.  we may  want to eventually change this (gj)

      dimension a(n),b(n),c(n),d(n),x(n)
      dimension alfa(nz),beta(nz)

      nm1=n-1

! apply the boundary condition at x(1)
! alfa(1) and beta(1) determined from b(1),c(1),d(1),a(1)=0.

      dst     = d(1)
      rb      = 1. / b(1)
      alfa(1) = -c(1) * rb
      beta(1) =   dst * rb

! calculate the alfas and betas of k on forward sweep

      do k=2,nm1
        ast     =  a(k)
        z       =  1. / ( b(k) + ast * alfa(k-1) )
        dst     =  d(k)
        alfa(k) = -c(k) * z
        beta(k) =  ( dst - ast * beta(k-1) ) * z
      enddo

! apply the boundary condition at x(n)
! x(n) determined from a(n),b(n),d(n),c(n)=0.

      x(n) = ( d(n) - a(n) *beta(nm1) ) / ( b(n) + a(n) * alfa(nm1) )

! calculate x of k from the alfas and betas on backward sweep

      do i=2,n
        k    = n + 1 - i
        x(k) = x(k+1) * alfa(k) + beta(k)
      enddo

      return
      end

*******************************************
*******************************************

!            msistim

*******************************************
*******************************************


      subroutine msistim ( iyr,iday,hr,glong,iyd,secut )

! msistim calculates time parameters for the 
! nrlmsise00 neutral atmosphere model.

! the arguments are defined as follows:

!   iyr    the julian year
!   iday   the day of the year
!   hr     the local time in hours
!   glong  the geocentric longitude in degrees east
!   iyd    the year and day in the form yydd
!   secut  the universal time in seconds

      iyd    = 1000 * mod(iyr,100) + iday
      hrut   = hr - glong /15.

      do while ( hrut .lt. 0.  )
        hrut = hrut + 24.
      enddo

      do while ( hrut. ge. 24. )
        hrut = hrut - 24.
      enddo

      secut  = hrut * 3600.

      return
      end


*******************************************
*******************************************

!            zenith

*******************************************
*******************************************

       subroutine zenith (hrut,nfl,nll)

       include 'param3_mpi-1.98.inc'
       include 'com3_mpi-1.98.inc' 

! geometric variables


! bdec: magnetic declination angle
! sdec: solar zenith angle
! cx:  cos of the zenith angle

       do i = 1,nz
         hrl   = mod(hrut + glons(i,nfl,nll) / 15.,24.)
         if ( lcr ) then
!            glons0(i,j,nll) = glons(i,j,nll)+15.*hrinit
!            hrl = mod(hrinit+glons0(i,nfl,nll) / 15.,24.)
            hrl = mod(glons(i,nfl,nll) / 15.,24.)
         endif
         sdec          = rtod * asin (  sin (2.*pie*(day-dayve)/sidyr)
     .                                * sin (solinc/rtod)             )
         cossdec       = cos ( po180 * sdec )
         sinsdec       = sin ( po180 * sdec )
         clat          = cos ( po180 * glats(i,nfl,nll) )
         slat          = sin ( po180 * glats(i,nfl,nll) )
         cx(i,nfl,nll) =   slat * sinsdec 
     .                   - clat * cossdec * cos ( 15.0*po180*hrl )
! MS: Since we will be taking acos of this value in photprod, make
! sure that the absolute value does not minutely exceed 1 because of
! round-off error.

        if (abs(abs(cx(i,nfl,nll))-1.) .lt. 1.e-6) 
     .     cx(i,nfl,nll) = sign(1.,cx(i,nfl,nll))
       enddo

       return
       end


*******************************************
*******************************************

!            mag_zenith

*******************************************
*******************************************

       subroutine mag_zenith (hrut,nfl,nll)

       include 'param3_mpi-1.98.inc'
       include 'com3_mpi-1.98.inc' 

! base zenith on fixed magnetic grid

! geometric variables

! bdec: magnetic declination angle
! sdec: solar zenith angle
! cx:  cos of the zenith angle

       do i = 1,nz
         blonsij = blons(i,nfl,nll)
!!!         glonsij = glons0(i,nfl,nll) -
!!!     .             (hrut - hrinit) * 15.
!!!         hrl   = mod(hrut + blonsij / 15.,24.)
         hrl   = mod(blonsij / 15.,24.)
!         hrl   = mod(hrut + glons(i,nfl,nll) / 15.,24.)
!         call magdec ( glats(i,nfl,nll),glons(i,nfl,nll),bdec )
!         cosbdec(i,nfl,nll) = cos ( po180 * bdec )
!         sinbdec(i,nfl,nll) = sin ( po180 * bdec )
         sdec          = rtod * asin (  sin (2.*pie*(day-dayve)/sidyr)
     .                                * sin (solinc/rtod)             )
         cossdec       = cos ( po180 * sdec )
         sinsdec       = sin ( po180 * sdec )
         clat          = cos ( po180 * glats(i,nfl,nll) )
         slat          = sin ( po180 * glats(i,nfl,nll) )
         cx(i,nfl,nll) =   slat * sinsdec 
     .                   - clat * cossdec * cos ( 15.0*po180*hrl )
! MS: Since we will be taking acos of this value in photprod, make
! sure that the absolute value does not minutely exceed 1 because of
! round-off error.

        if (abs(abs(cx(i,nfl,nll))-1.) .lt. 1.e-6) 
     .     cx(i,nfl,nll) = sign(1.,cx(i,nfl,nll))
       enddo

       return
       end

*******************************************
*******************************************

!      magdec

*******************************************
*******************************************

! this routine uses a table to calculate the magnetic declination
! (decpt) for use with a neutral wind model. the table was constructed 
! using " the earth's magnetic field" by robert merrill and
! michael mcelhinny, academic press p 18. (uah library # qc816.m47) 
! the declination at the specified location (degrees) is obtained
! by bilinear interpolation of the table values
! written by phil richards, august 1995 
! (modified by j huba june 1998)
! this code is not accurate near the magnetic poles
! glat, glong, and decpt are in degrees

       subroutine magdec ( glat,glong,decpt )
       implicit none

! nlong = # longs, nlat =  # lats. i,j,k are array indices

       integer nlong,nlat,i,j,k
       parameter (nlong = 13,nlat = 11) 
       real glat,glong,decpt
       real along(nlong),alat(nlat),decls(13,11),rlat,rlong

! fill up the longitude, latitude, and magnetic declination arrays

       data along /0,30,60,90,120,150,180,210,240,270,300,330,360/
       data alat  /-90,-80,-60,-40,-20,0,20,40,60,80,90/

! the declinations are loaded in latitude slices. the first nlong
! values are for the first latitude in the alat array

       data decls / -20,-45,-67,-95,-125,170,110,80,60,40,20,0,-20,
     .  -20,-40,-63,-82,-120,140,82,65,52,35,16,-2,-20, -20,-35,-55,-70,
     .  -60,35,45,40,38,30,10,-7.5,-20, -27,-28,-42,-35,-5,13,20,21,21,
     .  20,0,-21,-27, -21,-15,-22,-10,2,8,13,13,12,11,-7,-25,-21, -10,
     .  -2,-4,-3.5,2,6,12,9,8,6,-10,-20,-10, -5,1,0,-2,-2,1,10,12,11,5,
     .  -14,-16,-5, -5,3,5,2,-7,-5,7,17,17,2,-20,-16,-5, -7,7,16,10,-12,
     .  -11,6,25,30,-7.5,-36,-15,-7, -10,10,25,20,-11,-11,7,30,42,-40,
     .  -52,-33,-10, -12,10,27,25,-7,-7,7,32,50,-70,-60,-35,-12/

! find the index in the longitude array

       i = 1
       do while ( glong .ge. along(i) .and. i .lt. nlong )
         i = i + 1
       enddo
       j = i - 1

! find the index in the latitude array

       i = 1
       do while ( glat .ge. alat(i) .and. i .lt. nlat )
         i = i + 1
       enddo
       k = i - 1

! now do the bilinear fit (press et al. num. rec. page 96)

       rlat  = ( glat  - alat(k)  ) / ( alat (k+1) - alat(k)  )
       rlong = ( glong - along(j) ) / ( along(j+1) - along(j) )
       decpt = ( 1 - rlat ) * ( 1 - rlong ) * decls(j,k) 
     .         + rlong * ( 1 - rlat ) * decls(j+1,k) 
     .         + rlong * rlat  * decls(j+1,k+1) 
     .         + ( 1 - rlong ) * rlat * decls(j,k+1)

       return
       end


*******************************************
*******************************************

!            f1026

*******************************************
*******************************************

! subroutine to calculate the nighttime flux of
! lyman beta (1026) (note: line = 1)

      subroutine sf1026 ( f,line,nfl,nll )

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc' 

      real f(nz,nf,nl,91)

      imax = 1

! determine f for the 4 known values of theta

      do i = 1,nz
        if ( alts(i,nfl,nll) .lt. zaltnt(line,1) ) then
          do k = 1,4
            f( i,nfl,nll,int(thetant(line,k))+1-90 ) = 1.
          enddo
        elseif ( zaltnt(line,1) .le. alts(i,nfl,nll) .and.
     .           alts(i,nfl,nll) .le. zaltnt(line,2)       ) then
          f( i,nfl,nll,int(thetant(line,1))+1-90 ) =
     .       1.4e8 * tanh ( (alts(i,nfl,nll) - 90.) / 50. )
          f( i,nfl,nll,int(thetant(line,2))+1-90 ) =
     .       3.8e7 * tanh ( (alts(i,nfl,nll) - 90.) / 50. )
          f( i,nfl,nll,int(thetant(line,3))+1-90 ) =
     .       1.4e7 * tanh ( (alts(i,nfl,nll) - 93.) / 55. )
          f( i,nfl,nll,int(thetant(line,4))+1-90 ) =
     .       9.2e6 * tanh ( (alts(i,nfl,nll) - 94.) / 55. )
          imax = i
        else
          do k = 1,4
            f( i,nfl,nll,   int(thetant(line,k))+1-90 ) = 
     .      f( imax,nfl,nll,int(thetant(line,k))+1-90 )
          enddo
        endif           
      enddo

      do k = 1,4
        do i = 1,nz
          f( i,nfl,nll,int(thetant(line,k))+1-90 ) = 
     .    amax1 ( 1., f( i,nfl,nll,int(thetant(line,k))+1-90 ) ) 
        enddo
      enddo

! now interpolate to all valuse of theta (90 - 180)
 
      do k = 1,91
        k90 = 90 + k - 1
        ji  = 1
        ki  = int(thetant(line,1))
        do j = 1,4
          if ( k90 .gt. int(thetant(line,j)) ) then
            ji = j
            ki = int(thetant(line,ji))
          endif
        enddo
        jip1 = ji + 1
        kip1 = int(thetant(line,jip1))
        delk = float (   int(thetant(line,jip1)) 
     .                 - int(thetant(line,ji  )) )
        do i = 1,nz
          flog =   alog10(f(i,nfl,nll,ki+1-90)) 
     .           + (k90 - ki) / delk 
     .                        * (  alog10(f(i,nfl,nll,kip1+1-90)) 
     .                           - alog10(f(i,nfl,nll,ki  +1-90)) ) 
          f(i,nfl,nll,k) = 10 ** flog
        enddo
      enddo

      return
      end

*******************************************
*******************************************

!            f584

*******************************************
*******************************************

! subroutine to calculate the nighttime flux of
! he i (584) (note: line = 2)

      subroutine sf584 ( f,line,nfl,nll )

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc' 

      real f(nz,nf,nl,91)

      imax = 1

! determine f for the 4 known values of theta

      do i = 1,nz
        if ( alts(i,nfl,nll) .lt. zaltnt(line,1) ) then
          do k = 1,4
            f( i,nfl,nll,int(thetant(line,k))+1-90 ) = 1.
          enddo
        elseif ( zaltnt(line,1) .le. alts(i,nfl,nll) .and.
     .           alts(i,nfl,nll) .le. zaltnt(line,2)       ) then
          f( i,nfl,nll,int(thetant(line,1))+1-90 ) =
     .       1.85e5 * ( alts(i,nfl,nll) - 170. ) ** 1.20        
          f( i,nfl,nll,int(thetant(line,2))+1-90 ) =
     .       2.60e4 * ( alts(i,nfl,nll) - 170. ) ** 1.25        
          f( i,nfl,nll,int(thetant(line,3))+1-90 ) =
     .       2.60e3 * ( alts(i,nfl,nll) - 170. ) ** 1.20        
          f( i,nfl,nll,int(thetant(line,4))+1-90 ) =
     .       2.60e2 * ( alts(i,nfl,nll) - 170. ) ** 1.20       
          imax = i
        else
          do k = 1,4
            f( i   ,nfl,nll,int(thetant(line,k))+1-90 ) = 
     .      f( imax,nfl,nll,int(thetant(line,k))+1-90 )
          enddo
        endif           
      enddo

      do k = 1,4
        do i = 1,nz
          f( i,nfl,nll,int(thetant(line,k))+1-90 ) = 
     .    amax1 ( 1., f( i,nfl,nll,int(thetant(line,k))+1-90 ) ) 
        enddo
      enddo

! now interpolate to all valuse of theta (90 - 180)
! set f(i,nfl,nll,theta=180) = 1. 

      do k = 1,91
        k90 = 90 + k - 1
        ji  = 1
        ki  = int(thetant(line,1))
        do j = 1,4
          if ( k90 .gt. int(thetant(line,j)) ) then
            ji = j
            ki = int(thetant(line,ji))
          endif
        enddo
        if ( ji .ne. 4 ) then
          jip1 = ji + 1
          kip1 = int(thetant(line,jip1))
          delk = float (   int(thetant(line,jip1)) 
     .                   - int(thetant(line,ji  )) )
          do i = 1,nz
            flog =   alog10(f(i,nfl,nll,ki+1-90)) 
     .             + (k90 - ki) / delk 
     .                          * (  alog10(f(i,nfl,nll,kip1+1-90)) 
     .                             - alog10(f(i,nfl,nll,ki  +1-90)) ) 
            f(i,nfl,nll,k) = 10 ** flog
          enddo
        else
          delk = float (   180 
     .                   - int(thetant(line,ji  )) )
          do i = 1,nz
            flog =   alog10(f(i,nfl,nll,ki+1-90)) 
     .             + (k90 - ki) / delk 
     .                          * (  alog10(1.) 
     .                             - alog10(f(i,nfl,nll,ki  +1-90)) ) 
            f(i,nfl,nll,k) = 10 ** flog
          enddo
        endif
      enddo

      return
      end

*******************************************
*******************************************

!            f304

*******************************************
*******************************************

! subroutine to calculate the nighttime flux of
! he ii (304) (note: line = 3)

      subroutine sf304 ( f,line,nfl,nll )

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc' 

      real f(nz,nf,nl,91)

      imax = 1

! determine f for the 4 known values of theta

      do i = 1,nz
        if ( alts(i,nfl,nll) .lt. zaltnt(line,1) ) then
          do k = 1,4
            f( i,nfl,nll,int(thetant(line,k))+1-90 ) = 1.
          enddo
        elseif ( zaltnt(line,1) .le. alts(i,nfl,nll) .and.
     .           alts(i,nfl,nll) .le. zaltnt(line,2)       ) then
          f( i,nfl,nll,int(thetant(line,1))+1-90 ) =
     .       3.8e6 * tanh ( (alts(i,nfl,nll) - 138.) / 80. )
          f( i,nfl,nll,int(thetant(line,2))+1-90 ) =
     .       3.0e6 * tanh ( (alts(i,nfl,nll) - 138.) / 80. )
          f( i,nfl,nll,int(thetant(line,3))+1-90 ) =
     .       2.5e6 * tanh ( (alts(i,nfl,nll) - 138.) / 80. )
          f( i,nfl,nll,int(thetant(line,4))+1-90 ) =
     .       2.5e6 * tanh ( (alts(i,nfl,nll) - 138.) / 80. )
          imax = i
        else
          do k = 1,4
            f( i,   nfl,nll,int(thetant(line,k))+1-90 ) = 
     .      f( imax,nfl,nll,int(thetant(line,k))+1-90 )
          enddo
        endif           
      enddo

      do k = 1,4
        do i = 1,nz
          f( i,nfl,nll,int(thetant(line,k))+1-90 ) = 
     .    amax1 ( 1., f( i,nfl,nll,int(thetant(line,k))+1-90 ) ) 
        enddo
      enddo

! now interpolate to all valuse of theta (90 - 180)
! set f(i,nfl,nll,theta=180) = 1. 

      do k = 1,91
        k90 = 90 + k - 1
        ji  = 1
        ki  = int(thetant(line,1))
        do j = 1,4
          if ( k90 .gt. int(thetant(line,j)) ) then
            ji = j
            ki = int(thetant(line,ji))
          endif
        enddo
        if ( ji .ne. 4 ) then
          jip1 = ji + 1
          kip1 = int(thetant(line,jip1))
          delk = float (   int(thetant(line,jip1)) 
     .                   - int(thetant(line,ji  )) )
          do i = 1,nz
            flog =   alog10(f(i,nfl,nll,ki+1-90)) 
     .             + (k90 - ki) / delk 
     .                          * (  alog10(f(i,nfl,nll,kip1+1-90)) 
     .                             - alog10(f(i,nfl,nll,ki  +1-90)) ) 
            f(i,nfl,nll,k) = 10 ** flog
          enddo
        else
          delk = float (   180 
     .                   - int(thetant(line,ji  )) )
          do i = 1,nz
            flog =   alog10(f(i,nfl,nll,ki+1-90)) 
     .             + (k90 - ki) / delk 
     .                          * (  alog10(1.) 
     .                             - alog10(f(i,nfl,nll,ki  +1-90)) ) 
            f(i,nfl,nll,k) = 10 ** flog
          enddo
        endif
      enddo

      return
      end

*******************************************
*******************************************

!            f1216

*******************************************
*******************************************

!     subroutine to calculate the nighttime flux of
!     lyman alpha (1216) (note: line = 4)

      subroutine sf1216 ( f,line,nfl,nll )

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc' 

      real f(nz,nf,nl,91)

      imax = 1

! determine f for the 4 known values of theta

      do i = 1,nz
        if ( alts(i,nfl,nll) .lt. zaltnt(line,1) ) then
          do k = 1,4
            f( i,nfl,nll,int(thetant(line,k))+1-90 ) = 1.
          enddo
        elseif ( zaltnt(line,1) .le. alts(i,nfl,nll) .and.
     .           alts(i,nfl,nll) .le. zaltnt(line,2)       ) then
          f( i,nfl,nll,int(thetant(line,1))+1-90 ) =
     .       1.2e10 * tanh ( (alts(i,nfl,nll) - 80.) / 50. ) + 3.e9
          f( i,nfl,nll,int(thetant(line,2))+1-90 ) =
     .       4.0e9  * tanh ( (alts(i,nfl,nll) - 80.) / 50. ) + 1.e9
          f( i,nfl,nll,int(thetant(line,3))+1-90 ) =
     .       2.0e9  * tanh ( (alts(i,nfl,nll) - 65.) / 50. ) + 1.e8
          f( i,nfl,nll,int(thetant(line,4))+1-90 ) =
     .       1.5e9  * tanh ( (alts(i,nfl,nll) - 75.) / 50. ) + 1.e8
          imax = i
        else
          do k = 1,4
            f( i,   nfl,nll,int(thetant(line,k))+1-90 ) = 
     .      f( imax,nfl,nll,int(thetant(line,k))+1-90 )
          enddo
        endif           
      enddo

      do k = 1,4
        do i = 1,nz
          f( i,nfl,nll,int(thetant(line,k))+1-90 ) = 
     .    amax1 ( 1., f( i,nfl,nll,int(thetant(line,k))+1-90 ) ) 
        enddo
      enddo

! now interpolate to all valuse of theta (90 - 180)
 
      do k = 1,91
        k90 = 90 + k - 1
        ji  = 1
        ki  = int(thetant(line,1))
        do j = 1,4
          if ( k90 .gt. int(thetant(line,j)) ) then
            ji = j
            ki = int(thetant(line,ji))
          endif
        enddo
        jip1 = ji + 1
        kip1 = int(thetant(line,jip1))
        delk = float (   int(thetant(line,jip1)) 
     .                 - int(thetant(line,ji  )) )
        do i = 1,nz
          flog =   alog10(f(i,nfl,nll,ki+1-90)) 
     .           + (k90 - ki) / delk 
     .                        * (  alog10(f(i,nfl,nll,kip1+1-90)) 
     .                           - alog10(f(i,nfl,nll,ki  +1-90)) ) 
          f(i,nfl,nll,k) = 10 ** flog
        enddo
      enddo

      return
      end

*******************************************
*******************************************

!             open_u

*******************************************
*******************************************

      subroutine open_u
      character(len=14),parameter :: plotdir='IM/plotsSAMI3/'
! open output files (unformatted, except time.dat)

      open ( unit=70, file=plotdir//'time.dat' ,form='formatted'   )  
      open ( unit=71, file=plotdir//'deniu.dat',form='unformatted' )
      open ( unit=72, file=plotdir//'tiu.dat'  ,form='unformatted' )
      open ( unit=73, file=plotdir//'vsiu.dat' ,form='unformatted' )
      open ( unit=75, file=plotdir//'teu.dat'  ,form='unformatted' )
      open ( unit=78, file=plotdir//'vnu.dat'  ,form='unformatted' )
      open ( unit=92, file=plotdir//'dennu.dat',form='unformatted' )
      open ( unit=93, file=plotdir//'hipcu.dat',form='unformatted' )
      open ( unit=94, file=plotdir//'hihcu.dat',form='unformatted' )
      open ( unit=95, file=plotdir//'sigmapu.dat',form='unformatted' )
      open ( unit=96, file=plotdir//'sigmahu.dat',form='unformatted' )
      open ( unit=97, file=plotdir//'sigmapicu.dat',form='unformatted' )
      open ( unit=98, file=plotdir//'sigmahicu.dat',form='unformatted' )
      open ( unit=711, file=plotdir//'deniu1.dat',form='unformatted' )
      open ( unit=712, file=plotdir//'deniu2.dat',form='unformatted' )
      open ( unit=713, file=plotdir//'deniu3.dat',form='unformatted' )
      open ( unit=714, file=plotdir//'deniu4.dat',form='unformatted' )
      open ( unit=715, file=plotdir//'deniu5.dat',form='unformatted' )
      open ( unit=716, file=plotdir//'deniu6.dat',form='unformatted' )
      open ( unit=717, file=plotdir//'deniu7.dat',form='unformatted' )
      open ( unit=1718, file=plotdir//'deneu.dat',form='unformatted' )
      open ( unit=811, file=plotdir//'tiu1.dat',form='unformatted' )
      open ( unit=812, file=plotdir//'tiu2.dat',form='unformatted' )
      open ( unit=813, file=plotdir//'tiu3.dat',form='unformatted' )
      open ( unit=814, file=plotdir//'tiu4.dat',form='unformatted' )
      open ( unit=815, file=plotdir//'tiu5.dat',form='unformatted' )
      open ( unit=816, file=plotdir//'tiu6.dat',form='unformatted' )
      open ( unit=817, file=plotdir//'tiu7.dat',form='unformatted' )
      open ( unit=911, file=plotdir//'vsiu1.dat',form='unformatted' )
      open ( unit=912, file=plotdir//'vsiu2.dat',form='unformatted' )
      open ( unit=913, file=plotdir//'vsiu3.dat',form='unformatted' )
      open ( unit=914, file=plotdir//'vsiu4.dat',form='unformatted' )
      open ( unit=915, file=plotdir//'vsiu5.dat',form='unformatted' )
      open ( unit=916, file=plotdir//'vsiu6.dat',form='unformatted' )
      open ( unit=917, file=plotdir//'vsiu7.dat',form='unformatted' )
      open ( unit=1711, file=plotdir//'dennu1.dat',form='unformatted' )
      open ( unit=1712, file=plotdir//'dennu2.dat',form='unformatted' )
      open ( unit=1713, file=plotdir//'dennu3.dat',form='unformatted' )
      open ( unit=1714, file=plotdir//'dennu4.dat',form='unformatted' )
      open ( unit=1715, file=plotdir//'dennu5.dat',form='unformatted' )
      open ( unit=569,  file=plotdir//'rhsegv.dat',form='unformatted' )


      open ( unit=196, file=plotdir//'vnqu.dat',form='unformatted' )
      open ( unit=197, file=plotdir//'vnpu.dat',form='unformatted' )
      open ( unit=198, file=plotdir//'vnphiu.dat',form='unformatted' )
      open ( unit=201, file=plotdir//'jpu.dat',form='unformatted' )
      open ( unit=202, file=plotdir//'jphiu.dat',form='unformatted' )

      open ( unit=491, file=plotdir//'hipcpu.dat',form='unformatted' )
      open ( unit=492, file=plotdir//'hipcphiu.dat',form='unformatted' )
      open ( unit=493, file=plotdir//'hihcmu.dat',form='unformatted' )
      open ( unit=494, file=plotdir//'hidpvu.dat',form='unformatted' )
      open ( unit=495, file=plotdir//'hidphivu.dat',form='unformatted' )
      open ( unit=496, file=plotdir//'hidpgu.dat',form='unformatted' )
      open ( unit=497, file=plotdir//'hidphigu.dat',form='unformatted' )
      open ( unit=498, file=plotdir//'phiu.dat'    ,form='unformatted' )

! diagnostic files (unformatted)

      open ( unit=81, file=plotdir//'t1u.dat'  ,form='unformatted' )
      open ( unit=82, file=plotdir//'t2u.dat'  ,form='unformatted' )
      open ( unit=83, file=plotdir//'t3u.dat'  ,form='unformatted' )
      open ( unit=84, file=plotdir//'u1u.dat'  ,form='unformatted' )
      open ( unit=85, file=plotdir//'u2u.dat'  ,form='unformatted' )
      open ( unit=86, file=plotdir//'u3u.dat'  ,form='unformatted' )
      open ( unit=87, file=plotdir//'u4u.dat'  ,form='unformatted' )
      open ( unit=88, file=plotdir//'u5u.dat'  ,form='unformatted' )
      open ( unit=384, file=plotdir//'u1pu.dat'  ,form='unformatted' )
      open ( unit=385, file=plotdir//'u2su.dat'  ,form='unformatted' )
      open ( unit=386, file=plotdir//'u3hu.dat'  ,form='unformatted' )

      return
      end

*******************************************
*******************************************

!             open_f

*******************************************
*******************************************

      subroutine open_f
      character(len=14),parameter :: plotdir='IM/plotsSAMI3/'
! open output files (formatted)

      open ( unit=70, file=plotdir//'time.dat'      ,form='formatted' )  
      open ( unit=71, file=plotdir//'denif.dat'     ,form='formatted' )
      open ( unit=72, file=plotdir//'tif.dat'       ,form='formatted' )
      open ( unit=73, file=plotdir//'vsif.dat'      ,form='formatted' )
      open ( unit=75, file=plotdir//'tef.dat'       ,form='formatted' )
      open ( unit=78, file=plotdir//'vnf.dat'       ,form='formatted' )
      open ( unit=92, file=plotdir//'dennf.dat'     ,form='formatted' )

! diagnostic files (formatted)

      open ( unit=81, file=plotdir//'t1f.dat'  ,form='formatted' )
      open ( unit=82, file=plotdir//'t2f.dat'  ,form='formatted' )
      open ( unit=83, file=plotdir//'t3f.dat'  ,form='formatted' )
      open ( unit=84, file=plotdir//'u1f.dat'  ,form='formatted' )
      open ( unit=85, file=plotdir//'u2f.dat'  ,form='formatted' )
      open ( unit=86, file=plotdir//'u3f.dat'  ,form='formatted' )
      open ( unit=87, file=plotdir//'u4f.dat'  ,form='formatted' )
      open ( unit=88, file=plotdir//'u5f.dat'  ,form='formatted' )

      return
      end


*******************************************
*******************************************

!             output

*******************************************
*******************************************

!       subroutine output ( hr,ntm,istep,phi )
       subroutine output ( hr,ntm,istep,phi,denit,dennt,vsit,sumvsit,
     &                     tit,ut,vt,vpit,net,tet,tnt,u1t,
     &                     u2t,u3t,u4t,vnqt,vnpt,vnphit,jpt,jphit,
     &                     u1pt,u2st,u3ht,sigmapict,sigmahict,
     &                     sigmapt,sigmaht )

       include 'param3_mpi-1.98.inc'
       include 'com3_mpi-1.98.inc'

       real denit(nz,nf,nlt,nion)
       real dennt(nz,nf,nlt,nion)
       real vsit(nz,nf,nlt,nion)
       real sumvsit(nz,nf,nlt,nion)
       real tet(nz,nf,nlt),tit(nz,nf,nlt,nion),tnt(nz,nf,nlt)
       real ut(nz,nf,nlt),vt(nz,nf,nlt),vpit(nz,nf,nlt)
       real u1t(nz,nf,nlt),u2t(nz,nf,nlt),u3t(nz,nf,nlt),u4t(nz,nf,nlt)
       real vnqt(nz,nf,nlt),vnpt(nz,nf,nlt),vnphit(nz,nf,nlt)
       real jpt(nz,nf,nlt),jphit(nz,nf,nlt)
       real phi(nnx,nny)
       real u1pt(nz,nf,nlt),u2st(nz,nf,nlt),u3ht(nz,nf,nlt)
       real sigmapict(nz,nf,nlt),sigmahict(nz,nf,nlt)
       real sigmapt(nz,nf,nlt),sigmaht(nz,nf,nlt)

       real, dimension(:,:,:), allocatable :: denit1,denit2,denit3,
     .      denit4,denit5,denit6,denit7,denet
       real, dimension(:,:,:), allocatable :: dennt1,dennt2,dennt3,
     .      dennt4,dennt5,dennt6,dennt7
       real, dimension(:,:,:), allocatable :: tit1,tit2,tit3,
     .      tit4,tit5,tit6,tit7
       real, dimension(:,:,:), allocatable ::vsit1,vsit2,vsit3,
     .      vsit4,vsit5,vsit6,vsit7


       hr24   = mod (hr,24.)
       totsec = hr24 * 3600.
       thr    = totsec / 3600.
       nthr   = int(thr)
       tmin   = ( thr - nthr ) * 60.
       ntmin  = int(mod(tmin,60.))
       tsec   = ( tmin - ntmin ) * 60.
       ntsec  = int(tsec)

!       print *,'istep = ',istep,' ntm = ',ntm
!       print *,' hr = ',hr,' dt = ',dt

       write (70,100) ntm,nthr,ntmin,ntsec,hr

       allocate 
     .       (denit1(nz,nf,nlt),denit2(nz,nf,nlt),denit3(nz,nf,nlt),
     .        denit4(nz,nf,nlt),denit5(nz,nf,nlt),denit6(nz,nf,nlt),
     .        denit7(nz,nf,nlt),denet(nz,nf,nlt))


       do k = 1,nlt
         do j = 1,nf
           do i = 1,nz
             denit1(i,j,k) = denit(i,j,k,1)
             denit2(i,j,k) = denit(i,j,k,2)
             denit3(i,j,k) = denit(i,j,k,3)
             denit4(i,j,k) = denit(i,j,k,4)
             denit5(i,j,k) = denit(i,j,k,5)
             denit6(i,j,k) = denit(i,j,k,6)
             denit7(i,j,k) = denit(i,j,k,7)
           enddo
         enddo
       enddo

       open (144,file='denit_cg.rst',form='unformatted')
       write(144) denit1,denit2,denit3,denit4,denit5,denit6,denit7
       close(144)

       do k = 1,nlt
         do j = 1,nf
           do i = 1,nz
             denet(i,j,k) = 0
             do ni = nion1,nion2
               denet(i,j,k) = denet(i,j,k) + denit(i,j,k,ni)
             enddo
           enddo
         enddo
       enddo



       if ( .not. fmtout ) then
         write(1718) denet
         write(711) denit1
         write(712) denit2
         write(713) denit3
         write(714) denit4
         write(715) denit5
!         write(716) denit6
!         write(717) denit7
       endif

       deallocate (denit1,denit2,denit3,
     .           denit4,denit5,denit6,denit7,denet)


!       print *,'allocate denn'
       allocate
     .       (dennt1(nz,nf,nlt),dennt2(nz,nf,nlt),dennt3(nz,nf,nlt),
     .        dennt4(nz,nf,nlt),dennt5(nz,nf,nlt),dennt6(nz,nf,nlt),
     .        dennt7(nz,nf,nlt))

       do k = 1,nlt
         do j = 1,nf
           do i = 1,nz
             dennt1(i,j,k) = dennt(i,j,k,1)
             dennt2(i,j,k) = dennt(i,j,k,2)
             dennt3(i,j,k) = dennt(i,j,k,3)
             dennt4(i,j,k) = dennt(i,j,k,4)
             dennt5(i,j,k) = dennt(i,j,k,5)
             dennt6(i,j,k) = dennt(i,j,k,6)
             dennt7(i,j,k) = dennt(i,j,k,7)
           enddo
         enddo
       enddo


       if ( .not. fmtout ) then
         write(1711) dennt1
         write(1712) dennt2
         write(1713) dennt3
         write(1714) dennt4
         write(1715) dennt5
       endif

       deallocate (dennt1,dennt2,dennt3,
     .                  dennt4,dennt5,dennt6,dennt7)
!       print *,'deallocate denn'

       allocate 
     .       (tit1(nz,nf,nlt),tit2(nz,nf,nlt),tit3(nz,nf,nlt),
     .        tit4(nz,nf,nlt),tit5(nz,nf,nlt),tit6(nz,nf,nlt),
     .        tit7(nz,nf,nlt))

       do k = 1,nlt
         do j = 1,nf
           do i = 1,nz
             tit1(i,j,k) = tit(i,j,k,1)
             tit2(i,j,k) = tit(i,j,k,2)
             tit3(i,j,k) = tit(i,j,k,3)
             tit4(i,j,k) = tit(i,j,k,4)
             tit5(i,j,k) = tit(i,j,k,5)
             tit6(i,j,k) = tit(i,j,k,6)
             tit7(i,j,k) = tit(i,j,k,7)
           enddo
         enddo
       enddo

       if ( .not. fmtout ) then
         write(811) tit1
         write(812) tit2
         write(813) tit3
!         write(814) tit4
!         write(815) tit5
!         write(816) tit6
!         write(817) tit7

         write(75) tet

       endif
       deallocate (tit1,tit2,tit3,tit4,tit5,tit6,tit7)

       allocate
     .      (vsit1(nz,nf,nlt),vsit2(nz,nf,nlt),vsit3(nz,nf,nlt),
     .       vsit4(nz,nf,nlt),vsit5(nz,nf,nlt),vsit6(nz,nf,nlt),
     .       vsit7(nz,nf,nlt))


       do k = 1,nlt
         do j = 1,nf
           do i = 1,nz
             vsit1(i,j,k) = vsit(i,j,k,1)
             vsit2(i,j,k) = vsit(i,j,k,2)
             vsit3(i,j,k) = vsit(i,j,k,3)
             vsit4(i,j,k) = vsit(i,j,k,4)
             vsit5(i,j,k) = vsit(i,j,k,5)
             vsit6(i,j,k) = vsit(i,j,k,6)
             vsit7(i,j,k) = vsit(i,j,k,7)
           enddo
         enddo
       enddo

       if ( .not. fmtout ) then

         write(911) vsit1
         write(912) vsit2
         write(913) vsit3
!         write(914) vsit4
         write(915) vsit5
!         write(916) vsit6
!         write(917) vsit7

       endif

       deallocate (vsit1,vsit2,vsit3,
     .             vsit4,vsit5,vsit6,vsit7)

       if ( .not. fmtout ) then

!         write(78) vnt
!         write(81) t1t
!         write(82) t2t
!         write(83) t3t
         write(384) u1pt
         write(385) u2st
         write(386) u3ht
         write(84) u1t
         write(85) u2t
         write(86) u3t
         write(87) u4t
!         write(88) u5t
!         write(93) hipct
!         write(94) hihct
         write(95) sigmapt
         write(96) sigmaht
         write(97) sigmapict
         write(98) sigmahict

         write(196) vnqt
         write(197) vnpt
         write(198) vnphit
         write(201) jpt
         write(202) jphit

         write(491) hipcpt
         write(492) hipcphit
         write(493) hihcmt
         write(494) hidpvt
         write(495) hidphivt
         write(496) hidpgt
         write(497) hidphigt
         write(498) phi

!         write(569) cxxe,cyye,cxe,cye,rhse
!         write(569) rhseg,rhsev
!        close(569)

       endif

 100   format(1x,4i6,1p1e14.4)
 101   format(1x,1p10e16.6)

       return
       end


*******************************************
*******************************************

!             courant

*******************************************
*******************************************

       subroutine courant

       include 'param3_mpi-1.98.inc'
       include 'com3_mpi-1.98.inc'

! parallel motion

       dtnew = 1.e6
       do k = nion1,nion2
         do l = 1,nl
           do j = 1,nf
             do i = 1,nz
               dt1 = dels(i,j,l) / amax1(1.,abs(vsi(i,j,l,k)))
               if ( dt1 .le. dtnew ) then
                 dtnew = dt1
                 i0    = i
                 j0    = j
                 k0    = k
                 l0    = l
               endif
             enddo
           enddo
         enddo
       enddo
!         print *,'parallel',taskid,dtnew,i0,j0,l0,k0

! perpendicular motion

       do k = 1,nl
         do j = 1,nf
           do i = 1,nz
             dts = xdels(i,j,k) / amax1(1.,abs(vexbs(i,j,k)))
             dtp = xdelp(i,j,k) / amax1(1.,abs(vexbp(i,j,k)))
             dth = xdelh(i,j,k) / amax1(1.,abs(vexbh(i,j,k)))
             dt1 = amin1 ( dts,dtp,dth )
             if ( dt1 .le. dtnew ) then
                 dtnew = dt1
                 i0    = i
                 j0    = j
                 k0    = k
             endif
           enddo
         enddo
       enddo
!         print *,'perpendicular',taskid,dtnew,i0,j0,k0,
!     .   vexbs(i0,j0,k0),vexbp(i0,j0,k0),vexbh(i0,j0,k0)

       if ( dtnew .le. .01 ) then
         print *,' Time step too small'
         stop
       elseif ( dtnew .ge. 5e4 ) then
         print *,' Time step too big: dtnew',dtnew
         stop
       endif

       dtnew = 1.5 * dtnew
       if ( dtnew/dt .le. 1.0  ) dt = amin1(dt0,dtnew   )
       if ( dtnew/dt .ge. 1.2  ) dt = amin1(dt0,dt * 1.2)

!       print *,'in courant',dt0,dtnew,dt

       return
       end

*******************************************
*******************************************

!             EXB

*******************************************
*******************************************

       subroutine exb(hrut,phi,phialt,philon)

       include 'param3_mpi-1.98.inc'
       include 'com3_mpi-1.98.inc'

       real denic(nz,nf,nl,nion)
       real tic(nz,nf,nl,nion)
       real tec(nz,nf,nl)
       real fluxnp(nz,nfp1,nl,nion),fluxtp(nz,nfp1,nl,nion)
       real fluxtep(nz,nfp1,nl)
       real fluxns(nz,nf,nl,nion),fluxts(nz,nf,nl,nion)
       real fluxtes(nz,nf,nl)
       real fluxnh(nz,nf,nl,nion),fluxth(nz,nf,nl,nion)
       real fluxteh(nz,nf,nl)
       real vexb_p(nzp1,nfp1,nlp1),vexb_h(nzp1,nfp1,nlp1)
       real param(2)
       real phi(nnx,nny),phialt(nnx,nny),philon(nnx,nny)

! define the e x b drift

      param(1) = day
      param(2) = f10p7
      nzh      = nz / 2 

! note: modification of vexb because of field line variation
!       uses cos^3/sqrt(1.+3sin^2) instead of
!       uses sin^3/sqrt(1.+3cos^2) because
!       blats = 0 at the magnetic equator 
!       (as opposed to pi/2 in spherical coordinates)

      call vexb_phi(phi,phialt,philon)
      call current_jp_jphi

! here we add in exb drift from the potential

! reduce e x b velocities below alt_crit
! with 20 km decay (dela)

! and kill above some high altitude

!      alt_crit_high   =  12000. * re
      alt_crit_high   =  pcrit * re
      dela_high     = 2. * re
      dela = 20.

!     vexbp

      do k = 1,nl
        do j = 1,nfp1
          do i = 1,nz
            vexbp(i,j,k) = vexbp_phi(i,j,k)
            if ( baltp(i,j,k) .lt. alt_crit ) then
              arg0 = ( alt_crit - baltp(i,j,k) ) / dela
              fac = exp(-arg0*arg0)
              vexbp(i,j,k) = vexbp(i,j,k) * fac
            endif
            if ( baltp(i,j,k) .gt. alt_crit_high ) then
              arg0 = ( abs(alt_crit_high - baltp(i,j,k)) ) / dela_high
              fac = exp(-arg0*arg0)
!              print *,'pre vexbp',arg0,fac,vexbp(i,j,k)
              vexbp(i,j,k) = vexbp(i,j,k) * fac
!              print *,'post vexbp',vexbp(i,j,k)
            endif
          enddo
        enddo
      enddo

!     vexbs

      do k = 1,nl
        do j = 1,nf
          do i = 1,nzp1
            vexbs(i,j,k) = vexbs_phi(i,j,k)
            if ( baltp(i,j,k) .lt. alt_crit ) then
              arg0 = ( alt_crit - baltp(i,j,k) ) / dela
              fac = exp(-arg0*arg0)
              vexbs(i,j,k) = vexbs(i,j,k) * fac
            endif
            if ( baltp(i,j,k) .gt. alt_crit_high ) then
              arg0 = ( abs(alt_crit_high - baltp(i,j,k)) ) / dela_high
              fac = exp(-arg0*arg0)
              vexbs(i,j,k) = vexbs(i,j,k) * fac
            endif
          enddo
        enddo
      enddo

!     vexbh

      do k = 1,nlp1
        do j = 1,nf
          do i = 1,nz
            vexbh(i,j,k) = vexbh_phi(i,j,k)
            if ( baltp(i,j,k) .lt. alt_crit ) then
              arg0 = ( alt_crit - baltp(i,j,k) ) / dela
              fac = exp(-arg0*arg0)
              vexbh(i,j,k) = vexbh(i,j,k) * fac
            endif
            if ( baltp(i,j,k) .gt. alt_crit_high ) then
              arg0 = ( abs(alt_crit_high - baltp(i,j,k)) ) / dela_high
              fac = exp(-arg0*arg0)
              vexbh(i,j,k) = vexbh(i,j,k) * fac
            endif
          enddo
        enddo
      enddo

! limit e x b velocities

      do k = 1,nl
        do j = 1,nfp1
          do i = 1,nz
            if (vexbp(i,j,k) .gt. 0 .) 
     .        vexbp(i,j,k) = amin1(vexbp(i,j,k),vexb_max)
            if (vexbp(i,j,k) .lt. 0 .) 
     .        vexbp(i,j,k) = amax1(vexbp(i,j,k),-vexb_max)
          enddo
        enddo
      enddo


      do k = 1,nl
        do j = 1,nf
          do i = 1,nzp1
            if (vexbs(i,j,k) .gt. 0 .) 
     .        vexbs(i,j,k) = amin1(vexbs(i,j,k),vexb_max)
            if (vexbs(i,j,k) .lt. 0 .) 
     .        vexbs(i,j,k) = amax1(vexbs(i,j,k),-vexb_max)
          enddo
        enddo
      enddo

      do k = 1,nlp1
        do j = 1,nf
          do i = 1,nz
            if (vexbh(i,j,k) .gt. 0 .) 
     .        vexbh(i,j,k) = amin1(vexbh(i,j,k),vexb_max)
            if (vexbh(i,j,k) .lt. 0 .) 
     .        vexbh(i,j,k) = amax1(vexbh(i,j,k),-vexb_max)
          enddo
        enddo
      enddo



! output e x b velocities

      do k = 1,nl
        do j = 1,nf
          do i = 1,nz
            u1p(i,j,k) = vexbp(i,j,k)
            u2s(i,j,k) = vexbs(i,j,k)
            u3h(i,j,k) = vexbh(i,j,k)
          enddo
        enddo
      enddo

! calculate conserved particle number: denic
! and 'conserved' temperature: tic,tec

       do ni = nion1,nion2
         do k = 1,nl
           do j = 1,nf
             do i = 1,nz
               denic(i,j,k,ni) = deni(i,j,k,ni) * vol(i,j,k)
               tic(i,j,k,ni)   = ti(i,j,k,ni) * vol(i,j,k)
             enddo
           enddo
         enddo
       enddo

       do k = 1,nl
         do j = 1,nf
           do i = 1,nz
             tec(i,j,k)   = te(i,j,k) * vol(i,j,k)
           enddo
         enddo
       enddo

! calculate flux in p-direction at interface
! NOTE: neutral flux condition at outer boundary (JH 11/29/07)
! altered to consider NS pole densities

       do ni = nion1,nion2
         do k = 1,nl
           do j = 2,nf
             do i = 1,nz
               if ( vexbp(i,j,k) .ge. 0 ) then
                 fluxnp(i,j,k,ni) = deni(i,j-1,k,ni) * vexbp(i,j,k)
                 fluxtp(i,j,k,ni) = ti(i,j-1,k,ni)   * vexbp(i,j,k)
               else
                 fluxnp(i,j,k,ni) = deni(i,j,k,ni) * vexbp(i,j,k)
                 fluxtp(i,j,k,ni) = ti(i,j,k,ni)   * vexbp(i,j,k)
!                 if ( j .eq. nf ) then 
!                   fluxnp(i,j,k,ni) = fluxnp(i,j-1,k,ni) *
!     .                                areap(i,j-1,k)/areap(i,j,k) 
!                   fluxtp(i,j,k,ni) = fluxtp(i,j-1,k,ni) *
!     .                                areap(i,j-1,k)/areap(i,j,k) 
!                 endif
               endif 
             enddo
           enddo
         enddo
       enddo

       do k = 1,nl
         do j = 2,nf
           do i = 1,nz
             if ( vexbp(i,j,k) .ge. 0 ) then
               fluxtep(i,j,k) = te(i,j-1,k) * vexbp(i,j,k)
             else
               fluxtep(i,j,k) = te(i,j,k)   * vexbp(i,j,k)
!               if ( j .eq. nf ) then 
!                 fluxtep(i,j,k) = fluxtep(i,j-1,k) *
!     .                                areap(i,j-1,k)/areap(i,j,k) 
!               endif 
             endif
           enddo
         enddo
       enddo

!      flux at nfp1 (near magnetic north/south poles)

       do ni = nion1,nion2
         do k = 1,nl
           j = nfp1
             do i = 1,nz
               if ( vexbp(i,j,k) .ge. 0 ) then
                 fluxnp(i,j,k,ni) = deni(i,j-1,k,ni) * vexbp(i,j,k)
                 fluxtp(i,j,k,ni) = ti(i,j-1,k,ni)   * vexbp(i,j,k)
               else
!!                 if ( altp (i,j,k) .lt. alt_crit_avg) then
                   fluxnp(i,j,k,ni) = deni_mnp(i,ni) * vexbp(i,j,k)
                   fluxtp(i,j,k,ni) = ti_mnp(i,ni)   * vexbp(i,j,k)
!!                 else
!!                   fluxnp(i,j,k,ni) = deni(i,j-1,k,ni) * vexbp(i,j,k)
!!                   fluxtp(i,j,k,ni) = ti(i,j-1,k,ni)   * vexbp(i,j,k)
!!                 endif
               endif
             enddo
         enddo
       enddo

       do k = 1,nl
         j = nfp1
           do i = 1,nz
             if ( vexbp(i,j,k) .ge. 0 ) then
               fluxtep(i,j,k) = te(i,j-1,k) * vexbp(i,j,k)
             else
!!               if ( altp (i,j,k) .lt. alt_crit_avg) then
                 fluxtep(i,j,k) = te_mnp(i)   * vexbp(i,j,k)
!!               else
!!                 fluxtep(i,j,k) = te(i,j-1,k) * vexbp(i,j,k)
!!               endif

               fluxtep(i,j,k) = te_mnp(i)   * vexbp(i,j,k)

             endif
         enddo
       enddo

! calculate flux in s-direction at interface

       do ni = nion1,nion2
         do k = 1,nl
           do j = 1,nf
             do i = 2,nz
               if ( vexbs(i,j,k) .ge. 0 ) then
                 fluxns(i,j,k,ni) = deni(i-1,j,k,ni) * vexbs(i,j,k)
                 fluxts(i,j,k,ni) = ti(i-1,j,k,ni)   * vexbs(i,j,k)
               else
                 fluxns(i,j,k,ni) = deni(i,j,k,ni) * vexbs(i,j,k)
                 fluxts(i,j,k,ni) = ti(i,j,k,ni)   * vexbs(i,j,k)
               endif
             enddo
           enddo
         enddo
       enddo

       do k = 1,nl
         do j = 1,nf
           do i = 2,nz
             if ( vexbs(i,j,k) .ge. 0 ) then
               fluxtes(i,j,k) = te(i-1,j,k) * vexbs(i,j,k)
             else
               fluxtes(i,j,k) = te(i,j,k)   * vexbs(i,j,k)
             endif
           enddo
         enddo
       enddo

! calculate flux in h-direction at interface (k > 1)

       do ni = nion1,nion2
         do k = 2,nl
           do j = 1,nf
             do i = 1,nz
               if ( vexbh(i,j,k) .ge. 0 ) then
                 fluxnh(i,j,k,ni) = deni(i,j,k-1,ni) * vexbh(i,j,k)
                 fluxth(i,j,k,ni) = ti(i,j,k-1,ni)   * vexbh(i,j,k)
               else
                 fluxnh(i,j,k,ni) = deni(i,j,k,ni) * vexbh(i,j,k)
                 fluxth(i,j,k,ni) = ti(i,j,k,ni)   * vexbh(i,j,k)
               endif
             enddo
           enddo
         enddo
       enddo

       do k = 2,nl
         do j = 1,nf
           do i = 1,nz
             if ( vexbh(i,j,k) .ge. 0 ) then
               fluxteh(i,j,k) = te(i,j,k-1) * vexbh(i,j,k)
             else
               fluxteh(i,j,k) = te(i,j,k)   * vexbh(i,j,k)
             endif
           enddo
         enddo
       enddo

!      calculate flux in h-direction at interface (k = 1)                    
!      (invoke periodic boundary condition)                                    
                                                                               
       do ni = nion1,nion2
         do j = 1,nf
           do i = 1,nz
             if ( vexbh(i,j,1) .ge. 0 ) then
               fluxnh(i,j,1,ni) = deni(i,j,nl,ni) * vexbh(i,j,1) 
               fluxth(i,j,1,ni) = ti(i,j,nl,ni)   * vexbh(i,j,1)
             else
               fluxnh(i,j,1,ni) = deni(i,j,1,ni) * vexbh(i,j,1)
               fluxth(i,j,1,ni) = ti(i,j,1,ni)   * vexbh(i,j,1)
             endif
           enddo
         enddo
       enddo

       do j = 1,nf
         do i = 1,nz
           if ( vexbh(i,j,1) .ge. 0 ) then
             fluxteh(i,j,1) = te(i,j,nl) * vexbh(i,j,1)
           else
             fluxteh(i,j,1) = te(i,j,1)  * vexbh(i,j,1)
           endif 
         enddo
       enddo

! update total particle number and density
! and temperatures 
! NOTE: the temperature update is an approximation
!       (probably better than no update but, strictly
!       speaking, not exactly correct)

       do ni = nion1,nion2
         do k = 1,nlm1
           do j = 2,nf
             do i = 2,nzm1
               denic(i,j,k,ni) = denic(i,j,k,ni) 
     .                   + dt * ( areap(i,j,k)   * fluxnp(i,j,k,ni) -
     .                            areap(i,j+1,k) * fluxnp(i,j+1,k,ni) )
     .                   + dt * ( areas(i,j,k)   * fluxns(i,j,k,ni) -
     .                            areas(i+1,j,k) * fluxns(i+1,j,k,ni) )
     .                   + dt * ( areah(i,j,k)   * fluxnh(i,j,k,ni) -
     .                            areah(i,j,k+1) * fluxnh(i,j,k+1,ni) )
               deni(i,j,k,ni)  = denic(i,j,k,ni) / vol(i,j,k)

! brazen fix
                deni(i,j,k,ni)  = amax1(deni(i,j,k,ni),denmin)

               tic(i,j,k,ni) = tic(i,j,k,ni) 
     .                   + dt * ( areap(i,j,k)   * fluxtp(i,j,k,ni) -
     .                            areap(i,j+1,k) * fluxtp(i,j+1,k,ni) )
     .                   + dt * ( areas(i,j,k)   * fluxts(i,j,k,ni) -
     .                            areas(i+1,j,k) * fluxts(i+1,j,k,ni) )
     .                   + dt * ( areah(i,j,k)   * fluxth(i,j,k,ni) -
     .                            areah(i,j,k+1) * fluxth(i,j,k+1,ni) )
               ti(i,j,k,ni)  = tic(i,j,k,ni) / vol(i,j,k)
! brazen fix
               ti(i,j,k,ni)  = amax1(ti(i,j,k,ni),200.)
!               if (ti(i,j,k,ni) .lt. 0.) print *,'i,j,k,ni',
!     .                      i,j,k,ni,ti(i,j,k,ni)
               if (isnan(ti(i,j,k,ni))) then
                 print *,'Ti fixed',i,j,k,ni
                 ti(i,j,k,ni) = 200.
               endif
             enddo
           enddo
         enddo
       enddo

       do k = 1,nlm1
         do j = 2,nf
           do i = 2,nzm1
             tec(i,j,k) = tec(i,j,k) 
     .                     + dt * ( areap(i,j,k)   * fluxtep(i,j,k) -
     .                              areap(i,j+1,k) * fluxtep(i,j+1,k) )
     .                     + dt * ( areas(i,j,k)   * fluxtes(i,j,k) -
     .                              areas(i+1,j,k) * fluxtes(i+1,j,k) )
     .                     + dt * ( areah(i,j,k)   * fluxteh(i,j,k) -
     .                              areah(i,j,k+1) * fluxteh(i,j,k+1) )
             te(i,j,k)  = tec(i,j,k) / vol(i,j,k)
! brazen fix
               te(i,j,k)  = amax1(te(i,j,k),200.)
               if (te(i,j,k) .lt. 0.) print *,'i,j,k',
     .                      i,j,k,te(i,j,k)
           enddo
         enddo
       enddo

!      for k = nl

       do ni = nion1,nion2
         do j = 2,nf
           do i = 2,nzm1
             k                = nl
             deni0            = deni(i,j,k,ni)
             ti0              = ti(i,j,k,ni)
             denic(i,j,nl,ni) = denic(i,j,nl,ni) 
     .                 + dt * ( areap(i,j,nl)   * fluxnp(i,j,nl,ni) -
     .                          areap(i,j+1,nl) * fluxnp(i,j+1,nl,ni) )
     .                 + dt * ( areas(i,j,nl)   * fluxns(i,j,nl,ni) -
     .                          areas(i+1,j,nl) * fluxns(i+1,j,nl,ni) )
     .                 + dt * ( areah(i,j,nl)   * fluxnh(i,j,nl,ni) -
     .                          areah(i,j,1)    * fluxnh(i,j,1,ni) )
              deni(i,j,nl,ni)  = denic(i,j,nl,ni) / vol(i,j,nl)
              if ( deni(i,j,nl,ni) .le. 0. ) 
     .             deni(i,j,nl,ni) = deni0
              tic(i,j,nl,ni) = tic(i,j,nl,ni) 
     .                 + dt * ( areap(i,j,nl)   * fluxtp(i,j,nl,ni) -
     .                          areap(i,j+1,nl) * fluxtp(i,j+1,nl,ni) )
     .                 + dt * ( areas(i,j,nl)   * fluxts(i,j,nl,ni) -
     .                          areas(i+1,j,nl) * fluxts(i+1,j,nl,ni) )
     .                 + dt * ( areah(i,j,nl)   * fluxth(i,j,nl,ni) -
     .                          areah(i,j,1)    * fluxth(i,j,1,ni) )
             ti(i,j,nl,ni)  = tic(i,j,nl,ni) / vol(i,j,nl)
             if ( ti(i,j,nl,ni) .le. 0. ) 
     .            ti(i,j,nl,ni) = ti0
            enddo
         enddo
       enddo

       do j = 2,nf
         do i = 2,nzm1
           te0 = te(i,j,nl)
           tec(i,j,nl) = tec(i,j,nl) 
     .                   + dt * ( areap(i,j,nl)   * fluxtep(i,j,nl) -
     .                            areap(i,j+1,nl) * fluxtep(i,j+1,nl) )
     .                   + dt * ( areas(i,j,nl)   * fluxtes(i,j,nl) -
     .                            areas(i+1,j,nl) * fluxtes(i+1,j,nl) )
     .                   + dt * ( areah(i,j,nl)   * fluxteh(i,j,nl) -
     .                            areah(i,j,1)    * fluxteh(i,j,1) )
           te(i,j,nl)  = tec(i,j,nl) / vol(i,j,nl)
           if ( te(i,j,nl) .le. 0. ) 
     .          te(i,j,nl) = te0
         enddo
       enddo



! fill cells at j = 1 and nf with j = 2 and nfm1

       do ni = nion1,nion2
         do k = 1,nl
           do i = 2,nzm1
             deni(i,1,k,ni)  = deni(i,2,k,ni)
!             deni(i,nf,k,ni) = deni(i,nfm1,k,ni)
             ti(i,1,k,ni)    = ti(i,2,k,ni)
!             ti(i,nf,k,ni)   = ti(i,nfm1,k,ni)
           enddo
         enddo
       enddo

       do k = 1,nl
         do i = 2,nzm1
           te(i,1,k)    = te(i,2,k)
!           te(i,nf,k)   = te(i,nfm1,k)
         enddo
       enddo

! fill cells at i = 1 and nz with i = 2 and nzm1

       do ni = nion1,nion2
         do k = 2,nlm1
           do j = 1,nf
             deni(1,j,k,ni)  = deni(2,j,k,ni)
             deni(nz,j,k,ni) = deni(nzm1,j,k,ni)
             ti(1,j,k,ni)    = ti(2,j,k,ni)
             ti(nz,j,k,ni)   = ti(nzm1,j,k,ni)
           enddo
         enddo
       enddo

       do k = 2,nlm1
         do j = 1,nf
           te(1,j,k)    = te(2,j,k)
           te(nz,j,k)   = te(nzm1,j,k)
         enddo
       enddo

!   fix at j = nf (interpolate) 
 
      do ni = nion1,nion2 
        do k = 2,nlm1 
        j = nf 
          do i = 1,nz 
            slope  = (blatp(i,j,k) - blatp(i,j-1,k)  ) /
     &            (90.          - blatp(i,j-1,k)) 
            deni(i,j,k,ni) = deni(i,j-1,k,ni) +  
     &           slope * (deni_mnp(i,ni) - deni(i,j-1,k,ni))  
            deni(i,j,k,ni) = max(deni(i,j,k,ni),denmin)           
          enddo 
        enddo 
      enddo 
 

       return
       end


*******************************************
*******************************************

!             vdrift_model

*******************************************
*******************************************

C       ************************************************************
C       ************************************************************

	subroutine vdrift_model(xt,xl,param,y,z,fejer,ver,veh)

C       ************************************************************

C       ************************************************************
C       SUBROUTINE CALCULATES EQUATORIAL VERTICAL DRIFT AS DESCRIBED 
C       IN SCHERLIESS AND FEJER, JGR, 104, 6829-6842, 1999
C       ************************************************************

C       INPUT:   XT: SOLAR LOCAL TIME
C                XL: GEOGRAPHIC LONGITUDE (+ EAST)
C               
C             PARAM: 2-DIM ARRAY (DOY,F10.7CM)
C                    DOY     :Day of Year has to run from 1 to 365 (366)
C                    F10.7cm : F10.7cm solar flux
C             
C       OUTPUT:   Y: EQUATORIAL VERTICAL DRIFT

C       ************************************************************

!  z: longitudinal drift
!  fejer: logical variable (true: use fejer/scherliess model;
!                           false: use sinusoidal model)
!  ver: max velocity for 'radial' e x b drift for sinusoidal model
!  veh: max velocity for 'longitudinal' e x b drift for sinusoidal model


        include 'param3_mpi-1.98.inc'

        real ver
        logical fejer

        real param(2),coeff(624),funct(6)
        real coeff1(312),coeff2(312)
	real xt,xl,y
	real bspl4,bspl4_time,bspl4_long

	integer i,j,ind,il,kk
	integer index_t/13/,dim_t/78/
	integer index_l/8/,dim_l/48/
	integer index/104/,dim/624/
	integer nfunc/6/

        data coeff1/
     *  -10.80592, -9.63722,-11.52666, -0.05716, -0.06288,  0.03564,
     *   -5.80962, -7.86988, -8.50888, -0.05194, -0.05798, -0.00138,
     *    2.09876,-19.99896, -5.11393, -0.05370, -0.06585,  0.03171,
     *  -10.22653, -3.62499,-14.85924, -0.04023, -0.01190, -0.09656,
     *   -4.85180,-26.26264, -6.20501, -0.05342, -0.05174,  0.02419,
     *  -13.98936,-18.10416, -9.30503, -0.01969, -0.03132, -0.01984,
     *  -18.36633,-24.44898,-16.69001,  0.02033, -0.03414, -0.02062,
     *  -20.27621,-16.95623,-36.58234,  0.01445, -0.02044, -0.08297,
     *    1.44450,  5.53004,  4.55166, -0.02356, -0.04267,  0.05023,
     *    5.50589,  7.05381,  1.94387, -0.03147, -0.03548,  0.01166,
     *    3.24165, 10.05002,  4.26218, -0.03419, -0.02651,  0.07456,
     *    7.02218,  0.06708,-11.31012, -0.03252, -0.01021, -0.09008,
     *   -3.47588, -2.82534, -4.17668, -0.03719, -0.01519,  0.06507,
     *   -4.02607,-11.19563,-10.52923, -0.00592, -0.01286, -0.00477,
     *  -11.47478, -9.57758,-10.36887,  0.04555, -0.02249,  0.00528,
     *  -14.19283,  7.86422, -8.76821,  0.05758, -0.02398, -0.04075,
     *   14.58890, 36.63322, 27.57497,  0.01358, -0.02316,  0.04723,
     *   12.53122, 29.38367, 21.40356, -0.00071, -0.00553,  0.01484,
     *   18.64421, 26.27327, 18.32704,  0.00578,  0.03349,  0.11249,
     *    4.53014,  6.15099,  7.41935, -0.02860, -0.00395, -0.08394,
     *   14.29422,  9.77569,  2.85689, -0.00107,  0.04263,  0.10739,
     *    7.17246,  4.40242, -1.00794,  0.00089,  0.01436,  0.00626,
     *    7.75487,  5.01928,  4.36908,  0.03952, -0.00614,  0.03039,
     *   10.25556,  8.82631, 24.21745,  0.05492, -0.02968,  0.00177,
     *   21.86648, 24.03218, 39.82008,  0.00490, -0.01281, -0.01715,
     *   19.18547, 23.97403, 34.44242,  0.01978,  0.01564, -0.02434,
     *   26.30614, 14.22662, 31.16844,  0.06495,  0.19590,  0.05631,
     *   21.09354, 25.56253, 29.91629, -0.04397, -0.08079, -0.07903,
     *   28.30202, 16.80567, 38.63945,  0.05864,  0.16407,  0.07622,
     *   22.68528, 25.91119, 40.45979, -0.03185, -0.01039, -0.01206,
     *   31.98703, 24.46271, 38.13028, -0.08738, -0.00280,  0.01322,
     *   46.67387, 16.80171, 22.77190, -0.13643, -0.05277, -0.01982,
     *   13.87476, 20.52521,  5.22899,  0.00485, -0.04357,  0.09970,
     *   21.46928, 13.55871, 10.23772, -0.04457,  0.01307,  0.06589,
     *   16.18181, 16.02960,  9.28661, -0.01225,  0.14623, -0.01570,
     *   18.16289, -1.58230, 14.54986, -0.00375, -0.00087,  0.04991,
     *   10.00292, 11.82653,  0.44417, -0.00768,  0.15940, -0.01775,
     *   12.15362,  5.65843, -1.94855, -0.00689,  0.03851,  0.04851,
     *   -1.25167,  9.05439,  0.74164,  0.01065,  0.03153,  0.02433,
     *  -15.46799, 18.23132, 27.45320,  0.00899, -0.00017,  0.03385,
     *    2.70396, -0.87077,  6.11476, -0.00081,  0.05167, -0.08932,
     *    3.21321, -1.06622,  5.43623,  0.01942,  0.05449, -0.03084,
     *   17.79267, -3.44694,  7.10702,  0.04734, -0.00945,  0.11516,
     *    0.46435,  6.78467,  4.27231, -0.02122,  0.10922, -0.03331,
     *   15.31708,  1.70927,  7.99584,  0.07462,  0.07515,  0.08934,
     *    4.19893,  6.01231,  8.04861,  0.04023,  0.14767, -0.04308,
     *    9.97541,  5.99412,  5.93588,  0.06611,  0.12144, -0.02124,
     *   13.02837, 10.29950, -4.86200,  0.04521,  0.10715, -0.05465,
     *    5.26779,  7.09019,  1.76617,  0.09339,  0.22256,  0.09222,
     *    9.17810,  5.27558,  5.45022,  0.14749,  0.11616,  0.10418,
     *    9.26391,  4.19982, 12.66250,  0.11334,  0.02532,  0.18919,
     *   13.18695,  6.06564, 11.87835,  0.26347,  0.02858,  0.14801/
 
         data coeff2/
     *   10.08476,  6.14899, 17.62618,  0.09331,  0.08832,  0.28208,
     *   10.75302,  7.09244, 13.90643,  0.09556,  0.16652,  0.22751,
     *    6.70338, 11.97698, 18.51413,  0.15873,  0.18936,  0.15705,
     *    5.68102, 23.81606, 20.65174,  0.19930,  0.15645,  0.08151,
     *   29.61644,  5.49433, 48.90934,  0.70710,  0.40791,  0.26325,
     *   17.11994, 19.65380, 44.88810,  0.45510,  0.41689,  0.22398,
     *    8.45700, 34.54442, 27.25364,  0.40867,  0.37223,  0.22374,
     *   -2.30305, 32.00660, 47.75799,  0.02178,  0.43626,  0.30187,
     *    8.98134, 33.01820, 33.09674,  0.33703,  0.33242,  0.41156,
     *   14.27619, 20.70858, 50.10005,  0.30115,  0.32570,  0.45061,
     *   14.44685, 16.14272, 45.40065,  0.37552,  0.31419,  0.30129,
     *    6.19718, 18.89559, 28.24927,  0.08864,  0.41627,  0.19993,
     *    7.70847, -2.36281,-21.41381,  0.13766,  0.05113, -0.11631,
     *   -9.07236,  3.76797,-20.49962,  0.03343,  0.08630,  0.00188,
     *   -8.58113,  5.06009, -6.23262,  0.04967,  0.03334,  0.24214,
     *  -27.85742,  8.34615,-27.72532, -0.08935,  0.15905, -0.03655,
     *    2.77234,  0.14626, -4.01786,  0.22338, -0.04478,  0.18650,
     *    5.61364, -3.82235,-16.72282,  0.26456, -0.03119, -0.08376,
     *   13.35847, -6.11518,-16.50327,  0.28957, -0.01345, -0.19223,
     *   -5.37290, -0.09562,-27.27889,  0.00266,  0.22823, -0.35585,
     *  -15.29676,-18.36622,-24.62948, -0.31299, -0.23832, -0.08463,
     *  -23.37099,-13.69954,-26.71177, -0.19654, -0.18522, -0.20679,
     *  -26.33762,-15.96657,-42.51953, -0.13575, -0.00329, -0.28355,
     *  -25.42140,-14.14291,-21.91748, -0.20960, -0.19176, -0.32593,
     *  -23.36042,-23.89895,-46.05270, -0.10336,  0.03030, -0.21839,
     *  -19.46259,-21.27918,-32.38143, -0.17673, -0.15484, -0.11226,
     *  -19.06169,-21.13240,-34.01677, -0.25497, -0.16878, -0.11004,
     *  -18.39463,-16.11516,-19.55804, -0.19834, -0.23271, -0.25699,
     *  -19.93482,-17.56433,-18.58818,  0.06508, -0.18075,  0.02796,
     *  -23.64078,-18.77269,-22.77715, -0.02456, -0.12238,  0.02959,
     *  -12.44508,-21.06941,-19.36011,  0.02746, -0.16329,  0.19792,
     *  -26.34187,-19.78854,-24.06651, -0.07299, -0.03082, -0.03535,
     *  -10.71667,-26.04401,-16.59048,  0.02850, -0.09680,  0.15143,
     *  -18.40481,-23.37770,-16.31450, -0.03989, -0.00729, -0.01688,
     *   -9.68886,-20.59304,-18.46657,  0.01092, -0.07901,  0.03422,
     *   -0.06685,-19.24590,-29.35494,  0.12265, -0.24792,  0.05978,
     *  -15.32341, -9.07320,-13.76101, -0.17018, -0.15122, -0.06144,
     *  -14.68939,-14.82251,-13.65846, -0.11173, -0.14410, -0.07133,
     *  -18.38628,-18.94631,-19.00893, -0.08062, -0.14481, -0.12949,
     *  -16.15328,-17.40999,-14.08705, -0.08485, -0.06896, -0.11583,
     *  -14.50295,-16.91671,-25.25793, -0.06814, -0.13727, -0.12213,
     *  -10.92188,-14.10852,-24.43877, -0.09375, -0.11638, -0.09053,
     *  -11.64716,-14.92020,-19.99063, -0.14792, -0.08681, -0.12085,
     *  -24.09766,-16.14519, -8.05683, -0.24065, -0.05877, -0.23726,
     *  -25.18396,-15.02034,-15.50531, -0.12236, -0.09610, -0.00529,
     *  -15.27905,-19.36708,-12.94046, -0.08571, -0.09560, -0.03544,
     *   -7.48927,-16.00753,-13.02842, -0.07862, -0.10110, -0.05807,
     *  -13.06383,-27.98698,-18.80004, -0.05875, -0.03737, -0.11214,
     *  -13.67370,-16.44925,-16.12632, -0.07228, -0.09322, -0.05652,
     *  -22.61245,-21.24717,-18.09933, -0.05197, -0.07477, -0.05235,
     *  -27.09189,-21.85181,-20.34676, -0.05123, -0.05683, -0.07214,
     *  -27.09561,-22.76383,-25.41151, -0.10272, -0.02058, -0.16720/

        do i = 1,312
          coeff(i)     = coeff1(i)
          coeff(i+312) = coeff2(i)
        enddo

        xt = mod(xt,24.)

!  sinusoidal e x b model

!  longitudinal drift

        z = -veh * sin ( 2 * pie * ( xt - 3. ) / 24. )

!  'radial' drift

!       sinusoid e x b model
!       JK: put in pre-reversal enhancement in addition to sinusoid

        vpre = 40.0
        dpre = 1.0

        if ( .not. fejer ) then
          y = ver * sin ( 2 * pie * ( xt - 7. ) / 24. )
     .        + vpre*exp(-((xt-19)/dpre)**2)
          return
        endif

!  fejer-scherliess e x b model

	call g(param,funct,xl)

C       **********************************
	y=0.
C       **********************************
	do i=1,index_t
	  do il=1,index_l
	    kk=index_l*(i-1)+il
	    do j=1,nfunc
	       ind=nfunc*(kk-1)+j
	       bspl4=bspl4_time(i,xt)*bspl4_long(il,xl)
               y=y+bspl4*funct(j)*coeff(ind)
	    end do
          end do
         end do

	end

c       *************************************************
c       *************************************************
        real function bspl4_time(i,x1)
c       *************************************************
	implicit none 

	integer i,order/4/,j,k
	real t_t(0:39)
	real x,b(20,20),x1

        data t_t/
     *          0.00,2.75,4.75,5.50,6.25,
     *          7.25,10.00,14.00,17.25,18.00,
     *          18.75,19.75,21.00,24.00,26.75,
     *          28.75,29.50,30.25,31.25,34.00,
     *          38.00,41.25,42.00,42.75,43.75,
     *          45.00,48.00,50.75,52.75,53.50,
     *          54.25,55.25,58.00,62.00,65.25,
     *          66.00,66.75,67.75,69.00,72.00/
       
	x=x1
        if(i.ge.0) then
          if (x.lt.t_t(i-0)) then
	     x=x+24
	  end if
	end if
	do j=i,i+order-1
	   if(x.ge.t_t(j).and.x.lt.t_t(j+1)) then
	       b(j,1)=1
	   else
	       b(j,1)=0
	   end if
	end do

	do j=2,order
	     do k=i,i+order-j
		b(k,j) = ( x - t_t(k) ) / ( t_t(k+j-1) - t_t(k) ) * 
     .                   b(k,j-1)
		b(k,j) = b(k,j) + 
     .                   ( t_t(k+j)-x ) / ( t_t(k+j) - t_t(k+1) ) *
     .                    b(k+1,j-1)
	     end do
	end do

	bspl4_time=b(i,order)
	end


c       *************************************************
c       *************************************************
        real function bspl4_long(i,x1)
c       *************************************************
	implicit none 

	integer i,order/4/,j,k
	real t_l(0:24)
	real x,b(20,20),x1

        data t_l/
     *          0,10,100,190,200,250,280,310,
     *          360,370,460,550,560,610,640,670,
     *          720,730,820,910,920,970,1000,1030,1080/
       
	x=x1
        if(i.ge.0) then
          if (x.lt.t_l(i-0)) then
	     x=x+360
	  end if
	end if
	do j=i,i+order-1
	   if(x.ge.t_l(j).and.x.lt.t_l(j+1)) then
	       b(j,1)=1
	   else
	       b(j,1)=0
	   end if
	end do

	do j=2,order
	     do k=i,i+order-j
		b(k,j)=(x-t_l(k))/(t_l(k+j-1)-t_l(k))*b(k,j-1)
		b(k,j)=b(k,j)+(t_l(k+j)-x)/(t_l(k+j)-t_l(k+1))*
     .                 b(k+1,j-1)
	     end do
	end do

	bspl4_long=b(i,order)
	end

c       *************************************************
c       *************************************************
        subroutine g(param,funct,x)
c       *************************************************
	implicit none

        integer i
	real param(2),funct(6)
	real x,a,sigma,gauss,flux,cflux

c       *************************************************
	flux=param(2)
        if(param(2).le.75)  flux=75.
        if(param(2).ge.230) flux=230.
	cflux=flux

	a=0.
        if((param(1).ge.120).and.(param(1).le.240)) a=170.
        if((param(1).ge.120).and.(param(1).le.240)) sigma=60
        if((param(1).le.60).or.(param(1).ge.300)) a=170.
        if((param(1).le.60).or.(param(1).ge.300)) sigma=40

	if((flux.le.95).and.(a.ne.0)) then
	 gauss=exp(-0.5*((x-a)**2)/sigma**2)
         cflux=gauss*95.+(1-gauss)*flux
        end if
c       *************************************************

c       *************************************************
        do i=1,6
         funct(i)=0.
        end do
c       *************************************************

c       *************************************************
        if((param(1).ge.135).and.(param(1).le.230)) funct(1)=1
        if((param(1).le.45).or.(param(1).ge.320)) funct(2)=1
        if((param(1).gt.75).and.(param(1).lt.105)) funct(3)=1
        if((param(1).gt.260).and.(param(1).lt.290)) funct(3)=1
c       *************************************************

        if((param(1).ge.45).and.(param(1).le.75)) then  ! W-E
	 funct(2)=1.-(param(1)-45.)/30.
	 funct(3)=1-funct(2)
        end if
        if((param(1).ge.105).and.(param(1).le.135)) then  ! E-S
	 funct(3)=1.-(param(1)-105.)/30.
	 funct(1)=1-funct(3)
        end if
        if((param(1).ge.230).and.(param(1).le.260)) then  ! S-E
	 funct(1)=1.-(param(1)-230.)/30.
	 funct(3)=1-funct(1)
        end if
        if((param(1).ge.290).and.(param(1).le.320)) then  ! E-W
	 funct(3)=1.-(param(1)-290.)/30.
	 funct(2)=1-funct(3)
        end if

c       *************************************************
        funct(4)=(cflux-140)*funct(1)
        funct(5)=(cflux-140)*funct(2)
        funct(6)=(flux-140)*funct(3)
c       *************************************************

	end


! *********************
!
!     smoothz
!
! *********************

      subroutine smoothz(finout,ncomp) 
 
      include 'param3_mpi-1.98.inc' 
 
      dimension finout(nz), tempz(nz) 
 
c 
c This is the binomial filter (in x space) as described in  
c Birdsall appendix C. 
c We have the choice of a compensating filter or not. 
c if ncomp=0, no compensation, else compensation 
c 
 
c do smoothz in the z direction 
 
       do i = 1,nz 
          ip1 = i +1 
          if(i .eq. nz) ip1 = 1 
          im1 = i -1 
          if(i .eq. 1) im1 = nz 
          tempz(i) = .25*(finout(im1) +2.*finout(i)  
     &                   +finout(ip1)) 
       enddo 
       do i = 1,nz 
          finout(i) = tempz(i) 
       enddo 
 
       if ( ncomp .ne. 0 ) then 
 
c put in compensator  
c the alogrithm below is equivalent to  
c fftmp(i)=(1./16.)*(-ff0(i-2)+4.*ff0(i-1)+10.*ff0(i)+4.*ff0(i+1)-ff0(i+2)) 
 
c do compensation in the z direction 
 
       const = sqrt(1.4571072) 
       do i = 1,nz  
          ip1 = i +1 
          if(i .eq. nz) ip1 = 1 
          finout(i) = const*(finout(i) -.171573*finout(ip1)) 
       enddo 
       do i = nz,1,-1 
          im1 = i -1 
          if(i .eq. 1) im1 = nz 
          finout(i) = const*(finout(i) -.171573*finout(im1)) 
       enddo 
 
      endif 
 
      return 
      end 


!     ********************************************
!
!     VEXB_PHI
!
!     ********************************************

      subroutine vexb_phi (phi_sami,phialt,philon)
      use ModCoupleCimi,ONLY:DoCoupleCimi,PotCimiOnSamiGrid_C
      
      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc' 

      real phi_sami(nnx,nny),phi(nnx,nny),phialt(nnx,nny),philon(nnx,nny)

      if (DoCoupleCimi) then
         phi=phi_sami+PotCimiOnSamiGrid_C
      else
         phi=phi_sami
      endif
      
!      p face

       call phi_nfp1_nlp1(bradss,blonss,phihp,phi,phialt,philon)

!       open(1423,file='phiexb.dat',form='unformatted')
!       write(1423) phihp


       do k = 1,nl
         do j = 1,nfp1
           do i = 1,nz
             sinangx =  bdirpy(i,j,k) * ehpz(i,j,k) -
     .                  bdirpz(i,j,k) * ehpy(i,j,k)
             sinangy = -bdirpx(i,j,k) * ehpz(i,j,k) +
     .                  bdirpz(i,j,k) * ehpx(i,j,k)
             sinangz =  bdirpx(i,j,k) * ehpy(i,j,k) -
     .                  bdirpy(i,j,k) * ehpx(i,j,k)
             sinang  =  sqrt( sinangx * sinangx +
     .                        sinangy * sinangy +
     .                        sinangz * sinangz  )
!        if (j.eq.10 .and. k.eq.4.and.taskid.eq.4)
!     .   print *,'1',i,ehpx(i,j,k),ehpy(i,j,k),ehpz(i,j,k)

!              sinang     = 1.

              ehp(i,j,k) = -1. *
     .               ( phihp(j,k+1) - phihp(j,k) )/ delhp(i,j,k) /
     .                                              sinang

!             ehp(i,j,k) = -1.e-2/3. *
!     .               ( phihp(j,k+1) - phihp(j,k) )/ delhp(i,j,k) /
!     .                                              sinang
!             ehp(i,j,k) = -1.e-2/3. *
!     .               ( phi(j,k+1) - phi(j,k) )/ delhp(i,j,k) /
!     .                                              sinang
!       if(i.eq.11 .and.k.eq.4)print *,'ehp',sinang

           enddo
         enddo
       enddo

       do k = 1,nl
         do j = 1,nfp1
           do i = 1,nz
             vexbp_phi(i,j,k) =  ehp(i,j,k) * sol / bmag / bmpf(i,j,k) 
     .                            * tvexb0
           enddo
         enddo
       enddo

!       if ( taskid .eq. 1 ) then
!         i = nz/2
!         k = nl/2
!         do j = 1,nfp1
!           print *,'phihp,ehp,vexbp_phi',
!     .              i,j,k,phihp(j,k),ehp(i,j,k),vexbp_phi(i,j,k)
!         enddo
!       endif


!      print *, ' '   

!      h face

       do k = 1,nlp1
         do j = 1,nf
           do i = 1,nz
             sinangx =  bdirhy(i,j,k) * ephz(i,j,k) -
     .                  bdirhz(i,j,k) * ephy(i,j,k)
             sinangy = -bdirhx(i,j,k) * ephz(i,j,k) +
     .                  bdirhz(i,j,k) * ephx(i,j,k)
             sinangz =  bdirhx(i,j,k) * ephy(i,j,k) -
     .                  bdirhy(i,j,k) * ephx(i,j,k)
             sinang  =  sqrt( sinangx * sinangx +
     .                        sinangy * sinangy +
     .                        sinangz * sinangz  )

!        if (j.eq.100 .and. k.eq.4 .and. taskid.eq.1)
!     .   print *,'2',i,sinang,delph(i,j,k),sinang*delph(i,j,k)

!             sinang     = 1.
 
             eph(i,j,k) = -1. *
     .               ( phihp(j+1,k) - phihp(j,k) )/ 
     .                 delph(i,j,k)/sinang

           nzh = nz / 2

!           if (taskid.eq.1 .and. k.eq.2 .and. i.eq.nzh-1) 
!     .  print *,'sami3',j,phi(j+1,k),phi(j,k),hden,eph(i,j,k)
!           if(i.eq.nzh.and.j.eq.48.and.taskid.eq.6)
!     .  print *,taskid,hden,delph(i,j,k),sinang,phihp(j,k),bmhf(i,j,k)
!       if(i.eq.11 .and.k.eq.4)print *,'eph',sinang
!             eph(i,j,k) = -1.e-2/3. *
!     .               ( phihp(j+1,k) - phihp(j,k) )/ delph(i,j,k) /
!     .                                              sinang
!             eph(i,j,k) = -1.e-2/3. *
!     .               ( phi(j+1,k) - phi(j,k) )/ delph(i,j,k) /
!     .                                          sinang
           enddo
         enddo
       enddo

!      trying to smooth eph 
!      hopefully to help with region transitions

!      try  without it (JH 12/15/11)

!!       do k = 1,nlp1
!!         do j = 3,nf-2
!!           do i = 1,nz
!!             eph(i,j,k) = -1. *
!!     .             (phihp(j+3,k) - phihp(j-2,k))/ 
!!     .             (delph(i,j+1,k) + delph(i,j,k) + delph(i,j-1,k)+
!!     .              delph(i,j+2,k) + delph(i,j-2,k))
!!           enddo
!!         enddo
!!       enddo

       do k = 1,nlp1
         do j = 1,nf-1
           do i = 1,nz
             vexbh_phi(i,j,k) =  -eph(i,j,k) * sol / bmag / bmhf(i,j,k)
     .                            * tvexb0
!             vexbh_phi(i,j,k) =  -eph(i,j,k) * sol * bmhf(i,j,k) / bmag
!             vexbh_phi(i,j,k) =  0.
           enddo
         enddo
       enddo

       do k = 1,nlp1
           do i = 1,nz
             vexbh_phi(i,nf,k) = vexbh_phi(i,nf-1,k)
     .                            * tvexb0
           enddo
       enddo

!       do k = 1,nl
!         do j = 1,nf-1
!           do i = 1,nz
!             u1(i,j,k) =  eph(i,j,k)
!             u2(i,j,k) =  delph(i,j,k)
!             u3(i,j,k) =  phihp(j,k)
!             u4(i,j,k) =   PotCimiOnSamiGrid_C(k,j)
!           enddo
!         enddo
!       enddo

!       do k = 1,nl
!         do j = 1,nf
!           do i = 1,nz
!             u2(i,j,k) =  delph(i,j,k)
!           enddo
!         enddo
!       enddo

!      s face

       call phi_nfp1_nl(bradsh,blonsh,phish,phi,phialt,philon)
       call phi_nf_nlp1  (bradsp,blonsp,phisp,phi,phialt,philon)

       do k = 1,nl
         do j = 1,nf
           do i = 1,nzp1
             eps(i,j,k) = -1. * 
     .               ( phish(j+1,k) - phish(j,k) )/delps(i,j,k)
!             eps(i,j,k) = -1.e-2/3. * 
!     .               ( phish(j+1,k) - phish(j,k) )/delps(i,j,k)
           enddo
         enddo
       enddo

       do k = 1,nl
         do j = 1,nf
           do i = 1,nzp1
             ehs(i,j,k) = -1. * 
     .               ( phisp(j,k+1) - phisp(j,k) )/delhs(i,j,k)
!             ehs(i,j,k) = -1.e-2/3. * 
!     .               ( phisp(j,k+1) - phisp(j,k) )/delhs(i,j,k)
!             ehs(i,j,k) = -1.e-2/3. * 
!     .               ( phi(j,k+1) - phi(j,k) )/delhs(i,j,k)
           enddo
         enddo
       enddo

!       do k = 1,nl
!         do j = 1,nf
!           do i = 1,nz
!             u3(i,j,k) =  delps(i,j,k)
!             u4(i,j,k) =  eps(i,j,k)
!           enddo
!         enddo
!       enddo


       do k = 1,nl
         do j = 1,nf
           do i = 1,nz
!             vps =  eps(i,j,k) * sol / bmag / bmsf(i,j,k) *
!     .             ( (epsy(i,j,k) * bdirsz(i,j,k) - 
!     .                epsz(i,j,k) * bdirsy(i,j,k) ) * xnorms(i,j,k) -
!     .               (epsx(i,j,k) * bdirsz(i,j,k) -
!     .                epsz(i,j,k) * bdirsx(i,j,k) ) * ynorms(i,j,k) +
!     .               (epsx(i,j,k) * bdirsy(i,j,k) -
!     .                epsy(i,j,k) * bdirsx(i,j,k) ) * znorms(i,j,k) )
!             vsh = ehs(i,j,k) * sol / bmag / bmsf(i,j,k) *
!     .             ( (ehsy(i,j,k) * bdirsz(i,j,k) - 
!     .                ehsz(i,j,k) * bdirsy(i,j,k) ) * xnorms(i,j,k) -
!     .               (ehsx(i,j,k) * bdirsz(i,j,k) -
!     .                ehsz(i,j,k) * bdirsx(i,j,k) ) * ynorms(i,j,k) +
!     .               (ehsx(i,j,k) * bdirsy(i,j,k) -
!     .                ehsy(i,j,k) * bdirsx(i,j,k) ) * znorms(i,j,k) )

!             vps = 0.

!      need to check this

             vps = vexbp_phi(i,j,k) *
     .            ( vpsnx(i,j,k) * xnorms(i,j,k) +   
     .              vpsny(i,j,k) * ynorms(i,j,k) +   
     .              vpsnz(i,j,k) * znorms(i,j,k)   )  

             vph = vexbh_phi(i,j,k) *
     .            ( vhsnx(i,j,k) * xnorms(i,j,k) +   
     .              vhsny(i,j,k) * ynorms(i,j,k) +   
     .              vhsnz(i,j,k) * znorms(i,j,k)   )  
            
             vexbs_phi(i,j,k) = (vps + vph)
     .                            * tvexb0


!             if (abs(blats(i,j,k)) .gt. 72. ) vexbs_phi(i,j,k) = 0.

!             vexbs_phi(i,j,k) = 0.

           enddo
         enddo
       enddo

       do k = 1,nl
         do j = 1,nf
           vexbs_phi(nzp1,j,k) = vexbs_phi(nz,j,k)
         enddo
       enddo

      return
      end


!     ********************************************
!
!     current_jp_jphi
!
!     ********************************************

      subroutine current_jp_jphi

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc' 

!     jp and jphi currents
!     gravity not included

      do k = 1,nl
        do j = 1,nf
          do i = 1,nz
            ep        = eps(i,j,k)
            ephi      = ehs(i,j,k)
!            ep        = eph(i,j,k)
!            ephi      = ehp(i,j,k)
            boverc    = bmag * bms(i,j,k) / sol
            jp(i,j,k) = sigmap(i,j,k) * 
     .                  ( ep + boverc * vnphi(i,j,k) )
     .                + sigmah(i,j,k) *
     .                  ( -ephi + boverc * vnp(i,j,k) )
            jphi(i,j,k) = sigmap(i,j,k) * 
     .                    ( ephi - boverc * vnp(i,j,k) )
     .                  + sigmah(i,j,k) *
     .                    ( ep + boverc * vnphi(i,j,k) )
          enddo
        enddo
      enddo


      return
      end


!     ********************************************
!
!     PHI_NFP1_NLP1
!
!     ********************************************


! takes a variable in one grid (xu,zu):(nf,nl)
! and interpolates to another grid (x0,z0):(nx,ny)
! the nf,nf grid need not be orthogonal

      subroutine phi_nfp1_nlp1(brad,blon,phiout,phiin,phialt1,philon1)

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc' 

       real brad(nfp1,nlp1),blon(nfp1,nlp1),phiout(nfp1,nlp1)
       real phiin(nnx,nny),phialt1(nnx,nny),philon1(nnx,nny)

       do j = 1,nlp1
         do i = 1,nfp1
            phiout(i,j) = 0.
         enddo
       enddo

!       if ( taskid .eq. 24) then
!          do i =1,nnx
!            print *,taskid,i,philon1(i,40)
!          enddo
!          do i = 1,nlp1
!            print *,taskid,blon(41,i)
!          enddo
!       endif

       if ( taskid .ne. 1 .and. taskid .ne. numworkers ) then
         do k = 1,nlp1
           do j = nf,2,-1 
               kk = (taskid-1)*(nl-2) + (k-1)
               jj = j -1 
               phiout(j,k) = phiin(kk,jj)
           enddo
         enddo
       endif

       if ( taskid .eq. 1 ) then
         do k = 1,nlp1
           do j = nf,2,-1 
               kk = k - 1
               if ( k .eq. 1 ) kk = nnx - 1
               jj = j - 1 
               phiout(j,k) = phiin(kk,jj) 
           enddo
         enddo
       endif

       if ( taskid .eq. numworkers ) then
         do k = 1,nlp1
           do j = nf,2,-1 
               kk = (taskid - 1)*(nl-2) + (k-1)
               if ( k .eq. nlp1   ) kk = 2
               if ( k .eq. nlp1-1 ) kk = 1
               jj = j -1 
               phiout(j,k) = phiin(kk,jj)
           enddo
         enddo
       endif

       do k = 1,nlp1
         phiout(1,k)     = phiout(2,k)  
!         phiout(nf+1,k)  =  phiout(nf,k) * 
!     .                      (blatp(nz-1,nf+1,k) - 90.) /
!     .                      (blatp(nz-1,nf,k)   - 90.)

         slope           = (blatp(nz-1,nf+1,k) - blatp(nz-1,nf-1,k))/
     .                     (blatp(nz-1,nf,k)   - blatp(nz-1,nf-1,k))
         phiout(nf+1,k)  = phiout(nf-1,k) +
     .                     (phiout(nf,k) - phiout(nf-1,k)) * slope
 

       enddo

!       if ( taskid .eq. 1 ) then
!         open (9011,file='phi_inout.dat',form='unformatted')
!         write(9011) phiin,phiout
!         close(9011)
!       endif


       return
       end

!     ********************************************
!
!     PHI_NFP1_NL
!
!     ********************************************


! takes a variable in one grid (xu,zu):(nf,nl)
! and interpolates to another grid (x0,z0):(nx,ny)
! the nf,nf grid need not be orthogonal

      subroutine phi_nfp1_nl(brad,blon,phiout,phiin,phialt1,philon1)

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc' 

       real phiout(nfp1,nl)
       real phiin(nnx,nny)

       do j = 1,nl
         do i = 1,nfp1
            phiout(i,j) = 0.
         enddo
       enddo

       if ( taskid .ne. 1 .and. taskid .ne. numworkers ) then
         do k = 1,nl
           do j = nf,2,-1 
               kk = (taskid-1)*(nl-2) + (k-1)
               jj = j -1 
               phiout(j,k) = phiin(kk,jj)
           enddo
         enddo
       endif

       if ( taskid .eq. 1 ) then
         do k = 1,nl
           do j = nf,2,-1 
               kk = k - 1
               if ( k .eq. 1 ) kk = nnx - 1
               jj = j - 1 
               phiout(j,k) = phiin(kk,jj) 
           enddo
         enddo
       endif

       if ( taskid .eq. numworkers ) then
         do k = 1,nl
           do j = nf,2,-1 
               kk = (taskid - 1)*(nl-2) + (k-1)
               if ( k .eq. nl   ) kk = 1
               jj = j -1 
               phiout(j,k) = phiin(kk,jj)
           enddo
         enddo
       endif

       do k = 1,nl
         phiout(1,k)     = phiout(2,k)  
!         phiout(nf+1,k)  =  phiout(nf,k) * 
!     .                      (blatp(nz-1,nf+1,k) - 90.) /
!     .                      (blatp(nz-1,nf,k)   - 90.)

         slope           = (blatp(nz-1,nf+1,k) - blatp(nz-1,nf-1,k))/
     .                     (blatp(nz-1,nf,k)   - blatp(nz-1,nf-1,k))
         phiout(nf+1,k)  = phiout(nf-1,k) +
     .                     (phiout(nf,k) - phiout(nf-1,k)) * slope
 

       enddo

       return
       end



!     ********************************************
!
!     PHI_NF_NLP1
!
!     ********************************************


      subroutine phi_nf_nlp1(brad,blon,phiout,phiin,phialt1,philon1)

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc' 

      real brad(nf,nlp1),blon(nf,nlp1),phiout(nf,nlp1)
      real phiin(nnx,nny),phialt1(nnx,nny),philon1(nnx,nny)

       do j = 1,nlp1
         do i = 1,nf
            phiout(i,j) = 0.
         enddo
       enddo

!       if ( taskid .eq. 24) then
!          do i =1,nnx
!            print *,taskid,i,philon1(i,40)
!          enddo
!          do i = 1,nlp1
!            print *,taskid,blon(41,i)
!          enddo
!       endif

       if ( taskid .ne. 1 .and. taskid .ne. numworkers ) then
         do k = 1,nlp1
           do j = nf-1,2,-1 
               kk = (taskid-1)*(nl-2) + (k-1)
               jj = j - 1 
               phiout(j,k) = phiin(kk,jj)
           enddo
         enddo
       endif

       if ( taskid .eq. 1 ) then
         do k = 1,nlp1 
           do j = nf-1,2,-1 
               kk = k - 1
               if ( k .eq. 1 ) kk = nnx - 1
               jj = j - 1 
               phiout(j,k) = phiin(kk,jj) 
           enddo
         enddo
       endif

       if ( taskid .eq. numworkers ) then
         do k = 1,nlp1
           do j = nf-1,2,-1 
               kk = (taskid - 1)*(nl-2) + (k-1)
               if ( k .eq. nlp1   ) kk = 2
               if ( k .eq. nlp1-1 ) kk = 1
               jj = j -1 
               phiout(j,k) = phiin(kk,jj)
           enddo
         enddo
       endif

       do k = 1,nlp1
         phiout(nf,k)  = phiout(nf-1,k) * 
     .                   (blatp(nz-1,nf,k)   - 90.) /
     .                   (blatp(nz-1,nf-1,k) - 90.)
         phiout(1,k)   = phiout(2,k)
       enddo

       return
       end

!       **********************************
!       **********************************

!             POTPPHI

!       **********************************
!       **********************************

        subroutine potpphi(phi,phialt,philon,dphi,hrut,p_crit)

        include 'param3_mpi-1.98.inc'
        include 'com3_mpi-1.98.inc' 

!        parameter ( nyextra = 10, nnyt = nny + nyextra )

        real hipcp_pot(nnx,nnyt),hipcphi_pot(nnx,nnyt)
        real hidphig_pot(nnx,nnyt),hidphiv_pot(nnx,nnyt)
        real hidpg_pot(nnx,nnyt),hidpv_pot(nnx,nnyt)
        real hihcm_pot(nnx,nnyt)
        real hipc_pot(nnx,nnyt)
        real hihc_pot(nnx,nnyt)
        real hidv_pot(nnx,nnyt)

        real phi(nnx,nny),philon(nnx,nny)
!        real dpreal(nnyt),preal(nnyt)
        real dpreal(nnyt),preal(nnyt)
        real dbang(nnx,nnyt),blang(nnx)
        real*8 dphireal(nnx+1)
        real*8 dxij,dxip1j
        real*8 f11_lb(nnx+1),f11_ub(nnx+1)

        real*8 dphi(nnx+1,nnyt),si(nnx+1,nnyt),sih(nnx+1,nny)
        real*8 a1(nnx+1,nnyt),a2(nnx+1,nnyt),a3(nnx+1,nnyt)
        real*8 a4(nnx+1,nnyt),a5(nnx+1,nnyt)
        real*8 dphi0(nnx+1,nnyt)

        real ylonp(nnx),ylatp(nny)
        real zigm11(nnx,nny),zigm22(nnx,nny),zigm2(nnx,nny)
        real rim1(nnx,nny),rim2(nnx,nny)

        real vexbpphi(nnx,nny),vexbhphi(nnx,nny)
        real p_crit(nnx-1),phivs(nnx,nnyt)

        real phi_weimer(nnx,nny)

        data idpinter / 1 /
        data ipcrit   / 1 /

        if ( idpinter .eq. 1 ) then

          nzh  = nz / 2 
      
          do j = 1,nny
            jj         = j + nyextra
            dpreal(jj) = ppt(nzh,j+1,1) - ppt(nzh,j,1)
            preal(jj)  = ppt(nzh,j,1)
!            print *,'j,jj,preal',j,jj,preal(jj)
          enddo

          do j = nyextra,1,-1
            dpreal(j) = dpreal(j+1) * 1.4
!            tot = tot + dpreal(j)
          enddo

          do i = 1,nnx+1
           dphireal(i) = (blonp0t(i+1)-blonp0t(i))*pie/180.
          enddo

! define blang
! dbang (new dy) 

          do i = 1,nnx-1
            do j = 1,nny
              jj          = j + nyextra
!              dbang(i,jj) = ( blatpt(nzp1,j+1,i) -
!     .                        blatpt(nzp1,j,i)    ) * pie / 180.
              blang(i)    = -blonpt(nz/2+1,j,i) * pie / 180.
            enddo
          enddo

          blang(nnx) = blang(1)

!          do j = 1,nny
!            jj            = j + nyextra
!            dbang(nnx,jj) = dbang(1,jj)
!          enddo

!          do i = 1,nnx-1
!            do j = nyextra,1,-1
!              dbang(i,j) = ( blatpt(nzp1,1,i) * pie / 180.0 ) / nyextra
!            enddo
!          enddo

!          do j = nyextra,1,-1
!            dbang(nnx,j) = dbang(1,j)
!          enddo

!        initialize dphi and dphi0

         do j = 1,nnyt
           do i = 1,nnx+1
             dphi (i,j) = 0.
             dphi0(i,j) = 0.
           enddo
         enddo

         if (restart) then
           open(1232,file='dphi.rst',form='unformatted')
           read(1232) dphi0
           close(1232)
         endif

         idpinter = 0

        endif


!       set up conductances and driver for potential equation
!       zero-gradient in phi (x); zero-gradient in p (y)
!       note: transpose variables

        do j = 1,nny
          jj   = j + nyextra
          do i = 2,nnx-1
            hipcp_pot(i,jj)    = 0.25 * ( hipcpt(j,i-1)   + 
     .                                   hipcpt(j,i)     +
     .                                   hipcpt(j+1,i-1) + 
     .                                   hipcpt(j+1,i)     )
            hihcm_pot(i,jj)    = 0.25 * ( hihcmt(j,i-1)   + 
     .                                   hihcmt(j,i)     +
     .                                   hihcmt(j+1,i-1) + 
     .                                   hihcmt(j+1,i)     ) 
            hipcphi_pot(i,jj)  = 0.25 * ( hipcphit(j,i-1)   + 
     .                                   hipcphit(j,i)     +
     .                                   hipcphit(j+1,i-1) + 
     .                                   hipcphit(j+1,i)     )
            hidphig_pot(i,jj)  = 0.25 * ( hidphigt(j,i-1)   + 
     .                                   hidphigt(j,i)     +
     .                                   hidphigt(j+1,i-1) + 
     .                                   hidphigt(j+1,i)     )
            hidpg_pot(i,jj)    = 0.25 * ( hidpgt(j,i-1)   + 
     .                                   hidpgt(j,i)     +
     .                                   hidpgt(j+1,i-1) + 
     .                                   hidpgt(j+1,i)     )
            hidphiv_pot(i,jj)  = 0.25 * ( hidphivt(j,i-1)   + 
     .                                   hidphivt(j,i)     +
     .                                   hidphivt(j+1,i-1) + 
     .                                   hidphivt(j+1,i)     )
            hidpv_pot(i,jj)    = 0.25 * ( hidpvt(j,i-1)   + 
     .                                   hidpvt(j,i)     +
     .                                   hidpvt(j+1,i-1) + 
     .                                   hidpvt(j+1,i)     )
            hipc_pot(i,jj)     = 0.25 * ( hipct(j,i-1)   + 
     .                                   hipct(j,i)     +
     .                                   hipct(j+1,i-1) + 
     .                                   hipct(j+1,i)     )
            hihc_pot(i,jj)     = 0.25 * ( hihct(j,i-1)   + 
     .                                   hihct(j,i)     +
     .                                   hihct(j+1,i-1) + 
     .                                   hihct(j+1,i)     ) 
            hidv_pot(i,jj)     = 0.25 * ( hidvt(j,i-1)   + 
     .                                   hidvt(j,i)     +
     .                                   hidvt(j+1,i-1) + 
     .                                   hidvt(j+1,i)     )
          enddo
        enddo

        do j = nyextra+1,nnyt
            hipcp_pot(1,j)   = 0.5 * (hipcp_pot(2,j)    +
     .                                hipcp_pot(nnx-1,j) )
            hihcm_pot(1,j)   = 0.5 * (hihcm_pot(2,j)    +
     .                                hihcm_pot(nnx-1,j) )
            hipcphi_pot(1,j) = 0.5 * (hipcphi_pot(2,j)  +
     .                                hipcphi_pot(nnx-1,j) )
            hidphig_pot(1,j) = 0.5 * ( hidphig_pot(2,j) +
     .                                hidphig_pot(nnx-1,j) )
            hidpg_pot(1,j)   = 0.5 * ( hidpg_pot(2,j) +
     .                                hidpg_pot(nnx-1,j) )
            hidphiv_pot(1,j) = 0.5 * ( hidphiv_pot(2,j) +
     .                                hidphiv_pot(nnx-1,j) )
            hidpv_pot(1,j)   = 0.5 * ( hidpv_pot(2,j) +
     .                                hidpv_pot(nnx-1,j) )
            hipc_pot(1,j)    = 0.5 * (hipc_pot(2,j)    +
     .                                hipc_pot(nnx-1,j) )
            hihc_pot(1,j)    = 0.5 * (hihc_pot(2,j)    +
     .                                hihc_pot(nnx-1,j) )
            hidv_pot(1,j)    = 0.5 * (hidv_pot(2,j)    +
     .                                hidv_pot(nnx-1,j) )
        enddo

        do j = nyextra+1,nnyt
            hipcp_pot(nnx,j)   = hipcp_pot(1,j)
            hihcm_pot(nnx,j)   = hihcm_pot(1,j)
            hipcphi_pot(nnx,j) = hipcphi_pot(1,j)
            hidphig_pot(nnx,j) = hidphig_pot(1,j)
            hidpg_pot(nnx,j)   = hidpg_pot(1,j)
            hidphiv_pot(nnx,j) = hidphiv_pot(1,j)
            hidpv_pot(nnx,j)   = hidpv_pot(1,j)
            hipc_pot(nnx,j)    = hipc_pot(1,j)
            hihc_pot(nnx,j)    = hihc_pot(1,j)
            hidv_pot(nnx,j)    = hidv_pot(1,j)
        enddo

        do j = nyextra,1,-1
          do i = 1,nnx
            hipcp_pot(i,j)   = hipcp_pot(i,j+1)   * .02
            hihcm_pot(i,j)   = hihcm_pot(i,j+1)   * .02
            hipcphi_pot(i,j) = hipcphi_pot(i,j+1) * .02
            hidphig_pot(i,j) = hidphig_pot(i,j+1) * .01
            hidpg_pot(i,j)   = hidpg_pot(i,j+1)   * .01
            hidphiv_pot(i,j) = hidphiv_pot(i,j+1) * .01
            hidpv_pot(i,j)   = hidpv_pot(i,j+1)   * .01
            hipc_pot(i,j)    = hipc_pot(i,j+1)    * .02
            hihc_pot(i,j)    = hihc_pot(i,j+1)    * .02
            hidv_pot(i,j)    = hidv_pot(i,j+1)    * .01
          enddo
        enddo

        k   = 0
        k0  = 10

        do while ( k .le. k0 )

!     div j = 0 differencing (for Pedersen and Hall via iteration)

      do j = 2,nnyt-1
        do i = 1,nnx

          im1  = i - 1
          ip1  = i + 1
          jm1  = j - 1
          jp1  = j + 1

          dxij   = dphireal(i)
          dxip1j = dphireal(ip1)
 
!          dyij    = dbang(i,j)
!          dyijp1  = dbang(i,jp1)

          dyij    = dpreal(j)
          dyijp1  = dpreal(jp1)

          if ( i .eq. 1   ) im1 = nnx - 1
          if ( i .eq. nnx ) ip1 = 2

          delx    = dxij + dxip1j
          dely    = dyij + dyijp1
 
          delx_inv  = 1. / delx
          dely_inv  = 1. / dely
          delxy_inv = delx_inv * dely_inv

          pcphimx  = (0.5 * ( hipcphi_pot(im1,j) 
     .               + hipcphi_pot(i,j)   ))
          pcphipx  = (0.5 * ( hipcphi_pot(i,j)   
     .               + hipcphi_pot(ip1,j) ))

          pcpmy  = (0.5 * ( hipcp_pot(i,jm1) + hipcp_pot(i,j)   ))
          pcppy  = (0.5 * ( hipcp_pot(i,j)   + hipcp_pot(i,jp1) ))

          if (hall) then

            hcmjp = 0.5 * ( hihcm_pot(i,jp1) + hihcm_pot(i,j)   ) 
            hcmjm = 0.5 * ( hihcm_pot(i,j)   + hihcm_pot(i,jm1) ) 
            dhcy  = hcmjp - hcmjm

            hcmip = 0.5 * ( hihcm_pot(i,j)   + hihcm_pot(ip1,j) ) 
            hcmim = 0.5 * ( hihcm_pot(im1,j) + hihcm_pot(i,j)   ) 
            dhcx  = hcmip - hcmim

          else

            dhcy   = 0.
            dhcx   = 0.
            hcmjp  = 0.
            hcmjm  = 0.
            hcmip  = 0.
            hcmim  = 0.

          endif

          fphivmx   = 0.5 * ( hidphiv_pot(im1,j) + hidphiv_pot(i,j)   )
          fphivpx   = 0.5 * ( hidphiv_pot(i,j)   + hidphiv_pot(ip1,j) )
          dfphivx   = fphivpx - fphivmx

          fphigmx   = 0.5 * ( hidphig_pot(im1,j) + hidphig_pot(i,j)   )
          fphigpx   = 0.5 * ( hidphig_pot(i,j)   + hidphig_pot(ip1,j) )
          dfphigx   = fphigpx - fphigmx


          fpvmy   = 0.5 * ( hidpv_pot(i,jm1) + hidpv_pot(i,j)   )
          fpvpy   = 0.5 * ( hidpv_pot(i,j)   + hidpv_pot(i,jp1) )
          dfpvy   = fpvpy - fpvmy

          fpgmy   = 0.5 * ( hidpg_pot(i,jm1) + hidpg_pot(i,j)   )
          fpgpy   = 0.5 * ( hidpg_pot(i,j)   + hidpg_pot(i,jp1) )
          dfpgy   = fpgpy - fpgmy

          a11    = 2. * delx_inv * pcphimx / dxij 
          a12    = delxy_inv * dhcy

          a1(i,j) = a11 + a12

          a21     = 2. * dely_inv * pcpmy / dyij 
          a22     = delxy_inv * dhcx

          a2(i,j) = a21 - a22

          a41     = 2. * dely_inv * pcppy / dyijp1
          a42     = delxy_inv * dhcx

          a4(i,j) = a41 + a42

          a51     = 2. * delx_inv * pcphipx / dxip1j
          a52     = delxy_inv * dhcy

          a5(i,j) = a51 - a52

          a3(i,j) = -a51 - a11 - a41 - a21 

          sip1jp1 = -delxy_inv * (hcmip - hcmjp) * dphi0(ip1,jp1)

          sip1jm1 =  delxy_inv * (hcmip - hcmjm) * dphi0(ip1,jm1)

          sim1jp1 =  delxy_inv * (hcmim - hcmjp) * dphi0(im1,jp1)

          sim1jm1 = -delxy_inv * (hcmim - hcmjm) * dphi0(im1,jm1)

          si(i,j) =   
     .                2. * dfphigx * delx_inv 
     .              + 2. * dfphivx * delx_inv 
     .              - 2. * dfpgy   * dely_inv  
     .              + 2. * dfpvy   * dely_inv  
     .              + sip1jp1 + sip1jm1 + sim1jp1 + sim1jm1

        enddo
      enddo

!     zero gradient in y

      do i = 1,nnx
        a1(i,1)   = a1(i,2)
        a2(i,1)   = a2(i,2)
        a3(i,1)   = a3(i,2)
        a4(i,1)   = a4(i,2)
        a5(i,1)   = a5(i,2)
        si(i,1)   = si(i,2)

        dpp          = dpreal(nnyt) / dpreal(nnyt-1)
!        dpp          = dbang(i,nnyt) / dbang(i,nnyt-1)

        a1(i,nnyt)   = dpp*(a1(i,nnyt-1)-a1(i,nnyt-2))+a1(i,nnyt-1)
        a2(i,nnyt)   = dpp*(a2(i,nnyt-1)-a2(i,nnyt-2))+a2(i,nnyt-1)
        a3(i,nnyt)   = dpp*(a3(i,nnyt-1)-a3(i,nnyt-2))+a3(i,nnyt-1)
        a4(i,nnyt)   = dpp*(a4(i,nnyt-1)-a4(i,nnyt-2))+a4(i,nnyt-1)
        a5(i,nnyt)   = dpp*(a5(i,nnyt-1)-a5(i,nnyt-2))+a5(i,nnyt-1)
        si(i,nnyt)   = dpp*(si(i,nnyt-1)-si(i,nnyt-2))+si(i,nnyt-1)

      enddo

      do j = 1,nnyt 
        a1(nnx+1,j) = a1(2,j)
        a2(nnx+1,j) = a2(2,j)
        a3(nnx+1,j) = a3(2,j)
        a4(nnx+1,j) = a4(2,j)
        a5(nnx+1,j) = a5(2,j)
        si(nnx+1,j) = si(2,j)
      enddo

      do i = 2,nnx
         f11_lb(i)  = 0.
      enddo

      f11_lb(1)     = f11_lb(nnx)
      f11_lb(nnx+1) = f11_lb(2)

      phi90 = 0.

      do i = 2,nnx
!       f11_ub(i)  = 0.

       dbangi90     = (blatpt(nz-1,nf-2,1)-90.) 
       dbangif      = (blatpt(nz-1,nf-2,1) - blatpt(nz-1,nf-1,1))
       dbangf90     = (blatpt(nz-1,nf-1,1)-90.) 

       dphi(i,nnyt) = ( dphi(i,nnyt-1) * dbangf90 +
     .                  phi90          * dbangif  ) / dbangi90 

      enddo


      f11_ub(1)     = f11_ub(nnx)
      f11_ub(nnx+1) = f11_ub(2)

      if ( lmadala ) then
        call madala_sevp(a1,a2,a5,a4,a3,si,dphi,f11_lb,f11_ub)
      else
        do j = 1,nnyt
          do i = 1,nnx+1
            dphi(i,j) = 0.
          enddo
        enddo
      endif

      i00 = nnx/2
      j00 = nnyt/2

       err0 = abs(dphi0(i00,j00)-dphi(i00,j00))/
     .        abs(dphi0(i00,j00)+1.e-5)
       k = k + 1
       print *,k,err0
 
       tol_phi  = 2.e-2

       if ( err0 .le. tol_phi ) k0 = -10

      do j = 1,nnyt
        do i = 1,nnx+1
          dphi0(i,j) = dphi(i,j)
        enddo
      enddo

      enddo

!     parameters for volland/stern/macilwane potential

      fkp  = 6.

      akp_coef_i = 11.25
      akp_coef_f = 45.00

      if ( hrut .lt. storm_ti ) akp_coef = akp_coef_i
      if ( hrut .gt. storm_tf ) akp_coef = akp_coef_f
      if ( hrut .ge. storm_ti .and. hrut .le. storm_tf )
     .     akp_coef = akp_coef_i +
     .                (akp_coef_f - akp_coef_i) *
     .                (hrut - storm_ti)/(storm_tf - storm_ti)

      akp  = akp_coef / ( (1.-0.159*fkp+0.0093*fkp*fkp) ** 3 )

!      akp  = 45./ ( (1.-0.159*fkp+0.0093*fkp*fkp) ** 3 )
!      akp  = 11.25/ ( (1.-0.159*fkp+0.0093*fkp*fkp) ** 3 )

!        print *,'fkp,akp',fkp,akp

!     only for non-corotating code

!     rotate phi so that potential is 
!     aligned with local time midnight/noon 
!       angut  in degrees
!       angutr in radians

      angrot = 360. - plon * 180. / pie
      hr24   = mod (hrut,24.) 
      angut  = hr24 * 2. * pie / 24. * rtod - angrot 
      if ( angut .gt. 360. ) angut = angut - 360. 
      if ( angut .lt. 0.   ) angut = angut + 360. 
      angutr = angut * pie / 180.

      if ( lcr ) then
        angutr  = 0.
        angut   = 0.                
      endif

      if ( lweimer ) then
        call weimer(phi_weimer,angut,hrut)
      else 
        do j = nyextra+1,nnyt
          do i = 1,nnx
            jj = j - nyextra
            phi_weimer(i,jj) = 0.
          enddo
        enddo
      endif

      do j = nyextra+1,nnyt

        do i = 1,nnx
          jj = j - nyextra
          phivs(i,jj)   = 0.
          if ( lvs )
     .      phivs(i,jj)   = -akp * preal(j) * preal(j) * 
     .                      sin(blang(i)-angutr)
     .                    / 300.
          phicorot    = 0. 
          if ( lcr ) 
     .     phicorot    = -92.e3 / 300. / preal(j) * bmag /.31
           phi(i,jj)   = dphi(i,j) + phicorot + phivs(i,jj) +
     .                   phi_weimer(i,jj)
        enddo
      enddo

! find p_crit (i.e., open/closed field line boundary)
! here, boundary between corotation and convective potential
! along dusk longitude

! not used 

!      if ( ipcrit .eq. 1 ) then

!        nnx_dusk = 3*nlt/4 + 1

!        do j = 2,nny
!          if (phi(nnx_dusk,j) .ge. phi(nnx_dusk,j-1)) j_crit = j
!        enddo

!        phi_crit = phi(nnx_dusk,j_crit)
!        print *,'phi_crit',phi_crit

!        do j = 2,nny
!          do i = 0,(nnx-1)/2
!            if ( phi(i,j) .le. phi_crit ) then
!              p_crit(i)     = preal(j+nyextra)
!            endif
!          enddo
!        enddo

!        do j = nny,1,-1
!          do i = (nnx-1)/2+1,nnx-1
!            if ( phi(i,j) .ge. phi_crit ) then
!              p_crit(i)     = preal(j+nyextra)
!            endif
!          enddo
!        enddo

!        ipcrit = 0

!       endif

!      do mm = 1,12
!        call smoothy(phi)
!      enddo

      print *,phi(40,35)

! calculate grid that phi is on

!      do i = 1,nnx-1
!        ylonp(i) = blonpt(nz/2,1,i)
!      enddo
!      ylonp(nnx) = 360.

!      do j = 1,nny
!        ylatp(j) = blatpt(nz-1,j,1)
!      enddo

!      do j = nyextra + 1, nny + nyextra
!        jj = j - nyextra
!        do i = 1,nnx
!          zigm22(i,jj) = hipcp_pot(i,j)
!          zigm2(i,jj)  = hihcm_pot(i,j)
!          zigm11(i,jj) = hipcphi_pot(i,j)
!          rim1(i,jj)   = hidphiv_pot(i,j)
!          rim2(i,jj)   = hidpv_pot(i,j)
!        enddo
!      enddo

!      open(201,file='lonlatphiu.dat',form='unformatted')
!      write(201) ylonp,ylatp,phi,zigm22,zigm2,zigm11,rim1,rim2
!      close(201)

      return

      end


********************************************************
********************************************************
********************************************************
********************************************************
********************************************************

        subroutine weimer(phi_weimer,angut,hrut)

        include 'param3_mpi-1.98.inc'
        include 'com3_mpi-1.98.inc' 

        real*8 phi_weimer_real(nfp1,nlt+1)
        real phi_weimer_interp(nfp1,nlt+1)
        real phi_weimer(nnx,nny)
        real*8 hrutw1,hrutw2

        data iread_weimer / 1 /

!       read in weimer at first time

        if ( iread_weimer .eq. 1 ) then
          open(810,file='phi_weimer.inp',form='unformatted')
          read(810) hrutw1
          if ( restart ) then
            open (881,file='nweimer.rst',form='formatted')
            read(881,*) nweimer
            print *,'nweimer',nweimer
            do i = 1,nweimer-1
              read(810) phi_weimer_real
              read(810) hrutw2
            enddo
          else
            read(810) phi_weimer_real
            read(810) hrutw2
            print *,'hrutw2 = ',hrut,hrutw2
            nweimer       = 1
            print *,'nweimer',nweimer
          endif
          iread_weimer  = 0
        endif

!       read in weimer for subsequent time steps

        if ( hrut .ge. hrutw2 ) then
          read(810) phi_weimer_real
          read(810) hrutw2
          print *,'hrutw2 = ',hrut,hrutw2
          nweimer         = nweimer + 1
        endif

!       save nweimer for restart

        open (881,file='nweimer.rst',form='formatted')
        write(881,*) nweimer
        close(881)        

!       requires sami3 and weimer potential to
!       have the same latitude
!       (from sami3_rcm interp_rcm subroutine)

        do j = 1,nny
          do i = 1,nnx
             phi_weimer(i,j) = 0.
          enddo
        enddo

        dlon = 360. / float(nlt)

        do k = 1,nlt+1
          do j = nfp1,1,-1
            thlon = mod(blonp0t(k+1) + angut,360.)
            if ( thlon .lt. 0. ) thlon = thlon + 360.
            nlon = thlon/dlon + 1
            dnlon = thlon/dlon - nlon + 1
            phi_weimer_interp(j,k) = 
     &                 + dnlon * phi_weimer_real(j,nlon+1)
     &                 + (1. - dnlon) * phi_weimer_real(j,nlon)
!            if (j.eq.121) then
!              print *,'k121',k,thlon,nlon,dnlon,
!     .                  phi_weimer_real(j,nlon),phi_weimer_interp(j,k)
!            endif
!            if (j.eq.122) then
!              print *,'k122',k,thlon,nlon,dnlon,
!     .                  phi_weimer_real(j,nlon),phi_weimer_interp(j,k)
!            endif
          enddo
        enddo

!       here nnx = nlt + 1
!            nny = nf  - 1

        do j = 1,nny
          do i = 1,nnx
             phi_weimer(i,j) = phi_weimer_interp(j,i)
          enddo
        enddo

!       try to fix weimer potential problem with
!       simple linear interpolation

          j    = nny - 1
          do i = 1,nnx
             phi_weimer(i,j) = 0.5 * ( phi_weimer(i,j-1) +
     .                                 phi_weimer(i,j+1)   )
          enddo

        return
        end

********************************************************
********************************************************
********************************************************
********************************************************
********************************************************

        subroutine smoothx(f) 
 
        include 'param3_mpi-1.98.inc' 

        parameter ( nnxp2 = nnx + 2, nnyp2 = nny + 2 ) 
        parameter ( nnxp1 = nnx + 1, nnyp1 = nny + 1 ) 
         
        real f(nnx,nny),f0(nnx+2,nny+2) 
 
        u12 = 1. 
 
          do j = 1,nny 
            do i = 1,nnx 
              f0(i+1,j+1) = f(i,j) 
            enddo 
          enddo 
 
!       zero-gradient in x 
 
!        do j = 2,nnyp1 
!          f0(1,j)     = f0(2,j) 
!          f0(nnx+2,j)  = f0(nnx+1,j) 
!        enddo 

!       periodic in x 
 
        do j = 2,nnyp1 
          f0(1,j)     = f0(nnx,j) 
          f0(nnx+2,j)  = f0(3,j) 
        enddo 
 
!       zero gradient in y 
 
        do i = 1,nnxp2 
          f0(i,1)    = f0(i,2) 
          f0(i,nnyp2) = f0(i,nnyp1) 
        enddo       
 
!       sweep in x (1/2) 
    
          do j = 2,nnyp1 
            do i = 1,nnxp1 
              f0 (i,j) = 0.5 * ( f0(i,j) + u12*f0(i+1,j) ) 
            enddo 
          enddo 
 
          do j = 2,nnyp1 
            do i = nnxp2,2,-1 
              f0 (i,j) = 0.5 * ( f0(i,j) + u12*f0(i-1,j) ) 
            enddo 
          enddo 
 
!       now get f 
 
          do j = 1,nny 
            do i = 1,nnx 
              f(i,j)  = f0(i+1,j+1) 
            enddo 
          enddo 
 
        return  
        end 

********************************************************
********************************************************
********************************************************
********************************************************
********************************************************


        subroutine smoothy(f) 
 
        include 'param3_mpi-1.98.inc' 

        parameter ( nnxp2 = nnx + 2, nnyp2 = nny + 2 ) 
        parameter ( nnxp1 = nnx + 1, nnyp1 = nny + 1 ) 
         
        real f(nnx,nny),f0(nnx+2,nny+2) 
 
        u12 = 1. 
 
        do j = 1,nny 
          do i = 1,nnx 
            f0(i+1,j+1) = f(i,j) 
          enddo 
        enddo 
        
!       zero-gradient in x 
 
!        do j = 2,nnyp1 
!          f0(1,j)     = f0(2,j) 
!          f0(nnx+2,j)  = f0(nnx+1,j) 
!        enddo 
 
!       periodic in x 
 
        do j = 2,nnyp1 
          f0(1,j)     = f0(nnx,j) 
          f0(nnx+2,j)  = f0(3,j) 
        enddo 

!       zero gradient in y 
 
        do i = 1,nnxp2 
          f0(i,1)    = f0(i,2) 
          f0(i,nnyp2) = f0(i,nnyp1) 
        enddo       
 
!       sweep in y (1/2) 
 
          do j = 1,nnyp1 
            do i = 1,nnxp2 
              f0 (i,j) = 0.5 * ( f0(i,j) + u12*f0(i,j+1) ) 
            enddo 
          enddo 
 
          do j = nnyp2,2,-1 
            do i = 1,nnxp2 
              f0 (i,j) = 0.5 * ( f0(i,j) + u12*f0(i,j-1) ) 
            enddo 
          enddo 
 
!       now get f 
 
          do j = 1,nny 
            do i = 1,nnx 
              f(i,j)  = f0(i+1,j+1) 
            enddo 
          enddo 
 
        return 
        end 
 


! All variables made lower case (JH) 4/26/05

C*********************************************************************C
C*                                                                   *C
C*  chapman.for                                                      *C
C*                                                                   *C
C*  Written by:  David L. Huestis, Molecular Physics Laboratory      *C
C*                                                                   *C
C*  Copyright (c) 2000  SRI International                            *C
C*  All Rights Reserved                                              *C
C*                                                                   *C
C*  This software is provided on an as is basis; without any         *C
C*  warranty; without the implied warranty of merchantability or     *C
C*  fitness for a particular purpose.                                *C
C*                                                                   *C
C*********************************************************************C
C*
C*	To calculate the Chapman Function, Ch(X,chi0), the column 
C*	depth of an exponential atmosphere integrated along a line 
C*	from a given point to the sun, divided by the column depth for 
C*	a vertical sun.
C*
C*  USAGE:
C*
C*	  z = altitude above the surface
C*	  R = radius of the planet
C*	  H = atmospheric scale height
C*
C*	  X = (R+z)/H
C*	  chi0 = solar zenith angle (in degrees)
C*
C*	  implicit real*4(a-h,o-z)
C*	  depth = atm_chapman(X,chi0)	! analytical
C*	  depth = atm_chap_num(X,chi0)	! numerical (chi0 .le. 90)
C*
C*	  implicit real*8(a-h,o-z)
C*	  depth = atm8_chapman(X,chi0)	! analytical
C*	  depth = atm8_chap_num(X,chi0)	! numerical (chi0 .le. 90)
C*
C*  PERFORMANCE:
C*
C*	Compiled and linked using Microsoft FORTRAN 5.1, and executed 
C*	in MS-DOS mode under Windows 95 on a 160 MHz PC.
C*
C*    TIMING (in microseconds, typical)
C*
C*	  120	atm_chapman and atm8_chapman for X .lt. 36
C*	   25	atm_chapman and atm8_chapman for X .ge. 36
C*	  500	atm_chap_num
C*	 5000	atm8_chap_num
C*
C*    ACCURACY (maximum relative error, 0.le.chi0.le.90, 1.le.X.le.820)
C*
C*	6.0E-7	atm_chapman and atm8_chapman for X .lt. 60
C*	1.5E-7	atm_chapman and atm8_chapman for X .ge. 60
C*	6.0E-8	atm_chap_num
C*	1.E-15	atm8_chap_num (convergence test)
C*
C*    CODING
C*
C*	No claims are made that the code is optimized for speed, 
C*	accuracy, or compactness.  The principal objectives were 
C*
C*	  (1) Robustness with respect to argument values
C*	  (2) Rigorous mathematical derivation and error control
C*	  (3) Maximal use of "well known" mathematical functions
C*	  (4) Ease of readability and mapping of theory to coding
C*
C*	The real*8 accuracy could be improved with more accurate 
C*	representations of E1(), erfc(), I0(), I1(), K0(), K1().
C*
C*	In the course of development, many representations and 
C*	approximations of the Chapman Function were attempted that 
C*	failed to be robustly extendable to machine-precision.
C*
C*  INTERNET ACCESS:
C*
C*	Source: http://www-mpl.sri.com/software/chapman/chapman.html
C*	Author: mailto:david.huestis@sri.com
C*	        http://www-mpl.sri.com/bios/Huestis-DL.html
C*
C*  EDIT HISTORY:
C*
C*	01/22/2000 DLH	First complete documentation
C*
C*	01/15/2000 DLH	First complete version of chapman.for
C*
C**********************************************************************
C*
C*  THEORY:
C*
C*    INTRODUCTION
C*
C*	    This computer code models the absorption of solar radiation 
C*	by an atmosphere that depends exponentionally on altitude.  In 
C*	specific we calculate the effective column depth of a species 
C*	of local density, n(z), from a point at a given altitude, z0, 
C*	to the sun at a given solar zenith angle, chi0.  Following Rees 
C*	[Re89, Section 2.2] we write the column depth for chi0 .le. 90 
C*	degrees as
C*
C*   (A)  N(z0,chi0) = int{z=z0,infinity} 
C*	     [ n(z)/sqrt( 1 - ( sin(chi0) * (R+z0) / (R+z) ) **2 ) dz ]
C*
C*	where R is the radius of the solid planet (e.g. Earth).  For 
C*	chi0 .gt. 90 degrees we write
C*
C*	  N(z0,chi0) = 2*N(zs,90) - N(z0,180-chi0)
C*
C*	where zs = (R+z0)*sin(chi0)-R is the tangent height.
C*
C*	    For an exponential atmosphere, with
C*
C*	  n(z) = n(z0) * exp(-(z-z0)/H)
C*
C*	with a constant scale height, H, the column depth can be 
C*	represented by the Chapman function, Ch(X,chi0), named after 
C*	the author of the first quantitative mathematical investigation 
C*	[Ch31b] trough the relation
C*
C*	  N(z0,chi0) = H * n(z0) * Ch(X,chi0)
C*
C*	where X = (R+z0)/H is a dimensionless measure of the radius 
C*	of curvature, with values from about 300 to 1300 on Earth.
C*
C*
C*    APPROACH
C*
C*	    We provide function entry points for very stable and 
C*	reasonably efficient evaluation of Ch(X,chi0) with full 
C*	single-precision accuracy (.le. 6.0E-7 relative) for a wide 
C*	range of parameters.  A 15-digit-accurate double precision 
C*	numerical integration routine is also provided.
C*
C*	    Below we will develop (1) a compact asymptotic expansion of 
C*	good accuracy for moderately large values of X (.gt. 36) and all 
C*	values of chi0, (2) an efficient numerical integral for 
C*	all values of X and chi0, and (3) an explicit analytical 
C*	representation, valid for all values of X and chi0, based 
C*	the differential equation satisfied by Ch(X,chi0).
C*
C*	    All three of these represent new research results as well 
C*	as significant computational improvements over the previous 
C*	litearture, much of which is cited below.
C*
C*
C*    CHANGES OF THE VARIABLE OF INTEGRATION
C*
C*	Substituting y = (R+z)/(R+z0) - 1 we find
C*
C*   (B)  Ch(X,chi0) = X * int{y=0,infinity}
C*	     [ exp(-X*y) / sqrt( 1 - ( sin(chi0) / (1+y) )**2 ) dy ]
C*
C*	The futher substitutions s = (1+y)/sin(chi0), s0 = 1/sin(chi0) 
C*	give
C*
C*   (C)  Ch(X,chi0) = X*sin(chi0) * int{s=s0,infinity}
C*	     [ exp(X*(1-sin(chi0)*s)) * s / sqrt(s**2-1) ds ]
C*
C*	From this equation we can establish that
C*
C*	  Ch(X,90) = X*exp(X)*K1(X)
C*
C*	[AS64, Equations 9.6.23 and 9.6.27].  If we now substitute
C*	s = 1/sin(lambda) we obtain
C*
C*   (D)  Ch(X,chi0) = X*sin(chi0) * int{lambda=0,chi0} 
C*	    [ exp(X*(1-sin(chi0)*csc(lambda))) * csc(lambda)**2 dlambda]
C*
C*	which is the same as Chapman's original formulation [Ch31b, p486,
C*	eqn (10)].  If we first expand the square root in (B)
C*
C*	  1/sqrt(1-q) = 1 + q/( sqrt(1-q)*(1+sqrt(1-q)) )
C*
C*	with q = ( sin(chi0) / (1+y) )**2 = sin(lambda)**2, we obtain 
C*	a new form of (D) without numerical sigularities and simple 
C*	convergence to Ch(0,chi0) = Ch(X,0) = 1
C*
C*   (E)  Ch(X,chi0) = 1 + X*sin(chi0) * int{lambda=0,chi0} 
C*	    [ exp(X*(1-sin(chi0)*csc(lambda))) 
C*		/ (1 + cos(lambda) ) dlambda ]
C*
C*	Alternatively, we may substitute t**2 = y + t0**2, 
C*	into Equation (B), with t0**2 = 1-sin(chi0), finding
C*
C*   (F)  Ch(X,chi0) = X * int{s=t0,infinity} 
C*	    [ exp(-X*(t**2-t0**2)) * f(t,chi0) dt ]
C* 
C*	where
C*
C*	  f(t,chi0) = (t**2 + sin(chi0)) / sqrt(t**2+2*sin(chi0))
C*
C*	  f(t,chi0) = (t**2-t0**2+1)/sqrt(t**2-t0**2+1+sin(chi0))
C*
C*	    Below we will use Equation (F) above to develop a
C*	compact asymptotic expansion of good accuracy for moderately 
C*	large values of X (.gt. 36) and all values of chi0, Equation (E) 
C*	to develop an efficient numerical integral for Ch(X,chi0) for 
C*	all values of X and chi0, and Equation (C) to derive an explicit 
C*	analytical representation, valid for all values of X and chi0,  
C*	based on the differential equation satisfied by Ch(X,chi0).
C*
C*    atm_chapman(X,chi0) and atm8_chapman(X,chi0)
C*
C*	These routines return real*4 and real*8 values of Ch(X,chi0)
C*	selecting the asymptotic expansion or differential equation 
C*	approaches, depending on the value of X.  These routines also 
C*	handle the case of chi0 .gt. 90 degrees.
C*
C*    atm_chap_num(X,chi0) and atm8_chap_num(X,chi0)
C*
C*	These routines return real*4 and real*8 values of Ch(X,chi0)
C*	evaluated numerically.  They are both more accurate than the 
C*	corresponding atm*_chapman() functions, but take significantly 
C*	more CPU time.
C*
C*
C*    ASYMPTOTIC EXPANSION
C*
C*	From Equation (F) we expand, with t0**2 = 1-sin(chi0), 
C*
C*	  f(t,chi0) = sum{n=0,3} [ C(n,chi0) * (t**2-t0**2)**n ]
C*
C*	The function atm8_chap_asy(X,chi0) evaluates integrals of the 
C*	form
C*
C*	  int{t=t0,infinity} [exp(-X*(t**2-t0**2))*(t**2-t0**2)**n dt]
C*
C*	in terms of incomplete gamma functions, and sums them to 
C*	compute Ch(X,chi0).  For large values of X, this results in an 
C*	asymptotic expansion in negative powers of X, with coefficients 
C*	that are stable for all values of chi0.
C*
C*	In contrast, the asymptotic expansions of Chapman [Ch31b, 
C*	p488, Equation (22) and p490, Equation (38)], Hulburt [He39], 
C*	and Swider [Sw64, p777, Equation (43)] use negative powers of 
C*	X*cos(chi0)**2 or X*sin(chi0), and are accurate only for 
C*	small values or large values of chi0, respectively.
C*
C*	Taking only the first term in the present expansion gives the 
C*	simple formula
C*
C*	  Ch(X,chi0) = sqrt(pi*X/(1+sin(chi0))) * exp(X*(1-sin(chi0)))
C*		* erfc( sqrt(X*(1-sin(chi0))) )
C*
C*	This is slightly more accurate than the semiempirical 
C*	formula of Fitzmaurice [Fi64, Equation (3)], and sightly less 
C*	accurate than that of Swider [Sw64, p780, Equation (52), 
C*	corrected in SG69].
C*
C*
C*    NUMERICAL INTEGRATION
C*
C*	We are integrating
C*
C*   (E)  Ch(X,chi0) = 1 + X*sin(chi0) * int{lambda=0,chi0} 
C*	    [ exp(X*(1-sin(chi0)*csc(lambda))) 
C*		/ ( 1 + cos(lambda) ) dlambda ]
C*
C*	The integrand is numerically very smooth, and rapidly varying 
C*	only near lambda = 0.  For X .ne. 0 we choose the lower limit 
C*	of numerical integration such that the integrand is 
C*	exponentially small, 7.0E-13 (3.0E-20 for real*8).  The domain 
C*	of integration is divided into 64 equal intervals (6000 for 
C*	real*8), and integrated numerically using the 9-point closed 
C*	Newton-Cotes formula from Hildebrand [Hi56a, page 75, Equation
C*	(3.5.17)].
C*
C*
C*    INHOMOGENOUS DIFFERENTIAL EQUATION
C*
C*	    The function atm8_chap_deq(X,chi0) calculates Ch(X,chi0), 
C*	based on Equation (C) above, using the inhomogeneous 
C*	Bessel's equation as described below.  Consider the function 
C*
C*	  Z(Q) = int{s=s0,infinity} [ exp(-Q*s) / sqrt(s**2-1) ds ]
C*
C*	Differentiating with respect to Q we find that 
C*
C*	  Ch(X,chi0) = - Q * exp(X) * d/dQ [ Z(Q) ]
C*
C*	with Q = X*sin(chi0), s0 = 1/sin(chi0).  Differentiating 
C*	inside the integral, we find that
C*
C*	  Z"(Q) + Z'(Q)/Q - Z(Q) = sqrt(s0**2-1) * exp(-Q*s0) / Q
C*
C*	giving us an inhomogeneous modified Bessel's equation of order 
C*	zero.  Following Rabenstein [Ra66, pp43-45,149] the solution 
C*	of this equation can be written as
C*
C*	  Z(Q) = A*I0(Q) + B*K0(Q) - sqrt(s0**2-1) 
C*	         * int{t=Q,infinity} [ exp(-t*s0) 
C*		   * ( I0(Q)*K0(t) - I0(t)*K0(Q) ) dt ] 
C*
C*	with coefficients A and B to be determined by matching 
C*	boundary conditions.
C*
C*	    Differentiating with respect to Q we obtain
C*
C*	  Ch(X,chi0) = X*sin(chi0)*exp(X)*( 
C*		- A*I1(X*sin(chi0)) + B*K1(X*sin(chi0)) 
C*		+ cos(chi0) * int{y=X,infinity} [ exp(-y) 
C*		  * ( I1(X*sin(chi0))*K0(y*sin(chi0))
C*		    + K1(X*sin(chi0))*I0(y*sin(chi0)) ) dy ] )
C*
C*	Applying the boundary condition Ch(X,0) = 1 requires that 
C*	B = 0.  Similarly, the requirement that Ch(X,chi0) approach 
C*	the finite value of sec(chi0) as X approaches infinity [Ch31b, 
C*	p486, Equation (12)] implies A = 0.  Thus we have
C*
C*	  Ch(X,chi0) = X*sin(chi0)*cos(chi0)*exp(X)*
C*		int{y=X,infinity} [ exp(-y) 
C*		  * ( I1(X*sin(chi0))*K0(y*sin(chi0))
C*		    + K1(X*sin(chi0))*I0(y*sin(chi0)) ) dy ]
C*
C*	The function atm8_chap_deq(X,chi0) evaluates this expression.
C*	Since explicit approximations are available for I1(z) and K1(z),
C*	the remaining challenge is evaluation of the integrals
C*
C*	  int{y=X,infinity} [ exp(-y) I0(y*sin(chi0)) dy ]
C*
C*	and
C*
C*	  int{y=X,infinity} [ exp(-y) K0(y*sin(chi0)) dy ]
C*
C*	which are accomplished by term-by-term integration of ascending
C*	and descending power series expansions of I0(z) and K0(z).
C*
C*  REFERENCES:
C*
C*	AS64	M. Abramowitz and I. A. Stegun, "Handbook of 
C*		Mathematical Functions," NBS AMS 55 (USGPO, 
C*		Washington, DC, June 1964, 9th printing, November 1970).
C*
C*	Ch31b	S. Chapman, "The Absorption and Dissociative or
C*		Ionizing Effect of Monochromatic Radiation in an
C*		Atmosphere on a Rotating Earth: Part II. Grazing
C*		Incidence," Proc. Phys. Soc. (London), _43_, 483-501 
C*		(1931).
C*
C*	Fi64	J. A. Fitzmaurice, "Simplfication of the Chapman
C*		Function for Atmospheric Attenuation," Appl. Opt. _3_,
C*		640 (1964).
C*
C*	Hi56a	F. B. Hildebrand, "Introduction to Numerical
C*		Analysis," (McGraw-Hill, New York, 1956).
C*
C*	Hu39	E. O. Hulburt, "The E Region of the Ionosphere," 
C*		Phys. Rev. _55_, 639-645 (1939).
C*
C*	PFT86	W. H. Press, B. P. Flannery, S. A. Teukolsky, and 
C*		W. T. Vetterling, "Numerical Recipes," (Cambridge, 
C*		1986).
C*
C*	Ra66	A. L. Rabenstein, "Introduction to Ordinary
C*		Differential Equations," (Academic, NY, 1966).
C*
C*	Re89	M. H. Rees, "Physics and Chemistry of the Upper
C*		Atmosphere," (Cambridge, 1989).
C*
C*	SG69	W. Swider, Jr., and M. E. Gardner, "On the Accuracy 
C*		of Chapman Function Approximations," Appl. Opt. _8_,
C*		725 (1969).
C*
C*	Sw64	W. Swider, Jr., "The Determination of the Optical 
C*		Depth at Large Solar Zenith Angles," Planet. Space 
C*		Sci. _12_, 761-782 (1964).
C
C  ####################################################################
C
C	Chapman function calculated by various methods
C
C	  Ch(X,chi0) = atm_chapman(X,chi0)   : real*4 entry
C	  Ch(X,chi0) = atm8_chapman(X,chi0)  : real*8 entry
C
C	Internal service routines - user should not call, except for
C	testing.
C
C	  Ch(X,chi0) = atm8_chap_asy(X,chi0) : asymptotic expansion
C	  Ch(X,chi0) = atm8_chap_deq(X,chi0) : differential equation
C	  Ch(X,chi0) = atm_chap_num(X,chi0)  : real*4 numerical integral
C	  Ch(X,chi0) = atm8_chap_num(X,chi0) : real*8 numerical integral
C
C  ####################################################################

C  ====================================================================
C
C	These are the entries for the user to call.
C
C	chi0 can range from 0 to 180 in degrees.  For chi0 .gt. 90, the 
C	product X*(1-sin(chi0)) must not be too large, otherwise we 
C	will get an exponential overflow.
C
C	For chi0 .le. 90 degrees, X can range from 0 to thousands 
C	without overflow.
C
C  ====================================================================

	real*8 function atm_chapman( x, chi0 )
	real*8 atm8_chapman
	atm_chapman = atm8_chapman( dble(x), dble(chi0) )
	return
	end

c  ====================================================================

	real*8 function atm8_chapman( x, chi0 )
	implicit real*8(a-h,o-z)
	parameter (rad=57.2957795130823208768d0)

	if( (x .le. 0) .or. (chi0 .le. 0) .or. (chi0 .ge. 180) ) then
	  atm8_chapman = 1
	  return
	end if

	if( chi0 .gt. 90 ) then
	  chi = 180 - chi0
	else
	  chi = chi0
	end if

	if( x .lt. 36 ) then
	  atm8_chapman = atm8_chap_deq(x,chi)
	else
	  atm8_chapman = atm8_chap_asy(x,chi)
	end if

	if( chi0 .gt. 90 ) then
	  atm8_chapman = 2*exp(x*2*sin((90-chi)/(2*rad))**2)
     *		* atm8_chap_xk1(x*sin(chi/rad)) - atm8_chapman
	end if

	return
	end

c  ====================================================================
c
c	this chapman function routine calculates
c
c	  ch(x,chi0) = atm8_chap_asy(x,chi0)
c		     = sum{n=0,3} [c(n) * int{t=t0,infinity} 
c			[ exp(-x*(t**2-t0**2) * (t**2-t0**2)**n dy ] ]
c
c	with t0**2 = 1 - sin(chi0)
c
c  ====================================================================

	real*8 function atm8_chap_asy( x, chi0 )
	implicit real*8(a-h,o-z)
	parameter (rad=57.2957795130823208768d0)
	dimension c(0:3), xi(0:3), dn(0:3)
	common/atm8_chap_cm/fn(0:3)

	if( (x .le. 0) .or. (chi0 .le. 0) ) then
	  do i=0,3
	    fn(i) = 1
	  end do
	  go to 900
	end if

	sinchi = sin(chi0/rad)
	s1 = 1 + sinchi
	rx = sqrt(x)
	y0 = rx * sqrt( 2*sin( (90-chi0)/(2*rad) )**2 )

	c(0) = 1/sqrt(s1)
	fact = c(0)/s1
	c(1) = fact * (0.5d0+sinchi)
	fact = fact/s1
	c(2) = - fact * (0.125d0+0.5d0*sinchi)
	fact = fact/s1
	c(3) = fact * (0.0625d0+0.375d0*sinchi)

	call atm8_chap_gd3( y0, dn )
	fact = 2*rx
	do n=0,3
	  xi(n) = fact * dn(n)
	  fact = fact/x
	end do

	fn(0) = c(0) * xi(0)
	do i=1,3
	  fn(i) = fn(i-1) + c(i)*xi(i)
	end do

900	atm8_chap_asy = fn(3)
	return
	end

c  ====================================================================
c
c	this chapman function routine calculates
c
c	  ch(x,chi0) = atm8_chap_deq(x,chi0)
c		     = x * sin(chi0) * cos(chi0) * exp(x*sin(chi0))
c		       * int{y=x,infinity} [ exp(-y)*( 
c			 i1(x*sin(chi0))*k0(y*sin(chi0)) 
c			 + k1(x*sin(chi0))*i0(y*sin(chi0)) ) dy ]
c
c  ====================================================================

	real*8 function atm8_chap_deq( x, chi0 )
	implicit real*8(a-h,o-z)
	parameter (rad=57.2957795130823208768d0)
	common/atm8_chap_cm/xi1,xk1,yi0,yk0

	if( (x .le. 0) .or. (chi0 .le. 0) ) go to 800
	alpha = x * sin(chi0/rad)

c  --------------------------------------------------------------------
c
c	this code fragment calculates
c
c	  yi0 = exp(x*(1-sin(chi0))) * cos(chi0) * 
c		int{y=x,infinity} [ exp(-y) * i0(y*sin(chi0)) dy ]
c
c  --------------------------------------------------------------------

	yi0 = atm8_chap_yi0( x, chi0 )

c  --------------------------------------------------------------------
c
c	this code fragment calculates
c
c	  yk0 = exp(x*(1+sin(chi0))) x * sin(chi0) * cos(chi0) * 
c		int{y=x,infinity} [ exp(-y) * k0(y*sin(chi0)) dy ]
c
c  --------------------------------------------------------------------

	yk0 = atm8_chap_yk0( x, chi0 )

c  --------------------------------------------------------------------
c
c	this code fragment calculates
c
c	  xi1 = exp(-x*sin(chi0)) * i1(x*sin(chi0))
c
c  --------------------------------------------------------------------

	xi1 = atm8_chap_xi1( alpha )

c  --------------------------------------------------------------------
c
c	this code fragment calculates
c
c	  xk1 = x*sin(chi0) * exp(x*sin(chi0)) * k1(x*sin(chi0))
c
c  --------------------------------------------------------------------

	xk1 = atm8_chap_xk1( alpha )

c  --------------------------------------------------------------------
c
c	combine the terms
c
c  --------------------------------------------------------------------

	atm8_chap_deq = xi1*yk0 + xk1*yi0
	go to 900

800	atm8_chap_deq = 1
900	return
	end

c  ====================================================================
c
c	this chapman function routine calculates
c
c	  ch(x,chi0) = atm_chap_num(x,chi0) = numerical integral
c
c  ====================================================================

	real*4 function atm_chap_num(x,chi0)
	implicit real*8(a-h,o-z)
	real*4 x, chi0
	parameter (rad=57.2957795130823208768d0)
	parameter (n=65,nfact=8)
	dimension factor(0:nfact)
	data factor/14175.0d0, 23552.0d0, -3712.0d0, 41984.0d0,
     *	  -18160.0d0, 41984.0d0, -3712.0d0, 23552.0d0, 7912.0d0/

	if( (chi0 .le. 0) .or. (chi0 .gt. 90) .or. (x .le. 0) ) then
	  atm_chap_num = 1
	  return
	end if

	x8 = x
	chi0rad = chi0/rad
	sinchi = sin(chi0rad)

	alpha0 = asin( (x8/(x8+28)) * sinchi )
	delta = (chi0rad - alpha0)/(n-1)

	sum = 0

	do i=1,n
	  alpha = -(i-1)*delta + chi0rad

	  if( (i .eq. 1) .or. (x .le. 0) ) then
	    f = 1/(1+cos(alpha))
	  else if( alpha .le. 0 ) then
	    f = 0
	  else
	    f = exp(-x8*(sinchi/sin(alpha)-1) ) /(1+cos(alpha))
	  end if

	  if( (i.eq.1) .or. (i.eq.n) ) then
	    fact = factor(nfact)/2
	  else
	    fact = factor( mod(i-2,nfact)+1 )
	  end if

	  sum = sum + fact*f
	end do

	atm_chap_num = 1 + x8*sinchi*sum*delta/factor(0)
	return
	end

c  ====================================================================
c
c	this chapman function routine calculates
c
c	  ch(x,chi0) = atm8_chap_num(x,chi0) = numerical integral
c
c  ====================================================================

	real*8 function atm8_chap_num(x,chi0)
	implicit real*8(a-h,o-z)
	parameter (rad=57.2957795130823208768d0)
	parameter (n=601,nfact=8)
	dimension factor(0:nfact)
	data factor/14175.0d0, 23552.0d0, -3712.0d0, 41984.0d0,
     *	  -18160.0d0, 41984.0d0, -3712.0d0, 23552.0d0, 7912.0d0/

	if( (chi0 .le. 0) .or. (chi0 .gt. 90) .or. (x .le. 0) ) then
	  atm8_chap_num = 1
	  return
	end if

	chi0rad = chi0/rad
	sinchi = sin(chi0rad)

	alpha0 = asin( (x/(x+45)) * sinchi )
	delta = (chi0rad - alpha0)/(n-1)

	sum = 0

	do i=1,n
	  alpha = -(i-1)*delta + chi0rad

	  if( (i .eq. 1) .or. (x .le. 0) ) then
	    f = 1/(1+cos(alpha))
	  else if( alpha .le. 0 ) then
	    f = 0
	  else
	    f = exp(-x*(sinchi/sin(alpha)-1) ) /(1+cos(alpha))
	  end if

	  if( (i.eq.1) .or. (i.eq.n) ) then
	    fact = factor(nfact)/2
	  else
	    fact = factor( mod(i-2,nfact)+1 )
	  end if

	  sum = sum + fact*f
	end do

	atm8_chap_num = 1 + x*sinchi*sum*delta/factor(0)
	return
	end

c  ####################################################################
c
c	the following "bessel integral" routines return various 
c	combinations of integrals of bessel functions, powers, 
c	and exponentials, involving trigonometric functions of chi0.
c
c	for small values of z = x*sin(chi0) we expand
c
c	  i0(z) = sum{n=0,6} [ ai0(n) * z**(2*n) ]
c	  k0(z) = -log(z)*i0(z) + sum{n=0,6} [ ak0(n) * z**(2*n) ]
c
c	for large values of z we expand in reciprocal powers
c
c	  i0(z) = exp(z) * sum{n=0,8} [ bi0(n) * z**(-n-0.5) ]
c	  k0(z) = exp(-z) * sum{n=0,6} [ bk0(n) * z**(-n-0.5) ]
c
c	the expansion coefficients are calculated from those given 
c	by abramowitz and stegun [as64, pp378-9, section 9.8] and
c	press et al. [pft86, pp177-8, bessi0.for, bessk0.for].
c
c	for small values of x*sin(chi0) we break the integral
c	into two parts (with f(z) = i0(z) or k0(z)):
c
c	  int{y=x,infinity} [ exp(-y) * f(y*sin(chi0)) dy ]
c
c	    = int{y=x,x1} [ exp(-y) * f(y*sin(chi0)) dy ]
c	      + int{y=x1,infinity} [ exp(-y) * f(y*sin(chi0)) dy ]
c
c	where x1 = 3.75/sin(chi0) for i0 and 2/sin(chi0) for k0.
c
c	in the range y=x,x1 we integrate the term-by-term using
c
c	  int{z=a,b} [ exp(-z) * z**(2*n) dz ]
c	    = gamma(2*n+1,a) - gamma(2*n+1,b)
c
c	and a similar but more complicated formula for
c
c	  int{z=a,b} [ log(z) * exp(-z) * z**(2*n) dz ]
c
c	in the range y=x1,infinity we use
c
c	  int{z=b,infinity} [ exp(-z) * z**(-n-0.5) dz]
c	    = gamma(-n+0.5,b)
c
c  ####################################################################

c  ====================================================================
c
c	this bessel integral routine calculates
c
c	  yi0 = exp(x*(1-sin(chi0))) * cos(chi0) * 
c		int{y=x,infinity} [ exp(-y) * i0(y*sin(chi0)) dy ]
c
c  ====================================================================

	real*8 function atm8_chap_yi0( x, chi0 )
	implicit real*8(a-h,o-z)
	parameter (rad=57.2957795130823208768d0)
	dimension qbeta(0:8), gg(0:6)
	dimension ai0(0:6), bi0(0:8)

        data ai0/ 1.0000000d+00, 2.4999985d-01, 1.5625190d-02,
     *      4.3393973d-04, 6.8012343d-06, 6.5601736d-08,
     *      5.9239791d-10/
        data bi0/ 3.9894228d-01, 4.9822200d-02, 3.1685484d-02,
     *     -8.3090918d-02, 1.8119815d+00,-1.5259477d+01,
     *      7.3292025d+01,-1.7182223d+02, 1.5344533d+02/

	theta = (90-chi0)/(2*rad)
	sint = sin(theta)
	cost = cos(theta)
	sinchi = sin(chi0/rad)
	coschi = cos(chi0/rad)
	sc1m = 2*sint**2	! = (1-sinchi)

	alpha = x * sinchi

	if( alpha .le. 0 ) then
	  atm8_chap_yi0 = 1
	else if( alpha .lt. 3.75d0 ) then
	  x1 = 3.75d0/sinchi
	  call atm8_chap_gg06( x, x1, gg )
	  if( x .le. 1 ) then
	    rho = 1
	  else
	    rho = 1/x
	  end if
	  f = (sinchi/rho)**2
	  sum = ai0(6)*gg(6)
	  do i=5,0,-1
	    sum = sum*f + ai0(i)*gg(i)
c	    write(*,1900)i,sum,gg(i)
c1900	format(i5,1p5d14.6)
	  end do
	  call atm8_chap_gq85( x1*sc1m, qbeta )
	  sum2 = bi0(8) * qbeta(8)
	  do n=7,0,-1
	    sum2 = sum2/3.75d0 + bi0(n)*qbeta(n)
	  end do
	  atm8_chap_yi0 = exp(-alpha)*coschi*sum 
     *		+ exp((x-x1)*sc1m)*sum2*cost*sqrt(2/sinchi)
	else
	  call atm8_chap_gq85( x*sc1m, qbeta )
	  sum = bi0(8) * qbeta(8)
	  do n=7,0,-1
	    sum = sum/alpha + bi0(n)*qbeta(n)
	  end do
	  atm8_chap_yi0 = sum * cost * sqrt( 2 / sinchi )
	end if
	return
	end

c  ====================================================================
c
c	this bessel integral routine calculates
c
c	  yk0 = exp(x*(1+sin(chi0))) x * sin(chi0) * cos(chi0) * 
c		int{y=x,infinity} [ exp(-y) * k0(y*sin(chi0)) dy ]
c
c  ====================================================================

	real*8 function atm8_chap_yk0( x, chi0 )
	implicit real*8(a-h,o-z)
	parameter (rad=57.2957795130823208768d0)
	dimension ai0(0:6), ak0(0:6), bk0(0:6)
	dimension gf(0:6), gg(0:6), qgamma(0:8)
	
        data ai0/ 1.0000000d+00, 2.4999985d-01, 1.5625190d-02,
     *      4.3393973d-04, 6.8012343d-06, 6.5601736d-08,
     *      5.9239791d-10/
        data ak0/ 1.1593152d-01, 2.7898274d-01, 2.5249154d-02,
     *      8.4587629d-04, 1.4975897d-05, 1.5045213d-07,
     *      2.2172596d-09/
        data bk0/ 1.2533141d+00,-1.5664716d-01, 8.7582720d-02,
     *     -8.4995680d-02, 9.4059520d-02,-8.0492800d-02,
     *      3.4053120d-02/

	theta = (90-chi0)/(2*rad)
	sint = sin(theta)
	cost = cos(theta)
	sinchi = sin(chi0/rad)
	sc1 = 1+sinchi
	coschi = sin(2*theta)

	alpha = x * sinchi
	gamma = x * sc1

	if( alpha .le. 0 ) then
	  atm8_chap_yk0 = 0
	else if( alpha .lt. 2 ) then
	  x1 = 2/sinchi
	  call atm8_chap_gfg06( x, x1, gf, gg )
	  if( x .le. 1 ) then
	    rho = 1
	  else
	    rho = 1/x
	  end if
	  sl = log(sinchi)
	  f = (sinchi/rho)**2
	  sum = -ai0(6)*gf(6) + (-sl*ai0(6)+ak0(6))*gg(6)
	  do i=5,0,-1
	    sum = sum*f - ai0(i)*gf(i) + (-sl*ai0(i)+ak0(i))*gg(i)
c	    write(*,1900)i,sum,gf(i),gg(i)
c1900	format(i5,1p5d14.6)
	  end do
	  call atm8_chap_gq85( x1*sc1, qgamma )
	  sum2 = bk0(6)*qgamma(6)
	  do i=5,0,-1
	    sum2 = sum2*0.5d0 + bk0(i)*qgamma(i)
c	    write(*,1900)i,sum2,bk0(i),qgamma(i)
	  end do
	  sum = sum + exp(x-x1-2)*sum2/sqrt(sinchi*sc1)
	  atm8_chap_yk0 = sum * exp(alpha) * alpha * coschi
	else
	  call atm8_chap_gq85( gamma, qgamma )
	  sum = bk0(6) * qgamma(6)
	  do i=5,0,-1
	    sum = sum/alpha + bk0(i)*qgamma(i)
	  end do
	  atm8_chap_yk0 = sum * sint * sqrt( 2 * sinchi ) * x
	end if

	return
	end

c  ####################################################################
c
c	the following "pure math" routines return various combinations
c	of bessel functions, powers, and exponentials.
c
c  ####################################################################

c  ====================================================================
c
c	this bessel function math routine returns
c
c	  xi1 = exp(-|z|) * i1(z)
c
c	following press et al [pft86, page 178, bessi1.for] and 
c	abrahamson and stegun [as64, page 378, 9.8.3, 9.8.4].
c
c  ====================================================================

	real*8 function atm8_chap_xi1( z )
	implicit real*8(a-h,o-z)
        dimension ai1(0:6), bi1(0:8)

        data ai1/ 5.00000000d-01, 6.2499978d-02, 2.6041897d-03,
     *      5.4244512d-05, 6.7986797d-07, 5.4830314d-09,
     *      4.1909957d-11/
        data bi1/ 3.98942280d-01,-1.4955090d-01,-5.0908781d-02,
     *      8.6379434d-02,-2.0399403d+00, 1.6929962d+01,
     *     -8.0516146d+01, 1.8642422d+02,-1.6427082d+02/

	if( z .lt. 0 ) then
	  az = -z
	else if( z .eq. 0 ) then
	  atm8_chap_xi1 = 0
	  return
	else
	  az = z
	end if
	if( az .lt. 3.75d0 ) then
	  z2 = z*z
	  sum = ai1(6)
	  do i=5,0,-1
	    sum = sum*z2 + ai1(i)
	  end do
	  atm8_chap_xi1 = z*exp(-az) * sum
	else
	  sum = bi1(8)
	  do i=7,0,-1
	    sum = sum/az + bi1(i)
	  end do
	  atm8_chap_xi1 = sum*sqrt(az)/z
	end if
	return
	end

c  ====================================================================
c
c	this bessel function math routine returns
c
c	  xk1 = z * exp(+z) * k1(z)
c
c	following press et al [pft86, page 179, bessk1.for] and 
c	abrahamson and stegun [as64, page 379, 9.8.7, 9.8.8].
c
c  ====================================================================

	real*8 function atm8_chap_xk1( z )
	implicit real*8(a-h,o-z)
        dimension ak1(0:6), bk1(0:6)

        data ak1/ 1.00000000d+00, 3.8607860d-02,-4.2049112d-02,
     *     -2.8370152d-03,-7.4976641d-05,-1.0781641d-06,
     *     -1.1440430d-08/
        data bk1/ 1.25331414d+00, 4.6997238d-01,-1.4622480d-01,
     *      1.2034144d-01,-1.2485648d-01, 1.0419648d-01,
     *     -4.3676800d-02/

	if( z .le. 0 ) then
	  atm8_chap_xk1 = 1
	else if( z .lt. 2 ) then
	  xz = exp(z)
	  z2 = z*z
	  sum = ak1(6)
	  do i=5,0,-1
	    sum = sum*z2 + ak1(i)
	  end do
	  atm8_chap_xk1 = xz * ( sum 
     *		+ z*log(z/2)*atm8_chap_xi1(z)*xz )
	else
	  sum = bk1(6)
	  do i=5,0,-1
	    sum = sum/z + bk1(i)
	  end do
	  atm8_chap_xk1 = sum*sqrt(z)
	end if

	return
	end

c  ####################################################################
c
c	the following "pure math" routines return various combinations
c	of the error function, powers, and exponentials.
c
c  ####################################################################

c  ====================================================================
c
c	this error function math routine returns
c
c	  xerfc(x) = exp(x**2)*erfc(x)
c
c	following press et al. [pft86, p164, erfcc.for]
c
c  ====================================================================

	real*8 function atm8_chap_xerfc(x)
	implicit real*8(a-h,o-z)
        t=1.0d0/(1.0d0+0.5d0*x)
	atm8_chap_xerfc =
     *	  t*exp( -1.26551223d0 +t*(1.00002368d0 +t*( .37409196d0
     *       +t*(  .09678418d0 +t*(-.18628806d0 +t*( .27886807d0
     *	     +t*(-1.13520398d0 +t*(1.48851587d0 +t*(-.82215223d0
     *	     +t*   .17087277d0) ))))))))
        return
        end

c  ####################################################################
c
c	the following "pure math" routines return various combinations
c	of exponential integrals, powers, and exponentials.
c
c  ####################################################################

c  ====================================================================
c
c	this exponential math routine evaluates
c
c	  zxe1(x) = x*exp(x) int{y=1,infinity} [ exp(-x*y)/y dy ]
c
c	following abramowitz and stegun [as64, p229;231, equations
c	5.1.11 and 5.1.56]
c
c  ====================================================================

	real*8 function atm8_chap_zxe1(x)
	implicit real*8(a-h,o-z)
	parameter (gamma = 0.5772156649015328606d0)
	dimension ae1(0:4), be1(0:4), cein(1:10)

	data ae1/1.0d0, 8.5733287401d0, 18.0590169730d0,
     *	    8.6347608925d0, 0.2677737343d0 /
	data be1/1.0d0, 9.5733223454d0, 25.6329561486d0,
     *	    21.0996530827d0, 3.9584969228d0/
        data cein/ 1.00000000d+00,-2.50000000d-01, 5.55555556d-02,
     *    -1.0416666667d-02, 1.6666666667d-03,-2.3148148148d-04,
     *     2.8344671202d-05,-3.1001984127d-06, 3.0619243582d-07,
     *    -2.7557319224d-08/

	if( x .le. 0 ) then
	  atm8_chap_zxe1 = 0
	else if( x .le. 1 ) then
	  sum = cein(10)
	  do i=9,1,-1
	    sum = sum*x + cein(i)
	  end do
	  atm8_chap_zxe1 = x*exp(x)*( x * sum - log(x) - gamma )
	else
	  top = ae1(4)
	  bot = be1(4)
	  do i=3,0,-1
	    top = top/x + ae1(i)
	    bot = bot/x + be1(i)
	  end do
	  atm8_chap_zxe1 = top/bot
	end if
	return
	end

c  ####################################################################
c
c	the following "pure math" routines return various combinations
c	of incomplete gamma functions, powers, and exponentials.
c
c  ####################################################################

c  ====================================================================
c
c	this gamma function math routine calculates
c
c	dn(n) = int{t=z,infinity}
c		[ exp( -(t**2-z**2) ) * (t**2-z**2)**n dt ]
c
c  ====================================================================

	subroutine atm8_chap_gd3( z, dn )
	implicit real*8(a-h,o-z)
	parameter (rpi=1.7724538509055160273d0)
	dimension dn(0:3), xg(0:3)

	if( z .le. 0 ) then
	  dn(0) = rpi/2
	  do i=1,3
	    dn(i) = (i-0.5d0)*dn(i-1)
	  end do
	  return
	end if

	z2 = z*z
	if( z .ge. 7 ) r = 1/z2

	if( z .lt. 14 ) then
	  z4 = z2*z2
	  xg(0) = rpi * atm8_chap_xerfc(z)
	  xg(1) = 0.5d0*xg(0) + z
	  xg(2) = 1.5d0*xg(1) + z*z2
	  dn(0) = 0.5d0*xg(0)
	  dn(1) = 0.5d0*(xg(1)-z2*xg(0))
	  dn(2) = 0.5d0*(xg(2)-2*z2*xg(1)+z4*xg(0))
	else
	  dn(0) = ( 1 + r*(-0.5d0 +r*(0.75d0 +r*(-1.875d0
     *		+r*6.5625d0) ) ) )/(2*z)
	  dn(1) = ( 1 + r*(-1.0d0 +r*(2.25d0 +r*(-7.5d0
     *		+r*32.8125d0) ) ) )/(2*z)
	  dn(2) = ( 2 + r*(-3.0d0 +r*(9.00d0 +r*(-37.5d0
     *		+r*196.875d0) ) ) )/(2*z)
	end if

	if( z .lt. 7 ) then
	  z6 = z4*z2
	  xg(3) = 2.5d0*xg(2) + z*z4
	  dn(3) = 0.5d0*(xg(3)-3*z2*xg(2)+3*z4*xg(1)-z6*xg(0))
	else
	  dn(3) = ( 6 + r*(-12.0d0 +r*(45.0d0 +r*(-225.0d0
     *		+r*1378.125d0) ) ) )/(2*z)
	end if

	return
	end

c  ====================================================================
c
c	this gamma function math routine calculates
c
c	  gf06(n) = g(n,x) * int{y=x,z} [log(y) * exp(-y) * y**(2*n) dy]
c
c	and
c
c	  gg06(n) = g(n,x) * int{y=x,z} [ exp(-y) * y**(2*n) dy ]
c	          = g(n,x) * ( gamma(2*n+1,x) - gamma(2*n+1,z) )
c
c	for n=0,6, with g(n,x) = exp(x) * max(1,x)**(-2*n)
c
c  ====================================================================

	subroutine atm8_chap_gfg06( x, z, gf06, gg06 )
	implicit real*8 (a-h,o-z)
	parameter (gamma = 0.5772156649015328606d0)
	dimension gf06(0:6), gg06(0:6)
	dimension gh13x(13), gh13z(13), rgn(13), delta(13)
	call atm8_chap_gh13( x, x, gh13x )
	call atm8_chap_gh13( x, z, gh13z )
	if( x .le. 1 ) then
	  rho = 1
	else
	  rho = 1/x
	end if

	delta(1) = 0
	delta(2) = ( gh13x(1) - gh13z(1) ) * rho
	rgn(1) = 1
	rgn(2) = rho
	do n=2,12
	  delta(n+1) = rho*( n*delta(n) + gh13x(n) - gh13z(n) )
	  rgn(n+1) = (n*rho)*rgn(n)

	end do

	if( x .gt. 0 ) then
	  xe1_x = atm8_chap_zxe1(x)/x
	  xlog = log(x)
	end if
	if( z .gt. 0 ) then
	  xe1_z = exp(x-z)*atm8_chap_zxe1(z)/z
	  zlog = log(z)
	end if

	do k=0,6
	  n = 2*k+1
	  if( x .le. 0 ) then
	    gf06(k) = -gamma*rgn(n) + delta(n)
	  else
	    gf06(k) = xlog*gh13x(n) + rgn(n)*xe1_x + delta(n)
	  end if
	  if( z .le. 0 ) then
	    gf06(k) = gf06(k) + gamma*rgn(n)
	  else
	    gf06(k) = gf06(k) - (zlog*gh13z(n) + rgn(n)*xe1_z)
	  end if
	  gg06(k) = gh13x(n) - gh13z(n)
	end do

	return
	end

c  ====================================================================
c
c	this gamma function math routine calculates
c
c	  gg06(n) = g(n,x) * int{y=x,z} [ exp(-y) * y**(2*n) dy ]
c	          = g(n,x) * ( gamma(2*n+1,x) - gamma(2*n+1,z) )
c
c	for n=0,6, with g(n,x) = exp(x) * max(1,x)**(-2*n)
c
c  ====================================================================

	subroutine atm8_chap_gg06( x, z, gg06 )
	implicit real*8 (a-h,o-z)
	dimension gg06(0:6), gh13x(13), gh13z(13)
	call atm8_chap_gh13( x, x, gh13x )
	call atm8_chap_gh13( x, z, gh13z )
	do n=0,6
	  gg06(n) = gh13x(2*n+1) - gh13z(2*n+1)
	end do
	return
	end

c  ====================================================================
c
c	this gamma function math routine calculates
c
c	  gh13(n) = f(n,x) * int{y=z,infinity} [exp(-y) * y**(n-1) dy]
c	          = f(n,x) * gamma(n,z)
c
c	for n=1,13, with f(n,x) = exp(x) * max(1,x)**(-n+1)
c
c  ====================================================================

	subroutine atm8_chap_gh13( x, z, gh13 )
	implicit real*8 (a-h,o-z)
	dimension gh13(13), tab(12)

	if( z .le. 0 ) then
	  gh13(1) = 1
	  do n=1,12
	    gh13(n+1) = n*gh13(n)
	  end do
	  return
	end if

	if( x .le. 1 ) then
	  rho = 1
	else
	  rho = 1/x
	end if
	rhoz = rho * z
	exz = exp(x-z)
	tab(12) = exp( (x-z) + 12*log(rhoz) )
	do n=11,1,-1
	  tab(n) = tab(n+1)/rhoz
	end do
	gh13(1) = exz
	do n=1,12
	  gh13(n+1) = rho*n*gh13(n) + tab(n)
	end do
	return
	end

c  ====================================================================
c
c	this gamma function math subroutine calculates
c
c	  qn(x) = x**n * exp(x) * gamma(-n+0.5,x), n=0,8
c	    = x**n * exp(x) * int{y=x,infinity} [exp(-y)*y**(-n-0.5)dy]
c
c	for x .lt. 2 we first calculate
c
c	  q0(x) = sqrt(pi)*exp(x)*erfc(sqrt(x)) = exp(x)*gamma(0.5,x)
c
c	and use upward recursion.  else, we first calculate
c
c	  q8(x) = x**8 * exp(x) * gamma(-7.5,x)
c
c	following press et al. [pft86, pp162-63, gcf.for] and then
c	recur downward.  also see abramowitz and stegun [as64, 6.5].
c
c  ====================================================================

	subroutine atm8_chap_gq85( x, qn )
	implicit real*8(a-h,o-z)
	parameter (rpi=1.7724538509055160273d0)
	parameter (itmax=100,eps=3.0d-9)
	dimension qn(0:8)

	if( x .le. 0 ) then
	  qn(0) = rpi
	  do i=1,8
	    qn(i) = 0
	  end do
	  return
	end if

	rx = sqrt(x)

	if( x .lt. 2 ) then
	  qn(0) = rpi * atm8_chap_xerfc( rx )
	  do n=1,8
	    qn(n) = ( -rx*qn(n-1) + 1 ) * rx / ( n - 0.5d0 )
	  end do
	else
          gold=0.0d0
	  a0=1.0d0
	  a1=x
	  b0=0.0d0
	  b1=1.0d0
	  fac=1.0d0
	  do 11 n=1,itmax
	    an= (n)
	    ana=an + 7.5d0
	    a0=(a1+a0*ana)*fac
	    b0=(b1+b0*ana)*fac
	    anf=an*fac
	    a1=x*a0+anf*a1
	    b1=x*b0+anf*b1
	    fac=1./a1
            g=b1*fac
	    test = g*eps
	    del = g - gold
	    if( test .lt. 0 ) test = - test
	    if( (del .ge. -test) .and. (del .le. test) ) go to 12
	    gold=g
11        continue
12	  qn(8) = g * rx
	  do n=8,1,-1
	    qn(n-1) = ( (-n+0.5d0)*qn(n)/rx + 1 ) / rx
	  end do
	end if

	return
	end


*******************************************
*******************************************
*          splinenr
*    (from numerical recipes) 
*******************************************
*******************************************

      subroutine splinenr(x,y,n,yp1,ypn,y2)
      parameter (nmax=200)
      dimension x(n),y(n),y2(n),u(nmax)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      end

*******************************************
*******************************************
*          splintnr
*    (from numerical recipes) 
*******************************************
*******************************************

      subroutine splintnr(xa,ya,y2a,n,x,y)
      dimension xa(n),ya(n),y2a(n)
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input. from splintnr'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     *      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      end


