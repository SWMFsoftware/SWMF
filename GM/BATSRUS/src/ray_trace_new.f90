!^CFG COPYRIGHT UM
!^CFG FILE RAYTRACE

! The main subroutines in this file are
! subroutine ray_trace_accurate     - trace all rays starting from 3D MHD grid
! subroutine integrate_ray_accurate - integrate rays starting from 2D IM grid
! subroutine write_plot_line        - extract lines into plot file(s)

subroutine ray_trace_accurate

  ! Trace field lines from cell centers to the outer or inner boundaries

  use ModProcMH
  use ModRaytrace
  use CON_ray_trace, ONLY: ray_init
  use ModMain
  use ModPhysics,    ONLY: rBody
  use ModAdvance,    ONLY: State_VGB, Bx_, By_, Bz_, &
       B0xCell_BLK,B0yCell_BLK,B0zCell_BLK
  use ModGeometry,   ONLY: x_BLK,y_BLK,z_BLK,r_BLK,true_cell,XyzMax_D,XyzMin_D

  use ModMpi
  implicit none

  ! Arguments
  ! remember last call and the last grid number
  integer :: n_last=-1, iLastGrid=-1, iLastDecomposition=-1

  ! Indices corresponding to the starting point and directon of the ray
  integer :: i, j, k, iBlock, iRay

  ! Testing and timing
  logical :: oktest, oktest_me, oktime, oktime_me

  !----------------------------------------------------------------------------

  call set_oktest('ray_trace',oktest,oktest_me)

  if(oktest_me)then
     write(*,*)'GM ray_trace: n_last,n_step         =',n_last,n_step
     write(*,*)'GM ray_trace: iLastGrid,iNewGrid    =',iLastGrid,iNewGrid
     write(*,*)'GM ray_trace: iLastDecomp,iNewDecomp=',&
          iLastDecomposition,iNewDecomposition
  end if

  if(n_last == n_step .and. iLastGrid == iNewGrid .and. &
       iLastDecomposition == iNewDecomposition) &
       RETURN

  ! Initialize constants
  R_raytrace   = rBody
  R2_raytrace  = R_raytrace**2
  RayLengthMax = 4*sum(XyzMax_D - XyzMin_D)

  NameTask        = 'trace'
  NameVectorField = 'B'

  ! (Re)initialize CON_ray_trace
  call ray_init(iComm)

  ! Remember this call
  n_last=n_step; iLastGrid = iNewGrid; iLastDecomposition = iNewDecomposition

  call set_oktest('time_ray_trace',oktime,oktime_me)
  call timing_start('ray_trace')
  if(oktime)call timing_reset('ray_pass',2)

  ! Copy magnetic field into Bxyz_DGB
  Bxyz_DGB(:,:,:,:,1:nBlock) = State_VGB(Bx_:Bz_,:,:,:,1:nBlock)

  ! Fill in all ghost cells (faces+edges+corners) without monotone restrict
  call message_pass_cells8(.false.,.false.,.false.,3,Bxyz_DGB)

  ! Add B0
  Bxyz_DGB(1,:,:,:,1:nBlock) = Bxyz_DGB(1,:,:,:,1:nBlock) &
       + B0xCell_BLK(:,:,:,1:nBlock)
  Bxyz_DGB(2,:,:,:,1:nBlock) = Bxyz_DGB(2,:,:,:,1:nBlock) &
       + B0yCell_BLK(:,:,:,1:nBlock)
  Bxyz_DGB(3,:,:,:,1:nBlock) = Bxyz_DGB(3,:,:,:,1:nBlock) &
       + B0zCell_BLK(:,:,:,1:nBlock)

  ! Initial values
  ray=NORAY

  if(oktest_me)write(*,*)'rayface normalized B'
  if(oktime.and.iProc==0)then
     write(*,'(a)',ADVANCE='NO') 'setup and normalization:'
     call timing_show('ray_trace',1)
  end if

  ! This loop order seems to give optimal speed
  do k = 1, nK; do j = 1, nJ; do i = 1, nI
     do iBlock = 1, nBlock
        if(unusedBLK(iBlock))CYCLE

        oktest_ray = oktest .and. &
             x_BLK(i,j,k,iBlock)==xTest .and. &
             y_BLK(i,j,k,iBlock)==yTest .and. &
             z_BLK(i,j,k,iBlock)==zTest

        do iRay=1,2

           ! Short cut for inner and false cells
           if(R_BLK(i,j,k,iBlock) < rIonosphere .or. &
                .not.true_cell(i,j,k,iBlock))then
              ray(:,:,i,j,k,iBlock)=BODYRAY
              if(oktest_ray)write(*,*)'Shortcut BODYRAY iProc,iRay=',&
                   iProc,iRay
              CYCLE
           end if

           if(oktest_ray)write(*,*)'calling follow_ray iProc,iRay=',&
                iProc,iRay

           ! Follow ray in direction iRay
           call follow_ray(iRay, (/i,j,k,iBlock/), &
                (/ x_BLK(i,j,k,iBlock), y_BLK(i,j,k,iBlock), &
                z_BLK(i,j,k,iBlock) /) )

        end do             ! iRay
     end do          ! iBlock
  end do; end do; end do  ! i, j, k

  ! Do remaining rays passed from other PE-s
  call finish_ray

  ! Convert x, y, z to latitude and longitude, and status
  do iBlock=1,nBlock
     if(unusedBLK(iBlock)) CYCLE
     do k=1,nK; do j=1,nJ; do i=1,nI

        call xyz_to_latlonstatus(ray(:,:,i,j,k,iBlock))

        if(ray(3,1,i,j,k,iBlock)==-2.) write(*,*) &
             'Loop ray found at iProc,iBlock,i,j,k,ray=',&
             iProc,iBlock,i,j,k,ray(:,:,i,j,k,iBlock)

        if(ray(3,1,i,j,k,iBlock)==-3.) write(*,*) &
             'Strange ray found at iProc,iBlock,i,j,k,ray=',&
             iProc,iBlock,i,j,k,ray(:,:,i,j,k,iBlock)
     end do; end do; end do
  end do

  if(oktest_me)write(*,*)'ray lat, lon, status=',&
       ray(:,:,iTest,jTest,kTest,BlkTest)

  call timing_stop('ray_trace')

  if(oktime.and.iProc==0)then
     write(*,'(a)',ADVANCE='NO') 'Total ray tracing time:'
     call timing_show('ray_trace',1)
  end if
  call barrier_mpi
  if(oktest_me)write(*,*)'ray_trace completed.'

end subroutine ray_trace_accurate

!===========================================================================
subroutine finish_ray

  ! This subroutine is a simple interface for the last call to follow_ray
  call follow_ray(0, (/0, 0, 0, 0/), (/ 0., 0., 0. /))

end subroutine finish_ray

!===========================================================================

subroutine follow_ray(iRayIn,i_D,XyzIn_D)

  !DESCRIPTION:
  ! Follow ray in direction iRayIn (1 is parallel with the field, 
  !                                 2 is anti-parallel, 
  !                                 0 means that no ray is passed 
  ! Always follow rays received from other PE-s.
  !
  ! The passed ray is identified by the four dimensional index array i\_D.
  ! The meaning of i\_D d depends on the context: 
  !  3 cell + 1 block index for 3D ray tracing
  !  1 latitude + 1 longitude index for ray integration
  !  1 linear index for ray extraction.
  !
  ! The rays are followed until the ray hits the outer or inner 
  ! boundary of the computational domain. The results are saved into
  ! arrays defined in ModRayTrace or into files based on 
  ! ModRaytrace::NameTask:
  !
  ! If NameTask='trace', follow the ray, and save the final position into
  !    ModRayTrace::ray(:,iRayIn,i_D(1),i_D(2),i_D(3),i_D(4)) on the 
  !    processor that started the ray trace.
  !
  ! If NameTask='integrate', do integration along the ray and
  !    save the integrals into ModRayTrace::RayIntegral_VII(i_D(1),i_D(2))
  !
  ! If NameTask='extract' or 'extractfile' or 'extractfiles', 
  !    extract data along the ray and save it to either 
  !    an array ???, one file, or multiple files (one for each line).
  !    In this case the rays are indexed with i_D(1).
  !
  !EOP

  use ModRayTrace
  use CON_ray_trace

  use ModMain,     ONLY: iTest, jTest, kTest, BlkTest, ProcTest
  use ModGeometry, ONLY: x_BLK,y_BLK,z_BLK,XyzStart_BLK,Dx_BLK
  use ModProcMH
  use ModKind

  use ModMpi
  implicit none

  !INPUT ARGUMENTS:
  integer, intent(in) :: iRayIn     ! ray direction, 0 if no ray is passed
  integer, intent(in) :: i_D(4)     ! general index array for starting position
  real,    intent(in) :: XyzIn_D(3) ! coordinates of starting position

  !LOCAL VARIABLES:
  ! Cell, block and PE indexes for initial position and ray direction
  integer :: iStart, jStart, kStart, iBlockStart, iProcStart, iRay
  integer :: iStart_D(4)

  ! Current position of the ray
  integer :: iBlockRay
  real    :: XyzRay_D(3)

  ! Current length of ray
  real    :: RayLength

  ! Is the ray trace done
  logical :: DoneRay

  ! Shall we get ray from other PE-s
  logical :: DoGet

  ! Did we get rays from other PE-s
  logical :: IsFound

  ! Is the ray parallel with the vector field
  logical :: IsParallel

  integer, parameter :: MaxCount = 1000
  integer :: iFace, iCount, jProc, jBlock, i, j, k

  logical :: DoneAll
  integer :: iCountRay = 0, DiCountRay = 2000

  real(Real8_) :: CpuTimeNow

  character(len=*), parameter :: NameSub='follow_ray'

  logical :: DoTest, DoTestMe

  !-----------------------------------------------------------------------

  call set_oktest(NameSub, DoTest, DoTestMe)

  if(iRayIn /= 0)then

     ! Store starting indexes and ray direction
     iStart = i_D(1); jStart = i_D(2); kStart = i_D(3); 
     iBlockStart = i_D(4); iProcStart = iProc
     iRay   = iRayIn

     iStart_D = i_D
     if(DoTest)call set_oktest_ray

     ! Current position and length
     iBlockRay = i_D(4)
     XyzRay_D  = XyzIn_D
     RayLength = 0.0

     if(oktest_ray)write(*,'(a,6i4,3es12.4)')&
          'Local ray at iProc,i_D,iRay,XyzIn_D=',iProc,i_D,iRay,XyzIn_D

  end if

  ! If iRayIn==0 there are no more local rays to follow so get from other PEs
  DoGet = iRayIn==0
  IsFound = .true.

  RAYS: do

     if(DoGet)then
        GETRAY: do
           call ray_get(IsFound,iProcStart,iStart_D,XyzRay_D,RayLength,&
                IsParallel,DoneRay)

           if(IsFound)then
              if(DoTest)call set_oktest_ray

              if(IsParallel)then
                 iRay=1
              else
                 iRay=2
              end if
              if(oktest_ray)write(*,*)'Recv ray iProc,iRay,Done,XyzRay_D=',&
                   iProc,iRay,DoneRay,XyzRay_D

              if(DoneRay)then
                 if(NameTask /= 'trace')then
                    write(*,*)NameSub,' WARNING ',&
                         'received DoneRay=T for task ',trim(NameTask),' !!!'
                    CYCLE GETRAY
                 end if

                 ! Store the result into the ModRayTrace::ray
                 iStart      = iStart_D(1)
                 jStart      = iStart_D(2)
                 kStart      = iStart_D(3)
                 iBlockStart = iStart_D(4)

                 ray(:,iRay,iStart,jStart,kStart,iBlockStart)=XyzRay_D

                 if(oktest_ray)write(*,*)&
                      'Storing recv ray iProc,iRay,i,j,k,iBlock,ray=',&
                      iProc,iRay,iStart,jStart,kStart,iBlockStart,XyzRay_D

                 ! Get another ray from the others
                 CYCLE GETRAY
              else
                 ! Find block for the received ray
                 call xyz_to_peblk(XyzRay_D(1),XyzRay_D(2),XyzRay_D(3),&
                      jProc,iBlockRay,.false.,i,j,k)
                 if(jProc /= iProc)call stop_mpi(&
                      'GM_ERROR in ray_trace: Recvd ray is not in this PE')

                 if(oktest_ray)write(*,*)'Block for recv ray iProc,iBlock=',&
                      iProc,iBlockRay
              end if
           end if
           EXIT GETRAY
        end do GETRAY
     end if ! DoGet

     if(IsFound)then
        call follow_this_ray
        DoGet = .true.
     else
        if(iRayIn>0)then
           ! Stop working on received rays if there are no more
           ! but there are still local rays
           EXIT RAYS
        else
           ! Try to get more rays from others and check if everyone is done
           call ray_exchange(.true.,DoneAll)
           if(DoneAll)then
              EXIT RAYS
           else
              CYCLE RAYS
           end if
        end if
     end if

     iCountRay = iCountRay + 1

     if(iRayIn>0)then
        ! If there are still local rays, exchange only occasionally
        CpuTimeNow = MPI_WTIME()
        if(CpuTimeNow - CpuTimeStartRay > DtExchangeRay)then
           ! This PE is not done yet, so pass .false.
           call ray_exchange(.false., DoneAll)
           CpuTimeStartRay = CpuTimeNow
        end if
     end if

  end do RAYS

contains

  !=========================================================================
  subroutine follow_this_ray

    ! Initialize integrals for this segment
    if(NameTask=='integrate')RayIntegral_V = 0.0

    ! Follow the ray through the local blocks
    BLOCK: do iCount = 1, MaxCount

       if(iCount < MaxCount)then
          call follow_ray_block(iStart_D,iRay,iBlockRay,XyzRay_D,&
               RayLength,iFace)
       else
          write(*,*)NameSub,' WARNING ray passed through more than MaxCount=',&
               MaxCount,' blocks:'
          write(*,*)NameSub,'    iStart_D    =',iStart_D
          write(*,*)NameSub,'    XyzRay_D    =',XyzRay_D
          write(*,*)NameSub,'    XyzStart_BLK=',XyzStart_BLK(:,iBlockRay)
          iFace = ray_loop_
       end if

       select case(iFace)
       case(ray_block_)

          ! Find the new PE and block for the current position

          call xyz_to_peblk(XyzRay_D(1),XyzRay_D(2),XyzRay_D(3),&
               jProc,jBlock,.false.,i,j,k)

          if(jProc /= iProc)then
             ! Send ray to the next processor and return from here
             if(oktest_ray)write(*,*)'Sending ray iProc,jProc,iRay,Xyz=',&
                  iProc,jProc,iRay,XyzRay_D

             ! Add partial results to the integrals. 
             ! Pass .false., because this is not the final position
             if(NameTask=='integrate')call store_integral(.false.)

             call ray_put(iProcStart,iStart_D,jProc,XyzRay_D,RayLength,&
                  iRay==1,.false.)
             RETURN
          elseif(jBlock /= iBlockRay)then
             ! Continue the same ray in the next block
             iBlockRay = jBlock
             if(oktest_ray)write(*,'(a,3i4,3es12.4)')&
                  'Continuing ray iProc,jBlock,iRay,Xyz=',&
                  iProc,jBlock,iRay,XyzRay_D
             CYCLE BLOCK
          else
             write(*,*)'ERROR for follow_this_ray, iProc=',iProc
             write(*,*)'ERROR iBlockRay=jBlock=',jBlock
             write(*,*)'ERROR for iStart_D    =',iStart_D
             write(*,*)'ERROR for XyzRay_D    =',XyzRay_D
             write(*,*)'XyzStart_BLK, Dx_BLK  =',XyzStart_BLK(:,jBlock),&
                  Dx_BLK(jBlock)
             call stop_mpi(&
                  'GM_ERROR in follow_ray: continues in same BLOCK')
          end if
       case(ray_open_)
          ! The ray reached the outer boundary (or expected to do so)
          XyzRay_D = OPENRAY
          if(oktest_ray)write(*,*)&
               'follow_ray finished with OPENRAY, iProc,iRay=',iProc,iRay

       case(ray_loop_)
          ! The ray did not hit the wall of the block
          XyzRay_D = LOOPRAY
          if(oktest_ray)write(*,*)&
               'follow_ray finished with LOOPRAY, iProc,iRay=',iProc,iRay

       case(ray_body_)
          ! The ray hit a body
          XyzRay_D = BODYRAY
          if(oktest_ray)write(*,*)&
               'follow_ray finished with BODYRAY, iProc,iRay=',iProc,iRay

       case(ray_iono_)
          ! The ray hit the ionosphere 
          if(oktest_ray)write(*,'(a,2i4,3es12.4)')&
               'follow_this_ray finished on the ionosphere '// &
               'at iProc,iRay,Xyz=',iProc,iRay,XyzRay_D

       case default
          write(*,*)'Impossible value for iface=',iFace,&
               ' at XyzRay_D,iBlockRay=',XyzRay_D,iBlockRay
          call stop_mpi('GM_ERROR in follow_ray: impossible iFace value')
       end select

       ! Store integrals and the final position
       if(NameTask == 'integrate')call store_integral(.true.)

       ! Nothing more to do if not tracing
       if(NameTask /= 'trace') EXIT BLOCK

       ! For tracing either store results or send them back to starting PE
       if(iProcStart == iProc)then

          ! Store the result into the ModRayTrace::ray
          iStart      = iStart_D(1)
          jStart      = iStart_D(2)
          kStart      = iStart_D(3)
          iBlockStart = iStart_D(4)

          ray(:,iRay,iStart,jStart,kStart,iBlockStart)=XyzRay_D

          if(oktest_ray)write(*,*) &
               'Storing into iProc,iBlock,i,j,k,iRay,Xyz=',&
               iProc,iBlockStart,iStart,jStart,kStart,iRay,XyzRay_D

       else
          ! Send back result to iProcStart. 
          call ray_put(iProcStart,iStart_D,iProc,XyzRay_D,RayLength,&
               iRay==1,.true.)

          if(oktest_ray)write(*,*) &
               'Send result iProc,iProcStart,iRay,Xyz=',&
               iProc,iProcStart,iRay,XyzRay_D

       end if
       EXIT BLOCK

    end do BLOCK

  end subroutine follow_this_ray

  !===========================================================================
  subroutine store_integral(DoneRay)

    ! Store integrals of this ray into the 

    logical, intent(in) :: DoneRay

    integer :: iLat, iLon

    iLat = iStart_D(1)
    iLon = iStart_D(2)

    RayIntegral_VII(1:pInvB_,iLat,iLon) = &
         RayIntegral_VII(1:pInvB_,iLat,iLon) + RayIntegral_V

    if(DoneRay)then
       RayIntegral_VII(xEnd_:zEnd_,iLat,iLon) = XyzRay_D
       RayIntegral_VII(Length_,iLat,iLon)     = RayLength
    end if

  end subroutine store_integral

  subroutine set_oktest_ray

    select case(NameTask)
    case('extract','extractfile','extractfiles')
       oktest_ray = DoTest .and. iStart_D(1) == iTest
    case('integrate')
       oktest_ray = DoTest .and. all(iStart_D(1:2) == (/iTest,jTest/))
    case('trace')
       oktest_ray = DoTest .and. iProcStart == ProcTest .and. &
            all(iStart_D == (/iTest,jTest,kTest,BlkTest/))
    case default
       call stop_mpi(NameSub//' ERROR: invalid NameTask='//NameTask)
    end select

  end subroutine set_oktest_ray

end subroutine follow_ray

!==========================================================================
subroutine follow_ray_block(iStart_D,iRay,iBlock,Xyz_D,Length,iFace)

  !DESCRIPTION:
  ! Follow ray identified by index array iStart_D, 
  ! starting at initial position Xyz_D inside block iBlock,
  ! in direction iRay until we hit the wall of the block or the ionosphere. 
  ! Return Xyz_D with the final position. 
  ! Integrate and/or extract values if required.
  ! Also return Length increased by the length of the ray in this block.
  !
  ! Return iFace = 1..6 if the ray hit the east,west,south,north,bot,top walls
  ! Return ray_iono_    if the ray hit the ionosphere
  ! Return ray_loop_    if the ray did not hit anything 
  ! Return ray_body_    if the ray goes into or is inside a body
  ! Return ray_open_    if the ray goes outside the computational box
  !EOP

  use ModRayTrace
  use ModProcMH
  use ModNumConst, ONLY: cTiny, cOne
  use ModMain, ONLY: TypeCoordSystem, nI, nJ, nK, nBLK
  use ModGeometry, ONLY: XyzStart_BLK, XyzMax_D, XyzMin_D, &
       x_BLK,y_BLK,z_BLK, Dx_BLK, Dy_BLK, Dz_BLK, rMin_BLK
  use ModAdvance, ONLY : State_VGB, Bx_, By_, Bz_

  implicit none

  ! Arguments

  integer, intent(in) :: iStart_D(4)
  integer, intent(in) :: iRay
  integer, intent(in) :: iBlock
  real, intent(inout) :: Xyz_D(3)
  real, intent(inout) :: Length
  integer, intent(out):: iFace

  ! Local variables

  ! Block size
  real :: Dxyz_D(3)

  ! Initial and mid point normalized coordinates and direction of B field
  real, dimension(3) :: x_ini, x_mid, b_ini, b_mid

  ! True interpolated magnetic field, and true location
  real, dimension(3) :: b_D, xx_ini

  ! Radial distance from origin
  real :: r, r_ini

  ! dx is the difference between 1st and 2nd order RK to estimate accuracy
  ! dx_opt is the required accuracy, dx_rel=dx/dx_opt
  real :: dx_rel, dx_opt

  ! Ray length, maximum length in block, step size, next step size 
  real :: l, lmax, dl, dLength, dl_next

  ! Fraction of the last step inside the ionosphere
  real :: Fraction

  ! Step size limits
  real, parameter :: dl_max=1.0, dl_min=0.05, dl_tiny=1.e-6

  ! counter for ray integration
  integer :: nsegment 

  ! True if Rmin_BLK < R_raytrace
  logical :: DoCheckInnerBc

  ! True if the block already containes open rays
  logical :: DoCheckOpen

  ! Counter for entering follow_ray_iono
  integer :: n_iono

  ! Control volume limits in local coordinates
  real, dimension(3) :: xmin, xmax

  ! Current position of ray in normalized and physical coordinates
  real, dimension(3) :: x, xx

  ! Radial distance and square of it: r2=sum(xx**2)
  real :: r2

  ! Cell indices corresponding to current or final x position
  integer :: i1,j1,k1,i2,j2,k2

  ! Distance between x and i1,j1,k1, and i2,j2,k2
  real :: dx1, dy1, dz1, dx2, dy2, dz2

  ! dl/B in physical units
  real :: InvBDl, RhoP_V(2)

  ! Debugging
  logical :: okdebug=.false.
  !--------------------------------------------------------------------------

  if(oktest_ray)write(*,'(a,3i4,3es12.4)')&
       'Starting follow_ray_block: me,iBlock,iRay,Xyz_D=',&
       iProc,iBlock,iRay,Xyz_D

  ! Convert initial position to block coordinates
  Dxyz_D = (/Dx_BLK(iBlock), Dy_BLK(iBlock), Dz_BLK(iBlock)/)
  x      = (Xyz_D - XyzStart_BLK(:,iBlock))/Dxyz_D + 1.0

  ! Set flag if checking on the ionosphere is necessary
  DoCheckInnerBc=Rmin_BLK(iBlock)<R_raytrace

  ! Set flag if checking for open rays is useful
  DoCheckOpen = .false.
!!!! any(ray(1,iRay,1:nI,1:nJ,1:nK,iBlock)==OPENRAY)

  ! Set the boundaries of the control volume in block coordinates
  ! We go out to the first ghost cell centers for sake of speed and to avoid
  ! problems at the boundaries
  xmin=(/   0.0,   0.0,   0.0/)
  xmax=(/nI+1.0,nJ+1.0,nK+1.0/)

  ! Go out to the block interface at the edges of the computational domain
  where(XyzStart_BLK(:,iBlock)+Dxyz_D*(xmax-1.0) > XyzMax_D)xmax = xmax - 0.5
  where(XyzStart_BLK(:,iBlock)+Dxyz_D*(xmin-1.0) < XyzMin_D)xmin = xmin + 0.5

  ! Initial value
  dl_next=sign(dl_max,1.5-iRay)

  ! Accuracy in terms of x in normalized coordinates
  dx_opt   = 0.01

  ! Length and maximum length of ray within control volume
  l        = 0
  lmax     = 10*maxval(xmax-xmin)
  nsegment = 0
  n_iono   = 0

  ! Integration loop
  FOLLOW: do
     ! Check if we are inside the ionosphere
     if(DoCheckInnerBc)then
        ! Convert x to real coordinates xx
        xx = XyzStart_BLK(:,iBlock) + Dxyz_D * (x-1.)
        r2 = sum(xx**2)

        if(r2<=R2_raytrace)then

           if(NameVectorField /= 'B')then
              Xyz_D = xx
              iFace=ray_iono_
              EXIT FOLLOW
           end if

           ! Try mapping down to rIonosphere if we haven't tried yet
           if(n_iono<1)then
              if(.not.follow_ray_iono())then
                 ! We did not hit the surface of the ionosphere
                 ! continue the integration
                 n_iono=n_iono+1
              else
                 if(oktest_ray)write(*,'(a,3i4,6es12.4)')&
                      'Inside R_raytrace at me,iBlock,nsegment,x,xx=',&
                      iProc,iBlock,nsegment,x,xx

                 r=sqrt(r2)
                 xx_ini = XyzStart_BLK(:,iBlock) + Dxyz_D*(x_ini-1.0)
                 r_ini=sqrt(sum(xx_ini**2))

                 ! The fraction of the last step inside body is estimated from 
                 ! the radii.
                 Fraction = (R_raytrace - r) / (r_ini - r)

                 ! Reduce ray length
                 Length = Length - Fraction * dLength

                 ! Recalculate position
                 x = x - Fraction*(x-x_ini)
                 
                 if(NameTask == 'integrate')then
                    ! Reduce integrals with the fraction of the last step 
                    if(oktest_ray)write(*,'(a,4es12.4)')&
                         'Before reduction InvBdl, RayIntegral_V=', InvBdl, &
                         RayIntegral_V(InvB_),RayIntegral_V(RhoInvB_:pInvB_)

                    ! Recalculate dLength/abs(B)
                    InvBDl = Fraction * InvBDl

                    ! Reduce field line volume
                    RayIntegral_V(InvB_) = RayIntegral_V(InvB_) - InvBDl

                    ! Reduce density and pressure integrals
                    RayIntegral_V(RhoInvB_:pInvB_) = &
                         RayIntegral_V(RhoInvB_:pInvB_) - InvBDl * RhoP_V

                    if(oktest_ray)then
                       write(*,'(a,4es12.4)')&
                            'After  reduction InvBdl, RayIntegral_V=',InvBdl, &
                            RayIntegral_V(InvB_),RayIntegral_V(RhoInvB_:pInvB_)

                       write(*,*)'Reduction at InvBDl,RhoP_V   =',InvBDl,RhoP_V
                       write(*,*)'Reduction r_ini,r,R_raytrace =',&
                            r_ini,r,R_raytrace
                    end if
                    
                 end if

                 ! Exit integration loop (xx was set by follow_ray_iono)
                 Xyz_D=xx
                 iFace=ray_iono_
                 EXIT FOLLOW
              end if
           end if
        end if
     end if
     ! Integrate with 2nd order scheme
     dl    = dl_next
     x_ini = x

     ! Half step
     call interpolate_b(x_ini, b_D, b_ini)
     x_mid = x_ini + 0.5*dl*b_ini

     ! Extract ray values using the interpolation coefficients around x_ini
     if(NameTask(1:7)=='extract')call ray_extract(.false., x_ini)

     STEP: do
        ! Full step
        b_mid = b_ini ! In case interpolation would gives zero vector
        call interpolate_b(x_mid, b_D, b_mid)

        ! Calculate the difference between 1st and 2nd order integration
        ! and take ratio relative to dx_opt
        dx_rel = abs(dl) * maxval(abs(b_mid-b_ini)) / dx_opt

        if(oktest_ray.and.okdebug)&
             write(*,*)'me,iBlock,x_mid,b_mid,dx_rel=', &
             iProc,iBlock,x_mid,b_mid,dx_rel

        ! Make sure that dl does not change more than a factor of 2 or 0.5
        dx_rel = max(0.5, min(2., dx_rel))

        if(dx_rel > 1.)then
           ! Not accurate enough, decrease dl if possible

           if(abs(dl) <= dl_min + dl_tiny)then
              ! Cannot reduce dl further
              dl_next=dl
              EXIT STEP
           end if

           dl = sign(max(dl_min,abs(dl)/(dx_rel+0.001)),dl)

           ! New mid point using the reduced dl
           x_mid = x_ini + 0.5*dl*b_ini

           if(oktest_ray.and.okdebug)&
                write(*,*)'new decreased dl: me,iBlock,dl=', &
                iProc,iBlock,dl
        else
           ! Too accurate, increase dl if possible
           if(abs(dl) < dl_max - dl_tiny)then

              dl_next = sign(min(dl_max,abs(dl)/sqrt(dx_rel)),dl)

              if(oktest_ray.and.okdebug)&
                   write(*,*)'new increased dl_next: me,iBlock,dl_next=', &
                   iProc,iBlock,dl_next

           end if

           EXIT STEP
        end if
     end do STEP

     ! Update position after the full step
     x = x_ini + b_mid*dl

     ! Update number of segments and the length in block coordinates
     nsegment = nsegment + 1
     l        = l + abs(dl)

     ! Step size in MH units  !!! Use simpler formula for cubic cells ???
     dLength = abs(dl)*sqrt( sum((b_mid*Dxyz_D)**2) )

     ! Update ray length
     Length  = Length + dLength

     if(NameTask == 'integrate')then

        ! Interpolate density and pressure
        ! Use the last indexes and distances already set in interpolate_b
        RhoP_V = dx1*(   dy1*(   dz1*Extra_VGB(:,i2,j2,k2,iBlock)   &
             +                   dz2*Extra_VGB(:,i2,j2,k1,iBlock))  &
             +           dy2*(   dz1*Extra_VGB(:,i2,j1,k2,iBlock)   &
             +                   dz2*Extra_VGB(:,i2,j1,k1,iBlock))) &
             +   dx2*(   dy1*(   dz1*Extra_VGB(:,i1,j2,k2,iBlock)   &
             +                   dz2*Extra_VGB(:,i1,j2,k1,iBlock))  &
             +           dy2*(   dz1*Extra_VGB(:,i1,j1,k2,iBlock)   &
             +                   dz2*Extra_VGB(:,i1,j1,k1,iBlock)))

        ! Calculate physical step size divided by physical field strength
        InvBDl = dLength / sqrt( sum(b_D**2) )

        ! Intgrate field line volume = \int dl/B
        RayIntegral_V(InvB_) = RayIntegral_V(InvB_) + InvBDl

        ! Integrate density and pressure = \int Rho dl/B and \int P dl/B
        RayIntegral_V(RhoInvB_:pInvB_) = RayIntegral_V(RhoInvB_:pInvB_) + &
             InvBDl * RhoP_V

        ! Check if we crossed the Z=0 plane
        xx_ini = XyzStart_BLK(:,iBlock) + Dxyz_D*(x_ini - 1.)
        xx     = XyzStart_BLK(:,iBlock) + Dxyz_D*(x     - 1.)

        ! What about multiple crossings !!!
        ! Should we check in which direction the equator is crossed !!!
        ! We could check for south-to-north crossings here !!!
        if(xx(3)*xx_ini(3)<=0)then

           ! Interpolate x and y
           dz1 = abs(xx_ini(3))/(abs(xx(3))+abs(xx_ini(3))); dz2 = 1.0 - dz1

           RayIntegral_V(Z0x_:Z0y_) = dz2*xx_ini(1:2) + dz1*xx(1:2)

           ! Assign Z0b_ as the middle point value of the magnetic field
           RayIntegral_V(Z0b_) = sqrt(sum(b_D**2))
           if(oktest_ray)then
              write(*,'(a,3es12.4)') &
                   'Found z=0 crossing at xx_ini=',xx_ini
              write(*,'(a,3es12.4)') &
                   'Found z=0 crossing at xx    =',xx
              write(*,'(a,3es12.4)')&
                   'RayIntegral_V(Z0x_:Z0b_)    =',RayIntegral_V(Z0x_:Z0b_)

              !write(*,'(a,2es12.4)')'Weights =',dz1,dz2
              !write(*,'(a,3es12.4)')'b_D = ',b_D
           end if
        end if
     end if

     if(oktest_ray.and.okdebug)&
          write(*,*)'me,iBlock,nsegment,l,x=', &
          iProc,iBlock,nsegment,l,x

     if(DoCheckOpen)then
        if(all(ray(1,iRay,i1:i2,j1:j2,k1:k2,iBlock)==OPENRAY))then
           iFace = ray_open_
           nOpen=nOpen+1
           EXIT FOLLOW
        end if
     end if

     ! Check if the ray hit the wall of the control volume
     if(any(x<xmin) .or. any(x>xmax))then

        ! Convert back to real coordinates
        Xyz_D = XyzStart_BLK(:,iBlock) + Dxyz_D*(x-1)

        if(any(Xyz_D < XyzMin_D) .or. any(Xyz_D > XyzMax_D))then
           iFace = ray_open_
        else
           iFace = ray_block_
        end if

        EXIT FOLLOW
     end if

     ! Check if we have integrated for too long
     if( l > lmax .or. Length > RayLengthMax )then
        ! Seems to be a closed loop within a block
        if(oktest_ray) &
             write(*,*)'CLOSED LOOP at me,iBlock,x,xx=',&
             iProc,iBlock,x,XyzStart_BLK(:,iBlock)+Dxyz_D*(x-1.0)

        iFace=ray_loop_
        EXIT FOLLOW
     end if

  end do FOLLOW

  ! Extract last point if ray is done. 
  ! The interpolation coefficients are not known.
  if(iFace /= ray_block_ .and. NameTask(1:7) == 'extract') &
       call ray_extract(.true., x)

  if(oktest_ray) then
     write(*,'(a,4i4)')&
          'Finished follow_ray_block at me,iBlock,nsegment,iFace=',&
          iProc,iBlock,nsegment,iFace
     write(*,'(a,i4,6es12.4)')&
          'Finished follow_ray_block at me,x,xx=',&
          iProc,x,XyzStart_BLK(:,iBlock)+Dxyz_D*(x-1.)
  end if

contains
  !========================================================================

  subroutine interpolate_b(x_D,b_D,Dir_D)

    ! Interpolate the magnetic field at normalized location x_D 
    ! and return the result in b_D. 
    ! The direction of b_D (normalized to a unit vector) is returned 
    ! in Dir_D if the magnitude of b_D is not (almost) zero.

    real, intent(in)   :: x_D(3)      ! location
    real, intent(out)  :: b_D(3)      ! interpolated magnetic field
    real, intent(inout):: Dir_D(3)    ! direction vector

    !LOCAL VARIABLES:
    real :: AbsB

    character (len=*), parameter :: NameSub='interpolate_b'

    !-------------------------------------------------------------------------

    ! Determine cell indices corresponding to location x_D
    i1=floor(x_D(1)); i2=i1+1
    j1=floor(x_D(2)); j2=j1+1
    k1=floor(x_D(3)); k2=k1+1

    ! Distance relative to the cell centers

    dx1 = x_D(1) - i1; dx2 = cOne - dx1
    dy1 = x_D(2) - j1; dy2 = cOne - dy1
    dz1 = x_D(3) - k1; dz2 = cOne - dz1

    ! Interpolate the magnetic field

    b_D = dx1*(   dy1*(   dz1*Bxyz_DGB(:,i2,j2,k2,iBlock)   &
         +                dz2*Bxyz_DGB(:,i2,j2,k1,iBlock))  &
         +        dy2*(   dz1*Bxyz_DGB(:,i2,j1,k2,iBlock)   &
         +                dz2*Bxyz_DGB(:,i2,j1,k1,iBlock))) &
         +dx2*(   dy1*(   dz1*Bxyz_DGB(:,i1,j2,k2,iBlock)   &
         +                dz2*Bxyz_DGB(:,i1,j2,k1,iBlock))  &
         +        dy2*(   dz1*Bxyz_DGB(:,i1,j1,k2,iBlock)   &
         +                dz2*Bxyz_DGB(:,i1,j1,k1,iBlock)))

    ! Stretch according to normalized coordinates
    Dir_D = b_D/Dxyz_D

    ! Normalize to unity
    AbsB = sqrt(sum(Dir_D**2))

    if(AbsB > cTiny)then
       Dir_D = Dir_D/AbsB
    else
       ! This is actually a big problem, and should not happen
       write(*,*)NameSub,' ERROR at iProc,iBlock,Xyz_D = ',iProc,iBlock,&
            XyzStart_BLK(:,iBlock)+Dxyz_D*(x_D-1.0)
       write(*,*)NameSub,' ERROR b_D=',b_D,' abs(B)=',AbsB
    end if

  end subroutine interpolate_b

  !===========================================================================

  logical function follow_ray_iono()

    ! Follow ray inside ionosphere starting from xx which is given in
    ! real coordinates and use analytic mapping.
    ! On return xx contains the final coordinates.
    ! Return true if it was successfully integrated down to rIonosphere,
    ! return false if the ray exited R_raytrace or too many integration 
    ! steps were done

    use CON_planet_field, ONLY: map_planet_field
    use CON_planet,       ONLY: get_planet
    use ModMain, ONLY: Time_Simulation

    integer :: iHemisphere
    real    :: x_D(3), DipoleStrength=0.0
    !---------------------------------------------------------------------
    if(DipoleStrength==0)call get_planet(DipoleStrengthOut=DipoleStrength)

    call map_planet_field(Time_Simulation, xx, TypeCoordSystem//' NORM', &
         rIonosphere, x_D, iHemisphere)

    if(iHemisphere==0)then
       write(*,*)'iHemisphere==0 for xx=',xx
       write(*,*)'iBlock, iRay=',iBlock,iRay
       call stop_mpi('ERROR in follow_ray_iono')
    end if

    if(iHemisphere*DipoleStrength*sign(1.0,1.5-iRay) < 0.0)then
       xx = x_D
       follow_ray_iono = .true.
    else
       follow_ray_iono = .false.
    end if

  end function follow_ray_iono

  !=========================================================================

  subroutine ray_extract(DoInterpolate, x_D)

    use ModIO, ONLY: iUnitLine_I
    logical, intent(in) :: DoInterpolate
    real, intent(in) :: x_D(3)
    real :: Xyz_D(3)
    character(len=*), parameter :: NameSub='ray_extract'
    !----------------------------------------------------------------------
    !!! if(DoInterpolate) calculating coefficients to be implemented

    ! Convert x_D to real coordinates
    Xyz_D = XyzStart_BLK(:,iBlock) + Dxyz_D*(x_D - 1.)

    select case(NameTask)
    case('extractfiles')
       write(iUnitLine_I(iStart_D(1)),'(100es18.10)') Length, Xyz_D
    case('extractfile')
       write(iUnitLine_I(1),'(100es18.10)') Length, Xyz_D, real(iStart_D(1))
    case('extract')
       call stop_mpi(NameSub//' ERROR task=extract is not implemented yet')
    case default
       call stop_mpi(NameSub//' ERROR invalid task='//NameTask)
    end select

  end subroutine ray_extract

end subroutine follow_ray_block

!============================================================================

subroutine ray_trace_sorted

  ! This subroutine is an experiment to sort blocks and cells such that
  ! open field lines can be found very fast. It works well for simple problems,
  ! but it does not seem to improve the performance for realistic grids

  use ModMain, ONLY: MaxBlock, nBlock, nI, nJ, nK, unusedBLK
  use ModPhysics, ONLY: SW_Bx, SW_By, SW_Bz
  use ModGeometry, ONLY: XyzMin_D, XyzMax_D, XyzStart_BLK
  use ModSort, ONLY: sort_quick

  implicit none

  integer :: iStart, iEnd, iStride, jStart, jEnd, jStride, &
       kStart, kEnd, kStride

  real    :: Weight_D(3)                 ! weights for the directions
  real    :: SortFunc_B(MaxBlock)        ! sorting function
  integer :: iBlockSorted_B(MaxBlock)    ! sorted block inxdexes

  ! index order for sorted blocks
  integer :: iSort, iSortStart, iSortEnd, iSortStride

  ! Indices corresponding to the starting point and directon of the ray
  integer :: i, j, k, iBlock, iRay

  !-------------------------------------------------------------------------

  ! Sort blocks according to the direction of the solar wind magnetic field
  ! so that open rays are found fast from already calculated ray values.

  ! Weight X, Y and Z according to the SW_Bx, SW_By, SW_Bz components
  ! The Y and Z directions are preferred to X (usually SW_Bx=0 anyways).
  Weight_D(1) = sign(1.0,SW_Bx)
  ! Select Y or Z direction to be the slowest changing value
  ! to maximize overlap
  if(abs(SW_By) > abs(SW_Bz))then
     Weight_D(2) = sign(100.0,SW_By)
     Weight_D(3) = sign( 10.0,SW_Bz)
  else
     Weight_D(2) = sign( 10.0,SW_By)
     Weight_D(3) = sign(100.0,SW_Bz)
  end if

  do iBlock=1,nBlock
     if(unusedBLK(iBlock))then
        SortFunc_B(iBlock) = -10000.0
     else
        SortFunc_B(iBlock) = sum(Weight_D*&
             (XyzStart_BLK(:,iBlock) - XyzMin_D)/(XyzMax_D - XyzMin_D))
     end if
  end do

  call sort_quick(nBlock,SortFunc_B,iBlockSorted_B)

  call barrier_mpi
  !CpuTimeStartRay = MPI_WTIME()

  ! Assign face ray values to cell centers

  !nOpen = 0
  do iRay=1,2

     if(iRay==1)then
        iSortStart=nBlock; iSortEnd=1; iSortStride=-1
     else
        iSortStart=1; iSortEnd=nBlock; iSortStride=1
     end if

     if(iRay==1 .eqv. SW_Bx >= 0.0)then
        iStart = nI; iEnd=1; iStride=-1
     else
        iStart = 1; iEnd=nK; iStride= 1
     end if

     if(iRay==1 .eqv. SW_By >= 0.0)then
        jStart = nJ; jEnd=1; jStride=-1
     else
        jStart = 1; jEnd=nJ; jStride= 1
     end if

     if(iRay==1 .eqv. SW_Bz >= 0.0)then
        kStart = nK; kEnd=1; kStride=-1
     else
        kStart = 1; kEnd=nK; kStride= 1
     end if

     do iSort = iSortStart, iSortEnd, iSortStride
        iBlock = iBlockSorted_B(iSort)

        do k = kStart, kEnd, kStride
           do j = jStart, jEnd, jStride
              do i = iStart, iEnd, iStride

              end do
           end do
        end do
     end do

  end do

end subroutine ray_trace_sorted

!============================================================================

subroutine integrate_ray_accurate(nLat, nLon, Lat_I, Lon_I, Radius)

  use CON_ray_trace, ONLY: ray_init
  use CON_planet_field, ONLY: map_planet_field
  use ModRaytrace
  use ModMain, ONLY: nBlock, Time_Simulation, iTest, jTest
  use ModAdvance,    ONLY: State_VGB, Rho_, p_, Bx_, By_, Bz_, &
       B0xCell_BLK,B0yCell_BLK,B0zCell_BLK
  use ModPhysics, ONLY: rBody
  use ModProcMH
  use ModMpi
  use ModNumConst,       ONLY: cDegToRad, cTiny
  use ModCoordTransform, ONLY: sph_to_xyz
  use ModUtilities,      ONLY: check_allocate
  use ModGeometry,       ONLY: XyzMax_D, XyzMin_D
  implicit none

  !INPUT ARGUMENTS:
  integer, intent(in):: nLat, nLon
  real,    intent(in):: Lat_I(nLat), Lon_I(nLon), Radius

  real    :: Theta, Phi, Lat, Lon, XyzIono_D(3), Xyz_D(3)
  integer :: iLat, iLon, iHemisphere
  integer :: iProcFound, iBlockFound, i, j, k

  integer :: iError
  logical :: DoTest, DoTestMe
  character(len=*), parameter :: NameSub = 'integrate_ray_accurate'
  !-------------------------------------------------------------------------

  call set_oktest(NameSub, DoTest, DoTestMe)

  if(DoTest)write(*,*)NameSub,' starting on iProc=',iProc

  call timing_start('integrate_ray')

  oktest_ray = .false.

  ! Initialize R_raytrace, R2_raytrace
  R_raytrace   = rBody
  R2_raytrace  = R_raytrace**2
  RayLengthMax = 2*sum(XyzMax_D - XyzMin_D)

  NameTask        = 'integrate'
  NameVectorField = 'B'

  ! (Re)initialize CON_ray_trace
  call ray_init(iComm)

  ! Copy magnetic field into Bxyz_DGB
  Bxyz_DGB(:,:,:,:,1:nBlock) = State_VGB(Bx_:Bz_,:,:,:,1:nBlock)

  ! Fill in all ghost cells (faces+edges+corners) without monotone restrict
  call message_pass_cells8(.false.,.false.,.false.,3,Bxyz_DGB)

  ! Add B0 for faster interpolation
  Bxyz_DGB(1,:,:,:,1:nBlock) = Bxyz_DGB(1,:,:,:,1:nBlock) &
       + B0xCell_BLK(:,:,:,1:nBlock)
  Bxyz_DGB(2,:,:,:,1:nBlock) = Bxyz_DGB(2,:,:,:,1:nBlock) &
       + B0yCell_BLK(:,:,:,1:nBlock)
  Bxyz_DGB(3,:,:,:,1:nBlock) = Bxyz_DGB(3,:,:,:,1:nBlock) &
       + B0zCell_BLK(:,:,:,1:nBlock)

  ! Copy density and pressure into Extra_VGB
  Extra_VGB(1,:,:,:,1:nBlock) = State_VGB(rho_,:,:,:,1:nBlock)
  Extra_VGB(2,:,:,:,1:nBlock) = State_VGB(p_  ,:,:,:,1:nBlock)

  ! Fill in all ghost cells (faces+edges+corners) without monotone restrict
  call message_pass_cells8(.false.,.false.,.false.,2,Extra_VGB)

  ! Initialize storage for the integrals
  allocate(&
       RayIntegral_VII(nRayIntegral,nLat,nLon), &
       RayResult_VII(nRayIntegral,nLat,nLon), STAT=iError)
  call check_allocate(iError,NameSub//' RayIntegral_VII,RayResult_VII')
  RayIntegral_VII = 0.0
  RayResult_VII   = 0.0

  ! Integrate rays starting from the latitude-longitude pairs defined
  ! by the arrays Lat_I, Lon_I
  do iLat = 1, nLat

     Lat = Lat_I(iLat)
     Theta = cDegToRad*(90.0 - Lat)     

     do iLon = 1, nLon

        Lon = Lon_I(iLon)
        Phi = cDegToRad*Lon

        ! Convert to SMG Cartesian coordinates on the surface of the ionosphere
        call sph_to_xyz(Radius, Theta, Phi, XyzIono_D)

        ! Map from the ionosphere to rBody
        call map_planet_field(time_simulation, XyzIono_D, 'SMG NORM', &
             rBody+cTiny, Xyz_D, iHemisphere)

        ! Check if the mapping is on the north hemisphere
        if(iHemisphere == 0)then
           !           write(*,*)NameSub,' point did not map to rBody, ',&
                !                'implement analytic integrals here! Lat, Lon=', Lat, Lon
           CYCLE
        elseif(iHemisphere == -1)then
           write(*,*)NameSub,' point mapped to the southern hemisphere! ',&
                'Lat, Lon=',Lat,Lon
           CYCLE
        end if

        ! Convert SMG to GSM (Note: these are identical for ideal axes)
        !! Xyz_D = matmul(GmIm_DD,Xyz_D)

        ! Find processor and block for the location
        call xyz_to_peblk(Xyz_D(1), Xyz_D(2), Xyz_D(3), &
             iProcFound, iBlockFound, .true., i, j, k)

        ! If location is on this PE, follow and integrate ray
        if(iProc == iProcFound)then

           if(iLat==iTest .and. iLon==jTest)then
              write(*,'(a,2i3,a,i3,a,i4)') &
                   'start of ray iLat, iLon=',iLat, iLon,&
                   ' found on iProc=',iProc,' iBlock=',iBlockFound
              write(*,'(a,2i4,2es12.4)')'iLon, iLat, Lon, Lat=',&
                   iLon, iLat, Lon, Lat
              write(*,'(a,3es12.4)')'XyzIono_D=',XyzIono_D
              write(*,'(a,3es12.4)')'Xyz_D    =',Xyz_D
           end if

           call follow_ray(2, (/iLat, iLon, 0, iBlockFound/), Xyz_D)

        end if
     end do
  end do

  ! Do remaining rays obtained from other PE-s
  call finish_ray

  if(DoTest.and.iTest<=nLat.and.jTest<=nLon) &
       write(*,*)NameSub,' iProc, RayIntegral_VII=',&
       iProc, RayIntegral_VII(:,iTest,jTest)

  ! Add up local integrals onto the root PE
  call MPI_reduce(RayIntegral_VII, RayResult_VII, nLat*nLon*nRayIntegral, &
       MPI_REAL, MPI_SUM, 0, iComm, iError)

  call timing_stop('integrate_ray')

end subroutine integrate_ray_accurate

!============================================================================

subroutine test_ray_integral

  use ModRayTrace, ONLY: RayResult_VII, InvB_, xEnd_, zEnd_, Length_, &
       nRayIntegral, xyz_to_latlon
  use ModMain,     ONLY: iTest, jTest
  use ModProcMH,   ONLY: iProc
  use ModIoUnit,   ONLY: UNITTMP_
  use ModNumConst, ONLY: cTiny
  implicit none

  integer, parameter :: nLat=50, nLon=50
  real :: Lat_I(nLat), Lon_I(nLon), Lat, Lon
  integer :: iLat, iLon
  integer :: iError
  character(len=*), parameter :: NameSub='test_ray_integral'
  !-------------------------------------------------------------------------

  write(*,*)NameSub,' starting on iProc=',iProc

  ! Initialize the spherical grid
  do iLat = 1, nLat
     Lat_I(iLat) = 50.0 + 40.0*(iLat-0.5)/nLat
  end do
  do iLon = 1, nLon
     Lon_I(iLon) = 360.0*(iLon-0.5)/nLon
  end do

  ! Integrate all points on the spherical grid
  call integrate_ray_accurate(nLat,nLon,Lat_I,Lon_I,1.0)

  ! Write out results into a file
  if(iProc==0)then

     ! Take logarithm of field line volume for better plotting ?
     RayResult_VII(InvB_,:,:)=alog10(RayResult_VII(InvB_,:,:)+cTiny)

     open(UNITTMP_,file='test_ray_integral.dat')
     write(UNITTMP_,"(a79)")'test-ray-integral_var22'
     write(UNITTMP_,"(i7,1pe13.5,3i3)")0, 0.0, 2, 1, nRayIntegral
     write(UNITTMP_,"(3i4)")nLat, nLon
     write(UNITTMP_,"(100(1pe13.5))")0.0
     write(UNITTMP_,"(a79)")&
          'Lon Lat Bvol Z0x Z0y Z0b Rho P LatEnd LonEnd Zend Length Param'

     do iLat=1,nLat
        Lat = Lat_I(iLat)
        do iLon=1,nLon
           Lon = Lon_I(iLon)

           call xyz_to_latlon(RayResult_VII(xEnd_:zEnd_,iLat,iLon))

           if(iLat == iTest .and. iLon == jTest)then
              write(*,'(a,a)')'iLon iLat Lon Lat ',&
                   'Bvol Z0x Z0y Z0b Rho P LatEnd LonEnd Zend Length'
              write(*,'(2i4,100(1es12.4))') iLon, iLat, Lon, Lat, &
                   RayResult_VII(:,iLat,iLon)
           end if

           write(UNITTMP_,"(100(1pe18.10))")Lon,Lat,RayResult_VII(:,iLat,iLon)
        end do
     end do
     close(UNITTMP_)
  end if

  ! Deallocate buffers ???
  ! deallocate(RayIntegral_VII, RayResult_VII)

  ! Clean up CON_ray_trace ???
  !call clean_ray

  call timing_show('integrate_ray',1)

  write(*,*)NameSub,' finished on iProc=',iProc
  call mpi_finalize(iError)
  stop

end subroutine test_ray_integral

!==============================================================================

subroutine write_plot_line(iFile)

  use ModProcMH,   ONLY: iProc, iComm
  use ModRayTrace, ONLY: NameTask, NameVectorField, &
       oktest_ray, R_Raytrace, R2_Raytrace, RayLengthMax, Bxyz_DGB
  use CON_ray_trace, ONLY: ray_init
  use ModAdvance,  ONLY: State_VGB, RhoUx_, RhoUy_, RhoUz_, Bx_, By_, Bz_, &
       B0xCell_BLK, B0yCell_BLK, B0zCell_BLK
  use ModPhysics,  ONLY: rBody
  use ModIO,       ONLY: NamePlotDir, plot_type, plot_form, Plot_,iUnitLine_I,&
       NameLine_I, nLine_I, XyzStartLine_DII, IsParallelLine_II, IsSingleLine_I
  use ModMain,     ONLY: n_step, time_simulation, nI, nJ, nK, nBlock, unusedBLK
  use ModGeometry, ONLY: XyzMax_D, XyzMin_D, Dx_BLK, Dy_BLK, Dz_BLK
  use ModIoUnit,   ONLY: io_unit_new, UnitTmp_
  use ModUtilities,ONLY: flush_unit

  implicit none

  integer, intent(in) :: iFile ! The file index of the plot file

  character(len=100) :: NameFile, NameStart
  integer            :: nHeaderFile

  real    :: Xyz_D(3), Dx2Inv, Dy2Inv, Dz2Inv
  integer :: iPlotFile, iProcFound, iBlockFound, iLine, iRay, iError

  integer :: i, j, k, iBlock

  character(len=*), parameter :: NameSub = 'write_plot_line'
  logical :: DoTest, DoTestMe
  !-------------------------------------------------------------------------
  call set_oktest(NameSub, DoTest, DoTestMe)

  iPlotFile = iFile - Plot_

  ! Initialize R_raytrace, R2_raytrace
  oktest_ray = .false.
  R_raytrace   = rBody
  R2_raytrace  = R_raytrace**2
  RayLengthMax = 2*sum(XyzMax_D - XyzMin_D)

  if(IsSingleLine_I(iPlotFile))then
     NameTask    = 'extractfiles'
     nHeaderFile = nLine_I(iPlotFile)
  else
     NameTask    = 'extractfile'
     nHeaderFile = 1
  end if

  ! (Re)initialize CON_ray_trace
  call ray_init(iComm)

  ! Fill in all ghost cells (faces+edges+corners) without monotone restrict
  !call message_pass_cells8(.false.,.false.,.false.,8,State_VGB)

  NameVectorField = NameLine_I(iPlotFile)
  select case(NameVectorField)
  case('B')
     ! Store B1+B0 for faster interpolation
     Bxyz_DGB(1,:,:,:,1:nBlock) = State_VGB(Bx_,:,:,:,1:nBlock) &
          + B0xCell_BLK(:,:,:,1:nBlock)
     Bxyz_DGB(2,:,:,:,1:nBlock) = State_VGB(By_,:,:,:,1:nBlock) &
          + B0yCell_BLK(:,:,:,1:nBlock)
     Bxyz_DGB(3,:,:,:,1:nBlock) = State_VGB(Bz_,:,:,:,1:nBlock) &
          + B0zCell_BLK(:,:,:,1:nBlock)
  case('U')
     ! Store momentum field (same as velocity field after normalization)
     Bxyz_DGB(1,:,:,:,1:nBlock) = State_VGB(RhoUx_,:,:,:,1:nBlock)
     Bxyz_DGB(2,:,:,:,1:nBlock) = State_VGB(RhoUy_,:,:,:,1:nBlock)
     Bxyz_DGB(3,:,:,:,1:nBlock) = State_VGB(RhoUz_,:,:,:,1:nBlock)
  case('J')
     ! Store current
     do iBlock=1,nBlock;
        if(unusedBLK(iBlock)) CYCLE
        Dx2Inv = 0.5/Dx_BLK(iBlock)
        Dy2Inv = 0.5/Dy_BLK(iBlock)
        Dz2Inv = 0.5/Dz_BLK(iBlock)

        do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
           Bxyz_DGB(1,i,j,k,iBlock) = &
                (State_VGB(Bz_,i,j+1,k,iBlock)-State_VGB(Bz_,i,j-1,k,iBlock)) &
                *Dy2Inv - &
                (State_VGB(By_,i,j,k+1,iBlock)-State_VGB(By_,i,j,k-1,iBlock)) &
                *Dz2Inv
           Bxyz_DGB(2,i,j,k,iBlock) = &
                (State_VGB(Bx_,i,j,k+1,iBlock)-State_VGB(Bx_,i,j,k-1,iBlock)) &
                *Dz2Inv - &
                (State_VGB(Bz_,i+1,j,k,iBlock)-State_VGB(Bz_,i-1,j,k,iBlock)) &
                *Dx2Inv
           Bxyz_DGB(3,i,j,k,iBlock) = &
                (State_VGB(By_,i+1,j,k,iBlock)-State_VGB(By_,i-1,j,k,iBlock)) &
                *Dx2Inv - &
                (State_VGB(Bx_,i,j+1,k,iBlock)-State_VGB(Bx_,i,j-1,k,iBlock)) &
                *Dy2Inv
        end do; end do; end do
     end do
  case default
     write(*,*)NameSub,' WARNING: for iFile=',iFile,' invalid NameLine=',&
          trim(NameVectorField),'!!!'
     RETURN
  end select

  write(NameStart,'(a,i1,a)') &
       trim(NamePlotDir)//trim(plot_type(iFile))//'_',iPlotFile, &
       '_'//NameVectorField

  ! Write header information
  do iLine = 1, nHeaderFile

     NameFile = NameStart
     if(NameTask == 'extractfiles' .and. nLine_I(iPlotFile)>1)then
        write(NameFile,'(a,i1)') trim(NameFile),iLine
     else
        NameFile = NameStart
     end if
     write(NameFile,'(a,i7.7,a)') trim(NameFile) // '_n',n_step,'.l'

     if(iProc==0)then
        open(UnitTmp_,file=NameFile)
        select case(plot_form(iFile))
        case('idl')
           if(NameTask == 'extractfiles')then
              write(UnitTmp_,'(a79)') NameVectorField//' line_var11'
              write(UnitTmp_,'(i7,1pe13.5,3i3)') n_step,time_simulation,1,1,3
              write(UnitTmp_,'(a6)')             '__nn__'
              write(UnitTmp_,'(es13.5)')         real(iLine)
              write(UnitTmp_,'(a79)')            'Length x y z iLine'
           else
              write(UnitTmp_,'(a79)') NameVectorField//' lines_var11'
              write(UnitTmp_,'(i7,1pe13.5,3i3)') n_step,time_simulation,1,1,4
              write(UnitTmp_,'(a6)')             '__nn__'
              write(UnitTmp_,'(es13.5)')         real(nLine_I(iPlotFile))
              write(UnitTmp_,'(a79)')            'Length x y z Index nLine'
           end if
        case('tec')
           if(NameTask == 'extractfiles')then
              write(UnitTmp_,'(a)')'TITLE ="'//NameVectorField//' line"'
              write(UnitTmp_,'(a)')'VARIABLES="Length","x","y","z"'
              write(UnitTmp_,'(a,i2.2,a)')'ZONE T="'// &
                   NameVectorField//' line ',iLine,'", '//'I=__nn__'
           else
              write(UnitTmp_,'(a)')'TITLE ="'//NameVectorField//' lines"'
              write(UnitTmp_,'(a)')'VARIABLES="Length","x","y","z","Index"'
              write(UnitTmp_,'(a,i2.2,a)')'ZONE T="'// &
                   NameVectorField//' ',nLine_I(iPlotFile),' lines", '// &
                   'I=__nn__'
           end if
        end select
        close(UnitTmp_)
     end if

     ! Open output files
     write(NameFile,'(a,i4.4,a)') NameFile(1:len_trim(NameFile)-2) // &
          '_pe',iProc,'.lst'
     if(DoTest)write(*,*)NameSub,': saving results into file ',trim(NameFile)

     iUnitLine_I(iLine) = io_unit_new();
     open(iUnitLine_I(iLine), file = NameFile)

  end do

  do iLine = 1, nLine_I(iPlotFile)
     Xyz_D = XyzStartLine_DII(:, iLine, iPlotFile)

     call xyz_to_peblk(Xyz_D(1), Xyz_D(2), Xyz_D(3), &
          iProcFound, iBlockFound, .true., i, j, k)

     if(iProc == iProcFound)then
        if(DoTest)write(*,*)NameSub,' follows ray ',iLine,&
             ' from iProc,iBlock,i,j,k=',iProcFound, iBlockFound, i, j, k
        if(IsParallelLine_II(iLine,iPlotFile))then
           iRay = 1
        else
           iRay = 2
        end if
        call follow_ray(iRay, (/iLine, 0, 0, iBlockFound/), Xyz_D)
     end if
  end do

  ! Do remaining rays obtained from other PE-s
  call finish_ray

  do iLine=1,nHeaderFile
     close(iUnitLine_I(iLine))
  end do

end subroutine write_plot_line

