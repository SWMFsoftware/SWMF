!^CFG COPYRIGHT UM
!^CFG FILE RAYTRACE

subroutine ray_trace_accurate

  ! Trace field lines from cell centers to the outer or inner boundaries

  use ModProcMH
  use ModRaytrace
  use CON_ray_trace, ONLY: ray_init
  use ModMain
  use ModNumConst,   ONLY: cRadToDeg, cTiny, cZero, cOne
  use ModPhysics,    ONLY: rBody
  use ModAdvance,    ONLY: State_VGB, Bx_, By_, Bz_, &
       B0xCell_BLK,B0yCell_BLK,B0zCell_BLK
  use ModGeometry,   ONLY: x_BLK,y_BLK,z_BLK,r_BLK,true_cell

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

  if(  n_last             == n_step   .and. &
       iLastGrid          == iNewGrid .and. &
       iLastDecomposition == iNewDecomposition) RETURN

  if(n_last == -1)then
     ! Initialize constants for the first time
     R_raytrace = rBody
     if(iProc==0)write(*,*)'Setting R_raytrace=',R_raytrace
     R2_raytrace = R_raytrace**2

     ! Initialize CON_ray_trace
     call ray_init(iComm)
  end if

  ! Remember this call
  n_last=n_step; iLastGrid = iNewGrid; iLastDecomposition = iNewDecomposition

  call set_oktest('time_ray_trace',oktime,oktime_me)
  call timing_start('ray_trace')
  if(oktime)call timing_reset('ray_pass',2)

  ! Copy magnetic field into Bxyz_DGB
  Bxyz_DGB(:,:,:,:,1:nBlock) = &
       State_VGB(Bx_:Bz_,:,:,:,1:nBlock)

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

        oktest_ray = &
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
  call follow_ray(0,(/0,0,0,0/),(/0.0,0.0,0.0/))

!!! write(*,*)'iProc=',iProc,': open ray interpolations=',nOpen

  ! Convert x, y, z to latitude and longitude, and status
  do iBlock=1,nBlock
     if(unusedBLK(iBlock)) CYCLE
     do k=1,nK; do j=1,nJ; do i=1,nI
        call convert_ray(ray(:,:,i,j,k,iBlock))
     end do; end do; end do
  end do

  call timing_stop('ray_trace')

  if(oktime.and.iProc==0)then
     write(*,'(a)',ADVANCE='NO') 'Total ray tracing time:'
     call timing_show('ray_trace',1)
  end if
  call barrier_mpi
  if(oktest_me)write(*,*)'ray_trace completed.'

contains

  !=========================================================================

  subroutine convert_ray(Ray_DI)

    real, intent(inout) :: Ray_DI(3,2)
    real :: RayTmp_DI(3,2)
    logical :: DoTest

    DoTest = .false.

    ! Store input ray values into RayTmp_DI
    RayTmp_DI = Ray_DI

    do iRay=1,2

       if(DoTest)then
          write(*,*)'iProc,iBlock,i,j,k,iRay  =',iProc,iBlock,i,j,k,iRay
          write(*,*)'iBlock,iRay,Ray_DI(:,iRay),CLOSED=',&
               iBlock,iRay,Ray_DI(:,iRay),CLOSEDRAY
       end if

       ! Check if this direction is closed or not
       if(RayTmp_DI(1,iRay)>CLOSEDRAY)then
          ! Make sure that asin will work, -1<=RayTmp_DI(3,iRay)<=1
          RayTmp_DI(3,iRay) = max(-cOne+cTiny, RayTmp_DI(3,iRay))
          RayTmp_DI(3,iRay) = min( cOne-cTiny, RayTmp_DI(3,iRay))

          ! Calculate  -90 < theta=asin(z)  <  90
          Ray_DI(1,iRay) = cRadToDeg * asin(RayTmp_DI(3,iRay))

          if(DoTest)write(*,*)'iBlock,iRay,theta=',iBlock,iRay,Ray_DI(1,iRay)

          ! Calculate -180 < phi=atan2(y,z) < 180
          if(  abs(RayTmp_DI(1,iRay))<cTiny .and. &
               abs(RayTmp_DI(2,iRay))<cTiny) &
               RayTmp_DI(1,iRay)=1.

          Ray_DI(2,iRay) = &
               cRadToDeg * atan2(RayTmp_DI(2,iRay),RayTmp_DI(1,iRay))

          if(DoTest)write(*,*)'iBlock,iRay,phi  =',iBlock,iRay,Ray_DI(2,iRay)


          ! Get rid of negative phi angles
          if(Ray_DI(2,iRay) < cZero) Ray_DI(2,iRay) = Ray_DI(2,iRay)+360.

          if(DoTest)write(*,*)'iBlock,iRay,phi+ =',iBlock,iRay,Ray_DI(2,iRay)

       else
          ! Impossible values
          Ray_DI(1,iRay) = -100.
          Ray_DI(2,iRay) = -200.

          if(DoTest)write(*,*)'iBlock,iRay,Ray_DI- =',&
               iBlock,iRay,Ray_DI(1:2,iRay)
       endif

    end do ! iRay

    ! Calculate and store ray status in ray(3,1...)

    if(RayTmp_DI(1,1)>CLOSEDRAY .and. RayTmp_DI(1,2)>CLOSEDRAY)then
       Ray_DI(3,1)=3.      ! Fully closed
    elseif(RayTmp_DI(1,1)>CLOSEDRAY .and. RayTmp_DI(1,2)==OPENRAY)then
       Ray_DI(3,1)=2.      ! Half closed in positive direction
    elseif(RayTmp_DI(1,2)>CLOSEDRAY .and. RayTmp_DI(1,1)==OPENRAY)then
       Ray_DI(3,1)=1.      ! Half closed in negative direction
    elseif(RayTmp_DI(1,1)==OPENRAY .and. RayTmp_DI(1,2)==OPENRAY) then
       Ray_DI(3,1)=0.      ! Fully open
    elseif(RayTmp_DI(1,1)==BODYRAY)then
       Ray_DI(3,1)=-1.     ! Cells inside body
    elseif(RayTmp_DI(1,1)==LOOPRAY .and.  RayTmp_DI(1,2)==LOOPRAY) then
       Ray_DI(3,1)=-2.     ! Loop ray within block
       write(*,*)'Loop ray found at iProc,iBlock,i,j,k,ray=',&
            iProc,iBlock,i,j,k,Ray_DI
    else
       Ray_DI(3,1)=-3.     ! Strange status
       !        if(  x_BLK(i,j,k,iBlock)==xTest.and.&
       !             y_BLK(i,j,k,iBlock)==yTest.and.&
       !             z_BLK(i,j,k,iBlock)==zTest) &
       write(*,*)'Strange ray found at iProc,iBlock,i,j,k,ray=',&
            iProc,iBlock,i,j,k,RayTmp_DI
    end if

    if(DoTest)write(*,*)'iBlock,iRay,status=',iBlock,iRay,Ray_DI(3,1)

  end subroutine convert_ray

end subroutine ray_trace_accurate

!===========================================================================

subroutine follow_ray(iRayIn,i_D,XyzIn_D)

  ! Follow ray starting at initial cell iIn,jIn,kIn in block iBlockIn
  ! in direction iRayIn until we hit the outer or inner boundary of
  ! the computational domain. The final position is stored into
  ! ModRayTrace::ray(:,iRayIn,iIn,jIn,kIn,iBlock).

  use ModRayTrace
  use CON_ray_trace

  use ModMain, ONLY: xTest, yTest, zTest, iTest, jTest
  use ModGeometry, ONLY: x_BLK,y_BLK,z_BLK,XyzStart_BLK,Dx_BLK
  use ModProcMH
  use ModKind

  use ModMpi
  implicit none

  !INPUT ARGUMENTS:
  integer, intent(in) :: iRayIn        ! ray direction, 0 if no ray is passed
  integer, intent(in) :: i_D(4)        ! cell/block index for starting position
  real,    intent(in) :: XyzIn_D(3)    ! coordinates of starting position

  !LOCAL VARIABLES:
  ! Cell, block and PE indexes for initial position and ray direction
  integer :: iStart, jStart, kStart, iBlockStart, iProcStart, iRay
  real    :: XyzStart_D(3)

  ! Current position of the ray
  integer :: iBlockRay
  real    :: XyzRay_D(3)

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

  !-----------------------------------------------------------------------

  if(iRayIn /= 0)then

     ! Store starting indexes and ray direction
     iStart = i_D(1); jStart = i_D(2); kStart = i_D(3); 
     iBlockStart = i_D(4); iProcStart = iProc
     iRay   = iRayIn

     if(DoIntegrate)then
        ! Store initial RCM grid indexes as starting position
        XyzStart_D = (/ real(iStart), real(jStart), 0.0 /)
        oktest_ray = iStart == iTest .and. jStart == jTest
     else
        ! Store starting Cartesian position
        XyzStart_D = XyzIn_D
        oktest_ray = all(XyzStart_D==(/xTest,yTest,zTest/))
     end if

     ! Current position
     iBlockRay = iBlockStart
     XyzRay_D  = XyzIn_D

     if(oktest_ray)write(*,'(a,6i4,3es12.4)')&
          'Local ray at iProc,i_D,iRay,XyzIn_D=',iProc,i_D,iRay,XyzIn_D

  end if

  ! If iRayIn==0 there are no more local rays to follow so get from other PEs
  DoGet = iRayIn==0
  IsFound = .true.

  RAYS: do

     if(DoGet)then
        GETRAY: do
           call ray_get(IsFound,iProcStart,XyzStart_D,XyzRay_D,IsParallel,&
                DoneRay)

           if(IsFound)then
              if(DoIntegrate)then
                 oktest_ray = all(nint(XyzStart_D(1:2)) == (/iTest,jTest/))
              else
                 oktest_ray = all(XyzStart_D==(/xTest,yTest,zTest/))
              end if
              if(IsParallel)then
                 iRay=1
              else
                 iRay=2
              end if
              if(oktest_ray)write(*,*)'Recv ray iProc,iRay,Done,XyzRay_D=',&
                   iProc,iRay,DoneRay,XyzRay_D

              if(DoneRay)then

                 if(DoIntegrate)then
                    ! store integrals and whether this is a closed ray
                    call store_integral(XyzRay_D(1) > CLOSEDRAY)

                    CYCLE GETRAY
                 end if

                 ! Find the starting cell and check if it is local
                 call xyz_to_peblk(XyzStart_D(1),XyzStart_D(2),XyzStart_D(3),&
                      iProcStart,iBlockStart,.true.,iStart,jStart,kStart)

                 if(iProcStart /= iProc)call stop_mpi(&
                      'GM_ERROR in ray_trace: Done ray started from other PE')

                 ! Store the result into the ModRayTrace::ray
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
    if(DoIntegrate)RayIntegral_V = 0.0

    ! Follow the ray through the local blocks
    BLOCK: do iCount = 1, MaxCount

       if(iCount == MaxCount)then
          write(*,*)'XyzStart_D  =',XyzStart_D
          write(*,*)'XyzRay_D    =',XyzRay_D
          write(*,*)'XyzStart_BLK=',XyzStart_BLK(:,iBlockRay)
          call stop_mpi('follow_ray passed through more than MaxCount blocks')
       end if

       call follow_ray_block(iRay,iBlockRay,XyzRay_D,iFace)

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
             ! Pass .true., because this may well be a closed ray
             if(DoIntegrate)call store_integral(.true.)

             call ray_put(iProcStart,XyzStart_D,jProc,XyzRay_D,iRay==1,&
                  .false.)
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
             write(*,*)'ERROR for XyzStart_D  =',XyzStart_D
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

       ! The ray tracing of this ray is done if we got here
       if(DoIntegrate)then
          ! Store integrals and the information of being a closed ray
          call store_integral(iFace == ray_iono_)
          EXIT BLOCK
       end if

       if(iProcStart == iProc)then

          ! If went through other PE-s, find starting cell/block from position
          if(DoGet)then
             call xyz_to_peblk(&
                  XyzStart_D(1),XyzStart_D(2),XyzStart_D(3),&
                  jProc,iBlockStart,.true.,iStart,jStart,kStart)
             if(jProc /= iProc)call stop_mpi(&
                  'GM_ERROR: in follow_ray: wrong PE for start position')

             if(oktest_ray)write(*,*) &
                  'Found start pos iProc,iBlock,i,j,k,iRay=',&
                  jProc,iBlockStart,iStart,jStart,kStart,iRay

             if(oktest_ray)write(*,*) &
                  'Found start pos iProc,iBlock,i,j,k,iRay=',&
                  jProc,iBlockStart,iStart,jStart,kStart,iRay
          end if

          ! Store result locally
          ray(:,iRay,iStart,jStart,kStart,iBlockStart)=XyzRay_D

          if(oktest_ray)write(*,*) &
               'Storing into iProc,iBlock,i,j,k,iRay,Xyz=',&
               iProc,iBlockStart,iStart,jStart,kStart,iRay,XyzRay_D

       else
          ! Send back result to iProcStart. 
          call ray_put(iProcStart,XyzStart_D,iProc,XyzRay_D,iRay==1,.true.)

          if(oktest_ray)write(*,*) &
               'Send result iProc,iProcStart,iRay,Xyz=',&
               iProc,iProcStart,iRay,XyzRay_D

       end if
       EXIT BLOCK

    end do BLOCK

  end subroutine follow_this_ray

  !===========================================================================
  subroutine store_integral(IsClosedRay)

    ! Store integrals of this ray into the 

    logical, intent(in) :: IsClosedRay

    integer :: iLat, iLon

    iLat = nint(XyzStart_D(1))
    iLon = nint(XyzStart_D(2))

    RayIntegral_VII(1:nRayIntegral,iLat,iLon) = &
         RayIntegral_VII(1:nRayIntegral,iLat,iLon) + RayIntegral_V

    if(.not.IsClosedRay) RayIntegral_VII(ClosedRay_,iLat,iLon) = -1.0

  end subroutine store_integral

end subroutine follow_ray

!==========================================================================
subroutine follow_ray_block(iRay,iBlock,Xyz_D,iFace)

  ! Follow ray starting at initial position Xyz_D inside block iBlock
  ! in direction iRay until we hit the wall of the block or the ionosphere. 
  ! Return Xyz_D with the final position.
  ! Return iFace = 1..6 if the ray hit the east,west,south,north,bot,top walls
  ! Return ray_iono_    if the ray hit the ionosphere
  ! Return ray_loop_    if the ray did not hit anything 
  ! Return ray_body_    if the ray goes into or is inside a body
  ! Return ray_open_    if the ray goes outside the computational box

  use ModRayTrace
  use ModProcMH
  use ModNumConst, ONLY: cTiny, cOne
  use ModMain, ONLY: TypeCoordSystem, nI, nJ, nK, nBLK
  use ModGeometry, ONLY: XyzStart_BLK, XyzMax_D, XyzMin_D, &
       x_BLK,y_BLK,z_BLK, Dx_BLK, Dy_BLK, Dz_BLK, rMin_BLK
  use ModAdvance, ONLY : State_VGB, Bx_, By_, Bz_

  implicit none

  ! Arguments

  integer, intent(in) :: iRay
  integer, intent(in) :: iBlock
  real, intent(inout) :: Xyz_D(3)
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
  real :: l, lmax, dl, dl_next

  ! Step size limits
  real, parameter :: dl_max=1.0, dl_min=0.05, dl_tiny=1.e-6

  ! counter for ray integration
  integer :: nsegment 

  ! True if Rmin_BLK < R_raytrace
  logical :: check_inside

  ! True if the block already containes open rays
  logical :: DoCheckOpen

  ! Counter for entering follow_ray_iono
  integer :: n_iono

  ! Control volume limits in local coordinates
  real, dimension(3), parameter :: &
       xmin=(/   0.0,   0.0,   0.0/),&
       xmax=(/nI+1.0,nJ+1.0,nK+1.0/)

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

  ! Set flag if checking on the ionosphere is necessary
  check_inside=Rmin_BLK(iBlock)<R_raytrace

  ! Set flag if checking for open rays is useful
  DoCheckOpen = .false.
  !!!! any(ray(1,iRay,1:nI,1:nJ,1:nK,iBlock)==OPENRAY)

  ! Initial value
  dl_next=sign(dl_max,1.5-iRay)

  ! Accuracy in terms of x in normalized coordinates
  dx_opt=0.01

  ! Length and maximum length of ray within control volume
  l=0
  lmax=10*maxval(xmax-xmin)
  nsegment=0
  n_iono=0

  ! Convert initial position to block coordinates
  Dxyz_D = (/Dx_BLK(iBlock), Dy_BLK(iBlock), Dz_BLK(iBlock)/)
  x      = (Xyz_D - XyzStart_BLK(:,iBlock))/Dxyz_D + 1.0

  ! Integration loop
  FOLLOW: do
     ! Check if we are inside the ionosphere
     if(check_inside)then
        ! Convert x to real coordinates xx

        xx = XyzStart_BLK(:,iBlock) + Dxyz_D * (x-1.)

        r2=sum(xx**2)

        if(r2<=R2_raytrace)then

           if(oktest_ray)write(*,'(a,3i4,6es12.4)')&
                'Inside R_raytrace at me,iBlock,nsegment,x,xx=',&
                iProc,iBlock,nsegment,x,xx

           if(r2 <= rIonosphere2)then
              write(*,*)'WARNING r<rIonosphere should not happen!'
              if(nsegment==0)then
                 iFace=ray_body_
                 if(oktest_ray)write(*,*)&
                      'Initial point inside rIonosphere at me,iBlock,xx=',&
                      iProc,iBlock,xx
              else
                 r=sqrt(r2)
                 xx_ini = XyzStart_BLK(:,iBlock)+Dxyz_D*(x_ini-1.0)
                 r_ini=sqrt(sum(xx_ini**2))
                 ! Interpolate to the surface linearly along last segment
                 xx=(xx*(r_ini-rIonosphere)+xx_ini*(rIonosphere-r))/(r_ini-r)
                 ! Normalize xx in radial direction
                 xx=rIonosphere*xx/sqrt(sum(xx**2))
                 Xyz_D=xx
                 iFace=ray_iono_
              end if
              EXIT FOLLOW
           end if

           if(DoIntegrate)then
              r=sqrt(r2)
              xx_ini = XyzStart_BLK(:,iBlock)+Dxyz_D*(x_ini-1.0)
              r_ini=sqrt(sum(xx_ini**2))

              if(oktest_ray)write(*,'(a,4es12.4)')&
                   'Before reduction InvBdl, RayIntegral_V=', InvBdl, &
                   RayIntegral_V(InvB_),RayIntegral_V(RhoInvB_:pInvB_)

              ! Reduce integrals with the fraction of the last step which is
              ! inside rBody. This fraction is estimated from the radii.
              InvBDl = InvBDl * (R_raytrace - r) / (r_ini - r)

              ! Reduce field line volume
              RayIntegral_V(InvB_) = RayIntegral_V(InvB_) - InvBDl

              ! Reduce density and pressure integrals
              RayIntegral_V(RhoInvB_:pInvB_) = RayIntegral_V(RhoInvB_:pInvB_) &
                   - InvBDl * RhoP_V

              if(oktest_ray)then
                 write(*,'(a,4es12.4)')&
                      'After  reduction InvBdl, RayIntegral_V=',InvBdl, &
                      RayIntegral_V(InvB_),RayIntegral_V(RhoInvB_:pInvB_)

                 write(*,*)'Reduction at InvBDl,RhoP_V   =',InvBDl,RhoP_V
                 write(*,*)'Reduction r_ini,r,R_raytrace =',r_ini,r,R_raytrace
              end if

           end if

           ! Try mapping down to rIonosphere if we haven't tried yet
           if(n_iono<1)then
              if(follow_ray_iono())then
                 Xyz_D=xx
                 iFace=ray_iono_
                 EXIT FOLLOW
              else
                 ! We did not hit the surface of the ionosphere
                 ! continue the integration
                 n_iono=n_iono+1
              end if
           end if
        end if
     end if

     ! Integrate with 2nd order scheme
     dl=dl_next
     x_ini=x

     ! Half step
     call interpolate_b(x_ini,b_D,b_ini)
     x_mid=x_ini+0.5*dl*b_ini

     STEP: do
        ! Full step
        call interpolate_b(x_mid,b_D,b_mid)

        ! Calculate the difference between 1st and 2nd order integration
        ! and take ratio relative to dx_opt
        dx_rel=abs(dl)*maxval(abs(b_mid-b_ini))/dx_opt

        if(oktest_ray.and.okdebug)&
             write(*,*)'me,iBlock,x_mid,b_mid,dx_rel=', &
             iProc,iBlock,x_mid,b_mid,dx_rel

        ! Make sure that dl does not change more than a factor of 2 or 0.5
        dx_rel=max(0.5,min(2.,dx_rel))

        if(dx_rel>1.)then
           ! Not accurate enough, decrease dl if possible

           if(abs(dl) <= dl_min + dl_tiny)then
              ! Cannot reduce dl further
              dl_next=dl
              EXIT STEP
           end if

           dl = sign(max(dl_min,abs(dl)/(dx_rel+0.001)),dl)

           ! New mid point using the reduced dl
           x_mid=x_ini+0.5*dl*b_ini

           if(oktest_ray.and.okdebug)&
                write(*,*)'new decreased dl: me,iBlock,dl=', &
                iProc,iBlock,dl
        else
           ! Too accurate, increase dl if possible
           if(abs(dl)<dl_max-dl_tiny)then

              dl_next = sign(min(dl_max,abs(dl)/sqrt(dx_rel)),dl)

              if(oktest_ray.and.okdebug)&
                   write(*,*)'new increased dl_next: me,iBlock,dl_next=', &
                   iProc,iBlock,dl_next

           end if

           EXIT STEP
        end if
     end do STEP

     x=x_ini+b_mid*dl

     if(DoIntegrate)then

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
        InvBDl =  abs(dl) * sqrt( sum((b_mid*Dxyz_D)**2)/sum(b_D**2) )

        ! Intgrate field line volume = \int dl/B
        RayIntegral_V(InvB_) = RayIntegral_V(InvB_) + InvBDl

        ! Integrate density and pressure = \int Rho dl/B and \int P dl/B
        RayIntegral_V(RhoInvB_:pInvB_) = RayIntegral_V(RhoInvB_:pInvB_) + &
             InvBDl * RhoP_V

        ! Check if we crossed the Z=0 plane
        xx_ini = XyzStart_BLK(:,iBlock)+Dxyz_D*(x_ini-1.)
        xx     = XyzStart_BLK(:,iBlock)+Dxyz_D*(x    -1.)

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

     nsegment=nsegment+1
     l=l+abs(dl)

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
     if(l>lmax)then
        ! Seems to be a closed loop within a block
        if(oktest_ray) &
             write(*,*)'CLOSED LOOP at me,iBlock,x,xx=',&
             iProc,iBlock,x,XyzStart_BLK(:,iBlock)+Dxyz_D*(x-1.0)

        iFace=ray_loop_
        EXIT FOLLOW
     end if

  end do FOLLOW

  if(oktest_ray) then
     write(*,'(a,4i4)')&
          'Finished follow_ray_block at me,iBlock,nsegment,iFace=',&
          iProc,iBlock,nsegment,iFace
     write(*,'(a,i4,6es12.4)')&
          'Finished follow_ray_block at me,x,xx=',&
          iProc,x,XyzStart_BLK(:,iBlock)+Dxyz_D*(x-1.)
  end if

contains
  !===========================================================================

  subroutine interpolate_b(x_D,b_D,Dir_D)

    ! Obtain normalized bb field at normalized location x_D and put it into qb

    real, intent(in) :: x_D(3)      ! location
    real, intent(out):: b_D(3)      ! interpolated magnetic field
    real, intent(out):: Dir_D(3)    ! direction vector

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
       write(*,*)NameSub,' ERROR at iProc,iBlock,Xyz_D = ',iProc,iBlock,Xyz_D
       call stop_mpi(NameSub//' ERROR: magnetic field is zero')
       !!! Dir_D = 0.0
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

end subroutine follow_ray_block

!============================================================================

subroutine ray_trace_sorted

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

subroutine test_ray_integral

  use CON_ray_trace, ONLY: ray_init
  use ModRaytrace
  use ModMain, ONLY: nBlock
  use ModAdvance,    ONLY: State_VGB, Rho_, p_, Bx_, By_, Bz_, &
       B0xCell_BLK,B0yCell_BLK,B0zCell_BLK
  use ModPhysics, ONLY: rBody
  use ModProcMH
  use ModMpi
  use ModNumConst, ONLY: cDegToRad
  use ModCoordTransform, ONLY: sph_to_xyz
  use ModIoUnit, ONLY: UNITTMP_
  implicit none

  integer, parameter :: nLat=50, nLon=50

  real :: Result_VII(0:nRayIntegral, nLat, nLon)
  real :: Theta, Phi, Lat, Lon, Xyz_D(3)
  integer :: iProcFound, iBlockFound, i, j, k, iLat, iLon

  integer :: iError
  !-------------------------------------------------------------------------

  write(*,*)'Starting test_ray_integral iProc=',iProc

  oktest_ray = .true.

  ! Initialize R_raytrace, R2_raytrace
  R_raytrace = rBody
  if(iProc==0)write(*,*)'Setting R_raytrace=',R_raytrace
  R2_raytrace = R_raytrace**2

  ! Initialize CON_ray_trace
  call ray_init(iComm)

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

  ! Copy density and pressure into Extra_VGB
  Extra_VGB(1,:,:,:,1:nBlock) = State_VGB(rho_,:,:,:,1:nBlock)
  Extra_VGB(2,:,:,:,1:nBlock) = State_VGB(p_  ,:,:,:,1:nBlock)

  ! Fill in all ghost cells (faces+edges+corners) without monotone restrict
  call message_pass_cells8(.false.,.false.,.false.,2,Extra_VGB)

  ! Initialize storage for the integrals
  allocate(RayIntegral_VII(0:nRayIntegral,nLat,nLon))
  RayIntegral_VII = 0.0
  DoIntegrate = .true.

  ! Test cell in RCM grid and corresponding XYZ location
  do iLat = 1, nLat
     Lat = 20.0 + (50.0*iLat)/nLat
     Theta = cDegToRad*(90.0 - Lat)     
     do iLon = 1, nLon
        Lon = (360.0*iLon)/nLon
        Phi = cDegToRad*Lon
        call sph_to_xyz(rBody+0.1, Theta, Phi, Xyz_D)

        ! Find location
        call xyz_to_peblk(Xyz_D(1), Xyz_D(2), Xyz_D(3), &
             iProcFound, iBlockFound, .true., i, j, k)

        ! If location is on this PE, follow and integrate ray
        if(iProc == iProcFound)then
!           write(*,'(a,2i3,a,i3,a,i4)') &
!                'start of ray iLat, iLon=',iLat, iLon,&
!                ' found on iProc=',iProc,' iBlock=',iBlockFound
!           write(*,'(a,5es12.4)')'Lon, Lat, Xyz_D=',Lon, Lat, Xyz_D
           call follow_ray(2, (/iLat, iLon, 0, iBlockFound/), Xyz_D)
        end if
     end do
  end do

  ! Do remaining rays obtained from other PE-s
  call follow_ray(0, (/0, 0, 0, 0/), (/ 0., 0., 0. /))

  ! write(*,*)'iProc, RayIntegral_VII=',iProc, RayIntegral_VII

  ! Add up local integrals onto the root PE
  call MPI_reduce(RayIntegral_VII, Result_VII, nLat*nLon*(nRayIntegral+1), &
       MPI_REAL, MPI_SUM, 0, iComm, iError)

  ! Write out results
  if(iProc==0)then

     ! Take logarithm of field line volume for better plotting ?
     Result_VII(InvB_,:,:)=alog10(Result_VII(InvB_,:,:))

     open(UNITTMP_,file='test_ray_integral.dat')
     write(UNITTMP_,"(a79)")'test-ray-integral_var22'
     write(UNITTMP_,"(i7,1pe13.5,3i3)")0, 0.0, 2, 1, nRayIntegral+1
     write(UNITTMP_,"(3i4)")nLat, nLon
     write(UNITTMP_,"(100(1pe13.5))")0.0
     write(UNITTMP_,"(a79)")'Lon Lat Closed Bvol Z0x Z0y Z0b Rho P nothing'

!     write(*,*)'iLon iLat Lon Lat Closed Bvol Z0x Z0y Z0b Rho P'

     do iLat=1,nLat
        do iLon=1,nLon
           Lat = 20.0 + (50.0*iLat)/nLat
           Lon = (360.0 * iLon)/nLon

!           write(*,'(2i4,100(1es12.4))') iLon, iLat, Lat, Lon, &
!                Result_VII(:,iLat,iLon)

           write(UNITTMP_,"(100(1pe18.10))")Lon, Lat, Result_VII(:,iLat,iLon)
        end do
     end do
     close(UNITTMP_)
  end if

  ! Deallocate buffers ???
  deallocate(RayIntegral_VII)

  ! Clean up CON_ray_trace ???
  !call clean_ray

  write(*,*)'Finished test_ray_integral iProc=',iProc
  call mpi_finalize(iError)
  stop

end subroutine test_ray_integral
