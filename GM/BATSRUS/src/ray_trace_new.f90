!^CFG COPYRIGHT UM
!^CFG FILE RAYTRACE
subroutine OPTION_RAYTRACING(on,name)

  logical, intent(out) :: on
  character (len=40), intent(out) :: name

  on  =.true.
  name='RAY TRACING Toth 2.0'

end subroutine OPTION_RAYTRACING

!=============================================================================

subroutine ray_trace

  ! Trace field lines from cell centers to the outer or inner boundaries

  use ModProcMH
  use ModRaytrace
  use ModSort
  use CON_ray_trace, ONLY : ray_init
  use ModMain
  use ModNumConst, ONLY: cRadToDeg, cTiny, cZero, cOne
  use ModPhysics, ONLY: rBody, SW_Bx, SW_By, SW_Bz
  use ModAdvance, ONLY : State_VGB, Bx_, By_, Bz_, &
       B0xCell_BLK,B0yCell_BLK,B0zCell_BLK
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,r_BLK,true_cell,&
       XyzMin_D, XyzMax_D, XyzStart_BLK

  use ModMpi
  implicit none

  ! Arguments
  ! remember last call and the last grid number
  integer :: n_last=-1, iLastGrid=-1, iLastDecomposition=-1

  ! Indices corresponding to the starting point and directon of the ray
  integer :: i, j, k, iBlock, iRay

  integer :: iStart, iEnd, iStride, jStart, jEnd, jStride, &
       kStart, kEnd, kStride

  real    :: Weight_D(3)                 ! weights for the directions
  real    :: SortFunc_B(MaxBlock)        ! sorting function
  integer :: iBlockSorted_B(MaxBlock)    ! sorted block inxdexes

  ! index order for sorted blocks
  integer :: iSort, iSortStart, iSortEnd, iSortStride

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

  ! Fill in all ghost cells (faces+edges+corners) without monotone restrict
  call message_pass_cells8(.false.,.false.,.false.,3,&
       State_VGB(Bx_:Bz_,:,:,:,:))

  Bx_GB = State_VGB(Bx_,:,:,:,:) + B0xCell_BLK
  By_GB = State_VGB(By_,:,:,:,:) + B0yCell_BLK
  Bz_GB = State_VGB(Bz_,:,:,:,:) + B0zCell_BLK

  ! Initial values
  ray=NORAY

  if(oktest_me)write(*,*)'rayface normalized B'
  if(oktime.and.iProc==0)then
     write(*,'(a)',ADVANCE='NO') 'setup and normalization:'
     call timing_show('ray_trace',1)
  end if

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
  !!!
  !!! write(*,*)'Weight_D=',Weight_D

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
  TimeStartRay = MPI_WTIME()

  ! Assign face ray values to cell centers

  nOpen = 0
!  do iRay=1,2

  iRay=2
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

!        do k = kStart, kEnd, kStride
!           do j = jStart, jEnd, jStride
!              do i = iStart, iEnd, iStride

     do k = 1, nK; do j = 1, nJ; do i = 1, nI

!     do iSort = iSortStart, iSortEnd, iSortStride
!        iBlock = iBlockSorted_B(iSort)

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
                 call follow_ray(iRay, i, j, k, iBlock)

              end do ! i
           end do    ! j 
        end do       ! k
     end do          ! iBlock
  end do             ! iRay

  ! Do remaining rays passed from other PE-s
  call follow_ray(0,0,0,0,0)

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
      !!x_BLK(i,j,k,iBlock)==-2.125 .and. z_BLK(i,j,k,iBlock)==5.125 &
      !!     .and. abs(y_BLK(i,j,k,iBlock))==0.125

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
        write(*,*)'Loop ray found at iProc,iBlock,i,j,k,ray=',iProc,iBlock,i,j,k,Ray_DI
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

end subroutine ray_trace

!===========================================================================

subroutine follow_ray(iRayIn,iIn,jIn,kIn,iBlockIn)

  ! Follow ray starting at initial cell iIn,jIn,kIn in block iBlockIn
  ! in direction iRayIn until we hit the outer or inner boundary of
  ! the computational domain. The final position is stored into
  ! ModRayTrace::ray(:,iRayIn,iIn,jIn,kIn,iBlock).

  use ModRayTrace
  use CON_ray_trace

  use ModMain, ONLY: xTest, yTest, zTest
  use ModGeometry, ONLY: x_BLK,y_BLK,z_BLK,XyzStart_BLK,Dx_BLK
  use ModProcMH
  use ModKind

  use ModMpi
  implicit none


  integer, intent(in) :: iRayIn        ! ray direction, 0 if no ray is passed
  integer, intent(in) :: iIn, jIn, kIn ! cell index for starting position
  integer, intent(in) :: iBlockIn      ! block index for starting position

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

  real(Real8_) :: TimeNow

  !-----------------------------------------------------------------------

  if(iRayIn /= 0)then

     ! Store starting indexes and ray direction
     iStart = iIn; jStart = jIn; kStart = kIn; 
     iBlockStart = iBlockIn; iProcStart = iProc
     iRay   = iRayIn

     ! Store starting position
     XyzStart_D(1) = x_BLK(iIn,jIn,kIn,iBlockIn)
     XyzStart_D(2) = y_BLK(iIn,jIn,kIn,iBlockIn)
     XyzStart_D(3) = z_BLK(iIn,jIn,kIn,iBlockIn)

     ! Current position
     iBlockRay = iBlockStart
     XyzRay_D  = XyzStart_D

     oktest_ray = all(XyzStart_D==(/xTest,yTest,zTest/))
     if(oktest_ray)write(*,*)'Local ray at iProc,iBlock,i,j,k,iRay=',&
          iProc,iBlockRay,iIn,jIn,kIn,iRay

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
              oktest_ray = all(XyzStart_D==(/xTest,yTest,zTest/))
              if(IsParallel)then
                 iRay=1
              else
                 iRay=2
              end if
              if(oktest_ray)write(*,*)'Recv ray iProc,iRay,Done,XyzRay_D=',&
                   iProc,iRay,DoneRay,XyzRay_D

              if(DoneRay)then
                 ! Find the starting cell and check if it is local
                 call xyz_to_peblk(XyzStart_D(1),XyzStart_D(2),XyzStart_D(3),&
                      iProcStart,iBlockStart,.true.,iStart,jStart,kStart)

                 if(iProcStart /= iProc)call stop_mpi(&
                      'GM_ERROR in ray_trace: Done ray started from other PE')

                 ! Store the result into the ModRayTrace::ray
                 ray(:,iRay,iStart,jStart,kStart,iBlockStart)=XyzRay_D

                 if(oktest_ray)write(*,*)'Storing recv ray iProc,iRay,i,j,k,iBlock,ray=',&
                      iProc,iRay,iStart,jStart,kStart,iBlockStart,XyzRay_D

                 ! Get another ray from the others
                 CYCLE GETRAY
              else
                 ! Find block for the received ray
                 call xyz_to_peblk(XyzRay_D(1),XyzRay_D(2),XyzRay_D(3),&
                      jProc,iBlockRay,.false.,i,j,k)
                 if(jProc /= iProc)call stop_mpi(&
                      'GM_ERROR in ray_trace: Recvd ray is not in this PE')

                 if(oktest_ray)write(*,*)'Block for recv ray iProc,iBlock=',iProc,iBlockRay
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
           !TimeNow=MPI_WTIME()
           !write(*,*)'Final ray_exchange, iProc, time=',iProc,&
           !     TimeNow-TimeStartRay
           !TimeStartRay = TimeNow
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
        ! This PE is not done yet, so pass .false.

        !if(mod(iCountRay,DiCountRay)==0) then
        TimeNow=MPI_WTIME()
        if(TimeNow-TimeStartRay>DtExchangeRay)then
           call ray_exchange(iRayIn==0 .and. .not. IsFound, DoneAll)
           !write(*,*)'ray_exchange, iProc, time=',iProc,TimeNow-TimeStartRay
           TimeStartRay = TimeNow
           if(DoneAll) EXIT RAYS
        end if
     end if

  end do RAYS

contains

  !=========================================================================
  subroutine follow_this_ray

    ! Integrate ray
    BLOCK: do iCount = 1, MaxCount

       if(iCount == MaxCount)then
          write(*,*)'XyzStart_D  =',XyzStart_D
          write(*,*)'XyzRay_D    =',XyzRay_D
          write(*,*)'XyzStart_BLK=',XyzStart_BLK(:,iBlockRay)
          call stop_mpi('follow_ray integrated in more than MaxCount blocks')
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
             call ray_put(iProcStart,XyzStart_D,jProc,XyzRay_D,iRay==1,&
                  .false.)
             RETURN
          elseif(jBlock /= iBlockRay)then
             ! Continue the same ray in the next block
             iBlockRay = jBlock
             if(oktest_ray)write(*,*)'Continuing ray iProc,jBlock,iRay,Xyz=',&
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
          if(oktest_ray)write(*,*)&
               'follow_ray finished on the ionosphere at iProc,iRay,Xyz=',&
               iProc,iRay,XyzRay_D

       case default
          write(*,*)'Impossible value for iface=',iFace,&
               ' at XyzRay_D,iBlockRay=',XyzRay_D,iBlockRay
          call stop_mpi('GM_ERROR in follow_ray: impossible iFace value')
       end select

       ! The ray tracing of this ray is done if we got here
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

  ! Initial and mid point coordinates and bb field
  real, dimension(3) :: x_ini, x_mid, b_ini, b_mid, xx_ini
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
  !!!     xmin=(/   0.5,   0.5,   0.5/),&
  !!!     xmax=(/nI+0.5,nJ+0.5,nK+0.5/)
  !!! Maybe extend to 0 0 0 nI+1 nJ+1 nK+1 (but it did not help?!)

  ! Current position of ray in normalized and physical coordinates
  real, dimension(3) :: x, xx

  ! Radial distance and square of it: r2=sum(xx**2)
  real :: r2

  ! Cell indices corresponding to current or final x position
  integer :: i1,j1,k1,i2,j2,k2

  ! Distance between x and i1,j1,k1, and i2,j2,k2
  real :: dx1, dy1, dz1, dx2, dy2, dz2

  ! Debugging
  logical :: okdebug=.false.
  !--------------------------------------------------------------------------

  if(oktest_ray)&
       write(*,*)'follow_ray: me,iBlock,Xyz_D,iRay=',iProc,iBlock,Xyz_D,iRay

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

           if(oktest_ray)write(*,*)&
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
                 xx=(xx*(r_ini-rIonosphere)+xx_ini*(rIonosphere-r)) &
                      /(r_ini-r)
                 ! Normalize xx in radial direction
                 xx=rIonosphere*xx/sqrt(sum(xx**2))
                 Xyz_D=xx
                 iFace=ray_iono_
              end if
              EXIT FOLLOW
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
     call interpolate_bb(x_ini,b_ini)
     x_mid=x_ini+0.5*dl*b_ini

     STEP: do
        ! Full step
        call interpolate_bb(x_mid,b_mid)

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

           if(abs(dl)<=dl_min+dl_tiny)then
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
        if(oktest_ray)then
           write(*,*)'CLOSED LOOP at me,iBlock,x,xx=',&
                iProc,iBlock,x,&
                x_BLK(1,1,1,iBlock)+dx_BLK(iBlock)*(x(1)-1.),&
                y_BLK(1,1,1,iBlock)+dy_BLK(iBlock)*(x(2)-1.),&
                z_BLK(1,1,1,iBlock)+dz_BLK(iBlock)*(x(3)-1.)
        end if

        iFace=ray_loop_
        EXIT FOLLOW
     end if

  end do FOLLOW

  if(oktest_ray) &
       write(*,*)'Finished follow_ray at me,iBlock,nsegment,iFace,x,xx=',&
       iProc,iBlock,nsegment,iFace,x,&
       x_BLK(1,1,1,iBlock)+dx_BLK(iBlock)*(x(1)-1.),&
       y_BLK(1,1,1,iBlock)+dy_BLK(iBlock)*(x(2)-1.),&
       z_BLK(1,1,1,iBlock)+dz_BLK(iBlock)*(x(3)-1.)

contains
  !===========================================================================

  subroutine interpolate_bb(qx,qb)

    ! Obtain normalized bb field at normalized location qx and put it into qb

    real, intent(in) :: qx(3)
    real, intent(out):: qb(3)
    real :: qbD

    !-------------------------------------------------------------------------

    ! Determine cell indices corresponding to location qx
    i1=floor(qx(1)); i2=i1+1
    j1=floor(qx(2)); j2=j1+1
    k1=floor(qx(3)); k2=k1+1

    ! Distance relative to the cell centers

    dx1 = qx(1) - i1; dx2 = cOne - dx1
    dy1 = qx(2) - j1; dy2 = cOne - dy1
    dz1 = qx(3) - k1; dz2 = cOne - dz1

    ! Add in interpolated B1 values and take aspect ratios into account

    qb(1) = interpolate_bb1(Bx_GB)/dx_BLK(iBlock)
    qb(2) = interpolate_bb1(By_GB)/dy_BLK(iBlock)
    qb(3) = interpolate_bb1(Bz_GB)/dz_BLK(iBlock)

    ! Normalize
    qbD = sqrt(qb(1)**2 + qb(2)**2 + qb(3)**2)

    if(qbD>cTiny)then
       qb=qb/qbD
    else
       qb=0.
    end if

  end subroutine interpolate_bb

  !===========================================================================

  subroutine interpolate_bb0(qx,qb)

    ! Obtain normalized bb field at normalized location qx and put it into qb

    real, intent(in) :: qx(3)
    real, intent(out):: qb(3)
    real :: qbD

    !-------------------------------------------------------------------------

    ! Get B0 values for location

    xx(1)=x_BLK(1,1,1,iBlock)+dx_BLK(iBlock)*(qx(1)-1.)
    xx(2)=y_BLK(1,1,1,iBlock)+dy_BLK(iBlock)*(qx(2)-1.)
    xx(3)=z_BLK(1,1,1,iBlock)+dz_BLK(iBlock)*(qx(3)-1.)

    call get_b0(xx(1),xx(2),xx(3),qb)

    ! Determine cell indices corresponding to location qx
    i1=floor(qx(1)); i2=i1+1
    j1=floor(qx(2)); j2=j1+1
    k1=floor(qx(3)); k2=k1+1

    ! Distance relative to the cell centers

    dx1=qx(1)-i1; dx2=1.-dx1
    dy1=qx(2)-j1; dy2=1.-dy1
    dz1=qx(3)-k1; dz2=1.-dz1

    ! Add in interpolated B1 values and take aspect ratios into account

    qb(1)=(qb(1)+interpolate_bb1(State_VGB(Bx_,:,:,:,:)))/dx_BLK(iBlock)
    qb(2)=(qb(2)+interpolate_bb1(State_VGB(By_,:,:,:,:)))/dy_BLK(iBlock)
    qb(3)=(qb(3)+interpolate_bb1(State_VGB(Bz_,:,:,:,:)))/dz_BLK(iBlock)
!!!    end if

    ! Normalize
    qbD=sqrt(qb(1)**2 + qb(2)**2 + qb(3)**2)

    if(qbD>cTiny)then
       qb=qb/qbD
    else
       qb=0.
    end if

  end subroutine interpolate_bb0

  !===========================================================================

  real function interpolate_bb1(qbb)

    real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), &
         intent(in):: qbb

    !-------------------------------------------------------------------------

    ! Bilinear interpolation in 3D

    interpolate_bb1=&
         dx1*(   dy1*(   dz1*qbb(i2,j2,k2,iBlock)+&
         dz2*qbb(i2,j2,k1,iBlock))+&
         dy2*(   dz1*qbb(i2,j1,k2,iBlock)+&
         dz2*qbb(i2,j1,k1,iBlock)))+&
         dx2*(   dy1*(   dz1*qbb(i1,j2,k2,iBlock)+&
         dz2*qbb(i1,j2,k1,iBlock))+&
         dy2*(   dz1*qbb(i1,j1,k2,iBlock)+&
         dz2*qbb(i1,j1,k1,iBlock)))

  end function interpolate_bb1

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

!=============================================================================

subroutine convfaces2latlon

  call stop_mpi('convfaces2latlon to be written')

end subroutine convfaces2latlon

!=============================================================================

subroutine integrate_ray

  call stop_mpi('integrate_ray to be written')

end subroutine integrate_ray
