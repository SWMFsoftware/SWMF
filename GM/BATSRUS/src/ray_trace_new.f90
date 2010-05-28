!^CFG COPYRIGHT UM
!^CFG FILE RAYTRACE

! The main subroutines in this file are
! subroutine ray_trace_accurate        - trace all rays starting from 3D MHD grid
! subroutine integrate_ray_accurate    - integrate rays starting from 2D IM grid
! subroutine integrate_ray_accurate_1d - integrate rays starting from 1D list of points
! subroutine write_plot_line           - extract lines into plot file(s)

subroutine ray_trace_accurate

  ! Trace field lines from cell centers to the outer or inner boundaries

  use ModProcMH
  use ModRaytrace
  use CON_ray_trace,  ONLY: ray_init
  use ModMain
  use ModAdvance,     ONLY: State_VGB, Bx_, Bz_, B0_DGB
  use ModGeometry,    ONLY: x_BLK, y_BLK, z_BLK, r_BLK, true_cell, &
       XyzMax_D, XyzMin_D, x1, x2, y1, y2, z1, z2, UseCovariant
  use ModMessagePass, ONLY: message_pass_dir
  use ModMpi

  implicit none

  ! Indices corresponding to the starting point and directon of the ray
  integer :: i, j, k, iBlock, iRay

  ! Testing and timing
  logical :: okTest, okTestMe, okTime, okTimeMe

  !----------------------------------------------------------------------------

  call set_oktest('ray_trace',okTest,okTestMe)

  ! Initialize constants
  if(UseOldMethodOfRayTrace .and. .not.UseCovariant)then
     RayLengthMax = 4.*sum(XyzMax_D - XyzMin_D)
  else
     RayLengthMax = 4.*(abs(x2-x1) + abs(y2-y1) + abs(z2-z1))
  end if

  DoTraceRay     = .true.
  DoMapRay       = .false.
  DoIntegrateRay = .false.
  DoExtractRay   = .false.
  nRay_D         = (/ nI, nJ, nK, nBlock /)
  NameVectorField = 'B'

  ! (Re)initialize CON_ray_trace
  call ray_init(iComm)

  call set_oktest('time_ray_trace',okTime,okTimeMe)
  if(okTime)call timing_reset('ray_pass',2)

  ! Fill in all ghost cells
  call message_pass_dir(iDirMin=1, iDirMax=3, Width=2, SendCorners=.true., &
       ProlongOrder=1, nVar=nVar, Sol_VGB=State_VGB)

  ! Copy magnetic field into Bxyz_DGB
  do iBlock = 1, nBlock
     if(unusedBLK(iBlock))CYCLE
     Bxyz_DGB(:,:,:,:,iBlock) = State_VGB(Bx_:Bz_,:,:,:,iBlock)
     ! Add B0
     if(UseB0) Bxyz_DGB(:,:,:,:,iBlock) = &
          Bxyz_DGB(:,:,:,:,iBlock) + B0_DGB(:,:,:,:,iBlock)
  end do

  ! Initial values
  ray=NORAY

  if(okTestMe)write(*,*)'rayface normalized B'
  if(okTime.and.iProc==0)then
     write(*,'(a)',ADVANCE='NO') 'setup and normalization:'
     call timing_show('ray_trace',1)
  end if

  ! This loop order seems to give optimal speed
  CpuTimeStartRay = MPI_WTIME();
  do k = 1, nK; do j = 1, nJ; do i = 1, nI
     do iBlock = 1, nBlock
        if(unusedBLK(iBlock))CYCLE

        oktest_ray = okTest .and. &
             x_BLK(i,j,k,iBlock)==xTest .and. &
             y_BLK(i,j,k,iBlock)==yTest .and. &
             z_BLK(i,j,k,iBlock)==zTest

        do iRay=1,2

           ! Short cut for inner and false cells
           if(R_BLK(i,j,k,iBlock) < rIonosphere .or. &
                .not.true_cell(i,j,k,iBlock))then
              ray(:,:,i,j,k,iBlock)=BODYRAY
              if(oktest_ray)write(*,*)'Shortcut BODYRAY iProc,iRay=',iProc,iRay
              CYCLE
           end if

           if(oktest_ray)write(*,*)'calling follow_ray iProc,iRay=',iProc,iRay

           ! Follow ray in direction iRay
           call follow_ray(iRay, (/i,j,k,iBlock/), (/ x_BLK(i,j,k,iBlock), &
                y_BLK(i,j,k,iBlock), z_BLK(i,j,k,iBlock) /) )

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

        ! Comment out these statements as they spew thousands of lines with IM coupling
        !if(ray(3,1,i,j,k,iBlock)==-2.) write(*,*) &
        !     'Loop ray found at iProc,iBlock,i,j,k,ray=',&
        !     iProc,iBlock,i,j,k,ray(:,:,i,j,k,iBlock)

        !if(ray(3,1,i,j,k,iBlock)==-3.) write(*,*) &
        !     'Strange ray found at iProc,iBlock,i,j,k,ray=',&
        !     iProc,iBlock,i,j,k,ray(:,:,i,j,k,iBlock)
     end do; end do; end do
  end do

  if(okTestMe)write(*,*)'ray lat, lon, status=',&
       ray(:,:,iTest,jTest,kTest,BlkTest)

  ! Return ghost cells to values before trace started
  call exchange_messages

  if(okTime.and.iProc==0)then
     write(*,'(a)',ADVANCE='NO') 'Total ray tracing time:'
     call timing_show('ray_trace',1)
  end if

  if(okTestMe)write(*,*)'ray_trace completed.'

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
  ! arrays defined in ModRayTrace or into files based on the logicals 
  ! in ModRaytrace (more than one of these can be true):
  !
  ! If DoTraceRay, follow the ray from cell centers of the 3D AMR grid, 
  !    and save the final position into
  !    ModRayTrace::ray(:,iRayIn,i_D(1),i_D(2),i_D(3),i_D(4)) on the 
  !    processor that started the ray trace.
  !
  ! If DoMapRay, map the rays down to the ionosphere, save spherical
  !    coordinates (in SMG) into
  !    ModRayTrace::RayMap_DSII(4,i_D(1),i_D(2),i_D(3))
  !
  ! If DoIntegrateRay, do integration along the rays and
  !    save the integrals into ModRayTrace::RayIntegral_VII(i_D(1),i_D(2))
  !
  ! If DoExtractRay, extract data along the rays, collect and sort it
  !    In this case the rays are indexed with i_D(1).
  !
  !EOP

  use ModRayTrace
  use CON_ray_trace, ONLY: ray_exchange, ray_get, ray_put

  use ModMain,     ONLY: iTest, jTest, kTest, BlkTest, ProcTest,DoMultiFluidIMCoupling
  use ModGeometry, ONLY: XyzStart_BLK, Dx_BLK
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
  integer :: iCountRay = 0

  real(Real8_) :: CpuTimeNow

  character(len=*), parameter :: NameSub='follow_ray'

  logical :: DoTest = .false., DoTestMe = .false.
  !-----------------------------------------------------------------------

  ! call set_oktest(NameSub, DoTest, DoTestMe)

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
                 if(.not.DoTraceRay)then
                    write(*,*)NameSub,' WARNING ',&
                         'received DoneRay=T for DoTraceRay = .false. !'
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
    if(DoIntegrateRay)RayIntegral_V = 0.0

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
             if(DoIntegrateRay)call store_integral(.false.)

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
             write(*,*)'ERROR iBlockRay=jBlock=',iBlockRay,jBlock
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
       if(DoIntegrateRay)call store_integral(.true.)

       if(DoMapRay)then
          if(.not.allocated(RayMapLocal_DSII))then
             if(allocated(RayMap_DSII)) deallocate(RayMap_DSII)
             allocate(RayMap_DSII(3,nRay_D(1),nRay_D(2),nRay_D(3)))
             allocate(RayMapLocal_DSII(3,nRay_D(1),nRay_D(2),nRay_D(3)))
             RayMapLocal_DSII = 0.0
          end if
          RayMapLocal_DSII(:,iStart_D(1),iStart_D(2),iStart_D(3)) = XyzRay_D
       end if

       ! Nothing more to do if not tracing
       if(.not.DoTraceRay) EXIT BLOCK

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

    RayIntegral_VII(InvB_:pInvB_,iLat,iLon) = &
         RayIntegral_VII(InvB_:pInvB_,iLat,iLon) + RayIntegral_V(InvB_:pInvB_)

    if(DoMultiFluidIMCoupling)then
       RayIntegral_VII(HpRhoInvB_:OppInvB_,iLat,iLon) = &
            RayIntegral_VII(HpRhoInvB_:OppInvB_,iLat,iLon) + &
            RayIntegral_V(HpRhoInvB_:OppInvB_)
    endif

    if(DoneRay)then
       RayIntegral_VII(xEnd_:zEnd_,iLat,iLon) = XyzRay_D
       RayIntegral_VII(Length_,iLat,iLon)     = RayLength
    end if
  end subroutine store_integral

  subroutine set_oktest_ray

    if(DoIntegrateRay)then
       ! Test the ray starting from a given Lat-Lon grid point
       oktest_ray = DoTest .and. all(iStart_D(1:2) == (/iLatTest,iLonTest/))
    else if(DoTraceRay)then
       ! Test the ray starting from a given grid cell
       oktest_ray = DoTest .and. iProcStart == ProcTest .and. &
            all(iStart_D == (/iTest,jTest,kTest,BlkTest/))
    else
       ! Check the ray indexed in line plot files.
       oktest_ray = DoTest .and. iStart_D(1) == iTest
    end if

  end subroutine set_oktest_ray

end subroutine follow_ray

!==========================================================================
subroutine follow_ray_block(iStart_D,iRay,iBlock,XyzInOut_D,Length,iFace)

  !DESCRIPTION:
  ! Follow ray identified by index array iStart_D, 
  ! starting at initial position XyzInOut_D inside block iBlock,
  ! in direction iRay until we hit the wall of the block or the ionosphere. 
  ! Return XyzInOut_D with the final position. 
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
  use ModNumConst, ONLY: cTiny
  use ModMain, ONLY: TypeCoordSystem, nI, nJ, nK
  use ModGeometry, ONLY: XyzStart_BLK, XyzMax_D, XyzMin_D,  UseCovariant, &
       Dx_BLK, Dy_BLK, Dz_BLK, rMin_BLK, x_BLK,y_BLK,z_BLK, x1,x2,y1,y2,z1,z2
  use CON_planet, ONLY: DipoleStrength
  use ModMain,    ONLY: DoMultiFluidIMCoupling
  use ModMultiFLuid
  implicit none

  ! Arguments

  integer, intent(in) :: iStart_D(4)
  integer, intent(in) :: iRay
  integer, intent(in) :: iBlock
  real, intent(inout) :: XyzInOut_D(3)
  real, intent(inout) :: Length
  integer, intent(out):: iFace

  ! Local variables

  ! Block size
  real :: Dxyz_D(3)

  ! initial/mid/current points of: IJK and XYZ coordinates
  real, dimension(3) :: IjkIni_D, IjkMid_D, IjkCur_D, XyzIni_D, XyzMid_D, XyzCur_D

  ! General coordinates and reference Ijk
  real, dimension(3) :: Gen_D, Ijk_D

  ! Direction of B field, true interpolated field
  real, dimension(3) :: bNormIni_D, bNormMid_D, b_D

  ! SM coordinates
  real, dimension(3) :: XyzSMIni_D, XyzSMCur_D

  ! Radial distance from origin and square
  real :: rCur, r2Cur, rIni

  ! dx is the difference between 1st and 2nd order RK to estimate accuracy
  ! dxOpt is the required accuracy, dxRel=dx/dxOpt
  real :: dxRel, dxOpt

  ! Ray step size, step length, next step size 
  real :: dl, dlp, dLength, dlNext

  ! Fraction of the last step inside the ionosphere
  real :: Fraction

  ! Step size limits
  real :: dlMax, dlMin, dlTiny

  ! counter for ray integration
  integer :: nSegment
  integer :: nSegmentMax=10*(nI+nJ+nK)

  ! True if Rmin_BLK < R_raytrace
  logical :: DoCheckInnerBc

  ! True if the block already containes open rays
  logical :: DoCheckOpen

  ! Counter for entering follow_ray_iono
  integer :: nIono

  ! Control volume limits in local coordinates
  real, dimension(3) :: xmin, xmax

  ! Cell indices corresponding to current or final Ijk position
  integer :: i1,j1,k1,i2,j2,k2

  ! Distance between Ijk and i1,j1,k1, and i2,j2,k2
  real :: dx1, dy1, dz1, dx2, dy2, dz2

  ! dl/B in physical units
  real :: InvBDl, RhoP_V(2), RhoPMulti_V(4)

  ! Debugging
  logical :: okdebug=.false.

  logical :: IsWall

  integer :: iIonSecond
  !--------------------------------------------------------------------------

  iIonSecond = min(IonFirst_+1, IonLast_)

  if(oktest_ray)write(*,'(a,3i4,3es12.4)')&
       'Starting follow_ray_block: me,iBlock,iRay,XyzInOut_D=',&
       iProc,iBlock,iRay,XyzInOut_D

  ! Store local block deltas
  Dxyz_D  = (/Dx_BLK(iBlock), Dy_BLK(iBlock), Dz_BLK(iBlock)/)

  ! Convert initial position to block coordinates
  XyzCur_D = XyzInOut_D
  call xyz_to_ijk(XyzCur_D,IjkCur_D,iBlock, &
       XyzCur_D,XyzStart_BLK(:,iBlock),Dxyz_D)

  ! Set flag if checking on the ionosphere is necessary
  if(UseOldMethodOfRayTrace .and. .not.UseCovariant)then
     DoCheckInnerBc = Rmin_BLK(iBlock) < R_raytrace + sum(Dxyz_D)
  else
     DoCheckInnerBc = Rmin_BLK(iBlock) < 1.2*R_raytrace
  end if

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
  if(UseCovariant)then
     xmin(2)=0.0;  xmax(2)=nJ+1.0
     xmin(3)=0.0;  xmax(3)=nK+1.0
  end if

  ! Step size limits
  if(UseOldMethodOfRayTrace .and. .not.UseCovariant)then
     dlMax = 1.0
     dlMin = 0.05
     dlTiny= 1.e-6
  else
     dlMax =( abs(x_BLK(nI,nJ,nK,iBlock)-x_BLK(1,1,1,iBlock)) &
          +   abs(y_BLK(nI,nJ,nK,iBlock)-y_BLK(1,1,1,iBlock)) &
          +   abs(z_BLK(nI,nJ,nK,iBlock)-z_BLK(1,1,1,iBlock)) ) &
          /(nI + nJ + nK - 3)
     dlMin = dlMax*0.05
     dlTiny= dlMax*1.e-6
  end if

  ! Initial value
  dlNext=sign(dlMax,1.5-iRay)

  ! Accuracy in terms of a kind of normalized coordinates
  dxOpt = 0.01*dlMax

  ! Reference Ijk
  Ijk_D = (/ nI/2, nJ/2, nK/2 /)

  ! Length and maximum length of ray within control volume
  nSegment = 0
  nIono    = 0

  IsWall=.false.

  ! Integration loop
  FOLLOW: do

     ! Integrate with 2nd order scheme
     dl    = dlNext
     IjkIni_D = IjkCur_D
     XyzIni_D = XyzCur_D

     ! Half step
     call interpolate_b(IjkIni_D, b_D, bNormIni_D)
     if(UseOldMethodOfRayTrace .and. .not.UseCovariant)then
        IjkMid_D = IjkIni_D + 0.5*dl*bNormIni_D
        XyzMid_D = XyzStart_BLK(:,iBlock) + Dxyz_D*(IjkMid_D - 1.)
     else
        HALF: do
           ! Try a half step in XYZ space (and get IJK from it)
           XyzMid_D = XyzIni_D + 0.5*dl*bNormIni_D
           call xyz_to_ijk(XyzMid_D, IjkMid_D, iBlock, &
                XyzIni_D, XyzStart_BLK(:,iBlock), Dxyz_D)

           ! Check if it stepped too far, cut step if needed
           if(any(IjkMid_D<(xmin-0.5)) .or. any(IjkMid_D>(xmax+0.5)))then
              ! Step too far, reduce and try again
              dl = 0.5*dl

              if(abs(dl) < dlMin)then
                 ! Cannot reduce dl further
                 dl = 0.0
                 ! Obtain a point outside the block by mirroring the block
                 ! center Ijk_D to the starting location of this step IjkIni_D
                 IjkMid_D = 2*IjkIni_D - Ijk_D
                 ! Reduce length of Ijk_D --> IjkMid_D vector to end 
                 ! something like a 10th of a cell outside the block
                 dlp = 1.1*(1.-maxval(max(xmin-IjkMid_D,IjkMid_D-xmax) &
                      /(abs(IjkMid_D-IjkIni_D)+dlTiny)))
                 IjkMid_D=IjkIni_D+dlp*(IjkMid_D-IjkIni_D)

                 ! Make sure that IjkMid_D is just outside the control volume
                 IjkMid_D=max(xmin-.1,IjkMid_D)
                 IjkMid_D=min(xmax+.1,IjkMid_D)
                 call interpolate_xyz(IjkMid_D,XyzMid_D)
                 call interpolate_b(IjkMid_D, b_D, bNormMid_D)
                 IjkCur_D=IjkMid_D; XyzCur_D=XyzMid_D

                 ! We exited the block and have a good location to continued from
                 IsWall=.true.
                 EXIT HALF
              end if
           else
              !Step was OK, continue
              EXIT HALF
           end if
        end do HALF
     end if

     ! Extract ray values using around IjkIni_D
     if(DoExtractRay)call ray_extract(IjkIni_D,XyzIni_D)

     STEP: do
        if(IsWall)EXIT STEP

        ! Full step
        bNormMid_D = bNormIni_D ! In case interpolation would give zero vector
        call interpolate_b(IjkMid_D, b_D, bNormMid_D)

        ! Calculate the difference between 1st and 2nd order integration
        ! and take ratio relative to dxOpt
        dxRel = abs(dl) * maxval(abs(bNormMid_D-bNormIni_D)) / dxOpt

        if(oktest_ray.and.okdebug)&
             write(*,*)'me,iBlock,IjkMid_D,bNormMid_D,dxRel=', &
             iProc,iBlock,IjkMid_D,bNormMid_D,dxRel

        ! Make sure that dl does not change more than a factor of 2 or 0.5
        dxRel = max(0.5, min(2., dxRel))

        if(dxRel > 1.)then
           ! Not accurate enough, decrease dl if possible

           if(abs(dl) <= dlMin + dlTiny)then
              ! Cannot reduce dl further
              dlNext=dl
              EXIT STEP
           end if

           dl = sign(max(dlMin,abs(dl)/(dxRel+0.001)),dl)

           ! New mid point using the reduced dl
           if(UseOldMethodOfRayTrace .and. .not.UseCovariant)then
              IjkMid_D = IjkIni_D + 0.5*dl*bNormIni_D
              XyzMid_D = XyzStart_BLK(:,iBlock) + Dxyz_D*(IjkMid_D - 1.)
           else
              HALF2: do
                 ! Try new half step in XYZ space (and get IJK from it)
                 XyzMid_D = XyzIni_D + 0.5*dl*bNormIni_D
                 call xyz_to_ijk(XyzMid_D,IjkMid_D,iBlock, &
                      XyzIni_D,XyzStart_BLK(:,iBlock),Dxyz_D)

                 ! Check if it stepped too far, cut step if needed
                 if(any(IjkMid_D<(xmin-0.5)) .or. any(IjkMid_D>(xmax+0.5)))then
                    ! Step too far, reduce and try again
                    dl=0.5*dl

                    if(abs(dl)<dlMin)then
                       ! Cannot reduce dl further
                       dl=0.
                       ! Obtain a point outside the block by mirroring the block
                       ! center Ijk_D to the starting location of this step IjkIni_D
                       IjkMid_D=2.*IjkIni_D-Ijk_D
                       ! Reduce length of Ijk_D --> IjkMid_D vector to end 
                       ! something like a 10th of a cell outside the block
                       dlp = 1.1*(1.-maxval(max(xmin-IjkMid_D,IjkMid_D-xmax) &
                            /(abs(IjkMid_D-IjkIni_D)+dlTiny)))
                       IjkMid_D=IjkIni_D+dlp*(IjkMid_D-IjkIni_D)

                       ! Make sure that IjkMid_D is just outside the control volume
                       IjkMid_D=max(xmin-.1,IjkMid_D)
                       IjkMid_D=min(xmax+.1,IjkMid_D)
                       call interpolate_xyz(IjkMid_D,XyzMid_D)
                       call interpolate_b(IjkMid_D, b_D, bNormMid_D)
                       IjkCur_D=IjkMid_D; XyzCur_D=XyzMid_D

                       ! We exited the block and have a good location to continued from
                       IsWall=.true.
                       EXIT HALF2
                    end if
                 else
                    !Step was OK, continue
                    EXIT HALF2
                 end if
              end do HALF2
           end if

           if(oktest_ray.and.okdebug) write(*,*) &
                'new decreased dl: me,iBlock,dl=',iProc,iBlock,dl
        else
           ! Too accurate, increase dl if possible
           if(abs(dl) < dlMax - dlTiny)then
              dlNext = sign(min(dlMax,abs(dl)/sqrt(dxRel)),dl)

              if(oktest_ray.and.okdebug) write(*,*) &
                   'new increased dlNext: me,iBlock,dlNext=',iProc,iBlock,dlNext
           end if
           EXIT STEP
        end if
     end do STEP

     ! Update position after the full step
     if(.not.IsWall)then
        if(UseOldMethodOfRayTrace .and. .not.UseCovariant)then
           IjkCur_D = IjkIni_D + bNormMid_D*dl
           XyzCur_D = XyzStart_BLK(:,iBlock) + Dxyz_D*(IjkCur_D - 1.)
        else
           XyzCur_D = XyzIni_D + dl*bNormMid_D
           call xyz_to_ijk(XyzCur_D,IjkCur_D,iBlock, &
                XyzIni_D,XyzStart_BLK(:,iBlock),Dxyz_D)

           ! Check if it stepped too far, use midpoint if it did
           if(any(IjkCur_D<(xmin-0.5)) .or. any(IjkCur_D>(xmax+0.5)))then
              IjkCur_D=IjkMid_D; XyzCur_D=XyzMid_D
           end if
        end if
     end if  !! .not.IsWall

     ! Update number of segments
     nSegment = nSegment + 1

     ! Step size in MH units  !!! Use simpler formula for cubic cells ???
     if(UseOldMethodOfRayTrace .and. .not.UseCovariant)then
        dLength = abs(dl)*sqrt( sum((bNormMid_D*Dxyz_D)**2) )
     else
        dLength = sqrt( (abs(XyzCur_D(1)-XyzIni_D(1)))**2 &
             +          (abs(XyzCur_D(2)-XyzIni_D(2)))**2 &
             +          (abs(XyzCur_D(3)-XyzIni_D(3)))**2 )
     end if

     ! Update ray length
     Length  = Length + dLength

     if(DoIntegrateRay)then

        ! Interpolate density and pressure
        ! Use the last indexes and distances already set in interpolate_b
        RhoP_V = &
             +dx1*(dy1*(dz1*Extra_VGB(:,i2,j2,k2,iBlock)   &
             +          dz2*Extra_VGB(:,i2,j2,k1,iBlock))  &
             +     dy2*(dz1*Extra_VGB(:,i2,j1,k2,iBlock)   &
             +          dz2*Extra_VGB(:,i2,j1,k1,iBlock))) &
             +dx2*(dy1*(dz1*Extra_VGB(:,i1,j2,k2,iBlock)   &
             +          dz2*Extra_VGB(:,i1,j2,k1,iBlock))  &
             +     dy2*(dz1*Extra_VGB(:,i1,j1,k2,iBlock)   &
             +          dz2*Extra_VGB(:,i1,j1,k1,iBlock)))

        ! Calculate physical step size divided by physical field strength
        InvBDl = dLength / sqrt( sum(b_D**2) )

        ! Intgrate field line volume = \int dl/B
        RayIntegral_V(InvB_) = RayIntegral_V(InvB_) + InvBDl

        ! Integrate density and pressure = \int Rho dl/B and \int P dl/B
        RayIntegral_V(RhoInvB_:pInvB_) = RayIntegral_V(RhoInvB_:pInvB_) + &
             InvBDl * RhoP_V

        ! Integrate density and pressure for multifluid
        if(DoMultiFluidIMCoupling) then
           RhoPMulti_V = dx1*(   dy1*(   dz1*ExtraMulti_VGB(:,i2,j2,k2,iBlock)   &
                +                         dz2*ExtraMulti_VGB(:,i2,j2,k1,iBlock))  &
                +                 dy2*(   dz1*ExtraMulti_VGB(:,i2,j1,k2,iBlock)   &
                +                         dz2*ExtraMulti_VGB(:,i2,j1,k1,iBlock))) &
                +         dx2*(   dy1*(   dz1*ExtraMulti_VGB(:,i1,j2,k2,iBlock)   &
                +                         dz2*ExtraMulti_VGB(:,i1,j2,k1,iBlock))  &
                +                 dy2*(   dz1*ExtraMulti_VGB(:,i1,j1,k2,iBlock)   &
                +                         dz2*ExtraMulti_VGB(:,i1,j1,k1,iBlock)))
           
            RayIntegral_V(HpRhoInvB_:OppInvB_) = &
                 RayIntegral_V(HpRhoInvB_:OppInvB_) + InvBDl * RhoPMulti_V
 
       endif

        ! Check if we crossed the Z=0 plane in the SM coord system
        ! Convert GM position into SM frame using the transposed GmSm_DD
        XyzSMIni_D = matmul(XyzIni_D, GmSm_DD)
        XyzSMCur_D = matmul(XyzCur_D, GmSm_DD)

        ! Check if we have crossed the magnetic equator in the SM frame
        if(XyzSMCur_D(3)*XyzSMIni_D(3)<=0)then

           ! Crossing the magnetic equator in opposite direction 
           ! is not accepted !!!

           if(DipoleStrength*(iRay-1.5)<0)then
              if(XyzSMIni_D(3) <= 0 .and. XyzSMCur_D(3) >= 0)then

                 ! This write is necessary to avoid incorrect 
                 ! optimization by the ifort 8.070 compiler
                 write(*,'(a)',ADVANCE='NO') ''

                 iFace = ray_loop_
                 EXIT FOLLOW
              end if
           else
              if(XyzSMIni_D(3) >= 0 .and. XyzSMCur_D(3) <= 0)then

                 ! This write is necessary to avoid incorrect
                 ! optimization by the ifort 8.070 compiler
                 write(*,'(a)',ADVANCE='NO') ''

                 iFace = ray_loop_
                 EXIT FOLLOW
              end if
           end if

           ! Interpolate x and y
           dz1 = abs(XyzSMIni_D(3))/(abs(XyzSMCur_D(3))+abs(XyzSMIni_D(3)))
           dz2 = 1.0 - dz1

           RayIntegral_V(Z0x_:Z0y_) = dz2*XyzSMIni_D(1:2) + dz1*XyzSMCur_D(1:2)

           ! Assign Z0b_ as the middle point value of the magnetic field
           RayIntegral_V(Z0b_) = sqrt(sum(b_D**2))
           if(oktest_ray)then
              write(*,'(a,3es12.4)') &
                   'Found z=0 crossing at XyzSMIni_D=',XyzSMIni_D
              write(*,'(a,3es12.4)') &
                   'Found z=0 crossing at XyzSMCur_D=',XyzSMCur_D
              write(*,'(a,3es12.4)')&
                   'RayIntegral_V(Z0x_:Z0b_)=',RayIntegral_V(Z0x_:Z0b_)

              !write(*,'(a,2es12.4)')'Weights =',dz1,dz2
              !write(*,'(a,3es12.4)')'b_D = ',b_D
           end if
        end if
     end if

     if(oktest_ray.and.okdebug)&
          write(*,*)'me,iBlock,nSegment,IjkCur_D=', &
          iProc,iBlock,nSegment,IjkCur_D

     if(DoCheckOpen)then
        if(all(ray(1,iRay,i1:i2,j1:j2,k1:k2,iBlock)==OPENRAY))then
           nOpen=nOpen+1
           iFace = ray_open_
           EXIT FOLLOW
        end if
     end if

     ! Check if we got inside the ionosphere
     if(DoCheckInnerBc)then
        r2Cur = sum(XyzCur_D**2)

        if(r2Cur<=R2_raytrace)then

           if(NameVectorField /= 'B')then
              XyzInOut_D = XyzCur_D
              iFace=ray_iono_
              EXIT FOLLOW
           end if

           ! Try mapping down to rIonosphere if we haven't tried yet (a lot)
           if(nIono<5)then
              if(.not.follow_ray_iono())then
                 ! We did not hit the surface of the ionosphere
                 ! continue the integration
                 nIono=nIono+1
              else
                 if(oktest_ray)write(*,'(a,3i4,6es12.4)')&
                      'Inside R_raytrace at me,iBlock,nSegment,IjkCur_D,XyzCur_D=',&
                      iProc,iBlock,nSegment,IjkCur_D,XyzCur_D

                 rCur=sqrt(r2Cur)
                 rIni=sqrt(sum(XyzIni_D**2))

                 ! The fraction of the last step inside body is estimated from 
                 ! the radii.
                 Fraction = (R_raytrace - rCur) / (rIni - rCur)

                 ! Reduce ray length
                 Length = Length - Fraction * dLength

                 ! Recalculate position
                 IjkCur_D = IjkCur_D - Fraction*(IjkCur_D-IjkIni_D)
                 call interpolate_xyz(IjkCur_D,XyzCur_D)

                 if(DoIntegrateRay)then
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
                    ! If Multifluid
                    if(DoMultiFluidIMCoupling)then
                       RayIntegral_V(HpRhoInvB_:OppInvB_) = &
                            RayIntegral_V(HpRhoInvB_:OppInvB_) - InvBDl * RhoPMulti_V
                    endif

                    if(oktest_ray)then
                       write(*,'(a,4es12.4)')&
                            'After  reduction InvBdl, RayIntegral_V=',InvBdl, &
                            RayIntegral_V(InvB_),RayIntegral_V(RhoInvB_:pInvB_)

                       write(*,*)'Reduction at InvBDl,RhoP_V   =',InvBDl,RhoP_V
                       write(*,*)'Reduction rIni,rCur,R_raytrace =',&
                            rIni,rCur,R_raytrace
                    end if

                 end if

                 ! Exit integration loop (XyzInOut_D is set by follow_ray_iono)
                 iFace=ray_iono_
                 EXIT FOLLOW
              end if
           end if
        end if
     end if

     ! Check if the ray hit the wall of the control volume
     if(any(IjkCur_D<xmin) .or. any(IjkCur_D>xmax))then
        ! Compute generalized coords without pole or edge wrapping
        call xyz_to_gen(XyzCur_D,Gen_D)

        if(any(Gen_D < XyzMin_D) .or. any(Gen_D > XyzMax_D))then
           iFace = ray_open_
        else
           iFace = ray_block_
        end if

        XyzInOut_D = XyzCur_D
        EXIT FOLLOW
     end if

     if(UseCovariant)then
        ! Can also hit wall if spherical before reaching xmin,xmax
        if(  XyzCur_D(1)<x1 .or. XyzCur_D(2)<y1 .or. XyzCur_D(3)<z1 .or. &
             XyzCur_D(1)>x2 .or. XyzCur_D(2)>y2 .or. XyzCur_D(3)>z2 )then

           XyzInOut_D = XyzCur_D
           iFace = ray_open_
           EXIT FOLLOW
        end if
     end if

     ! Check if we have integrated for too long
     if( nSegment > nSegmentMax .or. Length > RayLengthMax )then
        ! Seems to be a closed loop within a block
        if(oktest_ray) write(*,*)'CLOSED LOOP at me,iBlock,IjkCur_D,XyzCur_D=', &
             iProc,iBlock,IjkCur_D,XyzCur_D

        iFace=ray_loop_
        EXIT FOLLOW
     end if

  end do FOLLOW

  ! Extract last point if ray is done. 
  if(iFace /= ray_block_ .and. DoExtractRay)call ray_extract(IjkCur_D,XyzCur_D)

  if(oktest_ray) then
     write(*,'(a,4i4)')&
          'Finished follow_ray_block at me,iBlock,nSegment,iFace=',&
          iProc,iBlock,nSegment,iFace
     write(*,'(a,i4,9es12.4)')&
          'Finished follow_ray_block at me,IjkCur_D,XyzCur_D,XyzInOut_D=',&
          iProc,IjkCur_D,XyzCur_D,XyzInOut_D
  end if

contains
  !========================================================================

  subroutine interpolate_b(IjkIn_D,b_D,bNorm_D)

    ! Interpolate the magnetic field at normalized location IjkIn_D 
    ! and return the result in b_D. 
    ! The direction of b_D (normalized to a unit vector) is returned 
    ! in bNorm_D if the magnitude of b_D is not (almost) zero.

    real, intent(in)   :: IjkIn_D(3)  ! location
    real, intent(out)  :: b_D(3)      ! interpolated magnetic field
    real, intent(inout):: bNorm_D(3)  ! unit magnetic field vector

    !LOCAL VARIABLES:
    real :: AbsB, Dir0_D(3)

    !-------------------------------------------------------------------------

    ! Determine cell indices corresponding to location IjkIn_D
    i1=floor(IjkIn_D(1)); i2=i1+1
    j1=floor(IjkIn_D(2)); j2=j1+1
    k1=floor(IjkIn_D(3)); k2=k1+1

    ! Distance relative to the cell centers
    dx1 = IjkIn_D(1) - i1; dx2 = 1.0 - dx1
    dy1 = IjkIn_D(2) - j1; dy2 = 1.0 - dy1
    dz1 = IjkIn_D(3) - k1; dz2 = 1.0 - dz1

    ! Interpolate the magnetic field
    b_D = dx1*(dy1*(dz1*Bxyz_DGB(:,i2,j2,k2,iBlock)+dz2*Bxyz_DGB(:,i2,j2,k1,iBlock)) &
         +     dy2*(dz1*Bxyz_DGB(:,i2,j1,k2,iBlock)+dz2*Bxyz_DGB(:,i2,j1,k1,iBlock))) &
         +dx2*(dy1*(dz1*Bxyz_DGB(:,i1,j2,k2,iBlock)+dz2*Bxyz_DGB(:,i1,j2,k1,iBlock))  &
         +     dy2*(dz1*Bxyz_DGB(:,i1,j1,k2,iBlock)+dz2*Bxyz_DGB(:,i1,j1,k1,iBlock)))

    ! Set bNorm_D only if the magnetic field is not very small. 
    ! Otherwise continue in the previous direction.
    if(.not.UseOldMethodOfRayTrace .or. UseCovariant)then
       AbsB = sqrt(sum(b_D**2))
       if(AbsB > cTiny) bNorm_D = b_D/AbsB
       RETURN
    end if

    ! Stretch according to normalized coordinates
    Dir0_D = b_D/Dxyz_D
    AbsB = sqrt(sum(Dir0_D**2))
    if(AbsB > cTiny)bNorm_D = Dir0_D/AbsB

  end subroutine interpolate_b

  !========================================================================

  subroutine interpolate_xyz(IjkIn_D,XyzOut_D)

    ! Interpolate X/Y/Z at normalized location IjkIn_D 
    ! and return the result in XyzOut_D.

    real, intent(in)   :: IjkIn_D(3)  ! Ijk location
    real, intent(out)  :: XyzOut_D(3) ! Xyz location

    !-------------------------------------------------------------------------

    ! Determine cell indices corresponding to location IjkIn_D
    i1=floor(IjkIn_D(1)); i2=i1+1
    j1=floor(IjkIn_D(2)); j2=j1+1
    k1=floor(IjkIn_D(3)); k2=k1+1

    ! Distance relative to the cell centers
    dx1 = IjkIn_D(1) - i1; dx2 = 1.0 - dx1
    dy1 = IjkIn_D(2) - j1; dy2 = 1.0 - dy1
    dz1 = IjkIn_D(3) - k1; dz2 = 1.0 - dz1

    ! Interpolate the magnetic field
    XyzOut_D(1) = &
         +dx1*(dy1*(dz1*x_BLK(i2,j2,k2,iBlock)+dz2*x_BLK(i2,j2,k1,iBlock))  &
         +     dy2*(dz1*x_BLK(i2,j1,k2,iBlock)+dz2*x_BLK(i2,j1,k1,iBlock))) &
         +dx2*(dy1*(dz1*x_BLK(i1,j2,k2,iBlock)+dz2*x_BLK(i1,j2,k1,iBlock))  &
         +     dy2*(dz1*x_BLK(i1,j1,k2,iBlock)+dz2*x_BLK(i1,j1,k1,iBlock)))
    XyzOut_D(2) = &
         +dx1*(dy1*(dz1*y_BLK(i2,j2,k2,iBlock)+dz2*y_BLK(i2,j2,k1,iBlock))  &
         +     dy2*(dz1*y_BLK(i2,j1,k2,iBlock)+dz2*y_BLK(i2,j1,k1,iBlock))) &
         +dx2*(dy1*(dz1*y_BLK(i1,j2,k2,iBlock)+dz2*y_BLK(i1,j2,k1,iBlock))  &
         +     dy2*(dz1*y_BLK(i1,j1,k2,iBlock)+dz2*y_BLK(i1,j1,k1,iBlock)))
    XyzOut_D(3) = &
         +dx1*(dy1*(dz1*z_BLK(i2,j2,k2,iBlock)+dz2*z_BLK(i2,j2,k1,iBlock))  &
         +     dy2*(dz1*z_BLK(i2,j1,k2,iBlock)+dz2*z_BLK(i2,j1,k1,iBlock))) &
         +dx2*(dy1*(dz1*z_BLK(i1,j2,k2,iBlock)+dz2*z_BLK(i1,j2,k1,iBlock))  &
         +     dy2*(dz1*z_BLK(i1,j1,k2,iBlock)+dz2*z_BLK(i1,j1,k1,iBlock)))

  end subroutine interpolate_xyz

  !===========================================================================

  logical function follow_ray_iono()

    ! Follow ray inside ionosphere starting from XyzCur_D which is given in
    ! real coordinates and use analytic mapping.
    ! On return XyzInOut_D contains the final coordinates.
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

    call map_planet_field(Time_Simulation, XyzCur_D, TypeCoordSystem//' NORM',&
         rIonosphere, x_D, iHemisphere)

    if(iHemisphere==0)then
       write(*,*)'iHemisphere==0 for XyzCur_D=',XyzCur_D
       write(*,*)'iBlock, iRay=',iBlock,iRay
       call stop_mpi('ERROR in follow_ray_iono')
    end if

    if(iHemisphere*DipoleStrength*sign(1.0,1.5-iRay) < 0.0)then
       XyzInOut_D = x_D
       follow_ray_iono = .true.
    else
       follow_ray_iono = .false.
    end if

  end function follow_ray_iono

  !=========================================================================

  subroutine ray_extract(x_D,Xyz_D)

    use CON_line_extract, ONLY: line_put
    use ModPhysics, ONLY: No2Si_V, UnitX_, UnitRho_, UnitU_, UnitP_, UnitB_
    use ModAdvance, ONLY: State_VGB, nVar, &
         Rho_, RhoUx_, RhoUz_, Ux_, Uz_, p_, Bx_, Bz_
    use ModMain, ONLY: UseB0
    use ModRaytrace, ONLY: DoExtractBGradB1, bGradB1_DGB
    use ModInterpolate, ONLY: trilinear

    real, intent(in) :: x_D(3),Xyz_D(3)

    real    :: State_V(nVar), B0_D(3), PlotVar_V(4+nVar+3)
    integer :: n, iLine
    !----------------------------------------------------------------------

    PlotVar_V(1)   = Length
    PlotVar_V(2:4) = Xyz_D

    if(DoExtractUnitSi) PlotVar_V(1:4) = PlotVar_V(1:4)*No2Si_V(UnitX_)

    if(DoExtractState)then

       ! Determine cell indices corresponding to location x_D
       i1=floor(x_D(1)); i2=i1+1
       j1=floor(x_D(2)); j2=j1+1
       k1=floor(x_D(3)); k2=k1+1

       ! Distance relative to the cell centers
       dx1 = x_D(1) - i1; dx2 = 1.0 - dx1
       dy1 = x_D(2) - j1; dy2 = 1.0 - dy1
       dz1 = x_D(3) - k1; dz2 = 1.0 - dz1

       ! Interpolate state to x_D
       State_V = &
            +dx1*(dy1*(dz1*State_VGB(:,i2,j2,k2,iBlock)   &
            +          dz2*State_VGB(:,i2,j2,k1,iBlock))  &
            +     dy2*(dz1*State_VGB(:,i2,j1,k2,iBlock)   &
            +          dz2*State_VGB(:,i2,j1,k1,iBlock))) &
            +dx2*(dy1*(dz1*State_VGB(:,i1,j2,k2,iBlock)   &
            +          dz2*State_VGB(:,i1,j2,k1,iBlock))  &
            +     dy2*(dz1*State_VGB(:,i1,j1,k2,iBlock)   &
            +          dz2*State_VGB(:,i1,j1,k1,iBlock)))

       ! Convert momentum to velocity
       State_V(Ux_:Uz_) = State_V(RhoUx_:RhoUz_)/State_V(Rho_)
       if(DoMultiFluidIMCoupling)then
          State_V(iUx_I(IonFirst_:iIonSecond)) = &
               State_V(iRhoUx_I(IonFirst_:iIonSecond))/ &
               State_V(iRho_I(IonFirst_:iIonSecond))
          State_V(iUy_I(IonFirst_:iIonSecond)) = &
               State_V(iRhoUy_I(IonFirst_:iIonSecond))/ &
               State_V(iRho_I(IonFirst_:iIonSecond))
          State_V(iUz_I(IonFirst_:iIonSecond)) = &
               State_V(iRhoUz_I(IonFirst_:iIonSecond))/ &
               State_V(iRho_I(IonFirst_:iIonSecond))
       end if

       ! Add B0 to the magnetic field
       if(UseB0)then
          call get_b0(Xyz_D(1),Xyz_D(2),Xyz_D(3),B0_D)
          State_V(Bx_:Bz_) = State_V(Bx_:Bz_) + B0_D
       end if

       ! Convert to SI units if required
       if(DoExtractUnitSi)then
          State_V(Rho_)    = State_V(Rho_)    * No2Si_V(UnitRho_)
          State_V(Ux_:Uz_) = State_V(Ux_:Uz_) * No2Si_V(UnitU_)
          State_V(Bx_:Bz_) = State_V(Bx_:Bz_) * No2Si_V(UnitB_)
          State_V(p_)      = State_V(p_)      * No2Si_V(UnitP_)
          if( DoMultiFluidIMCoupling)then
             State_V(iRho_I(IonFirst_:iIonSecond)) = &
                  State_V(iRho_I(IonFirst_:iIonSecond)) * No2Si_V(UnitRho_)
             State_V(iUx_I(IonFirst_:iIonSecond)) = &
                  State_V(iUx_I(IonFirst_:iIonSecond)) * No2Si_V(UnitU_)
             State_V(iUy_I(IonFirst_:iIonSecond)) = &
                  State_V(iUy_I(IonFirst_:iIonSecond)) * No2Si_V(UnitU_)
             State_V(iUz_I(IonFirst_:iIonSecond)) = &
                  State_V(iUz_I(IonFirst_:iIonSecond)) * No2Si_V(UnitU_)
             State_V(iP_I(IonFirst_:iIonSecond))      =  &
                  State_V(iP_I(IonFirst_:iIonSecond))  * No2Si_V(UnitP_)
          end if
       end if

       PlotVar_V(5:4+nVar) = State_V

       n = 4 + nVar

       if(DoExtractBGradB1)then

          n = n + 3

          ! Interpolate b.grad B1 into the last 3 elements
          PlotVar_V(n-2:n) = &
               trilinear(bGradB1_DGB(:,:,:,:,iBlock), &
               3, 0, nI+1, 0, nJ+1, 0, nK+1, x_D, DoExtrapolate=.false.)

          if(DoExtractUnitSi) PlotVar_V(n-2:n) = &
               PlotVar_V(n-2:n) * No2Si_V(UnitB_)/No2Si_V(UnitX_)

       end if

    else
       n = 4
    end if

    ! get a unique line index based on starting indexes
    ! ignore index 4 if nRay_D(4) is zero
    iLine = &
         ((max(0, min(nRay_D(4), iStart_D(4)-1))*nRay_D(3) &
         + max(0,iStart_D(3)-1) )*nRay_D(2) &
         + max(0,iStart_D(2)-1) )*nRay_D(1) &
         + iStart_D(1)

    if(iLine < 0)then
       write(*,*)'iLine=',iLine
       write(*,*)'nRay_D  =',nRay_D
       write(*,*)'iStart_D=',iStart_D
       call stop_mpi('DEBUG')
    end if
    call line_put(iLine,n,PlotVar_V(1:n))

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
  use ModRayTrace, ONLY: CpuTimeStartRay
  use ModMpi, ONLY: MPI_WTIME

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


  ! Assign face ray values to cell centers

  !nOpen = 0
  CpuTimeStartRay = MPI_WTIME()
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

subroutine integrate_ray_accurate(nLat, nLon, Lat_I, Lon_I, Radius, NameVar)

  use CON_ray_trace, ONLY: ray_init
  use CON_planet_field, ONLY: map_planet_field
  use CON_axes, ONLY: transform_matrix
  use ModRaytrace
  use ModMain,    ONLY: nBlock, Time_Simulation, TypeCoordSystem, UseB0, DoMultiFluidIMCoupling
  use ModPhysics, ONLY: rBody
  use ModAdvance, ONLY: nVar, State_VGB, Rho_, p_, Bx_, Bz_, B0_DGB
  use ModProcMH
  use ModMpi
  use ModMessagePass,    ONLY: message_pass_dir
  use ModNumConst,       ONLY: cDegToRad, cTiny
  use ModCoordTransform, ONLY: sph_to_xyz
  use ModUtilities,      ONLY: check_allocate
  use ModGeometry,       ONLY: XyzMax_D, XyzMin_D, x1, x2, y1, y2, z1, z2, &
       UseCovariant
  use CON_line_extract,  ONLY: line_init, line_collect
  use CON_planet,        ONLY: DipoleStrength
  use ModMultiFluid
  implicit none

  !INPUT ARGUMENTS:
  integer, intent(in):: nLat, nLon
  real,    intent(in):: Lat_I(nLat), Lon_I(nLon), Radius
  character(len=*), intent(in):: NameVar


  !DESCRIPTION:
  ! Lat_I(nLat) and Lon_I(nLon) are the coordinates of a 2D spherical 
  ! grid in the SM(G) coordinate system in degrees. The 2D grid is 
  ! at radius Radius given in units of planet radii.
  ! NameVar lists the variables that need to be extracted and/or integrated.
  ! The subroutine can calculate the integral of various quantities 
  ! and/or extract state variables along the field lines starting from the 2D 
  ! spherical grid.

  real    :: Theta, Phi, Lat, Lon, XyzIono_D(3), Xyz_D(3)
  integer :: iLat, iLon, iHemisphere, iRay
  integer :: iProcFound, iBlockFound, i, j, k
  integer :: nStateVar, iIonSecond
  integer :: iError
  logical :: DoTest, DoTestMe
  character(len=*), parameter :: NameSub = 'integrate_ray_accurate'
  !-------------------------------------------------------------------------

  iIonSecond = min(IonFirst_+1, IonLast_)

  call set_oktest(NameSub, DoTest, DoTestMe)

  if(DoTest)write(*,*)NameSub,' starting on iProc=',iProc,&
       ' with nLat, nLon, Radius=',nLat,nLon,Radius

  iLatTest = 23; iLonTest = 3

  call timing_start('integrate_ray')

  oktest_ray = .false.

  ! Initialize some basic variables
  R_raytrace      = rBody
  R2_raytrace     = R_raytrace**2
  if(UseOldMethodOfRayTrace .and. .not.UseCovariant)then
     RayLengthMax = 4.*sum(XyzMax_D - XyzMin_D)
  else
     RayLengthMax = 4.*(abs(x2-x1) + abs(y2-y1) + abs(z2-z1))
  end if

  DoIntegrateRay = index(NameVar, 'InvB') > 0 .or. index(NameVar, 'Z0') > 0
  DoExtractRay   = index(NameVar, '_I') > 0
  DoTraceRay     = .false.
  DoMapRay       = .false.

  if(DoTestMe)write(*,*)NameSub,' DoIntegrateRay,DoExtractRay,DoTraceRay=',&
       DoIntegrateRay,DoExtractRay,DoTraceRay

  if(DoExtractRay)then
     nRay_D  = (/ nLat, nLon, 0, 0 /)
     DoExtractState = .true.
     DoExtractUnitSi= .true.
     nStateVar = 4 + nVar
     call line_init(nStateVar)
  end if

  NameVectorField = 'B'

  ! (Re)initialize CON_ray_trace
  call ray_init(iComm)

  ! Copy magnetic field into Bxyz_DGB
  Bxyz_DGB(:,:,:,:,1:nBlock) = State_VGB(Bx_:Bz_,:,:,:,1:nBlock)

  ! Fill in all ghost cells
  call message_pass_dir(iDirMin=1, iDirMax=3, Width=2, SendCorners=.true., &
       ProlongOrder=2, nVar=3, Sol_VGB=Bxyz_DGB)

  ! Add B0 for faster interpolation
  if(UseB0) Bxyz_DGB(1:3,:,:,:,1:nBlock) = &
       Bxyz_DGB(1:3,:,:,:,1:nBlock) + B0_DGB(:,:,:,:,1:nBlock)

  if(DoIntegrateRay)then
     ! Copy density and pressure into Extra_VGB
     Extra_VGB(1,:,:,:,1:nBlock) = State_VGB(rho_,:,:,:,1:nBlock)
     Extra_VGB(2,:,:,:,1:nBlock) = State_VGB(p_  ,:,:,:,1:nBlock)

     ! Fill in all ghost cells (faces+edges+corners) without monotone restrict
     call message_pass_dir(iDirMin=1, iDirMax=3, Width=2, SendCorners=.true., &
          ProlongOrder=2, nVar=2, Sol_VGB=Extra_VGB)

     allocate(&
          RayIntegral_VII(nRayIntegral,nLat,nLon), &
          RayResult_VII(nRayIntegral,nLat,nLon), STAT=iError)
     call check_allocate(iError,NameSub//' RayIntegral_VII,RayResult_VII')
     RayIntegral_VII = 0.0
     RayResult_VII   = 0.0

     if(DoMultiFluidIMCoupling) then 
        ! Copy density and pressure into Extra_VGB
        ExtraMulti_VGB(1:2,:,:,:,1:nBlock) = &
             State_VGB(iRho_I(IonFirst_:iIonSecond),:,:,:,1:nBlock)
        ExtraMulti_VGB(3:4,:,:,:,1:nBlock) = &
             State_VGB(iP_I(IonFirst_:iIonSecond),:,:,:,1:nBlock)
        
        ! Fill in all ghost cells (faces+edges+co5Brners) without monotone restrict
        call message_pass_dir(iDirMin=1, iDirMax=3, Width=2, SendCorners=.true., &
          ProlongOrder=2, nVar=4, Sol_VGB=ExtraMulti_VGB)
     endif
  end if

  ! Transformation matrix between the SM and GM coordinates
  GmSm_DD = transform_matrix(time_simulation,'SMG',TypeCoordSystem)

  ! Integrate rays starting from the latitude-longitude pairs defined
  ! by the arrays Lat_I, Lon_I
  CpuTimeStartRay = MPI_WTIME()
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

        ! Figure out direction of tracing outward
        if(iHemisphere*DipoleStrength>0)then
           iRay = 1
        else
           iRay = 2
        end if

        ! Check if the mapping is on the north hemisphere
        if(iHemisphere == 0)then
           !  write(*,*)NameSub,' point did not map to rBody, ',&
           !   'implement analytic integrals here! Lat, Lon=', Lat, Lon
           CYCLE
        end if

        ! Convert SM position to GM (Note: these are identical for ideal axes)
        Xyz_D = matmul(GmSm_DD,Xyz_D)

        ! Find processor and block for the location
        call xyz_to_peblk(Xyz_D(1), Xyz_D(2), Xyz_D(3), &
             iProcFound, iBlockFound, .true., i, j, k)

        ! If location is on this PE, follow and integrate ray
        if(iProc == iProcFound)then

           if(DoTest .and. iLat==iLatTest .and. iLon==iLonTest)then
              write(*,'(a,2i3,a,i3,a,i4)') &
                   'start of ray iLat, iLon=',iLat, iLon,&
                   ' found on iProc=',iProc,' iBlock=',iBlockFound
              write(*,'(a,2i4,2es12.4)')'iLon, iLat, Lon, Lat=',&
                   iLon, iLat, Lon, Lat
              write(*,'(a,3es12.4)')'XyzIono_D=',XyzIono_D
              write(*,'(a,3es12.4)')'Xyz_D    =',Xyz_D
           end if

           call follow_ray(iRay, (/iLat, iLon, 0, iBlockFound/), Xyz_D)

        end if
     end do
  end do

  ! Do remaining rays obtained from other PE-s
  call finish_ray

  if(DoTest.and.iLatTest<=nLat.and.iLonTest<=nLon) &
       write(*,*)NameSub,' iProc, RayIntegral_VII=',&
       iProc, RayIntegral_VII(:,iLatTest,iLonTest)

  if(DoIntegrateRay) call MPI_reduce( &
       RayIntegral_VII, RayResult_VII, nLat*nLon*nRayIntegral, &
       MPI_REAL, MPI_SUM, 0, iComm, iError)

  if(DoExtractRay) call line_collect(iComm,0)

  call timing_stop('integrate_ray')

end subroutine integrate_ray_accurate

!============================================================================

subroutine integrate_ray_accurate_1d(nPts, XyzPt_DI, NameVar)

  use CON_ray_trace,     ONLY: ray_init
  use CON_planet_field,  ONLY: map_planet_field
  use CON_axes,          ONLY: transform_matrix
  use CON_line_extract,  ONLY: line_init, line_collect
  use CON_planet,        ONLY: DipoleStrength
  use ModRaytrace
  use ModMain,           ONLY: nBlock, Time_Simulation, TypeCoordSystem, &
       UseB0, DoMultiFluidIMCoupling
  use ModPhysics,        ONLY: rBody
  use ModAdvance,        ONLY: nVar, State_VGB, Rho_, p_, Bx_, Bz_, B0_DGB
  use ModProcMH
  use ModMpi
  use ModMessagePass,    ONLY: message_pass_dir
  use ModNumConst,       ONLY: cDegToRad, cTiny
  use ModCoordTransform, ONLY: sph_to_xyz
  use ModUtilities,      ONLY: check_allocate
  use ModGeometry,       ONLY: XyzMax_D, XyzMin_D, x1, x2, y1, y2, z1, z2, &
       UseCovariant
  use ModMultiFluid
  implicit none

  !INPUT ARGUMENTS:
  integer, intent(in):: nPts
  real,    intent(in):: XyzPt_DI(3,nPts)
  character(len=*), intent(in):: NameVar

  !DESCRIPTION:
  ! A 1D list of points is sent in with x,y,z values in GM coordinates.
  ! NameVar lists the variables that need to be extracted and/or integrated.
  ! The subroutine can calculate the integral of various quantities 
  ! and/or extract state variables along the field lines starting from the
  ! points sent in.

  real    :: Xyz_D(3)
  integer :: iPt, iRay
  integer :: iProcFound, iBlockFound, i, j, k
  integer :: nStateVar, iIonSecond
  integer :: iError
  logical :: DoTest, DoTestMe
  character(len=*), parameter :: NameSub = 'integrate_ray_accurate_1d'
  !-------------------------------------------------------------------------

  iIonSecond = min(IonFirst_+1, IonLast_)

  call set_oktest(NameSub, DoTest, DoTestMe)

  if(DoTest)write(*,*)NameSub,' starting on iProc=',iProc,' with nPts=',nPts

  call timing_start('integrate_ray_1d')

  oktest_ray = .false.

  ! Initialize some basic variables
  R_raytrace      = rBody
  R2_raytrace     = R_raytrace**2
  if(UseOldMethodOfRayTrace .and. .not.UseCovariant)then
     RayLengthMax = 4.*sum(XyzMax_D - XyzMin_D)
  else
     RayLengthMax = 4.*(abs(x2-x1) + abs(y2-y1) + abs(z2-z1))
  end if

  DoIntegrateRay = index(NameVar, 'InvB') > 0 .or. index(NameVar, 'Z0') > 0
  DoExtractRay   = index(NameVar, '_I') > 0
  DoTraceRay     = .false.
  DoMapRay       = .false.

  if(DoTestMe)write(*,*)NameSub,' DoIntegrateRay,DoExtractRay,DoTraceRay=',&
       DoIntegrateRay,DoExtractRay,DoTraceRay

  if(DoExtractRay)then
     nRay_D  = (/ 2, nPts, 0, 0 /)
     DoExtractState = .true.
     DoExtractUnitSi= .true.
     nStateVar = 4 + nVar
     call line_init(nStateVar)
  end if

  NameVectorField = 'B'

  ! (Re)initialize CON_ray_trace
  call ray_init(iComm)

  ! Copy magnetic field into Bxyz_DGB
  Bxyz_DGB(:,:,:,:,1:nBlock) = State_VGB(Bx_:Bz_,:,:,:,1:nBlock)

  ! Fill in all ghost cells
  call message_pass_dir(iDirMin=1, iDirMax=3, Width=2, SendCorners=.true., &
       ProlongOrder=2, nVar=3, Sol_VGB=Bxyz_DGB)

  ! Add B0 for faster interpolation
  if(UseB0) Bxyz_DGB(1:3,:,:,:,1:nBlock) = &
       Bxyz_DGB(1:3,:,:,:,1:nBlock) + B0_DGB(:,:,:,:,1:nBlock)

  if(DoIntegrateRay)then
     ! Copy density and pressure into Extra_VGB
     Extra_VGB(1,:,:,:,1:nBlock) = State_VGB(rho_,:,:,:,1:nBlock)
     Extra_VGB(2,:,:,:,1:nBlock) = State_VGB(p_  ,:,:,:,1:nBlock)

     ! Fill in all ghost cells (faces+edges+corners) without monotone restrict
     call message_pass_dir(iDirMin=1, iDirMax=3, Width=2, SendCorners=.true., &
          ProlongOrder=2, nVar=2, Sol_VGB=Extra_VGB)

     allocate(&
          RayIntegral_VII(nRayIntegral,nRay_D(1),nRay_D(2)), &
          RayResult_VII(nRayIntegral,nRay_D(1),nRay_D(2)), STAT=iError)
     call check_allocate(iError,NameSub//' RayIntegral_VII,RayResult_VII')
     RayIntegral_VII = 0.0
     RayResult_VII   = 0.0

     if(DoMultiFluidIMCoupling) then 
        ! Copy density and pressure into Extra_VGB
        ExtraMulti_VGB(1:2,:,:,:,1:nBlock) = &
             State_VGB(iRho_I(IonFirst_:iIonSecond),:,:,:,1:nBlock)
        ExtraMulti_VGB(3:4,:,:,:,1:nBlock) = &
             State_VGB(iP_I(IonFirst_:iIonSecond),:,:,:,1:nBlock)
        
        ! Fill in all ghost cells (faces+edges+co5Brners) without monotone restrict
        call message_pass_dir(iDirMin=1, iDirMax=3, Width=2, SendCorners=.true., &
          ProlongOrder=2, nVar=4, Sol_VGB=ExtraMulti_VGB)
     endif
  end if

  ! Transformation matrix between the SM and GM coordinates
  GmSm_DD = transform_matrix(time_simulation,'SMG',TypeCoordSystem)

  ! Integrate rays
  CpuTimeStartRay = MPI_WTIME()
  do iPt = 1, nPts
     Xyz_D=XyzPt_DI(:,iPt)

     ! Find processor and block for the location
     call xyz_to_peblk(Xyz_D(1), Xyz_D(2), Xyz_D(3), &
          iProcFound, iBlockFound, .true., i, j, k)

     ! If location is on this PE, follow and integrate ray
     if(iProc == iProcFound)then
        call follow_ray(1, (/1, iPt, 0, iBlockFound/), Xyz_D)
        call follow_ray(2, (/2, iPt, 0, iBlockFound/), Xyz_D)
     end if
  end do

  ! Do remaining rays obtained from other PE-s
  call finish_ray

  if(DoIntegrateRay) call MPI_reduce( RayIntegral_VII, RayResult_VII, &
       size(RayIntegral_VII), MPI_REAL, MPI_SUM, 0, iComm, iError)

  if(DoExtractRay) call line_collect(iComm,0)

  call timing_stop('integrate_ray_1d')

end subroutine integrate_ray_accurate_1d

!============================================================================

subroutine plot_ray_equator(iFile)

  use ModMain, ONLY: n_step, time_accurate, Time_Simulation, TypeCoordSystem
  use ModIo,   ONLY: StringDateOrTime, NamePlotDir, plot_range
  use ModAdvance, ONLY: nVar, Ux_, Uz_, Bx_, Bz_
  use ModProcMH,  ONLY: iProc
  use ModNumConst,       ONLY: cTwoPi
  use ModIoUnit,         ONLY: UnitTmp_
  use ModPlotFile,       ONLY: save_plot_file
  use ModRayTrace,       ONLY: RayMap_DSII
  use CON_line_extract,  ONLY: line_get, line_clean
  use CON_axes,          ONLY: transform_matrix

  implicit none

  !INPUT ARGUMENTS:
  integer, intent(in):: iFile

  !DESCRIPTION:
  ! Follow field lines starting from a 2D polar grid on the 
  ! magnetic equatorial plane in the SM(G) coordinate system.
  ! The grid parameters are given by plot_rang(1:4, iFile)
  ! The subroutine extracts coordinates and state variables
  ! along the field lines going in both directions 
  ! starting from the 2D polar grid.

  integer :: nRadius, nLon
  real    :: rMin, rMax
  integer :: iR, iLon
  integer :: iPoint, nPoint, nVarOut
  real, allocatable :: Radius_I(:), Longitude_I(:), PlotVar_VI(:,:)
  real    :: SmGm_DD(3,3)
  real    :: PlotVar_V(0:4+nVar)

  character(len=100) :: NameFile, NameFileEnd=''

  logical :: DoTest, DoTestMe
  character(len=*), parameter :: NameSub = 'plot_ray_equator'
  !-------------------------------------------------------------------------

  call set_oktest(NameSub, DoTest, DoTestMe)

  ! Extract grid info from plot_range (see MH_set_parameters for plot_type eqr)
  nRadius = plot_range(1, iFile)
  nLon    = plot_range(2, iFile)
  rMin    = plot_range(3, iFile)
  rMax    = plot_range(4, iFile)

  allocate(Radius_I(nRadius), Longitude_I(nLon))
  do iR = 1, nRadius
     Radius_I(iR) = rMin + (iR-1)*(rMax - rMin)/(nRadius - 1)
  end do
  do iLon = 1, nLon
     Longitude_I(iLon) = (iLon-1)*cTwoPi / (nLon - 1)
  end do

  call trace_ray_equator(nRadius, nLon, Radius_I, Longitude_I, .false.)

  deallocate(Radius_I, Longitude_I)

  if(iProc/=0) RETURN

  if(time_accurate) NameFileEnd = "_t"//StringDateOrTime
  write(NameFileEnd,'(a,i7.7,a)') trim(NameFileEnd) // '_n',n_step, '.out'

  ! Transformation matrix between the SM and GM coordinates
  SmGm_DD = transform_matrix(time_simulation,TypeCoordSystem,'SMG')

  ! Set the file name

  call line_get(nVarOut, nPoint)
  if(nVarOut /= 4+nVar)call stop_mpi(NameSub//': nVarOut error')
  allocate(PlotVar_VI(0:nVarOut, nPoint))
  call line_get(nVarOut, nPoint, PlotVar_VI, DoSort=.true.)

  NameFile = trim(NamePlotDir)//"eqr"//NameFileEnd

  open(UnitTmp_, FILE=NameFile, STATUS="replace")
  write(UnitTmp_, *) 'nRadius, nLon, nPoint=',nRadius, nLon, nPoint
  write(UnitTmp_, *) 'iLine l x y z rho ux uy uz bx by bz p'

  do iPoint = 1, nPoint
     ! Convert vectors to SM coordinates
     PlotVar_V = PlotVar_VI(:, iPoint)
     PlotVar_V(2:4) = matmul(SmGm_DD,PlotVar_V(2:4))
     PlotVar_V(4+Ux_:4+Uz_) = matmul(SmGm_DD,PlotVar_V(4+Ux_:4+Uz_))
     PlotVar_V(4+Bx_:4+Bz_) = matmul(SmGm_DD,PlotVar_V(4+Bx_:4+Bz_))
     ! Save into file
     write(UnitTmp_, *) PlotVar_V
  end do
  close(UnitTmp_)
  deallocate(PlotVar_VI)
  call line_clean

  ! Now save the mapping files
  NameFile = trim(NamePlotDir)//"map_north_"//NameFileEnd
  call save_plot_file( &
       NameFile, &
       StringHeaderIn = 'Mapping to northern ionosphere', &
       TimeIn       = time_simulation, &
       nStepIn      = n_step, &
       NameVarIn    = 'r Lon rIono ThetaIono PhiIono', &
       CoordMinIn_D = (/rMin,   0.0/), &
       CoordMaxIn_D = (/rMax, 360.0/), &
       VarIn_VII  = RayMap_DSII(:,1,:,:))

  NameFile = trim(NamePlotDir)//"map_south_"//NameFileEnd
  call save_plot_file( &
       NameFile, &
       StringHeaderIn = 'Mapping to southern ionosphere', &
       TimeIn       = time_simulation, &
       nStepIn      = n_step, &
       NameVarIn    = 'r Lon rIono ThetaIono PhiIono', &
       CoordMinIn_D = (/rMin,   0.0/), &
       CoordMaxIn_D = (/rMax, 360.0/), &
       VarIn_VII  = RayMap_DSII(:,2,:,:))

  deallocate(RayMap_DSII)

end subroutine plot_ray_equator

!============================================================================

subroutine trace_ray_equator(nRadius, nLon, Radius_I, Longitude_I, &
     DoMessagePass)

  use ModMain, ONLY: x_, y_, z_, nI, nJ, nK, UnusedBlk
  use CON_ray_trace, ONLY: ray_init
  use CON_axes, ONLY: transform_matrix
  use ModRaytrace, ONLY: oktest_ray, R_raytrace, R2_raytrace, RayLengthMax, &
       DoIntegrateRay, DoExtractRay, DoTraceRay, DoMapRay, &
       DoExtractState, DoExtractUnitSi, DoExtractBGradB1, bGradB1_DGB, &
       RayMap_DSII, RayMapLocal_DSII, &
       NameVectorField, Bxyz_DGB, nRay_D, CpuTimeStartRay, GmSm_DD, CLOSEDRAY
  use ModMain,    ONLY: nBlock, Time_Simulation, TypeCoordSystem, UseB0
  use ModPhysics, ONLY: rBody
  use ModAdvance, ONLY: nVar, State_VGB, Bx_, Bz_, B0_DGB
  use ModProcMH,  ONLY: iProc, iComm
  use ModMpi
  use ModGeometry,       ONLY: x1, x2, y1, y2, z1, z2, dx_BLK, dy_BLK, dz_BLK
  use CON_line_extract,  ONLY: line_init, line_collect, line_clean
  use ModMessagePass,    ONLY: message_pass_dir
  use ModCoordTransform, ONLY: xyz_to_sph

  implicit none

  !INPUT ARGUMENTS:
  integer, intent(in):: nRadius, nLon
  real,    intent(in):: Radius_I(nRadius), Longitude_I(nLon)
  logical, intent(in):: DoMessagePass

  !DESCRIPTION:
  ! Follow field lines starting from a 2D polar grid on the 
  ! magnetic equatorial plane in the SM(G) coordinate system.
  ! The grid parameters are given by the arguments.
  ! The subroutine extracts coordinates and state variables
  ! along the field lines going in both directions 
  ! starting from the 2D equatorial grid.
  ! Fill in ghost cells if DoMessagePass is true.

  integer :: iR, iLon, iSide
  integer :: iProcFound, iBlockFound, iBlock, i, j, k, iError
  real    :: r, Phi, XyzSm_D(3), Xyz_D(3), b_D(3)

  integer :: nStateVar

  logical :: DoTest, DoTestMe
  character(len=*), parameter :: NameSub = 'trace_ray_equator'
  !-------------------------------------------------------------------------

  call set_oktest(NameSub, DoTest, DoTestMe)

  ! Extract grid info from plot_range (see MH_set_parameters for plot_type eqr)
  if(DoTest)then
     write(*,*)NameSub,' starting on iProc=',iProc,&
          ' with nRadius, nLon=', nRadius, nLon
     write(*,*)NameSub,' Radius_I   =',Radius_I
     write(*,*)NameSub,' Longitude_I=',Longitude_I
  end if

  call timing_start(NameSub)

  ! Fill in all ghost cells
  if(DoMessagePass)  call message_pass_dir(iDirMin=1, iDirMax=3, Width=2, &
       SendCorners=.true., ProlongOrder=2, nVar=nVar, Sol_VGB=State_VGB)

  oktest_ray = .false.

  ! Initialize some basic variables
  R_raytrace   = rBody
  R2_raytrace  = R_raytrace**2
  RayLengthMax = 2.*(abs(x2-x1) + abs(y2-y1) + abs(z2-z1))

  DoIntegrateRay = .false.
  DoExtractRay   = .true.
  DoTraceRay     = .false.
  DoMapRay       = .true.

  if(DoExtractBGradB1)then
     allocate(bGradB1_DGB(3,0:nI+1,0:nJ+1,0:nK+1,nBlock))
     do iBlock = 1, nBlock
        if(UnusedBlk(iBlock)) CYCLE
        do k = 0, nK+1; do j = 0, nJ+1; do i = 0, nI+1; 
           b_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
           if(UseB0) b_D = b_D +  B0_DGB(:,i,j,k,iBlock)
           b_D = b_D/sqrt(max(1e-30,sum(b_D**2)))
        
           bGradB1_DGB(1,i,j,k,iBlock) = 0.5*sum(b_D *  &
                ( State_VGB(Bx_:Bz_,i+1,j,k,iBlock)     &
                - State_VGB(Bx_:Bz_,i-1,j,k,iBlock))) / &
                dx_BLK(iBlock)

           bGradB1_DGB(2,i,j,k,iBlock) = 0.5*sum(b_D *  &
                ( State_VGB(Bx_:Bz_,i,j+1,k,iBlock)     & 
                - State_VGB(Bx_:Bz_,i,j-1,k,iBlock))) / &
                dy_BLK(iBlock)

           bGradB1_DGB(3,i,j,k,iBlock) = 0.5*sum(b_D * &
                ( State_VGB(Bx_:Bz_,i,j,k+1,iBlock) &
                - State_VGB(Bx_:Bz_,i,j,k+1,iBlock))) / &
                dz_BLK(iBlock)

        end do; end do; end do
     end do
  end if

  nRay_D  = (/ 2, nRadius, nLon, 0 /)
  DoExtractState = .true.
  DoExtractUnitSi= .true.
  nStateVar = 4 + nVar
  if(DoExtractBGradB1) nStateVar = nStateVar + 3
  call line_init(nStateVar)

  NameVectorField = 'B'

  ! (Re)initialize CON_ray_trace
  call ray_init(iComm)

  ! Copy magnetic field into Bxyz_DGB
  Bxyz_DGB(:,:,:,:,1:nBlock) = State_VGB(Bx_:Bz_,:,:,:,1:nBlock)

  ! Add B0 for faster interpolation
  if(UseB0) Bxyz_DGB(1:3,:,:,:,1:nBlock) = &
       Bxyz_DGB(1:3,:,:,:,1:nBlock) + B0_DGB(:,:,:,:,1:nBlock)

  ! Transformation matrix between the SM and GM coordinates
  GmSm_DD = transform_matrix(time_simulation,'SMG',TypeCoordSystem)

  ! Integrate rays starting from the latitude-longitude pairs defined
  ! by the arrays Lat_I, Lon_I
  CpuTimeStartRay = MPI_WTIME()
  do iR = 1, nRadius

     r = Radius_I(iR)

     if(r < rBody*1.0001) CYCLE

     do iLon = 1, nLon

        Phi = Longitude_I(iLon)

        ! Convert polar coordinates to Cartesian coordinates in SM
        XyzSm_D(x_) = r*cos(Phi)
        XyzSm_D(y_) = r*sin(Phi)
        XyzSm_D(z_) = 0.0
        
        ! Convert SM position to GM (Note: these are identical for ideal axes)
        Xyz_D = matmul(GmSm_DD,XyzSm_D)

        ! Find processor and block for the location
        call xyz_to_peblk(Xyz_D(1), Xyz_D(2), Xyz_D(3), &
             iProcFound, iBlockFound, .true., i, j, k)

        ! If location is on this PE, follow and integrate ray
        if(iProc == iProcFound)then

           call follow_ray(1, (/1, iR, iLon, iBlockFound/), Xyz_D)
           call follow_ray(2, (/2, iR, iLon, iBlockFound/), Xyz_D)

        end if
     end do
  end do

  ! Do remaining rays obtained from other PE-s
  call finish_ray

  ! Collect all rays onto processor 0
  call line_collect(iComm,0)

  ! Clean data except on processor 0
  if(iProc /= 0)call line_clean

  ! Some procs never have their RayMap arrays allocated.
  if(.not.allocated(RayMap_DSII)) then
     allocate(RayMap_DSII(3,nRay_D(1),nRay_D(2),nRay_D(3)))
     allocate(RayMapLocal_DSII(3,nRay_D(1),nRay_D(2),nRay_D(3)))
     RayMapLocal_DSII = 0.0
     RayMap_DSII = 0.0
  end if

  ! Collect the ray mapping info to processor 0  
  call MPI_reduce(RayMapLocal_DSII, RayMap_DSII, size(RayMap_DSII), MPI_REAL, &
       MPI_SUM, 0, iComm, iError)
  deallocate(RayMapLocal_DSII)

  if(iProc == 0)then
     do iLon = 1, nLon; do iR = 1, nRadius; do iSide = 1, 2
        if(RayMap_DSII(1,iSide,iR,iLon) < CLOSEDRAY) CYCLE
        Xyz_D   = RayMap_DSII(:,iSide,iR,iLon)
        XyzSm_D = matmul(Xyz_D,GmSm_DD)
        call xyz_to_sph(XyzSm_D, RayMap_DSII(:,iSide,iR,iLon))
     end do; end do; end do
  else
     deallocate(RayMap_DSII)
  end if

  if(DoExtractBGradB1) deallocate(bGradB1_DGB)

  if(DoMessagePass)call exchange_messages

  call timing_stop(NameSub)

end subroutine trace_ray_equator

!============================================================================

subroutine test_ray_integral

  use ModRayTrace, ONLY: RayResult_VII, InvB_, xEnd_, zEnd_, &
       nRayIntegral, xyz_to_latlon, iLatTest, iLonTest
  use ModProcMH,   ONLY: iProc
  use ModIoUnit,   ONLY: UNITTMP_
  use ModNumConst, ONLY: cTiny
  use ModMain,     ONLY: DoMultiFluidIMCoupling
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
  call integrate_ray_accurate(nLat,nLon,Lat_I,Lon_I,1.0, &
       'InvB,RhoInvB,pInvB,Z0x,Z0y,Z0b')

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

           if(iLat == iLatTest .and. iLon == iLonTest)then
              if(.not. DoMultiFluidIMCoupling)then
                 write(*,'(a,a)')'iLon iLat Lon Lat ',&
                      'Bvol Z0x Z0y Z0b Rho P LatEnd LonEnd Zend Length'
              else
                 write(*,'(a,a)')'iLon iLat Lon Lat ',&
                      'Bvol Z0x Z0y Z0b Rho P LatEnd LonEnd Zend Length HpRho OpRho HpP OpP'
              endif
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

subroutine ray_lines(nLine, IsParallel_I, Xyz_DI)

  ! Extract nLine ray lines parallel or anti_parallel according to
  ! IsParallel_I(nLine), starting from positions Xyz_DI(3,nLine).
  ! The results are stored by CON_line_extract.

  use ModProcMH,   ONLY: iProc, iComm
  use ModRayTrace, ONLY: DoTraceRay, DoMapRay, DoIntegrateRay, DoExtractRay, &
       CpuTimeStartRay, oktest_ray, &
       nRay_D, NameVectorField, R_Raytrace, R2_Raytrace, RayLengthMax, Bxyz_DGB
  use CON_ray_trace, ONLY: ray_init
  use ModAdvance,  ONLY: State_VGB, RhoUx_, RhoUy_, RhoUz_, Bx_, By_, Bz_, B0_DGB
  use ModMain,     ONLY: nI, nJ, nK, nBlock, unusedBLK, UseB0
  use ModPhysics,  ONLY: rBody
  use ModGeometry, ONLY: XyzMax_D, XyzMin_D, Dx_BLK, Dy_BLK, Dz_BLK
  use ModMpi,      ONLY: MPI_WTIME

  implicit none

  !INPUT ARGUMENTS:
  integer, intent(in) :: nLine
  logical, intent(in) :: IsParallel_I(nLine)
  real,    intent(in) :: Xyz_DI(3, nLine)

  !EOP
  real    :: Xyz_D(3), Dx2Inv, Dy2Inv, Dz2Inv
  integer :: iProcFound, iBlockFound, iLine, iRay

  integer :: i, j, k, iBlock

  character(len=*), parameter :: NameSub = 'ray_lines'
  logical :: DoTest, DoTestMe
  !-------------------------------------------------------------------------
  call set_oktest(NameSub, DoTest, DoTestMe)

  ! Initialize R_raytrace, R2_raytrace
  oktest_ray = .false.
  R_raytrace   = rBody
  R2_raytrace  = R_raytrace**2
  RayLengthMax = 2*sum(XyzMax_D - XyzMin_D)

  DoTraceRay     = .false.
  DoMapRay       = .false.
  DoIntegrateRay = .false.
  DoExtractRay   = .true.
  nRay_D = (/ nLine, 0, 0, 0 /)

  ! (Re)initialize CON_ray_trace
  call ray_init(iComm)

  select case(NameVectorField)
  case('B')
     ! Store B1+B0 for faster interpolation
     if(UseB0)then
        Bxyz_DGB(1:3,:,:,:,1:nBlock) = State_VGB(Bx_:Bz_,:,:,:,1:nBlock) &
             + B0_DGB(:,:,:,:,1:nBlock)
     else
        Bxyz_DGB(1:3,:,:,:,1:nBlock) = State_VGB(Bx_:Bz_,:,:,:,1:nBlock)
     end if
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
     call stop_mpi(NameSub//': invalid NameVectorField='//NameVectorField)
  end select

  ! Start extracting rays
  CpuTimeStartRay = MPI_WTIME()
  do iLine = 1, nLine
     Xyz_D = Xyz_DI(:,iLine)

     call xyz_to_peblk(Xyz_D(1), Xyz_D(2), Xyz_D(3), &
          iProcFound, iBlockFound, .true., i, j, k)

     if(iProc == iProcFound)then
        if(DoTest)write(*,*)NameSub,' follows ray ',iLine,&
             ' from iProc,iBlock,i,j,k=',iProcFound, iBlockFound, i, j, k
        if(IsParallel_I(iLine))then
           iRay = 1
        else
           iRay = 2
        end if
        call follow_ray(iRay, (/iLine, 0, 0, iBlockFound/), Xyz_D)
     end if
  end do

  ! Do remaining rays obtained from other PE-s
  call finish_ray

end subroutine ray_lines

!==============================================================================

subroutine write_plot_line(iFile)

  use ModProcMH,   ONLY: iComm, iProc
  use ModRayTrace, ONLY: NameVectorField, DoExtractState, DoExtractUnitSi
  use ModVarIndexes,ONLY: nVar, NamePrimitiveVar, NamePrimitiveVarTec
  use ModIO,       ONLY: StringDateOrTime, &
       NamePlotDir, plot_type, plot_form, plot_dimensional, Plot_, &
       NameLine_I, nLine_I, XyzStartLine_DII, IsParallelLine_II, IsSingleLine_I
  use ModMain,     ONLY: n_step, time_accurate, time_simulation
  use ModIoUnit,   ONLY: UnitTmp_
  use CON_line_extract, ONLY: line_init, line_collect, line_get, line_clean

  implicit none

  integer, intent(in) :: iFile ! The file index of the plot file

  character(len=100) :: NameFile, NameStart, NameVar, StringTitle
  integer            :: nLineFile, nStateVar, nPlotVar
  integer            :: iPoint, nPoint, iPointNext, nPoint1

  real, pointer :: PlotVar_VI(:,:)

  integer :: iPlotFile, iLine, nLine, nVarOut

  logical :: IsSingleLine, IsIdl

  character(len=*), parameter :: NameSub = 'write_plot_line'
  logical :: DoTest, DoTestMe
  !-------------------------------------------------------------------------
  call set_oktest(NameSub, DoTest, DoTestMe)

  ! Set the global ModRaytrace variables for this plot file
  iPlotFile      = iFile - Plot_
  select case(NameLine_I(iPlotFile))
  case('A', 'B')
     NameVectorField = 'B'
  case('U','J')
     NameVectorField = NameLine_I(iPlotFile)
  case default
     write(*,*) NameSub//' WARNING invalid NameVectorField='// &
          NameVectorField//' for iPlotFile=',iPlotFile
     RETURN
  end select
  DoExtractState = index(plot_type(iFile),'pos')<1
  DoExtractUnitSi= plot_dimensional(iFile)

  ! Set the number lines and variables to be extracted
  nLine     = nLine_I(iPlotFile)
  nStateVar = 4
  if(DoExtractState) nStateVar = nStateVar + nVar

  ! Initialize CON_line_extract
  call line_init(nStateVar)

  ! Obtain the line data
  call ray_lines(nLine, IsParallelLine_II(1:nLine,iPlotFile), &
       XyzStartLine_DII(:,1:nLine,iPlotFile))

  ! Collect lines from all PE-s to Proc 0
  call line_collect(iComm,0)

  if(iProc==0)then
     call line_get(nVarOut, nPoint)
     if(nVarOut /= nStateVar)call stop_mpi(NameSub//': nVarOut error')
     allocate(PlotVar_VI(0:nVarOut, nPoint))
     call line_get(nVarOut, nPoint, PlotVar_VI, DoSort=.true.)
  end if

  call line_clean

  ! Only iProc 0 works on writing the plot files
  if(iProc /= 0) RETURN

  ! Write the result into 1 or more plot files from processor 0

  IsSingleLine = IsSingleLine_I(iPlotFile)

  if(IsSingleLine)then
     nLineFile = nLine
  else
     nLineFile = 1
  end if

  if(iPlotFile < 10)then
     write(NameStart,'(a,i1,a)') &
          trim(NamePlotDir)//trim(plot_type(iFile))//'_',iPlotFile
  else
     write(NameStart,'(a,i2,a)') &
          trim(NamePlotDir)//trim(plot_type(iFile))//'_',iPlotFile
  end if
  NameStart = trim(NameStart)//'_'//NameLine_I(iPlotFile)

  if(time_accurate)call get_time_string

  ! Set the title
  if(IsSingleLine)then
     StringTitle = NameVectorField//' line'
  else
     StringTitle = NameVectorField//' lines'
  end if

  ! Add the string describing the units
  if(DoExtractUnitSi)then
     StringTitle = trim(StringTitle)//" in SI units"
  else
     StringTitle = trim(StringTitle)//" in normalized units"
  end if

  ! The Length is used as a coordinate in the IDL file, so it is not a plot var
  nPlotVar = nStateVar - 1
  ! Add 1 for the Index array if it is needed in the plot file
  if(.not. IsSingleLine)nPlotVar = nPlotVar + 1

  ! Set the name of the variables
  select case(plot_form(iFile))
  case('idl')
     IsIdl = .true.
     NameVar = 'Length x y z'
     if(DoExtractState)NameVar = trim(NameVar)//' '//NamePrimitiveVar
     if(IsSingleLine)then
        NameVar = trim(NameVar)//' iLine'
     else
        NameVar = trim(NameVar)//' Index nLine'
     end if
  case('tec')
     IsIdl = .false.
     NameVar = '"X", "Y", "Z"'
     if(DoExtractState)NameVar = trim(NameVar)//' ,'//NamePrimitiveVarTec
     if(.not.IsSingleLine)NameVar = trim(NameVar)//', "Index"'
     NameVar = trim(NameVar)//', "Length"'
  case default
     call CON_stop(NameSub//' ERROR invalid plot form='//plot_form(iFile))
  end select

  ! Write out plot files
  ! If IsSingleLine is true write a new file for every line,
  ! otherwise write a single file

  iPointNext = 1
  do iLine = 1, nLineFile

     ! Set the file name
     NameFile = NameStart
     if(IsSingleLine .and. nLine > 1)then
        if(nLine < 10)then
           write(NameFile,'(a,i1)') trim(NameFile),iLine
        else
           write(NameFile,'(a,i2)') trim(NameFile),iLine
        end if
     end if
     if(time_accurate) NameFile = trim(NameFile)// "_t"//StringDateOrTime
     write(NameFile,'(a,i7.7,a)') trim(NameFile) // '_n',n_step

     if(IsIdl)then
        NameFile = trim(NameFile) // '.out'
     else
        NameFile = trim(NameFile) // '.dat'
     end if

     ! Figure out the number of points for this ray
     if(IsSingleLine) nPoint1 = count(nint(PlotVar_VI(0,1:nPoint))==iLine)

     open(UnitTmp_,file=NameFile)
     if(IsIdl)then
        write(UnitTmp_,'(a79)') trim(StringTitle)//'_var11'
        write(UnitTmp_,'(i7,1pe13.5,3i3)') &
             n_step,time_simulation,1,1,nPlotVar
        if(IsSingleLine)then
           write(UnitTmp_,'(i6)') nPoint1
           write(UnitTmp_,'(es13.5)') real(iLine)
        else
           write(UnitTmp_,'(i6)') nPoint
           write(UnitTmp_,'(es13.5)') real(nLine)
        end if
        write(UnitTmp_,'(a79)') NameVar
     else
        write(UnitTmp_,'(a)')'TITLE ="'//trim(StringTitle)//'"'
        write(UnitTmp_,'(a)')'VARIABLES='//trim(NameVar)
        if(IsSingleLine)then
           write(UnitTmp_,'(a,i2.2,a,i6)')'ZONE T="'// &
                NameVectorField//' line ',iLine,'", '//'I=',nPoint1
        else
           write(UnitTmp_,'(a,i2.2,a,i6)')'ZONE T="'// &
                NameVectorField//' ',nLine,' lines", '//'I=',nPoint
        end if
     end if

     ! Write out data
     if(IsSingleLine)then
        ! Write out the part corresponding to this line
        do iPoint = iPointNext, iPointNext + nPoint1 - 1
           if(IsIdl)then
              ! Write Length as the first variable: the 1D coordinate
              write(UnitTmp_,'(50es18.10)') PlotVar_VI(1:nStateVar,iPoint)
           else
              ! Write Length as the last variable, so that 
              ! x,y,z can be used as 3D coordinates
              write(UnitTmp_,'(50es18.10)') PlotVar_VI(2:nStateVar,iPoint),&
                   PlotVar_VI(1,iPoint)
           end if
        end do
        iPointNext = iPointNext + nPoint1
     else
        do iPoint = 1, nPoint
           if(IsIdl)then
              ! Write Index as the last variable
              write(UnitTmp_, '(50es18.10)') &
                   PlotVar_VI(1:nStateVar, iPoint), PlotVar_VI(0,iPoint)
           else
              ! Write Index and Length as the last 2 variables
              write(UnitTmp_, '(50es18.10)') &
                   PlotVar_VI(2:nStateVar, iPoint), PlotVar_VI(0:1,iPoint)
           end if
        end do
     end if
     close(UnitTmp_)
  end do

  deallocate(PlotVar_VI)

end subroutine write_plot_line

!============================================================================

subroutine xyz_to_ijk(XyzIn_D,IjkOut_D,iBlock,XyzRef_D,GenRef_D,dGen_D)
  use ModNumConst,  ONLY: cPi, cHalfPi, cTwoPi
  use ModCovariant, ONLY: TypeGeometry
  use ModMain,      ONLY: Phi_,Theta_,x_,y_,z_, nJ
  implicit none

  integer          ,intent(in)  :: iBlock
  real,dimension(3),intent(in)  :: XyzIn_D,XyzRef_D,GenRef_D,dGen_D
  real,dimension(3),intent(out) :: IjkOut_D
  real,dimension(3)             :: Gen_D

  select case(TypeGeometry)
  case('cartesian')
     Gen_D=XyzIn_D

  case('spherical','spherical_lnr','spherical_genr')   
     call xyz_to_gen(XyzIn_D,Gen_D)

     ! Did I cross the pole?
     if( (XyzIn_D(x_)*XyzRef_D(x_) + XyzIn_D(y_)*XyzRef_D(y_)) < 0.)then
        Gen_D(Phi_)=Gen_D(Phi_) + GenRef_D(Phi_) &
             - modulo((cPi+GenRef_D(Phi_)),cTwoPi)
        if(XyzIn_D(z_)>0.)then
           Gen_D(Theta_)=Gen_D(Theta_)+2.*(+cHalfPi-Gen_D(Theta_))
        else
           Gen_D(Theta_)=Gen_D(Theta_)+2.*(-cHalfPi-Gen_D(Theta_))
        end if
     end if

     ! Did I cross periodic boundary?
     if    ((+Gen_D(Phi_)-GenRef_D(Phi_))/dGen_D(Phi_) > 2*nJ)then
        Gen_D(Phi_)=Gen_D(Phi_)-cTwoPi
     elseif((-Gen_D(Phi_)+GenRef_D(Phi_))/dGen_D(Phi_) > 1*nJ)then 
        Gen_D(Phi_)=Gen_D(Phi_)+cTwoPi
     end if

  case('axial_torus')
     call stop_mpi('xyz_to_gen_ray2: cannot handle TypeGeometry='//TypeGeometry)

  case default
     call stop_mpi('xyz_to_gen_ray2: unknown TypeGeometry='//TypeGeometry)

  end select

  ! Gen_D is set, now compute Ijk
  IjkOut_D = (Gen_D - GenRef_D)/dGen_D + 1.

end subroutine xyz_to_ijk

!============================================================================

subroutine lcb_plot(iFile)
  use CON_line_extract,  ONLY: line_get, line_clean
  use CON_planet_field,  ONLY: map_planet_field
  use CON_axes,          ONLY: transform_matrix
  use ModIoUnit,         ONLY: UnitTmp_
  use ModAdvance,        ONLY: nVar, Ux_, Uz_, Bx_, Bz_
  use ModMain,           ONLY: Time_Simulation, TypeCoordSystem, time_accurate, n_step
  use ModNumConst,       ONLY: cRadToDeg, cDegToRad
  use ModProcMH,         ONLY: iProc, iComm
  use ModPhysics,        ONLY: Si2No_V, No2Si_V, UnitX_, UnitRho_, UnitP_, UnitB_, rBody
  use ModCoordTransform, ONLY: sph_to_xyz, xyz_to_sph
  use ModIO,             ONLY: StringDateOrTime, NamePlotDir, plot_range, plot_type
  use ModRaytrace,       ONLY: RayResult_VII, RayIntegral_VII, InvB_,RhoInvB_,pInvB_
  use ModMpi
  implicit none

  integer, intent(in) :: iFile

  character (len=*), parameter :: NameSub='lcb_plot'
  character (len=80) :: FileName,stmp
  integer, parameter :: nPts=11, nD=6
  integer :: i,j,k, nLine, iStart,iMid,iEnd, jStart,jMid,jEnd, iLon, nLon, iD, iLC
  integer :: iPoint, nPoint, nVarOut, iHemisphere, nFile, iError, nTP, iDirZ
  real :: PlotVar_V(0:4+nVar)
  real :: Radius, RadiusIono, Lon, zL,zU, zUs=40., xV,yV, Integrals(3)
  real :: XyzIono_D(3), Xyz_D(3)
  real :: Gsm2Smg_DD(3,3) = reshape( (/ 1.,0.,0.,  0.,1.,0.,  0.,0.,1. /), (/3,3/) )
  real :: Smg2Gsm_DD(3,3) = reshape( (/ 1.,0.,0.,  0.,1.,0.,  0.,0.,1. /), (/3,3/) )
  real, allocatable :: PlotVar_VI(:,:), XyzPt_DI(:,:), zPt_I(:)
  logical :: Map1,Map2, Odd, Skip, SaveIntegrals
  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest, DoTestMe)
  if(DoTest)write(*,*)NameSub,': starting'

  ! Extract grid info from plot_range (see MH_set_parameters for plot_type lcb)
  Radius = plot_range(1, iFile)
  nLon   = plot_range(2, iFile)

  SaveIntegrals=.false.
  if(index(plot_type(iFile),'int')>0) SaveIntegrals=.true.

  !Use a value of 1. for these plots.
  RadiusIono = 1.
  nTP=int( (rBody-RadiusIono)/.1 )

  if(.not.allocated(XyzPt_DI)) allocate(XyzPt_DI(3,nPts), zPt_I(nPts))

  ! Transformation matrix from default (GM) to SM coordinates
  Gsm2Smg_DD = transform_matrix(time_simulation,TypeCoordSystem,'SMG')
  Smg2Gsm_DD = transform_matrix(time_simulation,'SMG','GSM')

  if(iProc == 0)then
     FileName=trim(NamePlotDir)//'LCB-GM'
     if(iFile <  10) write(FileName, '(a,i1)') trim(FileName)//"_",iFile
     if(iFile >= 10) write(FileName, '(a,i2)') trim(FileName)//"_",iFile
     if(time_accurate)then
        call get_time_string
        FileName = trim(FileName) // "_t" // StringDateOrTime
     end if
     write(FileName,'(a,i7.7,a)') trim(FileName)//"_n",n_step,".dat"

     open( UnitTmp_, FILE=trim(FileName), STATUS="replace")
     write(UnitTmp_,'(a)')'TITLE="IE B traces (GM Coordinates)"'
     if(SaveIntegrals)then
        write(UnitTmp_,'(a)')'VARIABLES="X [R]", "Y [R]", "Z [R]", "1/B", "n/B", "p/B"'
     else
        write(UnitTmp_,'(a)')'VARIABLES="X [R]", "Y [R]", "Z [R]"'
     end if
  end if

  do iDirZ = -1,1,2
     !compute the last closed points on cylinder for positive and negative Z values

     do iLon=1,nLon
        Lon = (360./nLon)*(iLon-1)
        xV = Radius*cos(cDegToRad*Lon)
        yV = Radius*sin(cDegToRad*Lon)

        zL=0.;  zU=zUs*iDirZ
        Skip=.false.
        do iD=1,nD
           if(Skip) CYCLE

           iLC=-9

           ! Create in SM coords
           XyzPt_DI(1,:) = xV
           XyzPt_DI(2,:) = yV
           do i=1,nPts
              XyzPt_DI(3,i) = zL + ((i-1)*((zU-zL)/(nPts-1)))
           end do
           zPt_I = XyzPt_DI(3,:)

           ! Convert to GM coords
           do i=1,nPts
              XyzPt_DI(:,i) = matmul(Smg2Gsm_DD,XyzPt_DI(:,i))
           end do

           call integrate_ray_accurate_1d(nPts, XyzPt_DI, 'InvB,RhoInvB,pInvB,extract_I')

           if(iProc == 0)then

              Integrals = -1.
              
              call line_get(nVarOut, nPoint)
              if(nPoint>0)then
                 !PlotVar_VI variables = 'iLine l x y z rho ux uy uz bx by bz p'
                 allocate(PlotVar_VI(0:nVarOut, nPoint))
                 call line_get(nVarOut, nPoint, PlotVar_VI, DoSort=.true.)

                 k=0
                 do iPoint = 1, nPoint
                    nLine=PlotVar_VI(0,iPoint)
                    if(k == nLine) CYCLE
                    Odd=.true.;  if( (nLine/2)*2 == nLine )Odd=.false.

                    !\\
                    ! finish previous line
                    if(k/=0)then
                       if(Odd)then
                          iEnd = iPoint-1
                          Map2 = .false.
                          Xyz_D=PlotVar_VI(2:4,iEnd) * Si2No_V(UnitX_)
                          if(sqrt(Xyz_D(1)**2 + Xyz_D(2)**2 + Xyz_D(3)**2)<1.5*rBody)Map2 = .true.
                             
                          if(Map1 .and. Map2)then
                             iLC=k/2
                             jStart=iStart; jMid=iMid; jEnd=iEnd
                             Integrals(1) = sum(RayResult_VII(   InvB_,:,iLC)) / No2Si_V(UnitB_)
                             Integrals(2) = sum(RayResult_VII(RhoInvB_,:,iLC)) * No2Si_V(UnitRho_)
                             Integrals(3) = sum(RayResult_VII(  pInvB_,:,iLC)) * No2Si_V(UnitP_)
                          end if
                       else
                          iEnd = iPoint-1
                          Map1 = .false.
                          Xyz_D=PlotVar_VI(2:4,iEnd) * Si2No_V(UnitX_)
                          if(sqrt(Xyz_D(1)**2 + Xyz_D(2)**2 + Xyz_D(3)**2)<1.5*rBody)Map1 = .true.
                       end if
                    end if

                    !\\
                    ! start new line counters
                    k=nLine
                    if(Odd)then
                       iStart = iPoint
                    else
                       iMid = iPoint
                    end if
                 end do

                 !\\
                 ! finish last line
                 if(k/=0)then
                    iEnd = nPoint
                    Map2 = .false.
                    Xyz_D=PlotVar_VI(2:4,iEnd) * Si2No_V(UnitX_)
                    if(sqrt(Xyz_D(1)**2 + Xyz_D(2)**2 + Xyz_D(3)**2)<1.5*rBody)Map2 = .true.

                    if(Map1 .and. Map2)then
                       iLC=k/2
                       jStart=iStart; jMid=iMid; jEnd=iEnd
                       Integrals(1) = sum(RayResult_VII(   InvB_,:,iLC))
                       Integrals(2) = sum(RayResult_VII(RhoInvB_,:,iLC))
                       Integrals(3) = sum(RayResult_VII(  pInvB_,:,iLC))
                    end if
                 end if

                 !\\
                 ! write only last closed
                 if(iD == nD .and. iLC.ne.-9)then
                    j=(jEnd-jStart)+2*nTP
                    write(UnitTmp_,'(a,f7.2,a,a,i8,a)') 'ZONE T="LCB lon=',Lon,'"', &
                         ', I=',j,', J=1, K=1, ZONETYPE=ORDERED, DATAPACKING=POINT'
                    Xyz_D = PlotVar_VI(2:4,jMid-1) * Si2No_V(UnitX_)
                    do i=0,nTP-1
                       ! Map from the ionosphere to first point
                       call map_planet_field(time_simulation, Xyz_D, 'GSM NORM', &
                            RadiusIono+i*.1, XyzIono_D, iHemisphere)
                       if(SaveIntegrals)then
                          write(UnitTmp_, *) XyzIono_D,Integrals
                       else
                          write(UnitTmp_, *) XyzIono_D
                       end if
                    end do
                    do i=jMid-1,jStart+1,-1
                       PlotVar_V = PlotVar_VI(:, i)
                       Xyz_D = PlotVar_V(2:4) * Si2No_V(UnitX_)
                       if(SaveIntegrals)then
                          write(UnitTmp_, *) Xyz_D,Integrals
                       else
                          write(UnitTmp_, *) Xyz_D
                       end if
                    end do
                    do i=jMid,jEnd
                       PlotVar_V = PlotVar_VI(:, i)
                       Xyz_D = PlotVar_V(2:4) * Si2No_V(UnitX_)
                       if(SaveIntegrals)then
                          write(UnitTmp_, *) Xyz_D,Integrals
                       else
                          write(UnitTmp_, *) Xyz_D
                       end if
                    end do
                    do i=nTP-1,0,-1
                       ! Map from last point to the ionosphere
                       call map_planet_field(time_simulation, Xyz_D, 'GSM NORM', &
                            RadiusIono+i*.1, XyzIono_D, iHemisphere)
                       if(SaveIntegrals)then
                          write(UnitTmp_, *) XyzIono_D,Integrals
                       else
                          write(UnitTmp_, *) XyzIono_D
                       end if
                    end do
                 end if

                 deallocate(PlotVar_VI)
              end if
              call line_clean
           end if  !iProc==0

           if(allocated(RayIntegral_VII)) deallocate(RayIntegral_VII)
           if(allocated(RayResult_VII))   deallocate(RayResult_VII)

           ! set new zL and zU
           call MPI_Bcast(iLC,1,MPI_INTEGER,0,iComm,iError)
           if(iLC == -9)then
              Skip=.true.
           elseif(iLC == nPts)then
              zL=   zUs*iDirZ
              zU=2.*zUs*iDirZ
           else
              zL = zPt_I(iLC)
              zU = zPt_I(iLC+1)
           end if

        end do  !iD loop

     end do  !iLon loop

  end do  !iDirZ loop

  if(iProc == 0) close(UnitTmp_)

  if(DoTest)write(*,*)NameSub,': finished'
end subroutine lcb_plot

!**************************************************************

subroutine ieb_plot(iFile)
  use CON_line_extract,  ONLY: line_get, line_clean
  use CON_planet_field,  ONLY: map_planet_field
  use CON_axes,          ONLY: transform_matrix
  use ModIoUnit,         ONLY: UnitTmp_
  use ModAdvance,        ONLY: nVar, Ux_, Uz_, Bx_, Bz_
  use ModMain,           ONLY: Time_Simulation, TypeCoordSystem, time_accurate, n_step
  use ModNumConst,       ONLY: cRadToDeg, cDegToRad
  use ModProcMH,         ONLY: iProc, iComm
  use ModPhysics,        ONLY: Si2No_V, No2Si_V, UnitX_, rBody
  use ModCoordTransform, ONLY: sph_to_xyz, xyz_to_sph
  use ModIO,             ONLY: StringDateOrTime, NamePlotDir
  use ModRaytrace,       ONLY: RayResult_VII, RayIntegral_VII, rIonosphere
  implicit none

  integer, intent(in) :: iFile

  character (len=*), parameter :: NameSub='ieb_plot'
  character (len=80) :: FileName,stmp
  character (len=2) :: Coord
  character (len=1) :: NS
  integer :: i,j,k, nLat,nLon, nLine, nTP, iStart,iEnd, iLat,iLon, OC
  integer :: iPoint, nPoint, nVarOut, iHemisphere, nFile
  real :: PlotVar_V(0:4+nVar)
  real :: Radius, Lat,Lon, Theta,Phi, LonOC
  real :: XyzIono_D(3), RtpIono_D(3), Xyz_D(3)
  real :: Gsm2Smg_DD(3,3) = reshape( (/ 1.,0.,0.,  0.,1.,0.,  0.,0.,1. /), (/3,3/) )
  real :: Smg2Gsm_DD(3,3) = reshape( (/ 1.,0.,0.,  0.,1.,0.,  0.,0.,1. /), (/3,3/) )
  real, allocatable :: PlotVar_VI(:,:), IE_lat(:), IE_lon(:)
  logical :: MapDown
  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------
  call CON_set_do_test(NameSub,DoTest, DoTestMe)
  if(DoTest)write(*,*)NameSub,': starting'

  nLat=181
  nLon=36
  if(.not.allocated(IE_lat)) allocate(IE_lat(nLat), IE_lon(nLon))

  ! Load grid and convert to lat-lon in degrees
  do i=1,nLat
     IE_lat(i) = 90.-1.*(i-1)
  end do
  do i=1,nLon
     IE_lon(i) = 10.*(i-1)
  end do
  !  Radius = (6378.+100.)/6378.
  Radius = rIonosphere
  nTP=int( (rBody-Radius)/.1 )

  call integrate_ray_accurate(nLat, nLon, IE_lat, IE_lon, Radius, 'extract_I')

  if(iProc == 0)then

     ! Transformation matrix from default (GM) to SM coordinates
     Gsm2Smg_DD = transform_matrix(time_simulation,TypeCoordSystem,'SMG')
     Smg2Gsm_DD = transform_matrix(time_simulation,'SMG','GSM')

     call line_get(nVarOut, nPoint)
     if(nPoint>0)then
        !PlotVar_VI variables = 'iLine l x y z rho ux uy uz bx by bz p'
        allocate(PlotVar_VI(0:nVarOut, nPoint))
        call line_get(nVarOut, nPoint, PlotVar_VI, DoSort=.true.)

        do nFile=1,4

           if(nFile==1)then
              Coord = 'SM';  NS = 'N'
           elseif(nFile==2)then
              Coord = 'GM';  NS = 'N'
           elseif(nFile==3)then
              Coord = 'SM';  NS = 'S'
           elseif(nFile==4)then
              Coord = 'GM';  NS = 'S'
           end if
           FileName=trim(NamePlotDir)//'IEB-'//trim(Coord)//'-'//trim(NS)
           if(time_accurate)then
              call get_time_string
              FileName = trim(FileName) // "_t" // StringDateOrTime
           end if
           write(FileName,'(a,i7.7,a)') trim(FileName)//"_n", n_step,".dat"

           open( UnitTmp_, FILE=trim(FileName), STATUS="replace")
           if(Coord == 'GM')then
              write(UnitTmp_,'(a)')'TITLE="IE B traces (GM Coordinates)"'
           else
              write(UnitTmp_,'(a)')'TITLE="IE B traces (SM Coordinates)"'
           end if
           write(UnitTmp_,'(a)')'VARIABLES="X [R]", "Y [R]", "Z [R]", "Lat", "Lon", "OC"'

           k=0
           LonOC=-1.
           do iPoint = 1, nPoint
              nLine=PlotVar_VI(0,iPoint)
              if(k /= nLine)then
                 !\\
                 ! finish previous line
                 if(k/=0)then
                    iEnd = iPoint-1
                    MapDown = .false.
                    Xyz_D=PlotVar_VI(2:4,iEnd) * Si2No_V(UnitX_)
                    if(sqrt(Xyz_D(1)**2 + Xyz_D(2)**2 + Xyz_D(3)**2)<1.5*rBody)MapDown = .true.
                    j=(1+iEnd-iStart)+(nTP+1)
                    if(MapDown)j=j+(nTP+1)
                    OC=-1; if(MapDown)OC=2
                    if(MapDown .and. LonOC/=Lon)OC=1

!\
!                     write(UnitTmp_,'(a,2f7.2,a,a,i8,a)') 'ZONE T="IEB ll=',Lat,Lon,'"', &
!                          ', I=',j,', J=1, K=1, ZONETYPE=ORDERED, DATAPACKING=POINT'
!-
                    write(UnitTmp_,'(a,2f7.2,a,a,f7.2,a,i8,a)') 'ZONE T="IEB ll=',Lat,Lon,'"', &
                         ', STRANDID=1, SOLUTIONTIME=',Lon, &
                         ', I=',j,', J=1, K=1, ZONETYPE=ORDERED, DATAPACKING=POINT'
!/
                    write(stmp,'(f8.2)')Lat
                    write(UnitTmp_,'(a,a,a)') 'AUXDATA LAT="',trim(adjustl(stmp)),'"'
                    write(stmp,'(f8.2)')Lon
                    write(UnitTmp_,'(a,a,a)') 'AUXDATA LON="',trim(adjustl(stmp)),'"'

                    ! Convert to SMG Cartesian coordinates on the surface of the ionosphere
                    Theta = cDegToRad*(90.0 - Lat)     
                    Phi = cDegToRad*Lon
                    call sph_to_xyz(Radius, Theta, Phi, XyzIono_D)
                    Xyz_D=XyzIono_D
                    if(Coord == 'GM') Xyz_D = matmul(Smg2Gsm_DD,Xyz_D)
                    write(UnitTmp_, *) Xyz_D,Lat,Lon,OC
                    do i=1,nTP
                       ! Map from the ionosphere to rBody
                       call map_planet_field(time_simulation, XyzIono_D, 'SMG NORM', &
                            Radius+i*.1, Xyz_D, iHemisphere)
                       if(Coord == 'GM') Xyz_D = matmul(Smg2Gsm_DD,Xyz_D)
                       write(UnitTmp_, *) Xyz_D,Lat,Lon,OC
                    end do
                    do i=iStart,iEnd
                       ! Convert vectors to SM coordinates
                       PlotVar_V = PlotVar_VI(:, i)
                       PlotVar_V(2:4) = matmul(Gsm2Smg_DD,PlotVar_V(2:4))
                       PlotVar_V(2:4) = PlotVar_V(2:4) * Si2No_V(UnitX_)
                       Xyz_D = PlotVar_V(2:4)
                       if(Coord == 'GM') Xyz_D = matmul(Smg2Gsm_DD,Xyz_D)
                       write(UnitTmp_, *) Xyz_D,Lat,Lon,OC
                    end do
                    if(MapDown)then
                       Xyz_D=PlotVar_V(2:4)
                       do i=nTP,0,-1
                          ! Map from rBody to the ionosphere
                          call map_planet_field(time_simulation, Xyz_D, 'SMG NORM', &
                               Radius+i*.1, XyzIono_D, iHemisphere)
                          if(Coord == 'GM') XyzIono_D = matmul(Smg2Gsm_DD,XyzIono_D)
                          write(UnitTmp_, *) XyzIono_D,Lat,Lon,OC
                       end do
                    end if
                 end if

                 !\\
                 ! start new line counters
                 k=nLine
                 iStart = iPoint
                 iLon=1+((nLine-1)/nLat)
                 iLat=nLine-(iLon-1)*nLat
                 Lon=IE_lon(iLon)
                 Lat=IE_lat(iLat)
                 if(NS == 'N')then
                    if(Lat<0.)k=0
                 else
                    if(Lat>0.)k=0
                 end if
              end if
           end do

           !\\
           ! finish last line
           if(k/=0)then
              iEnd = nPoint
              MapDown = .false.
              Xyz_D=PlotVar_VI(2:4,iEnd) * Si2No_V(UnitX_)
              if(sqrt(Xyz_D(1)**2 + Xyz_D(2)**2 + Xyz_D(3)**2)<1.5*rBody)MapDown = .true.
              j=(1+iEnd-iStart)+(nTP+1)
              if(MapDown)j=j+(nTP+1)
              OC=-1; if(MapDown)OC=2
              if(MapDown .and. LonOC/=Lon)OC=1
!\
!              write(UnitTmp_,'(a,2f7.2,a,a,i8,a)') 'ZONE T="IEB ll=',Lat,Lon,'"', &
!                   ', I=',j,', J=1, K=1, ZONETYPE=ORDERED, DATAPACKING=POINT'
!-
              write(UnitTmp_,'(a,2f7.2,a,a,f7.2,a,i8,a)') 'ZONE T="IEB ll=',Lat,Lon,'"', &
                   ', STRANDID=1, SOLUTIONTIME=',Lon, &
                   ', I=',j,', J=1, K=1, ZONETYPE=ORDERED, DATAPACKING=POINT'
!/
              write(stmp,'(f8.2)')Lat
              write(UnitTmp_,'(a,a,a)') 'AUXDATA LAT="',trim(adjustl(stmp)),'"'
              write(stmp,'(f8.2)')Lon
              write(UnitTmp_,'(a,a,a)') 'AUXDATA LON="',trim(adjustl(stmp)),'"'

              ! Convert to SMG Cartesian coordinates on the surface of the ionosphere
              Theta = cDegToRad*(90.0 - Lat)     
              Phi = cDegToRad*Lon
              call sph_to_xyz(Radius, Theta, Phi, XyzIono_D)
              Xyz_D=XyzIono_D
              if(Coord == 'GM') Xyz_D = matmul(Smg2Gsm_DD,Xyz_D)
              write(UnitTmp_, *) Xyz_D,Lat,Lon,OC
              do i=1,nTP
                 ! Map from the ionosphere to rBody
                 call map_planet_field(time_simulation, XyzIono_D, 'SMG NORM', &
                      Radius+i*.1, Xyz_D, iHemisphere)
                 if(Coord == 'GM') Xyz_D = matmul(Smg2Gsm_DD,Xyz_D)
                 write(UnitTmp_, *) Xyz_D,Lat,Lon,OC
              end do
              do i=iStart,iEnd
                 ! Convert vectors to SM coordinates
                 PlotVar_V = PlotVar_VI(:, i)
                 PlotVar_V(2:4) = matmul(Gsm2Smg_DD,PlotVar_V(2:4))
                 PlotVar_V(2:4) = PlotVar_V(2:4) * Si2No_V(UnitX_)
                 Xyz_D = PlotVar_V(2:4)
                 if(Coord == 'GM') Xyz_D = matmul(Smg2Gsm_DD,Xyz_D)
                 write(UnitTmp_, *) Xyz_D,Lat,Lon,OC
              end do
              if(MapDown)then
                 Xyz_D=PlotVar_V(2:4)
                 do i=nTP,0,-1
                    ! Map from the ionosphere to rBody
                    call map_planet_field(time_simulation, Xyz_D, 'SMG NORM', &
                         Radius+i*.1, XyzIono_D, iHemisphere)
                    if(Coord == 'GM') XyzIono_D = matmul(Smg2Gsm_DD,XyzIono_D)
                    write(UnitTmp_, *) XyzIono_D,Lat,Lon,OC
                 end do
              end if
           end if

           close(UnitTmp_)
        end do

        deallocate(PlotVar_VI)
     end if
     call line_clean
  end if

  if(allocated(RayIntegral_VII)) deallocate(RayIntegral_VII)
  if(allocated(RayResult_VII))   deallocate(RayResult_VII)

  if(DoTest)write(*,*)NameSub,': finished'
end subroutine ieb_plot
