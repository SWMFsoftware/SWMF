!^CFG COPYRIGHT UM
!^CFG FILE RAYTRACE
!=============================================================================

subroutine ray_trace

  ! This parallel ray tracing algorithm was developed at the U.of M.
  ! by G. Toth and D. De Zeeuw. An overview of the scheme can be found in
  ! 
  ! D. L. De Zeeuw, S. Sazykin, R. A. Wolf, T. I. Gombosi,
  ! A. J. Ridley, G. T\'oth, 2004,\\
  ! Journal of Geophysical Research, 109, 12219,
  !
  ! Details of the algorithm are to be published later

  use ModMain,     ONLY: n_step, iNewGrid, iNewDecomposition, &
       time_simulation, TypeCoordSystem
  use CON_axes,    ONLY: transform_matrix
  use ModRaytrace, ONLY: init_mod_raytrace, GmSm_DD, &
       UseAccurateTrace, DnRaytrace, r_Raytrace, R2_Raytrace
  use ModPhysics,  ONLY: rBody
  implicit none

  ! remember last call and the last grid number
  integer :: n_last=-1, iLastGrid=-1, iLastDecomposition=-1

  logical :: oktest, oktest_me
  !-------------------------------------------------------------------------

  call set_oktest('ray_trace',oktest,oktest_me)

  if(oktest_me)then
     write(*,*)'GM ray_trace: n_last,n_step         =',n_last,n_step
     write(*,*)'GM ray_trace: iLastGrid,iNewGrid    =',iLastGrid,iNewGrid
     write(*,*)'GM ray_trace: iLastDecomp,iNewDecomp=',&
          iLastDecomposition,iNewDecomposition
  end if

  if(  n_last + DnRaytrace > n_step   .and. &
       iLastGrid          == iNewGrid .and. &
       iLastDecomposition == iNewDecomposition) RETURN

  ! Remember this call
  n_last=n_step; iLastGrid = iNewGrid; iLastDecomposition = iNewDecomposition

  call timing_start('ray_trace')

  call init_mod_raytrace

  ! Initialize R_raytrace, R2_raytrace to body radius
  R_raytrace  = rBody
  R2_raytrace = R_raytrace**2

  ! Transformation matrix between the SM(G) and GM coordinates
  GmSm_DD = transform_matrix(time_simulation,'SMG',TypeCoordSystem)

  if(UseAccurateTrace)then
     call ray_trace_accurate
  else
     call ray_trace_fast
  end if

  call timing_stop('ray_trace')

end subroutine ray_trace

!==============================================================================
subroutine ray_trace_fast

  use ModProcMH
  use ModMain
  use ModAdvance,  ONLY : Bx_, Bz_, State_VGB
  use ModParallel, ONLY : NOBLK, neiLEV
  use ModGeometry, ONLY : x_BLK, y_BLK, z_BLK, R_BLK, Rmin_BLK, &
       dx_BLK, dy_BLK, dz_BLK, true_cell
  use ModRaytrace
  use ModMpi
  implicit none

  ! Iteration parameters
  integer, parameter :: ray_iter_max=150
  integer :: ray_iter
  logical :: done_me, done
  real    :: dray_min

  real :: qqray(3)

  ! Minimum value of B for which integration of field lines makes any sense
  real, parameter :: smallB=1.e-8

  ! True if Rmin_BLK < R_raytrace
  logical :: check_inside

  ! Face index for the final point of the ray 
  integer :: iface

  ! Control volume limits in local coordinates
  real, dimension(3), parameter :: &
       xmin=(/   0.5,   0.5,   0.5/),&
       xmax=(/nI+0.5,nJ+0.5,nK+0.5/)

  ! Stride for ix
  integer :: i_stride

  ! Current position of ray in normalized and physical coordinates
  real, dimension(3) :: x, xx

  ! Radial distance and square of it: r2=sum(xx**2)
  real :: r2

  ! Cell indices corresponding to current or final x position
  integer :: i1,j1,k1,i2,j2,k2

  ! Distance between x and i1,j1,k1, and i2,j2,k2
  real :: dx1, dy1, dz1, dx2, dy2, dz2

  ! Weights for surface interpolation
  real :: weight(4)

  ! Cell indices
  integer :: i,j,k

  ! Indices corresponding to the starting point of the ray
  integer :: ix,iy,iz

  ! Current block and direction indices
  integer :: iBLK, iRay

  ! Testing and timing
  logical :: oktest, oktest_me, oktime, oktime_me
  integer :: loc(3)

  integer :: iError, iError1=-1
  !----------------------------------------------------------------------------

  call set_oktest('ray_trace',oktest,oktest_me)
  call set_oktest('time_ray_trace',oktime,oktime_me)
  if(oktime)call timing_reset('ray_pass',2)

  oktest_ray = .false.

  Bxyz_DGB = State_VGB(Bx_:Bz_,:,:,:,:)
  call message_pass_cells8(.false.,.false.,.false.,3,Bxyz_DGB)

  ! Initial values !!! Maybe LOOPRAY would be better??

  rayface=NORAY
  ray=NORAY

  do iBLK=1,nBlockMax
     if(unusedBLK(iBLK))then
        ! rayface in unused blocks is assigned to NORAY-1.
        rayface(:,:,:,:,:,iBLK)=NORAY-1.
        CYCLE
     end if
     ! Inner points of rayface should never be used, assign them to OPEN
     ! so that checking for blocks with fully open rays becomes easy
     rayface(:,:,2:nI,2:nJ,2:nK,iBLK)=OPENRAY
     
     ! Set rayface=OPENRAY at outer boundaries
     if(neiLEV(east_ ,iBLK)==NOBLK)rayface(:,:,   1,:,:,iBLK)=OPENRAY
     if(neiLEV(west_ ,iBLK)==NOBLK)rayface(:,:,nI+1,:,:,iBLK)=OPENRAY
     if(neiLEV(south_,iBLK)==NOBLK)rayface(:,:,:,   1,:,iBLK)=OPENRAY
     if(neiLEV(north_,iBLK)==NOBLK)rayface(:,:,:,nJ+1,:,iBLK)=OPENRAY
     if(neiLEV(bot_  ,iBLK)==NOBLK)rayface(:,:,:,:,   1,iBLK)=OPENRAY
     if(neiLEV(top_  ,iBLK)==NOBLK)rayface(:,:,:,:,nK+1,iBLK)=OPENRAY

  end do
  if(oktest_me)write(*,*)'ray_trace initialized ray and rayface arrays'

  ! Interpolate the B1 field to the nodes
  do iBLK=1, nBlockMax
     if(unusedBLK(iBLK))CYCLE

     do k=1,nK+1; do j=1,nJ+1; do i=1,nI+1;
        bb_x(i,j,k,iBLK)=sum(Bxyz_DGB(x_,i-1:i,j-1:j,k-1:k,iBLK))*0.125
        bb_y(i,j,k,iBLK)=sum(Bxyz_DGB(y_,i-1:i,j-1:j,k-1:k,iBLK))*0.125
        bb_z(i,j,k,iBLK)=sum(Bxyz_DGB(z_,i-1:i,j-1:j,k-1:k,iBLK))*0.125

        !!if(abs(bb_x(i,j,k,iBLK))<cTiny)bb_x(i,j,k,iBLK)=cTiny
        !!if(abs(bb_y(i,j,k,iBLK))<cTiny)bb_y(i,j,k,iBLK)=cTiny
        !!if(abs(bb_z(i,j,k,iBLK))<cTiny)bb_z(i,j,k,iBLK)=cTiny
     end do; end do; end do

  end do ! iBLK

  ! Average node values between shared faces
  call pass_and_average_nodes(.true.,bb_x)
  call pass_and_average_nodes(.true.,bb_y)
  call pass_and_average_nodes(.true.,bb_z)

  if(oktest_me)write(*,*)'rayface normalized B'
  if(oktime.and.iProc==0)then
     write(*,'(a)',ADVANCE='NO') 'setup and normalization:'
     call timing_show('ray_trace',1)
  end if

  if(oktest_me)write(*,*)'ray_trace starting iterations to obtain rayface'

  ! Iterate
  dray_min=rIonosphere*1.0e-6
  UsePreferredInterpolation = .false.
  ray_iter=0
  do

     if(oktest_me)write(*,*)'ray_iter=',ray_iter

     if(ray_iter>=ray_iter_max)exit

     ! Store rayface into ray so we can see if there is any change
     ray(:,:,:,:,:,1:nBlockMax)=rayface(:,:,:,:,:,1:nBlockMax)

     do iBLK=1,nBlockMax
        if(unusedBLK(iBLK))CYCLE

        ! Flag cells inside the ionosphere if necessary
        check_inside=Rmin_BLK(iBLK)<R_raytrace

        do iz=1,nK+1
           ! Exclude outer boundaries
           if(neiLEV(bot_,iBLK)==NOBLK.and.iz==   1)CYCLE
           if(neiLEV(top_,iBLK)==NOBLK.and.iz==nK+1)CYCLE
           do iy=1,nJ+1
              ! Exclude outer boundaries
              if(neiLEV(south_,iBLK)==NOBLK.and.iy==   1)CYCLE
              if(neiLEV(north_,iBLK)==NOBLK.and.iy==nJ+1)CYCLE

              ! Exclude inside points
              if(iz>1.and.iz<nK+1.and.iy>1.and.iy<nJ+1)then
                 ! iy and iz are inside, do endpoints only in ix
                 i_stride=nI
              else
                 ! iy or iz are on the surface, do all ix
                 i_stride=1
              end if

              do ix=1,nI+1,i_stride
                 ! Exclude outer boundaries
                 if(neiLEV(east_,iBLK)==NOBLK.and.ix==   1)CYCLE
                 if(neiLEV(west_,iBLK)==NOBLK.and.ix==nI+1)CYCLE

                 !oktest_ray = oktest_me .and. BLKtest==iBLK .and. &
                 !     ix==Itest.and.iy==Jtest.and.iz==Ktest

                 !oktest_ray = oktest .and. &
                 !     abs(x_BLK(ix,iy,iz,iBLK)-0.5*dx_BLK(iBLK)-xTest)<0.01 .and.&
                 !     abs(y_BLK(ix,iy,iz,iBLK)-0.5*dy_BLK(iBLK)-yTest)<0.01 .and.&
                 !     abs(z_BLK(ix,iy,iz,iBLK)-0.5*dz_BLK(iBLK)-zTest)<0.01

                 if(oktest_ray)write(*,*)'TESTING RAY: me,iBLK,ix,iy,iz,xx',&
                      iProc,iBLK,ix,iy,iz,&
                      x_BLK(ix,iy,iz,iBLK)-0.5*dx_BLK(iBLK),&
                      y_BLK(ix,iy,iz,iBLK)-0.5*dy_BLK(iBLK),&
                      z_BLK(ix,iy,iz,iBLK)-0.5*dz_BLK(iBLK)

                 ! Debug
                 !                  if(oktest_me) then
                 !                     oktest_ray = .true.
                 !                     write(*,*)'iBLK,ix,iy,iz=',iBLK,ix,iy,iz
                 !                  end if

                 if(ray_iter==0)then
                    do iray=1,2
                       ! Follow ray in direction iray
                       iface=follow_ray(.true.,ix-0.5,iy-0.5,iz-0.5)

                       ! Assign value to rayface
                       call assign_ray(.true.,rayface(:,iray,ix,iy,iz,iBLK))

                       ! Memorize ray integration results
                       rayend_ind(1,iray,ix,iy,iz,iBLK)=iface
                       if(iface>0)then
                          select case(iface)
                          case(east_,west_)
                             rayend_ind(2:3,iray,ix,iy,iz,iBLK)=(/j1,k1/)
                          case(south_,north_)
                             rayend_ind(2:3,iray,ix,iy,iz,iBLK)=(/i1,k1/)
                          case(top_,bot_)
                             rayend_ind(2:3,iray,ix,iy,iz,iBLK)=(/i1,j1/)
                          end select
                          rayend_pos(:,iray,ix,iy,iz,iBLK)=weight
                       end if
                    end do
!!$\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!!$                    call print_test(0)
!!$//////////////////////////////
                 else
                    do iray=1,2
                       ! Use stored values
                       iface=rayend_ind(1,iray,ix,iy,iz,iBLK)
                       if(iface>0)then
                          select case(iface)
                          case(east_)
                             i1=1; i2=1
                             j1=rayend_ind(2,iray,ix,iy,iz,iBLK); j2=j1+1
                             k1=rayend_ind(3,iray,ix,iy,iz,iBLK); k2=k1+1
                          case(west_)
                             i1=nI+1; i2=i1
                             j1=rayend_ind(2,iray,ix,iy,iz,iBLK); j2=j1+1
                             k1=rayend_ind(3,iray,ix,iy,iz,iBLK); k2=k1+1
                          case(south_)
                             j1=1; j2=1
                             i1=rayend_ind(2,iray,ix,iy,iz,iBLK); i2=i1+1
                             k1=rayend_ind(3,iray,ix,iy,iz,iBLK); k2=k1+1
                          case(north_)
                             j1=nJ+1; j2=nJ+1
                             i1=rayend_ind(2,iray,ix,iy,iz,iBLK); i2=i1+1
                             k1=rayend_ind(3,iray,ix,iy,iz,iBLK); k2=k1+1
                          case(bot_)
                             k1=1; k2=1
                             i1=rayend_ind(2,iray,ix,iy,iz,iBLK); i2=i1+1
                             j1=rayend_ind(3,iray,ix,iy,iz,iBLK); j2=j1+1
                          case(top_)
                             k1=nK+1; k2=k1
                             i1=rayend_ind(2,iray,ix,iy,iz,iBLK); i2=i1+1
                             j1=rayend_ind(3,iray,ix,iy,iz,iBLK); j2=j1+1
                          end select
                          !Debug
                          !if(oktest_ray)then
                          !   write(*,*)'before rayface_interpolate:'
                          !   write(*,*)'ix,iy,iz,iBLK=',ix,iy,iz,iBLK
                          !   write(*,*)'i1..k2=',i1,i2,j1,j2,k1,k2
                          !   write(*,*)'rayface=',&
                          !        rayface(:,iray,i1:i2,j1:j2,k1:k2,iBLK)
                          !   write(*,*)'weight=',&
                          !        rayend_pos(:,iray,ix,iy,iz,iBLK)
                          !end if

                          call rayface_interpolate(&
                               rayface(:,iray,i1:i2,j1:j2,k1:k2,iBLK),&
                               rayend_pos(:,iray,ix,iy,iz,iBLK),4,&
                               qqray)

                          rayface(:,iray,ix,iy,iz,iBLK)=qqray

                          !if(oktest_ray)then
                          !   write(*,*)'after rayface_interpolate:'
                          !   write(*,*)'ix,iy,iz,iBLK=',ix,iy,iz,iBLK
                          !   write(*,*)'i1..k2=',i1,i2,j1,j2,k1,k2
                          !   write(*,*)'rayface=',&
                          !        rayface(:,iray,i1:i2,j1:j2,k1:k2,iBLK)
                          !   write(*,*)'weight=',&
                          !        rayend_pos(:,iray,ix,iy,iz,iBLK)
                          !   write(*,*)'rayface=',&
                          !        rayface(:,iray,ix,iy,iz,iBLK)
                          !end if
                       end if
                    end do
                 end if ! ray_iter==0
              end do ! ix
           end do ! iy
        end do ! iz
     end do ! iBLK

     ! Exchange rayface information

     call timing_start('ray_pass')
     call ray_pass
     call timing_stop('ray_pass')

     !if(oktest_me)then
     !   write(*,*)'ray_trace finishing ray_pass'
     !   write(*,*)'rayface(:,1)=',rayface(:,1,itest,jtest,ktest,BLKtest)
     !   write(*,*)'rayface(:,2)=',rayface(:,2,itest,jtest,ktest,BLKtest)
     !end if

!!$\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!!$!     if(ray_iter==0)then
!!$        do iBLK=1,nBlockMax
!!$           if(unusedBLK(iBLK))CYCLE
!!$           do iz=1,nK+1; do iy=1,nJ+1; do ix=1,nI+1
!!$              call print_test(ray_iter)
!!$           end do; end do; end do
!!$        end do
!!$!     end if
!!$//////////////////////////////

     ray_iter=ray_iter+1

     if(oktime.and.iProc==0.and.ray_iter==1)then
        write(*,'(a)',ADVANCE='NO') 'first iteration:'
        call timing_show('ray_trace',1)
     end if

     ! Check if we are done by checking for significant changes in rayface
     done_me=all(abs(ray(:,:,:,:,:,1:nBlockMax)-&
          rayface(:,:,:,:,:,1:nBlockMax))<dray_min)

     call MPI_allreduce(done_me,done,1,MPI_LOGICAL,MPI_LAND,iComm,iError)

     if(Done)then
        Done_me = .true.
        do iBLK=1,nBlockMax
           if(unusedBLK(iBLK))CYCLE
           Done_me = all(rayface(1,:,:,:,:,iBLK) > LOOPRAY) !!! NORAY)
           if(.not.Done_me)EXIT
        end do
        call MPI_allreduce(Done_me,Done,1,MPI_LOGICAL,MPI_LAND,iComm,iError)
        if(Done) EXIT
        if(UsePreferredInterpolation)then
           if(iProc==0)call error_report('ray tracing, ray_iter=',&
                ray_iter+0.0,iError1,.true.)
           EXIT
        endif
        if(oktest_me)write(*,*)'Switching to UsePreferredInterpolation=.true.'
        UsePreferredInterpolation = .true.
     end if

  end do ! ray iteration

!!$\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!!$  do iBLK=1,nBlockMax
!!$     if(unusedBLK(iBLK))CYCLE
!!$     do iz=1,nK+1; do iy=1,nJ+1; do ix=1,nI+1
!!$        call print_test(999)
!!$     end do; end do; end do
!!$  end do
!!$//////////////////////////////

  ! Check for unassigned rayface in every used block
  if(oktest_me)then
     write(*,*)'ray_trace finished after ',ray_iter,' iterations'
     do iBLK=1,nBlockMax
        if(unusedBLK(iBLK))CYCLE
        do iray=1,2
           if(any(rayface(1,iray,1:nI,1:nJ,1:nK,iBLK)<BODYRAY))then
              loc=minloc(rayface(1,iray,1:nI,1:nJ,1:nK,iBLK))
              write(*,*)'LOOPRAYFACE: iray,me,loc,value,x,y,z=',&
                   iray,iProc,loc,iBLK,&
                   minval(rayface(1,iray,1:nI,1:nJ,1:nK,iBLK)),&
                   x_BLK(loc(1),loc(2),loc(3),iBLK)-0.5*dx_BLK(iBLK),&
                   y_BLK(loc(1),loc(2),loc(3),iBLK)-0.5*dy_BLK(iBLK),&
                   z_BLK(loc(1),loc(2),loc(3),iBLK)-0.5*dz_BLK(iBLK)
           end if
        end do
     end do
     if(index(Test_String,'ray_debugger')>0)call ray_debugger
  end if

  if(oktime.and.iProc==0)then
     write(*,'(i5,a)') ray_iter,' iterations:'
     call timing_show('ray_trace',1)
     call timing_show('ray_pass',2)
  end if

  if(oktest_me)write(*,*)'ray_trace starting cell center assignments'

  ! Assign face ray values to cell centers
  do iBLK=1,nBlockMax

     if(unusedBLK(iBLK))CYCLE

     ! Set flag if checking on the ionosphere is necessary
     check_inside=Rmin_BLK(iBLK)<R_raytrace

     do iray=1,2
        ! Some optimization for fully open blocks
        if(.not.check_inside)then
           if(all(rayface(1,iray,:,:,:,iBLK)==OPENRAY))then
              ray(:,iray,:,:,:,iBLK)=OPENRAY
              CYCLE
           end if
        end if

        do iz=1,nK; do iy=1,nJ; do ix=1,nI
!!$           oktest_ray = oktest_me .and. BLKtest==iBLK .and. &
!!$                ix==Itest.and.iy==Jtest.and.iz==Ktest

           if(oktest_ray)write(*,*)'TESTING'

           !Debug
           !write(*,*)'me,iBLK,ix,iy,iz=',iProc,iBLK,ix,iy,iz
           !oktest_ray = .true.

           ! Short cuts for inner and false cells
           if(R_BLK(ix,iy,iz,iBLK)<rIonosphere .or. &
                .not.true_cell(ix,iy,iz,iBLK))then
              ray(:,iray,ix,iy,iz,iBLK)=BODYRAY
              if(oktest_ray)write(*,*)'BODYRAY'
              CYCLE
           end if

           if(oktest_ray)write(*,*)'calling follow_ray'

           ! Follow ray in direction iray
           iface=follow_ray(.false.,real(ix),real(iy),real(iz))

           if(oktest_ray)write(*,*)'calling assign_ray'

           ! Assign value to ray
           call assign_ray(.false.,ray(:,iray,ix,iy,iz,iBLK))

        end do; end do; end do ! ix, iy, iz
     end do ! iray
  end do ! iBLK

  if(oktest_me)write(*,*)'ray_trace finished with ray=',&
          ray(:,:,itest,jtest,ktest,BLKtest)

  if(oktest)then
     ! Check for unassigned cell centers
     do iBLK=1,nBlockMax
        if(unusedBLK(iBLK))CYCLE
        do iray=1,2
           if(any(ray(1,iray,1:nI,1:nJ,1:nK,iBLK)<BODYRAY))then
              loc=minloc(ray(1,iray,1:nI,1:nJ,1:nK,iBLK))
              write(*,*)'LOOPRAY: iray,me,loc,value,x,y,z=',&
                   iray,iProc,loc,iBLK,&
                   minval(ray(1,iray,1:nI,1:nJ,1:nK,iBLK)),&
                   x_BLK(loc(1),loc(2),loc(3),iBLK),&
                   y_BLK(loc(1),loc(2),loc(3),iBLK),&
                   z_BLK(loc(1),loc(2),loc(3),iBLK)
           end if
        end do
     end do
  end if

  if(oktest_me)write(*,*)'ray_trace starting conversion to lat/lon'

  ! Convert x, y, z to latitude and longitude, and status
  do iBLK=1,nBlockMax
     if(unusedBLK(iBLK)) CYCLE
     do k=1,nK; do j=1,nJ; do i=1,nI
        call xyz_to_latlonstatus(ray(:,:,i,j,k,iBLK))
     end do; end do; end do
  end do

  if(oktime.and.iProc==0)then
     write(*,'(a)',ADVANCE='NO') 'Total ray tracing time:'
     call timing_show('ray_trace',1)
  end if
  call barrier_mpi
  if(oktest_me)write(*,*)'ray_trace completed.'

contains

  !===========================================================================

  subroutine ray_debugger

    ! Debug rayface values

    integer :: iPos_D(3),jX,kX,jY,kY,jZ,kZ

    do
       ! Read position
       write(*,'(a)',ADVANCE='NO')'Rayface x,y,z,iRay:'
       read(*,*)xTest,yTest,zTest,iRay
       if(xTest==0.0.and.yTest==0.0.and.zTest==0.0) EXIT

       ! Find position
       do iBLK=1,nBlockMax
          if(unusedBLK(iBLK))CYCLE
          do ix=1,nI+1
             if(abs(x_BLK(ix,1,1,iBLK)-0.5*dx_BLK(iBLK)-xTest)>0.01)CYCLE
             do iy=1,nJ+1
                if(abs(y_BLK(1,iy,1,iBLK)-0.5*dy_BLK(iBLK)-yTest)>0.01)CYCLE
                do iz=1,nK+1
                   if(abs(z_BLK(1,1,iz,iBLK)-0.5*dz_BLK(iBLK)-zTest)>0.01)&
                        CYCLE

                   ! Print information

                   write(*,*)'iProc,iBLK,ix,iy,iz=',iProc,iBLK,ix,iy,iz
                   write(*,*)' x,y,z_BLK(1,1,1),dx_BLK=',&
                        x_BLK(1,1,1,iBLK),y_BLK(1,1,1,iBLK),&
                        z_BLK(1,1,1,iBLK),dx_BLK(iBLK)

                   iPos_D=rayend_ind(:,iray,ix,iy,iz,iBLK)
                   if(iPos_D(1)>0)then
                      write(*,*)' rayface   =',rayface(:,iray,ix,iy,iz,iBLK)
                      write(*,*)' rayend_ind=',rayend_ind(:,iray,ix,iy,iz,iBLK)
                      write(*,*)' rayend_pos=',rayend_pos(:,iray,ix,iy,iz,iBLK)
                      select case(iPos_D(1))
                      case(1,2)
                         jX = 1+nI*(iPos_D(1)-1); jY = iPos_D(2); jZ = iPos_D(3)
                         kX = jX; kY = jY+1; kZ = jZ+1
                      case(3,4)
                         jY = 1+nJ*(iPos_D(1)-3); jX = iPos_D(2); jZ = iPos_D(3)
                         kX = jX+1; kY = jY; kZ = jZ+1
                      case(5,6)
                         jZ = 1+nK*(iPos_D(1)-5); jX = iPos_D(2); jY = iPos_D(3)
                         kX = jX+1; kY = jY+1; kZ = jZ
                      end select
                      write(*,*)' rayface(1,end)=',&
                           rayface(1,iray,jX:kX,jY:kY,jZ:kZ,iBLK)
                      write(*,*)' jX,kX,jY,kY,jZ,kZ=',jX,kX,jY,kY,jZ,kZ
                      write(*,*)' x,y,z(End)=',&
                           x_BLK(jx,jy,jz,iBLK)-0.5*dx_BLK(iBLK),&
                           y_BLK(jx,jy,jz,iBLK)-0.5*dy_BLK(iBLK),&
                           z_BLK(jx,jy,jz,iBLK)-0.5*dz_BLK(iBLK)
                   else
                      write(*,*)' rayend_ind=',iPos_D
                   end if
                end do
             end do
          end do
       end do
    end do

  end subroutine ray_debugger

  !===========================================================================

  subroutine print_test(inInt)
    integer, intent(in) :: inInt

    if(  abs(  4.0-(x_BLK(ix,iy,iz,iBLK)-0.5*dx_BLK(iBLK)))<0.01 .and. &
         abs(  6.0-(y_BLK(ix,iy,iz,iBLK)-0.5*dy_BLK(iBLK)))<0.01 .and. &
         abs( -6.5-(z_BLK(ix,iy,iz,iBLK)-0.5*dz_BLK(iBLK)))<0.01 )then
       write(*,'(i3,a,i3,a,i4,a,3i2,3(a,f9.2),a,6i3,a,6f10.3)') inInt, &
            ' DEBUG LOOPRAYFACE: PE=',iProc,' BLK=',iBLK,' loc=',ix,iy,iz,&
            ' x=',x_BLK(ix,iy,iz,iBLK)-0.5*dx_BLK(iBLK), &
            ' y=',y_BLK(ix,iy,iz,iBLK)-0.5*dy_BLK(iBLK), &
            ' z=',z_BLK(ix,iy,iz,iBLK)-0.5*dz_BLK(iBLK), &
            '   rayend_ind=',rayend_ind(:,:,ix,iy,iz,iBLK), &
            '   rayface=',rayface(:,:,ix,iy,iz,iBLK)
    end if
  end subroutine print_test

  !===========================================================================

  function follow_ray(surface_point,x_0,y_0,z_0) result(qface)

    ! Follow ray starting at initial position x_0,y_0,z_0 in direction iray
    ! until we hit the wall of the control volume or the ionosphere.
    ! Return 1,2,3,4,5,6 if the ray hit the east,west,south,north,bot,top walls
    ! Return ray_iono_   if the ray hit the ionosphere
    ! Return ray_loop_   if the ray did not hit anything 
    ! Return ray_out_    if the ray goes out of the box immediately
    ! Return ray_body_   if the ray goes into or is inside a body

    ! Arguments

    logical, intent(in):: surface_point
    real, intent(in)   :: x_0,y_0,z_0

    ! Result

    integer :: qface

    ! Local variables

    ! Initial and mid point coordinates and bb field
    real, dimension(3) :: x_ini, x_mid, b_ini, b_mid, xx_ini
    real :: r, r_ini

    ! dx is the difference between 1st and 2nd order RK to estimate accuracy
    ! dx_opt is the required accuracy, dx_rel=dx/dx_opt
    real :: dx_rel, dx_opt

    ! Ray length, max, step size, limits, next step size for backup to surface
    real :: l, lmax, dl, dl_max, dl_min, dl_next, dl_tiny, dl_back

    ! counter for ray integration
    integer :: nsegment 

    ! Counter for entering follow_ray_iono
    integer :: n_iono

    !--------------------------------------------------------------------------

    if(oktest_ray)&
         write(*,*)'follow_ray: me,iBLK,surface_point,x_0,y_0,z_0,iray=',&
         iProc,iBLK,surface_point,x_0,y_0,z_0,iray

    ! Step size limits
    dl_max=1.0
    dl_min=0.05
    dl_tiny=1.e-6

    ! Initial value
    dl_next=sign(dl_max,1.5-iray)

    ! Accuracy in terms of x in normalized coordinates
    dx_opt=0.01

    ! Length and maximum length of ray within control volume
    l=0
    lmax=10*maxval(xmax-xmin)
    nsegment=0
    n_iono=0

    ! Initial position
    x(1)=x_0
    x(2)=y_0
    x(3)=z_0

    ! Integration loop
    do
       ! Check if we are inside the ionosphere
       if(check_inside)then
          ! Convert x to real coordinates xx

          xx(1)=x_BLK(1,1,1,iBLK)+dx_BLK(iBLK)*(x(1)-1.)
          xx(2)=y_BLK(1,1,1,iBLK)+dy_BLK(iBLK)*(x(2)-1.)
          xx(3)=z_BLK(1,1,1,iBLK)+dz_BLK(iBLK)*(x(3)-1.)

          r2=sum(xx**2)

          if(r2<=R2_raytrace)then

             if(oktest_ray)write(*,*)&
                  'Inside R_raytrace at me,iBLK,nsegment,x,xx=',&
                  iProc,iBLK,nsegment,x,xx

             if(r2<=rIonosphere2)then
                if(nsegment==0)then
                   qface=ray_body_
                   if(oktest_ray)write(*,*)&
                        'Initial point inside rIonosphere at me,iBLK,xx=',&
                        iProc,iBLK,xx
                else
                   r=sqrt(r2)
                   xx_ini(1)=x_BLK(1,1,1,iBLK)+dx_BLK(iBLK)*(x_ini(1)-1.)
                   xx_ini(2)=y_BLK(1,1,1,iBLK)+dy_BLK(iBLK)*(x_ini(2)-1.)
                   xx_ini(3)=z_BLK(1,1,1,iBLK)+dz_BLK(iBLK)*(x_ini(3)-1.)
                   r_ini=sqrt(sum(xx_ini**2))
                   ! Interpolate to the surface linearly along last segment
                   xx=(xx*(r_ini-rIonosphere)+xx_ini*(rIonosphere-r)) &
                        /(r_ini-r)
                   ! Normalize xx in radial direction
                   xx=rIonosphere*xx/sqrt(sum(xx**2))
                   x=xx
                   qface=ray_iono_
                end if
                EXIT
             end if

             ! Try mapping down to rIonosphere if we haven't tried yet
             if(n_iono<1)then
                if(follow_ray_iono())then
                   x=xx
                   qface=ray_iono_
                   EXIT
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
       call interpolate_bb_node(x_ini,b_ini)
       x_mid=x_ini+0.5*dl*b_ini
          
       ! Check if the ray is pointing outwards
       if(nsegment==0.and.surface_point)then
          if(oktest_ray)write(*,*)'me,iBLK,x_ini,b_ini=', &
               iProc,iBLK,x_ini,b_ini

          if(any(x_mid<xmin) .or. any(x_mid>xmax))then
             qface=ray_out_
             if(oktest_ray)then
                write(*,*)'me,iBLK,x_mid=',iProc,iBLK,x_mid
                write(*,*)'ray points outwards: me,iBLK,dl,xx=', &
                     iProc,iBLK,dl,&
                     x_BLK(1,1,1,iBLK)+dx_BLK(iBLK)*(x_mid(1)-1.),&
                     y_BLK(1,1,1,iBLK)+dy_BLK(iBLK)*(x_mid(2)-1.),&
                     z_BLK(1,1,1,iBLK)+dz_BLK(iBLK)*(x_mid(3)-1.)
             end if
             RETURN
          end if
       end if

       do
          ! Full step
          call interpolate_bb_node(x_mid,b_mid)

          ! Calculate the difference between 1st and 2nd order integration
          ! and take ratio relative to dx_opt
          dx_rel=abs(dl)*maxval(abs(b_mid-b_ini))/dx_opt

          if(oktest_ray.and.okdebug)&
               write(*,*)'me,iBLK,x_mid,b_mid,dx_rel=', &
               iProc,iBLK,x_mid,b_mid,dx_rel

          ! Make sure that dl does not change more than a factor of 2 or 0.5
          dx_rel=max(0.5,min(2.,dx_rel))

          if(dx_rel>1.)then
             ! Not accurate enough, decrease dl if possible

             if(abs(dl)<=dl_min+dl_tiny)then
                ! Cannot reduce dl further
                dl_next=dl
                EXIT
             end if

             dl = sign(max(dl_min,abs(dl)/(dx_rel+0.001)),dl)

             ! New mid point using the reduced dl
             x_mid=x_ini+0.5*dl*b_ini
          
             if(oktest_ray.and.okdebug)&
                  write(*,*)'new decreased dl: me,iBLK,dl=', &
                  iProc,iBLK,dl
          else
             ! Too accurate, increase dl if possible
             if(abs(dl)<dl_max-dl_tiny)then

                dl_next = sign(min(dl_max,abs(dl)/sqrt(dx_rel)),dl)

                if(oktest_ray.and.okdebug)&
                     write(*,*)'new increased dl_next: me,iBLK,dl_next=', &
                     iProc,iBLK,dl_next

             end if

             EXIT
          end if
       end do

       x=x_ini+b_mid*dl

       nsegment=nsegment+1
       l=l+abs(dl)

       if(oktest_ray.and.okdebug)&
            write(*,*)'me,iBLK,nsegment,l,x=', &
            iProc,iBLK,nsegment,l,x

       ! Check if the ray hit the wall of the control volume
       if(any(x<xmin) .or. any(x>xmax))then

          ! Hit the wall, backup so that x is almost exactly on the wall
          ! just a little bit outside. Only if nsegment is more than 1!
          if(nsegment > 1)then
             dl_back = dl*maxval(max(xmin-x,x-xmax)/(abs(x-x_ini)+dl_tiny))
             x=x-dl_back*b_mid
          end if

          ! Find out which wall the ray hit
          if    (x(1)<=xmin(1))then; qface=1
          elseif(x(2)<=xmin(2))then; qface=3
          elseif(x(3)<=xmin(3))then; qface=5
          elseif(x(1)>=xmax(1))then; qface=2
          elseif(x(2)>=xmax(2))then; qface=4
          elseif(x(3)>=xmax(3))then; qface=6
          else
             write(*,*)'Error in follow_ray for me,iBLK,ix,iy,iz=',&
                  iProc,iBLK,ix,iy,iz
             write(*,*)'nsegment,x,dl,dl_back=',nsegment,x,dl,dl_back
             call stop_mpi('GM_follow_ray: Hit wall but which one?')
          end if

          ! Make sure that x is not outside the control volume
          x=max(xmin+dl_tiny,x)
          x=min(xmax-dl_tiny,x)

          EXIT
       end if

       ! Check if we have integrated for too long
       if(l>lmax)then
          ! Seems to be a closed loop within a block
          if(oktest_ray)then
             write(*,*)'CLOSED LOOP at me,iBLK,ix,iy,iz,x,xx=',&
                  iProc,iBLK,ix,iy,iz,x,&
                  x_BLK(1,1,1,iBLK)+dx_BLK(iBLK)*(x(1)-1.),&
                  y_BLK(1,1,1,iBLK)+dy_BLK(iBLK)*(x(2)-1.),&
                  z_BLK(1,1,1,iBLK)+dz_BLK(iBLK)*(x(3)-1.)
          end if

          qface=ray_loop_
          EXIT
       end if

    end do

    if(oktest_ray)write(*,*)'Finished follow_ray at me,iBLK,nsegment,qface,x,xx=',&
         iProc,iBLK,nsegment,qface,x,&
         x_BLK(1,1,1,iBLK)+dx_BLK(iBLK)*(x(1)-1.),&
         y_BLK(1,1,1,iBLK)+dy_BLK(iBLK)*(x(2)-1.),&
         z_BLK(1,1,1,iBLK)+dz_BLK(iBLK)*(x(3)-1.)

  end function follow_ray

  !===========================================================================

  subroutine interpolate_bb(qx,qb)

    ! Obtain normalized bb field at normalized location qx and put it into qb

    real, intent(in) :: qx(3)
    real, intent(out):: qb(3)

    !-------------------------------------------------------------------------

    ! Determine cell indices corresponding to location qx

    i1=floor(qx(1)); i2=i1+1
    j1=floor(qx(2)); j2=j1+1
    k1=floor(qx(3)); k2=k1+1

    if(i1<-1.or.i2>nI+2.or.j1<-1.or.j2>nJ+2.or.k1<-1.or.k2>nK+2)then
       write(*,*)'interpolate_bb: iProc, iBLK, qx=',iProc,iBLK, qx
       call stop_mpi('ERROR in interpolate_bb: location out of bounds')
    endif

    ! Distance relative to the cell centers

    dx1=qx(1)-i1; dx2=1.-dx1
    dy1=qx(2)-j1; dy2=1.-dy1
    dz1=qx(3)-k1; dz2=1.-dz1

    qb(1)=interpolate_bb1_node(bb_x)
    qb(2)=interpolate_bb1_node(bb_y)
    qb(3)=interpolate_bb1_node(bb_z)

  end subroutine interpolate_bb

  !===========================================================================

  real function interpolate_bb1(qbb)

    real, dimension(-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), &
         intent(in):: qbb

    !-------------------------------------------------------------------------

    ! Bilinear interpolation in 3D

    interpolate_bb1=&
         dx1*(   dy1*(   dz1*qbb(i2,j2,k2,iBLK)+&
                         dz2*qbb(i2,j2,k1,iBLK))+&
                 dy2*(   dz1*qbb(i2,j1,k2,iBLK)+&
                         dz2*qbb(i2,j1,k1,iBLK)))+&
         dx2*(   dy1*(   dz1*qbb(i1,j2,k2,iBLK)+&
                         dz2*qbb(i1,j2,k1,iBLK))+&
                 dy2*(   dz1*qbb(i1,j1,k2,iBLK)+&
                         dz2*qbb(i1,j1,k1,iBLK)))

  end function interpolate_bb1
  !==========================================================================
  function interpolate_bb_v(nVar,qbb)
    integer,intent(in)::nVar
    real, dimension(nVar,-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), &
         intent(in):: qbb
    real, dimension(nVar)::interpolate_bb_v

    !-------------------------------------------------------------------------

    ! Bilinear interpolation in 3D

    interpolate_bb_v=&
         dx1*(   dy1*(   dz1*qbb(:,i2,j2,k2,iBLK)+&
         dz2*qbb(:,i2,j2,k1,iBLK))+&
         dy2*(   dz1*qbb(:,i2,j1,k2,iBLK)+&
         dz2*qbb(:,i2,j1,k1,iBLK)))+&
         dx2*(   dy1*(   dz1*qbb(:,i1,j2,k2,iBLK)+&
         dz2*qbb(:,i1,j2,k1,iBLK))+&
         dy2*(   dz1*qbb(:,i1,j1,k2,iBLK)+&
         dz2*qbb(:,i1,j1,k1,iBLK)))

  end function interpolate_bb_v

  !===========================================================================

  subroutine interpolate_bb_node(qx,qb)

    ! Obtain normalized bb field at normalized location qx and put it into qb
    ! Interpolate B1 from nodes, take B0 from analytic expression

    real, intent(in) :: qx(3)
    real, intent(out):: qb(3)
    real :: qbD

    !-------------------------------------------------------------------------

    ! Determine cell indices corresponding to location qx

    i1=floor(qx(1)+0.5); i2=i1+1
    j1=floor(qx(2)+0.5); j2=j1+1
    k1=floor(qx(3)+0.5); k2=k1+1

    if(i1<0.or.i2>nI+2.or.j1<0.or.j2>nJ+2.or.k1<0.or.k2>nK+2)then
       write(*,*)'interpolate_bb_node: iProc, iBLK, qx=',iProc,iBLK, qx
       call stop_mpi('ERROR in interpolate_bb_node: location out of bounds')
    endif

    ! Get B0 values for location

    xx(1)=x_BLK(1,1,1,iBLK)+dx_BLK(iBLK)*(qx(1)-1.)
    xx(2)=y_BLK(1,1,1,iBLK)+dy_BLK(iBLK)*(qx(2)-1.)
    xx(3)=z_BLK(1,1,1,iBLK)+dz_BLK(iBLK)*(qx(3)-1.)
    if(UseB0)then
       call get_b0(xx(1),xx(2),xx(3),qb)
    else
       qb=0.00
    end if

    ! Make sure that the interpolation uses inside indexes only

    i1=max(1,i1)   ; j1=max(1,j1);    k1=max(1,k1)
    i2=min(nI+1,i2); j2=min(nJ+1,j2); k2=min(nK+1,k2)

    ! Distances relative to the nodes

    dx1=qx(1)+0.5-i1; dx2=1.-dx1
    dy1=qx(2)+0.5-j1; dy2=1.-dy1
    dz1=qx(3)+0.5-k1; dz2=1.-dz1

    ! Add in node interpolated B1 values and take aspect ratios into account

    qb(1)=(qb(1)+interpolate_bb1_node(bb_x))/dx_BLK(iBLK)
    qb(2)=(qb(2)+interpolate_bb1_node(bb_y))/dy_BLK(iBLK)
    qb(3)=(qb(3)+interpolate_bb1_node(bb_z))/dz_BLK(iBLK)

    ! Normalize
    qbD=sqrt(qb(1)**2 + qb(2)**2 + qb(3)**2)

    if(qbD>smallB)then
       qb=qb/qbD
    else
       qb=0.
    end if

  end subroutine interpolate_bb_node

  !===========================================================================

  real function interpolate_bb1_node(qbb)

    real, dimension(1:nI+1,1:nJ+1,1:nK+1,nBLK), &
         intent(in):: qbb

    !-------------------------------------------------------------------------

    ! Bilinear interpolation in 3D

    interpolate_bb1_node=&
         dx1*(   dy1*(   dz1*qbb(i2,j2,k2,iBLK)+&
                         dz2*qbb(i2,j2,k1,iBLK))+&
                 dy2*(   dz1*qbb(i2,j1,k2,iBLK)+&
                         dz2*qbb(i2,j1,k1,iBLK)))+&
         dx2*(   dy1*(   dz1*qbb(i1,j2,k2,iBLK)+&
                         dz2*qbb(i1,j2,k1,iBLK))+&
                 dy2*(   dz1*qbb(i1,j1,k2,iBLK)+&
                         dz2*qbb(i1,j1,k1,iBLK)))

  end function interpolate_bb1_node

  !===========================================================================

  logical function follow_ray_iono()

    ! Follow ray inside ionosphere starting from xx which is given in
    ! real coordinates and use analytic mapping.
    ! On return xx contains the final coordinates.
    ! Return true if it was successfully integrated down to rIonosphere,
    ! return false if the ray exited R_raytrace or too many integration 
    ! steps were done

    use ModMain,     ONLY: Time_Simulation
    use ModPhysics,  ONLY: DipoleStrengthSi ! only the sign of dipole is needed
    use CON_planet_field, ONLY: map_planet_field

    integer :: iHemisphere
    real    :: x_D(3)
    !---------------------------------------------------------------------
    call map_planet_field(Time_Simulation, xx, TypeCoordSystem//' NORM', &
         rIonosphere, x_D, iHemisphere)

    if(iHemisphere==0)then
       write(*,*)'iHemisphere==0 for xx=',xx
       write(*,*)'iBLK, iRay=',iBLK,iRay
       call stop_mpi('ERROR in follow_ray_iono')
    end if

    if(iHemisphere*DipoleStrengthSi*sign(1.0,1.5-iray) < 0.0)then
       xx = x_D
       follow_ray_iono = .true.
    else
       follow_ray_iono = .false.
    end if

  end function follow_ray_iono

  !===========================================================================
  subroutine evaluate_bb(qx,qb)

    ! Obtain normalized bb field at true location qx and put it into qb

    real, intent(in) :: qx(3)
    real, intent(out):: qb(3)

    !-------------------------------------------------------------------------

    ! Get B0
    call get_b0(qx(1),qx(2),qx(3),qb)

    ! Take aspect ratio of cells into account
    qb(1)=qb(1)/dx_BLK(iBLK)
    qb(2)=qb(2)/dy_BLK(iBLK)
    qb(3)=qb(3)/dz_BLK(iBLK)

    if(sum(abs(qb))==0.0)then
       write(*,*)'GM_ERROR in ray_trace::evaluate_bb: qb==0 at qx=',qx
       call stop_mpi('GM_ERROR in ray_trace::evaluate_bb')
    end if

    ! Normalize
    qb=qb/sqrt(sum(qb**2))

  end subroutine evaluate_bb

  !===========================================================================

  subroutine assign_ray(surface_point,qray)

    ! Assign value to qray(3) based on ray intersection 
    ! given by the global variables iface and position x(3)
    !
    ! iray is 1 if ray points in positive B direction and 2 otherwise
    !
    ! surface_point is true if the ray was started from the block face
    ! and false if it was started from a cell center

    logical, intent(in) :: surface_point
    ! Called with a segment of rayface array and it is used here to get qray
    real, intent(inout) :: qray(3)

    ! Temporary variable
    real :: qqray(3)

    ! Local variables

    ! Distances between x and the 4 grid points used for interpolation
    real :: d1,e1,d2,e2

    !-------------------------------------------------------------------------

    if(oktest_ray)write(*,*)&
         'assign_ray starting with surface_point, iray, iface=',&
         surface_point,iray,iface

    select case(iface)
    case(ray_out_)
       ! The ray points outward
       qray=OUTRAY
       if(oktest_ray)write(*,*)'assign_ray finished with qray=OUTRAY'
       return
    case(ray_loop_)
       ! The ray did not hit the wall of the block
       qray=LOOPRAY
       if(oktest_ray)write(*,*)'assign_ray finished with qray=LOOPRAY'
       return
    case(ray_body_)
       ! The ray hit a body
       qray=BODYRAY
       if(oktest_ray)write(*,*)'assign_ray finished with qray=BODYRAY'
       return       
    case(ray_iono_)
       ! The ray hit the ionosphere 
       qray=x
       if(oktest_ray)write(*,*)&
            'assign_ray finished with qray on ionosphere, qray=',qray
       return
    case(east_,west_)
       if(iface==east_)then
          i1=1
       else
          i1=nI+1
       endif
       i2=i1
       j1=floor(x(2)-xmin(2))+1; j2=j1+1
       k1=floor(x(3)-xmin(3))+1; k2=k1+1
       d1=x(2)-j1+0.5
       e1=x(3)-k1+0.5

    case(south_,north_)
       if(iface==south_)then
          j1=1
       else
          j1=nJ+1
       endif
       j2=j1
       i1=floor(x(1)-xmin(1))+1; i2=i1+1
       k1=floor(x(3)-xmin(3))+1; k2=k1+1
       d1=x(1)-i1+0.5
       e1=x(3)-k1+0.5

    case(bot_,top_)
       ! The ray hit the bot or top wall
       if(iface==bot_)then
          k1=1
       else
          k1=nK+1
       endif
       k2=k1
       i1=floor(x(1)-xmin(1))+1; i2=i1+1
       j1=floor(x(2)-xmin(2))+1; j2=j1+1
       d1=x(1)-i1+0.5
       e1=x(2)-j1+0.5

    case default
       write(*,*)'Impossible value for iface=',iface,' at ix,iy,iz,iBLK=',&
            ix,iy,iz,iBLK
       call stop_mpi('assign_ray')
    end select

    ! Calculate bilinear interpolation weights
    d2=1.-d1; e2=1.-e1
    weight(1)=d2*e2
    weight(2)=d1*e2
    weight(3)=d2*e1
    weight(4)=d1*e1

    if(oktest_ray)write(*,*)'weight=',weight

    ! Exclude the starting point if its among the 4 interpolated cells
    if(surface_point)then
       if((ix==i1.or.ix==i2).and.(iy==j1.or.iy==j2).and.(iz==k1.or.iz==k2))then
          select case(iface)
          case(east_,west_)
             weight(iy-j1+2*(iz-k1)+1)=0.
          case(south_,north_)
             weight(ix-i1+2*(iz-k1)+1)=0.
          case(bot_,top_)
             weight(ix-i1+2*(iy-j1)+1)=0.
          end select
          ! Normalize weights
          weight=weight/sum(weight)
          if(oktest_ray)write(*,*)'Excluded point: me,iBLK,ix,iy,iz,weight=',&
               iProc,iBLK,ix,iy,iz,weight
       end if
    end if

    if(oktest_ray)&
         write(*,*)'i1,j1,k1,i2,j2,k2,d1,e1=',i1,j1,k1,i2,j2,k2,d1,e1

    call rayface_interpolate(rayface(:,iray,i1:i2,j1:j2,k1:k2,iBLK),&
         weight,4,qqray)

    qray = qqray

    if(oktest_ray)write(*,*)'assign_ray finished qray=',qray

  end subroutine assign_ray

end subroutine ray_trace_fast

!=========================================================================
subroutine rayface_interpolate(qrayface,weight,nvalue,qray)

  ! Collect weights for qrayface values that differ less than dray_max
  ! and interpolate the values corresponding to the largest weight
  ! The result is returned in qray.
  ! Note that qray and qrayface may overlap, so their intent must be inout!

  use ModRaytrace
  implicit none

  integer, intent(in)    :: nvalue
  real,    intent(inout) :: qrayface(3,nvalue)
  real,    intent(in)    :: weight(nvalue)
  real,    intent(inout) :: qray(3)

  ! Local variables

  ! Cumulated weights corresponding to various kinds of qrayface values
  real :: qweight(4), weight_sum(4), ray_first(3,4), ray_sum(3,4)

  ! Difference between qrayface values, maximum for interpolation
  real :: dray, dray_max, ValueMax

  ! Number and indices of (cummulated) qrayface values, max location
  integer :: n, i, j, loc(1)

  !-----------------------------------------------------------------------

  if(.not.UsePreferredInterpolation .and. &
       maxval(weight)-minval(weight) > 0.0001)then
     qweight(1:nvalue)=weight
  else
     ValueMax = maxval(qrayface(1,:), MASK = weight > 0.0) !!! 0.01)

     if(ValueMax < CLOSEDRAY)then
!     if(ValueMax < OPENRAY-0.01)then
        qray = ValueMax
        RETURN
     end if

     where(qrayface(1,:)>=OPENRAY-0.01)
        qweight(1:nvalue)=weight
     elsewhere
        qweight(1:nvalue)=0.
     endwhere
  end if

  if(oktest_ray)then
     write(*,*)'rayface_interpolate'
     write(*,*)'qrayface(1,:)=',qrayface(1,:)
     write(*,*)'qrayface(2,:)=',qrayface(2,:)
     write(*,*)'qrayface(3,:)=',qrayface(3,:)
     write(*,*)'weight       =',weight
     write(*,*)'qweight      =',qweight
  end if

  ! Short cuts
  if(all(qrayface(1,:)==OPENRAY))then
     ! all surrounding rays are open
     qray=OPENRAY
     if(oktest_ray)write(*,*)'rayface_interpolate finished with fully OPENRAY'
     return
  end if

  if(all(qrayface(1,:)==NORAY))then
     ! all surrounding rays are unknown
     qray=NORAY
     if(oktest_ray)write(*,*)'rayface_interpolate finished with fully NORAY'
     return
  end if

  dray_max=0.2*rIonosphere
  n=0
  do j=1,nvalue
     i=1
     do
        if(i>n)then
           ! New type of ray
           n=i
           ray_first(:,i)=qrayface(:,j)
           weight_sum(i) =qweight(j)
           if(ray_first(1,i)>CLOSEDRAY)&
                ray_sum(:,i)=qweight(j)*qrayface(:,j)
           EXIT
        end if

        ! Calculate difference between qrayface(:,j) and ray_first(:,i)
        dray=sum(abs(qrayface(:,j)-ray_first(:,i)))

        if(dray<dray_max)then
           ! Same type of ray, cummulate it

           weight_sum(i)=weight_sum(i)+qweight(j)
           if(ray_first(1,i)>CLOSEDRAY)&
                ray_sum(:,i)=ray_sum(:,i)+qweight(j)*qrayface(:,j)
           EXIT
        end if
        ! Try next type
        i=i+1

        if(i>nvalue)call stop_mpi(&
             'Impossible value for i in rayface_interpolate')
     end do ! i
  end do ! j

  if(n==1)then
     ! Only one type of ray is interpolated
     if(ray_first(1,1)>CLOSEDRAY)then
        ! get result (weight_sum can be less than 1!)
        qray=ray_sum(:,1)/weight_sum(1)
     else
        ! identical rayface values, no need to average
        qray=ray_first(:,1)
     end if
  else
     ! Take the values corresponding to the largest cummulated weight
     loc=maxloc(weight_sum(1:n))
     i=loc(1)
     if(ray_first(1,i)>CLOSEDRAY)then
        ! take average 
        qray=ray_sum(:,i)/weight_sum(i)
     else
        ! identical rayface values, no need to average
        qray=ray_first(:,i)
     end if
  end if

  if(oktest_ray)then
     write(*,*)'rayface_interpolate: weight_sum=',weight_sum(1:n)
     write(*,*)'rayface_interpolate finished with qray=',qray
  end if

end subroutine rayface_interpolate

!=========================================================================
subroutine convFaces2LatLon(rayface_in)
  use ModRaytrace, ONLY: nI, nJ, nK, xyz_to_latlonstatus
  implicit none

  real, intent(inout), dimension(3,2,1:nI+1,1:nJ+1,1:nK+1):: rayface_in

  integer :: Di,i,j,k
  !---------------------------------------------------------------------
  do k=1,nK+1
     do j=1,nJ+1
        ! Exclude inside points
        if(k>1.and.k<nK+1.and.j>1.and.j<nJ+1)then
           ! j and k are inside, do endpoints only in i
           Di=nI
        else
           Di=1
        end if
        do i=1,nI+1,Di
           call xyz_to_latlonstatus(rayface_in(:,:,i,j,k))
        end do
     end do
  end do

end subroutine convFaces2LatLon

!=============================================================================

subroutine integrate_ray(dbg,iBLK,x_0,y_0,z_0,fvol,rvol,pvol)

  ! Follow ray starting at initial position x_0,y_0,z_0 in direction 1
  ! until we hit the wall of the control volume or the ionosphere.
  ! Return dS/B
  ! x_0, y_0, and z_0 sent in in real coordinates

  use ModProcMH
  use ModAdvance, ONLY : rho_, Bx_, Bz_, P_, State_VGB
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,Rmin_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModRaytrace
  implicit none

  ! Arguments

  integer, intent(in) :: iBLK
  real, intent(in)   :: x_0,y_0,z_0
  real, intent(out) :: fvol,rvol,pvol
  logical, intent(in) :: dbg

  ! Local variables

  logical :: end_inside

  ! True if Rmin_BLK < R_raytrace
  logical :: check_inside

  ! Minimum value of B for which integration of field lines makes any sense
  real, parameter :: smallB=1.e-8

  ! Control volume limits in local coordinates
  real, dimension(3), parameter :: &
       xmin=(/   0.5,   0.5,   0.5/),&
       xmax=(/nI+0.5,nJ+0.5,nK+0.5/)

  ! Current position of ray in normalized and physical coordinates
  real, dimension(3) :: x, xx

  ! Radial distance and square of it: r2=sum(xx**2)
  real :: r2

  ! Initial and mid point coordinates and bb field
  real, dimension(3) :: x_ini, x_mid, b_ini, b_mid

  real :: bNORM,rNORM,pNORM
  real :: local_fvol,local_rvol,local_pvol

  ! dx is the difference between 1st and 2nd order RK to estimate accuracy
  ! dx_opt is the required accuracy, dx_rel=dx/dx_opt
  real :: dx_rel, dx_opt

  ! Ray length, max, step size, limits, next step size for backup to surface
  real :: l, lmax, dl, dl_max, dl_min, dl_next, dl_tiny, dl_back

  ! Distance between x and i1,j1,k1, and i2,j2,k2
  real :: dx1, dy1, dz1, dx2, dy2, dz2

  ! Direction index
  integer :: iray

  ! Cell indices corresponding to current or final x position
  integer :: i1,j1,k1,i2,j2,k2

  ! counter for ray integration
  integer :: nsegment 

  real :: amount2add
  !--------------------------------------------------------------------------

  fvol=0.; rvol=0.; pvol=0.

  ! Set flag if checking on the ionosphere is necessary
  check_inside=Rmin_BLK(iBLK)<R_raytrace

  ! Step size limits
  dl_max=0.05
  dl_min=0.01
  dl_tiny=1.e-6

  do iray=1,2
     if(dbg)write(*,*)'Starting iray=',iray

     ! Initial value
     dl_next=sign(dl_max,1.5-iray)

     ! Accuracy in terms of x in normalized coordinates
     dx_opt=0.01

     ! Length and maximum length of ray within control volume
     l=0
     lmax=10*maxval(xmax-xmin)
     nsegment=0

     ! Initial position
     x(1)=1.+(x_0-x_BLK(1,1,1,iBLK))/dx_BLK(iBLK)
     x(2)=1.+(y_0-y_BLK(1,1,1,iBLK))/dy_BLK(iBLK)
     x(3)=1.+(z_0-z_BLK(1,1,1,iBLK))/dz_BLK(iBLK)

     end_inside=.false.
     local_fvol=0.; local_rvol=0.; local_pvol=0.

     if(dbg)write(*,*)'  initial values: dl_next=',dl_next,' dx_opt=',dx_opt, &
          ' lmax=',lmax,' x=',x,' check_inside=',check_inside
     if(iray==1 .or. check_inside)then
        call do_integration

        if(iray==1 .or. (iray==2 .and. end_inside)) then
           fvol=fvol+local_fvol
           rvol=rvol+local_rvol
           pvol=pvol+local_pvol
        end if
     end if

  end do


contains

  subroutine do_integration

    ! Integration loop
    do
       if(dbg)write(*,*)'  loop iter',nsegment,x

       ! Check if we are inside the ionosphere
       if(check_inside)then
          ! Convert x to real coordinates xx
          xx(1)=x_BLK(1,1,1,iBLK)+dx_BLK(iBLK)*(x(1)-1.)
          xx(2)=y_BLK(1,1,1,iBLK)+dy_BLK(iBLK)*(x(2)-1.)
          xx(3)=z_BLK(1,1,1,iBLK)+dz_BLK(iBLK)*(x(3)-1.)

          r2=sum(xx**2)

          if(r2<R2_raytrace)then
             if(dbg)write(*,*)'  inside raytrace limit, stopping.'

             if(nsegment==0)then
                !starting inside, leave
                EXIT
             end if

             end_inside=.true.

             EXIT
          end if
       end if

       ! Integrate with 2nd order scheme
       dl=dl_next
       x_ini=x

       ! Half step
       call interpolate_bbN(x_ini,b_ini,bNORM,rNORM,pNORM)

       do
          x_mid=x_ini+0.5*dl*b_ini

          ! Full step
          call interpolate_bbN(x_mid,b_mid,bNORM,rNORM,pNORM)

          ! Calculate the difference between 1st and 2nd order integration
          ! and take ratio relative to dx_opt
          dx_rel=abs(dl)*maxval(abs(b_mid-b_ini))/dx_opt

          ! Make sure that dl does not change more than a factor of 2 or 0.5
          dx_rel=max(0.5,min(2.,dx_rel))

          if(dx_rel>1.)then
             ! Not accurate enough, decrease dl if possible

             if(abs(dl)<=dl_min+dl_tiny)then
                ! Cannot reduce dl further
                dl_next=dl
                EXIT
             end if

             dl = sign(max(dl_min,abs(dl)/(dx_rel+0.001)),dl)

          else
             ! Too accurate, increase dl if possible
             if(abs(dl)<dl_max-dl_tiny)then
                dl_next = sign(min(dl_max,abs(dl)/sqrt(dx_rel)),dl)
             end if

             EXIT
          end if
       end do

       x=x_ini+b_mid*dl
       ! dl/|B| for a cubic cell (dx=dy=dz)
       ! amount2add = abs( dl * dx_BLK(iBLK))/bNORM

       ! dl/|B| for a cell with arbitrary aspect ratio
       amount2add = abs(dl) * sqrt( &
            (b_mid(1)*dx_BLK(iBLK))**2 + &
            (b_mid(2)*dy_BLK(iBLK))**2 + &
            (b_mid(3)*dz_BLK(iBLK))**2) / bNORM
       local_fvol=local_fvol+amount2add
       local_rvol=local_rvol+amount2add*rNORM
       local_pvol=local_pvol+amount2add*pNORM

       if(dbg)then
          write(*,*)'  take step:',dx_BLK(iBLK),dy_BLK(iBLK),dz_BLK(iBLK),dl
          write(*,*)'    b_mid=',b_mid,' bNORM=',bNORM,' b_ini=',b_ini
          write(*,*)'  take step and add:',amount2add,rNORM,pNORM
       end if

       nsegment=nsegment+1
       l=l+abs(dl)

       if(any(x<xmin) .or. any(x>xmax))then

          if(dbg)write(*,*)'  stepped out of box'

          ! Ray points outwards from surface
          if(nsegment==1)then
             if(dbg)write(*,*)'  zeroing volumes, only one segment completed.'

             local_fvol=0.; local_rvol=0.; local_pvol=0.
             return
          end if

          ! Hit the wall, backup so that x is almost exactly on the wall
          ! just a little bit outside
          dl_back = dl*maxval(max(xmin-x,x-xmax)/(abs(x-x_ini)+dl_tiny))
          x=x-dl_back*b_mid
          ! dl/|B| for arbitrary aspect ratio
          amount2add = abs(dl_back) * sqrt( &
               (dx_BLK(iBLK)*b_mid(1))**2 + &
               (dy_BLK(iBLK)*b_mid(2))**2 + &
               (dz_BLK(iBLK)*b_mid(3))**2) / bNORM
          local_fvol=local_fvol-amount2add
          local_rvol=local_rvol-amount2add*rNORM
          local_pvol=local_pvol-amount2add*pNORM

          if(dbg)write(*,*)'  adjust step and subtract:',amount2add,rNORM,pNORM

          EXIT
       end if

       ! Check if we have integrated for too long
       if(l>lmax)then
          if(dbg)write(*,*)'  too many segments, zero values and stop.'

          ! Seems to be a closed loop within a block
          local_fvol=0.; local_rvol=0.; local_pvol=0.
          EXIT
       end if

    end do

  end subroutine do_integration

  subroutine interpolate_bbN(qx,qb,qbD,qrD,qpD)

    ! Obtain normalized bb field at normalized location qx and put it into qb
    use ModAdvance,ONLY:nVar
    real, intent(in) :: qx(3)
    real, intent(out):: qb(3),qbD,qrD,qpD
    real,dimension(nVar)::Aux_V
    !-------------------------------------------------------------------------

    ! Get B0 values for location

    xx(1)=x_BLK(1,1,1,iBLK)+dx_BLK(iBLK)*(qx(1)-1.)
    xx(2)=y_BLK(1,1,1,iBLK)+dy_BLK(iBLK)*(qx(2)-1.)
    xx(3)=z_BLK(1,1,1,iBLK)+dz_BLK(iBLK)*(qx(3)-1.)

    call get_b0(xx(1),xx(2),xx(3),qb)

    ! Determine cell indices corresponding to location qx

    i1=floor(qx(1)); i2=i1+1
    j1=floor(qx(2)); j2=j1+1
    k1=floor(qx(3)); k2=k1+1

    if(i1<-1.or.i2>nI+2.or.j1<-1.or.j2>nJ+2.or.k1<-1.or.k2>nK+2)then
       write(*,*)'interpolate_bbN: iProc, i1,j1,k1=',iProc,i1,j1,k1
       write(*,*)'interpolate_bbN: iProc, iBLK, qx=',iProc,iBLK, qx
       write(*,*)'interpolate_bbN: iProc, x,y,z(1,1,1),dx=',&
            x_BLK(1,1,1,iBLK),y_BLK(1,1,1,iBLK),z_BLK(1,1,1,iBLK),&
            dx_BLK(iBLK)
       write(*,*)'interpolate_bbN: iProc, x_0,y_0,z_0=',x_0,y_0,z_0
       call stop_mpi('ERROR in interpolate_bbN: location out of bounds')
    endif

    ! Distance relative to the cell centers

    dx1=qx(1)-i1; dx2=1.-dx1
    dy1=qx(2)-j1; dy2=1.-dy1
    dz1=qx(3)-k1; dz2=1.-dz1

    ! Add in interpolated B1 values
    Aux_V=interpolate_bb_v(nVar,State_VGB)
    qb=qb+Aux_V(Bx_:Bz_)

    ! Get density and pressure

    qrD = Aux_V(rho_)
    qpD = Aux_V(P_)
    qbD = sqrt(qb(1)**2 + qb(2)**2 + qb(3)**2)

    ! Normalize
    if(qbD>smallB)then
       ! Take aspect ratio of cells into account
       qb(1)=qb(1)/dx_BLK(iBLK)
       qb(2)=qb(2)/dy_BLK(iBLK)
       qb(3)=qb(3)/dz_BLK(iBLK)
       qb=qb/sqrt(sum(qb**2))
    else
       qb=0.
    end if

  end subroutine interpolate_bbN
  !==========================================================================
  function interpolate_bb_v(nVar,qbb)
    integer,intent(in)::nVar
    real, dimension(nVar,-1:nI+2,-1:nJ+2,-1:nK+2,nBLK), &
         intent(in):: qbb
    real, dimension(nVar)::interpolate_bb_v

    !-------------------------------------------------------------------------

    ! Bilinear interpolation in 3D

    interpolate_bb_v=&
         dx1*(   dy1*(   dz1*qbb(:,i2,j2,k2,iBLK)+&
         dz2*qbb(:,i2,j2,k1,iBLK))+&
         dy2*(   dz1*qbb(:,i2,j1,k2,iBLK)+&
         dz2*qbb(:,i2,j1,k1,iBLK)))+&
         dx2*(   dy1*(   dz1*qbb(:,i1,j2,k2,iBLK)+&
         dz2*qbb(:,i1,j2,k1,iBLK))+&
         dy2*(   dz1*qbb(:,i1,j1,k2,iBLK)+&
         dz2*qbb(:,i1,j1,k1,iBLK)))

  end function interpolate_bb_v

end subroutine integrate_ray
