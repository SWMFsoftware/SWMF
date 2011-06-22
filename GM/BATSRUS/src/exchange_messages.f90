!^CFG COPYRIGHT UM
subroutine exchange_messages
  use ModProcMH
  use ModMain, ONLY : nI, nJ, nK, gcn, nBlockMax, unusedBLK, &
       TypeBc_I, time_loop, UseB, &
       UseConstrainB,&              !^CFG IF CONSTRAINB 
       UseProjection,&              !^CFG IF PROJECTION
       UseDivbDiffusion,&           !^CFG IF DIVBDIFFUSE
       time_simulation,nOrder,prolong_order,optimize_message_pass, UseBatl
  use ModVarIndexes
  use ModAdvance, ONLY : State_VGB,divB1_GB
  use ModMessagePass, ONLY: message_pass_dir
  use ModParallel, ONLY : UsePlotMessageOptions
  use ModGeometry, ONLY : far_field_BCs_BLK        
  use ModMPCells, ONLY : DoOneCoarserLayer
  use ModBoundaryCells,ONLY:SaveBoundaryCells
  use ModPhysics, ONLY : ShockSlope
  use ModFaceValue,ONLY: UseAccurateResChange
  use ModEnergy,   ONLY: calc_energy_ghost

  use BATL_lib, ONLY: message_pass_cell
  use ModBatlInterface, ONLY: UseBatlTest

  use ModMpi

  implicit none

  integer :: iBlock
  logical :: oktest, oktest_me, oktime, oktime_me
  logical :: DoRestrictFace, DoOneLayer, DoTwoCoarseLayers
  logical :: DoCorners, DoFaces

  integer :: nWidth, nCoarseLayer

  !---------------------------------------------------------------------------

!!$  !^CFG IF DEBUGGING BEGIN
!!$  call testmessage_pass_nodes
!!$  call time_message_passing
!!$  !^CFG END DEBUGGING

  DoRestrictFace = prolong_order==1
  if(UseConstrainB) DoRestrictFace = .false.   !^CFG IF CONSTRAINB

  DoTwoCoarseLayers = &
       nOrder==2 .and. prolong_order==1 .and. .not. DoOneCoarserLayer

  call set_oktest('exchange_messages',oktest,oktest_me)
  call set_oktest('time_exchange',oktime,oktime_me)

  call timing_start('exch_msgs')
  ! Ensure that energy and pressure are consistent and positive in real cells
  !if(prolong_order==2)then     !^CFG IF NOT PROJECTION
  do iBlock = 1, nBlockMax
     if (unusedBLK(iBlock)) CYCLE
     if (far_field_BCs_BLK(iBlock).and.prolong_order==2)&
          call set_outer_BCs(iBlock,time_simulation,.false.)        
     if(UseConstrainB)call correctP(iBlock)   !^CFG IF CONSTRAINB
     if(UseProjection)call correctP(iBlock)   !^CFG IF PROJECTION
  end do
  !end if                       !^CFG IF NOT PROJECTION
  if(oktest)write(*,*)'Checked negative P, me=',iProc

  if (UsePlotMessageOptions) then
     if(UseBatl)then
        if(UseBatlTest)then
           call message_pass_cell(nVar, State_VGB, nProlongOrderIn=1)
        else
           call message_pass_cell(nVar, State_VGB)
        end if
        call fix_boundary_ghost_cells(DoRestrictFace)
     else
        if(oktest)write(*,*)'calling message_pass with plot options'
        !                              Don't send just one layer
        !                                      Don't send faces only
        !                                               Don't monotone restrict
        call message_pass_cells8(.false.,.false.,.false.,nVar, State_VGB)
        if(UseB) call message_pass_cells(.false.,.false.,.false.,DivB1_GB)
        if(SaveBoundaryCells)call fix_boundary_ghost_cells(DoRestrictFace)
     end if
  elseif (optimize_message_pass=='all') then
     if(oktest)write(*,*)'calling message_pass with corners: me,type=',&
          iProc,optimize_message_pass
     ! If ShockSlope is not zero then even the first order scheme needs 
     ! two ghost cell layers to fill in the corner cells at the sheared BCs.
     DoOneLayer = nOrder == 1 .and. ShockSlope == 0.0

     if(UseBatl)then
        nWidth = 2;       if(DoOneLayer)        nWidth = 1
        nCoarseLayer = 1; if(DoTwoCoarseLayers) nCoarseLayer = 2
        call message_pass_cell(nVar, State_VGB, &
             nWidthIn=nWidth, nProlongOrderIn=1, &
             nCoarseLayerIn=nCoarseLayer, DoRestrictFaceIn = DoRestrictFace)
        call fix_boundary_ghost_cells(DoRestrictFace)
     else
        call message_pass_cells8(DoOneLayer, .false., DoRestrictFace, &
             nVar, State_VGB)
        if(SaveBoundaryCells)call fix_boundary_ghost_cells(DoRestrictFace)
     end if
  else
     if(oktest)write(*,*)'calling message_pass: me,type=',&
          iProc,optimize_message_pass

     select case(optimize_message_pass)
     case('max','dir','face','min')
        ! Pass corners
        call message_pass_dir(1,3,2,.true.,prolong_order,nVar,&
             Sol_VGB=State_VGB, DoTwoCoarseLayers=DoTwoCoarseLayers, &
             restrictface=DoRestrictFace)
     case('opt')
        ! Do not pass corners if not necessary
        DoCorners = nOrder == 2 .and. UseAccurateResChange
        call message_pass_dir(1,3,nOrder,DoCorners,prolong_order,nVar,&
             Sol_VGB=State_VGB, DoTwoCoarseLayers=DoTwoCoarseLayers, &
             restrictface=DoRestrictFace)
     case('allopt')
        ! Do not pass corners if not necessary
        DoFaces = .not.(nOrder == 2 .and. UseAccurateResChange)
        ! Pass one layer if possible
        DoOneLayer = nOrder == 1
        if(UseBatl)then
           nWidth = 2;       if(DoOneLayer)        nWidth = 1
           nCoarseLayer = 1; if(DoTwoCoarseLayers) nCoarseLayer = 2
           call message_pass_cell(nVar, State_VGB, &
                nWidthIn=nWidth, nProlongOrderIn=1, nCoarseLayerIn=nCoarseLayer,&
                DoSendCornerIn=.not.DoFaces, DoRestrictFaceIn=DoRestrictFace)
           call fix_boundary_ghost_cells(DoRestrictFace)
        else
           call message_pass_cells8(DoOneLayer,DoFaces,DoRestrictFace, &
                nVar, State_VGB)
           if(SaveBoundaryCells)call fix_boundary_ghost_cells(DoRestrictFace)
        end if
     case default
        call stop_mpi('Unknown optimize_message_pass='//optimize_message_pass)
     end select
  end if

  if(oktest)write(*,*)'Ensure that E and P consistent, me=',iProc

  do iBlock = 1, nBlockMax
     if (unusedBLK(iBlock)) CYCLE
     if (far_field_BCs_BLK(iBlock)) &                        
          call set_outer_BCs(iBlock,time_simulation,.false.) 
     if(time_loop.and. any(TypeBc_I=='buffergrid'))&
          call fill_in_from_buffer(iBlock)
     call calc_energy_ghost(iBlock)
  end do

  call timing_stop('exch_msgs')
  if(oktime)call timing_show('exch_msgs',1)

  if(oktest)write(*,*)'exchange_messages finished, me=',iProc

end subroutine exchange_messages

! ADJOINT SPECIFIC BEGIN
!============================================================================
subroutine exchange_messages_adjoint

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,gcn,nBlockMax,unusedBLK, &
       TypeBc_I,time_loop,&
       time_simulation,nOrder,prolong_order,optimize_message_pass, UseBatl
  use ModVarIndexes
  use ModAdvance, ONLY : State_VGB, divB1_GB
  use ModMessagePass, ONLY: message_pass_dir
  use ModParallel, ONLY : UsePlotMessageOptions
  use ModGeometry, ONLY : far_field_BCs_BLK        
  use ModMpi
  use ModMPCells, ONLY : DoOneCoarserLayer
  use ModBoundaryCells,ONLY:SaveBoundaryCells
  use ModPhysics, ONLY : ShockSlope
  use ModFaceValue,ONLY: UseAccurateResChange
  use ModEnergy,   ONLY: calc_energy_ghost_adjoint
  use BATL_lib, ONLY: message_pass_cell
  use ModAdjoint, ONLY : Adjoint_VGB

  implicit none

  integer :: iBlock
  logical :: oktest, oktest_me, oktime, oktime_me
  logical :: DoRestrictFace, DoOneLayer, DoTwoCoarseLayers
  logical :: DoCorners, DoFaces
  integer :: nWidth, nCoarseLayer

  character (len=*), parameter :: NameSub = 'exchange_messages_adjoint'
  !-------------------------------------------------------------------------
 
  do iBlock = 1, nBlockMax
     if (unusedBLK(iBlock)) CYCLE
     call calc_energy_ghost_adjoint(iBlock)
     if (far_field_BCs_BLK(iBlock)) &                        
          call set_outer_BCs_adjoint (iBlock,time_simulation,.false.) 
     if(time_loop.and. any(TypeBc_I=='buffergrid'))&
          call stop_mpi(NameSub // ' Not supported!')
  end do
  
  ! Need to modify to send ghost info to procs (reverse mode) in an additive way
  call stop_mpi(NameSub // ' needs reverse ghost message passing!')

end subroutine exchange_messages_adjoint
! ADJOINT SPECIFIC END



!============================================================================!
subroutine fill_in_from_buffer(iBlock)
  use ModGeometry,ONLY:R_BLK,x_BLK,y_BLK,z_BLK
  use ModMain,ONLY:nI,nJ,nK,gcn,rBuffMin,rBuffMax,nDim,x_,y_,z_
  use ModAdvance,ONLY:nVar,State_VGB,rho_,rhoUx_,rhoUz_,Ux_,Uz_
  use ModProcMH,ONLY:iProc
  implicit none
  integer,intent(in)::iBlock
  integer::i,j,k
  real::X_D(nDim)
  logical::DoWrite=.true.
  if(DoWrite)then
     DoWrite=.false.
     if(iProc==0)then
        write(*,*)'Fill in the cells near the inner boundary from the buffer'
     end if
  end if
  
  do k=1-gcn,nK+gcn
     do j=1-gcn,nJ+gcn
        do i=1-gcn,nI+gcn
           if(R_BLK(i,j,k,iBlock)>rBuffMax.or.R_BLK(i,j,k,iBlock)<rBuffMin)&
                CYCLE
           X_D(x_)=x_BLK(i,j,k,iBlock)
           X_D(y_)=y_BLK(i,j,k,iBlock)
           X_D(z_)=z_BLK(i,j,k,iBlock)
           !Get interpolated values from buffer grid:
           call get_from_spher_buffer_grid(&
                X_D,nVar,State_VGB(:,i,j,k,iBlock))
           !Transform primitive variables to conservative ones:
           State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock)=&
                State_VGB(Ux_:Uz_,i,j,k,iBlock)*&
                State_VGB(rho_   ,i,j,k,iBlock)
        end do
     end do
  end do
end subroutine fill_in_from_buffer
!^CFG IF DEBUGGING BEGIN
!============================================================================
! Test timing of various message passing options
subroutine time_message_passing
  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,gcn,nBlockMax,unusedBLK, &
       UseConstrainB,&              !^CFG IF CONSTRAINB 
       UseProjection,&              !^CFG IF PROJECTION
       time_simulation,nOrder,prolong_order,optimize_message_pass
  use ModVarIndexes
  use ModAdvance, ONLY : &
       State_VGB
  use ModMessagePass, ONLY: message_pass_dir
  use ModMpi
  implicit none

  integer :: iError
  real*8 :: time_this
  logical :: DoOneLayer, DoRestrictFace

  !---------------------------------------------------------------------------

  ! For first order, message pass cells can pass only one layer of ghost cells.
  DoOneLayer = nOrder==1

  DoRestrictFace = prolong_order==1
  if(UseConstrainB) DoRestrictFace = .false.   !^CFG IF CONSTRAINB

  if(iProc==0) &
       write(*,*)' Timing message passing options ...', &
       ' nOrder=',nOrder,' DoOneLayer=',DoOneLayer

!!!
  call message_pass_dir(1,3,2,.true.,prolong_order,nVar,&
       Sol_VGB=State_VGB, restrictface=DoRestrictFace)

  call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
  time_this=MPI_WTIME()
  call message_pass_dir(1,3,2,.true.,prolong_order,nVar,&
       Sol_VGB=State_VGB, restrictface=DoRestrictFace)

  call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
  if(iProc==0) &
       write(*,'(a,f8.5,a)')' dir-1,3,2,T  took',MPI_WTIME()-time_this,' sec'

!!!
  call message_pass_dir(1,3,nORDER,.false.,prolong_order,nVar,&
       Sol_VGB=State_VGB,restrictface=DoRestrictFace)

  call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
  time_this=MPI_WTIME()
  call message_pass_dir(1,3,nORDER,.false.,prolong_order,nVar,&
       Sol_VGB=State_VGB,restrictface=DoRestrictFace)

  call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
  if(iProc==0) &
       write(*,'(a,f8.5,a)')' dir-1,3,nORDER,F  took',MPI_WTIME()-time_this,' sec'

!!!
  call testmessage_pass_cells

!!!
  call message_pass_cells8(DoOneLayer, .false., DoRestrictFace,nVar, State_VGB)

  call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
  time_this=MPI_WTIME()
  call message_pass_cells8(DoOneLayer, .false., DoRestrictFace,nVar, State_VGB)

  call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
  if(iProc==0) &
       write(*,'(a,f8.5,a)')' 8state-DoOneLayer,F,T  took',MPI_WTIME()-time_this,' sec'

!!!
  call message_pass_cells8(DoOneLayer, .true., DoRestrictFace,nVar, State_VGB)

  call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------  
  time_this=MPI_WTIME()  
  call message_pass_cells8(DoOneLayer, .true., DoRestrictFace,nVar, State_VGB)

  call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
  if(iProc==0) &
       write(*,'(a,f8.5,a)')' 8state-DoOneLayer,T,T  took',MPI_WTIME()-time_this,' sec'

!!!
  call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
  call MPI_Finalize(iError)
  stop

end subroutine time_message_passing
!^CFG END DEBUGGING

!^CFG IF PROJECTION BEGIN
!============================================================================
subroutine correctP(iBlock)

  ! Make pressure and energy consistent and maintain thermal energy ratio 
  ! at a reasonable value (this also excludes negative pressure)

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,Itest,Jtest,Ktest,BLKtest
  use ModVarIndexes,ONLY:&
       rho_,rhoUx_,rhoUy_,rhoUz_,Bx_,By_,Bz_,P_,nVar
  use ModAdvance, ONLY : State_VGB, Energy_GBI
  use ModPhysics, ONLY : gm1, inv_gm1, Pratio_hi, Pratio_lo
  use ModGeometry, ONLY : x_BLK, y_BLK, z_BLK, true_cell
  implicit none

  integer, intent(in) :: iBlock

  integer :: i,j,k
  real :: inv_dratio, qp, qe, qth, qratio, qd, qde, qpmin, &
       qdesum, qdesumabs, qderelmax

  real, dimension(1:nI,1:nJ,1:nK) :: P_old

  integer :: ierror1=-1, ierror2=-1, ierror3=-1, loc(3)

  logical :: oktest, oktest_me
  !--------------------------------------------------------------------------

  if(iBlock==BLKtest)then
     call set_oktest('correctP',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if

  qpmin=1.
  qdesum=0.
  qdesumabs=0.
  qderelmax=0.

  P_old=State_VGB(P_,1:nI,1:nJ,1:nK,iBlock)

  inv_dratio=1./(Pratio_hi-Pratio_lo)

  do k=1,nK; do j=1,nJ; do i=1,nI

     if(.not.true_cell(i,j,k,iBlock))CYCLE

     ! Pressure and total energy
     qp=P_old(i,j,k)
     qe=Energy_GBI(i,j,k,iBlock,1)

     if(oktest_me.and.i==Itest.and.J==Jtest.and.K==Ktest)&
          write(*,*)'CorrectP at me,BLK,i,j,k=',&
          iProc,BLKtest,Itest,Jtest,Ktest, &
          ', initial P,E=',qp,qe

     ! Memorize smallest pressure
     qpmin=min(qp,qpmin)

     ! Thermal energy
     qth=inv_gm1*qp

     ! Deviation=extra total energy=qe-inv_gm1*qp-(rhoU**2/rho+B**2)/2
     qd=qE-qth                                                         &
          -0.5*(State_VGB(rhoUx_,i,j,k,iBlock)**2+                         &
          State_VGB(rhoUy_,i,j,k,iBlock)**2+                               &
          State_VGB(rhoUz_,i,j,k,iBlock)**2)/State_VGB(rho_,i,j,k,iBlock)      &
          -0.5*(State_VGB(Bx_,i,j,k,iBlock)**2+                            &
          State_VGB(By_,i,j,k,iBlock)**2+                                  &
          State_VGB(Bz_,i,j,k,iBlock)**2)

     ! Limited thermal/total energy ratio for correction
     qratio=min(Pratio_hi,max(Pratio_lo,min(qth,qth+qd)/qe))

     ! Total energy is modified by qde (=0 if qratio==Pratio_hi)
     qde=qd*(Pratio_hi-qratio)*inv_dratio

     ! Collect total energy change
     qdesum   =qdesum   -qde
     qdesumabs=qdesumabs+abs(qde)
     qderelmax=max(qderelmax,qde/qe)

     ! Pressure is modified
     State_VGB(P_,i,j,k,iBlock)=gm1*(qth+qd-qde)

     ! We should now have E=inv_gm1*P+(rhoU**2/rho+B**2)/2:
     !
     ! qp*inv_gm1+qd-qde + (rhoU**2/rho+B**2)/2 = qe-qde = E
     !
     ! Correct!

     if(oktest_me.and.i==Itest.and.J==Jtest.and.K==Ktest)then
        write(*,*)'qp,qth,qe,qd,qratio,qde=',qp,qth,qe,qd,qratio,qde
        write(*,*)'CorrectP, final P=',State_VGB(P_,i,j,k,iBlock)
     end if

  end do; end do; end do

  if(qpmin<0.)then
     if(ierror1==-1)then
        loc=minloc(P_old)
        write(*,*)'Negative P at me,iBLK,I,J,K,x,y,z,val',&
             iProc,iBlock,loc,&
             x_BLK(loc(1),loc(2),loc(3),iBlock),&
             y_BLK(loc(1),loc(2),loc(3),iBlock),&
             z_BLK(loc(1),loc(2),loc(3),iBlock),&
             P_old(loc(1),loc(2),loc(3))
     end if
     call error_report('Negative P in exchange msgs, min(P)', &
          qpmin,ierror1,.true.)
  end if
  if(qderelmax>1.0E-3)then
     call error_report('E change in exchange_msgs, dE',qdesum,ierror2,.false.)
     call error_report('|E| change in exchange_msgs, d|E|',qdesumabs,ierror3,&
          .false.)
  end if

  if(oktest_me)write(*,*)'CorrectP qpmin=',qpmin

end subroutine correctP
!^CFG END PROJECTION

