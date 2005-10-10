!^CFG COPYRIGHT UM
subroutine exchange_messages
  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,gcn,globalBLK,nBlockMax,unusedBLK, &
       UseConstrainB,&              !^CFG IF CONSTRAINB 
       UseProjection,&              !^CFG IF PROJECTION
       UseDivbDiffusion,&           !^CFG IF DIVBDIFFUSE
       time_simulation,nOrder,prolong_order,optimize_message_pass
  use ModVarIndexes
  use ModAdvance, ONLY : &
       State_VGB
  use ModInterface
  use ModParallel, ONLY : UsePlotMessageOptions
  use ModGeometry, ONLY : far_field_BCs_BLK            
  use ModMpi
  use ModMPCells, ONLY : DoOneCoarserLayer
  use ModBoundaryCells,ONLY:SaveBoundaryCells
  implicit none

  logical :: oktest, oktest_me, oktime, oktime_me
  logical :: DoRestrictFace, DoOneLayer, DoTwoCoarseLayers

  !---------------------------------------------------------------------------

!!$  !^CFG IF DEBUGGING BEGIN
!!$  call testmessage_pass_nodes
!!$  call time_message_passing
!!$  !^CFG END DEBUGGING

  ! For first order, message pass cells can pass only one layer of ghost cells.
  DoOneLayer = nOrder==1
  if(UseDivbDiffusion)DoOneLayer=.false.     !^CFG IF DIVBDIFFUSE

  DoRestrictFace = prolong_order==1
  if(UseConstrainB) DoRestrictFace = .false.   !^CFG IF CONSTRAINB

  DoTwoCoarseLayers = &
       nOrder==2 .and. prolong_order==1 .and. .not. DoOneCoarserLayer

  call set_oktest('exchange_messages',oktest,oktest_me)
  call set_oktest('time_exchange',oktime,oktime_me)

  call timing_start('exch_msgs')
  ! Ensure that energy and pressure are consistent and positive in real cells
  !if(prolong_order==2)then     !^CFG IF NOT PROJECTION
  do globalBLK = 1, nBlockMax
     if (unusedBLK(globalBLK)) CYCLE
     if (far_field_BCs_BLK(globalBLK).and.prolong_order==2)&
          call set_outer_BCs(globalBLK,time_simulation,.false.)        
     if(UseConstrainB)call correctP   !^CFG IF CONSTRAINB
     if(UseProjection)call correctP   !^CFG IF PROJECTION
  end do
  !end if                       !^CFG IF NOT PROJECTION
  if(oktest)write(*,*)'Checked negative P, me=',iProc

  if (UsePlotMessageOptions) then
     if(oktest)write(*,*)'calling message_pass with plot options'

     ! Send all faces, edges, and corners
     ! Don't do the monotone restriction
     ! Don't send just one layer
     call message_pass_cells_8state(.false.,.false.,.false.)

  elseif (optimize_message_pass=='all') then
     if(oktest)write(*,*)'calling message_pass with corners: me,type=',&
          iProc,optimize_message_pass
     call message_pass_cells_8state(DoOneLayer,.false.,DoRestrictFace)
     if(SaveBoundaryCells)call fix_boundary_ghost_cells(DoRestrictFace)
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
        ! Do not pass corners
        call message_pass_dir(1,3,nORDER,.false.,prolong_order,nVar,&
             Sol_VGB=State_VGB, DoTwoCoarseLayers=DoTwoCoarseLayers, &
             restrictface=DoRestrictFace)
     case('allopt')
        call message_pass_cells_8state(DoOneLayer,.true.,DoRestrictFace)
        if(SaveBoundaryCells)call fix_boundary_ghost_cells(DoRestrictFace)
      case default
        call stop_mpi('Unknown optimize_message_pass='//optimize_message_pass)
     end select
  end if

  if(oktest)write(*,*)'Ensure that E and P consistent, me=',iProc

  do globalBLK = 1, nBlockMax
     if (unusedBLK(globalBLK)) CYCLE
     if (far_field_BCs_BLK(globalBLK)) &                        
          call set_outer_BCs(globalBLK,time_simulation,.false.) 
     call correctE
  end do

  call timing_stop('exch_msgs')
  if(oktime)call timing_show('exch_msgs',1)

  if(oktest)write(*,*)'exchange_messages finished, me=',iProc

end subroutine exchange_messages

!^CFG IF DEBUGGING BEGIN
!============================================================================
! Test timing of various message passing options
subroutine time_message_passing
  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,gcn,globalBLK,nBlockMax,unusedBLK, &
       UseConstrainB,&              !^CFG IF CONSTRAINB 
       UseProjection,&              !^CFG IF PROJECTION
       UseDivbDiffusion,&           !^CFG IF DIVBDIFFUSE
       time_simulation,nOrder,prolong_order,optimize_message_pass
  use ModVarIndexes
  use ModAdvance, ONLY : &
       State_VGB
  use ModInterface
  use ModMpi
  implicit none

  integer :: iError
  real*8 :: time_this
  logical :: DoOneLayer, DoRestrictFace

  !---------------------------------------------------------------------------

  ! For first order, message pass cells can pass only one layer of ghost cells.
  DoOneLayer = nOrder==1
  if(UseDivbDiffusion)DoOneLayer=.false.     !^CFG IF DIVBDIFFUSE

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
  call message_pass_cells_8state(DoOneLayer, .false., DoRestrictFace)

  call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
  time_this=MPI_WTIME()
  call message_pass_cells_8state(DoOneLayer, .false., DoRestrictFace)

  call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
  if(iProc==0) &
       write(*,'(a,f8.5,a)')' 8state-DoOneLayer,F,T  took',MPI_WTIME()-time_this,' sec'

!!!
  call message_pass_cells_8state(DoOneLayer, .true., DoRestrictFace)

  call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------  
  time_this=MPI_WTIME()  
  call message_pass_cells_8state(DoOneLayer, .true., DoRestrictFace)

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
subroutine correctP

  ! Make pressure and energy consistent and maintain thermal energy ratio 
  ! at a reasonable value (this also excludes negative pressure)

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,Itest,Jtest,Ktest,BLKtest,globalBLK
  use ModVarIndexes,ONLY:&
       rho_,rhoUx_,rhoUy_,rhoUz_,Bx_,By_,Bz_,P_,nVar
  use ModAdvance, ONLY : &
       State_VGB, E_BLK
  use ModPhysics, ONLY : gm1, inv_gm1, Pratio_hi, Pratio_lo
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,true_cell
  implicit none

  integer :: i,j,k
  real :: inv_dratio, qp, qe, qth, qratio, qd, qde, qpmin, &
       qdesum, qdesumabs, qderelmax

  real, dimension(1:nI,1:nJ,1:nK) :: P_old

  integer :: ierror1=-1, ierror2=-1, ierror3=-1, loc(3)

  logical :: oktest, oktest_me
  !--------------------------------------------------------------------------

  if(globalBLK==BLKtest)then
     call set_oktest('correctP',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if

  qpmin=1.
  qdesum=0.
  qdesumabs=0.
  qderelmax=0.

  P_old=State_VGB(P_,1:nI,1:nJ,1:nK,globalBLK)

  inv_dratio=1./(Pratio_hi-Pratio_lo)

  do k=1,nK; do j=1,nJ; do i=1,nI

     if(.not.true_cell(i,j,k,globalBLK))CYCLE

     ! Pressure and total energy
     qp=P_old(i,j,k)
     qe=E_BLK(i,j,k,globalBLK)

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
          -0.5*(State_VGB(rhoUx_,i,j,k,globalBLK)**2+                         &
          State_VGB(rhoUy_,i,j,k,globalBLK)**2+                               &
          State_VGB(rhoUz_,i,j,k,globalBLK)**2)/State_VGB(rho_,i,j,k,globalBLK)      &
          -0.5*(State_VGB(Bx_,i,j,k,globalBLK)**2+                            &
          State_VGB(By_,i,j,k,globalBLK)**2+                                  &
          State_VGB(Bz_,i,j,k,globalBLK)**2)

     ! Limited thermal/total energy ratio for correction
     qratio=min(Pratio_hi,max(Pratio_lo,min(qth,qth+qd)/qe))

     ! Total energy is modified by qde (=0 if qratio==Pratio_hi)
     qde=qd*(Pratio_hi-qratio)*inv_dratio

     ! Collect total energy change
     qdesum   =qdesum   -qde
     qdesumabs=qdesumabs+abs(qde)
     qderelmax=max(qderelmax,qde/qe)

     ! Pressure is modified
     State_VGB(P_,i,j,k,globalBLK)=gm1*(qth+qd-qde)

     ! We should now have E=inv_gm1*P+(rhoU**2/rho+B**2)/2:
     !
     ! qp*inv_gm1+qd-qde + (rhoU**2/rho+B**2)/2 = qe-qde = E
     !
     ! Correct!

     if(oktest_me.and.i==Itest.and.J==Jtest.and.K==Ktest)then
        write(*,*)'qp,qth,qe,qd,qratio,qde=',qp,qth,qe,qd,qratio,qde
        write(*,*)'CorrectP, final P=',State_VGB(P_,i,j,k,globalBLK)
     end if

  end do; end do; end do

  if(qpmin<0.)then
     if(ierror1==-1)then
        loc=minloc(P_old)
        write(*,*)'Negative P at me,iBLK,I,J,K,x,y,z,val',&
             iProc,globalBLK,loc,&
             x_BLK(loc(1),loc(2),loc(3),globalBLK),&
             y_BLK(loc(1),loc(2),loc(3),globalBLK),&
             z_BLK(loc(1),loc(2),loc(3),globalBLK),&
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

!==============================================================================
subroutine correctE

  ! Correct total energy
  use ModVarIndexes
  use ModMain, ONLY : globalBLK,ndim,&
       nI,nJ,nK, gcn
  use ModAdvance, ONLY : State_VGB,E_BLK
  use ModPhysics, ONLY : inv_gm1
  implicit none

  integer::i,j,k

  do k=1-gcn,nK+gcn; do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
     if(State_VGB(rho_,i,j,k,globalBLK)<=0.0)cycle
     E_BLK(i,j,k,globalBLK)=inv_gm1*State_VGB(P_,i,j,k,globalBLK) &
          +0.5*(sum(State_VGB(rhoUx_:rhoU_+ndim,i,j,k,globalBLK)**2)/&
          State_VGB(rho_,i,j,k,globalBLK)    &
          +sum(State_VGB(Bx_:B_+ndim,i,j,k,globalBLK)**2) )
  end do; end do; end do

end subroutine correctE
