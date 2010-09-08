!^CFG COPYRIGHT UM
subroutine update_states(iStage,iBlock)
  use ModProcMH
  use ModMain
  use ModAdvance
  use ModGeometry, ONLY : x_BLK, y_BLK, z_BLK
  use ModUser, ONLY: user_update_states
  use ModMultiFluid, ONLY: select_fluid, iFluid, nFluid, iP
  use ModAdjoint, ONLY : AdjPreUpdate_, store_block_buffer, DoAdjoint  ! ADJOINT SPECIFIC
  implicit none

  integer, intent(in) :: iStage,iBlock
  integer :: i,j,k, iVar

  logical :: oktest, oktest_me
  character(len=*), parameter:: NameSub = 'update_states'
  !--------------------------------------------------------------------------
  if(iBlock==BLKtest .and. iProc==PROCtest)then
     call set_oktest(NameSub,oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  endif

  if(oktest_me)then
     write(*,*)NameSub,' dt=',time_BLK(iTest, jTest, kTest, iBlock)
     if(allocated(IsConserv_CB)) write(*,*)NameSub,' IsConserv=', &
          IsConserv_CB(iTest,jTest,kTest,iBlock)
     write(*,*)
     do iVar=1,nVar
        write(*,*)NameVar_V(iVar), '(TestCell)=',&
             State_VGB(iVar,Itest,Jtest,Ktest,BLKtest)
     end do
     do iFluid = 1, nFluid
        write(*,'(a,i2,a,g19.12)') &
             'E(',iFluid,')=',Energy_GBI(Itest,Jtest,Ktest,BLKtest,iFluid)
     end do
     write(*,*)'Fluxes and sources for ',NameVar_V(VarTest)
     write(*,*)'X fluxes L,R=',Flux_VX(VARtest,Itest,Jtest,Ktest),&
          Flux_VX(VARtest,Itest+1,Jtest,Ktest)
     write(*,*)'Y fluxes L,R=',Flux_VY(VARtest,Itest,Jtest,Ktest),&
          Flux_VY(VARtest,Itest,Jtest+1,Ktest)
     write(*,*)'Z fluxes L,R=',Flux_VZ(VARtest,Itest,Jtest,Ktest),&
          Flux_VZ(VARtest,Itest,Jtest,Ktest+1)
     write(*,*)'source=',Source_VC(VARtest,Itest,Jtest,Ktest)
  end if

  if (DoAdjoint) call store_block_buffer(iBlock,AdjPreUpdate_)   ! ADJOINT SPECIFIC

  if(UseUserUpdateStates)then
     call user_update_states(iStage,iBlock)
  else
     call update_states_MHD(iStage,iBlock)
  end if

  !^CFG IF DEBUGGING BEGIN
  if(index(test_string,'fixrho ')>0) &
       State_VGB(Rho_,1:nI,1:nJ,1:nK,iBlock)=StateOld_VCB(Rho_,:,:,:,iBlock)

  if(index(test_string,'fixrho ')>0)   State_VGB(Rho_,1:nI,1:nJ,1:nK,iBlock)&
       =                            StateOld_VCB(Rho_,:,:,:,iBlock)
  if(index(test_string,'fixrhoux ')>0) State_VGB(RhoUx_,1:nI,1:nJ,1:nK,iBlock)&
       =                            StateOld_VCB(RhoUx_,:,:,:,iBlock)
  if(index(test_string,'fixrhouy ')>0) State_VGB(RhoUy_,1:nI,1:nJ,1:nK,iBlock)&
       =                            StateOld_VCB(RhoUy_,:,:,:,iBlock)
  if(index(test_string,'fixrhouz ')>0) State_VGB(RhoUz_,1:nI,1:nJ,1:nK,iBlock)&
       =                            StateOld_VCB(RhoUz_,:,:,:,iBlock)
  if(index(test_string,'fixbx ')>0)    State_VGB(Bx_,1:nI,1:nJ,1:nK,iBlock)&
       =                            StateOld_VCB(Bx_,:,:,:,iBlock)
  if(index(test_string,'fixby ')>0)    State_VGB(By_,1:nI,1:nJ,1:nK,iBlock)&
       =                            StateOld_VCB(By_,:,:,:,iBlock)
  if(index(test_string,'fixbz ')>0)    State_VGB(Bz_,1:nI,1:nJ,1:nK,iBlock)&
       =                            StateOld_VCB(Bz_,:,:,:,iBlock)
  if(index(test_string,'fixe ')>0 .or. index(test_string,'fixp ')>0) then
     Energy_GBI(1:nI,1:nJ,1:nK,iBlock,:)=EnergyOld_CBI(:,:,:,iBlock,:)
     do iFluid=1, nFluid
        call select_fluid
        State_VGB(iP,1:nI,1:nJ,1:nK,iBlock)=StateOld_VCB(iP,:,:,:,iBlock)
     end do
  endif

  if(index(test_string,'fixothers')>0) then
     if(oktest_me)write(*,*)'Fix others!!!'
     do k=1,nK; do j=1,nJ; do i=1,nI
        if(  abs(x_BLK(i,j,k,iBlock)-Xtest)+ &
             abs(y_BLK(i,j,k,iBlock)-Ytest)+ &
             abs(z_BLK(i,j,k,iBlock)-Ztest) < 1.1 ) CYCLE

        State_VGB(:,i,j,k,iBlock)  = StateOld_VCB(:,i,j,k,iBlock)
        Energy_GBI(i,j,k,iBlock,:) = EnergyOld_CBI(i,j,k,iBlock,:)
     end do; end do; end do
  end if
  !^CFG END DEBUGGING

  if(oktest_me)then
     write(*,*)NameSub,' final:'
     do iVar=1,nVar
        write(*,*)NameVar_V(iVar),'(TestCell)  =',&
             State_VGB(iVar,Itest,Jtest,Ktest,BLKtest)
     end do
     do iFluid = 1, nFluid
        write(*,'(a,i2,a,g19.12)') &
             'E(',iFluid,')=',Energy_GBI(Itest,Jtest,Ktest,BLKtest,iFluid)
     end do
  end if

end subroutine update_states

! BEGIN ADJOINT SPECIFIC
!============================================================================
subroutine update_states_adjoint(iStage,iBlock) 
  use ModMain
  use ModUser, ONLY: user_update_states_adjoint
  use ModAdjoint, ONLY: AdjPreUpdate_, recall_block_buffer
  implicit none

  integer, intent(in) :: iStage,iBlock

  character(len=*), parameter:: NameSub = 'update_states_adjoint'
  !--------------------------------------------------------------------------

  if(UseUserUpdateStates)then
     call user_update_states_adjoint(iStage,iBlock)
  else
     call update_states_MHD_adjoint(iStage,iBlock)
  end if

  ! Restore stage-beginning state values from buffer
  call recall_block_buffer(iBlock,AdjPreUpdate_) 

end subroutine update_states_adjoint
! ADJOINT SPECIFIC END



!============================================================================
subroutine update_check(iStage)

  ! Check updated values for allowed change in density or pressure

  use ModProcMH
  use ModMain
  use ModImplicit, ONLY: UsePartImplicit !^CFG IF IMPLICIT
  use ModAdvance
  use ModPhysics
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,R_BLK,true_cell
  use ModNumConst, ONLY: cTiny
  use ModMpi
  use ModEnergy
  use ModMultiFluid, ONLY: IsMhd
  use ModMultiIon,   ONLY: DoRestrictMultiIon, IsMultiIon_CB

  implicit none

  integer,intent(in) :: iStage

  integer, parameter :: max_checks=25
  integer :: i,j,k, iVar, iBlock,num_checks, iError
  real :: time_fraction_rho, min_time_fraction_rho
  real :: time_fraction_p,   min_time_fraction_p
  real :: time_fraction, cell_time_fraction, report_tf, report_tf_all
  real, dimension(2) :: percent_chg_rho
  real, dimension(2) :: percent_chg_p
  real, dimension(4) :: PercentChangePE, PercentChangeMax

  real :: Value,B0_DC(3,nI,nJ,nK)
  integer :: i_D(3)
  logical :: update_check_done, IsNegative, DoStop
  logical :: oktest,  oktest_me
  logical :: oktest1, oktest1_me
  logical :: oktest2, oktest2_me
  logical :: oktest3, oktest3_me

  real    :: RhoChangeMax_I(2), pChangeMax_I(2)
  logical :: DoTest, DoTestMe
  character(len=*), parameter:: NameSub = 'update_check'
  character(len=*), parameter:: format="(a,i5,a,f6.1,a,3es11.3)"

  integer :: iError1=-1
  !-----------------------------------------------------------------------

  call set_oktest('fix_update',DoTest,DoTestMe)

  call set_oktest(NameSub, oktest,oktest_me)
  call set_oktest('convergence_history',oktest1,oktest1_me)
  call set_oktest('update_check_detail',oktest2,oktest2_me)
  call set_oktest('locations',oktest3,oktest3_me)

  ! Check for allowable percentage changed from update
  if(time_accurate) then
     !\\\
     !    TIME ACCURATE  ===================================================
     !///
     report_tf = 1.
     do num_checks = 1,max_checks
        percent_chg_rho = 0.1
        percent_chg_p   = 0.1
        do iBlock = 1, nBlockMax
           if (unusedBLK(iBlock)) CYCLE
           if (num_checks == 1) then
              do k=1,nK; do j=1,nJ; do i=1,nI
                 do iVar = 1, nVar
                    if (DefaultState_V(iVar) <= cTiny) CYCLE

                    ! For sake of backward compatibility
                    if(iVar == P_) CYCLE

                    ! Do not check multi-ion variables in regions that are
                    ! not truely multi-ion
                    if(UseMultiIon .and. iVar>p_ .and. DoRestrictMultiIon)then
                       if(.not.IsMultiIon_CB(i,j,k,iBlock)) CYCLE
                    end if

                    ! Do not check minor species if not necessary
                    if(UseMultiSpecies .and. &
                         iVar >= SpeciesFirst_ .and. iVar <= SpeciesLast_ &
                         .and. StateOld_VCB(iVar,i,j,k,iBlock) &
                         < SpeciesPercentCheck*0.01*&
                         StateOld_VCB(Rho_,i,j,k,iBlock)) CYCLE

                    percent_chg_rho(1) = &
                         max(percent_chg_rho(1), 100*abs( min(0.,&
                         (State_VGB(iVar,i,j,k,iBlock)- &
                         StateOld_VCB(iVar,i,j,k,iBlock)) &
                         /StateOld_VCB(iVar,i,j,k,iBlock) ) ) )
                    percent_chg_rho(2) = &
                         max(percent_chg_rho(2), 100*abs( max(0.,&
                         (State_VGB(iVar,i,j,k,iBlock)- &
                         StateOld_VCB(iVar,i,j,k,iBlock)) &
                         /StateOld_VCB(iVar,i,j,k,iBlock) ) ) )
                 end do
              end do; end do; end do
           end if
           percent_chg_p(1) = &
                max(percent_chg_p(1), 100. * abs( min( 0., minval( &
                (State_VGB(P_,1:nI,1:nJ,1:nK,iBlock)- &
                StateOld_VCB(P_,1:nI,1:nJ,1:nK,iBlock)) &
                /StateOld_VCB(P_,1:nI,1:nJ,1:nK,iBlock) ) ) ) )
           percent_chg_p(2) = &
                max(percent_chg_p(2), 100. * abs( max( 0., maxval( &
                (State_VGB(P_,1:nI,1:nJ,1:nK,iBlock)- &
                StateOld_VCB(P_,1:nI,1:nJ,1:nK,iBlock)) &
                /StateOld_VCB(P_,1:nI,1:nJ,1:nK,iBlock) ) ) ) )
        end do
        if(oktest)then
           ! Find location of maximum change
           call MPI_allreduce(percent_chg_rho, RhoChangeMax_I, 2, &
                MPI_REAL, MPI_MAX, iComm, iError)
           call MPI_allreduce(percent_chg_p, pChangeMax_I, 2, &
                MPI_REAL, MPI_MAX, iComm, iError)

           do iBlock = 1, nBlockMax
              if (unusedBLK(iBlock)) CYCLE
              do k=1,nK; do j=1,nJ; do i=1,nI
                 do iVar = 1, nVar
                    if (DefaultState_V(iVar) <= cTiny) CYCLE

                    if(UseMultiSpecies .and. &
                         iVar >= SpeciesFirst_ .and. iVar <= SpeciesLast_ &
                         .and. StateOld_VCB(iVar,i,j,k,iBlock) &
                         < SpeciesPercentCheck*0.01*&
                         StateOld_VCB(Rho_,i,j,k,iBlock)) CYCLE
                    
                    if(iVar == p_)then
                       if(pChangeMax_I(1) > percent_max_p(1) .and. &
                            1e-4 > abs(pChangeMax_I(1) - 100. * abs( &
                            (   State_VGB(P_,i,j,k,iBlock)- &
                            StateOld_VCB (P_,i,j,k,iBlock)) &
                            /StateOld_VCB(P_,i,j,k,iBlock) ))) &
                            write(*,format)NameSub//' nStep=',n_step,&
                            ' max p drop=',pChangeMax_I(1),&
                            '% at x,y,z=',&
                            x_BLK(i,j,k,iBlock),&
                            y_BLK(i,j,k,iBlock),&
                            z_BLK(i,j,k,iBlock)


                       if(pChangeMax_I(2) > percent_max_p(2) .and. &
                            1e-4 > abs(pChangeMax_I(2) - 100. * abs( &
                            (   State_VGB(P_,i,j,k,iBlock)- &
                            StateOld_VCB (P_,i,j,k,iBlock)) &
                            /StateOld_VCB(P_,i,j,k,iBlock)  ))) &
                            write(*,format)NameSub//' nStep=',n_step,&
                            ' max p increase=',&
                            pChangeMax_I(2), &
                            '% at x,y,z=',&
                            x_BLK(i,j,k,iBlock),&
                            y_BLK(i,j,k,iBlock),&
                            z_BLK(i,j,k,iBlock)
                            
                       CYCLE
                    end if
                    
                    if(RhoChangeMax_I(1) > percent_max_rho(1) .and. &
                         1e-4 > abs(RhoChangeMax_I(1) - 100*abs( &
                         (State_VGB(iVar,i,j,k,iBlock)- &
                         StateOld_VCB(iVar,i,j,k,iBlock)) &
                         /StateOld_VCB(iVar,i,j,k,iBlock) ))) &
                         write(*,format)NameSub//' nStep=',n_step,&
                         ' max '//trim(NameVar_V(iVar))//' drop=', &
                         RhoChangeMax_I(1), &
                         '% at x,y,z=',&
                         x_BLK(i,j,k,iBlock),&
                         y_BLK(i,j,k,iBlock),&
                         z_BLK(i,j,k,iBlock)

                    if(RhoChangeMax_I(2) > percent_max_rho(2) .and. &
                         1e-4 > abs(RhoChangeMax_I(2) - 100*abs( &
                         (State_VGB(iVar,i,j,k,iBlock)- &
                         StateOld_VCB(iVar,i,j,k,iBlock)) &
                         /StateOld_VCB(iVar,i,j,k,iBlock) ))) &
                         write(*,format)NameSub//' nStep=',n_step,&
                         ' max '//trim(NameVar_V(iVar))//' increase=',&
                         RhoChangeMax_I(2), &
                         '% at x,y,z=',&
                         x_BLK(i,j,k,iBlock),&
                         y_BLK(i,j,k,iBlock),&
                         z_BLK(i,j,k,iBlock)
                 end do
              end do; end do; end do
           end do
        end if !oktest
        time_fraction_rho = 1.0 / maxval(percent_chg_rho/percent_max_rho)
        call MPI_allreduce(time_fraction_rho, min_time_fraction_rho, 1, &
             MPI_REAL, MPI_MIN, iComm, iError)
        time_fraction_p   = 1.0 / maxval(percent_chg_p  /percent_max_p  )
        call MPI_allreduce(time_fraction_p, min_time_fraction_p, 1, &
             MPI_REAL, MPI_MIN, iComm, iError)
        if (min_time_fraction_rho >= 1. .and. min_time_fraction_p >= 1.) EXIT

        if (num_checks == 1) then
           time_fraction = 1.
           if (min_time_fraction_rho < 1.) &
                time_fraction = 0.9*min_time_fraction_rho
           if (min_time_fraction_p   < 1.) &
                time_fraction = min(time_fraction, 0.75)
        else
           time_fraction = 0.5
        end if
        dt = dt * time_fraction
        report_tf = report_tf*time_fraction
        do iBlock = 1, nBlockMax
           if (unusedBLK(iBlock)) CYCLE
           if(UseB0)then
              B0_DC=B0_DGB(:,1:nI,1:nJ,1:nK,iBlock)
           else
              B0_DC=0.00
           end if
           ! Fix the update in the cells
           do k=1,nK; do j=1,nJ; do i=1,nI
              call fix_update
           end do; end do; end do
        end do
     end do
     PercentChangePE(1:2) =  percent_chg_rho(1:2) - 0.1
     PercentChangePE(3:4) =  percent_chg_p(1:2) - 0.1

     !^CFG IF IMPLICIT BEGIN
     ! The part implicit scheme can get here if all blocks become explicit
     ! due to time step reductions. To be able to recover the time step,
     ! increase fixed time step if there was no time step reduction above.
     if(UsePartImplicit .and. dt == DtFixed) &
          DtFixed = min(DtFixedOrig, DtFixed*1.05)
     !^CFG END IMPLICIT

     if(oktest) then
        if (iProc == 0 .and. report_tf < 1.) &
             write(*,'(a,a,i6,a,f12.8,a,f12.8)') 'update_check TA:', &
             ' nStep=',n_step,'     dt reduction=',report_tf,' dt=',dt
     end if
  else
     !\\\
     !    LOCAL TIMESTEPPING
     !///
     report_tf = 1.
     PercentChangePE = cZero
     do iBlock = 1, nBlockMax
        if (unusedBLK(iBlock)) CYCLE
        if(UseB0)then
           B0_DC=B0_DGB(:,1:nI,1:nJ,1:nK,iBlock)
        else
           B0_DC=0.00
        end if
        do k=1,nK; do j=1,nJ; do i=1,nI
           cell_time_fraction = 1.
           do num_checks = 1, max_checks
              update_check_done = .true.
              percent_chg_rho = 0.1
              percent_chg_p   = 0.1
              if (num_checks == 1) then
                 do iVar = 1, nVar
                    if (DefaultState_V(iVar) <= cTiny) CYCLE

                    ! This is for backwards compatibility
                    if(iVar == P_) CYCLE

                    if(UseMultiSpecies .and. &
                         iVar >= SpeciesFirst_ .and. iVar <= SpeciesLast_ &
                         .and. StateOld_VCB(iVar,i,j,k,iBlock) < &
                         SpeciesPercentCheck*0.01*&
                         StateOld_VCB(Rho_,i,j,k,iBlock)) CYCLE

                    percent_chg_rho(1) = max(percent_chg_rho(1), &
                         0.1 + 100. * abs( min( 0., &
                         (State_VGB(iVar,i,j,k,iBlock)-&
                         StateOld_VCB(iVar,i,j,k,iBlock)) &
                         /StateOld_VCB(iVar,i,j,k,iBlock) ) ) )

                    percent_chg_rho(2) = max(percent_chg_rho(2), &
                         0.1 + 100. * abs( max( 0., &
                         (State_VGB(iVar,i,j,k,iBlock)-&
                         StateOld_VCB(iVar,i,j,k,iBlock)) &
                         /StateOld_VCB(iVar,i,j,k,iBlock) ) ) )
                 end do
              end if
              percent_chg_p(1) = 0.1 + 100. * abs( min( 0., &
                   (State_VGB(P_,i,j,k,iBlock)-&
                   StateOld_VCB(P_,i,j,k,iBlock)) &
                   /StateOld_VCB(P_,i,j,k,iBlock) ) )
              percent_chg_p(2) = 0.1 + 100. * abs( max( 0., &
                   (State_VGB(P_,i,j,k,iBlock)-&
                   StateOld_VCB(P_,i,j,k,iBlock)) &
                   /StateOld_VCB(P_,i,j,k,iBlock) ) )
              time_fraction_rho = 1.0 / maxval(percent_chg_rho/percent_max_rho)
              time_fraction_p   = 1.0 / maxval(percent_chg_p  /percent_max_p  )
              if (time_fraction_rho < 1. .or. time_fraction_p < 1.) then
                 if (num_checks == 1) then
                    time_fraction = 1.
                    if (time_fraction_rho < 1.) &
                         time_fraction = 0.9*time_fraction_rho
                    if (time_fraction_p   < 1.) &
                         time_fraction = min(time_fraction, 0.75)
                 else
                    time_fraction = 0.5
                 end if
                 update_check_done = .false.
                 cell_time_fraction = cell_time_fraction * time_fraction
                 if(oktest2) then
                    write(*,*) &
                         'update_check LT: changing cell value, PE=',iProc, &
                         ' BLK=',iBlock,' i,j,k=',i,' ',j,' ',k, &
                         '  time_fraction=',time_fraction
                    write(*,*) &
                         iProc,' ',iBlock,' ',i,' ',j,' ',k,' OLD:  ', &
                         NameVar_V(1),'=',StateOld_VCB(1,i,j,k,iBlock),'    ',&
                         NameVar_V(nVar),' ', StateOld_VCB(nVar,i,j,k,iBlock)
                    write(*,*) &
                         iProc,' ',iBlock,' ',i,' ',j,' ',k,' BAD:', &
                         ' ', NameVar_V(1),'=',State_VGB(1,i,j,k,iBlock), &
                         '   ', NameVar_V(nVar),State_VGB(nVar,i,j,k,iBlock)
                 end if
                 call fix_update
                 if(oktest2) then
                    write(*,*) &
                         iProc,' ',iBlock,' ',i,' ',j,' ',k,' NEW:', &
                         ' ',NameVar_V(Rho_),'=',State_VGB(Rho_,i,j,k,iBlock),&
                         '   ',NameVar_V(p_),'=', State_VGB(p_,i,j,k,iBlock)
                 end if
              end if
              if (update_check_done) EXIT
           end do
           PercentChangePE(1:2) = &
                max(percent_chg_rho(1:2)-0.1, PercentChangePE(1:2))
           PercentChangePE(3:4) = &
                max(percent_chg_p(1:2)-0.1, PercentChangePE(3:4))
           report_tf = min(report_tf, cell_time_fraction)
        end do; end do; end do
     end do
     call MPI_allreduce(report_tf, report_tf_all, 1, &
          MPI_REAL, MPI_MIN, iComm, iError)
     report_tf = report_tf_all
     if(oktest) then
        if (iProc == 0 .and. report_tf < 1.) &
             write(*,'(a,a,i6,a,f12.8)') 'update_check LT:', &
             ' nStep=',n_step,' max dt reduction=',report_tf
     end if
  end if

  if(oktest1.and.iStage==nStage) then
     call MPI_allreduce(PercentChangePE,  PercentChangeMax, 4, &
          MPI_REAL, MPI_MAX, iComm, iError)
     if(iProc==0) then
        write(*,*)'Maximum change in pressure on proc 0:',&
             - PercentChangeMax(3),' %,   ',  PercentChangeMax(4),' %'
        write(*,*) 'Maximum change in other positive variables:', &
             - PercentChangeMax(1),' %,   ',  PercentChangeMax(2),' %'
     end if
     if(oktest3) then
        do iBlock = 1,nBlockMax
           if(unusedBLK(iBlock))cycle
           do k=1,nK;do j=1,nJ; do i=1,nI
              if(.not.true_cell(i,j,k,iBlock))cycle
              if (abs(100. * abs( min( 0., &
                   (State_VGB(Rho_,i,j,k,iBlock)-&
                   StateOld_VCB(Rho_,i,j,k,iBlock)) &
                   /StateOld_VCB(Rho_,i,j,k,iBlock) ) )-&
                   PercentChangeMax(1))<cTiny*PercentChangeMax(1))&
                   write(*,*)'Maximum decrease in density at X Y Z=',&
                   x_BLK(i,j,k,iBlock), &
                   y_BLK(i,j,k,iBlock), &
                   z_BLK(i,j,k,iBlock),&
                   ': rho_old = ',StateOld_VCB(Rho_,i,j,k,iBlock),&
                   ' rho_new = ',State_VGB(Rho_,i,j,k,iBlock)
                 
              if (abs(100. * abs( max( 0., &
                   (State_VGB(Rho_,i,j,k,iBlock)-&
                   StateOld_VCB(Rho_,i,j,k,iBlock)) &
                   /StateOld_VCB(Rho_,i,j,k,iBlock) ) )-&
                   PercentChangeMax(2))<cTiny*PercentChangeMax(2))&
                   write(*,*)'Maximum increase in density at the point',&
                   x_BLK(i,j,k,iBlock), &
                   y_BLK(i,j,k,iBlock), &
                   z_BLK(i,j,k,iBlock),&
                   'is: rho_old = ',&
                   StateOld_VCB(Rho_,i,j,k,iBlock),&
                   'rho_new=',State_VGB(Rho_,i,j,k,iBlock)

              if (abs(100. * abs( min( 0., &
                   (State_VGB(p_,i,j,k,iBlock)-&
                   StateOld_VCB(p_,i,j,k,iBlock)) &
                   /StateOld_VCB(p_,i,j,k,iBlock) ) )-&
                   PercentChangeMax(3))<cTiny*PercentChangeMax(3))&
                   write(*,*)'Maximum decrease in',NameVar_V(p_), &
                   'at the point',&
                   x_BLK(i,j,k,iBlock), &
                   y_BLK(i,j,k,iBlock), &
                   z_BLK(i,j,k,iBlock),&
                   'is: valeu_old = ',StateOld_VCB(P_,i,j,k,iBlock),&
                   'value_new=',State_VGB(p_,i,j,k,iBlock)
              if (abs(100. * abs( max( 0., &
                   (State_VGB(p_,i,j,k,iBlock)-&
                   StateOld_VCB(p_,i,j,k,iBlock)) &
                   /StateOld_VCB(p_,i,j,k,iBlock) ) )-&
                   PercentChangeMax(4))<cTiny*PercentChangeMax(4))&
                   write(*,*)'Maximum increase in',NameVar_V(p_), &
                   'at the point',&
                   x_BLK(i,j,k,iBlock), &
                   y_BLK(i,j,k,iBlock), &
                   z_BLK(i,j,k,iBlock),&
                   'is: value_old = ',StateOld_VCB(p_,i,j,k,iBlock),&
                   'value_new=',State_VGB(p_,i,j,k,iBlock)
           end do; end do; end do
        end do
     end if
  end if

  if(iProc == 0 .and. report_tf<1.)&
       call error_report('Time step reduction, min(factor)',&
       report_tf,iError1,.true.)

  ! Check for positivity of variables
  IsNegative = .false.
  do iBlock = 1, nBlockMax
     if (unusedBLK(iBlock)) CYCLE
     do iVar = 1, nVar
        ! Ignore variables that do not have to be positive
        if(DefaultState_V(iVar) < cTiny) CYCLE

        ! Do not check species densities if check threshold is positive
        ! (i.e. minor species densities are allowed to go negative,
        !       and they are fixed to be zero in update_states)
        if(UseMultiSpecies .and. SpeciesPercentCheck > 0.0 &
           .and. iVar >= SpeciesFirst_ .and. iVar <= SpeciesLast_) CYCLE

        Value = minval(State_VGB(iVar,1:nI,1:nJ,1:nK,iBlock))

        if(Value < 0.0)then
           i_D = minloc(State_VGB(iVar,1:nI,1:nJ,1:nK,iBlock))
           i = i_D(1); j = i_D(2); k = i_D(3)
           write (*,'(a,3i3,2i5,i3,a,3f12.4,/,5x,a,a,es12.4)') &
                ' I J K iBlock iProc iVar=',i_D,iBlock,iProc,iVar, &
                ' X Y Z=', &
                x_BLK(i,j,k,iBlock), &
                y_BLK(i,j,k,iBlock), &
                z_BLK(i,j,k,iBlock), &
                ' Var='//trim(NameVar_V(iVar)), &
                ' Value=', State_VGB(iVar,i,j,k,iBlock)
           IsNegative = .true.
        end if
     end do
  end do
  if (IsNegative) then
     if(time_accurate) then
        write(*,'(a,i4,a,a,i6,a,f12.8,a,f12.8)') &
             'Negative updated value: PE=',iProc, &
             'update_check TA:',' nStep=',n_step, &
             '     dt reduction=',report_tf,' dt=',dt
     else
        write(*,'(a,i4,a,a,i6,a,f12.8)') &
             'Negative updated value: PE=',iProc, &
             'update_check LT:',' nStep=',n_step, &
             ' max dt reduction=',report_tf
     end if
  end if

  call MPI_allreduce(IsNegative, DoStop,1,MPI_LOGICAL,MPI_LOR,iComm,iError)
  if (DoStop) call stop_mpi('Stopping, negative density or pressure')

contains

  subroutine fix_update

    logical :: IsConserv

    real :: fullBx, fullBy, fullBz, fullBB, UdotBc2, rhoc2, gA2_Boris
    real :: rhoUx_Boris, rhoUy_Boris, rhoUz_Boris, E_Boris, &
         rhoUx_o_Boris, rhoUy_o_Boris, rhoUz_o_Boris, E_o_Boris

    character(len=*), parameter:: NameSub='fix_update'
    logical, parameter :: DoTestCell = .false.
    !-------------------------------------------------------------------------
    !DoTestCell = DoTestMe .and. iBlock==BlkTest .and. &
    !     i==iTest .and. j==jTest .and. k==kTest 

    if(allocated(IsConserv_CB))then
       IsConserv = IsConserv_CB(i,j,k,iBlock)
    else
       IsConserv = .not. UseNonConservative
    end if

    if(DoTestCell)then
       write(*,*)NameSub,' IsConserv    =',IsConserv
       write(*,*)NameSub,' initial state=',State_VGB(:,i,j,k,iBlock)
       write(*,*)NameSub,' old     state=',StateOld_VCB(:,i,j,k,iBlock)
    end if

    if(boris_correction .and. nFluid==1) then                  !^CFG IF BORISCORR BEGIN

       ! Convert old state
       fullBx = B0_DC(x_,i,j,k) + &
            StateOld_VCB(Bx_,i,j,k,iBlock)
       fullBy = B0_DC(y_,i,j,k) + &
            StateOld_VCB(By_,i,j,k,iBlock)
       fullBz = B0_DC(z_,i,j,k) + &
            StateOld_VCB(Bz_,i,j,k,iBlock)
       fullBB = fullBx**2 + fullBy**2 + fullBz**2
       rhoc2  = &
            StateOld_VCB(rho_,i,j,k,iBlock)*c2LIGHT
       UdotBc2= (StateOld_VCB(rhoUx_,i,j,k,iBlock)*fullBx + &
            StateOld_VCB(rhoUy_,i,j,k,iBlock)*fullBy + &
            StateOld_VCB(rhoUz_,i,j,k,iBlock)*fullBz)/rhoc2
       gA2_Boris=1.+fullBB/rhoc2


       ! rhoU_Boris = rhoU - ((U x B) x B)/c^2 
       !            = rhoU + (U B^2 - B U.B)/c^2
       !            = rhoU*(1+BB/(rho*c2)) - B UdotB/c^2
       rhoUx_o_Boris = StateOld_VCB(rhoUx_,i,j,k,iBlock)*ga2_Boris - &
            fullBx*UdotBc2
       rhoUy_o_Boris = StateOld_VCB(rhoUy_,i,j,k,iBlock)*ga2_Boris - &
            fullBy*UdotBc2
       rhoUz_o_Boris = StateOld_VCB(rhoUz_,i,j,k,iBlock)*ga2_Boris - &
            fullBz*UdotBc2

       if(IsConserv)then
          ! e_boris = e + 0.5/c^2 * (V x B)^2
          E_o_Boris = EnergyOld_CBI(i,j,k,iBlock,1) + (cHalf/c2LIGHT)*( &
               ((StateOld_VCB(rhoUy_,i,j,k,iBlock)*fullBz     &
               -StateOld_VCB(rhoUz_,i,j,k,iBlock)*fullBy)**2 &
               +(StateOld_VCB(rhoUx_,i,j,k,iBlock)*fullBz     &
               -StateOld_VCB(rhoUz_,i,j,k,iBlock)*fullBx)**2 &
               +(StateOld_VCB(rhoUx_,i,j,k,iBlock)*fullBy     &
               -StateOld_VCB(rhoUy_,i,j,k,iBlock)*fullBx)**2 &
               )/StateOld_VCB(rho_,i,j,k,iBlock)**2                 )
       end if

       ! Convert current state
       fullBx = B0_DC(x_,i,j,k) + State_VGB(Bx_,i,j,k,iBlock)
       fullBy = B0_DC(y_,i,j,k) + State_VGB(By_,i,j,k,iBlock)
       fullBz = B0_DC(z_,i,j,k) + State_VGB(Bz_,i,j,k,iBlock)
       fullBB = fullBx**2 + fullBy**2 + fullBz**2
       rhoc2  = State_VGB(rho_,i,j,k,iBlock)*c2LIGHT
       UdotBc2= (State_VGB(rhoUx_,i,j,k,iBlock)*fullBx + &
            State_VGB(rhoUy_,i,j,k,iBlock)*fullBy + &
            State_VGB(rhoUz_,i,j,k,iBlock)*fullBz)/rhoc2
       gA2_Boris=cOne+fullBB/rhoc2


       ! rhoU_Boris = rhoU - ((U x B) x B)/c^2 
       !            = rhoU + (U B^2 - B U.B)/c^2
       !            = rhoU*(1+BB/(rho*c2)) - B UdotB/c^2
       rhoUx_Boris = State_VGB(rhoUx_,i,j,k,iBlock)*ga2_Boris - fullBx*UdotBc2
       rhoUy_Boris = State_VGB(rhoUy_,i,j,k,iBlock)*ga2_Boris - fullBy*UdotBc2
       rhoUz_Boris = State_VGB(rhoUz_,i,j,k,iBlock)*ga2_Boris - fullBz*UdotBc2

       if(IsConserv)then
          ! e_boris = e + 0.5/c^2 * (V x B)^2
          E_Boris = Energy_GBI(i,j,k,iBlock,1) + (cHalf/c2LIGHT)*( &
               ((State_VGB(rhoUy_,i,j,k,iBlock)*fullBz     &
               -State_VGB(rhoUz_,i,j,k,iBlock)*fullBy)**2 &
               +(State_VGB(rhoUx_,i,j,k,iBlock)*fullBz     &
               -State_VGB(rhoUz_,i,j,k,iBlock)*fullBx)**2 &
               +(State_VGB(rhoUx_,i,j,k,iBlock)*fullBy     &
               -State_VGB(rhoUy_,i,j,k,iBlock)*fullBx)**2 &
               )/State_VGB(rho_,i,j,k,iBlock)**2                 )
       end if

       ! Interpolate
       !For possible extension to multifluid
       !State_VGB(iRho_I,i,j,k,iBlock) = &
       !     (   time_fraction) *   State_VGB(iRho_I,i,j,k,iBlock) + &
       !     (cOne-time_fraction) * StateOld_VCB(iRho_I,i,j,k,iBlock)
       State_VGB(rho_,i,j,k,iBlock) = &
            (   time_fraction) *   State_VGB(rho_,i,j,k,iBlock) + &
            (cOne-time_fraction) * StateOld_VCB(rho_,i,j,k,iBlock)
       rhoUx_Boris = &
            (   time_fraction) * rhoUx_Boris + &
            (cOne-time_fraction) * rhoUx_o_Boris
       rhoUy_Boris = &
            (   time_fraction) * rhoUy_Boris + &
            (cOne-time_fraction) * rhoUy_o_Boris
       rhoUz_Boris = &
            (   time_fraction) * rhoUz_Boris + &
            (cOne-time_fraction) * rhoUz_o_Boris
       State_VGB(Bx_,i,j,k,iBlock) = &
            (   time_fraction) *   State_VGB(Bx_,i,j,k,iBlock) + &
            (cOne-time_fraction) * StateOld_VCB(Bx_,i,j,k,iBlock)
       State_VGB(By_,i,j,k,iBlock) = &
            (   time_fraction) *   State_VGB(By_,i,j,k,iBlock) + &
            (cOne-time_fraction) * StateOld_VCB(By_,i,j,k,iBlock)
       State_VGB(Bz_,i,j,k,iBlock) = &
            (   time_fraction) *   State_VGB(Bz_,i,j,k,iBlock) + &
            (cOne-time_fraction) * StateOld_VCB(Bz_,i,j,k,iBlock)

       ! Convert Back
       fullBx = B0_DC(x_,i,j,k) + State_VGB(Bx_,i,j,k,iBlock)
       fullBy = B0_DC(y_,i,j,k) + State_VGB(By_,i,j,k,iBlock)
       fullBz = B0_DC(z_,i,j,k) + State_VGB(Bz_,i,j,k,iBlock)
       fullBB = fullBx**2 + fullBy**2 + fullBz**2
       rhoc2  = State_VGB(rho_,i,j,k,iBlock)*c2LIGHT
       UdotBc2= (rhoUx_Boris*fullBx + rhoUy_Boris*fullBy + rhoUz_Boris*fullBz)&
            /rhoc2
       gA2_Boris= cOne/(cOne+fullBB/rhoc2)

       ! rhoU = 1/(rho c^2 + B^2) * (I rho c^2 + B B) * rhoU_Boris
       !      = 1/[1+BB/(rho c^2)]* (rhoU_Boris + (rhoUBorisdotB/(rho c^2) * B)

       State_VGB(rhoUx_,i,j,k,iBlock) =  gA2_Boris * (rhoUx_Boris + UdotBc2*fullBx)
       State_VGB(rhoUy_,i,j,k,iBlock) =  gA2_Boris * (rhoUy_Boris + UdotBc2*fullBy)
       State_VGB(rhoUz_,i,j,k,iBlock) =  gA2_Boris * (rhoUz_Boris + UdotBc2*fullBz)

       if(IsConserv)then
          E_boris= &
               (   time_fraction) * E_Boris + &
               (cOne-time_fraction) * E_o_Boris

          ! E = E_boris - 0.5/c^2 * (V x B)^2
          Energy_GBI(i,j,k,iBlock,1) = E_Boris - (cHalf/c2LIGHT)*( &
               ((State_VGB(rhoUy_,i,j,k,iBlock)*fullBz     &
               -State_VGB(rhoUz_,i,j,k,iBlock)*fullBy)**2 &
               +(State_VGB(rhoUx_,i,j,k,iBlock)*fullBz     &
               -State_VGB(rhoUz_,i,j,k,iBlock)*fullBx)**2 &
               +(State_VGB(rhoUx_,i,j,k,iBlock)*fullBy     &
               -State_VGB(rhoUy_,i,j,k,iBlock)*fullBx)**2 &
               )/State_VGB(rho_,i,j,k,iBlock)**2               )


          if((nStage==1.and..not.time_accurate).or.&
               (nStage>1.and.iStage==1)) &
               Energy_GBI(i,j,k,iBlock,1) =  Energy_GBI(i,j,k,iBlock,1) - &
               (cHalf/time_fraction - cHalf)*&
               ((State_VGB(Bx_,i,j,k,iBlock) - &
               StateOld_VCB(Bx_,i,j,k,iBlock))**2 +&
               (State_VGB(By_,i,j,k,iBlock) - &
               StateOld_VCB(By_,i,j,k,iBlock))**2 +&
               (State_VGB(Bz_,i,j,k,iBlock) - &
               StateOld_VCB(Bz_,i,j,k,iBlock))**2 )

          call calc_pressure1_point(i,j,k,iBlock)

          ! For multifluid update all other energies and call calc_pressure_point
       else
          ! For possible extension to multifluid
          !State_VGB(iP_I,i,j,k,iBlock) = &
          !     (   time_fraction) *   State_VGB(iP_I,i,j,k,iBlock) + &
          !     (cOne-time_fraction) * StateOld_VCB(iP_I,i,j,k,iBlock)
          !call calc_energy_point(i,j,k,iBlock)

          State_VGB(p_,i,j,k,iBlock) = &
               (   time_fraction) *   State_VGB(p_,i,j,k,iBlock) + &
               (cOne-time_fraction) * StateOld_VCB(p_,i,j,k,iBlock)

          call calc_energy1_point(i,j,k,iBlock)
       end if
    else                      ! Non-Boris interpolation                
       !^CFG END BORISCORR
       State_VGB(1:nVar,i,j,k,iBlock) = &
            (   time_fraction) *   State_VGB(1:nVar,i,j,k,iBlock) + &
            (cOne-time_fraction) * StateOld_VCB(1:nVar,i,j,k,iBlock)
       if(IsConserv)then
          Energy_GBI(i,j,k,iBlock,:) = &
               (   time_fraction) *   Energy_GBI(i,j,k,iBlock,:) + &
               (cOne-time_fraction) * EnergyOld_CBI(i,j,k,iBlock,:)

          if(IsMhd .and. (nStage==1.and..not.time_accurate).or.&
               (nStage>1.and.iStage==1)) then
             Energy_GBI(i,j,k,iBlock,1) = &
                  Energy_GBI(i,j,k,iBlock,1) - &
                  (cHalf/time_fraction - cHalf)*&
                  sum((State_VGB(Bx_:Bz_,i,j,k,iBlock) - &
                  StateOld_VCB(Bx_:Bz_,i,j,k,iBlock))**2)
          end if

          call calc_pressure_point(i,j,k,iBlock)
       else
          call calc_energy_point(i,j,k,iBlock)
       end if
    end if                       !^CFG IF BORISCORR
    time_BLK(i,j,k,iBlock) = time_BLK(i,j,k,iBlock)*time_fraction

    if(DoTestCell)write(*,*)NameSub,' final state=',State_VGB(:,i,j,k,iBlock)

  end subroutine fix_update

end subroutine update_check

!=============================================================================
subroutine select_conservative

  ! Set the global variable IsConserv_CB

  use ModNumConst
  use ModMain
  use ModAdvance
  use ModGeometry
  implicit none

  integer :: iBlock, iCrit, i, j, k

  real :: UxDx_G( 0:nI+1, 1:nJ,   1:nK  )
  real :: UyDy_G( 1:nI,   0:nJ+1, 1:nK  )
  real :: UzDz_G( 1:nI,   1:nJ,   0:nK+1)
  real :: DivU_C( 1:nI,   1:nJ, 1:nK)
  real :: aFast_C(1:nI,   1:nJ, 1:nK)

  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------
  call set_oktest('select_conservative',DoTest,DoTestMe)

  call timing_start('nonconservative')

  if(DoTestMe)write(*,*)'select_conservative: starting with ',&
       'UseNonConservative, nConservCrit=',UseNonConservative, nConservCrit

  if(.not.allocated(IsConserv_CB))then
     allocate(IsConserv_CB(nI,nJ,nK,MaxBlock))
     if(DoTestMe)write(*,*)'select_conservative: allocated IsConserv_CB'
  end if

  if(nConservCrit < 1)then
     ! There are no criteria so use fully non-conservative
     IsConserv_CB = .false.
     if(DoTestMe)write(*,*)'select_conservative: set IsConserv_CB = F'
     RETURN
  endif

  if(any(TypeConservCrit_I == 'p' .or. TypeConservCrit_I == 'gradp' &
       .or. TypeConservCrit_I == 'jumpp') )then

     if(DoTestMe)write(*,*)'select_conservative: Apply physics based criteria'

     ! These all have to be true to use non-conservative, 
     ! so any of them can switch the original non-conservative to conservative
     IsConserv_CB = .false.

     do iBlock = 1, nBlock
        if( UnusedBlk(iBlock) ) CYCLE

        do iCrit = 1, nConservCrit
           select case(TypeConservCrit_I(iCrit))
           case('p')
              if(UseB0)then
                 do k=1,nK; do j=1,nJ; do i=1,nI
                    IsConserv_CB(i,j,k,iBlock) = &
                         IsConserv_CB(i,j,k,iBlock) .or. &
                         State_VGB(P_,i,j,k,iBlock) > pCoeffConserv * &
                         (Energy_GBI(i,j,k,iBlock,1) + 0.5 * &
                         ((State_VGB(Bx_,i,j,k,iBlock) &
                         + B0_DGB(x_,i,j,k,iBlock))**2 &
                         +(State_VGB(By_,i,j,k,iBlock) &
                         + B0_DGB(y_,i,j,k,iBlock))**2 &
                         +(State_VGB(Bz_,i,j,k,iBlock) &
                         + B0_DGB(z_,i,j,k,iBlock))**2 &
                         -State_VGB(Bx_,i,j,k,iBlock)**2 &
                         -State_VGB(By_,i,j,k,iBlock)**2 &
                         -State_VGB(Bz_,i,j,k,iBlock)**2 &
                         ))
                 end do;end do;end do
              else
                 do k=1,nK; do j=1,nJ; do i=1,nI
                    IsConserv_CB(i,j,k,iBlock) = &
                         IsConserv_CB(i,j,k,iBlock) .or. &
                         State_VGB(P_,i,j,k,iBlock) > pCoeffConserv * &
                         Energy_GBI(i,j,k,iBlock,1) 
                 end do;end do;end do
              end if
           case('gradp')
              ! Switch to conservative if gradient of pressure is large
              do k=1,nK; do j=1,nJ; do i=1,nI
                 IsConserv_CB(i,j,k,iBlock) = IsConserv_CB(i,j,k,iBlock) .or. &
                      (abs(State_VGB(P_,i+1,j,k,iBlock) &
                      -    State_VGB(P_,i-1,j,k,iBlock))  &
                      +abs(State_VGB(P_,i,j+1,k,iBlock) &
                      -    State_VGB(P_,i,j-1,k,iBlock))  &
                      +abs(State_VGB(P_,i,j,k+1,iBlock) &
                      -    State_VGB(P_,i,j,k-1,iBlock))) &
                      > GradPCoeffConserv * min(    &
                      State_VGB(P_,i,j,k,iBlock),   &
                      State_VGB(P_,i+1,j,k,iBlock), &
                      State_VGB(P_,i-1,j,k,iBlock), &
                      State_VGB(P_,i,j+1,k,iBlock), &
                      State_VGB(P_,i,j-1,k,iBlock), &
                      State_VGB(P_,i,j,k+1,iBlock), &
                      State_VGB(P_,i,j,k-1,iBlock))
              end do; end do; end do
           case('jumpp')
              ! Switch to conservative if pressure jump is large
              do k=1,nK; do j=1,nJ; do i=1,nI
                 IsConserv_CB(i,j,k,iBlock) = IsConserv_CB(i,j,k,iBlock) .or. &
                      maxval(State_VGB(P_,i-2:i+2,j,k,iBlock)) &
                      > GradPCoeffConserv* &
                      minval(State_VGB(P_,i-2:i+2,j,k,iBlock)) .or. &
                      nJ > 1 .and. &
                      maxval(State_VGB(P_,i,j-2:j+2,k,iBlock)) &
                      > GradPCoeffConserv* &
                      minval(State_VGB(P_,i:i,j-2:j+2,k,iBlock)) .or. &
                      nK > 1 .and. &
                      maxval(State_VGB(P_,i,j,k-2:k+2,iBlock)) &
                      > GradPCoeffConserv* &
                      minval(State_VGB(P_,i:i,j,k-2:k+2,iBlock))
              end do; end do; end do
           case default
              CYCLE
           end select

           if(DoTestMe.and.iBlock==BlkTest)&
                write(*,*)'select_conservative: TypeCrit, IsConserv=',&
                TypeConservCrit_I(iCrit), &
                IsConserv_CB(iTest,jTest,kTest,iBlock)
        end do
     end do
  else
     ! If there are no physics based criteria we start from 
     ! the assumption of conservative everywhere
     IsConserv_CB = .true.

     if(DoTestMe.and.iBlock==BlkTest)&
          write(*,*)'select_conservative: default IsConserv is true'
  endif

  do iBlock = 1, nBlock
     if( UnusedBlk(iBlock) ) CYCLE

     ! Apply geometry based criteria
     ! Any of these can switch from conservative to non-conservative
     do iCrit = 1, nConservCrit
        select case(TypeConservCrit_I(iCrit))
        case('r')
           ! Switch to non-conservative inside radius rConserv
           IsConserv_CB(:,:,:,iBlock) = IsConserv_CB(:,:,:,iBlock) .and. &
                R_BLK(1:nI,1:nJ,1:nK,iBlock) > rConserv
        case('parabola')
           ! Switch to non-conservative behind a parabola inside the bow shock
           IsConserv_CB(:,:,:,iBlock) = IsConserv_CB(:,:,:,iBlock) .and. &
                x_BLK(1:nI,1:nJ,1:nK,iBlock) > xParabolaConserv - &
                ( y_BLK(1:nI,1:nJ,1:nK,iBlock)**2 &
                + z_BLK(1:nI,1:nJ,1:nK,iBlock)**2 ) / yParabolaConserv
        case default
           CYCLE
        end select
        if(DoTestMe.and.iBlock==BlkTest)&
             write(*,*)'select_conservative: TypeCrit, IsConserv=',&
             TypeConservCrit_I(iCrit), IsConserv_CB(iTest,jTest,kTest,iBlock)
     end do
  end do

  call timing_stop('nonconservative')

end subroutine select_conservative
