!^CFG COPYRIGHT UM
subroutine calc_timestep
  use ModProcMH
  use ModMain
  use ModAdvance, ONLY : VdtFace_x,VdtFace_y,VdtFace_z,time_BLK
  use ModNumConst
  use ModGeometry, ONLY : true_cell,true_BLK,&
       iVolumeCounterBLK,iVolumeCounterI,&               !^CFG IF NOT CARTESIAN
       VolumeInverse_I
  implicit none

  logical :: oktest, oktest_me
  integer :: iVolumeCounter,i,j,k
  !--------------------------------------------------------------------------

  if(globalBLK==BLKtest)then
     call set_oktest('calc_timestep',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  endif

  !\
  ! Compute the allowable local time step 
  ! for each cell based on local face values
  !/
  iVolumeCounter=globalBLK
  iVolumeCounter=iVolumeCounterBLK*(iVolumeCounter-iVolumeCounterI) !^CFG IF NOT CARTESIAN
  do k=1,nK
     do j=1,nJ
        do i=1,nI
           iVolumeCounter=iVolumeCounter+iVolumeCounterI !^CFG IF NOT CARTESIAN
           time_BLK(i,j,k,globalBLK) = cOne /(VolumeInverse_I(iVolumeCounter)&
                *(max(VdtFace_x(i,j,k),VdtFace_x(i+1,j,k))+ &
                max(VdtFace_y(i,j,k),VdtFace_y(i,j+1,k))+ &
                max(VdtFace_z(i,j,k),VdtFace_z(i,j,k+1))))
        end do
     end do
  end do

  if(oktest_me)then
     write(*,*)'left  VdtFace_x,y,z=',&
          VdtFace_x(Itest,Jtest,Ktest),&
          VdtFace_y(Itest,Jtest,Ktest),&
          VdtFace_z(Itest,Jtest,Ktest)
     write(*,*)'right VdtFace_x,y,z=',&
          VdtFace_x(Itest+1,Jtest,Ktest),&
          VdtFace_y(Itest,Jtest+1,Ktest),&
          VdtFace_z(Itest,Jtest,Ktest+1)
     write(*,*)'time_BLK=',time_BLK(Itest,Jtest,Ktest,globalBLK)
  end if

  !\
  ! Compute maximum stable time step for this solution block 
  !/
  if(true_BLK(globalBLK)) then
     dt_BLK(globalBLK) = minval(time_BLK(:,:,:,globalBLK))
  else
     ! If the block has no true cells, set dt_BLK=1.0E20
     dt_BLK(globalBLK) = min(cHuge,&
          minval(time_BLK(:,:,:,globalBLK), &
          MASK=true_cell(1:nI,1:nJ,1:nK,globalBLK)))
  end if

  if(oktest_me)write(*,*)'dt_BLK, loc=',dt_BLK(globalBLK),&
       minloc(time_BLK(:,:,:,globalBLK),&
       MASK=true_cell(1:nI,1:nJ,1:nK,globalBLK))

  ! Reset time_BLK for fixed time step (but dt_BLK is kept!)
  if(UseDtFixed &
       .and..not.UsePartLocal &               !^CFG IF IMPLICIT
       ) time_BLK(:,:,:,globalBLK) = dt

  ! Set time step to zero inside body.
  if(.not.true_BLK(globalBLK)) then
     where (.not.true_cell(1:nI,1:nJ,1:nK,globalBLK))&
          time_BLK(:,:,:,globalBLK) = cZero
  end if

end subroutine calc_timestep

subroutine set_global_timestep
  use ModProcMH
  use ModMain
  use ModAdvance, ONLY : time_BLK,State_VGB,rho_,Bx_,By_,Bz_,P_
  use ModAdvance,ONLY:   B0xCell_BLK,B0yCell_BLK,B0zCell_BLK
  use ModGeometry, ONLY : true_cell,true_BLK,dx_BLK,XyzStart_BLK
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK
  use ModImplicit, ONLY: UsePartImplicit, implicitBLK     !^CFG IF IMPLICIT
  use ModPhysics,ONLY:UnitSI_x,UnitSI_U,UnitSI_t,UnitSI_B,UnitSI_rho,g
  use ModMpi
  implicit none

  integer :: iBlock
  integer :: iError,lPosition_D(3)
  real    :: dtMinPE,UAlfvMax

  logical :: oktest, oktest_me
  !--------------------------------------------------------------------------

  call set_oktest('calc_timestep',oktest,oktest_me)
  if(oktest_me)write(*,*)'Starting set_global_timestep'

  if(UseDtFixed)then
     dt = DtFixed
  else
     !\
     ! Impose global time step for time-accurate 
     ! calculations as required
     !/
     if(UsePartImplicit .or. UsePartLocal)then        !^CFG IF IMPLICIT BEGIN
        ! Implicit blocks are not taken into account for partially implicit run
        dtMinPE = minval(dt_BLK(1:nBlock),&
             MASK=.not.(unusedBLK(1:nBlock) .or. implicitBLK(1:nBlock)))
     else                                             !^CFG END IMPLICIT
        dtMinPE = minval(dt_BLK(1:nBlock), MASK=.not.unusedBLK(1:nBlock))
     end if                                           !^CFG IF IMPLICIT

     ! Set dt to minimum time step over all the PE-s
     call MPI_allreduce(dtMinPE, dt, 1, MPI_REAL, MPI_MIN, iComm, iError)

     if(oktest .and. dtMinPE==dt)then
        do iBlock=1,nBlock
           if(unusedBLK(iBlock)) CYCLE
           if(dt_BLK(iBlock)==dt)then
              write(*,*)'Time step dt=',dt,'=', dt*UnitSI_t,&
                   ' s  is comtrolled by block with PE, iBlock=', iProc, iBlock
              write(*,*)'The coordinates of (1,1,1) cell center are ',&
                   XyzStart_BLK(:,iBlock)
              write(*,*)'Block size, dimensionless and physical, are',&
                   dx_BLK(iBlock), ' , ', dx_BLK(iBlock)*UnitSI_x,' m'
              UAlfvMax=sqrt(&
                   maxval(&
                   ((State_VGB(Bx_,1:nI,1:nJ,1:nK,iBlock)+&
                     B0xCell_BLK(1:nI,1:nJ,1:nK,iBlock))**2+&
                    (State_VGB(By_,1:nI,1:nJ,1:nK,iBlock)+&
                     B0yCell_BLK(1:nI,1:nJ,1:nK,iBlock))**2+&
                    (State_VGB(Bz_,1:nI,1:nJ,1:nK,iBlock)+&
                     B0zCell_BLK(1:nI,1:nJ,1:nK,iBlock))**2+&
                     g*State_VGB(P_,1:nI,1:nJ,1:nK,iBlock))/&
                     State_VGB(rho_,1:nI,1:nJ,1:nK,iBlock),&
                     MASK=true_cell(1:nI,1:nJ,1:nK,iBlock)))
              lPosition_D=maxloc(((State_VGB(Bx_,1:nI,1:nJ,1:nK,iBlock)+&
                     B0xCell_BLK(1:nI,1:nJ,1:nK,iBlock))**2+&
                    (State_VGB(By_,1:nI,1:nJ,1:nK,iBlock)+&
                     B0yCell_BLK(1:nI,1:nJ,1:nK,iBlock))**2+&
                    (State_VGB(Bz_,1:nI,1:nJ,1:nK,iBlock)+&
                     B0zCell_BLK(1:nI,1:nJ,1:nK,iBlock))**2+&
                     g*State_VGB(P_,1:nI,1:nJ,1:nK,iBlock))/&
                     State_VGB(rho_,1:nI,1:nJ,1:nK,iBlock),&
                     MASK=true_cell(1:nI,1:nJ,1:nK,iBlock))
              write(*,*)'Maximal magnetosonic speed =',UAlfvMax*UnitSI_U,&
                     ' m/s is reached in the point with coords=',&
                     x_BLK(lPosition_D(1),lPosition_D(2),lPosition_D(3),iBlock),&
                     y_BLK(lPosition_D(1),lPosition_D(2),lPosition_D(3),iBlock),&
                     z_BLK(lPosition_D(1),lPosition_D(2),lPosition_D(3),iBlock)
              write(*,*)'Magnetic field in this point: B0:',&
                     B0xCell_BLK(lPosition_D(1),lPosition_D(2),lPosition_D(3),iBlock)*UnitSI_B,&
                     B0yCell_BLK(lPosition_D(1),lPosition_D(2),lPosition_D(3),iBlock)*UnitSI_B,&
                     B0zCell_BLK(lPosition_D(1),lPosition_D(2),lPosition_D(3),iBlock)*UnitSI_B,&
                     ' T,  B1:',&
                     State_VGB(Bx_,lPosition_D(1),lPosition_D(2),lPosition_D(3),iBlock)*UnitSI_B,&
                     State_VGB(By_,lPosition_D(1),lPosition_D(2),lPosition_D(3),iBlock)*UnitSI_B,&
                     State_VGB(Bz_,lPosition_D(1),lPosition_D(2),lPosition_D(3),iBlock)*UnitSI_B,&
                     ' T,  Density=',&
                     State_VGB(rho_,lPosition_D(1),lPosition_D(2),lPosition_D(3),iBlock)*UnitSI_rho
              EXIT
           end if
        end do
     end if

  end if

  do iBlock = 1,nBlock
     if (unusedBLK(iBlock)) CYCLE

     if(UsePartLocal)then                             !^CFG IF IMPLICIT BEGIN
        ! Set smaller of the stable and the global time steps
        time_BLK(:,:,:,iBlock)= &
             min(dt,time_BLK(:,:,:,iBlock))
     else                                             !^CFG END IMPLICIT
        time_BLK(:,:,:,iBlock) = dt
     end if                                           !^CFG IF IMPLICIT

     !\
     ! Reset time step to zero inside body.
     !/
     if(.not.true_BLK(iBlock))then
        where(.not.true_cell(1:nI,1:nJ,1:nK,iBlock)) &
             time_BLK(:,:,:,iBlock) = 0.00
     end if

  end do

  ! Set global time step to the actual time step used
  dt = cfl*dt

  if(oktest_me)write(*,*)'Finished set_global_timestep with dt=',dt

end subroutine set_global_timestep
