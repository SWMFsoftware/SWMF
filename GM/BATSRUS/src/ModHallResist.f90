!^CFG COPYRIGHT UM
module ModHallResist

  implicit none

  ! Logical for adding hall resistivity
  logical:: UseHallResist=.false.
  logical:: IsNewBlockHall=.true.
  ! Coefficient for taking whistler wave speed into account
  real:: HallCmaxFactor = 1.0

  ! Non-diagonal part (Hall) resistivity with an arbitrary factor
  real:: IonMassPerCharge, HallFactor

  ! Arrays for the implicit preconditioning
  real, allocatable :: HallJ_CD(:,:,:,:)
  real, allocatable :: BxPerRho_G(:,:,:),ByPerRho_G(:,:,:),BzPerRho_G(:,:,:)

  character(len=20) :: NameHallRegion ='all'
  real :: R1Hall=0.0, R2Hall=1000.0, HallWidth=10.0 

  save

contains
  !============================================================================
  subroutine init_hall_resist
    use ModSize,    ONLY: nI, nJ, nK, nDim
    use ModConst,   ONLY: cProtonMass, cElectronCharge, cMu
    use ModPhysics, ONLY: UnitSI_B, UnitSI_T, UnitSI_X, UnitSI_Rho

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='init_hall_resist'
    !--------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    IonMassPerCharge =HallFactor/cMu*(cProtonMass/cElectronCharge) &
         * UnitSI_B*UnitSI_T/(UnitSI_X**2 * UnitSI_Rho)

    if (DoTestMe) then
       write(*,*) ''
       write(*,*) '>>>>>>>>>>>>>>>>> HALL Resistivity Parameters <<<<<<<<<<'
       write(*,*)
       write(*,*) 'HallFactor       = ',HallFactor
       write(*,*) 'HallCmaxFactor   = ',HallCmaxFactor
       write(*,*) 'IonMassPerCharge = ',IonMassPerCharge
       ! Omega_Bi=B0/IonMassPerCharge'
       write(*,*)
       write(*,*) '>>>>>>>>>>>>>>>>>                       <<<<<<<<<<<<<<<<<'
       write(*,*) ''
    end if

    if(.not.allocated(HallJ_CD)) allocate(&
         HallJ_CD(nI,nJ,nK,nDim), &
         BxPerRho_G(0:nI+1,0:nJ+1,0:nK+1),&
         ByPerRho_G(0:nI+1,0:nJ+1,0:nK+1),&
         BzPerRho_G(0:nI+1,0:nJ+1,0:nK+1) )

  end subroutine init_hall_resist

  !============================================================================

  subroutine get_face_current(iDir, i, j, k, iBlock, Jx, Jy, Jz)

    use ModProcMH,  ONLY: iProc
    use ModAdvance, ONLY: State_VGB, B_, Bx_, By_, Bz_
    use ModMain,    ONLY: nBlock, nI, nJ, nK, x_, y_, z_, &
         iTest,jTest,kTest,VarTest,BlkTest,ProcTest, &
         n_step, nStage
    use ModGeometry,ONLY: Dx_BLK, Dy_BLK, Dz_BLK, &
         X_BLK, y_BLK, z_BLK
    use ModParallel, ONLY:neiLeast,neiLwest,neiLsouth, &
         neiLnorth,neiLtop,neiLbot,BlkNeighborLev
    implicit none

    real, parameter :: cTwoThird = 2.0/3.0, cFourThird = 4.0/3.0

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='get_face_current'

    integer, intent(in):: iDir, i, j, k, iBlock
    real, intent(out)  :: Jx, Jy, Jz

    integer :: iL, iR, jL, jR, kL, kR
    real :: InvDx, InvDy, InvDz


    integer :: i1, j1, k1, i2, j2, k2, iC, jC, kC, iSide, jSide, kSide
    real, save :: b_DG(3,-1:nI+2,-1:nJ+2,-1:nK+2)
    real       :: B1_DG(3,-1:nI+2,-1:nJ+2,-1:nK+2)

    real :: Bc_D(3) ! interpolated coarse cell B field

    real,parameter :: C1 = 8./15., F1 = 2./3., F2 = -1./5.
    real,parameter :: p0 = 5./32., m0 =-3./32., c0= 15./16.
    real :: Ax, Bx, Cx, Ay, By, Cy, Az, Bz, Cz
    integer:: ip, im, jp, jm, kp, km, iDim

    real :: jLeft_D(3), jRight_D(3)
    !-------------------------------------------------------------------------
    if(iProc==ProcTest.and.iBlock==BlkTest.and. &
             (i==iTest.or.i==iTest+1.and.iDir==x_) .and. &
             (j==jTest.or.j==jTest+1.and.iDir==y_) .and. &
             (k==kTest.or.k==kTest+1.and.iDir==z_))then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTestMe = .false.; DoTest =.false.
    end if

    if( IsNewBlockHall ) then
       IsNewBlockHall = .false.

       ! Copy State_VGB into local array and overwrite ghost cells
       B1_DG = State_VGB(Bx_:Bz_,:,:,:,iBlock)
       b_DG  = b1_DG

       ! Fix ghost edge and corner ghost cells
       do kSide = -1,1; do jSide = -1,1; do iSide = -1,1
          ! If the corner or edge is not coarser but the face is
          ! then average out the 8 fine cells so that the
          ! general interpolation formulas remain 2nd order
          ! accurate
          if( BlkNeighborLev(iSide,jSide,kSide,iBlock) /= 1 .and. ( &
               BlkNeighborLev(iSide,0,0,iBlock) == 1 .or. &
               BlkNeighborLev(0,jSide,0,iBlock) == 1 .or. &
               BlkNeighborLev(0,0,kSide,iBlock) == 1)) then

             if(iSide==0)then
                iL = 1; iR = nI
             elseif(iSide==1)then
                iL = nI+1; iR=nI+2
             else
                iL = -1; iR = 0
             end if
             if(jSide==0)then
                jL = 1; jR = nJ
             elseif(jSide==1)then
                jL = nJ+1; jR=nJ+2
             else
                jL = -1; jR = 0
             end if
             if(kSide==0)then
                kL = 1; kR = nK
             elseif(kSide==1)then
                kL = nK+1; kR=nK+2
             else
                kL = -1; kR = 0
             end if

             do k1=kL,kR,2; do j1=jL,jR,2; do i1=iL,iR,2; do iDim=1,3
                B1_DG(iDim,i1:i1+1,j1:j1+1,k1:k1+1)= &
                     0.125*sum(b_DG(iDim,i1:i1+1,j1:j1+1,k1:k1+1))
             end do; end do; end do; end do

          end if
       end do; end do; end do

       ! Do six faces
       if(NeiLeast(iBlock) == 1)then
          do k1=1, nK, 2; do j1=1, nJ, 2; do k2 = k1,k1+1; do j2 = j1,j1+1
             jp = 3*j2 - 2*j1 -1 ; jm = 4*j1 -3*j2 +2
             kp = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
             Bc_D = c0*B1_DG(:,0,j2,k2) &
                  + p0*B1_DG(:,0,jp,kp) &
                  + m0*B1_DG(:,0,jm,km)
             b_DG(:,0,j2,k2) = &
                  C1*Bc_D + F1*b_DG(:,1,j2,k2) + F2*b_DG(:,2,j2,k2)
             !b_DG(:,0,j2,k2) = min( max( b_DG(:,0,j2,k2), &
             !     min(Bc_D, b_DG(:,1,j2,k2),b_DG(:,2,j2,k2))), &
             !     max(Bc_D, b_DG(:,1,j2,k2),b_DG(:,2,j2,k2)))
          end do; end do; end do; end do
       end if

       if(NeiLwest(iBlock) == 1)then
          do k1=1, nK, 2; do j1=1, nJ, 2; do k2 = k1,k1+1; do j2 = j1,j1+1
             jp = 3*j2 - 2*j1 -1 ; jm = 4*j1 -3*j2 +2
             kp = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
             Bc_D = c0*B1_DG(:,nI+1,j2,k2) &
                  + p0*B1_DG(:,nI+1,jp,kp) &
                  + m0*B1_DG(:,nI+1,jm,km)
             b_DG(:,nI+1,j2,k2)= &
                  C1*Bc_D + F1*b_DG(:,nI,j2,k2) + F2*b_DG(:,nI-1,j2,k2)
             !b_DG(:,nI+1,j2,k2) = min( max( b_DG(:,nI+1,j2,k2), &
             !     min(Bc_D, b_DG(:,nI,j2,k2),b_DG(:,nI-1,j2,k2))), &
             !     max(Bc_D, b_DG(:,nI,j2,k2),b_DG(:,nI-1,j2,k2)))
          end do; end do; end do; end do
       end if

       if(NeiLsouth(iBlock) == 1)then
          do k1=1, nK, 2; do i1=1, nI, 2; do k2 = k1,k1+1; do i2 = i1,i1+1
             ip = 3*i2 - 2*i1 -1 ; im = 4*i1 -3*i2 +2
             kp = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
             Bc_D = c0*B1_DG(:,i2,0,k2) &
                  + p0*B1_DG(:,ip,0,kp) &
                  + m0*B1_DG(:,im,0,km)
             b_DG(:,i2,0,k2) = &
                  C1*Bc_D + F1*b_DG(:,i2,1,k2) + F2*b_DG(:,i2,2,k2)
             !b_DG(:,i2,0,k2) = min(max( b_DG(:,i2,0,k2), &
             !     min(Bc_D, b_DG(:,i2,1,k2),b_DG(:,i2,2,k2))), &
             !     max(Bc_D, b_DG(:,i2,1,k2),b_DG(:,i2,2,k2)))
          end do; end do; end do; end do
       end if

       if(NeiLnorth(iBlock) == 1)then
          do k1=1, nK, 2; do i1=1, nI, 2; do k2 = k1,k1+1; do i2 = i1,i1+1
             ip = 3*i2 - 2*i1 -1 ; im = 4*i1 -3*i2 +2
             kp = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
             Bc_D = c0*B1_DG(:,i2,nJ+1,k2) &
                  + p0*B1_DG(:,ip,nJ+1,kp) &
                  + m0*B1_DG(:,im,nJ+1,km)
             b_DG(:,i2,nJ+1,k2) = &
                  C1*Bc_D + F1*b_DG(:,i2,nJ,k2) + F2*b_DG(:,i2,nJ-1,k2)
             !b_DG(:,i2,nJ+1,k2) = min(max( b_DG(:,i2,nJ+1,k2), &
             !     min(Bc_D, b_DG(:,i2,nJ,k2),b_DG(:,i2,nJ-1,k2))), &
             !     max(Bc_D, b_DG(:,i2,nJ,k2),b_DG(:,i2,nJ-1,k2)))
          end do; end do; end do; end do
       end if

       if(NeiLbot(iBlock) == 1)then
          do j1=1, nJ, 2; do i1=1, nI, 2; do j2 = j1,j1+1; do i2 = i1,i1+1
             ip = 3*i2 - 2*i1 -1 ; im = 4*i1 -3*i2 +2
             jp = 3*j2 - 2*j1 -1 ; jm = 4*j1 -3*j2 +2
             Bc_D = c0*B1_DG(:,i2,j2,0) &
                  + p0*B1_DG(:,ip,jp,0) &
                  + m0*B1_DG(:,im,jm,0)
             b_DG(:,i2,j2,0) = &
                  C1*Bc_D + F1*b_DG(:,i2,j2,1) + F2*b_DG(:,i2,j2,2)
             !b_DG(:,i2,j2,0) = min(max( b_DG(:,i2,j2,0), &
             !     min(Bc_D, b_DG(:,i2,j2,1), b_DG(:,i2,j2,2))), &
             !     max(Bc_D, b_DG(:,i2,j2,1), b_DG(:,i2,j2,2)))
          end do; end do; end do; end do
       end if

       if(NeiLtop(iBlock) == 1)then
          do j1=1, nJ, 2; do i1=1, nI, 2; do j2 = j1,j1+1; do i2 = i1,i1+1
             ip = 3*i2 - 2*i1 -1 ; im = 4*i1 -3*i2 +2
             jp = 3*j2 - 2*j1 -1 ; jm = 4*j1 -3*j2 +2
             Bc_D = c0*B1_DG(:,i2,j2,nK+1) &
                  + p0*B1_DG(:,ip,jp,nK+1) &
                  + m0*B1_DG(:,im,jm,nK+1)
             b_DG(:,i2,j2,nK+1) = &
                  C1*Bc_D + F1*b_DG(:,i2,j2,nK) + F2*b_DG(:,i2,j2,nK-1)
             !b_DG(:,i2,j2,nK+1) = min(max( b_DG(:,i2,j2,nK+1), &
             !     min(Bc_D, b_DG(:,i2,j2,nK), b_DG(:,i2,j2,nK-1))), &
             !     max(Bc_D, b_DG(:,i2,j2,nK), b_DG(:,i2,j2,nK-1)))
          end do; end do; end do; end do
       end if

       ! Do 12 edges
       ! 4 X edges
       do kSide = -1,1,2; do jSide = -1,1,2
          if(BlkNeighborLev(0, jSide, kSide, iBlock) /= 1)CYCLE
          j1=1; if(jSide==1) j1=nJ; j2 = j1-jSide; jC = j1+jSide
          k1=1; if(kSide==1) k1=nK; k2 = k1-kSide; kC = k1+kSide
          do i1 = 1,nI,2; do i2 = i1, i1+1
             ip = 3*i2 - 2*i1 -1 ; im = 4*i1 -3*i2 +2
             Bc_D = c0*B1_DG(:,i2,jC,kC) &
                  + p0*B1_DG(:,ip,jC,kC) &
                  + m0*B1_DG(:,im,jC,kC)
             b_DG(:,i2,jC,kC) = &
                  C1* Bc_D + F1*b_DG(:,i2,j1,k1) + F2*b_DG(:,i2,j2,k2)
             !b_DG(:,i2,jC,kC) = min(max( b_DG(:,i2,jC,kC), &
             !     min(Bc_D, b_DG(:,i2,j1,k1), b_DG(:,i2,j2,k2))), &
             !     max(Bc_D, b_DG(:,i2,j1,k1), b_DG(:,i2,j2,k2)))
          end do;end do
       end do;end do
       ! 4 Y edges
       do kSide = -1,1,2; do iSide = -1,1,2
          if(BlkNeighborLev(iSide, 0, kSide, iBlock) /= 1)CYCLE
          i1=1; if(iSide==1) i1=nI; i2 = i1-iSide; iC = i1+iSide
          k1=1; if(kSide==1) k1=nK; k2 = k1-kSide; kC = k1+kSide
          do j1 = 1, nJ, 2; do j2 = j1, j1+1
             jp = 3*j2 - 2*j1 -1 ; jm = 4*j1 -3*j2 +2
             Bc_D = c0*B1_DG(:,iC,j2,kC) &
                  + p0*B1_DG(:,iC,jp,kC) &
                  + m0*B1_DG(:,iC,jm,kC)
             b_DG(:,iC,j2,kC) = &
                  C1*Bc_D + F1*b_DG(:,i1,j2,k1) + F2*b_DG(:,i2,j2,k2)
             !b_DG(:,iC,j2,kC) = min(max( b_DG(:,iC,j2,kC), &
             !     min(Bc_D, b_DG(:,i1,j2,k1), b_DG(:,i2,j2,k2))), &
             !     max(Bc_D, b_DG(:,i1,j2,k1), b_DG(:,i2,j2,k2)))
          end do;end do
       end do;end do
       ! 4 Z edges
       do jSide = -1,1,2; do iSide = -1,1,2
          if(BlkNeighborLev(iSide, jSide, 0, iBlock) /= 1)CYCLE
          i1=1; if(iSide==1) i1=nI; i2 = i1-iSide; iC = i1+iSide
          j1=1; if(jSide==1) j1=nJ; j2 = j1-jSide; jC = j1+jSide
          do k1 = 1, nK, 2 ; do k2 = k1, k1 + 1
             kp = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
             Bc_D = c0*B1_DG(:,iC,jC,k2) &
                  + p0*B1_DG(:,iC,jC,kp) &
                  + m0*B1_DG(:,iC,jC,km)
             b_DG(:,iC,jC,k2) = &
                  C1*Bc_D + F1*b_DG(:,i1,j1,k2) + F2*b_DG(:,i2,j2,k2)
             !b_DG(:,iC,jC,k2) = min(max( b_DG(:,iC,jC,k2), &
             !     min(Bc_D, b_DG(:,i1,j1,k2), b_DG(:,i2,j2,k2))), &
             !     max(Bc_D, b_DG(:,i1,j1,k2), b_DG(:,i2,j2,k2)))
          end do;end do         
       end do;end do

    end if

    InvDx = 1.0/dx_Blk(iBlock)
    InvDy = 1.0/dy_Blk(iBlock)
    InvDz = 1.0/dz_Blk(iBlock)

    ! Central difference with averaging in orthogonal direction
    iR = i+1; iL = i-1; 
    jR = j+1; jL = j-1; 
    kR = k+1; kL = k-1; 

    Ax = -0.25*InvDx; Bx = 0.0; Cx = +0.25*InvDx
    Ay = -0.25*InvDy; By = 0.0; Cy = +0.25*InvDy
    Az = -0.25*InvDz; Bz = 0.0; Cz = +0.25*InvDz

    if(i==1)then
       if(NeiLeast(iBlock)==-1 &
            .or. (iDir==y_ .and. &
            (j==1    .and. BlkNeighborLev(-1,-1,0, iBlock)==-1) .or. &
            (j==nJ+1 .and. BlkNeighborLev(-1, 1,0, iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. BlkNeighborLev(-1, 0,-1, iBlock)==-1) .or. &
            (k==nK+1 .and. BlkNeighborLev(-1, 0, 1, iBlock)==-1))&
            ) then
          iL = i+1; iR = i+2; Ax=InvDx; Bx=-0.75*InvDx; Cx=-0.25*InvDx
       end if
    elseif(i==nI)then
       if(NeiLwest(iBlock)==-1 &
            .or. (iDir==y_ .and. &
            (j==1    .and. BlkNeighborLev( 1,-1, 0, iBlock)==-1) .or. &
            (j==nJ+1 .and. BlkNeighborLev( 1, 1, 0, iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. BlkNeighborLev( 1, 0,-1, iBlock)==-1) .or. &
            (k==nK+1 .and. BlkNeighborLev( 1, 0, 1, iBlock)==-1))&
            ) then
          iL = i-1; iR = i-2; Ax=-InvDx; Bx=0.75*InvDx; Cx=0.25*InvDx
       end if
    end if

    if(j==1)then
       if(NeiLsouth(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. BlkNeighborLev(-1,-1,0, iBlock)==-1) .or. &
            (i==nI+1 .and. BlkNeighborLev( 1,-1,0, iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. BlkNeighborLev( 0,-1,-1, iBlock)==-1) .or. &
            (k==nK+1 .and. BlkNeighborLev( 0,-1, 1, iBlock)==-1))&
            )then
          jL = j+1; jR = j+2; Ay=InvDy; By=-0.75*InvDy; Cy=-0.25*InvDy
       end if
    elseif(j==nJ)then
       if(NeiLnorth(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. BlkNeighborLev(-1, 1,0, iBlock)==-1) .or. &
            (i==nI+1 .and. BlkNeighborLev( 1, 1,0, iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. BlkNeighborLev( 0, 1,-1, iBlock)==-1) .or. &
            (k==nK+1 .and. BlkNeighborLev( 0, 1, 1, iBlock)==-1))&
            )then
          jL = j-1; jR = j-2; Ay=-InvDy; By=0.75*InvDy; Cy=0.25*InvDy
       end if
    end if

    if(k==1)then
       if(NeiLbot(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. BlkNeighborLev(-1,0,-1, iBlock)==-1) .or. &
            (i==nI+1 .and. BlkNeighborLev( 1,0,-1, iBlock)==-1)) &
            .or. (iDir==y_ .and. &
            (j==1    .and. BlkNeighborLev( 0,-1,-1, iBlock)==-1) .or. &
            (j==nJ+1 .and. BlkNeighborLev( 0, 1,-1, iBlock)==-1))&
            )then
          kL = k+1; kR = k+2; Az=InvDz; Bz=-0.75*InvDz; Cz=-0.25*InvDz
       end if
    elseif(k==nK)then
       if(NeiLtop(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. BlkNeighborLev(-1,0, 1, iBlock)==-1) .or. &
            (i==nI+1 .and. BlkNeighborLev( 1,0, 1, iBlock)==-1)) &
            .or. (iDir==y_ .and. &
            (j==1    .and. BlkNeighborLev( 0,-1,1, iBlock)==-1) .or. &
            (j==nJ+1 .and. BlkNeighborLev( 0, 1,1, iBlock)==-1))&
            )then
          kL = k-1; kR = k-2; Az=-InvDz; Bz=0.75*InvDz; Cz=0.25*InvDz
       end if
    end if

    select case(iDir)
    case(x_)
       Jy = -InvDx* (b_DG(z_,i  ,j,k) - b_DG(z_,i-1,j,k)) &
            + Az*(b_DG(x_,i-1,j,kL)+b_DG(x_,i,j ,kL))     &
            + Bz*(b_DG(x_,i-1,j,k )+b_DG(x_,i,j ,k ))     &
            + Cz*(b_DG(x_,i-1,j,kR)+b_DG(x_,i,j ,kR))

       Jz = +InvDx* (b_DG(y_,i,j,k) - b_DG(y_,i-1,j,k)) &
            - Ay*(b_DG(x_,i-1,jL,k)+b_DG(x_,i,jL,k)) &
            - By*(b_DG(x_,i-1,j ,k)+b_DG(x_,i,j ,k)) &
            - Cy*(b_DG(x_,i-1,jR,k)+b_DG(x_,i,jR,k)) 

       Jx = + Ay*(b_DG(z_,i-1,jL,k)+b_DG(z_,i,jL,k )) &
            + By*(b_DG(z_,i-1,j ,k)+b_DG(z_,i,j ,k )) &
            + Cy*(b_DG(z_,i-1,jR,k)+b_DG(z_,i,jR,k )) &
            - Az*(b_DG(y_,i-1,j,kL)+b_DG(y_,i,j ,kL)) &
            - Bz*(b_DG(y_,i-1,j,k )+b_DG(y_,i,j ,k )) &
            - Cz*(b_DG(y_,i-1,j,kR)+b_DG(y_,i,j ,kR))

       !call get_current(i-1,j,k,iBlock,jLeft_D)

    case(y_)
       Jx = + InvDy*(b_DG(z_,i,j,k) - b_DG(z_,i,j-1,k)) &
            - Az*(b_DG(y_,i,j-1,kL) + b_DG(y_,i,j ,kL)) &
            - Bz*(b_DG(y_,i,j-1,k ) + b_DG(y_,i,j ,k )) &
            - Cz*(b_DG(y_,i,j-1,kR) + b_DG(y_,i,j ,kR))

       Jz = - InvDy*(b_DG(x_,i,j,k) - b_DG(x_,i,j-1,k)) &
            + Ax*(b_DG(y_,iL,j-1,k) + b_DG(y_,iL,j ,k)) &
            + Bx*(b_DG(y_,i ,j-1,k) + b_DG(y_,i ,j ,k)) &
            + Cx*(b_DG(y_,iR,j-1,k) + b_DG(y_,iR,j ,k))

       Jy = + Az*(b_DG(x_,i,j-1,kL) + b_DG(x_,i,j ,kL)) &
            + Bz*(b_DG(x_,i,j-1,k ) + b_DG(x_,i,j ,k )) &
            + Cz*(b_DG(x_,i,j-1,kR) + b_DG(x_,i,j ,kR)) &
            - Ax*(b_DG(z_,iL,j-1,k) + b_DG(z_,iL,j ,k)) &
            - Bx*(b_DG(z_,i ,j-1,k) + b_DG(z_,i ,j ,k)) &
            - Cx*(b_DG(z_,iR,j-1,k) + b_DG(z_,iR,j ,k))

       !call get_current(i,j-1,k,iBlock,jLeft_D)

    case(z_)

       Jx = -InvDz* (b_DG(y_,i,j,k) - b_DG(y_,i,j,k-1)) & 
            + Ay*(b_DG(z_,i,jL,k-1) + b_DG(z_,i,jL,k))  &
            + By*(b_DG(z_,i,j ,k-1) + b_DG(z_,i,j ,k))  &
            + Cy*(b_DG(z_,i,jR,k-1) + b_DG(z_,i,jR,k)) 

       Jy = +InvDz* (b_DG(x_,i,j,k) - b_DG(x_,i,j,k-1)) &
            - Ax*(b_DG(z_,iL,j,k-1) + b_DG(z_,iL,j,k))  &
            - Bx*(b_DG(z_,i ,j,k-1) + b_DG(z_,i ,j,k))  &
            - Cx*(b_DG(z_,iR,j,k-1) + b_DG(z_,iR,j,k)) 

       Jz = + Ax*(b_DG(y_,iL,j,k-1) + b_DG(y_,iL,j,k))  &
            + Bx*(b_DG(y_,i ,j,k-1) + b_DG(y_,i ,j,k))  &
            + Cx*(b_DG(y_,iR,j,k-1) + b_DG(y_,iR,j,k))  &
            - Ay*(b_DG(x_,i,jL,k-1) + b_DG(x_,i,jL,k))  &
            - By*(b_DG(x_,i,j ,k-1) + b_DG(x_,i,j ,k))  &
            - Cy*(b_DG(x_,i,jR,k-1) + b_DG(x_,i,jR,k)) 

       !call get_current(i,j,k-1,iBlock,jLeft_D)

    case default
       write(*,*)'Error in get_face_current: iDir=',iDir
       call stop_mpi('DEBUG')
    end select

    !call get_current(i,j,k,iBlock,jRight_D)
    !
    !Jx = min(Jx, max(jLeft_D(x_),jRight_D(x_)))
    !Jy = min(Jy, max(jLeft_D(y_),jRight_D(y_)))
    !Jz = min(Jz, max(jLeft_D(z_),jRight_D(z_)))
    !
    !Jx = max(Jx, min(jLeft_D(x_),jRight_D(x_)))
    !Jy = max(Jy, min(jLeft_D(y_),jRight_D(y_)))
    !Jz = max(Jz, min(jLeft_D(z_),jRight_D(z_)))

    if(DoTestMe)then
       write(*,*)NameSub,': Jx=',Jx
       write(*,*)NameSub,': Jy=',Jy
       write(*,*)NameSub,': Jz=',Jz
    end if

  end subroutine get_face_current

  !==========================================================================
  subroutine test_face_current

    use ModMain,     ONLY: nI, nJ, nK, nBlock, UnusedBlk, x_, y_, z_, &
         east_, west_, south_, north_, bot_, top_
    use ModAdvance,  ONLY: State_VGB, Bx_, By_, Bz_
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, dx_BLK, dy_BLK, dz_BLK
    use ModParallel, ONLY: NeiLev
    use ModFaceValue,ONLY: correct_monotone_restrict

    integer, parameter :: nTest = 2
    integer :: i,j,k,iBlock,iTest
    real :: Jx, Jy, Jz
    !------------------------------------------------------------------------

    write(*,*)'test_face_current starting !!!'

    do iTest = 1, nTest

       do iBlock = 1, nBlock
          if(UnusedBlk(iBlock)) CYCLE

          select case(iTest)
          case(1)
             State_VGB(Bx_,:,:,:,iBlock) = &
                  + 1*x_BLK(:,:,:,iBlock) &
                  + 2*y_BLK(:,:,:,iBlock) &
                  + 3*z_BLK(:,:,:,iBlock) &
                  + 4*x_BLK(:,:,:,iBlock)*y_BLK(:,:,:,iBlock) &
                  + 5*x_BLK(:,:,:,iBlock)*z_BLK(:,:,:,iBlock) &
                  + 6*y_BLK(:,:,:,iBlock)*z_BLK(:,:,:,iBlock)
             State_VGB(By_,:,:,:,iBlock) = &
                  + 10*x_BLK(:,:,:,iBlock) &
                  + 20*y_BLK(:,:,:,iBlock) &
                  + 30*z_BLK(:,:,:,iBlock) &
                  + 40*x_BLK(:,:,:,iBlock)*y_BLK(:,:,:,iBlock) &
                  + 50*x_BLK(:,:,:,iBlock)*z_BLK(:,:,:,iBlock) &
                  + 60*y_BLK(:,:,:,iBlock)*z_BLK(:,:,:,iBlock)
             State_VGB(Bz_,:,:,:,iBlock) = &
                  + 100*x_BLK(:,:,:,iBlock) &
                  + 200*y_BLK(:,:,:,iBlock) &
                  + 300*z_BLK(:,:,:,iBlock) &
                  + 400*x_BLK(:,:,:,iBlock)*y_BLK(:,:,:,iBlock) &
                  + 500*x_BLK(:,:,:,iBlock)*z_BLK(:,:,:,iBlock) &
                  + 600*y_BLK(:,:,:,iBlock)*z_BLK(:,:,:,iBlock)
          case(2)
             State_VGB(Bx_,:,:,:,iBlock) = 1.0 + &
                  0.01*x_BLK(:,:,:,iBlock)**2 + &
                  0.02*y_BLK(:,:,:,iBlock)**2 + &
                  0.03*z_BLK(:,:,:,iBlock)**2
             State_VGB(By_,:,:,:,iBlock) = 10.0 + &
                  0.1*x_BLK(:,:,:,iBlock)**2 + &
                  0.2*y_BLK(:,:,:,iBlock)**2 + &
                  0.3*z_BLK(:,:,:,iBlock)**2
             State_VGB(Bz_,:,:,:,iBlock) = 100.0 + &
                  1.0*x_BLK(:,:,:,iBlock)**2 + &
                  2.0*y_BLK(:,:,:,iBlock)**2 + &
                  3.0*z_BLK(:,:,:,iBlock)**2
          end select

       end do

       call message_pass_cells_8state(.false., .false., .true.)

       do iBlock = 1, nBlock
          if(UnusedBlk(iBlock)) CYCLE

          call correct_monotone_restrict(iBlock)

          do k=1, nK; do j=1,nJ; do i=1,nI+1
             call get_face_current(x_, i, j, k, iBlock, Jx, Jy, Jz)
             call check_error('x')
          end do; end do; end do

          do k=1, nK; do j=1,nJ+1; do i=1,nI
             call get_face_current(y_, i, j, k, iBlock, Jx, Jy, Jz)
             call check_error('y')
          end do; end do; end do

          do k=1, nK+1; do j=1,nJ; do i=1,nI
             call get_face_current(z_, i, j, k, iBlock, Jx, Jy, Jz)
             call check_error('z')
          end do; end do; end do

       end do

       write(*,*)'test_face_current: test ',iTest,' passed !!!'

    end do

    call stop_mpi('test_face_current succeeded !!! ')

  contains

    subroutine check_error(NameDir)
      character, intent(in):: NameDir
      real :: x, y, z, JxGood, JyGood, JzGood, Tolerance

      !--------------------------------------------------------------------
      ! Face center coordinates
      x = x_BLK(i,j,k,iBlock)
      y = y_BLK(i,j,k,iBlock)
      z = z_BLK(i,j,k,iBlock)

      select case(NameDir)
      case('x')
         if(i==1   .and.neiLEV(east_,iBlock)==-1) RETURN
         if(i==nI+1.and.neiLEV(west_,iBlock)==-1) RETURN
         x = x - 0.5*dx_BLK(iBlock)
      case('y')
         if(j==1   .and.neiLEV(south_,iBlock)==-1) RETURN
         if(j==nJ+1.and.neiLEV(north_,iBlock)==-1) RETURN
         y = y - 0.5*dy_BLK(iBlock)
      case('z')
         if(k==1   .and.neiLEV(bot_,iBlock)==-1) RETURN
         if(k==nK+1.and.neiLEV(top_,iBlock)==-1) RETURN
         z = z - 0.5*dz_BLK(iBlock)
      end select

      select case(iTest)
      case(1)
         JxGood = 200 + 400*x + 600*z -  30 -  50*x -  60*y
         JyGood =   3 +   5*x +   6*y - 100 - 400*y - 500*z
         JzGood =  10 +  40*y +  50*z -   2 -   4*x -   6*z

         Tolerance = 1e-6
      case(2)
         JxGood = 4.0 *y - 0.6 *z
         JyGood = 0.06*z - 2.0 *x
         JzGood = 0.2 *x - 0.04*y

         Tolerance = 5e-2
      end select

      if(       abs(Jx-JxGood) > Tolerance*abs(JxGood) + 1e-6  &
           .or. abs(Jy-JyGood) > Tolerance*abs(JyGood) + 1e-6  &
           .or. abs(Jz-JzGood) > Tolerance*abs(JzGood) + 1e-6) then

         write(*,*)'Face=',NameDir
         write(*,*)'iTest, i,j,k,iBlock=',iTest, i,j,k,iBlock
         write(*,*)'x,y,z=', &
              x_BLK(i,j,k,iBlock), y_BLK(i,j,k,iBlock), z_BLK(i,j,k,iBlock)
         write(*,*)'dx,dy,dz=',&
              dx_BLK(iBlock), dy_BLK(iBlock), dz_BLK(iBlock)
         write(*,*)'Bad  Jx,Jy,Jz =',Jx,Jy,Jz
         write(*,*)'Good Jx,Jy,Jz =',JxGood,JyGood,JzGood
         write(*,*)'NeiLev=',NeiLev(:,iBlock)
         call stop_mpi('Error')
      end if

    end subroutine check_error

  end subroutine test_face_current
  
  real function hall_factor(iDir, iFace, jFace, kFace , iBlock)
    use ModMain,     ONLY: nI, nJ, nK, nBlock
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, dx_BLK, dy_BLK, dz_BLK
    use ModPhysics, ONLY : Rbody
    integer, intent(in)::iDir, iFace, jFace, kFace, iBlock 
    real :: r,x,y,z

    !--------------------------------------------------------------
    select case(iDir)
    case(0)  !for cell center
       x = x_BLK(iFace,jFace,kFace,iBlock)
       y = y_BLK(iFace,jFace,kFace,iBlock)
       z = z_BLK(iFace,jFace,kFace,iBlock)       
    case(1)
       x = 0.5*sum(x_BLK(iFace-1:iFace,jFace,kFace,iBlock))
       y = 0.5*sum(y_BLK(iFace-1:iFace,jFace,kFace,iBlock))
       z = 0.5*sum(z_BLK(iFace-1:iFace,jFace,kFace,iBlock))
    case(2)
       x = 0.5*sum(x_BLK(iFace,jFace-1:jFace,kFace,iBlock))
       y = 0.5*sum(y_BLK(iFace,jFace-1:jFace,kFace,iBlock))
       z = 0.5*sum(z_BLK(iFace,jFace-1:jFace,kFace,iBlock))
    case(3)
       x = 0.5*sum(x_BLK(iFace,jFace,kFace-1:KFace,iBlock))
       y = 0.5*sum(y_BLK(iFace,jFace,kFace-1:KFace,iBlock))
       z = 0.5*sum(z_BLK(iFace,jFace,kFace-1:KFace,iBlock))
    end select
    
    select case(NameHallRegion)
    case("all")
       hall_factor = 1.0
    case("user")
       ! hall_factor= &
       !   user_hall_factor(x, y, z, iDir, iFace, jFace, kFace, iBlock)
    case("shell")
       hall_factor = 1.0
       r=sqrt(x*x+y*y+z*z)       
       if(r < R1Hall .or. r > R2Hall)then
          hall_factor=0.0
       else if(r > R2Hall-HallWidth)then
          hall_factor = (R2Hall-r)/HallWidth
       else if(r < R1Hall+HallWidth)then       
          hall_factor = (r-R1Hall)/HallWidth
       end if
    case default
    end select
  end function hall_factor

end module ModHallResist
