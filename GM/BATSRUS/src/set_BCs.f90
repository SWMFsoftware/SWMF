!^CFG COPYRIGHT UM
subroutine set_BCs(TimeBcIn,DoResChangeOnlyIn)
  use ModProcMH
  use ModMain
  use ModAdvance
  use ModNumConst
  use ModGeometry, ONLY:&
       IsBoundaryCell_GI , IsBoundaryBlock_IB,true_cell,MinBoundary,MaxBoundary
  use ModBoundaryCells
  use ModFaceBc, ONLY: TimeBc, DoResChangeOnly, iBlockBc, iBoundary

  implicit none

  real,     intent(in) :: TimeBcIn
  logical, intent (in) :: DoResChangeOnlyIn
  logical :: oktest, oktest_me
  character(len=*), parameter:: NameSub='set_bcs'
  !----------------------------------------------------------------------------

  call timing_start(NameSub)

  if(globalBLK==BLKtest.and.iProc==PROCtest)then
     call set_oktest(NameSub, oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  endif

  ! set variables in module
  TimeBc          = TimeBcIn
  DoResChangeOnly = DoResChangeOnlyIn
  iBlockBc        = GlobalBlk


  if(oktest_me)call write_face_state('Initial')

  call set_boundary_cells(iBlockBc)
  if(SaveBoundaryCells)then
     do iBoundary=MinBoundarySaved,MaxBoundarySaved
        IsBoundaryCell_GI(:,:,:,iBoundary)=&
             IsBoundaryCell_IGB(iBoundary,:,:,:,iBlockBc)
     end do
  end if
  
  !\
  ! Apply boundary conditions
  !/
  do iBoundary = MinBoundary, MaxBoundary
     if(IsBoundaryBlock_IB(iBoundary,globalBLK)) call set_face_BCs( &
          IsBoundaryCell_GI(:,:,:,iBoundary),&
          true_cell(:,:,:,globalBLK) )
  end do
  
  if(oktest_me)call write_face_state('Final')

  call timing_stop(NameSub)

contains

  subroutine write_face_state(String)

    character(len=*), intent(in):: String
    
    write(*,*)NameSub,' ',String,' face states:'
    write(*,*)'east  VarL_x,VarR_x=',&
         LeftState_VX(VarTest,  Itest, Jtest, Ktest),  &
         RightState_VX(VarTest, Itest, Jtest, Ktest)
    write(*,*)'west  VarL_x,VarR_x=',&
         LeftState_VX(VarTest,  Itest+1, Jtest, Ktest), &
         RightState_VX(VarTest, Itest+1, Jtest, Ktest)
    write(*,*)'south VarL_y,VarR_y=',&
         LeftState_VY(VarTest,  Itest, Jtest, Ktest),  &
         RightState_VY(VarTest, Itest, Jtest, Ktest)
    write(*,*)'north VarL_y,VarR_y=',&
         LeftState_VY(VarTest,  Itest, Jtest+1, Ktest), &
         RightState_VY(VarTest, Itest, Jtest+1, Ktest)
    write(*,*)'bot   VarL_z,VarR_z=',&
         LeftState_VZ(VarTest,  Itest, Jtest, Ktest), &
         RightState_VZ(VarTest, Itest, Jtest, Ktest)
    write(*,*)'top   VarL_z,VarR_z=',&
         LeftState_VZ(VarTest,  Itest, Jtest, Ktest+1), &
         RightState_VZ(VarTest, Itest, Jtest, Ktest+1)

  end subroutine write_face_state

end subroutine set_BCs

!=============================================================================

subroutine set_face_BCs(IsBodyCell_G,IsTrueCell_G)

  use ModMain
  use ModGeometry, ONLY : x_BLK, y_BLK, z_BLK
  use ModAdvance
  use ModParallel, ONLY : neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth
  use ModNumConst
  use ModFaceBc

  implicit none

  logical, dimension(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn), intent(in) :: &
       IsBodyCell_G, IsTrueCell_G 

  integer :: i,j,k

  character (len=*), parameter :: NameSub = 'set_face_BCs'
  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------
  if(iBlockBc==BLKtest.and.iProc==PROCtest)then
     call set_oktest(NameSub, DoTest, DoTestMe)
  else
     DoTest = .false.; DoTestMe = .false.
  end if

  TypeBc = TypeBc_I(iBoundary)

  !\
  ! Apply body BCs as required.
  !/                            
                   
  do k = kMinFaceX, kMaxFaceX
     do j = jMinFaceX, jMaxFaceX
        do i = 1, nIFace
           !\
           ! Apply BCs at X-direction faces as necessary.
           !/
           !====================================
           !         NUMBERING!
           !
           !     C     F     C  B  F     C
           !     +     |     +  !  |     +
           !
           !     i-1   i     i     i+1   i+1
           !
           !====================================
           if (IsTrueCell_G(i-1,j,k) .and. &
                IsBodyCell_G(i,j,k) .and. &
                (.not.DoResChangeOnly .or. &
                ((i == nIFace .and. neiLwest(iBlockBc)==+1) .or. &
                (i == 1       .and. neiLeast(iBlockBc)==+1)) )) then

              iSide = West_

              FaceCoords_D(x_) = 0.5*sum(x_BLK(i-1:i,j,k,iBlockBc))
              FaceCoords_D(y_) = 0.5*sum(y_BLK(i-1:i,j,k,iBlockBc))
              FaceCoords_D(z_) = 0.5*sum(z_BLK(i-1:i,j,k,iBlockBc))
              B0Face_D(x_) = B0xFace_x_BLK(i,j,k,iBlockBc)
              B0Face_D(y_) = B0yFace_x_BLK(i,j,k,iBlockBc)
              B0Face_D(z_) = B0zFace_x_BLK(i,j,k,iBlockBc)

              VarsTrueFace_V= LeftState_VX(:,i,j,k)

              call set_face_bc

              RightState_VX(:,i,j,k) = VarsGhostFace_V

           end if

           if (IsTrueCell_G(i,j,k) .and. &
                IsBodyCell_G(i-1,j,k)  .and. &
                (.not.DoResChangeOnly .or. &
                (i == 1         .and. neiLeast(iBlockBc)==+1) .or. &
                (i == nIFace    .and. neiLwest(iBlockBc)==+1)  )) then

              iSide = East_

              FaceCoords_D(x_) = 0.5*sum(x_BLK(i-1:i,j,k,iBlockBc))
              FaceCoords_D(y_) = 0.5*sum(y_BLK(i-1:i,j,k,iBlockBc))
              FaceCoords_D(z_) = 0.5*sum(z_BLK(i-1:i,j,k,iBlockBc))
              B0Face_D(x_)     = B0xFace_x_BLK(i,j,k,iBlockBc)
              B0Face_D(y_)     = B0yFace_x_BLK(i,j,k,iBlockBc)
              B0Face_D(z_)     = B0zFace_x_BLK(i,j,k,iBlockBc)
            
              VarsTrueFace_V = RightState_VX(:,i,j,k)

              call set_face_bc

              LeftState_VX(:,i,j,k) = VarsGhostFace_V
           end if
        end do !end i loop
     end do !end j loop
  end do !end k loop

  do k = kMinFaceY,kMaxFaceY
     do j = 1 , nJFace
        do i = iMinFaceY, iMaxFaceY
           !\
           ! Apply BCs at Y-direction faces as necessary.
           !/
           if (IsTrueCell_G(i,j-1,k) .and. &
                IsBodyCell_G(i,j,k)  .and. &
                ( .not.DoResChangeOnly .or. &
                (j == nJFace .and. neiLnorth(iBlockBc)==+1) .or. &
                (j == 1      .and. neiLsouth(iBlockBc)==+1) )) then

              iSide = North_

              FaceCoords_D(x_) = 0.5*sum(x_BLK(i,j-1:j,k,iBlockBc))
              FaceCoords_D(y_) = 0.5*sum(y_BLK(i,j-1:j,k,iBlockBc))
              FaceCoords_D(z_) = 0.5*sum(z_BLK(i,j-1:j,k,iBlockBc))
              B0Face_D(x_)     = B0xFace_y_BLK(i,j,k,iBlockBc)
              B0Face_D(y_)     = B0yFace_y_BLK(i,j,k,iBlockBc)
              B0Face_D(z_)     = B0zFace_y_BLK(i,j,k,iBlockBc)

              VarsTrueFace_V = LeftState_VY(:,i,j,k)
          
              call set_face_bc

              RightState_VY(:,i,j,k) = VarsGhostFace_V           
           end if

           if (IsTrueCell_G(i,j,k) .and. &
                IsBodyCell_G(i,j-1,k)  .and. &
                (.not.DoResChangeOnly .or. &
                (j ==1       .and. neiLsouth(iBlockBc)==+1) .or. &
                (j == nJFace .and. neiLnorth(iBlockBc)==+1) )) then

              iSide = South_

              FaceCoords_D(x_) = 0.5*sum(x_BLK(i,j-1:j,k,iBlockBc))
              FaceCoords_D(y_) = 0.5*sum(y_BLK(i,j-1:j,k,iBlockBc))
              FaceCoords_D(z_) = 0.5*sum(z_BLK(i,j-1:j,k,iBlockBc))
              B0Face_D(x_)     = B0xFace_y_BLK(i,j,k,iBlockBc)
              B0Face_D(y_)     = B0yFace_y_BLK(i,j,k,iBlockBc)
              B0Face_D(z_)     = B0zFace_y_BLK(i,j,k,iBlockBc)

              VarsTrueFace_V = RightState_VY(:,i,j,k)

              call set_face_bc

              LeftState_VY(:,i,j,k) = VarsGhostFace_V
           end if
        end do !end j loop
     end do !end i loop
  end do !end k loop

  do k = 1, nKFace
     do j = jMinFaceZ,jMaxFaceZ
        do i = iMinFaceZ,iMaxFaceZ
           !\
           ! Apply BCs at Z-direction faces as necessary.
           !/
           if (IsTrueCell_G(i,j,k-1) .and. &
                IsBodyCell_G(i,j,k) .and. &
                (.not.DoResChangeOnly .or. &
                (k == nKFace .and. neiLtop(iBlockBc)==+1) .or. &
                (k == 1       .and. neiLbot(iBlockBc)==+1)) ) then

              iSide = Top_

              FaceCoords_D(x_)= 0.5*sum(x_BLK(i,j,k-1:k,iBlockBc))
              FaceCoords_D(y_)= 0.5*sum(y_BLK(i,j,k-1:k,iBlockBc))
              FaceCoords_D(z_)= 0.5*sum(z_BLK(i,j,k-1:k,iBlockBc))
              B0Face_D(x_)    = B0xFace_z_BLK(i,j,k,iBlockBc)
              B0Face_D(y_)    = B0yFace_z_BLK(i,j,k,iBlockBc)
              B0Face_D(z_)    = B0zFace_z_BLK(i,j,k,iBlockBc)

              VarsTrueFace_V =  LeftState_VZ(:,i,j,k)
        
              call set_face_bc

              RightState_VZ(:,i,j,k) = VarsGhostFace_V
           end if

           if (IsTrueCell_G(i,j,k).and. &
                IsBodyCell_G(i,j,k-1).and. &
                (.not.DoResChangeOnly .or. &
                (k == 1         .and. neiLbot(iBlockBc)==+1) .or. &
                (k == nKFace .and. neiLtop(iBlockBc)==+1))  ) then

              iSide = Bot_

              FaceCoords_D(x_) = 0.5*sum(x_BLK(i,j,k-1:k,iBlockBc))
              FaceCoords_D(y_) = 0.5*sum(y_BLK(i,j,k-1:k,iBlockBc))
              FaceCoords_D(z_) = 0.5*sum(z_BLK(i,j,k-1:k,iBlockBc))
              B0Face_D(x_)     = B0xFace_z_BLK(i,j,k,iBlockBc)
              B0Face_D(y_)     = B0yFace_z_BLK(i,j,k,iBlockBc)
              B0Face_D(z_)     = B0zFace_z_BLK(i,j,k,iBlockBc)

              VarsTrueFace_V =  RightState_VZ(:,i,j,k)
         
              call set_face_bc

              LeftState_VZ(:,i,j,k) = VarsGhostFace_V
          
           end if
        end do !end i loop
     end do !end j loop
  end do !end k loop
 
contains

  subroutine set_face_bc

    use ModPhysics, ONLY : xBody2,yBody2,zBody2 !^CFG IF SECONDBODY
    use ModPhysics, ONLY : FaceState_VI
    use ModUser, ONLY: user_face_bcs
    
    implicit none

    real, dimension(1:3) :: uRot_D, uIono_D

    real:: FaceState_V(nVar)
    real:: PressureJumpLimit=0.0,DensityJumpLimit=0.1    !
    real ::BdotR,BRefl_D(nDim), BOrig_D(nDim)
    real:: BdotU,RInv
    real:: cosTheta,sinTheta,cosPhi,sinPhi
    real:: UrTrue,UtTrue,BpTrue,BrGhost,BtGhost,BpGhost

    !------------------------------------------------------------------------

    ! User defined boundary conditions
    if( index(TypeBc, 'user') > 0 .or. &
         (UseUserInnerBCs .and. iBoundary <= body1_) .or. &
         (UseUserOuterBCs .and. iBoundary >= east_ ) )then
       iFace = i; jFace = j; kFace = k
       call user_face_bcs(VarsGhostFace_V)
       return
    end if

    !^CFG IF SECONDBODY BEGIN
    if(iBoundary==body2_)then
       FaceCoords_D(x_)= FaceCoords_D(x_) - xBody2
       FaceCoords_D(y_)= FaceCoords_D(y_) - yBody2
       FaceCoords_D(z_)= FaceCoords_D(z_) - zBody2
    end if
    !^CFG END SECONDBODY

    FaceState_V = FaceState_VI(:,iBoundary)  
!    if(iBoundary==west_.and. &
     if(  &
         (TypeBc == 'vary'.or.TypeBc == 'inflow') ) then                          
       call get_solar_wind_point(&
            TimeBc,           &
            FaceCoords_D(x_), &
            FaceCoords_D(y_), &
            FaceCoords_D(z_), &
            FaceState_V)
    end if                                                        

    select case(TypeBc) 
    case('linetied','ionospherefloat')
       VarsGhostFace_V=VarsTrueFace_V
       VarsGhostFace_V(iRhoUx_I) = -VarsGhostFace_V(iRhoUx_I)
       VarsGhostFace_V(iRhoUy_I) = -VarsGhostFace_V(iRhoUy_I)
       VarsGhostFace_V(iRhoUz_I) = -VarsGhostFace_V(iRhoUz_I)
    case('float','outflow')
       VarsGhostFace_V=VarsTrueFace_V
    case('heliofloat')
       VarsGhostFace_V = VarsTrueFace_V
       RInv     = cOne/sqrt(dot_product(FaceCoords_D,FaceCoords_D))
       cosTheta = FaceCoords_D(z_)*RInv
       sinTheta = sqrt(FaceCoords_D(x_)**2+FaceCoords_D(y_)**2)*RInv
       cosPhi   = FaceCoords_D(x_)/&
            sqrt(FaceCoords_D(x_)**2+FaceCoords_D(y_)**2+cTolerance**2)
       sinPhi   = FaceCoords_D(y_)/&
            sqrt(FaceCoords_D(x_)**2+FaceCoords_D(y_)**2+cTolerance**2)
       BdotU    = dot_product(VarsTrueFace_V(Bx_:Bz_),        &
                              VarsTrueFace_V(Ux_:Uz_))/       &
                  (dot_product(VarsTrueFace_V(Ux_:Uz_),        &
                              VarsTrueFace_V(Ux_:Uz_))+cTolerance**2)
       UrTrue   = dot_product(VarsTrueFace_V(Ux_:Uz_),        &
                    FaceCoords_D(x_:z_))*RInv
       UtTrue   = ((VarsTrueFace_V(Ux_)*FaceCoords_D(x_)+     &
                    VarsTrueFace_V(Uy_)*FaceCoords_D(y_))*    &
                    FaceCoords_D(z_)-VarsTrueFace_V(Uz_)*     &
                   (FaceCoords_D(x_)**2+FaceCoords_D(y_)**2))/&
               sqrt(FaceCoords_D(x_)**2+FaceCoords_D(y_)**2+  &
                    cTolerance**2)*RInv
       BpTrue   =  (VarsTrueFace_V(By_)*FaceCoords_D(x_)-     &
                    VarsTrueFace_V(Bx_)*FaceCoords_D(y_))/    &
                  ((FaceCoords_D(x_)**2+FaceCoords_D(y_)**2+  &
                    cTolerance**2)*RInv)*sinTheta
       BrGhost  = UrTrue*BdotU; BtGhost = UtTrue*BdotU;
       BpGhost  = BpTrue
       VarsGhostFace_V(Bx_) = BrGhost*FaceCoords_D(x_)*RInv+  &
            BtGhost*cosTheta*cosPhi-BpGhost*sinPhi
       VarsGhostFace_V(By_) = BrGhost*FaceCoords_D(y_)*RInv+  &
            BtGhost*cosTheta*sinPhi+BpGhost*cosPhi
       VarsGhostFace_V(Bz_) = BrGhost*FaceCoords_D(z_)*RInv-  &
            BtGhost*sinTheta
    case('fixed','inflow','vary')
       VarsGhostFace_V    = FaceState_V
       VarsGhostFace_V(Bx_:Bz_)  =  VarsGhostFace_V(Bx_:Bz_) - B0Face_D
    case('reflect','reflectb')
       BOrig_D = VarsTrueFace_V(Bx_:Bz_)
       if(TypeBc == 'reflectb') BOrig_D = BOrig_D + B0Face_D

       select case(iBoundary)                                                 
       case(body1_,body2_)
          BdotR   = 2*dot_product(BOrig_D, FaceCoords_D)/sum(FaceCoords_D**2)
          BRefl_D = FaceCoords_D*BdotR
       case(east_,west_)  
          BRefl_D(X_) = 2*BOrig_D(x_)
          BRefl_D(Y_) = 0.0
          BRefl_D(Z_) = 0.0         
       case(south_,north_)                                                 
          BRefl_D(X_) = 0.0
          BRefl_D(Y_) = 2*BOrig_D(y_)
          BRefl_D(Z_) = 0.0 
       case(bot_,top_)                                                     
          BRefl_D(X_) = 0.0
          Brefl_D(Y_) = 0.0 
          BRefl_D(Z_) = 2*BOrig_D(z_)
       end select
       VarsGhostFace_V           = VarsTrueFace_V
       VarsGhostFace_V(iRhoUx_I) = -VarsGhostFace_V(iRhoUx_I)
       VarsGhostFace_V(iRhoUy_I) = -VarsGhostFace_V(iRhoUy_I)
       VarsGhostFace_V(iRhoUz_I) = -VarsGhostFace_V(iRhoUz_I)
       VarsGhostFace_V(Bx_:Bz_)  = VarsGhostFace_V(Bx_:Bz_) - BRefl_D
    case('fixedB1')
       VarsGhostFace_V  = FaceState_V
    case('ionosphere')
       where(DefaultState_V(1:nVar-1)>cTiny)
          VarsGhostFace_V(1:nVar-1)    = VarsTrueFace_V(1:nVar-1)+&
            sign(cOne,FaceState_V(1:nVar-1) - VarsTrueFace_V(1:nVar-1))*&
            min(abs(FaceState_V(1:nVar-1)-VarsTrueFace_V(1:nVar-1)),&
            DensityJumpLimit*VarsTrueFace_V(1:nVar-1))
       elsewhere
          VarsGhostFace_V(1:nVar-1)=VarsTrueFace_V(1:nVar-1)
       end where
       VarsGhostFace_V(iP_I)      = VarsTrueFace_V(iP_I) + &
            sign(cOne,FaceState_V(iP_I) - VarsTrueFace_V(iP_I))*&
            min(abs(FaceState_V(iP_I) - VarsTrueFace_V(iP_I)),&
            PressureJumpLimit*VarsTrueFace_V(iP_I))

       VarsGhostFace_V(iRhoUx_I) = -VarsGhostFace_V(iRhoUx_I)
       VarsGhostFace_V(iRhoUy_I) = -VarsGhostFace_V(iRhoUy_I)
       VarsGhostFace_V(iRhoUz_I) = -VarsGhostFace_V(iRhoUz_I)
    case('coronatoih')    !Only for nVar=8

       !Get interpolated values from buffer grid:
       call get_from_spher_buffer_grid(&
            FaceCoords_D,nVar,FaceState_V)

       VarsGhostFace_V = FaceState_V
    case('polarwind')
       ! Get density/ies and velocity from polarwind code
       call read_pw_buffer(FaceCoords_D,nVar,FaceState_V)
       VarsGhostFace_V = FaceState_V

       ! Apply floating conditions on P and B
       VarsGhostFace_V(iP_I)    = VarsTrueFace_V(iP_I)
       VarsGhostFace_V(Bx_:Bz_) = VarsTrueFace_V(Bx_:Bz_)
    case default
       call stop_mpi('Incorrect TypeBc_I='//TypeBc)
    end select
!^CFG IF IONOSPHERE BEGIN
    if (UseIonosphere .and. iBoundary == Body1_) then
       call calc_inner_bc_velocity(TimeBc, FaceCoords_D, &
            VarsTrueFace_V(Bx_:Bz_), B0Face_D, uIono_D)
       
       select case(TypeBc)
       case('reflect','linetied','polarwind','ionosphere','ionospherefloat')
          VarsGhostFace_V(iRhoUx_I) = 2*uIono_D(x_) + VarsGhostFace_V(iRhoUx_I)
          VarsGhostFace_V(iRhoUy_I) = 2*uIono_D(y_) + VarsGhostFace_V(iRhoUy_I)
          VarsGhostFace_V(iRhoUz_I) = 2*uIono_D(z_) + VarsGhostFace_V(iRhoUz_I)
       case default
          call stop_mpi('UseIonosphere is not compatible with TypeBc_I=' &
               //TypeBc)
       end select
    end if
!^CFG END IONOSPHERE
    if (UseRotatingBc .and. iBoundary==Body1_) then

       !\
       ! The program is called which calculates the cartesian corotation 
       ! velocity vector uRot_D as a function of the radius-vector "FaceCoords"
       !/
       call calc_corotation_velocities(FaceCoords_D, uRot_D)

       select case(TypeBc)
       case('reflect','linetied', &
            'ionosphere','ionospherefloat','polarwind')
          VarsGhostFace_V(iRhoUx_I) = 2*uRot_D(x_) + VarsGhostFace_V(iRhoUx_I)
          VarsGhostFace_V(iRhoUy_I) = 2*uRot_D(y_) + VarsGhostFace_V(iRhoUy_I)
          VarsGhostFace_V(iRhoUz_I) = 2*uRot_D(z_) + VarsGhostFace_V(iRhoUz_I)
       case default
          call stop_mpi('UseRotatingBc is not compatible with TypeBc_I=' &
               //TypeBc) 
       end select
    end if 
  end subroutine set_face_bc

end subroutine set_face_BCs
