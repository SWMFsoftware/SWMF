!^CFG COPYRIGHT UM
subroutine set_BCs(iter,time_now,DoResChangeOnly)
  use ModProcMH
  use ModMain
  use ModAdvance
  use ModNumConst
  use ModGeometry, ONLY:&
       IsBoundaryCell_GI , IsBoundaryBlock_IB,true_cell,MinBoundary,MaxBoundary
  implicit none

  integer, intent(in) :: iter
  real, intent(in) :: time_now
  logical, intent (in) :: DoResChangeOnly
  integer:: iBoundary
  logical :: oktest, oktest_me
  !----------------------------------------------------------------------------

  call timing_start('set_BCs')

  if(globalBLK==BLKtest.and.iProc==PROCtest)then
     call set_oktest('set_BCs',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  endif

  if(oktest_me.and. .not.DoResChangeOnly)write(*,*)'set_BCs initial Ux,Uy=',&
       LeftState_VX(Ux_,Itest,Jtest,Ktest),LeftState_VX(Uy_,Itest,Jtest,Ktest)

  call set_boundary_cells(globalBLK)
  !\
  ! Apply boundary conditions
  !/
  do iBoundary=MinBoundary,MaxBoundary
     if(IsBoundaryBlock_IB(iBoundary,globalBLK))&
          call set_face_BCs(iter,time_now,DoResChangeOnly,&
          IsBoundaryCell_GI(:,:,:,iBoundary),&
          true_cell(:,:,:,globalBLK),iBoundary)
  end do
  

  if(oktest_me.and. .not.DoResChangeOnly)write(*,*)'set_BCs final Ux,Uy=',&
       LeftState_VX(Ux_,Itest,Jtest,Ktest),&
       LeftState_VX(Uy_,Itest,Jtest,Ktest)

  if(oktest_me)then
     write(*,*)'east  PfaceL_x,PfaceR_x=',&
          LeftState_VX(P_,Itest,Jtest,Ktest),  RightState_VX(P_,Itest,Jtest,Ktest)
     write(*,*)'west  PfaceL_x,PfaceR_x=',&
          LeftState_VX(P_,Itest+1,Jtest,Ktest),RightState_VX(P_,Itest+1,Jtest,Ktest)
     write(*,*)'south PfaceL_y,PfaceR_y=',&
          LeftState_VY(P_,Itest,Jtest,Ktest),  RightState_VY(P_,Itest,Jtest,Ktest)
     write(*,*)'north PfaceL_y,PfaceR_y=',&
          LeftState_VY(P_,Itest,Jtest+1,Ktest),RightState_VY(P_,Itest,Jtest+1,Ktest)
     write(*,*)'bot   PfaceL_z,PfaceR_z=',&
          LeftState_VZ(P_,Itest,Jtest,Ktest),  RightState_VZ(P_,Itest,Jtest,Ktest)
     write(*,*)'top   PfaceL_z,PfaceR_z=',&
          LeftState_VZ(P_,Itest,Jtest,Ktest+1),RightState_VZ(P_,Itest,Jtest,Ktest+1)
  end if

  call timing_stop('set_BCs')

end subroutine set_BCs

subroutine set_face_BCs(iter,time_now,DoResChangeOnly,&
     IsBodyCell,IsTrueCell,iBoundary)

  use ModMain
  use ModGeometry, ONLY : x_BLK, y_BLK, z_BLK
  use ModAdvance
  use ModParallel, ONLY : neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth
  use ModNumConst
  implicit none

  integer, intent(in) :: iter, iBoundary
  real, intent(in) :: time_now
  logical, intent (in) :: DoResChangeOnly
  logical, dimension(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn), intent(in) :: &
       IsBodyCell, IsTrueCell 
  logical :: UseIonosphereHere, UseCorotationHere
  character (len=20) :: TypeBcHere
  
  real,dimension(nFaceValueVars)::VarsTrueFace_V,VarsGhostFace_V


  ! Negative iBoundary indicates which body we are computing for.
  ! Positive iBoundary numerates the sides for the outer BCs.
  ! Here iter and time now are the iteration number and physical time,
  ! which may be used for time-dependent BCs (not used now).
  ! 
  ! RHere is 3d array of the distance values apart from the inner boundary center,
  ! xBodyHere, yBodyHere,zBodyHere are the physical coordinates of the inner boundary center 
  ! UseIonesphereHere=.true. allows to call subroutine calc_innerBC_velocities
  ! UseCorotationHere = .true. allows to call calc_corotation_velocities

  integer :: i,j,k

  real,dimension(1:nDim):: FaceCoords_D, B0Face_D

  if(iBoundary==Body1_)then
     UseIonosphereHere=UseIonosphere
     UseCorotationHere=UseCorotation
  else
     UseIonosphereHere=.false.
     UseCorotationHere=.false.
  end if

  TypeBcHere=TypeBc_I(iBoundary)
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
           if (IsTrueCell(i-1,j,k) .and. &
                IsBodyCell(i,j,k) .and. &
                (.not.DoResChangeOnly .or. &
                ((i == nIFace .and. neiLwest(globalBLK)==+1) .or. &
                (i == 1       .and. neiLeast(globalBLK)==+1)) )) then

              FaceCoords_D(x_) = 0.50*(x_BLK(i-1,j,k,globalBLK)+x_BLK(i,j,k,globalBLK))
              FaceCoords_D(y_) = 0.50*(y_BLK(i-1,j,k,globalBLK)+y_BLK(i,j,k,globalBLK))
              FaceCoords_D(z_) = 0.50*(z_BLK(i-1,j,k,globalBLK)+z_BLK(i,j,k,globalBLK))
              B0Face_D(x_) = B0xFace_x_BLK(i,j,k,globalBLK)
              B0Face_D(y_) = B0yFace_x_BLK(i,j,k,globalBLK)
              B0Face_D(z_) = B0zFace_x_BLK(i,j,k,globalBLK)

              VarsTrueFace_V= LeftState_VX(:,i,j,k)

              call set_body_BCs(West_)             

              RightState_VX(:,i,j,k) = VarsGhostFace_V

           end if

           if (IsTrueCell(i,j,k) .and. &
                IsBodyCell(i-1,j,k)  .and. &
                (.not.DoResChangeOnly .or. &
                (i == 1         .and. neiLeast(globalBLK)==+1) .or. &
                (i == nIFace    .and. neiLwest(globalBLK)==+1)  )) then

              FaceCoords_D(x_) = 0.50*(x_BLK(i,j,k,globalBLK)+x_BLK(i-1,j,k,globalBLK))
              FaceCoords_D(y_) = 0.50*(y_BLK(i,j,k,globalBLK)+y_BLK(i-1,j,k,globalBLK))
              FaceCoords_D(z_) = 0.50*(z_BLK(i,j,k,globalBLK)+z_BLK(i-1,j,k,globalBLK))
              B0Face_D(x_)     = B0xFace_x_BLK(i,j,k,globalBLK)
              B0Face_D(y_)     = B0yFace_x_BLK(i,j,k,globalBLK)
              B0Face_D(z_)     = B0zFace_x_BLK(i,j,k,globalBLK)
            
              VarsTrueFace_V = RightState_VX(:,i,j,k)

              call set_body_BCs(East_)

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
           if (IsTrueCell(i,j-1,k) .and. &
                IsBodyCell(i,j,k)  .and. &
                ( .not.DoResChangeOnly .or. &
                (j == nJFace .and. neiLnorth(globalBLK)==+1) .or. &
                (j == 1      .and. neiLsouth(globalBLK)==+1) )) then
              FaceCoords_D(x_) = 0.50*(x_BLK(i,j-1,k,globalBLK)+x_BLK(i,j,k,globalBLK))
              FaceCoords_D(y_) = 0.50*(y_BLK(i,j-1,k,globalBLK)+y_BLK(i,j,k,globalBLK))
              FaceCoords_D(z_) = 0.50*(z_BLK(i,j-1,k,globalBLK)+z_BLK(i,j,k,globalBLK))
              B0Face_D(x_)     = B0xFace_y_BLK(i,j,k,globalBLK)
              B0Face_D(y_)     = B0yFace_y_BLK(i,j,k,globalBLK)
              B0Face_D(z_)     = B0zFace_y_BLK(i,j,k,globalBLK)

              VarsTrueFace_V = LeftState_VY(:,i,j,k)
          
              call set_body_BCs(North_)

              RightState_VY(:,i,j,k) = VarsGhostFace_V           
           end if

           if (IsTrueCell(i,j,k) .and. &
                IsBodyCell(i,j-1,k)  .and. &
                (.not.DoResChangeOnly .or. &
                (j ==1         .and. neiLsouth(globalBLK)==+1) .or. &
                (j == nJFace .and. neiLnorth(globalBLK)==+1) )) then
              FaceCoords_D(x_) = 0.50*(x_BLK(i,j-1,k,globalBLK)+x_BLK(i,j,k,globalBLK))
              FaceCoords_D(y_) = 0.50*(y_BLK(i,j-1,k,globalBLK)+y_BLK(i,j,k,globalBLK))
              FaceCoords_D(z_) = 0.50*(z_BLK(i,j-1,k,globalBLK)+z_BLK(i,j,k,globalBLK))
              B0Face_D(x_)     = B0xFace_y_BLK(i,j,k,globalBLK)
              B0Face_D(y_)     = B0yFace_y_BLK(i,j,k,globalBLK)
              B0Face_D(z_)     = B0zFace_y_BLK(i,j,k,globalBLK)

              VarsTrueFace_V = RightState_VY(:,i,j,k)

              call set_body_BCs(South_)

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
           if (IsTrueCell(i,j,k-1) .and. &
                IsBodyCell(i,j,k) .and. &
                (.not.DoResChangeOnly .or. &
                (k == nKFace .and. neiLtop(globalBLK)==+1) .or. &
                (k == 1       .and. neiLbot(globalBLK)==+1)) ) then
              FaceCoords_D(x_)= 0.50*(x_BLK(i,j,k,globalBLK)+x_BLK(i,j,k-1,globalBLK))
              FaceCoords_D(y_)= 0.50*(y_BLK(i,j,k,globalBLK)+y_BLK(i,j,k-1,globalBLK))
              FaceCoords_D(z_)= 0.50*(z_BLK(i,j,k,globalBLK)+z_BLK(i,j,k-1,globalBLK))
              B0Face_D(x_)    = B0xFace_z_BLK(i,j,k,globalBLK)
              B0Face_D(y_)    = B0yFace_z_BLK(i,j,k,globalBLK)
              B0Face_D(z_)    = B0zFace_z_BLK(i,j,k,globalBLK)

              VarsTrueFace_V =  LeftState_VZ(:,i,j,k)
        
              call set_body_BCs(Top_)

              RightState_VZ(:,i,j,k) = VarsGhostFace_V
           end if

           if (IsTrueCell(i,j,k).and. &
                IsBodyCell(i,j,k-1).and. &
                (.not.DoResChangeOnly .or. &
                (k == 1         .and. neiLbot(globalBLK)==+1) .or. &
                (k == nKFace .and. neiLtop(globalBLK)==+1))  ) then
              FaceCoords_D(x_) = 0.50*(x_BLK(i,j,k-1,globalBLK)+x_BLK(i,j,k,globalBLK))
              FaceCoords_D(y_) = 0.50*(y_BLK(i,j,k-1,globalBLK)+y_BLK(i,j,k,globalBLK))
              FaceCoords_D(z_) = 0.50*(z_BLK(i,j,k-1,globalBLK)+z_BLK(i,j,k,globalBLK))
              B0Face_D(x_)     = B0xFace_z_BLK(i,j,k,globalBLK)
              B0Face_D(y_)     = B0yFace_z_BLK(i,j,k,globalBLK)
              B0Face_D(z_)     = B0zFace_z_BLK(i,j,k,globalBLK)

              VarsTrueFace_V =  RightState_VZ(:,i,j,k)
         
              call set_body_BCs(Bot_)

              LeftState_VZ(:,i,j,k) = VarsGhostFace_V
          
           end if
        end do !end i loop
     end do !end j loop
  end do !end k loop
 
contains

  subroutine set_body_BCs(iSide)
    use ModPhysics, ONLY : xBody2,yBody2,zBody2 !^CFG IF SECONDBODY
    use ModPhysics, ONLY : &
         FaceState_VI,UnitSI_rho,UnitSI_U,UnitSI_B,UnitSI_p
    implicit none
    integer,intent(in)::iSide !is defined with respect to the TRUE CELL!!!
    real, dimension(1:3) :: v_phi, FaceIono_D

    real:: FaceState_V(nFaceValueVars)
    real:: PressureJumpLimit=0.0,DensityJumpLimit=0.1    !
    ! Theoretical upper limit for DensityJumpLimit is (1-cfl) for nOrder=nStage=1
    ! Practical estimate 0.1 (may be (1-cfl)/2 ???) for nOrder=nStage=2
    !
    ! ___________ATTENTION!!!!!___________
    ! A non-zero value for the DensityJumpLimit physically corresponds 
    ! to a partially penetrable surface through which matter can diffuse. The 
    ! diffusion rate is restricted in a non-linear fashion due to the use 
    ! of a small DensityJumpLimit.  The "diffusion mass velocity" is 
    ! automatically lower than Alfven velocity*DensityJumpLimit.  Even the use 
    ! of small DensityJumpLimit STRONGLY INFLUENCES the numerical solution 
    ! in cases when RhoFace and RhoBodyHere are strongly different  near some 
    ! part of the boundary.  Even when using a non-zero DensityJumpLimit, the
    ! outer density tends to RhoBodyHere after a huge number of iterations.
    !
    !___________ATTENTION!!!!!!___________
    ! A non-zero value for PressureJumpLimit is not recommended. A theoretical 
    ! upper limit is unknown and numerical experiments with a non-zero value are
    ! not robust. NO REAL PHYSICAL process corresponds to a non-zero 
    ! PressureJumpLimit which would be like "pressure diffusion".

    real ::BdotR,BRefl_D(nDim)
    real:: BdotU,RInv
    real:: cosTheta,sinTheta,cosPhi,sinPhi
    real:: UrTrue,UtTrue,BpTrue,BrGhost,BtGhost,BpGhost
!^CFG IF USERFILES BEGIN
!    For new ionosphere, multispecies, multifluids 
    if(index(TypeBcHere,'user')>0.or.UseUserInnerBCs)then
       call user_face_bcs(i,j,k,globalBLK,iSide,iBoundary,&
            iter,time_now, FaceCoords_D,&
            VarsTrueFace_V,VarsGhostFace_V,&
            B0Face_D,  UseIonosphereHere,UseCorotationHere)
       return
    end if
!^CFG END USERFILES

    !^CFG IF SECONDBODY BEGIN
    if(iBoundary==body2_)then
       FaceCoords_D(x_)= FaceCoords_D(x_)-xBody2
       FaceCoords_D(y_)= FaceCoords_D(y_)-yBody2
       FaceCoords_D(z_)= FaceCoords_D(z_)-zBody2
    end if
    !^CFG END SECONDBODY                                                                
    FaceState_V=FaceState_VI(:,iBoundary)  
    if(iBoundary==west_.and.TypeBcHere=='vary') then                          
       call get_solar_wind_point(&
            time_now,         &
            FaceCoords_D(y_), &
            FaceCoords_D(z_), &
            FaceState_V(rho_),&
            FaceState_V(Ux_ ),&
            FaceState_V(Uy_ ),&
            FaceState_V(Uz_ ),&
            FaceState_V(Bx_ ),&
            FaceState_V(By_ ),&
            FaceState_V(Bz_ ),&
            FaceState_V(P_))            
    end if                                                        

    select case(TypeBcHere) 
    case('linetied','ionospherefloat')
       VarsGhostFace_V=VarsTrueFace_V
       VarsGhostFace_V(Ux_:Uz_)=-VarsGhostFace_V(Ux_:Uz_)
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
                  dot_product(VarsTrueFace_V(Ux_:Uz_),        &
                              VarsTrueFace_V(Ux_:Uz_))
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
    case('reflect')   
       select case(iBoundary)                                                 
       case(body1_,body2_)                                                        
          BdotR  = dot_product(VarsTrueFace_V(Bx_:Bz_),FaceCoords_D)*&
               cTwo/dot_product(FaceCoords_D,FaceCoords_D)
          BRefl_D = FaceCoords_D*BdotR
       case(east_,west_)  
          BRefl_D(X_) = cTwo*VarsTrueFace_V(Bx_)
          BRefl_D(Y_) = cZero
          BRefl_D(Z_) = cZero         
       case(south_,north_)                                                 
          BRefl_D(X_) = cZero
          BRefl_D(Y_) = cTwo*VarsTrueFace_V(By_)
          BRefl_D(Z_) = cZero         
       case(bot_,top_)                                                     
          BRefl_D(X_) = cZero;
          Brefl_D(Y_) = cZero; 
          BRefl_D(Z_) = cTwo*VarsTrueFace_V(Bz_)         
       end select
       VarsGhostFace_V=VarsTrueFace_V
       VarsGhostFace_V(Ux_:Uz_)     = - VarsGhostFace_V(Ux_:Uz_)
       VarsGhostFace_V(Bx_:Bz_)     = VarsGhostFace_V(Bx_:Bz_) - BRefl_D
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
       VarsGhostFace_V(P_)      = VarsTrueFace_V(P_) + &
            sign(cOne,FaceState_V(P_) - VarsTrueFace_V(P_))*&
            min(abs(FaceState_V(P_) - VarsTrueFace_V(P_)),&
            PressureJumpLimit*VarsTrueFace_V(P_))

       VarsGhostFace_V(Ux_:Uz_)     = - VarsGhostFace_V(Ux_:Uz_)
    case('coronatoih')    !Only for nVar=8
       !Get interpolated values from buffer grid:
       call get_from_spher_buffer_grid(&
            FaceCoords_D,nVar,FaceState_V)
       !Transform to primitive variables
       FaceState_V(Ux_:Uz_)=&
            FaceState_V(rhoUx_:rhoUz_)/FaceState_V(rho_)
       !Convert from SI:
       FaceState_V(rho_)      =&
             FaceState_V(rho_)    /UnitSI_rho
       FaceState_V(Ux_:Uz_)   =&
             FaceState_V(Ux_:Uz_) /UnitSI_U
       FaceState_V(Bx_:Bz_)   =&
            FaceState_V(Bx_:Bz_)  /UnitSI_B
       FaceState_V(P_)        =&
             FaceState_V(P_)      /UnitSI_P
       VarsGhostFace_V = cTwo*FaceState_V - VarsTrueFace_V
    case default
       call stop_mpi('Incorrect TypeBc_I='//TypeBcHere)
    end select
!^CFG IF IONOSPHERE BEGIN
    if (UseIonosphereHere) then
       call calc_inner_BC_velocities(iter,time_now,&
            FaceCoords_D(x_),FaceCoords_D(y_),FaceCoords_D(z_), &
            VarsTrueFace_V(Bx_),VarsTrueFace_V(By_),VarsTrueFace_V(Bz_),&
            B0Face_D(x_),B0Face_D(y_),B0Face_D(z_),&
            FaceIono_D(x_),FaceIono_D(y_),FaceIono_D(z_))
       select case(TypeBcHere)
       case('reflect','linetied',&
            'ionosphere','ionospherefloat')
          VarsGhostFace_V(Ux_:Uz_) = cTwo*FaceIono_D + VarsGhostFace_V(Ux_:Uz_)
       case default
          call stop_mpi('UseIonosphere is not compatible with TypeBc_I='//TypeBcHere)
       end select
    end if
!^CFG END IONOSPHERE
    if (UseCorotationHere) then

       !\
       ! The program is called which calculates the cartesian corotation velocity vector
       ! v_phi as a function of the radius-vector "FaceCoords"
       !/
       call calc_corotation_velocities(iter,time_now,FaceCoords_D,v_phi)

       select case(TypeBcHere)
       case('reflect','linetied', &
            'ionosphere','ionospherefloat')
          VarsGhostFace_V(Ux_:Uz_) = VarsGhostFace_V(Ux_:Uz_) + cTwo*v_phi
       case default
          call stop_mpi('UseCorotation is not compatible with TypeBc_I='//TypeBcHere) 
       end select
    end if 
  end subroutine set_body_BCs

end subroutine set_face_BCs
