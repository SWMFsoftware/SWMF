!^CFG COPYRIGHT UM
subroutine set_BCs(TimeBcIn, DoResChangeOnlyIn)
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

subroutine set_face_BCs(IsBodyCell_G, IsTrueCell_G)

  use ModMain
  use ModGeometry, ONLY : x_BLK, y_BLK, z_BLK
  use ModAdvance
  use ModParallel, ONLY : neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth
  use ModNumConst
  use ModFaceBc
  use ModPhysics, ONLY: PolarRho_I, PolarU_I, PolarP_I, PolarTheta
  use ModSolarwind, ONLY: get_solar_wind_point
  use CON_axes, ONLY: transform_matrix

  implicit none

  logical, dimension(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn), intent(in) :: &
       IsBodyCell_G, IsTrueCell_G 

  integer :: i,j,k

  ! Variables used for polar wind boundary condition
  real :: GmToSmg_DD(3,3), CoordSm_D(3), Cos2PolarTheta, bUnit_D(3)
  logical :: IsPolarFace

  character (len=*), parameter :: NameSub = 'set_face_bcs'
  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------
  if(iBlockBc==BLKtest.and.iProc==PROCtest)then
     call set_oktest(NameSub, DoTest, DoTestMe)
  else
     DoTest = .false.; DoTestMe = .false.
  end if

  TypeBc = TypeBc_I(iBoundary)

  if(TypeBc == 'polarwind') then
     GmToSmg_DD = transform_matrix(Time_Simulation, TypeCoordSystem, 'SMG')
     Cos2PolarTheta = cos(PolarTheta)**2
  end if

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
              B0Face_D = B0_DX(:,i,j,k)
     
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
              B0Face_D = B0_DX(:,i,j,k)
            
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
              B0Face_D     = B0_DY(:,i,j,k)
              
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
              B0Face_D = B0_DY(:,i,j,k)

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
              B0Face_D = B0_DZ(:,i,j,k)

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
              B0Face_D = B0_DZ(:,i,j,k)
     
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
    use ModMultiFluid, ONLY: iUx_I, iUy_I, iUz_I
    
    implicit none

    real, parameter:: PressureJumpLimit=0.0, DensityJumpLimit=0.1

    real:: uRot_D(nDim), uIono_D(nDim)
    real:: FaceState_V(nVar)
    real:: bDotR, Brefl_D(nDim), Borig_D(nDim)
    real:: bDotU, rInv
    real:: CosTheta, SinTheta, CosPhi, SinPhi
    real:: UrTrue, UtTrue, BpTrue, BrGhost, BtGhost, BpGhost

    real:: bUnit_D(3)
    logical:: IsPolarFace
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

    ! Default fixed/initial state for this boundary
    FaceState_V = FaceState_VI(:, iBoundary)  

    select case(TypeBc) 
    case('linetied','ionospherefloat')
       VarsGhostFace_V        =  VarsTrueFace_V
       VarsGhostFace_V(iUx_I) = -VarsTrueFace_V(iUx_I)
       VarsGhostFace_V(iUy_I) = -VarsTrueFace_V(iUy_I)
       VarsGhostFace_V(iUz_I) = -VarsTrueFace_V(iUz_I)

    case('float','outflow')
       VarsGhostFace_V = VarsTrueFace_V

    case('heliofloat')
       VarsGhostFace_V = VarsTrueFace_V
       rInv     = 1.0/sqrt(sum(FaceCoords_D**2))
       CosTheta = FaceCoords_D(z_)*rInv
       SinTheta = sqrt(FaceCoords_D(x_)**2+FaceCoords_D(y_)**2)*RInv
       CosPhi   = FaceCoords_D(x_)/ &
            sqrt(FaceCoords_D(x_)**2+FaceCoords_D(y_)**2+cTolerance**2)
       SinPhi   = FaceCoords_D(y_)/ &
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

    case('fixedB1')
       VarsGhostFace_V = FaceState_V

    case('fixed')
       VarsGhostFace_V = FaceState_V
       VarsGhostFace_V(Bx_:Bz_) = VarsGhostFace_V(Bx_:Bz_) - B0Face_D

    case('inflow','vary')
       call get_solar_wind_point(TimeBc, FaceCoords_D(x_), VarsGhostFace_V)
       VarsGhostFace_V(Bx_:Bz_) = VarsGhostFace_V(Bx_:Bz_) - B0Face_D

    case('reflect','reflectb')
       ! reflect the full velocity vector and
       ! reflect the normal component of B1 (reflect) or full B (reflectb)

       Borig_D = VarsTrueFace_V(Bx_:Bz_)
       if(TypeBc == 'reflectb') Borig_D = Borig_D + B0Face_D

       select case(iBoundary)                                                 
       case(body1_, body2_)
          bDotR   = 2*dot_product(Borig_D, FaceCoords_D)/sum(FaceCoords_D**2)
          Brefl_D = FaceCoords_D*bDotR
       case(east_, west_)  
          Brefl_D = (/ 2*Borig_D(x_), 0.0, 0.0 /)
       case(south_, north_)                                                 
          Brefl_D = (/ 0.0, 2*Borig_D(y_), 0.0 /)
       case(bot_, top_)                                                     
          Brefl_D = (/ 0.0, 0.0, 2*Borig_D(z_) /)
       end select

       ! Apply floating condition on densities and pressures
       VarsGhostFace_V          =  VarsTrueFace_V

       ! Reflect all components of velocities
       VarsGhostFace_V(iUx_I)   = -VarsGhostFace_V(iUx_I)
       VarsGhostFace_V(iUy_I)   = -VarsGhostFace_V(iUy_I)
       VarsGhostFace_V(iUz_I)   = -VarsGhostFace_V(iUz_I)

       ! Reflect B1 or full B
       VarsGhostFace_V(Bx_:Bz_) =  VarsTrueFace_V(Bx_:Bz_) - BRefl_D

    case('ionosphere', 'polarwind')

       if(TypeBc == 'polarwind')then
          CoordSm_D = matmul(GmToSmg_DD, FaceCoords_D)
          IsPolarFace = CoordSm_D(z_)**2/sum(CoordSm_D**2) > Cos2PolarTheta
       else
          IsPolarFace = .false.
       end if

       if(IsPolarFace)then
          ! polarwind type conditions
          if(UsePw)then
             ! Get density/ies and velocity from polarwind code
             call read_pw_buffer(FaceCoords_D,nVar,FaceState_V)
             VarsGhostFace_V = FaceState_V

             ! Apply floating conditions on P and B
             VarsGhostFace_V(iP_I)    = VarsTrueFace_V(iP_I)
             VarsGhostFace_V(Bx_:Bz_) = VarsTrueFace_V(Bx_:Bz_)
          else
             ! Use variables set in the #POLARBOUNDARY command
             VarsGhostFace_V(iRho_I) = PolarRho_I

             ! Align flow with the magnetic field
             bUnit_D = B0Face_D / sqrt(sum(B0Face_D**2))
             ! Make sure it points outward
             if(sum(bUnit_D*FaceCoords_D) < 0.0) bUnit_D = -bUnit_D
             VarsGhostFace_V(iUx_I)  = PolarU_I*bUnit_D(x_)
             VarsGhostFace_V(iUy_I)  = PolarU_I*bUnit_D(y_)
             VarsGhostFace_V(iUz_I)  = PolarU_I*bUnit_D(z_)
             VarsGhostFace_V(iP_I)   = PolarP_I
          end if
       else
          ! Ionosphere type conditions
          where(DefaultState_V(1:nVar) > cTiny)
             ! Use body densities but limit jump
             VarsGhostFace_V = VarsTrueFace_V + &
                  sign(1.0, FaceState_V - VarsTrueFace_V)*   &
                  min( abs(FaceState_V - VarsTrueFace_V)     &
                  ,    DensityJumpLimit*VarsTrueFace_V   )
          elsewhere
             ! Apply floating
             VarsGhostFace_V = VarsTrueFace_V
          end where

          if(PressureJumpLimit > 0.0) then
             ! Use body pressures but limit jump
             VarsGhostFace_V(iP_I) = VarsTrueFace_V(iP_I) + &
                  sign(cOne,FaceState_V(iP_I) - VarsTrueFace_V(iP_I))*&
                  min(abs(FaceState_V(iP_I) - VarsTrueFace_V(iP_I)),&
                  PressureJumpLimit*VarsTrueFace_V(iP_I))
          else
             ! Use floating BC for pressure (correct for zero radial velocity)
             VarsGhostFace_V(iP_I) = VarsTrueFace_V(iP_I)
          end if

          ! Change sign for velocities (plasma frozen into dipole field)
          VarsGhostFace_V(iUx_I) = -VarsTrueFace_V(iUx_I)
          VarsGhostFace_V(iUy_I) = -VarsTrueFace_V(iUy_I)
          VarsGhostFace_V(iUz_I) = -VarsTrueFace_V(iUz_I)
       end if

    case('coronatoih')    !Only for nVar=8

       !Get interpolated values from buffer grid:
       call get_from_spher_buffer_grid(&
            FaceCoords_D,nVar,FaceState_V)

       VarsGhostFace_V = FaceState_V
    case default
       call stop_mpi('Incorrect TypeBc_I='//TypeBc)
    end select

    !^CFG IF IONOSPHERE BEGIN
    if (UseIe .and. iBoundary == Body1_) then
       ! Get the E x B / B^2 velocity
       call calc_inner_bc_velocity(TimeBc, FaceCoords_D, &
            VarsTrueFace_V(Bx_:Bz_) + B0Face_D, uIono_D)

       ! Subtract the radial component of the velocity (no outflow/inflow)
       uIono_D = uIono_D &
            - FaceCoords_D * sum(FaceCoords_D * uIono_D) / sum(FaceCoords_D**2)

       select case(TypeBc)
       case('reflect','linetied','polarwind','ionosphere','ionospherefloat')
          VarsGhostFace_V(iUx_I) = 2*uIono_D(x_) + VarsGhostFace_V(iUx_I)
          VarsGhostFace_V(iUy_I) = 2*uIono_D(y_) + VarsGhostFace_V(iUy_I)
          VarsGhostFace_V(iUz_I) = 2*uIono_D(z_) + VarsGhostFace_V(iUz_I)

       case default
          call stop_mpi('Coupling with IE is not compatible with TypeBc_I=' &
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
          VarsGhostFace_V(iUx_I) = 2*uRot_D(x_) + VarsGhostFace_V(iUx_I)
          VarsGhostFace_V(iUy_I) = 2*uRot_D(y_) + VarsGhostFace_V(iUy_I)
          VarsGhostFace_V(iUz_I) = 2*uRot_D(z_) + VarsGhostFace_V(iUz_I)
       case default
          call stop_mpi('UseRotatingBc is not compatible with TypeBc_I=' &
               //TypeBc) 
       end select
    end if 
  end subroutine set_face_bc

end subroutine set_face_BCs
