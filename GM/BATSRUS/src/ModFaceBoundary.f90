!^CFG COPYRIGHT UM
!==============================================================================
module ModFaceBoundary

  use ModVarIndexes, ONLY: nVar

  implicit none

  SAVE

  private ! except

  ! Public methods
  public :: set_face_boundary

  ! True if only boundaries at resolution changes are updated
  logical, public :: DoResChangeOnly

  ! The type and index of the boundary
  character(len=20), public :: TypeBc

  ! Negative iBoundary indicates which body we are computing for.
  ! Zero corresponds to the user defined extra boundary.
  ! Positive iBoundary indexes the outer boundaries from 1 to 2*nDim.
  integer, public :: iBoundary

  ! Index of the face
  integer, public :: iFace, jFace, kFace, iBlockBc

  ! The side of the cell defined with respect to the cell inside the domain
  integer, public :: iSide

  ! The values on the physical side and the ghost cell side of the boudary
  real, public :: VarsTrueFace_V(nVar), VarsGhostFace_V(nVar)

  ! The coordinates of the face center and the B0 field at that point
  real, public :: FaceCoords_D(3), B0Face_D(3)

  ! The time at which the (time dependent) boundary condition is calculated
  real, public :: TimeBc

contains

  !============================================================================
  subroutine set_face_boundary(iBlock, TimeBcIn, DoResChangeOnlyIn)

    use ModProcMH, ONLY: iProc
    use ModSize, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModMain, ONLY: ProcTest, BlkTest, iTest, jTest, kTest, VarTest
    use ModAdvance, ONLY: LeftState_VX, LeftState_VY, LeftState_VZ, &
         RightState_VX, RightState_VY, RightState_VZ
    use ModGeometry, ONLY: true_cell, MinBoundary, MaxBoundary
    use ModBoundaryCells, ONLY: iBoundary_GB

    integer, intent(in) :: iBlock
    real,    intent(in) :: TimeBcIn
    logical, intent(in) :: DoResChangeOnlyIn

    logical :: oktest, oktest_me
    character(len=*), parameter:: NameSub='set_face_boundary'
    logical, allocatable :: IsBodyCell_G(:,:,:)
    !--------------------------------------------------------------------------

    if(.not.allocated(IsBodyCell_G))&
         allocate(IsBodyCell_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK))

    call timing_start(NameSub)

    if(iBlock==BLKtest.and.iProc==PROCtest)then
       call set_oktest(NameSub, oktest,oktest_me)
    else
       oktest=.false.; oktest_me=.false.
    endif

    ! set variables in module
    TimeBc          = TimeBcIn
    DoResChangeOnly = DoResChangeOnlyIn
    iBlockBc        = iBlock

    if(oktest_me)call write_face_state('Initial')

    call set_boundary_cells(iBlockBc)

    IsBodyCell_G(:,:,:) = &
         iBoundary_GB(:,:,:,iBlockBc) >= MinBoundary .and. &
         iBoundary_GB(:,:,:,iBlockBc) <= MaxBoundary

    call set_face_bc(IsBodyCell_G, true_cell(:,:,:,iBlockBc) )

    if(oktest_me)call write_face_state('Final')

    call timing_stop(NameSub)

  contains

    subroutine write_face_state(String)

      character(len=*), intent(in):: String

      write(*,*) NameSub,' ',String,' face states:'
      write(*,*) 'VarL_x, VarR_x(iTest)  =',&
           LeftState_VX(VarTest,  Itest, Jtest, Ktest),  &
           RightState_VX(VarTest, Itest, Jtest, Ktest)
      write(*,*) 'VarL_x, VarR_x(iTest+1)=',&
           LeftState_VX(VarTest,  Itest+1, Jtest, Ktest), &
           RightState_VX(VarTest, Itest+1, Jtest, Ktest)
      write(*,*) 'VarL_y, VarR_y(jTest)  =',&
           LeftState_VY(VarTest,  Itest, Jtest, Ktest),  &
           RightState_VY(VarTest, Itest, Jtest, Ktest)
      write(*,*) 'VarL_y, VarR_y(jTest+1)=',&
           LeftState_VY(VarTest,  Itest, Jtest+1, Ktest), &
           RightState_VY(VarTest, Itest, Jtest+1, Ktest)
      write(*,*) 'VarL_z, VarR_z(kTest)  =',&
           LeftState_VZ(VarTest,  Itest, Jtest, Ktest), &
           RightState_VZ(VarTest, Itest, Jtest, Ktest)
      write(*,*) 'VarL_z, VarR_z(kTest+1)=',&
           LeftState_VZ(VarTest,  Itest, Jtest, Ktest+1), &
           RightState_VZ(VarTest, Itest, Jtest, Ktest+1)

    end subroutine write_face_state

  end subroutine set_face_boundary

  !============================================================================

  subroutine set_face_bc(IsBodyCell_G, IsTrueCell_G)

    use ModMain
    use ModProcMH, ONLY: iProc
    use ModB0, ONLY: B0_DX, B0_DY, B0_DZ
    use ModAdvance, ONLY: UseAnisoPressure, &
         LeftState_VX, LeftState_VY, LeftState_VZ, &
         RightState_VX, RightState_VY, RightState_VZ
    use ModParallel, ONLY: &
         neiLtop, neiLbot, neiLeast, neiLwest, neiLnorth, neiLsouth
    use ModNumConst
    use ModPhysics, ONLY: PolarRho_I, PolarU_I, PolarP_I, PolarTheta, &
         UseCpcpBc, Rho0Cpcp, RhoPerCpcp, &
         Io2No_V, No2Si_V, UnitRho_, UnitElectric_, UnitX_
    use ModSolarwind, ONLY: get_solar_wind_point
    use CON_axes, ONLY: transform_matrix
    use BATL_lib, ONLY: Xyz_DGB

    logical, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK), intent(in):: &
         IsBodyCell_G, IsTrueCell_G

    integer :: i,j,k

    ! Variables used for polar wind boundary condition
    real :: GmToSmg_DD(3,3), CoordSm_D(3), Cos2PolarTheta

    ! External function for ionosphere    !^CFG IF IONOSPHERE
    real, external :: logvar_ionosphere   !^CFG IF IONOSPHERE
    real:: RhoCpcp

    character (len=*), parameter :: NameSub = 'set_face_bc'
    logical :: DoTest, DoTestMe
    !--------------------------------------------------------------------------
    if(iBlockBc==BLKtest.and.iProc==PROCtest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if

    if(TypeBc_I(body1_) == 'polarwind') then
       GmToSmg_DD = transform_matrix(Time_Simulation, TypeCoordSystem, 'SMG')
       Cos2PolarTheta = cos(PolarTheta)**2
    end if

    !^CFG IF IONOSPHERE BEGIN
    ! Calculate inner BC density from cross polar cap potential if required
    ! Use KeV units for Cpcp and amu/cc for density.
    if(UseCpcpBc .and. UseIe) &
         RhoCpcp = Io2No_V(UnitRho_)*(Rho0Cpcp + RhoPerCpcp &
         * 0.5*(logvar_ionosphere('cpcpn') + logvar_ionosphere('cpcps')) &
         * (No2Si_V(UnitElectric_)*No2Si_V(UnitX_))/1000.0)
    !^CFG END IONOSPHERE

    !\
    ! Apply body boundary conditions as required.
    !/                            

    B0Face_D = 0.0

    do k = kMinFace, kMaxFace
       do j = jMinFace, jMaxFace
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

                iSide = 2

                FaceCoords_D = 0.5*(Xyz_DGB(:,i-1,j,k,iBlockBc) &
                     +              Xyz_DGB(:,i  ,j,k,iBlockBc))

                if(UseB0)B0Face_D = B0_DX(:,i,j,k)

                VarsTrueFace_V= LeftState_VX(:,i,j,k)

                call set_face(i-1,j,k,i,j,k)

                RightState_VX(:,i,j,k) = VarsGhostFace_V

             end if

             if (IsTrueCell_G(i,j,k) .and. &
                  IsBodyCell_G(i-1,j,k)  .and. &
                  (.not.DoResChangeOnly .or. &
                  (i == 1         .and. neiLeast(iBlockBc)==+1) .or. &
                  (i == nIFace    .and. neiLwest(iBlockBc)==+1)  )) then

                iSide = 1

                FaceCoords_D = 0.5*(Xyz_DGB(:,i-1,j,k,iBlockBc) &
                     +              Xyz_DGB(:,i  ,j,k,iBlockBc))

                if(UseB0)B0Face_D = B0_DX(:,i,j,k)

                VarsTrueFace_V = RightState_VX(:,i,j,k)

                call set_face(i,j,k,i-1,j,k)

                LeftState_VX(:,i,j,k) = VarsGhostFace_V
             end if
          end do !end i loop
       end do !end j loop
    end do !end k loop

    if(nDim == 1) RETURN
    do k = kMinFace, kMaxFace
       do j = 1 , nJFace
          do i = iMinFace, iMaxFace
             !\
             ! Apply BCs at Y-direction faces as necessary.
             !/
             if (IsTrueCell_G(i,j-1,k) .and. &
                  IsBodyCell_G(i,j,k)  .and. &
                  ( .not.DoResChangeOnly .or. &
                  (j == nJFace .and. neiLnorth(iBlockBc)==+1) .or. &
                  (j == 1      .and. neiLsouth(iBlockBc)==+1) )) then

                iSide = 4

                FaceCoords_D = 0.5*(Xyz_DGB(:,i,j-1,k,iBlockBc) &
                     +              Xyz_DGB(:,i,j  ,k,iBlockBc))

                if(UseB0)B0Face_D = B0_DY(:,i,j,k)

                VarsTrueFace_V = LeftState_VY(:,i,j,k)

                call set_face(i,j-1,k,i,j,k)

                RightState_VY(:,i,j,k) = VarsGhostFace_V           
             end if

             if (IsTrueCell_G(i,j,k) .and. &
                  IsBodyCell_G(i,j-1,k)  .and. &
                  (.not.DoResChangeOnly .or. &
                  (j ==1       .and. neiLsouth(iBlockBc)==+1) .or. &
                  (j == nJFace .and. neiLnorth(iBlockBc)==+1) )) then

                iSide = 3

                FaceCoords_D = 0.5*(Xyz_DGB(:,i,j-1,k,iBlockBc) &
                     +              Xyz_DGB(:,i,j  ,k,iBlockBc))

                if(UseB0)B0Face_D = B0_DY(:,i,j,k)

                VarsTrueFace_V = RightState_VY(:,i,j,k)

                call set_face(i,j,k,i,j-1,k)

                LeftState_VY(:,i,j,k) = VarsGhostFace_V
             end if
          end do !end j loop
       end do !end i loop
    end do !end k loop

    if(nDim == 2) RETURN
    do k = 1, nKFace
       do j = jMinFace, jMaxFace
          do i = iMinFace, iMaxFace
             !\
             ! Apply BCs at Z-direction faces as necessary.
             !/
             if (IsTrueCell_G(i,j,k-1) .and. &
                  IsBodyCell_G(i,j,k) .and. &
                  (.not.DoResChangeOnly .or. &
                  (k == nKFace .and. neiLtop(iBlockBc)==+1) .or. &
                  (k == 1       .and. neiLbot(iBlockBc)==+1)) ) then

                iSide = 6

                FaceCoords_D = 0.5*(Xyz_DGB(:,i,j,k-1,iBlockBc) &
                     +              Xyz_DGB(:,i,j,k  ,iBlockBc))

                if(UseB0)B0Face_D = B0_DZ(:,i,j,k)

                VarsTrueFace_V =  LeftState_VZ(:,i,j,k)

                call set_face(i,j,k-1,i,j,k)

                RightState_VZ(:,i,j,k) = VarsGhostFace_V
             end if

             if (IsTrueCell_G(i,j,k).and. &
                  IsBodyCell_G(i,j,k-1).and. &
                  (.not.DoResChangeOnly .or. &
                  (k == 1         .and. neiLbot(iBlockBc)==+1) .or. &
                  (k == nKFace .and. neiLtop(iBlockBc)==+1))  ) then

                iSide = 5

                FaceCoords_D = 0.5*(Xyz_DGB(:,i,j,k-1,iBlockBc) &
                     +              Xyz_DGB(:,i,j,k  ,iBlockBc))

                if(UseB0)B0Face_D = B0_DZ(:,i,j,k)

                VarsTrueFace_V =  RightState_VZ(:,i,j,k)

                call set_face(i,j,k,i,j,k-1)

                LeftState_VZ(:,i,j,k) = VarsGhostFace_V

             end if
          end do !end i loop
       end do !end j loop
    end do !end k loop

  contains

    subroutine set_face(iTrue, jTrue, kTrue, iGhost, jGhost, kGhost)

      use ModPhysics, ONLY : xBody2,yBody2,zBody2 !^CFG IF SECONDBODY
      use ModPhysics, ONLY : FaceState_VI,Si2No_V,No2Si_V,UnitX_,UnitN_, &
           UnitU_, UnitTemperature_, UnitJ_, UnitPoynting_,OrbitPeriod, &
           UseOutflowPressure, pOutflow
      use ModMain
      use ModMultiFluid
      use CON_planet_field, ONLY: get_planet_field, map_planet_field
      use ModConst,   ONLY: cElectronCharge, cBoltzmann,cProtonMass
      use ModPlanetConst, ONLY: Earth_, rPlanet_I
      use ModUtilities
      use ModBoundaryCells, ONLY: iBoundary_GB
      use ModNumConst, ONLY: cTwoPi

      ! indexes of the true and ghost cells on the two sides of the face
      integer, intent(in):: iTrue, jTrue, kTrue, iGhost, jGhost, kGhost

      real, parameter:: PressureJumpLimit=0.0, DensityJumpLimit=0.1
      real, parameter:: LatitudeCap = 55.0

      real:: uRot_D(MaxDim), uIono_D(MaxDim)
      real:: FaceState_V(nVar), State_V(Bx_:nVar+3)
      real:: bDotR, Brefl_D(MaxDim), Borig_D(MaxDim)
      real:: bDotU, rInv
      real:: CosTheta, SinTheta, CosPhi, SinPhi
      real:: UrTrue, UtTrue, BpTrue, BrGhost, BtGhost, BpGhost
      real:: Ub_V(2), b, b1,b4,bFace_D(3), JouleHeating, FluxIono, FluxPw
      real:: bUnit_D(3), GseToGeo_D(3), XyzMap_D(3), SmgFaceCoords_D(3), &
           GeoFaceCoords_D(3)
      logical:: IsPolarFace
      real:: SinLatitudeCap, zCap, eCap, ePar, &
           TheTmp,DtTmp,DaTmp, Cosx, Jlocal_D(3), Jpar
      integer:: iHemisphere
      integer :: iIonSecond
      !------------------------------------------------------------------------

      iBoundary = iBoundary_GB(iGhost,jGhost,kGhost,iBlockBc)
      TypeBc = TypeBc_I(iBoundary)

      ! User defined boundary conditions
      if( index(TypeBc, 'user') > 0 .or. &
           (UseUserInnerBCs .and. iBoundary <= body1_) .or. &
           (UseUserOuterBCs .and. iBoundary >= 1 ) )then
         iFace = i; jFace = j; kFace = k
         call user_set_face_boundary(VarsGhostFace_V)
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

      case('float')
         VarsGhostFace_V = VarsTrueFace_V

      case('outflow')
         VarsGhostFace_V = VarsTrueFace_V
         if(UseOutflowPressure) VarsGhostFace_V(p_) = pOutflow

      case('heliofloat')
         VarsGhostFace_V = VarsTrueFace_V
         rInv     = 1.0/sqrt(sum(FaceCoords_D**2))
         CosTheta = FaceCoords_D(z_)*rInv
         SinTheta = sqrt(FaceCoords_D(x_)**2+FaceCoords_D(y_)**2)*RInv
         CosPhi   = FaceCoords_D(x_)/ &
              sqrt(FaceCoords_D(x_)**2+FaceCoords_D(y_)**2+cTolerance**2)
         SinPhi   = FaceCoords_D(y_)/ &
              sqrt(FaceCoords_D(x_)**2+FaceCoords_D(y_)**2+cTolerance**2)
         BdotU    = sum(VarsTrueFace_V(Bx_:Bz_)*VarsTrueFace_V(Ux_:Uz_)) &
              /(sum(VarsTrueFace_V(Ux_:Uz_)*VarsTrueFace_V(Ux_:Uz_)) &
              + cTolerance**2)
         UrTrue   = sum(VarsTrueFace_V(Ux_:Uz_)*FaceCoords_D(x_:z_))*RInv
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
         call get_solar_wind_point(TimeBc, FaceCoords_D, VarsGhostFace_V)
         VarsGhostFace_V(Bx_:Bz_) = VarsGhostFace_V(Bx_:Bz_) - B0Face_D

      case('reflect','reflectb')
         ! reflect the full velocity vector and
         ! reflect the normal component of B1 (reflect) or full B (reflectb)

         Borig_D = VarsTrueFace_V(Bx_:Bz_)
         if(TypeBc == 'reflectb') Borig_D = Borig_D + B0Face_D

         select case(iBoundary)
         case(body1_, body2_)
            bDotR   = 2*sum(Borig_D*FaceCoords_D)/sum(FaceCoords_D**2)
            Brefl_D = FaceCoords_D*bDotR
         case(1, 2)  
            Brefl_D = (/ 2*Borig_D(x_), 0.0, 0.0 /)
         case(3, 4)                                                 
            Brefl_D = (/ 0.0, 2*Borig_D(y_), 0.0 /)
         case(5, 6)                                                     
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

      case('ionosphere', 'polarwind','ionosphereoutflow')

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

            ! Apply CPCP dependent density if required      !^CFG IF IONOSPHERE
            if(UseCpcpBc .and. UseIe) &                     !^CFG IF IONOSPHERE
                 VarsGhostFace_V(Rho_) = RhoCpcp            !^CFG IF IONOSPHERE

            if(PressureJumpLimit > 0.0) then
               ! Use body pressures but limit jump
               VarsGhostFace_V(iP_I) = VarsTrueFace_V(iP_I) + &
                    sign(1.0, FaceState_V(iP_I) - VarsTrueFace_V(iP_I))*&
                    min(abs(FaceState_V(iP_I) - VarsTrueFace_V(iP_I)),&
                    PressureJumpLimit*VarsTrueFace_V(iP_I))
               if(UseAnisoPressure) &
                    VarsGhostFace_V(Ppar_) = VarsTrueFace_V(Ppar_) + &
                    sign(1.0, FaceState_V(Ppar_) - VarsTrueFace_V(Ppar_))*&
                    min(abs(FaceState_V(Ppar_) - VarsTrueFace_V(Ppar_)),&
                    PressureJumpLimit*VarsTrueFace_V(Ppar_))
            else
               ! Use floating BC for pressure
               ! (correct for zero radial velocity)
               VarsGhostFace_V(iP_I) = VarsTrueFace_V(iP_I)
               if(UseAnisoPressure) &
                    VarsGhostFace_V(Ppar_) = VarsTrueFace_V(Ppar_)
            end if

            ! Change sign for velocities (plasma frozen into dipole field)
            VarsGhostFace_V(iUx_I) = -VarsTrueFace_V(iUx_I)
            VarsGhostFace_V(iUy_I) = -VarsTrueFace_V(iUy_I)
            VarsGhostFace_V(iUz_I) = -VarsTrueFace_V(iUz_I)

            !---------------------------------------------------
            ! Ionosphere outflow in multifluids  --- Yiqun 2008
            !---------------------------------------------------

            if(TypeBc == 'ionosphereoutflow')then      

               iIonSecond = min(IonFirst_ + 1, IonLast_)

               if (TypeCoordSystem /= 'SMG') then 
                  SmgFaceCoords_D = matmul(transform_matrix(TimeBc, &
                       TypeCoordSystem, 'SMG'), FaceCoords_D)
               else 
                  SmgFaceCoords_D = FaceCoords_D
               endif

               SinLatitudeCap = sin(LatitudeCap * cDegToRad)
               zCap = sqrt(sum(SmgFaceCoords_D**2))*SinLatitudeCap

               if(abs(SmgFaceCoords_D(z_)) > zCap)then
                  ! for the polar region
                  if(UseIe .and. UseMultiIon) then 

                     if (TypeCoordSystem /= 'GEO') then 
                        GeoFaceCoords_D = matmul(transform_matrix(TimeBc, &
                             TypeCoordSystem, 'GEO'), FaceCoords_D)
                     else 
                        GeoFaceCoords_D = FaceCoords_D
                     endif
                     GseToGeo_D = matmul(transform_matrix(TimeBc,'GSE','GEO'),&
                          (/0,0,1/))

                     ! For the cap region (refer to Tom Moore 2003?)
                     ! Get the Op flux from IE calculation
                     ! Get the Hp flux from fluxpw, 
                     ! which is constant for certain solar zenith angle
                     ! Fix the velocities(V), thermal energies
                     ! Get the densities(rho), thermal pressure(P)

                     ! get the magnetic field
                     call get_planet_field(TimeBc, FaceCoords_D,&
                          TypeCoordSystem//'NORM', bFace_D)
                     b =  sqrt(sum(bFace_D**2))
                     bUnit_D = bFace_D / B


                     ! get the magnetic field at 4000km 
                     call map_planet_field(TimeBc, FaceCoords_D, &
                          TypeCoordSystem//'NORM', &
                          (4000.0+rPlanet_I(Earth_))/rPlanet_I(Earth_), &
                          XyzMap_D, iHemisphere)
                     call get_planet_field(TimeBc, XyzMap_D, &
                          TypeCoordSystem//'NORM', bFace_D)

                     b4 =  sqrt(sum(bFace_D**2))

                     ! get the joule heating mapped from the ionosphere 
                     ! (already in nomalized unit)
                     call map_inner_bc_jouleheating(TimeBc, FaceCoords_D, &
                          JouleHeating)

                     ! get the O+ flux based on Strangeway's formula, 
                     ! and scale it
                     FluxIono = 2.142e7*(JouleHeating * No2Si_V(UnitPoynting_)&
                          * 1.0e3)**1.265 * 1.0e4 &
                          * Si2No_V(UnitU_) * Si2No_V(UnitN_) * (b4/b)**0.265

                     ! thermal energy = 0.1 + 1.6 * S^1.26 
                     ! (S is joule heating in mW/m^2 at inner boundary)
                     ! to specify the O+ temperature
                     eCap = 0.1 + 9.2 * &
                          ((b4/b)* JouleHeating*No2Si_V(UnitPoynting_) &
                          * 1.0e3)**0.35    !eV

                     ! Get the field aligned current at this location, 
                     ! so comes the parallel energy
                     call get_point_data(1.0, FaceCoords_D, 1, nBlock, Bx_, &
                          nVar+3, State_V)

                     Jlocal_D = State_V(nVar+1:nVar+3)
                     Jpar = sum(bUnit_D * Jlocal_D)  !in normalized unit

                     ! parallel energy 
                     ! (ePar = eV=e*(1500[V/mmA/m^2] * (J//-0.33)^2 [mmA/m^2]))
                     if(abs(Jpar*No2Si_V(UnitJ_))*1.0e6 > 0.33)then 
                        ePar = 1500 * (abs(Jpar*No2Si_V(UnitJ_))*1.0e6 &
                             - 0.33)**2 !eV
                     else
                        ePar = 0.
                     end if

                     ! Get the velocity along B, 
                     ! superpose the parallel velocity and the thermal velocity
                     Ub_V(1) = (sqrt(2 * (ePar + eCap) * cElectronCharge / &
                          (MassFluid_I(IonFirst_)*cProtonMass))) &
                          * Si2No_V(UnitU_)

                     Ub_V(2) = (sqrt(2 * (ePar + eCap) * cElectronCharge / &
                          (MassFluid_I(iIonSecond)*cProtonMass))) &
                          * Si2No_V(UnitU_)

                     ! .OR. Pick the constant velocities and thermal energy
                     ! Ub_V(2) = 10*Io2No_V(UnitU_)  !20km/s
                     ! Ub_V(1) = 20*Io2No_V(UnitU_) 

                     ! SZA x is determind by 
                     ! cosx = sin(the)sin(da)+cos(the)cos(da)cos(dt)
                     ! where, the is the latitude, da is solar declination( 
                     ! angle between solar ray and equatorial plane), 
                     ! dt is local time angle
                     TheTmp = asin(GeoFaceCoords_D(z_)/ &
                          sqrt(sum(GeoFaceCoords_D**2)))      !latitutde
                     DaTmp = acos(GseToGeo_D(z_))               !declination
                     DtTmp = acos(SmgFaceCoords_D(x_)/ &
                          sqrt(SmgFaceCoords_D(x_)**2 + &
                          SmgFaceCoords_D(y_)**2))        !local time angle

                     if(SmgFaceCoords_D(y_)<0.0) DtTmp =  cTwoPi - DtTmp
                     Cosx = sin(TheTmp)*sin(DaTmp) + &
                          cos(TheTmp)*cos(DaTmp)*cos(DtTmp)

                     ! get the magnetic field at 1000km
                     call map_planet_field(TimeBc, FaceCoords_D, &
                          TypeCoordSystem//'NORM', &
                          (1000.0+rPlanet_I(Earth_))/rPlanet_I(Earth_), &
                          XyzMap_d, iHemisphere)
                     call get_planet_field(TimeBc, XyzMap_D, &
                          TypeCoordSystem//'NORM', bFace_D)
                     b1 =  sqrt(sum(bFace_D**2))

                     ! get the Hp flux by mapping the flux at 1000km 
                     ! into the inner boudnary
                     if (acos(Cosx)*cRadToDeg < 90 .and. &
                          acos(Cosx)*cRadToDeg > 0.) then
                        FluxPw = 2.0e8 * 1.0e4 * (b/b1) &
                             * Si2No_V(UnitU_) * Si2No_V(UnitN_)
                     elseif(acos(Cosx)*cRadToDeg < 110) then
                        FluxPw = 2.0*10.**(8-(acos(Cosx)*cRadToDeg &
                             - 90.)/20.*2.5) * 1.0e4 * (b/b1) * &
                             Si2No_V(UnitU_) * Si2No_V(UnitN_)
                     else
                        FluxPw = 2.0*10.**5.5 * 1.0e4 * (b/b1) * &
                             Si2No_V(UnitU_) * Si2No_V(UnitN_)
                     endif

                     ! get the densities
                     VarsGhostFace_V(iRho_I(IonFirst_)) = FluxPw/Ub_V(1) *   &
                          MassFluid_I(IonFirst_)
                     VarsGhostFace_V(iRho_I(iIonSecond)) = FluxIono/Ub_V(2) * &
                          MassFluid_I(iIonSecond)     

                     ! Make sure it points outward
                     if(sum(bUnit_D*FaceCoords_D) < 0.0) bUnit_D = -bUnit_D

                     VarsGhostFace_V(iUx_I(IonFirst_)) = Ub_V(1) * bUnit_D(x_)
                     VarsGhostFace_V(iUy_I(IonFirst_)) = Ub_V(1) * bUnit_D(y_)
                     VarsGhostFace_V(iUz_I(IonFirst_)) = Ub_V(1) * bUnit_D(z_)

                     VarsGhostFace_V(iUx_I(iIonSecond)) = Ub_V(2) * bUnit_D(x_)
                     VarsGhostFace_V(iUy_I(iIonSecond)) = Ub_V(2) * bUnit_D(y_)
                     VarsGhostFace_V(iUz_I(iIonSecond)) = Ub_V(2) * bUnit_D(z_)

                     ! get the pressure
                     VarsGhostFace_V(iP_I(iIonSecond))   =  2./3. * eCap * &
                          cElectronCharge / cBoltzmann &
                          * Si2No_V(UnitTemperature_)  &
                          * VarsGhostFace_V(iRho_I(iIonSecond)) &
                          /MassFluid_I(iIonSecond)
                     VarsGhostFace_V(iP_I(IonFirst_))   =  2./3. * eCap * &
                          cElectronCharge / cBoltzmann & 
                          * Si2No_V(UnitTemperature_)  &
                          * VarsGhostFace_V(iRho_I(IonFirst_)) &
                          /MassFluid_I(IonFirst_)

                     ! for the 'all' fluid
                     VarsGhostFace_V(Rho_) = sum(VarsGhostFace_V( &
                          iRho_I(IonFirst_:iIonSecond)))
                     VarsGhostFace_V(iUx_I(1))  = sum(VarsGhostFace_V( &
                          iRho_I(IonFirst_:iIonSecond)) &
                          * VarsGhostFace_V(iUx_I(IonFirst_:iIonSecond)))&
                          /sum(VarsGhostFace_V(iRho_I(IonFirst_:iIonSecond))) 
                     VarsGhostFace_V(iUy_I(1))  = sum(VarsGhostFace_V( &
                          iRho_I(IonFirst_:iIonSecond)) &
                          * VarsGhostFace_V(iUy_I(IonFirst_:iIonSecond)))&
                          /sum(VarsGhostFace_V(iRho_I(IonFirst_:iIonSecond)))
                     VarsGhostFace_V(iUz_I(1))  = sum(VarsGhostFace_V( &
                          iRho_I(IonFirst_:iIonSecond)) &
                          * VarsGhostFace_V(iUz_I(IonFirst_:iIonSecond)))&
                          /sum(VarsGhostFace_V(iRho_I(IonFirst_:iIonSecond)))
                     VarsGhostFace_V(P_)        = sum(VarsGhostFace_V( &
                          iP_I(IonFirst_:iIonSecond)))

                  else
                     call stop_mpi('ionosphereoutflow should have IE '// &
                          'coupled and multifluids')
                  end if
               end if ! polar cap region
            end if ! ionosphereoutflow type of innerboundary

         end if

      case('buffergrid')
         ! REVISION: June  2011 - R. Oran - generalized.   
         ! Inner boundary conditions based on coupling to another BATSRUS
         ! component through a buffer grid. 

         ! Allows specifying boundary conditions for additional variables
         ! that were not passed through the buffer.
         ! Coupling flags for these options are set in CON_couple_all.f90.

         ! Note: the older case 'coronatoih' is redirected to here for
         ! backwards compatability.

         !Get interpolated values from buffer grid:
         call get_from_spher_buffer_grid(&
              FaceCoords_D,nVar,FaceState_V)

         VarsGhostFace_V = FaceState_V
         if(UseB0)VarsGhostFace_V(Bx_:Bz_)=VarsGhostFace_V(Bx_:Bz_) - B0Face_D

         if(.not. IsFullyCoupledFluid .and. nFluid > 1) then
            ! Only variable associated with the main MHD plasma are passed
            ! through the buffer grid. BC's for fluids must be specified
            ! somehow.

            if (DoOhNeutralBc) then   
               ! Get face BCs for neutrals in the outerheliosphere
               ! (based on M. Opher)
               ! Pop I is going through the inner BCs

               ! PopII leaves the domain at a supersonic velocity
               ! (50km/s while for their temperature 1.E5K their C_s=30km/s)
               ! For the transient case of inward flow, use a fraction of ions

               ! Pop III has the velocity and temperature of the ions at inner
               !boundary,  the density is taken to be a fraction of the ions

               VarsGhostFace_V(iRho_I(iFluid):iP_I(iFluid)) = &
                    VarsTrueFace_V(iRho_I(iFluid):iP_I(iFluid))

               do iFluid = 3,nFluid
                  if(sum(VarsTrueFace_V(iRhoUx_I(iFluid):iRhoUz_I(iFluid)) * &
                       FaceCoords_D) <= 0.0)then
                     VarsGhostFace_V(iRho_I(iFluid):iP_I(iFluid)) = &
                          VarsTrueFace_V(iRho_I(iFluid):iP_I(iFluid))
                  else
                     VarsGhostFace_V(iRho_I(iFluid)) = VarsGhostFace_V(Rho_)* &
                          RhoBcFactor_I(iFluid)
                     VarsGhostFace_V(iP_I(iFluid)) = VarsGhostFace_V(p_) *  &
                          RhoBcFactor_I(iFluid)
                     VarsGhostFace_V(iRhoUx_I(iFluid):iRhoUz_I(iFluid)) = &
                          VarsGhostFace_V(Ux_:Uz_) *uBcFactor_I(iFluid)
                  end if
               end do

            else 
               ! If this component is multyfluid and coupled to a single-fluid,
               ! stop with an error
               call CON_stop( &
                    NameSub//': BCs for multifluid must be specified.')
            end if
         end if

         !^CFG IF SECONDBODY BEGIN
      case('Body2Orbit')
         VarsGhostFace_V = FaceState_V
         VarsGhostFace_V(Bx_:Bz_) = VarsGhostFace_V(Bx_:Bz_) - B0Face_D

         ! Setting velocity BCs to be the second body orbital velocity: 
         VarsGhostFace_V(Ux_) = &
              -(cTwoPi*yBody2/OrbitPeriod)*No2Si_V(UnitX_)*Si2No_V(UnitU_) 
         VarsGhostFace_V(Uy_) = &
              (cTwoPi*xBody2/OrbitPeriod)*No2Si_V(UnitX_)*Si2No_V(UnitU_)
         VarsGhostFace_V(Uz_) =  0.0
         !^CFG END SECONDBODY

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
              - FaceCoords_D * sum(FaceCoords_D*uIono_D) / sum(FaceCoords_D**2)

         select case(TypeBc)
         case('reflect','linetied','polarwind','ionosphere', &
              'ionospherefloat', 'ionosphereoutflow')
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
         ! velocity vector uRot_D as a function of the radius-vector
         ! "FaceCoords"
         !/
         call calc_corotation_velocities(FaceCoords_D, uRot_D)

         select case(TypeBc)
         case('reflect','linetied', &
              'ionosphere','ionospherefloat','polarwind','ionosphereoutflow')
            VarsGhostFace_V(iUx_I) = 2*uRot_D(x_) + VarsGhostFace_V(iUx_I)
            VarsGhostFace_V(iUy_I) = 2*uRot_D(y_) + VarsGhostFace_V(iUy_I)
            VarsGhostFace_V(iUz_I) = 2*uRot_D(z_) + VarsGhostFace_V(iUz_I)
         case default
            call stop_mpi('UseRotatingBc is not compatible with TypeBc_I=' &
                 //TypeBc) 
         end select
      end if
    end subroutine set_face

  end subroutine set_face_bc

end module ModFaceBoundary
