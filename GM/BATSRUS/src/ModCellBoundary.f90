!^CFG COPYRIGHT UM
!==============================================================================
module ModCellBoundary

  use ModSetOuterBc

  use ModVarIndexes
  use ModAdvance, ONLY : State_VGB
  use BATL_size, ONLY: nI, nJ, nK, nGI, nGJ, nGK, nG
  implicit none

  SAVE

  private ! except

  ! Public methods
  public :: set_cell_boundary

  ! Local variables
  integer:: iMin, iMax, jMin, jMax, kMin, kMax, iSide

contains

  !============================================================================
  subroutine set_cell_boundary(iBlock, time_now, DoSetEnergy)

    ! Set ghost cells values in State_VGB.
    ! Set ghost cells in Energy_GBI if DoSetEnergy is true.

    use ModProcMH
    use ModMain
    use ModParallel, ONLY: NOBLK, NeiLev
    use ModGeometry, ONLY: far_field_BCs_BLK, MaxBoundary, XyzMin_D
    use BATL_lib, ONLY: IsRzGeometry, IsCylindricalAxis, IsRlonLat
    use ModPhysics
    use ModUser, ONLY: user_set_cell_boundary
    use ModMultiFluid, ONLY: iFluid, nFluid, iRhoUx_I, iRhoUy_I, iRhoUz_I
    use ModEnergy, ONLY: calc_energy

    integer, intent(in) :: iBlock
    real,    intent(in) :: time_now
    logical, intent(in) :: DoSetEnergy

    integer :: iStart, iLast, iVar

    ! Coefficient +1 or -1 for symmetric vs. anti-symmetric BCs
    real:: SymmCoeff_V(nVar)

    integer :: iGhost, jGhost, kGhost, iGhost2, jGhost2, kGhost2
    logical :: DoTest, DoTestMe, IsFound

    character(len=*), parameter :: NameSub = 'set_cell_boundary'
    !--------------------------------------------------------------------------
    if(iBlock==BLKtest.and.iProc==PROCtest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    endif

    if(.not.far_field_BCs_BLK(iBlock))then
       write(*,*) NameSub,' warning: iBlock=',iBlock,' is not far_field block'
       RETURN
    end if

    if(DoTestMe)write(*,*)NameSub,': iBlock, set_E, neiLEV=',&
         iBlock, DoSetEnergy,neiLEV(:,iBlock)

    if(DoTestMe)then
       Ighost =Itest; Jghost=Jtest;  Kghost=Ktest
       Ighost2=Itest; Jghost2=Jtest; Kghost2=Ktest
       select case(DimTest)
       case(x_)
          if(iTest== 1)then
             iGhost=0; iGhost2 = -1
          endif
          if(iTest==nI)then
             iGhost=nI+1; iGhost2 = nI+2
          endif
       case(y_)
          if(jTest== 1)then
             jGhost=0; jGhost2 = -1
          endif
          if(jTest==nJ)then
             jGhost=nJ+1; jGhost2 = nJ+2
          endif
       case(z_)
          if(kTest== 1)then
             kGhost=0; kGhost2 = -1
          endif
          if(kTest==nK)then
             kGhost=nK+1; kGhost2 = nK+2
          endif
       end select

       write(*,*)'iTest,  jTest,  kTest   =',iTest,  jTest,  kTest
       write(*,*)'iGhost, jGhost, kGhost  =',iGhost, jGhost, kGhost
       write(*,*)'iGhost2,jGhost2,kGhost2 =',iGhost2,jGhost2,kGhost2
       do iVar=1,nVar
          write(*,*)'initial',NameVar_V(iVar),   'cell,ghost,ghost2=',&
               State_VGB(iVar,Itest,Jtest,Ktest,iBlock),&
               State_VGB(iVar,Ighost,Jghost,Kghost,iBlock), &
               State_VGB(iVar,Ighost2,Jghost2,Kghost2,iBlock)
       end do
    end if

    iStart=max(MaxBoundary+1,1)
    ! Do not apply cell boundary conditions at the pole 
    ! This is either handled by message passing or supercell
    if(IsRLonLat) then
       iLast = 2
    else
       iLast = 6
    end if

    ! Do not work on ignored directions
    if(nK == 1 .and. iLast > 4) iLast = 4
    if(nJ == 1 .and. iLast > 2) iLast = 2

    do iSide = iStart, iLast

       ! Check if this side of the block is indeed an outer boundary
       if(neiLEV(iSide,iBlock)/=NOBLK) CYCLE

       ! Do not apply cell boundary conditions at the pole 
       ! This is either handled by message passing or supercell
       if(IsCylindricalAxis .and. iSide == 1) CYCLE

       ! Set index limits for the ghost cell range
       iMin = 1-nGI; iMax = nI + nGI
       jMin = 1-nGJ; jMax = nJ + nGJ
       kMin = 1-nGK; kMax = nK + nGK
       select case(iSide)
       case(1)
          iMax = 0
       case(2)
          iMin = nI + 1
       case(3)
          jMax = 0
       case(4)
          jMin = nJ + 1
       case(5)
          kMax = 0
       case(6)
          kMin = nK + 1
       end select

       select case(TypeBc_I(iSide))
       case('coupled')
          ! For SC-IH coupling the extra wave energy variable needs a BC
          if(NameThisComp == 'SC') call BC_cont(ScalarFirst_, ScalarLast_)
       case('periodic')
          call stop_mpi('The neighbors are not defined at periodic boundaries')
       case('float','outflow')
          call BC_cont(1, nVar)
          if(UseOutflowPressure .and. TypeBc_I(iSide) == 'outflow') &
               call  BC_fixed(p_,p_,(/ pOutflow /) )
          !^CFG IF IMPLICIT BEGIN
          if(UseRadDiffusion) &
               call set_radiation_outflow_bc(WaveFirst_, WaveLast_, iSide)
          !^CFG END IMPLICIT
       case('raeder')
          call BC_cont(1, nVar)
          if(iSide==4.or.iSide==3)then
             call BC_fixed(By_,By_,DefaultState_V(By_))
          elseif(iSide==5.or.iSide==6)then
             call BC_fixed(Bz_,Bz_,DefaultState_V(Bz_))
          end if
       case('reflect')
          ! Scalars are symmetric
          SymmCoeff_V = 1.0
          ! Normal vector components are mirror symmetric
          select case(iSide)
          case(1,2)
             SymmCoeff_V(iRhoUx_I(1:nFluid)) = -1.0
             if(UseB)SymmCoeff_V(Bx_) = -1.0
          case(3,4)
             SymmCoeff_V(iRhoUy_I(1:nFluid)) = -1.0
             if(UseB)SymmCoeff_V(By_) = -1.0
             ! For RZ geometry, mirror Z components too at the axis
             if(IsRzGeometry .and. XyzMin_D(2)==0.0 .and. iSide==3)then
                SymmCoeff_V(iRhoUz_I(1:nFluid)) = -1.0
                if(UseB)SymmCoeff_V(Bz_) = -1.0
             end if
          case(5,6)
             SymmCoeff_V(iRhoUz_I(1:nFluid)) = -1.0
             if(UseB)SymmCoeff_V(Bz_) = -1.0
          end select
          call BC_symm(1, nVar, SymmCoeff_V)
       case('linetied')
          ! Most variables float
          call BC_cont(1,nVar)

          ! The density and momentum use symmetric/antisymmetric BCs
          SymmCoeff_V = 1.0
          do iFluid=1, nFluid
             SymmCoeff_V(iRhoUx_I(iFluid):iRhoUz_I(iFluid)) = -1.0
             call BC_symm(iRho_I(iFluid), iRhoUz_I(iFluid), SymmCoeff_V)
          end do
       case('fixed','inflow','vary','ihbuffer')
          if(time_accurate &
               .and.(TypeBc_I(iSide)=='vary'.or.TypeBc_I(iSide)=='inflow'))then
             call BC_solar_wind(time_now)
          else if(TypeBc_I(iSide)=='ihbuffer'.and.time_loop)then
             call BC_solar_wind_buffer
          else
             call BC_fixed(1,nVar,CellState_VI(:,iSide))
             if(UseB0)call BC_fix_B0(Bx_,Bz_)
          end if
       case('fixedB1','fixedb1')
          call BC_fixed(1,nVar,CellState_VI(:,iSide))
       case('shear')
          call BC_shear
       case('none')
       case default
          IsFound=.false.
          if(UseUserOuterBcs .or. TypeBc_I(iSide) == 'user')&
               call user_set_cell_boundary(iBlock,iSide,TypeBc_I(iSide),IsFound)

          if(.not. IsFound) call stop_mpi( &
               NameSub // ': unknown TypeBc_I=' //TypeBc_I(iSide))
       end select

       if(DoSetEnergy) call calc_energy( &
            iMin, iMax, jMin, jMax, kMin, kMax, iBlock, 1, nFluid)

    end do

    if(DoTestMe)then
       do iVar=1,nVar
          write(*,*)'final',NameVar_V(iVar),'   cell,ghost,ghost2=',&
               State_VGB(iVar,Itest,Jtest,Ktest,iBlock),&
               State_VGB(iVar,Ighost,Jghost,Kghost,iBlock),&
               State_VGB(iVar,Ighost2,Jghost2,Kghost2,iBlock)
       end do
    end if

  contains

    !============================================================================
    subroutine BC_cont(iVarStart, iVarLast)

      ! Continuous: ghost = phys1
      integer, intent(in) :: iVarStart,iVarLast

      integer:: i, j, k
      !--------------------------------------------------------------------------

      select case(iSide)
      case(1)
         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            State_VGB(iVarStart:iVarLast,i,j,k,iBlock) = State_VGB(iVarStart:iVarLast,1,j,k,iBlock)
         end do; end do; end do
      case(2)
         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            State_VGB(iVarStart:iVarLast,i,j,k,iBlock) = State_VGB(iVarStart:iVarLast,nI,j,k,iBlock)
         end do; end do; end do
      case(3)
         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            State_VGB(iVarStart:iVarLast,i,j,k,iBlock) = State_VGB(iVarStart:iVarLast,i,1,k,iBlock)
         end do; end do; end do
      case(4)
         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            State_VGB(iVarStart:iVarLast,i,j,k,iBlock) = State_VGB(iVarStart:iVarLast,i,nJ,k,iBlock)
         end do; end do; end do
      case(5)
         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            State_VGB(iVarStart:iVarLast,i,j,k,iBlock) = State_VGB(iVarStart:iVarLast,i,j,1,iBlock)
         end do; end do; end do
      case(6)
         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            State_VGB(iVarStart:iVarLast,i,j,k,iBlock) = State_VGB(iVarStart:iVarLast,i,j,nK,iBlock)
         end do; end do; end do
      end select

    end subroutine BC_cont

    !============================================================================
    subroutine BC_shear

      ! Shear: ghost = phys(shifted)

      integer :: i, j, k, Dn
      !--------------------------------------------------------------------------
      ! For the corners and boundaries 5 and 6 fill with unsheared data first
      call BC_cont(1, nVar)

      ! If the shock is not tilted, there is nothing to do
      if( ShockSlope == 0.0) RETURN

      ! Shear according to ShockSlope
      if(ShockSlope < 0.0)then
         call stop_mpi('ShockSlope must be positive!')
      elseif(ShockSlope >= 1.0)then
         Dn = nint(ShockSlope)
         if(abs(Dn-ShockSlope)>cTiny)&
              call stop_mpi('ShockSlope > 1 should be a round number!')
         select case(iSide)
         case(1)
            ! Shift by (Dix,Diy) = (Dn,-1)
            do k = kMin, kMax; do j = jMin+1, jMax; do i = iMax, iMin, -1
               State_VGB(:,i,j,k,iBlock) = State_VGB(:,i+Dn,j-1,k,iBlock)
            end do; end do; end do
         case(2)
            ! Shift by (Dix,Diy) = (-Dn,+1)
            do k = kMin, kMax; do j = jMin, jMax-1; do i = iMin, iMax
               State_VGB(:,i,j,k,iBlock) = State_VGB(:,i-Dn,j+1,k,iBlock)
            end do; end do; end do
         case(3)
            ! Shift by (Dix,Diy) = (-Dn,+1)
            do k = kMin, kMax; do j = jMax, jMin, -1; do i = iMin+Dn, iMax
               State_VGB(:,i,j,k,iBlock) = State_VGB(:,i-Dn,j+1,k,iBlock)
            end do; end do; end do
         case(4)
            ! Shift by (Dix,Diy) = (+Dn,-1)
            do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax-Dn
               State_VGB(:,i,j,k,iBlock) = State_VGB(:,i+Dn,j-1,k,iBlock)
            end do; end do; end do
         end select
      else
         ! ShockSlope < 1
         Dn = nint(1/ShockSlope)
         if(abs(Dn - 1/ShockSlope) > 1e-6)call stop_mpi( &
              'ShockSlope < 1 should be the inverse of a round number!')
         select case(iSide)
         case(1)
            ! Shift by (Dix,Diy) = (1,-Dn)
            do k = kMin, kMax; do j = jMin+Dn, jMax; do i = iMax, iMin, -1
               State_VGB(:,i,j,k,iBlock) = State_VGB(:,i+1,j-Dn,k,iBlock)
            end do; end do; end do
         case(2)
            ! Shift by (Dix,Diy) = (-1,+Dn)
            do k = kMin, kMax; do j = jMin, jMax-Dn; do i = iMin, iMax
               State_VGB(:,i,j,k,iBlock) = State_VGB(:,i-1,j+Dn,k,iBlock)
            end do; end do; end do
         case(3)
            ! Shift by (Dix,Diy) = (-1,+Dn)
            do k = kMin, kMax; do j = jMax, jMin, -1; do i = iMin+1, iMax
               State_VGB(:,i,j,k,iBlock) = State_VGB(:,i-1,j+Dn,k,iBlock)
            end do; end do; end do
         case(4)
            ! Shift by (Dix,Diy) = (+1,-Dn)
            do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax-1
               State_VGB(:,i,j,k,iBlock) = State_VGB(:,i+1,j-Dn,k,iBlock)
            end do; end do; end do
         end select
      end if

    end subroutine BC_shear

    !============================================================================
    subroutine BC_symm(iVarStart, iVarLast, Coeff_V)

      ! Symmetry with optional sign change: ghost_i = Coeff*phys_i

      integer, intent(in) :: iVarStart, iVarLast
      real, intent(in):: Coeff_V(iVarStart:iVarLast)

      integer:: i, j, k
      !--------------------------------------------------------------------------
      select case(iSide)
      case(1)
         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            State_VGB(iVarStart:iVarLast,i,j,k,iBlock) = &
                 Coeff_V*State_VGB(iVarStart:iVarLast,1-i,j,k,iBlock)
         end do; end do; end do
      case(2)
         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            State_VGB(iVarStart:iVarLast,i,j,k,iBlock) = &
                 Coeff_V*State_VGB(iVarStart:iVarLast,2*nI+1-i,j,k,iBlock)
         end do; end do; end do
      case(3)
         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            State_VGB(iVarStart:iVarLast,i,j,k,iBlock) = &
                 Coeff_V*State_VGB(iVarStart:iVarLast,i,1-j,k,iBlock)
         end do; end do; end do
      case(4)
         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            State_VGB(iVarStart:iVarLast,i,j,k,iBlock) = &
                 Coeff_V*State_VGB(iVarStart:iVarLast,i,2*nJ+1-j,k,iBlock)
         end do; end do; end do
      case(5)
         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            State_VGB(iVarStart:iVarLast,i,j,k,iBlock) = &
                 Coeff_V*State_VGB(iVarStart:iVarLast,i,j,1-k,iBlock)
         end do; end do; end do
      case(6)
         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            State_VGB(iVarStart:iVarLast,i,j,k,iBlock) = &
                 Coeff_V*State_VGB(iVarStart:iVarLast,i,j,2*nK+1-k,iBlock)
         end do; end do; end do
      end select

    end subroutine BC_symm

    !============================================================================
    subroutine BC_fixed(iVarStart, iVarLast, State_V)

      ! ghost = State_V

      integer, intent(in) :: iVarStart, iVarLast
      real,    intent(in) :: State_V(iVarStart:iVarLast)

      integer:: i, j, k
      !--------------------------------------------------------------------------

      do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
         State_VGB(iVarStart:iVarLast,i,j,k,iBlock) = State_V
      end do; end do; end do

    end subroutine BC_fixed

    !============================================================================
    subroutine BC_fix_B0(iVarStart, iVarLast)

      use ModB0, ONLY: B0_DGB

      integer, intent(in) :: iVarStart,iVarLast

      ! Set B = B - B0 in ghost cells

      integer:: i, j, k
      !--------------------------------------------------------------------------
      do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
         State_VGB(iVarStart:iVarLast,i,j,k,iBlock) =  &
              State_VGB(iVarStart:iVarLast,i,j,k,iBlock) - B0_DGB(:,i,j,k,iBlock)
      end do; end do; end do

    end subroutine BC_fix_B0

    !============================================================================

    subroutine BC_solar_wind(time_now)

      use ModVarIndexes
      use ModGeometry,   ONLY: x2
      use ModB0,         ONLY: B0_DGB
      use ModMultiFluid, ONLY: &
           iRho_I, iUx_I, iUy_I, iUz_I, iRhoUx_I, iRhoUy_I, iRhoUz_I
      use ModSolarwind,  ONLY: get_solar_wind_point
      use ModMain,       ONLY: UseB0, y_, z_
      use BATL_lib,      ONLY: Xyz_DGB

      ! Current simulation time in seconds
      real, intent(in) :: time_now 

      ! index and location of a single point
      integer :: i, j, k
      real :: x, y, z
      ! Varying solar wind parameters
      real :: SolarWind_V(nVar)

      !logical :: DoTest, DoTestMe
      !character(len=*), parameter:: NameSub = 'BC_solar_wind'
      !--------------------------------------------------------------------------

      !call set_oktest(NameSub, DoTest, DoTestMe)

      do k = kMin, kMax
         z = Xyz_DGB(z_,1,1,k,iBlock)
         do j = jMin, jMax
            y = Xyz_DGB(y_,1,j,1,iBlock)
            do i = iMin, iMax

               ! x= Xyz_DGB(x_,i,j,k,iBlock) ! for cell based BC this would be best
               x=x2 ! for face based and for west side as the inflow
               call get_solar_wind_point(time_now, (/x, y, z/), SolarWind_V)

               State_VGB(:,i,j,k,iBlock) = SolarWind_V

               ! Convert velocities to momenta
               State_VGB(iRhoUx_I, i,j,k,iBlock) = &
                    State_VGB(iUx_I, i,j,k,iBlock)*State_VGB(iRho_I,i,j,k,iBlock)
               State_VGB(iRhoUy_I, i,j,k,iBlock) = &
                    State_VGB(iUy_I, i,j,k,iBlock)*State_VGB(iRho_I,i,j,k,iBlock)
               State_VGB(iRhoUz_I, i,j,k,iBlock) = &
                    State_VGB(iUz_I, i,j,k,iBlock)*State_VGB(iRho_I,i,j,k,iBlock)

               ! Subtract B0:   B1 = B - B0
               if(UseB0) State_VGB(Bx_:Bz_,i,j,k,iBlock)    = &
                    State_VGB(Bx_:Bz_,i,j,k,iBlock) - B0_DGB(:,i,j,k,iBlock)
            end do
         end do
      end do

    end subroutine BC_solar_wind

    !============================================================================

    subroutine BC_solar_wind_buffer

      use ModSize, ONLY: y_, z_
      use ModVarIndexes, ONLY: Bx_, Bz_
      use ModAdvance, ONLY : State_VGB, B0_DGB
      use BATL_lib, ONLY: Xyz_DGB

      ! index and location of a single point
      integer :: i, j, k
      real    :: y, z
      !--------------------------------------------------------------------------
      do k = kMin, kMax
         z = Xyz_DGB(z_,1,1,k,iBlock)
         do j = jMin, jMax
            y = Xyz_DGB(y_,1,j,1,iBlock)
            do i = iMin, iMax
               call read_ih_buffer(y,z,State_VGB(:,i,j,k,iBlk))
               ! Subtract B0
               State_VGB(Bx_:Bz_,i,j,k,iBlock) = State_VGB(Bx_:Bz_,i,j,k,iBlock) &
                    - B0_DGB(:,i,j,k,iBlock)
            end do
         end do
      end do

    end subroutine BC_solar_wind_buffer
    !============================================================================

    subroutine set_radiation_outflow_bc(iVarFirst, iVarLast, iSide)

      use ModAdvance,  ONLY: State_VGB, nI, nJ, nK, nWave
      use ModPhysics,  ONLY: Si2No_V, UnitX_
      use ModSize,     ONLY: x_, y_, z_
      use ModUser,     ONLY: user_material_properties
      use BATL_lib,    ONLY: CellSize_DB

      integer, intent(in) :: iVarFirst, iVarLast, iSide

      integer:: iVar, i, j, k, iWave
      integer:: i1G, j1G, k1G ! index for first ghost cell layer
      integer:: i2G, j2G, k2G ! index for last  ghost cell layer
      integer:: i1P, j1P, k1P ! index for first physical cell layer
      integer:: i2P, j2P, k2P ! index for second physical cell layer
      integer:: Di, Dj, Dk    ! +1 if last ghost cell layer has larger index 
      !                         than first ghost cell layer and -1 otherwise

      real   :: OpacityRosselandSi_W(nWave), Coef
      !--------------------------------------------------------------------------

      select case(iSide)
      case(1,2)
         if(iSide == 1)then
            i1P = 1; i2P = 2; i1G = 0; i2G = iMin; Di = -1
         else
            i1P = nI; i2P = nI-1; i1G = nI+1; i2G = iMax; Di = +1
         end if
         do k = 1, nK; do j = 1, nJ
            call user_material_properties(State_VGB(:,i1P,j,k,iBlock), &
                 i1P, j, k, iBlock, OpacityRosselandOut_W=OpacityRosselandSi_W)

            do iVar = iVarFirst, iVarLast
               iWave = iVar - iVarFirst + 1
               Coef = 2/sqrt( (3*OpacityRosselandSi_W(iWave) &
                    /          Si2No_V(UnitX_) * CellSize_DB(x_,iBlock))**2 &
                    + ((State_VGB(iVar,i2p,j,k,iBlock) &
                    -   State_VGB(iVar,i1p,j,k,iBlock)) &
                    /  State_VGB(iVar,i1p,j,k,iBlock))**2)
               State_VGB(iVar,i1G,j,k,iBlock) = State_VGB(iVar,i1P,j,k,iBlock) &
                    *(Coef - 0.5)/(Coef + 0.5)
               ! Extrapolate linearly to the rest of the ghost cells
               do i = i1G + Di, i2G, Di
                  State_VGB(iVar,i,j,k,iBlock) &
                       = 2*State_VGB(iVar,i-Di,j,k,iBlock) &
                       -   State_VGB(iVar,i-2*Di,j,k,iBlock)
               end do
            end do
         end do; end do

      case(3,4)
         if(iSide == 1)then
            j1P = 1; j2P = 2; j1G = 0; j2G = jMin; Dj = -1
         else
            j1P = nJ; j2P = nJ-1; j1G = nJ+1; j2G = jMax; Dj = +1
         end if
         do k = 1, nK; do i = 1, nI
            call user_material_properties(State_VGB(:,i,j1p,k,iBlock), &
                 i, j1p, k, iBlock, OpacityRosselandOut_W=OpacityRosselandSi_W)

            do iVar = iVarFirst, iVarLast
               iWave = iVar - iVarFirst + 1
               Coef = 2/sqrt( (3*OpacityRosselandSi_W(iWave) &
                    /          Si2No_V(UnitX_) * CellSize_DB(y_,iBlock))**2 &
                    + ((State_VGB(iVar,i,j2p,k,iBlock) &
                    -   State_VGB(iVar,i,j1p,k,iBlock)) &
                    /  State_VGB(iVar,i,j1p,k,iBlock))**2)
               State_VGB(iVar,i,j1g,k,iBlock) = State_VGB(iVar,i,j1p,k,iBlock) &
                    *(Coef - 0.5)/(Coef + 0.5)
               ! Extrapolate linearly to the rest of the ghost cells
               do j = j1G + Dj, j2G, Dj
                  State_VGB(iVar,i,j,k,iBlock) &
                       = 2*State_VGB(iVar,i,j-Dj,k,iBlock) &
                       -   State_VGB(iVar,i,j-2*Dj,k,iBlock)
               end do
            end do
         end do; end do

      case(5,6)
         if(iSide == 1)then
            k1P = 1; k2P = 2; k1G = 0; k2G = kMin; Dk = -1
         else
            k1P = nK; k2P = nK-1; k1G = nK+1; k2G = kMax; Dk = +1
         end if
         do j = 1, nJ; do i = 1, nI
            call user_material_properties(State_VGB(:,i,j,k1p,iBlock), &
                 i, j, k1p, iBlock, OpacityRosselandOut_W=OpacityRosselandSi_W)

            do iVar = iVarFirst, iVarLast
               iWave = iVar - iVarFirst + 1
               Coef = 2/sqrt( (3*OpacityRosselandSi_W(iWave) &
                    /          Si2No_V(UnitX_) * CellSize_DB(z_,iBlock))**2 &
                    + ((State_VGB(iVar,i,j,k2p,iBlock) &
                    -   State_VGB(iVar,i,j,k1p,iBlock)) &
                    /  State_VGB(iVar,i,j,k1p,iBlock))**2)
               State_VGB(iVar,i,j,k1g,iBlock)=State_VGB(iVar,i,j,k1p,iBlock) &
                    *(Coef - 0.5)/(Coef + 0.5)
               ! Extrapolate linearly to the rest of the ghost cells
               do k = k1G + Dk, k2G, Dk
                  State_VGB(iVar,i,j,k,iBlock) &
                       = 2*State_VGB(iVar,i,j,k-Dk,iBlock) &
                       -   State_VGB(iVar,i,j,k-2*Dk,iBlock)
               end do
            end do
         end do; end do
      end select

    end subroutine set_radiation_outflow_bc

  end subroutine set_cell_boundary


end module ModCellBoundary
