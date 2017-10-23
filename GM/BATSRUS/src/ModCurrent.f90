!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModCurrent

  use ModCoordTransform, ONLY: sph_to_xyz
  use CON_axes,          ONLY: transform_matrix

  implicit none

contains
  !==============================================================================

  subroutine get_point_data(WeightOldState, XyzIn_D, iBlockMin, iBlockMax, &
       iVarMin, iVarMax, StateCurrent_V)

    ! Interpolate the (new and/or old) state vector from iVarMin to iVarMax and 
    ! the current (if iVarMax=nVar+3) for input position 
    ! XyzIn_D given in Cartesian coordinates. The interpolated state 
    ! is second order accurate everywhere except where there is a 
    ! resolution change in more than one direction for the cell centers 
    ! surrounding the given position. In these exceptional cases the 
    ! interpolated state is first order accurate. The interpolation algorithm
    ! is based on trilinear interpolation, but it is generalized for
    ! trapezoidal hexahedrons.

    use ModProcMH, ONLY: iProc
    use ModVarIndexes, ONLY: nVar
    use ModMain, ONLY: nI, nJ, nK, nIJK_D, Unused_B, BlkTest, ProcTest
    use ModAdvance, ONLY: State_VGB, StateOld_VGB
    use ModParallel, ONLY: NeiLev
    use ModGeometry, ONLY: XyzStart_BLK
    use BATL_lib, ONLY: IsCartesianGrid, CellSize_DB, xyz_to_coord

    ! Weight of the old state in the interpolation
    real, intent(in)  :: WeightOldState

    ! Input position is in generalized coordinates
    real, intent(in)  :: XyzIn_D(3)

    ! Block index range (typically 1:nBlock or iBlock:iBlock)
    integer, intent(in) :: iBlockMin, iBlockMax

    ! Do we need to calculate currents
    integer, intent(in) :: iVarMin, iVarMax

    ! Weight and interpolated state at the input position
    real, intent(out) :: StateCurrent_V(0:iVarMax-iVarMin+1)

    ! Local variables

    ! Maximum index for state variable and number of state variables
    integer :: iStateMax, nState

    ! Position in generalized coordinates
    real :: Xyz_D(3)

    ! Cell size and buffer size for current block
    real,    dimension(3) :: Dxyz_D, DxyzInv_D, DxyzLo_D, DxyzHi_D

    ! Position of cell center to the lower index direction
    integer, dimension(3) :: IjkLo_D 

    ! Position of point and current cell center
    real :: x, y, z, xI, yJ, zK

    ! Bilinear weights
    real    :: WeightX,WeightY,WeightZ,WeightXyz

    ! Dimension, cell, block index and MPI error code
    integer :: iDim,i,j,k,iLo,jLo,kLo,iHi,jHi,kHi,iBlock

    ! Current at the cell center
    real:: Current_D(3)

    ! Testing
    logical :: DoTest, DoTestMe
    !----------------------------------------------------------------------------
    if(iProc == ProcTest .and. iBlockMin <= BlkTest .and. iBlockMax >= BlkTest)then
       call set_oktest('get_point_data', DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if

    ! Calculate maximum index and the number of state variables
    iStateMax = min(iVarMax, nVar)
    nState    = iStateMax - iVarMin + 1

    ! Convert to generalized coordinates if necessary
    if(IsCartesianGrid)then
       Xyz_D = XyzIn_D
    else
       call xyz_to_coord(XyzIn_D, Xyz_D)
    end if

    ! Set state and weight to zero, so MPI_reduce will add it up right
    StateCurrent_V = 0.0

    ! Loop through all blocks
    BLOCK: do iBlock = iBlockMin, iBlockMax
       if(Unused_B(iBlock)) CYCLE

       if(DoTestMe)write(*,*)'get_point_data called with XyzIn_D=',XyzIn_D

       ! Put cell size of current block into an array
       Dxyz_D = CellSize_DB(:,iBlock)

       ! Set buffer zone according to relative size of neighboring block
       do iDim = 1, 3
          ! Block at the lower index side
          select case(NeiLev(2*iDim-1,iBlock))
          case(1)
             DxyzLo_D(iDim) = 1.5*Dxyz_D(iDim)
          case(-1)
             DxyzLo_D(iDim) = 0.75*Dxyz_D(iDim)
          case default
             DxyzLo_D(iDim) = Dxyz_D(iDim)
          end select
          ! Check if point is inside the buffer zone on the lower side
          if(Xyz_D(iDim)<XyzStart_BLK(iDim,iBlock) - DxyzLo_D(iDim)) CYCLE BLOCK

          ! Block at the upper index side
          select case(NeiLev(2*iDim,iBlock))
          case(1)
             DxyzHi_D(iDim) = 1.5*Dxyz_D(iDim)
          case(-1)
             DxyzHi_D(iDim) = 0.75*Dxyz_D(iDim)
          case default
             DxyzHi_D(iDim) = Dxyz_D(iDim)
          end select
          ! Check if point is inside the buffer zone on the upper side
          if(Xyz_D(iDim) > XyzStart_BLK(iDim,iBlock) + &
               (nIJK_D(iDim)-1)*Dxyz_D(iDim) + DxyzHi_D(iDim)) CYCLE BLOCK
       end do

       ! Find closest cell center indexes towards the lower index direction
       IjkLo_D = floor((Xyz_D - XyzStart_BLK(:,iBlock))/Dxyz_D)+1

       ! Set the size of the box for bilinear interpolation

       ! At resolution change the box size is the sum
       ! average of the cell size of the neighboring blocks

       ! Also make sure that IjkLo_D is not out of bounds
       do iDim = 1,3
          if(IjkLo_D(iDim) < 1)then
             IjkLo_D(iDim)   = 0
             DxyzInv_D(iDim) = 1/DxyzLo_D(iDim)
          elseif(IjkLo_D(iDim) >= nIJK_D(iDim))then
             IjkLo_D(iDim)   = nIJK_D(iDim)
             DxyzInv_D(iDim) = 1/DxyzHi_D(iDim)
          else
             DxyzInv_D(iDim) = 1/Dxyz_D(iDim)
          end if
       end do

       if(DoTest)then
          write(*,*)'Point found at iProc,iBlock,iLo,jLo,kLo=',&
               iProc,iBlock,IjkLo_D
          write(*,*)'iProc, XyzStart_BLK,Dx=', iProc, &
               XyzStart_BLK(:,iBlock), Dxyz_D(1)
       end if

       ! Set the index range for the physical cells
       iLo = max(IjkLo_D(1),1)
       jLo = max(IjkLo_D(2),1)
       kLo = max(IjkLo_D(3),1)
       iHi = min(IjkLo_D(1)+1,nI)
       jHi = min(IjkLo_D(2)+1,nJ)
       kHi = min(IjkLo_D(3)+1,nK)

       ! Put the point position into scalars
       x = Xyz_D(1); y = Xyz_D(2); z = Xyz_D(3)

       ! Loop through the physical cells to add up their contribution
       do k = kLo, kHi
          zK = XyzStart_BLK(3,iBlock) + (k-1)*Dxyz_D(3)
          WeightZ = 1 - DxyzInv_D(3)*abs(z-zK)
          do j = jLo, jHi
             yJ = XyzStart_BLK(2,iBlock) + (j-1)*Dxyz_D(2)
             WeightY = 1 - DxyzInv_D(2)*abs(y-yJ)
             do i = iLo, iHi
                xI = XyzStart_BLK(1,iBlock) + (i-1)*Dxyz_D(1)
                WeightX = 1 - DxyzInv_D(1)*abs(x-xI)

                WeightXyz = WeightX*WeightY*WeightZ

                if(WeightXyz>0.0)then
                   ! Add up the weight
                   StateCurrent_V(0) = StateCurrent_V(0) + WeightXyz

                   ! Add contibutions from the old state if required
                   if(WeightOldState > 0.0) &
                        StateCurrent_V(1:nState) = StateCurrent_V(1:nState) + &
                        WeightXyz * WeightOldState * &
                        StateOld_VGB(iVarMin:iStateMax,i,j,k,iBlock)

                   ! Add contibutions from the current state if required
                   if(WeightOldState < 1.0) &
                        StateCurrent_V(1:nState) = StateCurrent_V(1:nState) + &
                        WeightXyz * (1 - WeightOldState) * &
                        State_VGB(iVarMin:iStateMax,i,j,k,iBlock)

                   ! The current is always based on the new state
                   if(iVarMax == nVar + 3)then
                      call get_current(i, j, k, iBlock, Current_D)
                      StateCurrent_V(nState+1:nState+3) = &
                           StateCurrent_V(nState+1:nState+3) &
                           + WeightXyz * Current_D
                   end if

                   if(DoTest)write(*,*)'Contribution iProc,i,j,k,WeightXyz=',&
                        iProc,i,j,k,WeightXyz
                end if
             end do
          end do
       end do
    end do BLOCK

  end subroutine get_point_data
  !============================================================================

  subroutine get_current(i, j, k, iBlock, Current_D, nOrderResChange, &
       DoIgnoreBody)

    ! Calculate the current in a cell of a block
    ! Avoid using ghost cells at resolution changes.
    ! Avoid using cells inside the body.
    ! If the optional argument nOrderReschange is present and equals 1,
    ! then use first order scheme at resolution changes,
    ! otherwise use second order scheme when possible.

    use ModAdvance,  ONLY: State_VGB, Bx_, By_, Bz_
    use ModGeometry, ONLY: True_Cell, true_BLK, r_BLK
    use BATL_lib, ONLY: IsCartesianGrid, IsRzGeometry, Xyz_DGB, CellSize_DB
    use ModParallel, ONLY: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot
    use ModSize,     ONLY: nI, nJ, nK, x_, y_, z_
    use ModB0, ONLY: UseCurlB0, rCurrentFreeB0, set_b0_source, CurlB0_DC

    integer, intent(in) :: i, j, k, iBlock
    real,    intent(out):: Current_D(3)

    integer, optional, intent(in):: nOrderResChange

    logical, optional, intent(in):: DoIgnoreBody

    logical:: UseFirstOrder
    integer:: iL, iR, jL, jR, kL, kR
    real   :: Ax, Bx, Cx, Ay, By, Cy, Az, Bz, Cz
    real   :: InvDx2, InvDy2, InvDz2

    integer:: iBlockLast = -1
    !--------------------------------------------------------------------------

    ! Exclude body cells
    if(.not.True_Cell(i,j,k,iBlock) .and. .not.present(DoIgnoreBody))then
       Current_D = 0.0

       RETURN
    endif

    UseFirstOrder = .false.
    if(present(nOrderReschange)) UseFirstOrder = nOrderResChange == 1

    InvDx2 = 0.5/CellSize_DB(x_,iBlock)
    InvDy2 = 0.5/CellSize_DB(y_,iBlock)
    InvDz2 = 0.5/CellSize_DB(z_,iBlock)

    ! Central difference
    iR = i+1; iL = i-1;
    Ax = -InvDx2; Bx = 0.0; Cx = +InvDx2
    ! Avoid the ghost cells at resolution changes by using
    ! second-order one-sided difference
    if(i==1 .and. abs(NeiLeast(iBlock))==1)then
       if(UseFirstOrder)then
          iL = i; Ax = -2*InvDx2; Cx = 2*InvDx2
       else
          iL = i+1; iR = i+2; Ax = 4*InvDx2; Bx =-3*InvDx2; Cx =-InvDx2
       end if
    elseif(i==nI .and. abs(NeiLwest(iBlock))==1)then
       if(UseFirstOrder)then
          iR = i; Ax = -2*InvDx2; Cx = 2*InvDx2
       else
          iL = i-1; iR = i-2; Ax =-4*InvDx2; Bx = 3*InvDx2; Cx = InvDx2
       end if
    end if

    ! For y direction
    if(nJ == 1)then
       ! 1D
       jR = j; jL = j
       Ay = 0.0; By = 0.0; Cy = 0.0
    else
       jR = j+1; jL = j-1;
       Ay = -InvDy2; By = 0.0; Cy = +InvDy2
       if(j==1 .and. abs(NeiLsouth(iBlock))==1)then
          if(UseFirstOrder)then
             jL = j; Ay = -2*InvDy2; Cy = 2*InvDy2
          else
             jL = j+1; jR = j+2; Ay = 4*InvDy2; By = -3*InvDy2; Cy = -InvDy2
          end if
       elseif(j==nJ .and. abs(NeiLnorth(iBlock))==1)then
          if(UseFirstOrder)then
             jR = j; Ay = -2*InvDy2; Cy = 2*InvDy2
          else
             jL = j-1; jR = j-2; Ay = -4*InvDy2; By = 3*InvDy2; Cy = InvDy2
          end if
       end if
    end if

    ! For z direction
    if(nK == 1)then
       ! 1D or 2D
       kR = k; kL = k
       Az = 0.0; Bz = 0.0; Cz = 0.0
    else
       kR = k+1; kL = k-1
       Az = -InvDz2; Bz = 0.0; Cz = +InvDz2
       if(k==1 .and. abs(NeiLbot(iBlock))==1)then
          if(UseFirstOrder)then
             kL = k; Az = -2*InvDz2; Cz = 2*InvDz2
          else
             kL = k+1; kR = k+2; Az = 4*InvDz2; Bz =-3*InvDz2; Cz =-InvDz2
          end if
       elseif(k==nK .and. abs(NeiLtop(iBlock))==1)then
          if(UseFirstOrder)then
             kR = k; Az = -2*InvDz2; Cz = 2*InvDz2
          else
             kL = k-1; kR = k-2; Az = -4*InvDz2; Bz = 3*InvDz2; Cz = InvDz2
          end if
       end if
    end if

    ! Use first-order one-sided difference near the body if needed.
    ! If even first-order fails, then set the current to zero and exit.
    if(.not.true_BLK(iBlock) .and. .not.present(DoIgnoreBody))then
       if(.not.True_Cell(iL,j,k,iBlock).and..not.True_Cell(iR,j,k,iBlock))then
          Current_D = 0.0
          RETURN
       elseif(.not.True_Cell(iL,j,k,iBlock))then
          Ax = 0.0
          if(iR==i+2)then
             Bx =-InvDx2; Cx = InvDx2
          elseif(iR==i-2)then
             Bx = InvDx2; Cx =-InvDx2
          else ! iR==i+1
             Bx =-2.0*InvDx2; Cx = 2.0*InvDx2
          end if
       elseif(.not.True_Cell(iR,j,k,iBlock))then
          Cx = 0.0
          if(iL==i+1)then
             Ax = 2.0*InvDx2; Bx =-2.0*InvDx2
          else ! iL==i-1
             Ax =-2.0*InvDx2; Bx = 2.0*InvDx2
          end if
       end if

       if(nJ > 1)then
          ! 2D or 3D
          if(  .not.True_Cell(i,jL,k,iBlock) .and. &
               .not.True_Cell(i,jR,k,iBlock))then
             Current_D = 0.0
             RETURN
          elseif(.not.True_Cell(i,jL,k,iBlock))then
             Ay = 0.0
             if(jR==j+2)then
                By =-InvDy2; Cy = InvDy2
             elseif(jR==j-2)then
                By = InvDy2; Cy =-InvDy2
             else ! jR==j+1
                By =-2.0*InvDy2; Cy = 2.0*InvDy2
             end if
          elseif(.not.True_Cell(i,jR,k,iBlock))then
             Cy = 0.0
             if(jL==j+1)then
                Ay = 2.0*InvDy2; By =-2.0*InvDy2
             else ! jL==j-1
                Ay =-2.0*InvDy2; By = 2.0*InvDy2
             end if
          end if
       end if

       if(nK > 1)then
          ! 3D
          if(  .not.True_Cell(i,j,kL,iBlock) .and. &
               .not.True_Cell(i,j,kR,iBlock))then
             Current_D = 0.0
             RETURN
          elseif(.not.True_Cell(i,j,kL,iBlock))then
             Az = 0.0
             if(kR==k+2)then
                Bz =-InvDz2; Cz = InvDz2
             elseif(kR==k-2)then
                Bz = InvDz2; Cz =-InvDz2
             else ! kR==k+1
                Bz =-2.0*InvDz2; Cz = 2.0*InvDz2
             end if
          elseif(.not.True_Cell(i,j,kR,iBlock))then
             Cz = 0.0
             if(kL==k+1)then
                Az = 2.0*InvDz2; Bz =-2.0*InvDz2
             else ! kL==k-1
                Az =-2.0*InvDz2; Bz = 2.0*InvDz2
             end if
          end if
       end if
    end if

    if(IsCartesianGrid)then
       call calc_cartesian_j
    else
       call calc_gencoord_j
    end if

    ! Add curl B0 if necessary
    if(UseCurlB0)then
       ! Curl B0 is zero inside rCurrentFreeB0
       if(r_BLK(i,j,k,iBlock) < rCurrentFreeB0) RETURN
       ! Curl B0 is only calculated for the physical cells
       if(i<1 .or. i>nI .or. j<1 .or. j>nJ .or. k<1 .or. k>nK)RETURN
       ! Optimize for multiple calls for the same block
       if(iBlock /= iBlockLast)call set_b0_source(iBlock)
       iBlockLast = iBlock
       ! Add curl(B0) to the current
       Current_D = Current_D + CurlB0_DC(:,i,j,k)
    end if

  contains
    !==========================================================================
    subroutine calc_cartesian_j

      !------------------------------------------------------------------------
      Current_D(x_) = &
           + Ay*State_VGB(Bz_,i,jL,k ,iBlock) &
           + By*State_VGB(Bz_,i,j ,k ,iBlock) &
           + Cy*State_VGB(Bz_,i,jR,k ,iBlock) &
           - Az*State_VGB(By_,i,j ,kL,iBlock) &
           - Bz*State_VGB(By_,i,j ,k ,iBlock) &
           - Cz*State_VGB(By_,i,j ,kR,iBlock)

      Current_D(y_) = &
           + Az*State_VGB(Bx_,i ,j,kL,iBlock) &
           + Bz*State_VGB(Bx_,i ,j,k ,iBlock) &
           + Cz*State_VGB(Bx_,i ,j,kR,iBlock) &
           - Ax*State_VGB(Bz_,iL,j,k ,iBlock) &
           - Bx*State_VGB(Bz_,i ,j,k ,iBlock) &
           - Cx*State_VGB(Bz_,iR,j,k ,iBlock)

      Current_D(z_) = &
           + Ax*State_VGB(By_,iL,j ,k,iBlock) &
           + Bx*State_VGB(By_,i ,j ,k,iBlock) &
           + Cx*State_VGB(By_,iR,j ,k,iBlock) &
           - Ay*State_VGB(Bx_,i ,jL,k,iBlock) &
           - By*State_VGB(Bx_,i ,j ,k,iBlock) &
           - Cy*State_VGB(Bx_,i ,jR,k,iBlock)

      ! Correct current for rz-geometry: Jz = Jz + Bphi/radius
      if(IsRzGeometry) Current_D(x_) = Current_D(x_) &
           + State_VGB(Bz_,i,j,k,iBlock)/Xyz_DGB(y_,i,j,k,iBlock)

    end subroutine calc_cartesian_j

    !==========================================================================

    subroutine calc_gencoord_j

      use ModCoordTransform, ONLY: inverse_matrix

      real :: DxyzDcoord_DD(3,3), DcoordDxyz_DD(3,3), DbDcoord_DD(3,3)
      !------------------------------------------------------------------------

      ! Get the dCartesian/dGencoord matrix with central difference
      DxyzDcoord_DD(:,1) = InvDx2 &
           *(Xyz_DGB(:,i+1,j,k,iBlock) - Xyz_DGB(:,i-1,j,k,iBlock))

      DxyzDcoord_DD(:,2) = InvDy2 &
           *(Xyz_DGB(:,i,j+1,k,iBlock) - Xyz_DGB(:,i,j-1,k,iBlock))

      if(nK > 1)then
         DxyzDcoord_DD(:,3) = InvDz2 &
              *(Xyz_DGB(:,i,j,k+1,iBlock) - Xyz_DGB(:,i,j,k-1,iBlock))
      else
         DxyzDcoord_DD(:,3) = (/ 0., 0., 1./)
      end if

      DcoordDxyz_DD = inverse_matrix(DxyzDcoord_DD, DoIgnoreSingular=.true.)

      ! Calculate the partial derivatives dB/dGencoord
      DbDcoord_DD(:,1) = &
           + Ax*State_VGB(Bx_:Bz_,iL,j,k,iBlock) &
           + Bx*State_VGB(Bx_:Bz_,i ,j,k,iBlock) &
           + Cx*State_VGB(Bx_:Bz_,iR,j,k,iBlock)

      DbDcoord_DD(:,2) = &
           + Ay*State_VGB(Bx_:Bz_,i,jL,k,iBlock) &
           + By*State_VGB(Bx_:Bz_,i,j ,k,iBlock) &
           + Cy*State_VGB(Bx_:Bz_,i,jR,k,iBlock)

      if(nK > 1)then
         DbDcoord_DD(:,3) = &
              + Az*State_VGB(Bx_:Bz_,i,j,kL,iBlock) &
              + Bz*State_VGB(Bx_:Bz_,i,j,k ,iBlock) &
              + Cz*State_VGB(Bx_:Bz_,i,j,kR,iBlock)
      else
         DbDcoord_DD(:,3) = 0.0
      end if

      ! Jx = Dbz/Dy - Dby/Dz = Dbz/Dcoord.Dcoord/Dy - DBy/Dcoord.Dccord/dz
      Current_D(x_) = &
           + sum(DbDcoord_DD(z_,:)*DcoordDxyz_DD(:,y_)) &
           - sum(DbDcoord_DD(y_,:)*DcoordDxyz_DD(:,z_))

      ! Jy = Dbx/Dz - Dbz/Dx
      Current_D(y_) = &
           + sum(DbDcoord_DD(x_,:)*DcoordDxyz_DD(:,z_)) &
           - sum(DbDcoord_DD(z_,:)*DcoordDxyz_DD(:,x_))

      ! Jz = Dby/Dx - Dbx/Dy
      Current_D(z_) = &
           + sum(DbDcoord_DD(y_,:)*DcoordDxyz_DD(:,x_)) &
           - sum(DbDcoord_DD(x_,:)*DcoordDxyz_DD(:,y_))

    end subroutine calc_gencoord_j

  end subroutine get_current

  !===========================================================================

  subroutine calc_field_aligned_current(nTheta, nPhi, rIn, Fac_II, bSm_DII, &
       LatBoundary, Theta_I, Phi_I)

    use ModVarIndexes,     ONLY: Bx_, Bz_, nVar
    use ModMain,           ONLY: Time_Simulation, TypeCoordSystem, nBlock
    use ModPhysics,        ONLY: rCurrents, UnitB_, Si2No_V
    use CON_planet_field,  ONLY: get_planet_field, map_planet_field
    use ModProcMH,         ONLY: iProc, iComm
    use ModNumConst,       ONLY: cHalfPi, cTwoPi, cPi
    use ModMpi

    ! Map the grid points from the rIn radius to rCurrents.
    ! Calculate the field aligned currents there, use the ratio of the
    ! magnetic field strength.
    ! The result is saved into Fac_II, bSm_DII.

    ! Size of spherical grid at rIn
    integer, intent(in) :: nTheta, nPhi

    ! Radius of spherical grid where FAC is needed
    real, intent(in) :: rIn

    ! Field aligned current at rIn
    real, intent(out):: Fac_II(nTheta,nPhi)

    ! Magnetic field at rIn in SM coordinates
    real, intent(out):: bSm_DII(3,nTheta,nPhi)

    ! Lowest latitude that maps up to rCurrents
    real, intent(out), optional :: LatBoundary

    ! Coordinate arrays allow non-uniform grid
    real, intent(in), optional:: Theta_I(nTheta)
    real, intent(in), optional:: Phi_I(nPhi)

    real, allocatable :: bCurrentLocal_VII(:,:,:), bCurrent_VII(:,:,:)

    integer :: i, j, iHemisphere, iError
    real    :: Phi, Theta, Xyz_D(3),XyzIn_D(3),B0_D(3)
    real    :: b_D(3), bRcurrents,Fac, Jr_D(3),bUnit_D(3)
    real    :: bIn_D(3), bIn
    real    :: GmSmg_DD(3,3)
    real    :: State_V(Bx_-1:nVar+3)
    real    :: dPhi, dTheta
    logical :: DoMap
    !-------------------------------------------------------------------------
    if(.not.allocated(bCurrentLocal_VII)) allocate( &
         bCurrentLocal_VII(0:6,nTheta,nPhi), &
         bCurrent_VII(0:6,nTheta,nPhi))

    bCurrentLocal_VII = 0.0
    bCurrent_VII      = 0.0

    Fac_II = 0.0

    GmSmg_DD = transform_matrix(Time_Simulation, 'SMG', TypeCoordSystem)

    if(present(LatBoundary)) LatBoundary = 100.0

    if (abs(rIn-rCurrents)<1.0e-3)then
       DoMap = .false.
    else
       DoMap = .true.
    end if

    dPhi = cTwoPi/(nPhi-1)
    dTheta = cPi /(nTheta-1)

    do j = 1, nPhi
       if(present(Phi_I))then
          Phi = Phi_I(j)
       else
          Phi = (j-1) * dPhi
       end if
       do i = 1, nTheta
          if(present(Theta_I))then
             Theta = Theta_I(i)
          else
             Theta = (i-1) * dTheta
          end if

          call sph_to_xyz(rIn,Theta,Phi, XyzIn_D)

          if (DoMap)then
             ! if not reach to the mapped position, then skip this line
             call map_planet_field(Time_Simulation, XyzIn_D, 'SMG NORM', &
                  rCurrents, Xyz_D, iHemisphere)

             if(iHemisphere == 0) then
                ! Assign weight 1, magnetic field of 1,0,0 and current 0,0,0 
                bCurrentLocal_VII(:,i,j) = &
                     (/1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
                CYCLE
             end if
          else
             Xyz_D = XyzIn_D
          end if

          ! Convert to GM coordinates
          Xyz_D = matmul(GmSmg_DD, Xyz_D)

          if(present(LatBoundary)) &
               LatBoundary = min( abs(Theta - cHalfPi), LatBoundary )

          ! Get the B0 field at the mapped position
          call get_planet_field(Time_Simulation, Xyz_D, &
               TypeCoordSystem//' NORM', B0_D)

          B0_D = B0_D*Si2No_V(UnitB_)

          ! Extract currents and magnetic field for this position
          call get_point_data(0.0, Xyz_D, 1, nBlock, Bx_, nVar+3, State_V)

          bCurrentLocal_VII(0,  i,j) = State_V(Bx_-1)        ! Weight
          bCurrentLocal_VII(1:3,i,j) = State_V(Bx_:Bz_) + &  ! B1 and B0
               State_V(Bx_-1)*B0_D
          bCurrentLocal_VII(4:6,i,j) = State_V(nVar+1:nVar+3) ! Currents

          if(.false. .and. i==6 .and. j==6)then
             write(*,*)'iHemispher=',iHemisphere
             write(*,*)'Phi,Theta=',Phi,Theta
             write(*,*)'XyzIn_D  =', XyzIn_D
             write(*,*)'Xyz_D    =',Xyz_D
             write(*,*)'rCurrents=',rCurrents, sqrt(sum(Xyz_D**2))
             write(*,*)'B0_D     =',B0_D
             write(*,*)'bCurrentLocal_VII =',bCurrentLocal_VII(:,i,j)
             call stop_mpi('DEBUG')
          end if
       end do
    end do

    call MPI_reduce(bCurrentLocal_VII,bCurrent_VII,nTheta*nPhi*7, &
         MPI_REAL,MPI_SUM,0,iComm,iError)

    !\
    ! Map the field aligned current to rIn sphere
    !/

    if(iProc==0)then

       do j = 1, nPhi
          if(present(Phi_I))then
             Phi = Phi_I(j)
          else
             Phi = (j-1) * dPhi
          end if
          do i = 1, nTheta
             if(present(Theta_I))then
                Theta = Theta_I(i)
             else
                Theta = (i-1) * dTheta
             end if

             ! Divide MHD values by the total weight if it exceeds 1.0 
             if(bCurrent_VII(0,i,j) > 1.0) bCurrent_VII(:,i,j) = &
                  bCurrent_VII(:,i,j) / bCurrent_VII(0,i,j)

             ! Extract magnetic field and current
             b_D  = bCurrent_VII(1:3,i,j)
             Jr_D = bCurrent_VII(4:6,i,j)

             ! The strength of the field
             bRcurrents = sqrt(sum(b_D**2))

             ! Convert b_D into a unit vector
             bUnit_D = b_D / bRcurrents
             ! get the field aligned current
             Fac = sum(bUnit_D*Jr_D)

             if(DoMap)then
                ! Calculate magnetic field strength at the rIn grid point
                call sph_to_xyz(rIn, Theta, Phi, XyzIn_D)

                ! Convert to GM coordinates
                XyzIn_D = matmul(GmSmg_DD, XyzIn_D)

                call get_planet_field(Time_Simulation, XyzIn_D, &
                     TypeCoordSystem//' NORM', bIn_D)

                ! Convert to normalized units and get magnitude
                bIn_D = bIn_D*Si2No_V(UnitB_)
                bIn   = sqrt(sum(bIn_D**2))

                ! Multiply by the ratio of the magnetic field strengths
                Fac = bIn / bRcurrents * Fac
             else
                bIn_D = b_D
             end if
             !store the field alinged current
             Fac_II(i,j) = Fac 
             !store the B field in SM coordinates !!

             bSm_DII(:,i,j) = matmul(bIn_D, GmSmg_DD)

          end do
       end do
    end if
    deallocate(bCurrentLocal_VII, bCurrent_VII)

  end subroutine calc_field_aligned_current

end module ModCurrent
