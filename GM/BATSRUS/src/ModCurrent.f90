!This code is a copyright protected software (c) 2002- University of Michigan      
module ModCurrent

  use ModCoordTransform, ONLY: sph_to_xyz
  use CON_axes,          ONLY: transform_matrix

  implicit none

contains
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
    use ModGeometry, ONLY: True_Cell, true_BLK
    use BATL_lib, ONLY: IsCartesianGrid, IsRzGeometry, Xyz_DGB, CellSize_DB
    use ModParallel, ONLY: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot
    use ModSize,     ONLY: nI, nJ, nK, x_, y_, z_

    integer, intent(in) :: i, j, k, iBlock
    real,    intent(out):: Current_D(3)

    integer, optional, intent(in):: nOrderResChange

    logical, optional, intent(in):: DoIgnoreBody

    logical:: UseFirstOrder
    integer:: iL, iR, jL, jR, kL, kR
    real   :: Ax, Bx, Cx, Ay, By, Cy, Az, Bz, Cz
    real   :: InvDx2, InvDy2, InvDz2
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
          if(.not.True_Cell(i,jL,k,iBlock).and..not.True_Cell(i,jR,k,iBlock))then
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
          if(.not.True_Cell(i,j,kL,iBlock).and..not.True_Cell(i,j,kR,iBlock))then
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

      DxyzDcoord_DD(:,3) = InvDz2 &
           *(Xyz_DGB(:,i,j,k+1,iBlock) - Xyz_DGB(:,i,j,k-1,iBlock))

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

      DbDcoord_DD(:,3) = &
           + Az*State_VGB(Bx_:Bz_,i,j,kL,iBlock) &
           + Bz*State_VGB(Bx_:Bz_,i,j,k ,iBlock) &
           + Cz*State_VGB(Bx_:Bz_,i,j,kR,iBlock)

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
    use ModMain,           ONLY: Time_Simulation, TypeCoordSystem, nBlock, &
         Test_String
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

    logical :: UseWrongB = .false.

    !-------------------------------------------------------------------------
    UseWrongB = index(Test_String,'FACBUG') > 0

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

          if(UseWrongB) B0_D = matmul(B0_D, GmSmg_DD)

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

             if(UseWrongB)then
                bSm_DII(:,i,j) = b_D
             else
                bSm_DII(:,i,j) = matmul(bIn_D, GmSmg_DD)
             end if

          end do
       end do
    end if
    deallocate(bCurrentLocal_VII, bCurrent_VII)

  end subroutine calc_field_aligned_current

end module ModCurrent
