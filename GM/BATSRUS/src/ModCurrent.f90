!^CFG COPYRIGHT UM      
module ModCurrent

  implicit none

contains
  !============================================================================  
  subroutine calc_field_aligned_current(nTheta, nPhi, rIn, Fac_II, b_DII)

    use ModVarIndexes,     ONLY: Bx_, Bz_, nVar
    use ModMain,           ONLY: Time_Simulation, TypeCoordSystem, nBlock
    use ModPhysics,        ONLY: rCurrents,UnitB_,Si2No_V, UnitJ_,No2Si_V,No2Io_V
    use CON_planet_field,  ONLY: get_planet_field, map_planet_field
    use ModProcMH,         ONLY: iProc, iComm
    use ModNumConst, ONLY: cHalfPi, cTwoPi, cPi
    use ModMpi
    use ModCoordTransform, ONLY: sph_to_xyz
    use CON_axes,          ONLY: transform_matrix
    

    ! Map the grid points from the rIn radius to rCurrents.
    ! Calculate the field aligned currents there, use the ratio of the
    ! magnetic field strength, and project to the radial direction.
    ! The result is saved into Fac_II, b_DII.

    integer, intent(in) :: nTheta, nPhi
    real, intent(in) :: rIn
    real, intent(out), dimension(nTheta,nPhi):: Fac_II
    real, intent(out), dimension(3,nTheta,nPhi):: b_DII
    real, allocatable :: bCurrentLocal_VII(:,:,:), bCurrent_VII(:,:,:)
    integer :: i, j, iHemisphere, iError
    real    :: Phi, Theta, Xyz_D(3),XyzIn_D(3),B0_D(3)
    real    :: b_D(3), b, bRcurrents,Fac, Jr_D(3),bUnit_D(3)
    real    :: GmSmg_DD(3,3),GsmtoSmg_DD(3,3)
    real    :: State_V(Bx_-1:nVar+3)
    real    :: LatBoundary, dPhi, dTheta
    logical :: DoMap
    !-------------------------------------------------------------------------

    if(.not.allocated(bCurrentLocal_VII))then
       allocate(bCurrentLocal_VII(0:6,nTheta,nPhi))
       allocate(bCurrent_VII(0:6,nTheta,nPhi))
    end if

    bCurrentLocal_VII = 0.0
    bCurrent_VII = 0.0

    GmSmg_DD = transform_matrix(Time_Simulation, 'SMG', TypeCoordSystem)
!    GsmtoSmg_DD = transform_matrix(TIme_Simulation,TypeCoordSystem,'SMG')

    LatBoundary = 100.0

    if (abs(rIn-rCurrents)<1.0e-3)then
       DoMap = .false.
    else
       DoMap = .true.
    end if

    dPhi = cTwoPi/(nPhi-1)
    dTheta = cPi /(nTheta-1)

    do j = 1, nPhi
       Phi = (j-1) * dPhi
       do i = 1, nTheta
          Theta = (i-1) * dTheta
          if (i==1)Theta = 1.0E-4

          call sph_to_xyz(rIn,Theta,Phi, XyzIn_D)
       
          if (DoMap)then
             ! if not reach to the mapped position, then skip this line
             call map_planet_field(Time_Simulation, XyzIn_D, 'SMG NORM', &
                  rCurrents, Xyz_D, iHemisphere)
             
             if(iHemisphere == 0) then
             ! Assign weight 1, magnetic field of 1,0,0 and current 0,0,0 
                bCurrentLocal_VII(:,i,j) = (/1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
                CYCLE
             end if
          else
             Xyz_D = XyzIn_D
          end if

          LatBoundary = min( abs(Theta - cHalfPi), LatBoundary )
          
          ! Get the B0 field at the mapped position                               
          call get_planet_field(Time_Simulation, Xyz_D,'SMG NORM',B0_D)
          B0_D = B0_D*Si2No_V(UnitB_)            

          ! Convert to GM coordinates
          Xyz_D = matmul(GmSmg_DD, Xyz_D)

          ! Extract currents and magnetic field for this position
          call get_point_data(0.0, Xyz_D, 1, nBlock, Bx_, nVar+3, State_V)

          bCurrentLocal_VII(0,  i,j) = State_V(Bx_-1)        ! Weight
          bCurrentLocal_VII(1:3,i,j) = State_V(Bx_:Bz_) + &  ! B1 and B0
               State_V(Bx_-1)*B0_D
          bCurrentLocal_VII(4:6,i,j) = State_V(nVar+1:nVar+3) ! Currents
          ! convert the vectors to SMG coordinates
!          bCurrentLocal_VII(1:3,i,j) = matmul(GsmtoSmg_DD,bCurrentLocal_VII(1:3,i,j))
           bCurrentLocal_VII(1:3,i,j) = matmul(bCurrentLocal_VII(1:3,i,j),GmSmg_DD)

          if(.false. .and. i==6 .and. j==6)then
             write(*,*)'iHemispher=',iHemisphere
             write(*,*)'Phi,Theta=',Phi,Theta
             write(*,*)'Xyz_D    =',Xyz_D
             write(*,*)'rCurrents=',rCurrents
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
          Phi = (j-1) * dPhi
          do i = 1, nTheta
            Theta = (i-1) * dTheta
             if(i==1)Theta = 1.0E-4
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
                ! Calculate magnetic field strength at the rIn grid points
                call sph_to_xyz(rIn,Theta,Phi, XyzIn_D)
                call get_planet_field(Time_Simulation, XyzIn_D,'SMG NORM',b_D)
                
                ! Convert to GM units and get magnitude
                b_D = b_D*Si2No_V(UnitB_)
                b   = sqrt(sum(b_D**2))
                
                ! Multiply by the ratio of the magnetic field strengths
                Fac = b / bRcurrents * Fac

             end if
            !store the field alinged current
             Fac_II(i,j) = Fac 
             !store the B field
             b_DII(:,i,j) = b_D
             
          end do
       end do
       ! Save the latitude boundary information to the equator.                     
       ! The reason to use Si2No_V(UnitJ_) is in order to be converted back
       ! in the coupler (ie.,IE). For other purpose (not for the coupler), 
       ! Please avoid this point !
       Fac_II(nTheta/2:nTheta/2+1,1) = LatBoundary * Si2No_V(UnitJ_)

    end if
    deallocate(bCurrentLocal_VII,bCurrent_VII)
    
  end subroutine calc_field_aligned_current

end module ModCurrent
