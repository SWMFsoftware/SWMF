subroutine GM_put_from_pw(Buffer_VI, nFieldLine, nVar, Name_V)

  use ModMain, ONLY: x_, y_
  use ModVarIndexes, ONLY: UseMultiSpecies, SpeciesFirst_, SpeciesLast_
  use ModUtilities, ONLY: lower_case
  use ModPhysics, ONLY: UnitSi_Rho, UnitSi_RhoU
  use ModPwGrid
  use ModNumConst, ONLY: cTwoPi
  use ModTriangulate,ONLY:calc_triangulation

  implicit none
  character (len=*),parameter :: NameSub='GM_put_from_pw'

  integer, parameter :: Theta_=1, Phi_=2
  
  integer, intent(in)           :: nVar, nFieldLine
  real, intent(in)              :: Buffer_VI(nVar, nFieldLine)
  character (len=*), intent(in) :: Name_V(nVar)


  logical, save :: DoInitialize = .true.

  integer :: iVar, i, nTriangle
  real    :: SinThetaOuter
  character (len=40):: NameVar
  logical :: DoTest, DoTestMe
  !----------------------------------------------------------------------------

  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  if(DoInitialize)then
     DoInitialize=.false.

     nVarPw  = nVar
     nLinePw = nFieldLine
     nPoint = nLinePw + nOuterPoint

     do iVar = 1, nVar
        NameVar = Name_V(iVar)
        call lower_case(NameVar)
        if(iRhoPwFirst == -1 .and. NameVar(1:1) == 'd') iRhoPwFirst = iVar
        if(NameVar(1:1) == 'd')                         iRhoPwLast  = iVar
        if(iUPwFirst == -1 .and. NameVar(1:1) == 'v')   iUPwFirst   = iVar
        if(NameVar(1:1) == 'v')                         iUPwLast    = iVar
     end do

     nSpeciesPw = iRhoPwLast - iRhoPwFirst + 1

     if(nSpeciesPw /= iUPwLast - iUPwFirst + 1) then
        write(*,*)NameSub,'iRhoPwFirst,iRhoPwLast,iUPwFirst,iUPwLast=',&
             iRhoPwFirst,iRhoPwLast,iUPwFirst,iUPwLast
        call stop_mpi(NameSub// &
             ' ERROR: iRhoPwLast - iRhoPwFirst /= iUPwLast - iUPwFirst !')
     end if

     allocate(CoordXyPw_DI(nCoord, nPoint), iNodeTriangle_II(3, nPoint))

     if(UseMultiSpecies)then

        if(SpeciesLast_-SpeciesFirst_+1 /= nSpeciesPw) then
           write(*,*)NameSub,' SpeciesFirst_, SpeciesLast_, nSpeciesPw=',&
                SpeciesFirst_, SpeciesLast_, nSpeciesPw
           call stop_mpi(NameSub// &
                ' ERROR: nSpeciesPw /= nSpecies in ModEquation')
        end if
        ! Total momentum, density and species densities
        allocate(StatePw_VI(nSpeciesPw + 2, nPoint))
     else
        ! Total momentum and density
        allocate(StatePw_VI(2, nPoint))
     end if
  end if

  if(nLinePw /= nFieldLine .or. nVarPw /= nVar)then
     write(*,*)NameSub,' ERROR: nLinePw=',nLinePw,'/= nFieldLine=',nFieldLine,&
          ' or nVarPw=',nVarPw,' /= nVar=',nVar
     call CON_stop(NameSub,' nFieldLine or nVar has changed')
  end if

  ! Convert to X, Y on a unit sphere
  CoordXyPw_DI(x_,1:nLinePw) =  &
         sin(Buffer_VI(Theta_,:)) * cos(Buffer_VI(Phi_,:))
  CoordXyPw_DI(y_,1:nLinePw) =  &
         sin(Buffer_VI(Theta_,:)) * sin(Buffer_VI(Phi_,:))

  StatePw_VI(RhoUb_,1:nLinePw)=sum(Buffer_VI(iRhoPwFirst:iRhoPwLast,:) &
       *Buffer_VI(iUPwFirst:iUPwLast,:), dim=1) &
       / UnitSI_rhoU

  StatePw_VI(RhoPw_,1:nLinePw)=sum(Buffer_VI(iRhoPwFirst:iRhoPwLast,:),dim=1) &
       / UnitSI_rho

  if(UseMultiSpecies) StatePw_VI(RhoPw_+1: RhoPw_+nSpeciesPw, 1:nLinePw) &
       = Buffer_VI(iRhoPwFirst:iRhoPwLast,:) / UnitSI_rho

  ! Set coordinates for the outer points
  SinThetaOuter = sin(maxval(Buffer_VI(Theta_,:)+dThetaOuter))

  do i = 1, nOuterPoint
     CoordXyPw_DI(x_,nLinePw+i) = SinThetaOuter * cos(i*cTwoPi/nOuterPoint)
     CoordXyPw_DI(y_,nLinePw+i) = SinThetaOuter * sin(i*cTwoPi/nOuterPoint)
  end do


  call calc_triangulation(nPoint, CoordXyPw_DI, iNodeTriangle_II,nTriangle)

  if(DoTestMe)then
     write(*,*)'!!! nVarPw, nLinePw, nSpeciesPw=',nVarPw, nLinePw, nSpeciesPw
     write(*,*)'!!! Buffer_VI   =',Buffer_VI
     write(*,*)'!!! CoordXyPw_DI=',CoordXyPw_DI
     write(*,*)'!!! StatePw_VI  =',StatePw_VI
  end if

end subroutine GM_put_from_pw

!==============================================================================

subroutine read_pw_buffer(CoordIn_D, nVarIn, State_V)

  use CON_coupler, ONLY: PW_, Grid_C
  use CON_axes, ONLY: transform_matrix
  use ModMain, ONLY: TypeCoordSystem, Time_Simulation, x_, y_, z_
  use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUz_, &
       SpeciesFirst_, SpeciesLast_, UseMultiSpecies
  use ModPwGrid
  use ModTriangulate, ONLY: find_triangle, triangle_area

  real, intent(in)    :: CoordIn_D(3)
  integer, intent(in) :: nVarIn
  real, intent(inout)   :: State_V(nVarIn)

  real       :: Node_DI(2,3)
  real       :: Triangle1_DI(2,3), Triangle2_DI(2,3), Triangle3_DI(2,3)
  logical    :: IsTriangleFound
  integer    :: node1,node2,node3
  real       :: Area1,Area2,Area3,Area

  ! A short time relative to the rotation period of Earth
  real, parameter :: dTimeMax = 600.0 ! [s]
  real    :: TimeSimLast = -1.0
  real    :: B0_D(3), XyzPw_D(3), Xy_D(2), PwGm_DD(3,3)
  integer :: iPoint
  logical :: DoInitialize = .true.
  character (len=3) :: NamePwCoord = '???'

  character (len=*), parameter :: NameSub = 'read_pw_buffer'
  !--------------------------------------------------------

  if(DoInitialize)then
     DoInitialize = .false.
     ! Fill in outer points with body values coming in via State_V
     do iPoint=nLinePw+1, nPoint
        StatePw_VI(RhoUb_, iPoint) = 0.0
        StatePw_VI(RhoPw_, iPoint) = State_V(Rho_)
        if(UseMultiSpecies) StatePw_VI(RhoPw_+1:nVarPw, iPoint) = &
             State_V(SpeciesFirst_:SpeciesLast_)
     end do

     NamePwCoord = Grid_C(PW_) % TypeCoord
  end if

  if(TypeCoordSystem == NamePwCoord)then
     XyzPw_D = CoordIn_D
  else
     if(abs(Time_Simulation - TimeSimLast) > dTimeMax)then
        ! Update GM-PW transformation if necessary
        PwGm_DD = transform_matrix(Time_Simulation, &
             NamePwCoord, TypeCoordSystem)
        TimeSimLast = Time_Simulation
     end if
     ! Convert from GM coordinates to buffer coordinates
     XyzPw_D    = matmul( PwGm_DD, CoordIn_D)
  end if

  ! Project to unit sphere and calculate X, Y coordinates
  ! Disregard the 3rd coordinate (Z) that determines the hemisphere !!!
  Xy_D(1:2) = XyzPw_D(1:2)/sqrt(sum(XyzPw_D**2))


  !Find triangle containing point
  !call find_triangle(nPoint, CoordXyPw_DI, iNodeTriangle_II, Xy_D, &
  !     node1, node2, node3, IsTriangleFound)

  Node_DI(:,1)=CoordXyPw_DI(:,node1)
  Node_DI(:,2)=CoordXyPw_DI(:,node2)
  Node_DI(:,3)=CoordXyPw_DI(:,node3)

  ! interpolate values
  if (IsTriangleFound) then
     Triangle1_DI(:,1)     = Xy_D(:)
     Triangle1_DI(:,2:3)   = Node_DI(:,2:3)

     Triangle2_DI(:,1)     = Xy_D(:)
     Triangle2_DI(:,1:3:2) = Node_DI(:,1:3:2)

     Triangle3_DI(:,1)     = Xy_D(:)
     Triangle3_DI(:,2:3)   = Node_DI(:,2:3)

     Area1 = triangle_area(Triangle1_DI)
     Area2 = triangle_area(Triangle2_DI)
     Area3 = triangle_area(Triangle3_DI)

     ! Normalize areas to weights
     Area  = Area1+Area2+Area3
     Area1 = Area1 / Area
     Area2 = Area2 / Area
     Area3 = Area3 / Area

     State_V(Rho_) = &
          Area1*StatePw_VI(RhoPw_,node1) + &
          Area2*StatePw_VI(RhoPw_,node2) + &
          Area3*StatePw_VI(RhoPw_,node3)

     ! Calculate field aligned momentum vector
     call get_b0(CoordIn_D(x_), CoordIn_D(y_), CoordIn_D(z_), B0_D)
     B0_D = B0_D / sqrt(sum(B0_D**2))
     State_V(RhoUx_:RhoUz_) = B0_D * ( &
          Area1*StatePw_VI(RhoUb_,node1) + &
          Area2*StatePw_VI(RhoUb_,node2) + &
          Area3*StatePw_VI(RhoUb_,node3))

     if(UseMultiSpecies) &
          State_V(SpeciesFirst_:SpeciesLast_) = &
          Area1*StatePw_VI(iRhoPwFirst:iRhoPwLast,node1) + &
          Area2*StatePw_VI(iRhoPwFirst:iRhoPwLast,node2) + &
          Area3*StatePw_VI(iRhoPwFirst:iRhoPwLast,node3)

  else
     write(*,*)NameSub, 'CoordIn_D:', CoordIn_D
     call stop_mpi(NameSub,' ERROR: point is outside of all triangles!')
  end if

end subroutine read_pw_buffer
