subroutine GM_put_from_pw(Buffer_VI, nVar, nFieldLine, Name_V)

  use CON_coupler, ONLY: PW_, Grid_C
  use ModMain, ONLY: x_, y_, TypeCoordSystem, Time_Simulation
  use ModVarIndexes, ONLY: UseMultiSpecies, SpeciesFirst_, SpeciesLast_
  use ModUtilities, ONLY: lower_case
  use ModPhysics, ONLY: Si2No_V, UnitRho_, UnitRhoU_, rCurrents
  use ModPwGrid
  use ModNumConst, ONLY: cTwoPi
  use ModTriangulate,ONLY:calc_triangulation
  use CON_axes, ONLY: transform_matrix
  use ModCoordTransform, ONLY: dir_to_xyz
  use CON_planet_field,  ONLY: map_planet_field

  implicit none
  character (len=*),parameter :: NameSub='GM_put_from_pw'

  integer, parameter :: Theta_=1, Phi_=2
  
  integer, intent(in)           :: nVar, nFieldLine
  real, intent(in)              :: Buffer_VI(nVar, nFieldLine)
  character (len=*), intent(in) :: Name_V(nVar)


  logical, save :: DoInitialize = .true.

  integer :: iVar, i, iLine, iHemisphere
  real    :: SinThetaOuter, GmPw_DD(3,3), Theta, Phi, XyzPw_D(3)
  character (len=40):: NameVar
  logical :: DoTest, DoTestMe
  !----------------------------------------------------------------------------

  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  ! If only 2 variables are passed, they are the coordinates colatitude 
  ! and longitude needed for the GM->PW pressure coupling.
  if(nVar == 2)then
     if(.not.allocated(CoordXyzPw_DI)) allocate(CoordXyzPw_DI(3, nFieldLine))

     ! Convert from spherical to Cartesian coordinates at radius rCurrents
     ! where the pressure is going to be taken from
     do iLine = 1, nFieldLine
        ! Map field line from ionosphere (taken at r=1) to rCurrents
        Theta = Buffer_VI(Theta_, iLine)
        Phi   = Buffer_VI(Phi_,   iLine)
        call dir_to_xyz(Theta, Phi, XyzPw_D)
        call map_planet_field(Time_Simulation, XyzPw_D, 'SMG NORM', &
             rCurrents,  CoordXyzPw_DI(:,iLine), iHemisphere)
        if(iHemisphere == 0) CoordXyzPw_DI(:,iLine) = 0.0
     end do

     ! Convert from PW to GM coordinates if necessary
     NamePwCoord = Grid_C(PW_) % TypeCoord
     if(TypeCoordSystem /= NamePwCoord) then
        GmPw_DD = &
             transform_matrix(Time_Simulation, NamePwCoord, TypeCoordSystem)
        do iLine = 1, nFieldLine
           CoordXyzPw_DI(:,iLine) = matmul( GmPw_DD, CoordXyzPw_DI(:,iLine))
        end do
     end if

     RETURN
  end if

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

     allocate(CoordXyPw_DI(nCoord, nPoint), iNodeTriangle_II(3, 2*nPoint))

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

  ! Convert to X, Y on a unit sphere (this will be used for triangulation)
  CoordXyPw_DI(x_,1:nLinePw) =  &
         sin(Buffer_VI(Theta_,:)) * cos(Buffer_VI(Phi_,:))
  CoordXyPw_DI(y_,1:nLinePw) =  &
         sin(Buffer_VI(Theta_,:)) * sin(Buffer_VI(Phi_,:))

  ! Total density in normalized units
  StatePw_VI(RhoPw_,1:nLinePw)=sum(Buffer_VI(iRhoPwFirst:iRhoPwLast,:),dim=1) &
       * Si2No_V(UnitRho_)

  ! Field aligned velocity = total moment/total density
  StatePw_VI(Ub_,1:nLinePw)=sum(Buffer_VI(iRhoPwFirst:iRhoPwLast,:) &
       *Buffer_VI(iUPwFirst:iUPwLast,:), dim=1) &
       * Si2No_V(UnitRhou_) / StatePw_VI(RhoPw_,1:nLinePw)

  if(UseMultiSpecies) StatePw_VI(RhoPw_+1: RhoPw_+nSpeciesPw, 1:nLinePw) &
       = Buffer_VI(iRhoPwFirst:iRhoPwLast,:) * Si2No_V(UnitRho_)

  ! Set coordinates for the outer points
  SinThetaOuter = sin(maxval(Buffer_VI(Theta_,:)+dThetaOuter))

  do i = 1, nOuterPoint
     CoordXyPw_DI(x_,nLinePw+i) = SinThetaOuter * cos(i*cTwoPi/nOuterPoint)
     CoordXyPw_DI(y_,nLinePw+i) = SinThetaOuter * sin(i*cTwoPi/nOuterPoint)
  end do

  call calc_triangulation(nPoint, CoordXyPw_DI, iNodeTriangle_II, nTriangle)

  if(DoTestMe)then
     write(*,*)NameSub,': nVarPw, nLinePw, nSpeciesPw=',&
          nVarPw, nLinePw, nSpeciesPw
     write(*,*)NameSub,': nPoint, nTriangle=',nPoint, nTriangle
     write(*,*)NameSub,': CoordXyPw_DI'
     do i=1,nPoint
        write(*,*) i, CoordXyPw_DI(:,i)
     end do
     write(*,*)NameSub,': iNodeTriangle_II'
     do i=1,nTriangle
        write(*,*) i, iNodeTriangle_II(:,i)
     end do
  end if

end subroutine GM_put_from_pw

!==============================================================================

subroutine read_pw_buffer(CoordIn_D, nVarIn, State_V)

  use CON_coupler, ONLY: PW_, Grid_C
  use CON_axes, ONLY: transform_matrix
  use ModMain, ONLY: TypeCoordSystem, Time_Simulation, x_, y_, z_
  use ModVarIndexes, ONLY: Rho_, Ux_, Uz_, &
       SpeciesFirst_, SpeciesLast_, UseMultiSpecies
  use ModPwGrid
  use ModTriangulate, ONLY: find_triangle, triangle_area
  use ModPhysics, ONLY: No2Io_V, UnitU_, UnitRho_
  use CON_planet_field,  ONLY: map_planet_field

  implicit none

  real, intent(in)    :: CoordIn_D(3)
  integer, intent(in) :: nVarIn
  real, intent(inout)   :: State_V(nVarIn)

  real       :: Node_DI(2,3)
  real       :: Triangle1_DI(2,3), Triangle2_DI(2,3), Triangle3_DI(2,3)
  logical    :: IsTriangleFound
  integer    :: iNode1,iNode2,iNode3
  real       :: Area1,Area2,Area3,Area

  ! A short time relative to the rotation period of Earth
  real, parameter :: dTimeMax = 600.0 ! [s]
  real    :: B0_D(3), XyzPw_D(3), Xyz_D(3), Xy_D(2)
  integer :: iPoint, iHemisphere
  logical :: DoInitialize = .true.

  real    :: TimeSimLast = -1.0 - dTimeMax
  real, save:: PwGm_DD(3,3)

  character (len=*), parameter :: NameSub = 'read_pw_buffer'
  
  logical :: DoTest, DoTestMe
  !--------------------------------------------------------

  if(DoInitialize)then
     DoInitialize = .false.
     ! Fill in outer points with body values coming in via State_V
     do iPoint=nLinePw+1, nPoint
        StatePw_VI(Ub_, iPoint) = 0.0
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
        ! Update GM to PW transformation if necessary
        PwGm_DD = &
             transform_matrix(Time_Simulation, TypeCoordSystem, NamePwCoord)
        TimeSimLast = Time_Simulation
     end if
     ! Convert from GM coordinates to buffer coordinates
     XyzPw_D = matmul( PwGm_DD, CoordIn_D)
  end if

  ! Map down to radius 1.0 where the PW fieldline positions are defined
  call map_planet_field(Time_Simulation, XyzPw_D, 'SMG NORM', &
       1.0, Xyz_D, iHemisphere)

  ! Project to unit sphere and calculate X, Y coordinates
  ! Disregard the 3rd coordinate (Z) that determines the hemisphere !!!
  Xy_D(1:2) = Xyz_D(1:2)

  !Find triangle containing point Xy_D

  call find_triangle(nPoint, nTriangle, Xy_D, CoordXyPw_DI, iNodeTriangle_II, &
       iNode1, iNode2, iNode3, IsTriangleFound)

  ! Point is not covered: leave the input state variables alone
  if (.not.IsTriangleFound) RETURN

  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  ! interpolate values

  Node_DI(:,1)=CoordXyPw_DI(:,iNode1)
  Node_DI(:,2)=CoordXyPw_DI(:,iNode2)
  Node_DI(:,3)=CoordXyPw_DI(:,iNode3)

  Triangle1_DI(:,1)     = Xy_D(:)
  Triangle1_DI(:,2:3)   = Node_DI(:,2:3)

  Triangle2_DI(:,1)     = Xy_D(:)
  Triangle2_DI(:,2:3)   = Node_DI(:,1:3:2)

  Triangle3_DI(:,1)     = Xy_D(:)
  Triangle3_DI(:,2:3)   = Node_DI(:,1:2)

  Area1 = triangle_area(Triangle1_DI)
  Area2 = triangle_area(Triangle2_DI)
  Area3 = triangle_area(Triangle3_DI)

  ! Normalize areas to weights
  Area  = Area1+Area2+Area3
  Area1 = Area1 / Area
  Area2 = Area2 / Area
  Area3 = Area3 / Area

  State_V(Rho_) = &
       Area1*StatePw_VI(RhoPw_,iNode1) + &
       Area2*StatePw_VI(RhoPw_,iNode2) + &
       Area3*StatePw_VI(RhoPw_,iNode3)

  ! Calculate field aligned momentum vector
  call get_b0(CoordIn_D(x_), CoordIn_D(y_), CoordIn_D(z_), B0_D)
  B0_D = B0_D / sqrt(sum(B0_D**2))
  ! Make sure unit vector is pointing outward
  if(sum(B0_D*CoordIn_D) < 0.)B0_D = -B0_D
  State_V(Ux_:Uz_) = B0_D * ( &
       Area1*StatePw_VI(Ub_,iNode1) + &
       Area2*StatePw_VI(Ub_,iNode2) + &
       Area3*StatePw_VI(Ub_,iNode3))

  if(UseMultiSpecies) &
       State_V(SpeciesFirst_:SpeciesLast_) = &
       Area1*StatePw_VI(iRhoPwFirst:iRhoPwLast,iNode1) + &
       Area2*StatePw_VI(iRhoPwFirst:iRhoPwLast,iNode2) + &
       Area3*StatePw_VI(iRhoPwFirst:iRhoPwLast,iNode3)

  if(DoTestMe)then
     write(*,*)NameSub,' finished with'
     write(*,*)'CoordIn_D   =',CoordIn_D
     write(*,*)'XyzPw_D     =',XyzPw_D
     write(*,*)'Xy_D        =',Xy_D
     write(*,*)'Node_DI(:,1)=',Node_DI(:,1)
     write(*,*)'Node_DI(:,2)=',Node_DI(:,2)
     write(*,*)'Node_DI(:,3)=',Node_DI(:,3)
     write(*,*)'Area1,2,3=',Area1,Area2,Area3
     write(*,*)'State_V(Rho_)   =',State_V(Rho_)*No2Io_V(UnitRho_)
     write(*,*)'B0_D            =',B0_D
     write(*,*)'State_V(Ux_:Uz_)=',State_V(Ux_:Uz_)*No2Io_V(UnitU_)
  end if

end subroutine read_pw_buffer
