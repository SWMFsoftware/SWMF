subroutine GM_put_from_pw(Buffer_VI, nVar, nFieldLine, Name_V)

  use CON_coupler, ONLY: PW_, Grid_C
  use ModMain, ONLY: x_, y_, TypeCoordSystem, Time_Simulation
  use ModVarIndexes, ONLY: UseMultiSpecies, SpeciesFirst_, SpeciesLast_
  use ModMultiFluid, ONLY: UseMultiIon, nIonFluid
  use ModPhysics, ONLY: Si2No_V, UnitRho_, UnitU_, rCurrents
  use ModPwGrid
  use ModPhysics, ONLY: rBody
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

  integer :: i, iLine, iHemisphere
  real    :: SinThetaOuter, GmPw_DD(3,3), Theta, Phi, XyzPw_D(3), Factor

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

     ! Subtract 2 coords and half of the remaining variables are velocities
     nSpeciesPw = (nVar-2)/2
     ! First 2 variables are coordinates, 
     ! then nSpeciesPw densities then nSpeciesPw velocities
     iRhoPwFirst = 3
     iRhoPwLast  = 2 + nSpeciesPw
     iUPwFirst   = iRhoPwLast + 1
     iUPwLast    = nVar

     if(UseMultiSpecies)then
        ! Species densities and total velocity
        iRhoGmFirst = 1
        iRhoGmLast  = SpeciesLast_-SpeciesFirst_+1
        iUGmFirst   = iRhoGmLast + 1
        iUGmLast    = iUGmFirst

        iRhoPwLast = min(iRhoPwLast, iRhoPwFirst + SpeciesLast_ - SpeciesFirst_)
     elseif(UseMultiIon)then
        ! Fluid densities and velocities
        iRhoGmFirst = 1
        iRhoGmLast  = nIonFluid
        iUGmFirst   = nIonFluid + 1
        iUGmLast    = 2*nIonFluid
        
        iRhoPwLast = min(iRhoPwLast, iRhoPwFirst + nIonFluid - 1)
        iUPwLast   = min(iUPwLast,   iUPwFirst   + nIonFluid - 1)
     else
        ! Total density and velocity
        iRhoGmFirst = 1
        iRhoGmLast  = 1
        iUGmFirst   = 2
        iUGmLast    = 2
     end if
     allocate(StateGm_VI(1:iUGmLast, nPoint))
     allocate(CoordXyPw_DI(nCoord, nPoint), iNodeTriangle_II(3, 2*nPoint))

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

  if(UseMultiIon)then
     StateGm_VI(iRhoGmFirst:iRhoGmLast, 1:nLinePw) &
          = Buffer_VI(iRhoPwFirst:iRhoPwLast,:)*Si2No_V(UnitRho_)

     StateGm_VI(iUGmFirst:iUGmLast, 1:nLinePw) &
          = Buffer_VI(iUPwFirst:iUPwLast,:) * Si2No_V(UnitU_)
  else
     ! Set total 
     if(UseMultiSpecies)then
        StateGm_VI(iRhoGmFirst:iRhoGmLast, 1:nLinePw) &
             = Buffer_VI(iRhoPwFirst:iRhoPwLast,:) * Si2No_V(UnitRho_)
     else
        ! Total density in normalized units
        StateGm_VI(iRhoGmFirst, 1:nLinePw) &
             = sum(Buffer_VI(iRhoPwFirst:iRhoPwLast,:),dim=1) * Si2No_V(UnitRho_)
     end if

     ! Field aligned velocity = total moment/total density
     StateGm_VI(iUGmFirst,1:nLinePw) &
          = sum(Buffer_VI(iRhoPwFirst:iRhoPwLast,:) &
          *     Buffer_VI(iUPwFirst:iUPwLast,:), dim=1) &
          / sum(Buffer_VI(iRhoPwFirst:iRhoPwLast,:), dim=1) &
          * Si2No_V(UnitU_) 
  end if

  ! Reduce densities according to a 1/r3 law based on magnetic field strength
  ! and assuming that PWOM solves up to R=2.25
  !Factor = (2.25/rBody)**3
  !StateGm_VI(iRhoGmFirst:iRhoGmLast,1:nLinePw) = &
  !     StateGm_VI(iRhoGmFirst:iRhoGmLast,1:nLinePw)*Factor

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
  use ModVarIndexes, ONLY: Rho_, Ux_, Uy_, Uz_, &
       SpeciesFirst_, SpeciesLast_, UseMultiSpecies, IsMhd
  use ModMultiFluid, ONLY: UseMultiIon, nIonFluid, &
       iRhoIon_I, iUxIon_I, iUyIon_I, iUzIon_I
      
  use ModPwGrid
  use ModTriangulate, ONLY: find_triangle
  use ModPhysics, ONLY: No2Io_V, UnitU_, UnitRho_
  use CON_planet_field,  ONLY: map_planet_field

  implicit none

  real, intent(in)    :: CoordIn_D(3)
  integer, intent(in) :: nVarIn
  real, intent(inout) :: State_V(nVarIn)

  logical    :: IsTriangleFound
  integer    :: iNode1, iNode2, iNode3, iIon, iUGm
  real       :: Area1, Area2, Area3

  ! A short time relative to the rotation period of Earth
  real, parameter :: dTimeMax = 600.0 ! [s]
  real    :: B0_D(3), XyzPw_D(3), Xyz_D(3), Xy_D(2)
  integer :: iPoint, iHemisphere
  logical :: DoInitialize = .true.

  real    :: TimeSimLast = -1.0 - dTimeMax
  real, save:: PwGm_DD(3,3)

  character (len=*), parameter :: NameSub = 'read_pw_buffer'
  
  logical, parameter :: DoTestMe = .false.
  !--------------------------------------------------------

  if(DoInitialize)then
     DoInitialize = .false.
     ! Fill in outer points with body values coming in via State_V
     do iPoint=nLinePw+1, nPoint
        StateGm_VI(iUGmFirst:iUGmLast, iPoint) = 0.0
        if(UseMultiSpecies)then
           StateGm_VI(iRhoGmFirst:iRhoGmLast, iPoint) = &
                State_V(SpeciesFirst_:SpeciesLast_)
        elseif(UseMultiIon)then
           StateGm_VI(iRhoGmFirst:iRhoGmLast, iPoint) = State_V(iRhoIon_I)
        else
           StateGm_VI(iRhoGmFirst, iPoint) = State_V(Rho_)
        end if
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

  call find_triangle(&
       nPoint, nTriangle, Xy_D, CoordXyPw_DI, &
       iNodeTriangle_II(:,1:nTriangle), &
       iNode1, iNode2, iNode3, Area1, Area2, Area3, IsTriangleFound)

  ! Point is not covered: leave the input state variables alone
  if (.not.IsTriangleFound) RETURN

  ! Calculate field aligned momentum vector
  call get_b0(CoordIn_D(x_), CoordIn_D(y_), CoordIn_D(z_), B0_D)
  B0_D = B0_D / sqrt(sum(B0_D**2))

  ! Make sure unit vector is pointing outward
  if(sum(B0_D*CoordIn_D) < 0.)B0_D = -B0_D

  ! interpolate values
  if(UseMultiIon)then
     State_V(iRhoIon_I) = &
          Area1*StateGm_VI(iRhoGmFirst:iRhoGmLast,iNode1) + &
          Area2*StateGm_VI(iRhoGmFirst:iRhoGmLast,iNode2) + &
          Area3*StateGm_VI(iRhoGmFirst:iRhoGmLast,iNode3)
     do iIon = 1, nIonFluid
        ! Calculate field aligned velocity vector per fluid
        iUGm = iUGmFirst + iIon - 1
        State_V(iUxIon_I(iIon):iUzIon_I(iIon)) = B0_D * &
             ( Area1*StateGm_VI(iUGm, iNode1) &
             + Area2*StateGm_VI(iUGm, iNode2) &
             + Area3*StateGm_VI(iUGm, iNode3))
     end do
     if(IsMhd)then
        ! Calculate total density and average velocity
        State_V(Rho_)= sum(State_V(iRhoIon_I))
        State_V(Ux_) = sum(State_V(iRhoIon_I)*State_V(iUxIon_I))/State_V(Rho_)
        State_V(Uy_) = sum(State_V(iRhoIon_I)*State_V(iUyIon_I))/State_V(Rho_)
        State_V(Uz_) = sum(State_V(iRhoIon_I)*State_V(iUzIon_I))/State_V(Rho_)
     end if
  else
     if(UseMultiSpecies) then
        ! Interpolate species densities and put total into Rho_
        State_V(SpeciesFirst_:SpeciesLast_) = &
             Area1*StateGm_VI(iRhoGmFirst:iRhoGmLast,iNode1) + &
             Area2*StateGm_VI(iRhoGmFirst:iRhoGmLast,iNode2) + &
             Area3*StateGm_VI(iRhoGmFirst:iRhoGmLast,iNode3)
        State_V(Rho_) = sum(State_V(SpeciesFirst_:SpeciesLast_))
     else
        ! For simple MHD the StateGm_VI contains the total density
        State_V(Rho_) = &
             Area1*StateGm_VI(iRhoGmFirst, iNode1) + &
             Area2*StateGm_VI(iRhoGmFirst, iNode2) + &
             Area3*StateGm_VI(iRhoGmFirst, iNode3)
     end if

     ! Field aligned velocity
     State_V(Ux_:Uz_) = B0_D * ( &
          Area1*StateGm_VI(iUGmFirst, iNode1) + &
          Area2*StateGm_VI(iUGmFirst, iNode2) + &
          Area3*StateGm_VI(iUGmFirst, iNode3))

  end if
  if(DoTestMe)then
     write(*,*)NameSub,' finished with'
     write(*,*)'CoordIn_D     =', CoordIn_D
     write(*,*)'XyzPw_D       =', XyzPw_D
     write(*,*)'Xy_D          =', Xy_D
     write(*,*)'Area1,2,3     =', Area1, Area2, Area3
     write(*,*)'State_V(Rho_) =', State_V(Rho_)*No2Io_V(UnitRho_)
     write(*,*)'B0_D          =', B0_D
     write(*,*)'State_V(Ux:Uz)=', State_V(Ux_:Uz_)*No2Io_V(UnitU_)
  end if

end subroutine read_pw_buffer
