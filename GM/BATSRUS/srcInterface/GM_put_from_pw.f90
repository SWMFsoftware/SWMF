subroutine GM_put_from_pw(Buffer_VI, nVar, nFieldLine, Name_V)

  use CON_coupler, ONLY: PW_, Grid_C
  use ModMain, ONLY: x_, y_, TypeCoordSystem, Time_Simulation
  use ModVarIndexes, ONLY: UseMultiSpecies, SpeciesFirst_, SpeciesLast_
  use ModMultiFluid, ONLY: UseMultiIon, nIonFluid
  use ModPhysics, ONLY: Si2No_V, UnitRho_, UnitU_, rCurrents
  use ModPwGrid
  use ModPhysics, ONLY: rBody
  use ModNumConst, ONLY: cTwoPi,cHalfPi
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
  real    :: SinThetaOuter, GmPw_DD(3,3), Theta, Phi, XyzPw_D(3), Factor, &
       SinThetaOuter1,SinThetaOuter2,tmp1_array(nFieldLine),tmp2_array(nFieldLine)
  integer :: Tmp_array(nFieldLine)
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

     Tmp_array = 0
     where(Buffer_VI(Theta_,:)<cHalfPi)
        Tmp_array = 1
     end where
     nLinePw1 = sum(Tmp_array)
     nPoint1 = nLinePw1 + nOuterPoint    

     nLinePw2 = nLinePw - nLinePw1
     nPoint2 = nLinePw2 + nOuterPoint

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

     allocate(StateGm1_VI(1:iUGmLast, nPoint), CoordXyPw1_DI(nCoord, nPoint), &
          iNodeTriangle1_II(3, 2*nPoint1))
     allocate(StateGm2_VI(1:iUGmLast, nPoint), CoordXyPw2_DI(nCoord, nPoint))
     if(nLinePw2 /=0)then
        allocate(iNodeTriangle2_II(3, 2*nPoint2))
     end if
  end if

  if(nLinePw1+nLinePw2 /= nFieldLine .or. nVarPw /= nVar)then
     write(*,*)NameSub,' ERROR: nLinePw=',nLinePw1,'/= nFieldLine=',nFieldLine,&
          ' or nVarPw=',nVarPw,' /= nVar=',nVar
     call CON_stop(NameSub,' nFieldLine or nVar has changed')
  end if

  ! Convert to X, Y on a unit sphere (this will be used for triangulation)

  where(Buffer_VI(Theta_,:) > cHalfPi)
     CoordXyPw2_DI(x_,1:nFieldline) =  &
          sin(Buffer_VI(Theta_,:)) * cos(Buffer_VI(Phi_,:))
     CoordXyPw2_DI(y_,1:nFieldline) =  &
          sin(Buffer_VI(Theta_,:)) * sin(Buffer_VI(Phi_,:))
  elsewhere
     CoordXyPw1_DI(x_,1:nFieldline) =  &
          sin(Buffer_VI(Theta_,:)) * cos(Buffer_VI(Phi_,:))
     CoordXyPw1_DI(y_,1:nFieldline) =  &
          sin(Buffer_VI(Theta_,:)) * sin(Buffer_VI(Phi_,:))
  end where

  if(UseMultiIon)then
     do i = iRhoGmFirst, iRhoGmLast
        where(Buffer_VI(Theta_,:)< cHalfPi)
           StateGm1_VI(i,1:nFieldline) &
                = Buffer_VI(i+iRhoPwFirst-1,:)*Si2No_V(UnitRho_)
        elsewhere
           StateGm2_VI(i,1:nFieldline) &
                = Buffer_VI(i+iRhoPwFirst-1,:)*Si2No_V(UnitRho_)
        end where
     end do
     
     do i = iUGmFirst, iUGmLast
        where(Buffer_VI(Theta_,:)< cHalfPi)
           StateGm1_VI(i, 1:nFieldline) &
                = Buffer_VI(i-iUGmFirst+iUPwFirst, :) * Si2No_V(UnitU_)
        elsewhere
           StateGm2_VI(i, 1:nFieldline) &
                = Buffer_VI(i-iUGmFirst+iUPwFirst, :) * Si2No_V(UnitU_)
        end where
     end do
  else
     ! Set total 
     if(UseMultiSpecies)then
        do i=iRhoGmFirst, iRhoGmLast
           where(Buffer_VI(Theta_,:)< cHalfPi)
              StateGm1_VI(i, 1:nFieldline) &
                   = Buffer_VI(i+iRhoGmFirst-1,:) * Si2No_V(UnitRho_)
           elsewhere
              StateGm2_VI(i, 1:nFieldline) &
                   = Buffer_VI(i+iRhoGmFirst-1,:) * Si2No_V(UnitRho_)
           end where
        end do
     else
        ! Total density in normalized units
        where(Buffer_VI(Theta_,:)< cHalfPi)
           StateGm1_VI(iRhoGmFirst, 1:nFieldline) &
                = sum(Buffer_VI(iRhoPwFirst:iRhoPwLast,:),dim=1) * Si2No_V(UnitRho_)
        elsewhere
           StateGm2_VI(iRhoGmFirst, 1:nFieldline) &
                = sum(Buffer_VI(iRhoPwFirst:iRhoPwLast,:),dim=1) * Si2No_V(UnitRho_)
        end where
     end if

     ! Field aligned velocity = total moment/total density
      where(Buffer_VI(Theta_,:)< cHalfPi)
         StateGm1_VI(iUGmFirst,1:nFieldline) &
              = sum(Buffer_VI(iRhoPwFirst:iRhoPwLast,:) &
              *     Buffer_VI(iUPwFirst:iUPwLast,:), dim=1) &
              / sum(Buffer_VI(iRhoPwFirst:iRhoPwLast,:), dim=1) &
              * Si2No_V(UnitU_) 
      elsewhere
         StateGm2_VI(iUGmFirst,1:nFieldline)&
              = sum(Buffer_VI(iRhoPwFirst:iRhoPwLast,:) &
              *     Buffer_VI(iUPwFirst:iUPwLast,:), dim=1) &
              / sum(Buffer_VI(iRhoPwFirst:iRhoPwLast,:), dim=1) &
              * Si2No_V(UnitU_)
      end where
  end if

  ! Reduce densities according to a 1/r3 law based on magnetic field strength
  ! and assuming that PWOM solves up to R=2.25
  !Factor = (2.25/rBody)**3
  !StateGm_VI(iRhoGmFirst:iRhoGmLast,1:nLinePw) = &
  !     StateGm_VI(iRhoGmFirst:iRhoGmLast,1:nLinePw)*Factor

  ! Set coordinates for the outer points
  tmp1_array = cTwoPi
  where(Buffer_VI(Theta_,:)>cHalfPi)
     tmp1_array = Buffer_VI(Theta_,:)
  end where
  SinThetaOuter2 = sin(minval(tmp1_array)-dThetaOuter)

  tmp1_array = 0.0
  where(Buffer_VI(Theta_,:)<cHalfPi)
     tmp1_array = Buffer_VI(Theta_,:)
  end where  
  SinThetaOuter1 = sin(maxval(tmp1_array)+dThetaOuter)

  do i = 1, nOuterPoint
     CoordXyPw1_DI(x_,nLinePw1+i) = &
          SinThetaOuter1 * cos(i*cTwoPi/nOuterPoint)
     CoordXyPw1_DI(y_,nLinePw1+i) = &
          SinThetaOuter1 * sin(i*cTwoPi/nOuterPoint)

     if(nLinePw2 /=0)then
        CoordXyPw2_DI(x_,nLinePw1+nLinePw2+i) = &
             SinThetaOuter2 * cos(i*cTwoPi/nOuterPoint)
        CoordXyPw2_DI(y_,nLinePw1+nLinePw2+i) = &
             SinThetaOuter2 * sin(i*cTwoPi/nOuterPoint)
     end if
  end do

  call calc_triangulation(nPoint1, CoordXyPw1_DI(:,1:nPoint1), &
       iNodeTriangle1_II, nTriangle1)
  
  if(nLinePw2 /=0)then
     call calc_triangulation(nPoint2, CoordXyPw2_DI(:,nLinePw1+1:nLinePw1+nPoint2),&
          iNodeTriangle2_II, nTriangle2)
  end if

  if(DoTestMe)then
     write(*,*)NameSub,': nVarPw, nLinePw, nSpeciesPw=',&
          nVarPw, nLinePw2, nSpeciesPw
     write(*,*)NameSub,': nPoint2, nTriangle2=',nPoint2, nTriangle2
     write(*,*)NameSub,': CoordXyPw_DI'
     do i=nLinePw1+1,nLinePw1+nPoint2
        write(*,*) i, CoordXyPw2_DI(:,i)
     end do
     do i=1,nPoint1
        write(*,*) i, CoordXyPw1_DI(:,i)
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
     do iPoint=nLinePw1+1, nPoint1
        StateGm1_VI(iUGmFirst:iUGmLast, iPoint) = 0.0
        if(UseMultiSpecies)then
           StateGm1_VI(iRhoGmFirst:iRhoGmLast, iPoint) = &
                State_V(SpeciesFirst_:SpeciesLast_)
        elseif(UseMultiIon)then
           StateGm1_VI(iRhoGmFirst:iRhoGmLast, iPoint) = State_V(iRhoIon_I)
        else
           StateGm1_VI(iRhoGmFirst, iPoint) = State_V(Rho_)
        end if
     end do
     
      if (nLinePw2 /=0) then                                                                                            
        do iPoint=nLinePw2+1, nPoint2
            StateGm2_VI(iUGmFirst:iUGmLast, nLinePw1+iPoint) = 0.0
            if(UseMultiSpecies)then
               StateGm2_VI(iRhoGmFirst:iRhoGmLast, nLinePw1+iPoint) = &
                    State_V(SpeciesFirst_:SpeciesLast_)
            elseif(UseMultiIon)then
               StateGm2_VI(iRhoGmFirst:iRhoGmLast, nLinePw1+iPoint) = State_V(iRhoIon_I)
            else
               StateGm2_VI(iRhoGmFirst, nLinePw1+iPoint) = State_V(Rho_)
            end if
        end do
     end if

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

  if(Xyz_D(3) < 0..and. nLinePw2 /=0)then
     call find_triangle(&
          nPoint2, nTriangle2, Xy_D, CoordXyPw2_DI(:,nLinePw1+1:nLinePw1+nPoint2), &
          iNodeTriangle2_II(:,1:nTriangle2), &
          iNode1, iNode2, iNode3, Area1, Area2, Area3, IsTriangleFound)
  else
     call find_triangle(&
          nPoint1, nTriangle1, Xy_D, CoordXyPw1_DI(:,1:nPoint1), &
          iNodeTriangle1_II(:,1:nTriangle1), &
          iNode1, iNode2, iNode3, Area1, Area2, Area3, IsTriangleFound)
  end if
  ! Point is not covered: leave the input state variables alone
  if (.not.IsTriangleFound) RETURN

  ! Calculate field aligned momentum vector
  call get_b0(CoordIn_D(x_), CoordIn_D(y_), CoordIn_D(z_), B0_D)
  B0_D = B0_D / sqrt(sum(B0_D**2))

  ! Make sure unit vector is pointing outward
  if(sum(B0_D*CoordIn_D) < 0.)B0_D = -B0_D

  ! Put into a temporary array StateGm_VI                                                                            
  if(Xyz_D(3) < 0 .and. nLinePw2 /=0) then
     allocate(StateGm_VI(1:iUGmLast, nPoint2))
     StateGm_VI = StateGm2_VI(:,nLinePw1+1:nLinePw1+nPoint2)
  else
     allocate(StateGm_VI(1:iUGmLast, nPoint1))
     StateGm_VI = StateGm1_VI(:,1:nPoint1)
  end if

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

  deallocate(StateGm_VI)

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
