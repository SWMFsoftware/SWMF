!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModPwGrid

  use ModNumConst,   ONLY: cDegToRad
  use ModMultiFluid, ONLY: nIonFluid
  
  implicit none
  save

  character (len=3) :: NamePwCoord = '???'
  integer, parameter :: nCoord=3

  ! StateGm_VI stores the PW variables that occur in GM too
  ! These are interpolated then to the GM grid
  real, allocatable :: CoordXyzPw_DI(:,:), CoordXyPw_DI(:,:), StateGm_VI(:,:),&
       CoordXyzPw1_DI(:,:), CoordXyzPw2_DI(:,:), &
       StateGm1_VI(:,:),StateGm2_VI(:,:)
  
  integer, allocatable :: list1_I(:),lptr1_I(:), lend1_I(:),&
       list2_I(:),lptr2_I(:), lend2_I(:)  

  integer, parameter :: nOuterPoint=12
  real,    parameter :: dThetaOuter=5.0*cDegToRad

  integer :: iCoupleFluid(nIonFluid) = -1

  integer :: nVarPw=-1, nSpeciesPw=-1, nLinePw=-1, nPoint=-1, nTriangle=-1
  integer :: iRhoPwFirst=-1, iRhoPwLast=-1, iUPwFirst=-1, iUPwLast=-1
  integer :: iRhoGmFirst=-1, iRhoGmLast=-1, iUGmFirst=-1, iUGmLast=-1
  integer :: nPoint1=-1, nPoint2=-1,nTriangle1=-1,nTriangle2=-1
  integer :: nLinePw1 = -1, nLinePw2 = -1

end module ModPwGrid

module GM_couple_pw

  implicit none

contains 

  subroutine GM_get_for_pw(nTotalLine, Buffer_I)

    use ModVarIndexes, ONLY: p_
    use ModPwGrid,     ONLY: CoordXyzPw_DI
    use ModMain,       ONLY: nBlock
    use ModProcMH,     ONLY: nProc, iComm
    use ModMpi
    use ModPhysics, ONLY: No2Si_V, UnitP_
    implicit none

    integer,intent(in)  :: nTotalLine
    real,   intent(out) :: Buffer_I(nTotalLine)
    real, allocatable   :: LocalBuffer_I(:)
    real                :: WeightAndP_I(2)
    integer             :: iLine,iError
    !----------------------------------------------------------------------------

    allocate(LocalBuffer_I(nTotalLine))

    !Get the pressure at each field line location.
    do iLine=1,nTotalLine
       call get_point_data(&
            0.0, CoordXyzPW_DI(:,iLine), 1, nBlock, p_, p_, WeightAndP_I)
       LocalBuffer_I(iLine) = WeightAndP_I(2)
    enddo

    Buffer_I = LocalBuffer_I

    ! Sum contributions from all PE-s to processor 0
    if(nProc > 1) call MPI_reduce(LocalBuffer_I, Buffer_I, nTotalLine,&
         MPI_REAL, MPI_SUM, 0, iComm, iError)

    ! Convert pressure to SI units
    Buffer_I = Buffer_I * No2Si_V(UnitP_)

    deallocate(LocalBuffer_I)

  end subroutine GM_get_for_pw
  subroutine GM_put_from_pw(Buffer_VI, nVar, nFieldLine, Name_V)

    use CON_coupler, ONLY: PW_, Grid_C
    use ModMain, ONLY: x_, y_,z_, TypeCoordSystem, Time_Simulation
    use ModVarIndexes, ONLY: UseMultiSpecies, SpeciesFirst_, &
         SpeciesLast_, NameVar_V
    use ModMultiFluid, ONLY: UseMultiIon, nIonFluid
    use ModPhysics, ONLY: Si2No_V, UnitRho_, UnitU_, rCurrents
    use ModPwGrid
    use ModNumConst, ONLY: cTwoPi,cHalfPi
    use ModTriangulateSpherical,ONLY:trmesh, trplot
    use CON_axes, ONLY: transform_matrix
    use ModCoordTransform, ONLY: dir_to_xyz
    use CON_planet_field,  ONLY: map_planet_field
    use ModIoUnit, ONLY: UnitTmp_

    implicit none
    character (len=*),parameter :: NameSub='GM_put_from_pw'

    integer, parameter :: Theta_=1, Phi_=2

    integer, intent(in)           :: nVar, nFieldLine
    real, intent(in)              :: Buffer_VI(nVar, nFieldLine)
    character (len=*), intent(in) :: Name_V(nVar)

    logical, save :: DoInitialize = .true.

    integer :: i, iLine, iHemisphere, j
    real    :: GmPw_DD(3,3), Theta, Phi, XyzPw_D(3), &
         SinThetaOuter1,SinThetaOuter2,CosThetaOuter1,CosThetaOuter2,&
         tmp1_array(nFieldLine)
    integer :: Tmp_array(nFieldLine)
    logical :: DoTest, DoTestMe
    integer :: iError=0 !error counter for triangulation subroutine
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

          ! Build iCoupleFluid, which connects PW fluids with GM fluids.
          ! Assumption: MHD is using STANDARD FLUID NAMES and PW only has
          ! Op, Hep, and Hp to share.
          j = 1
          do i=1, size(NameVar_V)
             if(DoTestMe)write(*,*) 'Checking var ', trim(NameVar_V(i)),'.'
             select case( trim(NameVar_V(i)) )
             case('HpRho')
                iCoupleFluid(j) = iRhoPwFirst
             case('OpRho')
                iCoupleFluid(j) = iRhoPwFirst + 1
             case('HepRho')
                iCoupleFluid(j) = iRhoPwFirst + 2
             end select
             ! Keep track of which ion species we are on by counting the
             ! number of variables that have the name <Species>Rho.
             if (index(NameVar_V(i), 'Rho')>1) j = j + 1
          end do

          if(DoTestMe)write(*,*) 'nIonFluid = ', nIonFluid
          if(DoTestMe)write(*,*) 'iCoupleFluid = (',iCoupleFluid,')'

       else
          ! Total density and velocity
          iRhoGmFirst = 1
          iRhoGmLast  = 1
          iUGmFirst   = 2
          iUGmLast    = 2
       end if

       allocate(StateGm1_VI(1:iUGmLast, nPoint), CoordXyzPw1_DI(nCoord, nPoint), &
            list1_I(6*(nPoint1-2)),lptr1_I(6*(nPoint1-2)), lend1_I(nPoint1))
       allocate(StateGm2_VI(1:iUGmLast, nPoint), CoordXyzPw2_DI(nCoord, nPoint))
       !Initialize Arrays to zero
       StateGm1_VI(:,:) =0.0
       StateGm2_VI(:,:) =0.0
       CoordXyzPw1_DI(:,:)=0.0
       CoordXyzPw2_DI(:,:)=0.0
       list1_I(:)=0
       !list2_I(:)=0
       lend1_I(:)=0
       !lend2_I(:)=0
       lptr1_I(:)=0
       !lptr2_I(:)=0
       if(nLinePw2 /=0)then
          allocate(list2_I(6*(nPoint1-2)),lptr2_I(6*(nPoint1-2)), &
               lend2_I(nPoint1))
       end if
    end if

    if(nLinePw1+nLinePw2 /= nFieldLine .or. nVarPw /= nVar)then
       write(*,*)NameSub,' ERROR: nLinePw=',nLinePw1,'/= nFieldLine=',nFieldLine,&
            ' or nVarPw=',nVarPw,' /= nVar=',nVar
       call CON_stop(NameSub,' nFieldLine or nVar has changed')
    end if

    ! Convert to X, Y on a unit sphere (this will be used for triangulation)

    where(Buffer_VI(Theta_,:) > cHalfPi)
       CoordXyzPw2_DI(x_,1:nFieldline) =  &
            sin(Buffer_VI(Theta_,:)) * cos(Buffer_VI(Phi_,:))
       CoordXyzPw2_DI(y_,1:nFieldline) =  &
            sin(Buffer_VI(Theta_,:)) * sin(Buffer_VI(Phi_,:))
       CoordXyzPw2_DI(z_,1:nFieldline) =  &
            cos(Buffer_VI(Theta_,:)) 
    elsewhere
       CoordXyzPw1_DI(x_,1:nFieldline) =  &
            sin(Buffer_VI(Theta_,:)) * cos(Buffer_VI(Phi_,:))
       CoordXyzPw1_DI(y_,1:nFieldline) =  &
            sin(Buffer_VI(Theta_,:)) * sin(Buffer_VI(Phi_,:))
       CoordXyzPw1_DI(z_,1:nFieldline) =  &
            cos(Buffer_VI(Theta_,:)) 
    end where

    if(UseMultiIon)then
       do i = iRhoGmFirst, iRhoGmLast
          ! Set values based on if fluid is shared by PW & GM.
          if (iCoupleFluid(i)==-1) then ! Fluid NOT in PW.
             where(Buffer_VI(Theta_,:)< cHalfPi) ! Northern Hemisphere.
                StateGm1_VI(i,           1:nFieldline) = 1.67E-24 * Si2No_V(UnitRho_)
                StateGm1_VI(i+nIonFluid, 1:nFieldLine) = 0.0
             elsewhere                           ! Southern Hemisphere
                StateGm2_VI(i,           1:nFieldline) = 1.67E-24 * Si2No_V(UnitRho_)
                StateGm2_VI(i+nIonFluid, 1:nFieldLine) = 0.0
             end where
          else                          ! Fluid IS in PW.
             where(Buffer_VI(Theta_,:)< cHalfPi) ! Northern Hemisphere.
                StateGm1_VI(i,           1:nFieldLine) = &
                     Buffer_VI(iCoupleFluid(i),            :) * Si2No_V(UnitRho_)
                StateGm1_VI(i+nIonFluid, 1:nFieldLine) = &
                     Buffer_VI(iCoupleFluid(i)+nSpeciesPw, :) * Si2No_V(UnitU_)
             elsewhere                           ! Southern Hemisphere
                StateGm2_VI(i,           1:nFieldLine) = &
                     Buffer_VI(iCoupleFluid(i),            :) * Si2No_V(UnitRho_)
                StateGm2_VI(i+nIonFluid, 1:nFieldLine) = &
                     Buffer_VI(iCoupleFluid(i)+nSpeciesPw, :) * Si2No_V(UnitU_)
             end where
          end if
       end do

       if(DoTestMe) then
          write(*,*)'Max Hp Density:', maxval(Buffer_VI(3, :))
          write(*,*)'Max Op Density:', maxval(Buffer_VI(4, :))
          write(*,*)'Max He Density:', maxval(Buffer_VI(5, :))
          write(*,*)'Max Sw Density:', maxval(StateGm1_VI(1, :))/Si2No_V(UnitRho_)

          write(*,*)'Max Hp Velocity:', maxval(Buffer_VI(6, :))
          write(*,*)'Max Op Velocity:', maxval(Buffer_VI(7, :))
          write(*,*)'Max He Velocity:', maxval(Buffer_VI(8, :))
          write(*,*)'Max Sw Velocity:', maxval(StateGm1_VI(4, :))/Si2No_V(UnitU_)
       end if

    else
       ! Set total 
       if(UseMultiSpecies)then
          do i=iRhoGmFirst, iRhoGmLast
             where(Buffer_VI(Theta_,:)< cHalfPi)
                StateGm1_VI(i, 1:nFieldline) &
                     = Buffer_VI(i+iRhoPwFirst-1,:) * Si2No_V(UnitRho_)
             elsewhere
                StateGm2_VI(i, 1:nFieldline) &
                     = Buffer_VI(i+iRhoPwFirst-1,:) * Si2No_V(UnitRho_)
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
    CosThetaOuter2 = cos(minval(tmp1_array)-dThetaOuter)

    tmp1_array = 0.0
    where(Buffer_VI(Theta_,:)<cHalfPi)
       tmp1_array = Buffer_VI(Theta_,:)
    end where
    SinThetaOuter1 = sin(maxval(tmp1_array)+dThetaOuter)
    CosThetaOuter1 = cos(maxval(tmp1_array)+dThetaOuter)

    do i = 1, nOuterPoint
       CoordXyzPw1_DI(x_,nLinePw1+i) = &
            SinThetaOuter1 * cos(i*cTwoPi/nOuterPoint)
       CoordXyzPw1_DI(y_,nLinePw1+i) = &
            SinThetaOuter1 * sin(i*cTwoPi/nOuterPoint)
       CoordXyzPw1_DI(z_,nLinePw1+i) = CosThetaOuter1
       if(nLinePw2 /=0)then
          CoordXyzPw2_DI(x_,nLinePw1+nLinePw2+i) = &
               SinThetaOuter2 * cos(i*cTwoPi/nOuterPoint)
          CoordXyzPw2_DI(y_,nLinePw1+nLinePw2+i) = &
               SinThetaOuter2 * sin(i*cTwoPi/nOuterPoint)
          CoordXyzPw2_DI(z_,nLinePw1+nLinePw2+i) = CosThetaOuter2
       end if
    end do


    !
    !  Create the triangulation.
    !
    call trmesh ( nPoint1, CoordXyzPw1_DI(x_,1:nPoint1), &
         CoordXyzPW1_DI(y_,1:nPoint1), CoordXyzPW1_DI(z_,1:nPoint1), &
         list1_I, lptr1_I, lend1_I, iError )
    if ( iError == -2 ) then
       write(*,*)NameSub, &
            ' WARNING: Error in TRMESH, First three nodes are collinear'
       call CON_stop(NameSub,' Problem With Triangulation')
    else if ( iError > 0 ) then
       write(*,*)NameSub, &
            ' ERROR: Error in TRMESH, Duplicate nodes encountered'
       call CON_stop(NameSub,' Problem With Triangulation')
    end if
    !


    if(nLinePw2 /=0)then
       call trmesh ( nPoint2, CoordXyzPw2_DI(x_,nLinePw1+1:nLinePw1+nPoint2), &
            CoordXyzPW2_DI(y_,nLinePw1+1:nLinePw1+nPoint2), &
            CoordXyzPW2_DI(z_,nLinePw1+1:nLinePw1+nPoint2), &
            list2_I, lptr2_I, lend2_I, iError )
       if ( iError == -2 ) then
          write(*,*)NameSub, &
               ' WARNING: Error in TRMESH, First three nodes are collinear'
          call CON_stop(NameSub,' Problem With Triangulation')
       else if ( iError > 0 ) then
          write(*,*)NameSub, &
               ' ERROR: Error in TRMESH, Duplicate nodes encountered'
          call CON_stop(NameSub,' Problem With Triangulation')
       end if
    end if

    if(DoTestMe)then
       write(*,*)NameSub,': nVarPw, nLinePw, nSpeciesPw=',&
            nVarPw, nLinePw2, nSpeciesPw
       write(*,*)NameSub,': nPoint2, nTriangle2=',nPoint2, nTriangle2
       write(*,*)NameSub,': CoordXyzPw_DI'
       do i=nLinePw1+1,nLinePw1+nPoint2
          write(*,*) i-nLinePw1, CoordXyzPw2_DI(:,i)
       end do
       do i=1,nPoint1
          write(*,*) i, CoordXyzPw1_DI(:,i)
       end do
       write(*,*) NameSub,': Writing plot of triangulation'
       open ( UnitTmp_, file = 'TestTriangulation.eps' )
       call trplot ( UnitTmp_, 7.5, 90.0, 0.0, 90.0, nPoint1, CoordXyzPw1_DI(x_,1:nPoint1), &
            CoordXyzPw1_DI(y_,1:nPoint1), CoordXyzPw1_DI(z_,1:nPoint1), list1_I, lptr1_I, &
            lend1_I, 'test1 triangulation',.true., iError )
       close(UnitTmp_)
    end if



  end subroutine GM_put_from_pw

 end module GM_couple_pw
