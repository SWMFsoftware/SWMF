module ModHeidiHydrogenGeo

  implicit none

contains

 !============================================================
    subroutine get_rairden_density(nPoint,RadialDistance_I, L, RhoH_I)
    use ModIoUnit,     ONLY : UnitTmp_
    use ModHeidiInput, ONLY: TypeHModel
    use ModHeidiIO,    ONLY: NameInputDirectory,NameRun
    use ModHeidiMain,  ONLY: LZ, Re
    use ModHeidiSize,  ONLY: nR
    use NeutralHydrogenModel

    use ModNumConst, ONLY:cPi
    
    integer, intent(in)  :: nPoint                 ! number of points along the field line
    real   , intent(in)  :: RadialDistance_I(nPoint), L
    real,    intent(out) :: RhoH_I(nPoint)          !interpolated density
    integer, parameter   :: nUniform =100           !number of points on new refined grid
    integer, parameter   :: nRairden = 82           !number of radial point in the Rairden hydrogen file
    character (len=100)  :: StringHeader
    real                 :: Rad_I(nUniform), Weight, RhoHUniform_I(nUniform)
    real                 :: RairdenDistance_I(nRairden),RairdenDensity
    real                 :: LnRairdenDensity_I(nRairden)
    real                 :: rMax, rMin, dR
    real                 :: LnRhoH,r
    integer              :: iRairden,i
    integer              :: iR,iPoint



    real :: HDensity(40), R_I(40)
    !------------------------------------------------------------------------------
   
    select case(TypeHModel)
    case('Rairden')

       RhoHUniform_I = 0.0
       !\
       ! Read the Rairden Geocorona Hydrogen density file.
       !/
       open(UNITTMP_,FILE=NameInputDirectory//'RairdenHydrogenGeocorona.dat',status='old')
       
       do i = 1, 4
          read (UnitTmp_,*) StringHeader
       end do
       
       do iRairden = 1, nRairden
          read(UnitTmp_,*) RairdenDistance_I(iRairden), RairdenDensity
          LnRairdenDensity_I(iRairden) = log(RairdenDensity*10.**6)    ! need the density in m^-3 NOT cm^-3!!!!
       end do
       close(UnitTmp_) 
       
       !\
       !Interpolate the hydrogen density to a new, well refined grid.
       !/
       rMax = RairdenDistance_I(nRairden)
       rMin = RairdenDistance_I(1)  
       dR = (rMax-rMin)/(nUniform-1)
       
       do i =1, nUniform
          Rad_I(i) = rMin + dR*(i-1)
          do iRairden = 1, nRairden
             if (RairdenDistance_I(iRairden)>Rad_I(i)) then
                iR = iRairden-1
                LnRhoH = LnRairdenDensity_I(iR)+ (Rad_I(i)-RairdenDistance_I(iR))*&
                     (LnRairdenDensity_I(iR+1) - LnRairdenDensity_I(iR))/&
                     (RairdenDistance_I(iR+1) - RairdenDistance_I(iR))
                RhoHUniform_I(i) = exp(LnRhoH)
                EXIT
             end if
          end do
       end do
       
       do iPoint =1, nPoint
          r = RadialDistance_I(iPoint)
          if (r>=rmax) then
             RhoH_I(iPoint) = RhoHUniform_I(nUniform)
          else
             i = 1 + abs((r - rMin))/dR
             Weight = (r-Rad_I(i))/dR
             RhoH_I(iPoint) = Weight*RhoHUniform_I(i+1) + (1-Weight)*RhoHUniform_I(i)
          end if
       end do
       


    case('Ostgaard')
      
          do iR = 1, nR
             call get_ostgaard_density(LZ(iR), RhoHUniform_I(iR))
          end do
          

          rMax = LZ(nR)
          rMin = LZ(1)
          dR = (rMax-rMin)/(nR-1)
          
          do iPoint =1, nPoint
             r = RadialDistance_I(iPoint)
             if (r>=rmax) then
                RhoH_I(iPoint) = RhoHUniform_I(nUniform)
          else
             i = 1 + abs((r - rMin))/dR
             Weight = (r-Rad_I(i))/dR
             RhoH_I(iPoint) = Weight*RhoHUniform_I(i+1) + (1-Weight)*RhoHUniform_I(i)
          end if
       end do
       
    case('Hodges')

       ! Include here the Hodges densities.

!!$       call get_hodges_density(cPi/2., 0.0, HDensity, R_I)
!!$
!!$       open(unit=3, file='Hodges.dat')
!!$       write(3,*) 'Header'
!!$       write(3,*) 'R, HDens, R/Re, H'
!!$       do i =1, 40
!!$          write(3,*) R_I(i), HDensity(i),  R_
!!$       end do
!!$
!!$       STOP

    end select


  end subroutine get_rairden_density

  !=================================================================================
  subroutine get_hydrogen_density(nPoint, L, bField_I, bMirror,iMirror_I,dLength_I,HDensity_I,AvgHDensity)

    integer, intent(in)  :: nPoint
    real,    intent(in)  :: L
    real,    intent(in)  :: bField_I(nPoint) 
    real,    intent(in)  :: bMirror  
    integer, intent(in)  :: iMirror_I(2)
    real,    intent(in)  :: dLength_I(nPoint-1)
    real,    intent(in)  :: HDensity_I(nPoint) 
    real,    intent(out) :: AvgHDensity
    real                 :: DeltaS1, DeltaS2, b1, b2, Coeff
    integer              :: iPoint, iFirst, iLast
    !-----------------------------------------------------------------------------

    iFirst = iMirror_I(1)
    iLast  = iMirror_I(2)

    AvgHDensity = 0.0

    if (iFirst > iLast) RETURN

    Coeff = sqrt(bMirror)

    DeltaS1 = abs((bMirror-bField_I(iFirst))*&
         (dLength_I(iFirst-1))/(bField_I(iFirst-1)-bField_I(iFirst)))

    AvgHDensity= AvgHDensity + HDensity_I(iFirst)*Coeff*2.*DeltaS1/(sqrt(bMirror-bField_I(iFirst)))

    do iPoint = iFirst, iLast-1
       b1 = bField_I(iPoint)
       b2 = bField_I(iPoint+1)

       AvgHDensity = AvgHDensity + HDensity_I(iPoint)*Coeff*2.*dLength_I(iPoint)/(b1 - b2) &
            *( sqrt(bMirror  - b2) - sqrt(bMirror  - b1) )

    end do


    DeltaS2 = abs((bMirror-bField_I(iLast))*(dLength_I(iLast))/(bField_I(iLast+1)-bField_I(iLast)))


    AvgHDensity = AvgHDensity + HDensity_I(iLast)* Coeff*2.*DeltaS2/(sqrt(bMirror-bField_I(iLast)))

  end subroutine get_hydrogen_density
  !===============================================================

  subroutine test_density
    use ModNumConst,   ONLY: cPi
    use ModHeidiBField, ONLY: initialize_b_field,find_mirror_points


    integer, parameter   :: nPoint = 11
    integer, parameter   :: nPitch = 1
    integer, parameter   :: nR = 1
    integer, parameter   :: nPhi =1
    real                 :: L_I(1)
    real                 :: Phi_I(nPhi)    
    real                 :: bFieldMagnitude_III(nPoint,nR,nPhi) 
    real                 :: GradBCrossB_VIII(3,nPoint,nR,nPhi)
    real                 :: GradB_VIII(3,nPoint,nR,nPhi)
    real                 :: Rho_I(nPoint)    
    real                 :: PitchAngle
    real                 :: dLength_III(nPoint-1,nR,nPhi)      ! Length interval between i and i+1  
    real                 :: Length_III(nPoint,nR,nPhi) 
    real                 :: RadialDistance_III(nPoint,nR,nPhi)
    real                 :: AvgHDensity(nR)
    real                 :: bMirror
    integer              :: iMirror_I(2)
    integer              :: iPoint
    real                 :: dBdt_III(nPoint,nR,nPhi)
    !-----------------------------------------------------------------------------

    L_I(1) = 2.0
    PitchAngle = cPi/10.
    
    call initialize_b_field(L_I, Phi_I, nPoint, nR, nPhi, bFieldMagnitude_III, &
         RadialDistance_III,Length_III, dLength_III,GradBCrossB_VIII,GradB_VIII,dBdt_III)

    call find_mirror_points (iPoint,  PitchAngle, bFieldMagnitude_III, &
         bMirror,iMirror_I)
   
    call get_rairden_density(nPoint,RadialDistance_III, L_I(1), Rho_I)

     
    call get_hydrogen_density(nPoint, L_I(1), bFieldMagnitude_III, bMirror,iMirror_I(1),dLength_III,Rho_I,AvgHDensity(1))

  end subroutine test_density

end module ModHeidiHydrogenGeo
