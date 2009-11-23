module ModHeidiHydrogenGeo

  implicit none

contains

  subroutine get_rairden_density(nPoint, nRadialPoint,L, RhoH_I)

    use ModIoUnit,  ONLY : UnitTmp_
    use ModHeidiIO, ONLY : NameInputDirectory,NameRun
    integer, intent(in)  :: nPoint                  ! number of points along the field line
    integer, intent(in)  :: nRadialPoint            ! number of points in the radial direction = NR
    real,    intent(in)  :: L                       ! x^2+y^2
    real,    intent(out) :: RhoH_I(nPoint)          !interpolated density
    integer, parameter   :: nUniform =100           !number of points on new refined grid
    integer, parameter   :: nRairden = 82           !number of radial point in the Rairden hydrogen file
    character (len=100)  :: StringHeader
    real                 :: Lat
    real                 :: RadialDistance_I(nPoint)
    real                 :: LatMax, LatMin, dLat
    real                 :: Rad_I(nUniform), Weight, RhoHUniform_I(nUniform)
    real                 :: RairdenDistance_I(nRairden),RairdenDensity
    real                 :: LnRairdenDensity_I(nRairden)
    real                 :: rMax, rMin, dR
    real                 :: LnRhoH,r
    integer              :: iRairden,i
    integer              :: iR,iStep
    
    !------------------------------------------------------------------------------

    ! Read the Rairden Geocorona Hydrogen density file.
    
    open(UNITTMP_,FILE=NameInputDirectory//'RairdenHydrogenGeocorona.dat',status='old')
    
    do i = 1, 4
       read (UnitTmp_,*) StringHeader
    end do

    do iRairden = 1, nRairden
       read(UnitTmp_,*) RairdenDistance_I(iRairden), RairdenDensity
       LnRairdenDensity_I(iRairden) = log(RairdenDensity*10.**6)    ! need the density in m^-3 NOT cm^-3!!!!
    end do
    close(UnitTmp_) 

    !interpolate the hydrogen density to a new, well refined grid.

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
   
     !write the new interpolated values to a file
!!$    open(unit=2, file='density.dat')
!!$    write (2,*)'values for rairden density and interpolated density'
!!$    write (2,*)' R, LnH, H'
!!$    do i = 1, nUniform
!!$       write(2,*) Rad_I(i),RhoHUniform_I(i)
!!$    end do
!!$    close(2)

    ! calculate distance to any point along the field line 
    LatMax = acos(sqrt(1./L))
    LatMin = -LatMax
    dLat   = (LatMax - LatMin)/(nPoint - 1)
    do iStep = 1, nPoint
       Lat = LatMin + (iStep-1)*dLat
       RadialDistance_I(iStep) = L*cos(Lat)*cos(Lat)
    end do
    
    do iStep =1, nPoint
       r = RadialDistance_I(iStep)
       i = 1 + abs((r - rMin))/dR
       Weight = (r-Rad_I(i))/dR
       RhoH_I(iStep) = Weight*RhoHUniform_I(i+1) + (1-Weight)*RhoHUniform_I(i)
    end do
  
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
    real                 :: Inv2L
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
    real                 :: Rho_II(nPoint)    
    real                 :: PitchAngle
    real                 :: dLength_III(nPoint-1,nR,nPhi)      ! Length interval between i and i+1  
    real                 :: Length_III(nPoint,nR,nPhi) 
    real                 :: RadialDistance_III(nPoint,nR,nPhi)
    real                 :: AvgHDensity(nR),HDensity(nR,nPoint)
    real                 :: bMirror
    integer              :: iMirror_I(2)
    real, parameter      :: Pi = 3.141592654   
    real                 :: Ds_I(nPoint)
    integer              :: iPoint
    !-----------------------------------------------------------------------------

    L_I(1) = 2.0
    PitchAngle = Pi/10.
   
    call initialize_b_field(L_I, Phi_I, nPoint, nR, nPhi, bFieldMagnitude_III, &
          RadialDistance_III,Length_III, dLength_III,GradBCrossB_VIII,GradB_VIII)
    
    call find_mirror_points (iPoint,  PitchAngle, bFieldMagnitude_III, &
                      bMirror,iMirror_I)
    
    call get_rairden_density(nPoint, nR,L_I(1), Rho_II)
    
    call get_hydrogen_density(nPoint, L_I(1), bFieldMagnitude_III, bMirror,iMirror_I(1),dLength_III,Rho_II,AvgHDensity(1))
       
  end subroutine test_density

end module ModHeidiHydrogenGeo
