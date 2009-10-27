module ModHeidiHydrogenGeo

  implicit none

contains

  subroutine get_rairden_density(nPoint, nRadialPoint,L_I, RhoH_II)

    use ModIoUnit, ONLY : UnitTmp_
    integer, intent(in)  :: nPoint                  ! number of points along the field line
    integer, intent(in)  :: nRadialPoint           ! number of points in the radial direction = NR
    real,    intent(in)  :: L_I(nRadialPoint)      ! x^2+y^2
    real,    intent(out) :: RhoH_II(nRadialPoint, nPoint)  !interpolated density
    integer, parameter   :: nUniform =100          !number of points on new refined grid
    integer, parameter   :: nRairden = 81          !number of radial point in the Rairden hydrogen file
    character (len=100)  :: StringHeader
    real                 :: Lat
    real                 :: RadialDistance_II(nRadialPoint, nPoint)
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

    open(Unit=UnitTmp_,file="RairdenHydrogenGeocorona.dat", status = 'old')

    do i = 1, 5
       read (UnitTmp_,*) StringHeader
    end do

    do iRairden = 1, nRairden
       read(UnitTmp_,*) RairdenDistance_I(iRairden), RairdenDensity
       LnRairdenDensity_I(iRairden) = log(RairdenDensity)
    end do
    close(UnitTmp_) 

    !interpolate the hydrogen density to a new, well refined grid.

    rMax = RairdenDistance_I(nRairden)
    rMin = RairdenDistance_I(1)  
    dR = (rMax-rMin)/(nUniform-1)

    do i =1, nUniform
       Rad_I(i) = rMin + dR*(i-1)
       do iRairden =1, nRairden
          if (RairdenDistance_I(iRairden)>=Rad_I(i)) then
             iR = iRairden-1
             LnRhoH = LnRairdenDensity_I(iR)+ (Rad_I(i)-RairdenDistance_I(iR))*&
                  (LnRairdenDensity_I(iR+1) - LnRairdenDensity_I(iR))/&
                  (RairdenDistance_I(iR+1) - RairdenDistance_I(iR))
             RhoHUniform_I(i) = exp(LnRhoH)
             EXIT
          end if
       end do
    end do
   
    ! write the new interpolated values to a file
    open(unit=2, file='density.dat')
    write (2,*)'values for rairden density and interpolated density'
    write (2,*)' R, LnH, H'
    do i = 1, nUniform
       write(2,*) Rad_I(i),RhoHUniform_I(i)
    end do
    close(2)

    ! calculate distance to any point along the field line
    do iR = 1, nRadialPoint
       LatMax = acos(sqrt(1./L_I(iR)))
       write(*,*) 'LatMax', LatMax
       LatMin = -LatMax
       dLat   = (LatMax - LatMin)/(nPoint - 1)
       do iStep = 1, nPoint
          Lat = LatMin + (iStep-1)*dLat
          RadialDistance_II(iR,iStep) = abs(L_I(iR)/cos(Lat))
       end do
    end do
    
    do iR = 1, nRadialPoint
       do iStep =1, nPoint
          r = RadialDistance_II(iR,iStep)
          i = 1 + (r - rMin)/dR
          Weight = (Rad_I(i) - r)/dR
                    
          write (*,*) 'iR, iStep', iR, iStep
          write(*,*) 'i,r,weight', i, r, Weight
          
          RhoH_II(iR,iStep) = Weight*RhoHUniform_I(i+1) + (1-Weight)*RhoHUniform_I(i)
          write(*,*) ' RhoH_II', RhoH_II
       end do
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
    real                 :: InvL
    real                 :: DeltaS1, DeltaS2, b1, b2, Coeff
    integer              :: iStep, iFirst, iLast
    !-----------------------------------------------------------------------------

    iFirst = iMirror_I(1)
    iLast  = iMirror_I(2)

    AvgHDensity = 0.0

    if (iFirst > iLast) RETURN

    InvL = 1.0/L
    Coeff = InvL/sqrt(bMirror)

    DeltaS1 = abs((bMirror-bField_I(iFirst))*(dLength_I(iFirst-1))/(bField_I(iFirst-1)-bField_I(iFirst)))

    AvgHDensity= AvgHDensity + HDensity_I(iFirst)*Coeff*(2./3.)*DeltaS1*sqrt(bMirror-bField_I(iFirst))

    do iStep = iFirst, iLast-1
       b1 = bField_I(iStep)
       b2 =  bField_I(iStep+1)
       AvgHDensity = AvgHDensity + HDensity_I(iStep)*Coeff*(2./3.)*dLength_I(iStep)/(b1 - b2) &
            *( sqrt(bMirror  - b2)**3 - sqrt(bMirror  - b1)**3 )  

    end do

    DeltaS2 = abs((bMirror-bField_I(iLast))*(dLength_I(iLast))/(bField_I(iLast+1)-bField_I(iLast)))
    AvgHDensity= AvgHDensity + HDensity_I(iLast)*Coeff*(2./3.)*DeltaS2*(sqrt(bMirror-bField_I(iLast)))

  end subroutine get_hydrogen_density
  !===============================================================

  subroutine test_density

    use ModHeidiBField, ONLY: initialize_b_field,find_mirror_points
 

    integer, parameter   :: nPoint = 11
    integer, parameter   :: nPitch = 1
    integer, parameter   :: nR = 1
    integer, parameter   :: nPhi =1
    real                 :: L_I(nR)
    real                 :: Phi_I(nPhi)    
    real                 :: bFieldMagnitude_III(nPoint,nR,nPhi) 
    real                 :: Rho_II(nPoint)    
    real                 :: PitchAngle_I(nPitch) 
    real                 :: dLength_III(nPoint-1,nR,nPhi)      ! Length interval between i and i+1  
    real                 :: Length_III(nPoint,nR,nPhi) 
    real                 :: RadialDistance_III(nPoint,nR,nPhi)
    real                 :: AvgHDensity(nR),HDensity(nR,nPoint)
    real                 :: bMirror_I(nPitch),bMirror
    integer              :: iMirror_I(2)
    real, parameter      :: Pi = 3.141592654   
    real                 :: Ds_I(nPoint)
    integer              :: iStep
    !-----------------------------------------------------------------------------

    L_I =2.0
    PitchAngle_I(1) = Pi/10.
    
    call initialize_b_field (L_I, Phi_I, nPoint, nR, nPhi, bFieldMagnitude_III, RadialDistance_III,Length_III, dLength_III)
!call initialize_b_field(L_I(1), iStep, bField_I, RadialDistance_I, Length_I, Ds_I)
            
    call find_mirror_points (iStep, nPitch, PitchAngle_I,bFieldMagnitude_III, bMirror_I,iMirror_I)
    
    bMirror = bMirror_I(1)
 
    call get_rairden_density(nPoint, nR,L_I(1), Rho_II)
    call get_hydrogen_density(nPoint, L_I(1), bFieldMagnitude_III, bMirror,iMirror_I(1),dLength_III,Rho_II,AvgHDensity(1))
       
  end subroutine test_density
  
end module ModHeidiHydrogenGeo
