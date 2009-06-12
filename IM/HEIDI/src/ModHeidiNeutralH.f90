module ModHeidiHydrogenGeo

  implicit none

contains

  subroutine get_rairden_density(nStep, nRadialPoint,L_I, HDensity_II)

    !use ModIoUnit,      ONLY : io_unit_new 

    integer, intent(in)  :: nStep                  ! number of points along the field line
    integer, intent(in)  :: nRadialPoint           ! number of points in the radial direction = NR
    real,    intent(in)  :: L_I(nRadialPoint)      ! x^2+y^2
    real,    intent(out) :: HDensity_II(nRadialPoint, nStep)  !interpolated density
    integer, parameter   :: n =100                 !number of points on new refined grid
    integer, parameter   :: nPoint = 81            !number of radial point in the Rairden hydrogen file
    real,    parameter   :: epsilon=1e-30    
    character (len=100)  :: header
    real                 :: RairdenDistance_I(nPoint),RairdenDensity_I(nPoint)
    real                 :: LnRairdenDensity_I(nPoint)
    real                 :: Lat_II(nRadialPoint,nStep)
    real                 :: RadialDistance_II(nRadialPoint, nStep)
    real                 :: LatMax_I(nRadialPoint), LatMin_I(nRadialPoint), dLat_I(nRadialPoint)
    real                 :: Rad_I(n), weight_I(n),H_I(n),LnH_I(n)
    integer              :: iPoint, i
    integer              :: iUnitRairden
    real                 :: rMax, rMin, dR
    integer              :: iRadialPoint,iStep,iR  
    !------------------------------------------------------------------------------

    ! Read the Rairden Geocorona Hydrogen density file.

    !iUnitRairden = io_unit_new()
    iUnitRairden = 11
    open(Unit=iUnitRairden,file="RairdenHydrogenGeocorona.dat", status = 'old')

    do i = 1, 5
       read (iUnitRairden,*) header
    end do

    do iPoint = 1, nPoint
       read(iUnitRairden,*) RairdenDistance_I(iPoint),RairdenDensity_I(iPoint)
       LnRairdenDensity_I(iPoint)=log(RairdenDensity_I(iPoint))
    end do
    close(iUnitRairden) 

    !interpolate the hydrogen density to new, well refined grid.

    rMax = RairdenDistance_I(81) 
    rMin = RairdenDistance_I(1)  
    dR = (rMax-rMin)/(n-1)

    Rad_I(1) = Rmin

    do i =1, n
       Rad_I(i+1) = Rad_I(i) + dR
       weight_I(i) = Rad_I(i)/dR
       do iPoint =1, nPoint
          if (RairdenDistance_I(iPoint)>=Rad_I(i)) then
             iR = iPoint
             LnH_I(i) = LnRairdenDensity_I(iR)+ (Rad_I(i)-RairdenDistance_I(iR))*&
                  (LnRairdenDensity_I(iR+1)-LnRairdenDensity_I(iR))/&
                  (RairdenDistance_I(iR+1)-RairdenDistance_I(iR))
             EXIT
          end if
       end do
       H_I(i) = exp(LnH_I(i))
    end do

    open(unit=2, file='density.dat')
    write (2,*)'values for rairden density and interpolated density'
    write (2,*)' R, LnH, H'

    do i=1,n
       write(2,*) Rad_I(i), LnH_I(i),H_I(i)
    end do
    close(2)

    do iRadialPoint = 1, nRadialPoint
       LatMax_I(iRadialPoint) = acos(sqrt(1./L_I(iRadialPoint)))
       LatMin_I(iRadialPoint) = -LatMax_I(iRadialPoint)
       dLat_I(iRadialPoint)   = (LatMax_I(iRadialPoint)-LatMin_I(iRadialPoint))/(nStep-1)
       do iStep =2, nStep
          Lat_II(iRadialPoint,1) = LatMin_I(iRadialPoint) 
          Lat_II(iRadialPoint,iStep) = Lat_II(iRadialPoint,iStep-1)+ dLat_I(iRadialPoint)
       end do
    end do

    ! calculate distance to any point along the field line

    do iRadialPoint = 1, nRadialPoint
       do iStep =1, nStep
          RadialDistance_II(iRadialPoint,iStep) = abs(L_I(iRadialPoint)/(sin(Lat_II(iRadialPoint,iStep))+epsilon))
       end do
    end do

    do iRadialPoint = 1, nRadialPoint
       do iStep =1, nStep
          do i = 1, n
             if (Rad_I(i)>=RadialDistance_II(iRadialPoint,iStep)) then
                iR = i
                EXIT
             end if
          end do
          !HDensity_II(iRadialPoint,iStep) = weight_I(i)*H_I(iR)-weight_I(i+1)*H_I(iR+1)
          HDensity_II(iRadialPoint,iStep) = H_I(iR)+(RadialDistance_II(iRadialPoint,nStep)-Rad_I(iR))*&
               (H_I(iR+1)-H_I(iR))/dR
       end do
    end do

  end subroutine get_rairden_density

  !=================================================================================
  subroutine get_hydrogen_density(nStep, L, bField_I, bMirror,iMirror_I,dLength_I,HDensity_I,AvgHDensity)

    integer, intent(in)  :: nStep
    real,    intent(in)  :: L
    real,    intent(in)  :: bField_I(nStep) 
    real,    intent(in)  :: bMirror  
    integer, intent(in)  :: iMirror_I(2)
    real,    intent(in)  :: dLength_I(nStep-1)
    real,    intent(in)  :: HDensity_I(nStep) 
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
    
    integer, parameter  :: nStep =11                  
    integer,parameter   :: nRadialPoint =1 
    real                :: L_I(nRadialPoint)
    real                :: HDensity_II(nRadialPoint, nStep)    

    integer, parameter   :: nPitch = 1
    real                 :: PitchAngle_I(nPitch) 
    real                 :: bField_I(nRadialPoint,nStep) 
    real                 :: Length_I(nRadialPoint,nStep) 
    real                 :: RadialDistance_I(nRadialPoint,nStep) 
    real                 :: dLength_I(nStep)
    real                 :: AvgHDensity(nRadialPoint),HDensity(nRadialPoint,nStep)
    real                 :: bMirror_I(nRadialPoint,nPitch),bMirror(nRadialPoint)
    integer              :: iMirror_II(nRadialPoint,2,nPitch)
    real, parameter      :: Pi = 3.141592654   
    !-----------------------------------------------------------------------------

    L_I =2.0
    PitchAngle_I(1) = Pi/10.
    
    call initialize_b_field(L, iStep, bField_I, RadialDistance_I, Length_I, Ds_I)
            
    call find_mirror_points (iStep, nPitch, PitchAngle_I, bField_I, bMirror_I,iMirror_II)
    
    bMirror = bMirror_I(1)
 
    call get_rairden_density(nStep, nRadialPoint,L_I, HDensity_II)
    call get_hydrogen_density(nStep, L, bField_I, bMirror,iMirror_I,dLength_I,HDensity_I,AvgHDensity)
       
  end subroutine test_density
  !===============================================================
end module ModHeidiHydrogenGeo
