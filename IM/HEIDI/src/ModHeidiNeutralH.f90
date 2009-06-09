module ModHeidiHydrogenGeo
  
  implicit none
  
contains

  subroutine get_rairden_density(nStep, nRadialPoint,L, HDensity)
    
    !use ModHeidiBField, ONLY: initialize_b_field
    !use ModIoUnit,      ONLY : io_unit_new 
   
    character (len=100) :: header
    integer, parameter  :: nPoint = 81            !number of radial point in the Rairden hydrogen file
    integer, intent(in) :: nStep                  ! number of points along the field line
    integer, intent(in) :: nRadialPoint           ! number of points in the radial direction = NR
    real, intent(in)    :: L(nRadialPoint)
    
    integer             :: iUnitRairden
    integer             :: i, iPoint, iRadialPoint,iStep  
    integer             :: iMin,iMax

    real                :: RairdenDistance(nPoint),RairdenDensity(nPoint)
    real                :: LnRairdenDensity(nPoint)
    real                :: Lat(nRadialPoint,nStep)
    
    real                :: RadialDistance(nRadialPoint, nStep)
    real                :: LnHDensity(nRadialPoint, nStep)     !interpolated ln(hydrogen density)
    real, intent(out)   :: HDensity(nRadialPoint, nStep)       !interpolated hydrogen density
    real                :: LatMax(nRadialPoint), LatMin(nRadialPoint), dLat(nRadialPoint)
    !------------------------------------------------------------------------------
    
    do iRadialPoint = 1, nRadialPoint
       LatMax(iRadialPoint) = acos(sqrt(1./L(iRadialPoint)))
       LatMin(iRadialPoint) = -LatMax(iRadialPoint)
       dLat(iRadialPoint)   = (LatMax(iRadialPoint)-LatMin(iRadialPoint))/(nStep-1)
       Lat(iRadialPoint,1)  = LatMin(iRadialPoint)
       do iStep =1, nStep
          Lat(iRadialPoint,iStep) = Lat(iRadialPoint,iStep) + dLat(iRadialPoint)
       end do
    end do


    !iUnitRairden = io_unit_new()
    ! Read the Rairden Geocorona Hydrogen density file

    iUnitRairden = 11
    open(Unit=iUnitRairden,file="RairdenHydrogenGeocorona.dat", status = 'old')
    
    do i = 1, 5
       read (iUnitRairden,*) header
    end do
    
    do iPoint = 1, nPoint
       read(iUnitRairden,*) RairdenDistance(iPoint),RairdenDensity(iPoint)
       LnRairdenDensity(iPoint)=log(RairdenDensity(iPoint))
       write(*,*)  RairdenDistance(iPoint),RairdenDensity(iPoint),LnRairdenDensity(iPoint)
    end do
    close(iUnitRairden)
    
    ! calculate distance from the center of the Earth to any point along the field line
    
    do iRadialPoint = 1, nRadialPoint
       do iStep =1, nStep
          RadialDistance(iRadialPoint,iStep) = L(iRadialPoint)/sin(Lat(iRadialPoint,iStep))
       end do
    end do
    
       
    do iRadialPoint = 1, nRadialPoint
       do iStep =1, nStep
          do iPoint = 1, nPoint
             
             if (RadialDistance(iRadialPoint,iStep)>=RairdenDistance(iPoint)) then
                iMin=iPoint-1
                !EXIT
             end if
             if (RadialDistance(iRadialPoint,iStep)<=RairdenDistance(iPoint))then   
                iMax=iPoint+1
                !EXIT
             end if
             LnHDensity(iRadialPoint,iStep) = LnRairdenDensity(iMin)+(RadialDistance(iRadialPoint,iStep)&
                  -RairdenDistance(iMin))*(LnRairdenDensity(iMax)-LnRairdenDensity(iMin))&
                  /(RairdenDistance(iMax)-RairdenDistance(iMin))
             
          end do
          HDensity(iRadialPoint,iStep) = exp(LnHDensity(iRadialPoint,iStep))
       end do
       
    end do
    !write(*,*)' HDensity', HDensity
  
  end subroutine get_rairden_density

!=================================================================================
  subroutine get_hydrogen_density(nStep, iMirror_I, bMirror, bField_I, dLength_I, L, AvgHDensity,HDensity)

    integer,intent(in)  :: nStep
    integer             :: iMirror_I(2)
    real                :: bMirror  
    real                :: dLength_I(nStep-1)
    real                :: bField_I(nStep) 
    real,intent(out)    :: AvgHDensity
    real,intent(in)     :: L
    real                :: InvL
    integer             :: iStep, iFirst, iLast
    real                :: DeltaS1, DeltaS2, b1, b2, Coeff
    real,intent(in)     :: HDensity(nStep) 
    !-----------------------------------------------------------------------------
    do iRadialPoint =1, nRadialPoint
  
    iFirst = iMirror_I(1)
    iLast  = iMirror_I(2)
    
    AvgHDensity = 0.0

    if (iFirst > iLast) RETURN

    InvL = 1.0/L
    Coeff = InvL/sqrt(bMirror)

    DeltaS1 = abs((bMirror-bField_I(iFirst))*(dLength_I(iFirst-1))/(bField_I(iFirst-1)-bField_I(iFirst)))


    AvgHDensity= AvgHDensity + HDensity(iFirst)*Coeff*(2./3.)*DeltaS1*sqrt(bMirror-bField_I(iFirst))

    do iStep = iFirst, iLast-1
       b1 = bField_I(iStep)
       b2 =  bField_I(iStep+1)
       AvgHDensity = AvgHDensity + HDensity(iStep)*Coeff*(2./3.)*dLength_I(iStep)/(b1 - b2) &
            *( sqrt(bMirror  - b2)**3 - sqrt(bMirror  - b1)**3 )  

    end do

    DeltaS2 = abs((bMirror-bField_I(iLast))*(dLength_I(iLast))/(bField_I(iLast+1)-bField_I(iLast)))
    AvgHDensity= AvgHDensity + HDensity(iLast)*Coeff*(2./3.)*DeltaS2*(sqrt(bMirror-bField_I(iLast)))


  end subroutine get_hydrogen_density
  !===============================================================
  
  subroutine test_density
    
    
    integer, parameter   :: nStep = 11,nRadialPoint=20
    real                 :: L(nRadialPoint)
    integer, parameter   :: nPitch = 1
    real                 :: PitchAngle_I(nPitch) 
    real                 :: bField_I(nstep) 
    real                 :: Length_I(nStep) 
    real                 :: RadialDistance_I(nStep) 
    real                 :: dLength_I(nStep)
    real                 :: AvgHDensity,HDensity(nStep)
    real                 :: bMirror_I(nPitch),bMirror
    integer              :: iMirror_II(2,nPitch)
    real, parameter      :: Pi = 3.141592654
    integer              :: iStep
    !----------------------------------------------------------------------------------
    L = 2.0
  
    do iStep =11, nStep
       call initialize_b_field(L, iStep, bField_I, RadialDistance_I, Length_I, dLength_I)
       PitchAngle_I(1) = Pi/10.
       call find_mirror_points (iStep, nPitch, PitchAngle_I, bField_I, bMirror_I,iMirror_II)
       bMirror = bMirror_I(1)
       call get_rairden_density(nStep, nRadialPoint,L, HDensity)
       call get_hydrogen_density(nStep, iMirror_II, bMirror, bField_I, dLength_I,L, AvgHDensity)
       
    end do
  
  end subroutine test_density
  
end module ModHeidiHydrogenGeo
