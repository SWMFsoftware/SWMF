module ModHeidiBACoefficients
  
  implicit none
  
contains
  
  subroutine get_mudot
    !Calculates the bounce-average rate of change of the pitch angle
    
    integer,parameter   :: nStep  =11  
    real                :: RadialDistance(nStep) 
    integer             :: iFirst, iLast, iStep
    real                :: dRdt(nStep), dBdt(nStep)
    real                :: bField_I(nStep) 
    real                :: bMirror  
    integer             :: iMirror_I(2)
    real                :: dLength_I(nStep-1)
    real                :: CoeffRadial, CoeffMagnetic
    real                :: InvR(nStep), InvB(nStep)
    real                :: DeltaS1, DeltaS2, b1, b2, Coeff
    real                :: FactorR(nStep), FactorB(nStep)
    
    !-----------------------------------------------------------------------------
    iFirst = iMirror_I(1)
    iLast  = iMirror_I(2)
    
    ! \ 
    !   Calculate the bounced average of FactorR = (1/R)*dR/dt
    !   Calculate the bounced average of FactorB = (1/B)*dB/dt
    ! / 
    
    
    do iStep =1, nStep
       InvR(iStep)    = 1.0/RadialDistance(iStep)
       InvB(iStep)    = 1.0/bField_I(iStep)
       FactorR(iStep) = InvR(iStep)*dRdt(iStep)
       FactorB(iStep) = InvB(iStep)*dBdt(iStep)
    end do
    
    Coeff = 1./sqrt(bMirror)
    
    CoeffRadial   = 0.0
    CoeffMagnetic = 0.0
    
    if (iFirst > iLast) RETURN
    
    DeltaS1 = abs((bMirror-bField_I(iFirst))*(dLength_I(iFirst-1))/(bField_I(iFirst-1)-bField_I(iFirst)))

    CoeffRadial   = CoeffRadial + FactorR(iFirst)*Coeff*(2./3.)*DeltaS1*sqrt(bMirror-bField_I(iFirst))
    CoeffMagnetic = CoeffMagnetic + FactorB(iFirst)*Coeff*(2./3.)*DeltaS1*sqrt(bMirror-bField_I(iFirst))

    do iStep = iFirst, iLast-1
       b1 =  bField_I(iStep)
       b2 =  bField_I(iStep+1)
       CoeffRadial = CoeffRadial + FactorR(iStep)*Coeff*(2./3.)*dLength_I(iStep)/(b1 - b2) &
            *( sqrt(bMirror  - b2)**3 - sqrt(bMirror  - b1)**3 ) 
       CoeffMagnetic = CoeffMagnetic + FactorB(iStep)*Coeff*(2./3.)*dLength_I(iStep)/(b1 - b2) &
            *( sqrt(bMirror  - b2)**3 - sqrt(bMirror  - b1)**3 ) 

    end do

    DeltaS2 = abs((bMirror-bField_I(iLast))*(dLength_I(iLast))/(bField_I(iLast+1)-bField_I(iLast)))
    CoeffRadial   = CoeffRadial + FactorR(iLast)*Coeff*(2./3.)*DeltaS2*(sqrt(bMirror-bField_I(iLast)))
    CoeffMagnetic = CoeffMagnetic + FactorB(iLast)*Coeff*(2./3.)*DeltaS2*(sqrt(bMirror-bField_I(iLast)))

    

  end subroutine get_mudot
  !===============================================================
  subroutine get_drift_velocity
    ! \ This subroutine calculates the drift velocity for an 
    !   arbitrary magnetic field. 
    ! /
    
    use  ModCoordTransform
    
    implicit none
    
    integer,parameter :: nR = 20,nLat = 20,nLon =24
    !integer   :: nR, nLat, nLon
    real      :: BField(3, nR, nLat, nLon) 
    real      :: RadialDist(nR)
    real      :: Lat(nLat)
    real      :: Lon(nLon)
    real      :: Bsquared(3, nR, nLat, nLon)
    real      :: gradB(3,nR, nLat, nLon)
    real      :: B1(nR,nLat,nLon), B2(nR,nLat,nLon), B3(nR,nLat,nLon)
    real      :: GradBCrossB(nR, nLat, nLon)
   
    real      :: gradB1(nR,nLat,nLon), gradB2(nR,nLat,nLon), gradB3(nR,nLat,nLon)
    integer   :: iR, iLat, iLon
    
   !-----------------------------------------------------------------------------
   ! Define the half square of the magnetic field.
   do iR = 1, nR
      do iLat =1, nLat
         do iLon =1, nLon
            Bsquared(:,iR,ilat,iLon) = 0.5*BField(:,iR,ilat,iLon)* &
                 BField(:,iR,ilat,iLon)
         end do
      end do
   end do

   ! Define the gradient of (B^2)/2

   B1 = BField(1,iR,iLat,iLon)
   B2 = BField(2,iR,iLat,iLon)
   B3 = BField(3,iR,iLat,iLon)


   do iR = 2, nR -1
      do iLat = 2, nLat-1
         do iLon = 2, nLon-1
            if (iR == 1 .or. iR == nR)         gradB1(iR, iLat, iLon) = 0.0
            if (iLat ==1 .or. iLat == nLat)    gradB2(iR, iLat, iLon) = 0.0
            if (iLon == 1 .or. iLon == nLon)   gradB3(iR, iLat, iLon) = 0.0
        
            
            gradB1(iR, iLat, iLon) = (B1(iR+1,iLat,iLon)-B1(iR-1,iLat,iLon))&
                 /(RadialDist(iR+1)-RadialDist(iR-1))
            
            gradB2(iR, iLat, iLon) = (1./RadialDist(iR))*(B2(iR,iLat+1,iLon)-B2(iR,iLat-1,iLon))&
                 /(Lat(iLat+1)-Lat(iLat-1))
            
            gradB3(iR, iLat, iLon) = (1./(RadialDist(iR)*sin(Lat(iLat))))&
                 *(B3(iR,iLat,iLon+1)-B3(iR,iLat,iLon-1))/(Lon(iLon+1)-Lon(iLon-1))
         
            gradB(1,iR, iLat, iLon) = gradB1(iR, iLat, iLon)
            gradB(2,iR, iLat, iLon) = gradB2(iR, iLat, iLon)
            gradB(3,iR, iLat, iLon) = gradB3(iR, iLat, iLon)
            
            ! Calculate the grad(B^2)xB
            
            GradBCrossB(iR, iLat, iLon) =  &
                 cross_product33(gradB1(iR, iLat, iLon),&
                 gradB2(iR, iLat, iLon),&
                 gradB3(iR, iLat, iLon))
            
            
         end do
      end do
   end do

   ! Calculate the grad(B^2)xB

   !function cross_product33(aX, aY, aZ, bX, bY, bZ) result(c_D)
   
    
  end subroutine get_drift_velocity


 !=============================================================== 



  subroutine test_coefficients

!!$    use ModHeidiBField, ONLY: initialize_b_field,find_mirror_points
!!$    
!!$    integer, parameter   :: nStep =11                  
!!$    integer,parameter    :: nRadialPoint =1 
!!$    real                 :: L_I(nRadialPoint)
!!$    real                 :: Rho_II(nStep)    
!!$    integer, parameter   :: nPitch = 1
!!$    real                 :: PitchAngle_I(nPitch) 
!!$    real                 :: bField_I(nRadialPoint,nStep) 
!!$    real                 :: Length_I(nRadialPoint,nStep) 
!!$    real                 :: RadialDistance_I(nRadialPoint,nStep) 
!!$    real                 :: dLength_I(nStep)
!!$    real                 :: bMirror_I(nPitch),bMirror
!!$    integer              :: iMirror_I(2)
!!$    real, parameter      :: Pi = 3.141592654   
!!$    real                 :: Ds_I(nStep)
!!$    integer              :: iStep
!!$    !-----------------------------------------------------------------------------
!!$
!!$    L_I =2.0
!!$    PitchAngle_I(1) = Pi/10.
!!$    
!!$    call initialize_b_field(L_I(1), iStep, bField_I, RadialDistance_I, Length_I, Ds_I)
!!$            
!!$    call find_mirror_points (iStep, nPitch, PitchAngle_I, bField_I, bMirror_I,iMirror_I)
!!$    
!!$    bMirror = bMirror_I(1)
!!$ 
!!$    
  end subroutine test_coefficients
  
end module ModHeidiBACoefficients
