module ModHeidiBACoefficients

  implicit none

contains

subroutine get_drift_velocity(nPoint,nR,nPhi)
  
  ! \ This subroutine calculates the drift velocity for a
  !   stretched dipole (TypeBFieldGrid = 'stretched2')
  ! /

  use ModCoordTransform, ONLY: cross_product
  use ModConst,          ONLY: cElectronCharge

  integer, intent(in)    :: nPoint        
  integer, intent(in)    :: nR        
  integer, intent(in)    :: nPhi      
  real                   :: bField_VIII(3,nPoint,nR,nPhi)
  real                   :: GradB2over2_VIII(3,nPoint,nR,nPhi)
  real                   :: GradBCrossB_VIII(3,nPoint,nR,nPhi)
  real                   :: bFieldMagnitude_III(nPoint,nR,nPhi)
  real                   :: b4_III(nPoint,nR,nPhi)
  real                   :: VDrift_VIII(3,nPoint,nR,nPhi)
  real                   :: Energy
  integer                :: iPoint, iR, iPhi
  !-----------------------------------------------------------------------------

  do iPhi = 1, nPhi
     do iR =1, nR
        do iPoint =1 ,nPoint
           GradBCrossB_VIII(:,iPoint,iR,iPhi) = &
                cross_product(GradB2over2_VIII(:,iPoint,iR,iPhi), &
                bField_VIII(:,iPoint,iR,iPhi))
           
           b4_III(iPoint,iR,iPhi)= bFieldMagnitude_III(iPoint,iR,iPhi)*&
                bFieldMagnitude_III(iPoint,iR,iPhi)*&
                bFieldMagnitude_III(iPoint,iR,iPhi)*&
                bFieldMagnitude_III(iPoint,iR,iPhi)
                   
           VDrift_VIII(1,iPoint,iR,iPhi) = -Energy*GradBCrossB_VIII(1,iPoint,iR,iPhi)/&
                (cElectronCharge*b4_III(iPoint,iR,iPhi))   ! Vr
           VDrift_VIII(2,iPoint,iR,iPhi) = -Energy*GradBCrossB_VIII(2,iPoint,iR,iPhi)/&
                (cElectronCharge*b4_III(iPoint,iR,iPhi))   ! VTheta
           VDrift_VIII(3,iPoint,iR,iPhi) = -Energy*GradBCrossB_VIII(3,iPoint,iR,iPhi)/&
                 (cElectronCharge*b4_III(iPoint,iR,iPhi))  ! VPhi
           
        end do
     end do
  end do

end subroutine get_drift_velocity

!===============================================================

subroutine get_bounced_drift(nPoint, L, bField_I, bMirror,b4_I,&
     GradBCrossB_I,iMirror_I,dLength_I,BouncedDrift)

    integer, intent(in)  :: nPoint
    real,    intent(in)  :: L
    real,    intent(in)  :: bField_I(nPoint) 
    real,    intent(in)  :: bMirror  
    integer, intent(in)  :: iMirror_I(2)
    real,    intent(in)  :: dLength_I(nPoint-1)
    real,    intent(in)  :: b4_I(nPoint),GradBCrossB_I(nPoint) 
    real,    intent(out) :: BouncedDrift
    real                 :: DeltaS1, DeltaS2, b1, b2, Coeff
    integer              :: iPoint, iFirst, iLast
    !-----------------------------------------------------------------------------

    iFirst = iMirror_I(1)
    iLast  = iMirror_I(2)

    BouncedDrift = 0.0

    if (iFirst > iLast) RETURN

    Coeff = 1/sqrt(bMirror)

    DeltaS1 = abs((bMirror-bField_I(iFirst))*(dLength_I(iFirst-1))/(bField_I(iFirst-1)-bField_I(iFirst)))

    BouncedDrift= BouncedDrift + (1./b4_I(iFirst))*GradBCrossB_I(iFirst)*Coeff*(2./3.)*DeltaS1*sqrt(bMirror-bField_I(iFirst))

    do iPoint = iFirst, iLast-1
       b1 = bField_I(iPoint)
       b2 = bField_I(iPoint+1)
       BouncedDrift = BouncedDrift + (1./b4_I(iFirst))*GradBCrossB_I(iPoint)*Coeff*(2./3.)*dLength_I(iPoint)/(b1 - b2) &
            *( sqrt(bMirror  - b2)**3 - sqrt(bMirror  - b1)**3 )  

    end do

    DeltaS2 = abs((bMirror-bField_I(iLast))*(dLength_I(iLast))/(bField_I(iLast+1)-bField_I(iLast)))
    BouncedDrift = BouncedDrift + (1./b4_I(iLast))*GradBCrossB_I(iLast)*Coeff*(2./3.)*DeltaS2*(sqrt(bMirror-bField_I(iLast)))



end subroutine get_bounced_drift

!===============================================================
subroutine get_radial_drift(nPoint, L, bField_I, bMirror,RadialDistance,&
     RadialDrift_I,iMirror_I,dLength_I,BouncedRadialDrift,BouncedInvRdRdt)

    
  integer, intent(in)  :: nPoint
  real,    intent(in)  :: L
  real,    intent(in)  :: bField_I(nPoint), RadialDistance(nPoint)
  real,    intent(in)  :: bMirror  
  integer, intent(in)  :: iMirror_I(2)
  real,    intent(in)  :: dLength_I(nPoint-1)
  real,    intent(in)  :: RadialDrift_I(nPoint)
  real,    intent(out) :: BouncedRadialDrift,BouncedInvRdRdt
  real                 :: DeltaS1, DeltaS2, b1, b2, Coeff
  integer              :: iPoint, iFirst, iLast
  !-----------------------------------------------------------------------------
  
  iFirst = iMirror_I(1)
  iLast  = iMirror_I(2)
  
  BouncedRadialDrift = 0.0
  
    if (iFirst > iLast) RETURN

    Coeff = 1/sqrt(bMirror)

    DeltaS1 = abs((bMirror-bField_I(iFirst))*(dLength_I(iFirst-1))/(bField_I(iFirst-1)-bField_I(iFirst)))

    BouncedRadialDrift = BouncedRadialDrift + RadialDrift_I(iFirst)*Coeff*(2./3.)*DeltaS1*sqrt(bMirror-bField_I(iFirst))
    BouncedInvRdRdt    =  BouncedInvRdRdt + &
         (1./RadialDistance(iFirst))*RadialDrift_I(iFirst)*Coeff*(2./3.)*DeltaS1*sqrt(bMirror-bField_I(iFirst))
    
    do iPoint = iFirst, iLast-1
       b1 = bField_I(iPoint)
       b2 = bField_I(iPoint+1)
       BouncedRadialDrift = BouncedRadialDrift +  RadialDrift_I(iPoint)*Coeff*(2./3.)*dLength_I(iPoint)/(b1 - b2) &
            *( sqrt(bMirror  - b2)**3 - sqrt(bMirror  - b1)**3 )  
       
       BouncedInvRdRdt =  BouncedInvRdRdt + &
            (1./RadialDistance(iPoint))*RadialDrift_I(iPoint)*Coeff*(2./3.)*dLength_I(iPoint)/(b1 - b2) &
            *( sqrt(bMirror  - b2)**3 - sqrt(bMirror  - b1)**3 )  


    end do

    DeltaS2 = abs((bMirror-bField_I(iLast))*(dLength_I(iLast))/(bField_I(iLast+1)-bField_I(iLast)))
    BouncedRadialDrift = BouncedRadialDrift + RadialDrift_I(iLast)*Coeff*(2./3.)*DeltaS2*(sqrt(bMirror-bField_I(iLast)))
    BouncedInvRdRdt    =  BouncedInvRdRdt + &
         (1./RadialDistance(iLast))*RadialDrift_I(iLast)*Coeff*(2./3.)*DeltaS2*(sqrt(bMirror-bField_I(iLast)))
             

end subroutine get_radial_drift

!===============================================================
!!$  subroutine test_coefficients

!!$    use ModHeidiBField, ONLY: initialize_b_field,find_mirror_points
!!$    
!!$    integer, parameter   :: nStep =11                  
!!$    integer,parameter    :: nRadialPoint =1 
!!$    real                 :: L_I(nRadialPoint)
!!$    real                 :: Rho_II(nStep)    
!!$    integer, parameter   :: nPitch = 1
!!$    real                 :: PitchAngle_I(nPitch) 
!!$    real                 :: bField_VIII_I(nRadialPoint,nStep) 
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
!!$    call initialize_b_field(L_I(1), iStep, bField_VIII_I, RadialDistance_I, Length_I, Ds_I)
!!$            
!!$    call find_mirror_points (iStep, nPitch, PitchAngle_I, bField_VIII_I, bMirror_I,iMirror_I)
!!$    
!!$    bMirror = bMirror_I(1)
!!$ 
!!$    
!!$  end subroutine test_coefficients

end module ModHeidiBACoefficients
