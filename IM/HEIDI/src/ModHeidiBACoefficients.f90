module ModHeidiBACoefficients

  implicit none

contains

subroutine get_drift_velocity(nPoint,nR,nPhi, Energy, bFieldMagnitude_III,VDrift_VIII,GradBCrossB_VIII)
  
  ! \ This subroutine calculates the drift velocity for a
  !   stretched dipole (TypeBFieldGrid = 'stretched2')
  ! /

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
  real                   :: Energy,b2,b4
  integer                :: iPoint, iR, iPhi
  !-----------------------------------------------------------------------------

  do iPhi = 1, nPhi
     do iR =1, nR
        do iPoint =1 ,nPoint

           b2 = bFieldMagnitude_III(iPoint,iR,iPhi)*bFieldMagnitude_III(iPoint,iR,iPhi)
           b4 = b2*b2
           b4_III(iPoint,iR,iPhi) = b4
                   
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

subroutine get_bounced_drift(nPoint, L, bField_I, bMirror,iMirror_I,&
     dLength_I,Drift_I,BouncedDrift,RadialDistance_I,BouncedInvRdRdt)
  
  integer, intent(in)  :: nPoint
  real,    intent(in)  :: L
  real,    intent(in)  :: bField_I(nPoint) 
  real,    intent(in)  :: bMirror  
  integer, intent(in)  :: iMirror_I(2)
  real,    intent(in)  :: dLength_I(nPoint-1)
  real,    intent(in)  :: Drift_I(nPoint)
  real,    intent(in)  :: RadialDistance_I(nPoint)
  real,    intent(out) :: BouncedDrift
  real, optional       :: BouncedInvRdRdt
  real                 :: Inv2L
  real                 :: DeltaS1, DeltaS2, b1, b2, Coeff
  integer              :: iPoint, iFirst, iLast

    !-----------------------------------------------------------------------------

    iFirst = iMirror_I(1)
    iLast  = iMirror_I(2)

    BouncedDrift    = 0.0
    BouncedInvRdRdt = 0.0

    if (iFirst > iLast) RETURN

    Coeff = sqrt(bMirror) 
    DeltaS1 = abs((bMirror-bField_I(iFirst))*&
         (dLength_I(iFirst-1))/(bField_I(iFirst-1)-bField_I(iFirst)))

   BouncedDrift    = BouncedDrift + Drift_I(iFirst)*Coeff*2.*DeltaS1/(sqrt(bMirror-bField_I(iFirst)))
   BouncedInvRdRdt =  BouncedInvRdRdt +&
        (1./RadialDistance_I(iFirst))*Drift_I(iFirst)*Coeff*2.*DeltaS1/(sqrt(bMirror-bField_I(iFirst)))
   
   do iPoint = iFirst, iLast-1
      b1 = bField_I(iPoint)
      b2 = bField_I(iPoint+1)
      
       BouncedDrift =  BouncedDrift + Drift_I(iPoint)*Coeff*2.*dLength_I(iPoint)/(b1 - b2) &
            *( sqrt(bMirror  - b2) - sqrt(bMirror  - b1) )

       
       BouncedInvRdRdt= BouncedInvRdRdt + &
            (1./RadialDistance_I(iPoint))*Drift_I(iPoint)*Coeff*2.*dLength_I(iPoint)/(b1 - b2) &
            *( sqrt(bMirror  - b2) - sqrt(bMirror  - b1) )
   
    end do

      
    DeltaS2 = abs((bMirror-bField_I(iLast))*(dLength_I(iLast))/(bField_I(iLast+1)-bField_I(iLast)))

   
   BouncedDrift  = BouncedDrift + Drift_I(iLast)* Coeff*2.*DeltaS2/(sqrt(bMirror-bField_I(iLast)))
   BouncedInvRdRdt = BouncedInvRdRdt + (1./RadialDistance_I(iLast))*Drift_I(iLast)* Coeff*2.*DeltaS2/(sqrt(bMirror-bField_I(iLast)))
 
end subroutine get_bounced_drift

!===============================================================

end module ModHeidiBACoefficients
