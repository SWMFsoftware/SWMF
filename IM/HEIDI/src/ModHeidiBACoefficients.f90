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
                   
           VDrift_VIII(1,iPoint,iR,iPhi) =  Energy*GradBCrossB_VIII(1,iPoint,iR,iPhi)/&
                (b4_III(iPoint,iR,iPhi))   ! Vr                           
           VDrift_VIII(2,iPoint,iR,iPhi) =  Energy*GradBCrossB_VIII(2,iPoint,iR,iPhi)/&
                (b4_III(iPoint,iR,iPhi))   ! VTheta                       
           VDrift_VIII(3,iPoint,iR,iPhi) =  Energy*GradBCrossB_VIII(3,iPoint,iR,iPhi)/&
                (b4_III(iPoint,iR,iPhi))  ! VPhi 

           
        end do
     end do
  end do


end subroutine get_drift_velocity

!===============================================================

subroutine get_bounced_drift(nPoint, L, bField_I, bMirror,iMirror_I,&
     dLength_I,Drift_I,BouncedDrift,RadialDistance_I)
  
  integer, intent(in)  :: nPoint
  real,    intent(in)  :: L
  real,    intent(in)  :: bField_I(nPoint) 
  real,    intent(in)  :: bMirror  
  integer, intent(in)  :: iMirror_I(2)
  real,    intent(in)  :: dLength_I(nPoint-1)
  real,    intent(in)  :: Drift_I(nPoint)
  real,    intent(in)  :: RadialDistance_I(nPoint)
  real,    intent(out) :: BouncedDrift
  real                 :: Inv2L
  real                 :: DeltaS1, DeltaS2, b1, b2, Coeff
  integer              :: iPoint, iFirst, iLast

  !-----------------------------------------------------------------------------
  
  iFirst = iMirror_I(1)
  iLast  = iMirror_I(2)
  
  BouncedDrift    = 0.0
  
  if (iFirst > iLast) RETURN
  
  Inv2L = 1./(2*L)
  Coeff = Inv2L*sqrt(bMirror) 
  DeltaS1 = abs((bMirror-bField_I(iFirst))*&
       (dLength_I(iFirst-1))/(bField_I(iFirst-1)-bField_I(iFirst)))
  
  BouncedDrift    = BouncedDrift + Drift_I(iFirst)*Coeff*2.*DeltaS1/(sqrt(bMirror-bField_I(iFirst)))
  
  do iPoint = iFirst, iLast-1
     b1 = bField_I(iPoint)
     b2 = bField_I(iPoint+1)
     
     BouncedDrift =  BouncedDrift + Drift_I(iPoint)*Coeff*2.*dLength_I(iPoint)/(b1 - b2) &
          *( sqrt(bMirror  - b2) - sqrt(bMirror  - b1) )

  end do
  
  
  DeltaS2 = abs((bMirror-bField_I(iLast))*(dLength_I(iLast))/(bField_I(iLast+1)-bField_I(iLast)))
  
  
  BouncedDrift  = BouncedDrift + Drift_I(iLast)* Coeff*2.*DeltaS2/(sqrt(bMirror-bField_I(iLast)))

end subroutine get_bounced_drift

!===============================================================

subroutine get_bounced_VGradB(nPoint, L, bField_I, GradB_I,bMirror,iMirror_I,&
     dLength_I,Drift_I,RadialDistance_I,BouncedVGradB)
  
  integer, intent(in)  :: nPoint
  real,    intent(in)  :: L
  real,    intent(in)  :: bField_I(nPoint) 
  real,    intent(in)  :: bMirror  
  integer, intent(in)  :: iMirror_I(2)
  real,    intent(in)  :: dLength_I(nPoint-1)
  real,    intent(in)  :: Drift_I(nPoint), GradB_I(nPoint)
  real,    intent(in)  :: RadialDistance_I(nPoint)
  real,    intent(out) :: BouncedVGradB
  real                 :: DeltaS1, DeltaS2, b1, b2, Coeff
  real                 :: Temp, TempFirst,TempLast
  integer              :: iPoint, iFirst, iLast

    !-----------------------------------------------------------------------------

    iFirst = iMirror_I(1)
    iLast  = iMirror_I(2)

    BouncedVGradB    = 0.0
    if (iFirst > iLast) RETURN
    
    Coeff = sqrt(bMirror) 
    DeltaS1 = abs((bMirror-bField_I(iFirst))*&
         (dLength_I(iFirst-1))/(bField_I(iFirst-1)-bField_I(iFirst)))
    
   TempFirst = (1./bField_I(iFirst)) * Drift_I(iFirst) * GradB_I(iFirst)
!    TempFirst = Drift_I(iFirst) * GradB_I(iFirst)
  BouncedVGradB    = BouncedVGradB + TempFirst*Coeff*2.*DeltaS1/(sqrt(bMirror-bField_I(iFirst)))
   
   
   do iPoint = iFirst, iLast-1
      b1 = bField_I(iPoint)
      b2 = bField_I(iPoint+1)
      Temp = (1./bField_I(iPoint)) * Drift_I(iPoint) * GradB_I(iPoint)
!      Temp =  Drift_I(iPoint) * GradB_I(iPoint)      

      BouncedVGradB  =  BouncedVGradB  + Temp *Coeff*2.*dLength_I(iPoint)/(b1 - b2) &
            *( sqrt(bMirror  - b2) - sqrt(bMirror  - b1) )

    end do

      
    DeltaS2 = abs((bMirror-bField_I(iLast))*(dLength_I(iLast))/(bField_I(iLast+1)-bField_I(iLast)))

   TempLast = (1./bField_I(iLast)) * Drift_I(iLast) * GradB_I(iLast)
!   TempLast = Drift_I(iLast) * GradB_I(iLast)
   BouncedVGradB   = BouncedVGradB + TempLast* Coeff*2.*DeltaS2/(sqrt(bMirror-bField_I(iLast)))
   

 end subroutine get_bounced_VGradB

!===============================================================

end module ModHeidiBACoefficients
